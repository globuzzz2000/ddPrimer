#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File preparation module for ddPrimer pipeline.

Contains functionality for:
1. VCF file validation and preparation (bgzip, indexing, normalization)
2. FASTA file indexing and validation
3. INFO/AF field addition to VCF files when missing
4. Chromosome name mapping and harmonization across files
5. Interactive file preparation with user consent

This module ensures all input files are properly formatted and compatible
before proceeding with the ddPrimer pipeline, minimizing downstream errors
and improving processing reliability.
"""

import os
import subprocess
import tempfile
import gzip
import logging
import shutil
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Set
from ..config import Config, FileError, ExternalToolError, SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)


class FilePreparator:
    """
    Handles comprehensive file preparation and validation for ddPrimer pipeline.
    
    This class validates and prepares VCF, FASTA, and GFF files to ensure
    compatibility and proper formatting. It performs automatic fixes when
    possible and requests user permission for file modifications.
    
    Attributes:
        temp_dir: Temporary directory for intermediate files
        chromosome_mapper: ChromosomeMapper instance for compatibility analysis
        prepared_files: Dictionary tracking prepared file paths
        
    Example:
        >>> preparator = FilePreparator()
        >>> result = preparator.prepare_files("input.vcf", "genome.fasta", "annotations.gff")
        >>> if result['success']:
        ...     # Use result['vcf_file'], result['fasta_file'], result['gff_file']
        ...     pass
    """
    
    def __init__(self):
        """
        Initialize file preparator.
        
        Creates temporary directory and sets up for file preparation and
        compatibility analysis between input files.
        
        Raises:
            FileError: If temporary directory creation fails
            ExternalToolError: If required tools are not available
        """
        logger.debug("=== FILE PREPARATOR INITIALIZATION DEBUG ===")
        
        try:
            self.temp_dir = Config.get_temp_dir()
            self.prepared_files = {}
            
            # Validate required dependencies
            self._validate_dependencies()
            
            logger.debug(f"Initialized FilePreparator with temp dir: {self.temp_dir}")
            logger.debug("=== END FILE PREPARATOR INITIALIZATION DEBUG ===")
            
        except Exception as e:
            error_msg = f"Failed to initialize FilePreparator: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END FILE PREPARATOR INITIALIZATION DEBUG ===")
            raise FileError(error_msg) from e
    
    def _validate_dependencies(self):
        """
        Validate that all required external tools are available.
        
        Raises:
            ExternalToolError: If required tools are missing
        """
        required_tools = [
            ('bcftools', 'bcftools --version'),
            ('bgzip', 'bgzip --version'),
            ('tabix', 'tabix --version'),
            ('samtools', 'samtools --version')
        ]
        
        missing_tools = []
        
        for tool_name, version_cmd in required_tools:
            try:
                result = subprocess.run(version_cmd.split(), 
                                      capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    logger.debug(f"{tool_name} is available")
                else:
                    missing_tools.append(tool_name)
            except (FileNotFoundError, subprocess.TimeoutExpired):
                missing_tools.append(tool_name)
        
        if missing_tools:
            error_msg = (
                f"Required tools not found: {', '.join(missing_tools)}\n"
                "Please install missing tools:\n"
                "  Ubuntu/Debian: sudo apt-get install bcftools bgzip tabix samtools\n"
                "  macOS: brew install bcftools htslib samtools\n"
                "  Conda: conda install -c bioconda bcftools htslib samtools"
            )
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name=missing_tools[0])
    
    def prepare_files(self, vcf_file: str, fasta_file: str, 
                     gff_file: Optional[str] = None) -> Dict:
        """
        Main entry point for file preparation workflow.
        
        Analyzes all input files for compatibility and formatting issues,
        then prepares corrected versions with user consent when needed.
        
        Args:
            vcf_file: Path to VCF file
            fasta_file: Path to FASTA file  
            gff_file: Optional path to GFF file
            
        Returns:
            Dictionary with preparation results and final file paths
            
        Raises:
            FileError: If input files cannot be accessed
            ExternalToolError: If external tools fail
        """
        logger.debug("=== FILE PREPARATION WORKFLOW DEBUG ===")
        logger.info("\nAnalyzing input files for compatibility and formatting...")
        
        try:
            # Validate input files exist
            self._validate_input_files(vcf_file, fasta_file, gff_file)
            
            # Analyze files and determine what needs to be prepared
            analysis = self._analyze_files(vcf_file, fasta_file, gff_file)
            
            # Check if any preparation is needed
            if not analysis['needs_preparation']:
                logger.info("All files are properly formatted and compatible!")
                return {
                    'success': True,
                    'vcf_file': vcf_file,
                    'fasta_file': fasta_file,
                    'gff_file': gff_file,
                    'changes_made': False,
                    'issues_found': analysis['issues']
                }
            
            # Report issues found
            self._report_issues(analysis['issues'])
            
            # Ask user for permission to create corrected files
            if not self._get_user_consent(analysis['issues']):
                logger.info(f"\nFile preparation canceled by user.")
                return {
                    'success': False,
                    'reason': 'User declined file preparation',
                    'issues_found': analysis['issues']
                }
            
            # Prepare files
            logger.info("\nCorrecting files...")
            prepared_files = self._prepare_corrected_files(vcf_file, fasta_file, gff_file, analysis)
            
            logger.info("File preparation completed successfully!")
            logger.debug("=== END FILE PREPARATION WORKFLOW DEBUG ===")
            
            return {
                'success': True,
                'vcf_file': prepared_files.get('vcf', vcf_file),
                'fasta_file': prepared_files.get('fasta', fasta_file),
                'gff_file': prepared_files.get('gff', gff_file),
                'changes_made': True,
                'issues_found': analysis['issues'],
                'prepared_files': prepared_files
            }
            
        except Exception as e:
            error_msg = f"Error in file preparation workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END FILE PREPARATION WORKFLOW DEBUG ===")
            raise FileError(error_msg) from e
    
    def _validate_input_files(self, vcf_file: str, fasta_file: str, 
                            gff_file: Optional[str] = None):
        """
        Validate that input files exist and are accessible.
        
        Args:
            vcf_file: Path to VCF file
            fasta_file: Path to FASTA file
            gff_file: Optional path to GFF file
            
        Raises:
            FileError: If files cannot be accessed
        """
        files_to_check = [
            (vcf_file, "VCF"),
            (fasta_file, "FASTA")
        ]
        
        if gff_file:
            files_to_check.append((gff_file, "GFF"))
        
        for file_path, file_type in files_to_check:
            if not os.path.exists(file_path):
                error_msg = f"{file_type} file not found: {file_path}"
                logger.error(error_msg)
                raise FileError(error_msg)
            
            if not os.access(file_path, os.R_OK):
                error_msg = f"{file_type} file not readable: {file_path}"
                logger.error(error_msg)
                raise FileError(error_msg)
    
    def _analyze_files(self, vcf_file: str, fasta_file: str, 
                      gff_file: Optional[str] = None) -> Dict:
        """
        Analyze all files to determine what preparation is needed.
        
        Args:
            vcf_file: Path to VCF file
            fasta_file: Path to FASTA file
            gff_file: Optional path to GFF file
            
        Returns:
            Dictionary containing analysis results and preparation requirements
        """
        logger.debug("=== FILE ANALYSIS DEBUG ===")
        
        issues = []
        
        # Analyze VCF file
        vcf_issues = self._analyze_vcf_file(vcf_file)
        issues.extend(vcf_issues)
        
        # Analyze FASTA file
        fasta_issues = self._analyze_fasta_file(fasta_file)
        issues.extend(fasta_issues)
        
        # Analyze GFF file if provided
        if gff_file:
            gff_issues = self._analyze_gff_file(gff_file)
            issues.extend(gff_issues)
        
        # Check chromosome compatibility
        compat_issues = self._analyze_chromosome_compatibility(vcf_file, fasta_file, gff_file)
        issues.extend(compat_issues)
        
        needs_preparation = len(issues) > 0
        
        logger.debug(f"Analysis complete: {len(issues)} issues found")
        logger.debug("=== END FILE ANALYSIS DEBUG ===")
        
        return {
            'needs_preparation': needs_preparation,
            'issues': issues
        }
    
    def _get_vcf_chromosomes(self, vcf_file: str) -> Dict[str, int]:
        """
        Extract chromosome names and variant counts from VCF file.
        
        Uses bcftools to efficiently extract chromosome information from
        VCF files, supporting both compressed and uncompressed formats.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            Dictionary mapping chromosome names to variant counts
            
        Raises:
            ExternalToolError: If bcftools execution fails
        """
        logger.debug(f"Extracting chromosome information from VCF: {vcf_file}")
        
        try:
            cmd = ['bcftools', 'query', '-f', '%CHROM\n', vcf_file]
            logger.debug(f"Running command: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                chromosomes = result.stdout.strip().split('\n')
                chrom_counts = defaultdict(int)
                
                for chrom in chromosomes:
                    if chrom and chrom.strip():  # Skip empty lines
                        chrom_counts[chrom.strip()] += 1
                
                logger.debug(f"Found {len(chrom_counts)} unique chromosomes in VCF")
                return dict(chrom_counts)
            else:
                error_msg = f"bcftools execution failed for {vcf_file}: {result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="bcftools")
                
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error reading VCF chromosomes from {vcf_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="bcftools") from e
    
    def _get_fasta_sequences(self, fasta_file: str) -> Dict[str, int]:
        """
        Extract sequence names and lengths from FASTA file.
        
        Parses FASTA files to extract sequence identifiers and calculate
        sequence lengths, supporting both compressed and uncompressed formats.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Dictionary mapping sequence names to lengths
            
        Raises:
            FileError: If FASTA file cannot be parsed
        """
        logger.debug(f"Extracting sequence information from FASTA: {fasta_file}")
        
        try:
            sequences = {}
            current_seq = None
            current_length = 0
            
            # Determine file opener based on extension
            opener = gzip.open if fasta_file.endswith('.gz') else open
            mode = 'rt' if fasta_file.endswith('.gz') else 'r'
                
            with opener(fasta_file, mode) as f:
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_seq is not None:
                            sequences[current_seq] = current_length
                        
                        # Start new sequence
                        seq_id = line[1:].split()[0]  # Take first part after >
                        current_seq = seq_id
                        current_length = 0
                    elif line:
                        # Count sequence characters
                        current_length += len(line)
                
                # Don't forget the last sequence
                if current_seq is not None:
                    sequences[current_seq] = current_length
            
            logger.debug(f"Found {len(sequences)} sequences in FASTA")
            return sequences
            
        except Exception as e:
            error_msg = f"Error reading FASTA sequences from {fasta_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
    
    def _extract_numeric_component(self, name: str) -> float:
        """
        Extract numeric component from chromosome name for intelligent sorting.
        
        Analyzes chromosome names to extract numeric components for proper
        sorting, handling various naming conventions including Chr1, Chromosome1,
        accession numbers, and special chromosomes.
        
        Args:
            name: Chromosome or sequence name
            
        Returns:
            Numeric component if found, large number for special/unrecognized names
        """
        if not name or not isinstance(name, str):
            return float('inf')
            
        name_upper = name.upper()
        
        # Pattern 1: Simple numbers (1, 2, 3, etc.)
        if name.isdigit():
            return int(name)
        
        # Pattern 2: Chr1, Chr2, etc.
        if name_upper.startswith('CHR') and len(name) > 3:
            chr_part = name[3:]
            if chr_part.isdigit():
                return int(chr_part)
        
        # Pattern 3: Chromosome1, Chromosome2, etc.
        if name_upper.startswith('CHROMOSOME') and len(name) > 10:
            chr_part = name[10:]
            if chr_part.isdigit():
                return int(chr_part)
        
        # Pattern 4: For accession numbers like CP002684.1, extract the main number
        base_name = name.split('.')[0] if '.' in name else name
        numbers = re.findall(r'\d+', base_name)
        if numbers:
            # Use the last/longest number found (usually the main identifier)
            return int(numbers[-1])
        
        # Pattern 5: Special chromosomes (X, Y, MT, etc.)
        special_chroms = {
            'X': 100, 'Y': 101, 'MT': 102, 'MITO': 102, 'MITOCHONDRIAL': 102,
            'CHLOROPLAST': 103, 'PLASTID': 103, 'PT': 103, 'CP': 103
        }
        
        for special, value in special_chroms.items():
            if special in name_upper:
                return value
        
        return float('inf')  # Put unrecognized names at end
    
    def _filter_nuclear_chromosomes(self, fasta_seqs: Dict[str, int]) -> Dict[str, int]:
        """
        Filter FASTA sequences to keep only likely nuclear chromosomes.
        
        Identifies and filters out organellar genomes (mitochondrial, chloroplast)
        and small sequences that are unlikely to be nuclear chromosomes.
        
        Args:
            fasta_seqs: Dictionary mapping sequence names to lengths
            
        Returns:
            Filtered dictionary with only nuclear chromosomes
        """
        nuclear_seqs = {}
        
        # Define size thresholds (organellar genomes are typically much smaller)
        MIN_NUCLEAR_SIZE = 1_000_000  # 1 Mb minimum for nuclear chromosomes
        
        for seq_name, length in fasta_seqs.items():
            seq_upper = seq_name.upper()
            
            # Check for organellar indicators in the name
            organellar_indicators = [
                'MT', 'MITO', 'MITOCHONDRIAL', 'MITOCHONDRION',
                'PT', 'PLASTID', 'CHLOROPLAST', 'CHLORO',
                'PLASMID', 'PLAS'
            ]
            
            is_organellar = any(indicator in seq_upper for indicator in organellar_indicators)
            is_too_small = length < MIN_NUCLEAR_SIZE
            
            if not is_organellar and not is_too_small:
                nuclear_seqs[seq_name] = length
                logger.debug(f"Including nuclear sequence: {seq_name} ({length:,} bp)")
            else:
                if is_organellar:
                    logger.debug(f"Excluding organellar sequence: {seq_name} ({length:,} bp)")
                elif is_too_small:
                    logger.debug(f"Excluding small sequence: {seq_name} ({length:,} bp)")
        
        logger.debug(f"Filtered to {len(nuclear_seqs)} nuclear sequences from {len(fasta_seqs)} total")
        return nuclear_seqs
    
    def _suggest_chromosome_mapping(self, vcf_chroms: Dict[str, int], 
                                   fasta_seqs: Dict[str, int]) -> Dict[str, str]:
        """
        Suggest intelligent mapping based on analysis of both files.
        
        Analyzes chromosome naming patterns and sequence characteristics
        to suggest automatic mappings between VCF chromosomes and FASTA
        sequences based on numerical ordering and sequence properties.
        
        Args:
            vcf_chroms: VCF chromosomes and their variant counts
            fasta_seqs: FASTA sequences and their lengths
            
        Returns:
            Suggested mapping dictionary from VCF chromosome to FASTA sequence
        """
        logger.debug("Analyzing files for automatic chromosome mapping")
        
        # Sort VCF chromosomes by numeric component, then alphabetically
        vcf_sorted = sorted(vcf_chroms.keys(), 
                           key=lambda x: (self._extract_numeric_component(x), x))
        
        # Filter FASTA sequences to exclude likely organellar genomes
        nuclear_fasta = self._filter_nuclear_chromosomes(fasta_seqs)
        
        # Sort nuclear FASTA sequences by numeric component, then alphabetically
        fasta_sorted = sorted(nuclear_fasta.keys(), 
                             key=lambda x: (self._extract_numeric_component(x), x))
        
        logger.debug(f"VCF chromosomes (sorted): {vcf_sorted}")
        logger.debug(f"Nuclear FASTA sequences (sorted): {fasta_sorted}")
        
        # Main strategy: Map chromosomes in order after sorting both by numeric component
        if len(vcf_sorted) <= len(fasta_sorted):
            # Take the first N nuclear FASTA sequences
            main_fasta = fasta_sorted[:len(vcf_sorted)]
            mapping = dict(zip(vcf_sorted, main_fasta))
            
            logger.debug(f"Suggested mapping: {mapping}")
            return mapping
        else:
            # VCF has more chromosomes than nuclear FASTA sequences - problematic
            logger.warning(f"VCF has {len(vcf_sorted)} chromosomes but only {len(fasta_sorted)} nuclear FASTA sequences")
            return {}
    
    def _check_chromosome_compatibility(self, vcf_file: str, fasta_file: str) -> Dict:
        """
        Check if VCF and FASTA files have compatible chromosome names.
        
        Performs comprehensive analysis of chromosome naming compatibility
        between VCF and FASTA files, providing detailed compatibility
        information and suggested mapping strategies.
        
        Args:
            vcf_file: Path to VCF file
            fasta_file: Path to FASTA file
            
        Returns:
            Dictionary containing analysis results with compatibility info and suggested actions
        """
        logger.debug("Checking chromosome compatibility between VCF and FASTA files")
        
        try:
            # Get chromosome information from both files
            vcf_chroms = self._get_vcf_chromosomes(vcf_file)
            fasta_seqs = self._get_fasta_sequences(fasta_file)
            
            # Check for exact matches
            vcf_set = set(vcf_chroms.keys())
            fasta_set = set(fasta_seqs.keys())
            
            exact_matches = vcf_set & fasta_set
            vcf_only = vcf_set - fasta_set
            fasta_only = fasta_set - vcf_set
            
            # Generate analysis report
            analysis = {
                'vcf_chromosomes': vcf_chroms,
                'fasta_sequences': fasta_seqs,
                'exact_matches': exact_matches,
                'vcf_only': vcf_only,
                'fasta_only': fasta_only,
                'compatible': len(exact_matches) > 0,
                'needs_mapping': len(vcf_only) > 0 or len(fasta_only) > 0
            }
            
            # Log the analysis
            logger.debug(f"VCF file: {len(vcf_chroms)} chromosomes")
            logger.debug(f"FASTA file: {len(fasta_seqs)} sequences")
            
            if exact_matches:
                logger.debug(f"{len(exact_matches)} chromosome(s) match exactly")
            
            if vcf_only:
                logger.debug(f"{len(vcf_only)} chromosome(s) only in VCF")
            
            if fasta_only:
                logger.debug(f"{len(fasta_only)} sequence(s) only in FASTA")
            
            # Provide recommendations
            if analysis['compatible'] and not analysis['needs_mapping']:
                logger.debug("Files are fully compatible - no chromosome renaming needed")
                analysis['action'] = 'none'
            elif analysis['needs_mapping']:
                logger.debug("Files need chromosome name mapping")
                analysis['action'] = 'mapping'
                # Generate automatic mapping suggestion
                analysis['suggested_mapping'] = self._suggest_chromosome_mapping(vcf_chroms, fasta_seqs)
            else:
                logger.debug("No compatible chromosomes found")
                analysis['action'] = 'error'
            
            return analysis
            
        except Exception as e:
            error_msg = f"Error in compatibility check: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
    
    def _analyze_vcf_file(self, vcf_file: str) -> List[Dict]:
        """
        Analyze VCF file for formatting and content issues.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            List of issue dictionaries
        """
        logger.debug(f"Analyzing VCF file: {vcf_file}")
        issues = []
        
        # Check if file is properly bgzipped
        if not self._is_properly_bgzipped(vcf_file):
            issues.append({
                'type': 'vcf_compression',
                'file': vcf_file,
                'description': 'VCF file needs proper bgzip compression',
                'action': 'Recompress with bgzip'
            })
        
        # Check if file is indexed
        if not self._is_vcf_indexed(vcf_file):
            issues.append({
                'type': 'vcf_index',
                'file': vcf_file,
                'description': 'VCF file needs tabix indexing',
                'action': 'Create tabix index'
            })
        
        # Check if INFO/AF field is present
        if not self._has_af_field(vcf_file):
            issues.append({
                'type': 'vcf_af_field',
                'file': vcf_file,
                'description': 'VCF file lacks INFO/AF field',
                'action': 'Add INFO/AF field using bcftools +fill-tags'
            })
        
        # Check if VCF is normalized
        if not self._is_vcf_normalized(vcf_file):
            issues.append({
                'type': 'vcf_normalization',
                'file': vcf_file,
                'description': 'VCF file needs normalization',
                'action': 'Normalize variants with bcftools norm'
            })
        
        logger.debug(f"VCF analysis found {len(issues)} issues")
        return issues
    
    def _analyze_fasta_file(self, fasta_file: str) -> List[Dict]:
        """
        Analyze FASTA file for indexing requirements.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            List of issue dictionaries
        """
        logger.debug(f"Analyzing FASTA file: {fasta_file}")
        issues = []
        
        # Check if file is indexed
        if not self._is_fasta_indexed(fasta_file):
            issues.append({
                'type': 'fasta_index',
                'file': fasta_file,
                'description': 'FASTA file needs samtools indexing',
                'action': 'Create samtools faidx index'
            })
        
        logger.debug(f"FASTA analysis found {len(issues)} issues")
        return issues
    
    def _analyze_gff_file(self, gff_file: str) -> List[Dict]:
        """
        Analyze GFF file for formatting issues.
        
        Args:
            gff_file: Path to GFF file
            
        Returns:
            List of issue dictionaries
        """
        logger.debug(f"Analyzing GFF file: {gff_file}")
        issues = []
        
        # Check if file is sorted and indexed (optional improvement)
        if not self._is_gff_indexed(gff_file):
            issues.append({
                'type': 'gff_index',
                'file': gff_file,
                'description': 'GFF file could benefit from tabix indexing',
                'action': 'Sort and index GFF file'
            })
        
        logger.debug(f"GFF analysis found {len(issues)} issues")
        return issues
    
    def _analyze_chromosome_compatibility(self, vcf_file: str, fasta_file: str,
                                        gff_file: Optional[str] = None) -> List[Dict]:
        """
        Analyze chromosome name compatibility between files.
        
        Args:
            vcf_file: Path to VCF file
            fasta_file: Path to FASTA file
            gff_file: Optional path to GFF file
            
        Returns:
            List of issue dictionaries
        """
        logger.debug("Analyzing chromosome compatibility")
        issues = []
        
        try:
            # Check VCF-FASTA compatibility
            compat_analysis = self._check_chromosome_compatibility(vcf_file, fasta_file)
            
            if compat_analysis['needs_mapping']:
                issues.append({
                    'type': 'chromosome_mapping',
                    'files': [vcf_file, fasta_file],
                    'description': 'Chromosome names differ between VCF and FASTA files',
                    'action': 'Rename chromosomes in VCF to match FASTA',
                    'mapping': compat_analysis.get('suggested_mapping', {})
                })
            
            # If GFF provided, check GFF-FASTA compatibility
            if gff_file:
                gff_chroms = self._get_gff_chromosomes(gff_file)
                fasta_chroms = set(compat_analysis['fasta_sequences'].keys())
                
                gff_only = set(gff_chroms) - fasta_chroms
                if gff_only:
                    issues.append({
                        'type': 'gff_chromosome_mapping',
                        'files': [gff_file, fasta_file],
                        'description': f'GFF contains chromosomes not in FASTA: {", ".join(sorted(gff_only))}',
                        'action': 'Filter GFF to match FASTA chromosomes'
                    })
        
        except Exception as e:
            logger.warning(f"Error analyzing chromosome compatibility: {str(e)}")
            # Don't raise error, just note the issue
            issues.append({
                'type': 'compatibility_check_failed',
                'description': f'Could not verify chromosome compatibility: {str(e)}',
                'action': 'Manual verification recommended'
            })
        
        logger.debug(f"Compatibility analysis found {len(issues)} issues")
        return issues
    
    def _is_properly_bgzipped(self, vcf_file: str) -> bool:
        """
        Check if VCF file is properly bgzip compressed.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            True if properly bgzipped, False otherwise
        """
        try:
            # Test with bcftools view
            cmd = ['bcftools', 'view', '-h', vcf_file]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            return result.returncode == 0
        except Exception:
            return False
    
    def _is_vcf_indexed(self, vcf_file: str) -> bool:
        """
        Check if VCF file has tabix index.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            True if indexed, False otherwise
        """
        # Check for .tbi or .csi index files
        index_files = [
            f"{vcf_file}.tbi",
            f"{vcf_file}.csi"
        ]
        
        return any(os.path.exists(idx) for idx in index_files)
    
    def _has_af_field(self, vcf_file: str) -> bool:
        """
        Check if VCF file has INFO/AF field.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            True if AF field is present, False otherwise
        """
        try:
            # Get header and check for AF field
            cmd = ['bcftools', 'view', '-h', vcf_file]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                header = result.stdout
                # Look for INFO AF field definition
                return '##INFO=<ID=AF,' in header
            return False
            
        except Exception:
            logger.debug("Error checking AF field presence")
            return False
    
    def _is_vcf_normalized(self, vcf_file: str) -> bool:
        """
        Check if VCF file appears to be normalized.
        
        This is a heuristic check - we assume if the file works with bcftools
        and doesn't have obvious normalization issues, it's probably normalized.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            True if appears normalized, False otherwise
        """
        try:
            # Sample a few variants to check for obvious normalization issues
            cmd = ['bcftools', 'view', '-H', vcf_file]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')[:100]  # Check first 100 variants
                
                for line in lines:
                    if not line.strip():
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) >= 5:
                        ref = fields[3]
                        alt = fields[4]
                        
                        # Check for obvious unnormalized patterns
                        # 1. Multi-character variants that could be left-aligned
                        if len(ref) > 1 and len(alt) > 1:
                            if ref[0] == alt[0]:  # Common prefix suggests not left-aligned
                                logger.debug(f"Potential normalization issue: REF={ref}, ALT={alt}")
                                return False
                
                return True
            
            return False
            
        except Exception:
            logger.debug("Error checking VCF normalization")
            # If we can't check, assume it needs normalization to be safe
            return False
    
    def _is_fasta_indexed(self, fasta_file: str) -> bool:
        """
        Check if FASTA file has samtools index.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            True if indexed, False otherwise
        """
        index_file = f"{fasta_file}.fai"
        return os.path.exists(index_file)
    
    def _is_gff_indexed(self, gff_file: str) -> bool:
        """
        Check if GFF file has tabix index.
        
        Args:
            gff_file: Path to GFF file
            
        Returns:
            True if indexed, False otherwise
        """
        # Check for .tbi or .csi index files
        index_files = [
            f"{gff_file}.tbi",
            f"{gff_file}.csi",
            f"{gff_file}.gz.tbi",
            f"{gff_file}.gz.csi"
        ]
        
        return any(os.path.exists(idx) for idx in index_files)
    
    def _get_gff_chromosomes(self, gff_file: str) -> Set[str]:
        """
        Extract chromosome names from GFF file.
        
        Args:
            gff_file: Path to GFF file
            
        Returns:
            Set of chromosome names found in GFF
        """
        chromosomes = set()
        
        try:
            opener = gzip.open if gff_file.endswith('.gz') else open
            mode = 'rt' if gff_file.endswith('.gz') else 'r'
            
            with opener(gff_file, mode) as f:
                for line_num, line in enumerate(f, 1):
                    if line_num > 1000:  # Sample first 1000 lines
                        break
                    
                    line = line.strip()
                    if line and not line.startswith('#'):
                        fields = line.split('\t')
                        if len(fields) >= 1:
                            chromosomes.add(fields[0])
            
        except Exception as e:
            logger.debug(f"Error reading GFF chromosomes: {str(e)}")
        
        return chromosomes
    
    def _report_issues(self, issues: List[Dict]):
        """
        Report found issues to user.
        
        Args:
            issues: List of issue dictionaries
        """
        logger.info(f"\nFound {len(issues)} issue(s) that need to be addressed:")
        
        for i, issue in enumerate(issues, 1):
            logger.info(f"\n{i}. {issue['description']}")

    
    def _get_user_consent(self, issues: List[Dict]) -> bool:
        """
        Ask user for consent to create corrected files.
        
        Args:
            issues: List of issues to be fixed
            
        Returns:
            True if user consents, False otherwise
        """
        logger.info("\nTo proceed with the pipeline, corrected copies of your files will be created.")
        
        while True:
            try:
                response = input("\n>>> Create corrected files? <<<\n[y/n]: ").strip().lower()
                if response in ['y', 'yes']:
                    return True
                elif response in ['n', 'no']:
                    return False
                else:
                    print("Please enter 'y' for yes or 'n' for no.")
            except (EOFError, KeyboardInterrupt):
                return False
    
    def _prepare_corrected_files(self, vcf_file: str, fasta_file: str,
                               gff_file: Optional[str], analysis: Dict) -> Dict:
        """
        Create corrected versions of files based on analysis.
        
        Args:
            vcf_file: Original VCF file path
            fasta_file: Original FASTA file path
            gff_file: Optional original GFF file path
            analysis: Analysis results with issues to fix
            
        Returns:
            Dictionary mapping file types to prepared file paths
        """
        prepared_files = {}
        
        # Prepare VCF file
        vcf_issues = [issue for issue in analysis['issues'] 
                     if issue.get('type', '').startswith('vcf') or 
                        issue.get('type') == 'chromosome_mapping']
        
        if vcf_issues:
            prepared_vcf = self._prepare_vcf_file(vcf_file, vcf_issues)
            if prepared_vcf:
                prepared_files['vcf'] = prepared_vcf
        
        # Prepare FASTA file
        fasta_issues = [issue for issue in analysis['issues'] 
                       if issue.get('type', '').startswith('fasta')]
        
        if fasta_issues:
            prepared_fasta = self._prepare_fasta_file(fasta_file, fasta_issues)
            if prepared_fasta:
                prepared_files['fasta'] = prepared_fasta
        
        # Prepare GFF file
        if gff_file:
            gff_issues = [issue for issue in analysis['issues'] 
                         if issue.get('type', '').startswith('gff')]
            
            if gff_issues:
                prepared_gff = self._prepare_gff_file(gff_file, gff_issues)
                if prepared_gff:
                    prepared_files['gff'] = prepared_gff
        
        return prepared_files
    
    def _prepare_vcf_file(self, vcf_file: str, issues: List[Dict]) -> Optional[str]:
        """
        Prepare corrected VCF file addressing all identified issues.
        
        Args:
            vcf_file: Original VCF file path
            issues: List of VCF-related issues to fix
            
        Returns:
            Path to prepared VCF file, or None if preparation failed
        """
        logger.debug("=== VCF PREPARATION DEBUG ===")
        logger.debug("Preparing VCF file...")
        
        try:
            # Create output path
            base_name = Path(vcf_file).stem
            if base_name.endswith('.vcf'):
                base_name = base_name[:-4]
            
            output_path = os.path.join(
                os.path.dirname(vcf_file),
                f"{base_name}_prepared.vcf.gz"
            )
            
            # Process through pipeline
            current_file = vcf_file
            temp_files = []
            
            # Step 1: Ensure proper compression
            compression_issues = [issue for issue in issues if issue['type'] == 'vcf_compression']
            if compression_issues:
                logger.debug("Applying bgzip compression")
                current_file = self._bgzip_vcf(current_file)
                temp_files.append(current_file)
            
            # Step 2: Add AF field if missing
            af_issues = [issue for issue in issues if issue['type'] == 'vcf_af_field']
            if af_issues:
                logger.debug("Adding INFO/AF field")
                current_file = self._add_af_field(current_file)
                temp_files.append(current_file)
            
            # Step 3: Apply chromosome mapping BEFORE normalization
            mapping_issues = [issue for issue in issues if issue['type'] == 'chromosome_mapping']
            if mapping_issues:
                logger.debug("Applying chromosome mapping")
                mapping = mapping_issues[0].get('mapping', {})
                if mapping:
                    current_file = self._apply_chromosome_mapping(current_file, mapping)
                    temp_files.append(current_file)
            
            # Step 4: Normalize VCF (after chromosome mapping)
            norm_issues = [issue for issue in issues if issue['type'] == 'vcf_normalization']
            if norm_issues:
                logger.debug("Normalizing VCF")
                current_file = self._normalize_vcf(current_file)
                temp_files.append(current_file)
            
            # Step 5: Copy to final location and index
            shutil.copy2(current_file, output_path)
            self._index_vcf(output_path)
            
            # Clean up temporary files
            for temp_file in temp_files:
                if temp_file != vcf_file and os.path.exists(temp_file):
                    try:
                        os.unlink(temp_file)
                    except OSError:
                        pass
            
            logger.info(f"Prepared VCF file: {output_path}")
            logger.debug("=== END VCF PREPARATION DEBUG ===")
            return output_path
            
        except Exception as e:
            error_msg = f"Error preparing VCF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END VCF PREPARATION DEBUG ===")
            return None
    
    def _prepare_fasta_file(self, fasta_file: str, issues: List[Dict]) -> Optional[str]:
        """
        Prepare FASTA file (mainly indexing).
        
        Args:
            fasta_file: Original FASTA file path
            issues: List of FASTA-related issues to fix
            
        Returns:
            Path to prepared FASTA file, or None if preparation failed
        """
        logger.debug("=== FASTA PREPARATION DEBUG ===")
        logger.debug("Preparing FASTA file...")
        
        try:
            # For FASTA, we usually just need to create index
            # Copy file is not necessary unless we need to modify it
            
            index_issues = [issue for issue in issues if issue['type'] == 'fasta_index']
            if index_issues:
                logger.debug("Creating FASTA index")
                self._index_fasta(fasta_file)
            
            logger.info("FASTA file preparation completed")
            logger.debug("=== END FASTA PREPARATION DEBUG ===")
            return fasta_file  # Return original file as no copy was needed
            
        except Exception as e:
            error_msg = f"Error preparing FASTA file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END FASTA PREPARATION DEBUG ===")
            return None
    
    def _prepare_gff_file(self, gff_file: str, issues: List[Dict]) -> Optional[str]:
        """
        Prepare GFF file (sorting and indexing).
        
        Args:
            gff_file: Original GFF file path
            issues: List of GFF-related issues to fix
            
        Returns:
            Path to prepared GFF file, or None if preparation failed
        """
        logger.debug("=== GFF PREPARATION DEBUG ===")
        logger.debug("Preparing GFF file...")
        
        try:
            # For GFF, we might need to sort and index
            index_issues = [issue for issue in issues if issue['type'] == 'gff_index']
            
            if index_issues:
                # Create sorted and compressed version
                base_name = Path(gff_file).stem
                if base_name.endswith('.gff'):
                    base_name = base_name[:-4]
                
                output_path = os.path.join(
                    os.path.dirname(gff_file),
                    f"{base_name}_prepared.gff.gz"
                )
                
                logger.debug("Sorting and compressing GFF file")
                self._sort_and_compress_gff(gff_file, output_path)
                self._index_gff(output_path)
                
                logger.info(f"Prepared GFF file: {output_path}")
                logger.debug("=== END GFF PREPARATION DEBUG ===")
                return output_path
            
            logger.info("GFF file preparation completed")
            logger.debug("=== END GFF PREPARATION DEBUG ===")
            return gff_file  # Return original file as no changes needed
            
        except Exception as e:
            error_msg = f"Error preparing GFF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END GFF PREPARATION DEBUG ===")
            return None
    
    def _bgzip_vcf(self, vcf_file: str) -> str:
        """
        Compress VCF file with bgzip.
        
        Args:
            vcf_file: Path to input VCF file
            
        Returns:
            Path to bgzipped VCF file
            
        Raises:
            ExternalToolError: If bgzip compression fails
        """
        try:
            # Create temporary output file
            temp_fd, temp_path = tempfile.mkstemp(
                suffix=".vcf.gz",
                prefix="ddprimer_bgzip_",
                dir=self.temp_dir
            )
            os.close(temp_fd)
            
            # Check if input is already compressed
            if vcf_file.endswith('.gz'):
                # Decompress and recompress with bgzip
                decompress_cmd = ['gunzip', '-c', vcf_file]
                compress_cmd = ['bgzip', '-c']
                
                decompress_proc = subprocess.Popen(decompress_cmd, stdout=subprocess.PIPE,
                                                 stderr=subprocess.PIPE)
                
                with open(temp_path, 'wb') as temp_file:
                    compress_proc = subprocess.Popen(compress_cmd, stdin=decompress_proc.stdout,
                                                   stdout=temp_file, stderr=subprocess.PIPE)
                    decompress_proc.stdout.close()
                    
                    compress_returncode = compress_proc.wait()
                    decompress_returncode = decompress_proc.wait()
                
                if decompress_returncode != 0 or compress_returncode != 0:
                    os.unlink(temp_path)
                    error_msg = "Failed to recompress VCF with bgzip"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="bgzip")
            else:
                # Compress uncompressed VCF
                compress_cmd = ['bgzip', '-c', vcf_file]
                
                with open(temp_path, 'wb') as temp_file:
                    result = subprocess.run(compress_cmd, stdout=temp_file,
                                          stderr=subprocess.PIPE)
                
                if result.returncode != 0:
                    os.unlink(temp_path)
                    error_msg = f"Failed to compress VCF with bgzip: {result.stderr.decode()}"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="bgzip")
            
            logger.debug(f"Successfully bgzipped VCF: {temp_path}")
            return temp_path
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error bgzipping VCF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="bgzip") from e
    
    def _add_af_field(self, vcf_file: str) -> str:
        """
        Add INFO/AF field to VCF file using bcftools +fill-tags.
        
        Args:
            vcf_file: Path to input VCF file
            
        Returns:
            Path to VCF file with AF field added
            
        Raises:
            ExternalToolError: If bcftools execution fails
        """
        try:
            # Create temporary output file
            temp_fd, temp_path = tempfile.mkstemp(
                suffix=".vcf.gz",
                prefix="ddprimer_af_",
                dir=self.temp_dir
            )
            os.close(temp_fd)
            
            # Use bcftools +fill-tags to add AF field
            cmd = [
                'bcftools', '+fill-tags', vcf_file,
                '-Oz', '-o', temp_path,
                '--', '-t', 'AF'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                os.unlink(temp_path)
                error_msg = f"bcftools +fill-tags failed: {result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="bcftools")
            
            logger.debug(f"Successfully added AF field: {temp_path}")
            return temp_path
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error adding AF field to VCF: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="bcftools") from e
    
    def _normalize_vcf(self, vcf_file: str) -> str:
        """
        Normalize VCF file using bcftools norm.
        
        Args:
            vcf_file: Path to input VCF file
            
        Returns:
            Path to normalized VCF file
            
        Raises:
            ExternalToolError: If bcftools normalization fails
        """
        try:
            # Create temporary output file
            temp_fd, temp_path = tempfile.mkstemp(
                suffix=".vcf.gz",
                prefix="ddprimer_norm_",
                dir=self.temp_dir
            )
            os.close(temp_fd)
            
            # Normalize with bcftools norm
            cmd = [
                'bcftools', 'norm', vcf_file,
                '-f', self._get_reference_file(),
                '-Oz', '-o', temp_path
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                os.unlink(temp_path)
                error_msg = f"bcftools norm failed: {result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="bcftools")
            
            logger.debug(f"Successfully normalized VCF: {temp_path}")
            return temp_path
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error normalizing VCF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="bcftools") from e
    
    def _apply_chromosome_mapping(self, vcf_file: str, mapping: Dict[str, str]) -> str:
        """
        Apply chromosome name mapping to VCF file.
        
        Args:
            vcf_file: Path to input VCF file
            mapping: Dictionary mapping VCF chromosome names to FASTA names
            
        Returns:
            Path to VCF file with renamed chromosomes
            
        Raises:
            ExternalToolError: If chromosome renaming fails
        """
        try:
            # Create temporary output file
            temp_fd, temp_path = tempfile.mkstemp(
                suffix=".vcf.gz",
                prefix="ddprimer_remap_",
                dir=self.temp_dir
            )
            os.close(temp_fd)
            
            # Create chromosome renaming file
            rename_fd, rename_file = tempfile.mkstemp(
                suffix=".txt",
                prefix="ddprimer_chrmap_",
                dir=self.temp_dir
            )
            
            try:
                with os.fdopen(rename_fd, 'w') as f:
                    for old_name, new_name in mapping.items():
                        f.write(f"{old_name}\t{new_name}\n")
                
                # Apply chromosome renaming
                cmd = [
                    'bcftools', 'annotate', vcf_file,
                    '--rename-chrs', rename_file,
                    '-Oz', '-o', temp_path
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    os.unlink(temp_path)
                    error_msg = f"bcftools annotate (rename-chrs) failed: {result.stderr}"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="bcftools")
                
                logger.debug(f"Successfully remapped chromosomes: {temp_path}")
                return temp_path
                
            finally:
                # Clean up rename file
                if os.path.exists(rename_file):
                    os.unlink(rename_file)
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error applying chromosome mapping: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="bcftools") from e
    
    def _index_vcf(self, vcf_file: str):
        """
        Create tabix index for VCF file.
        
        Args:
            vcf_file: Path to VCF file to index
            
        Raises:
            ExternalToolError: If indexing fails
        """
        try:
            cmd = ['tabix', '-p', 'vcf', vcf_file]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                error_msg = f"tabix indexing failed: {result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="tabix")
            
            logger.debug(f"Successfully indexed VCF: {vcf_file}")
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error indexing VCF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="tabix") from e
    
    def _index_fasta(self, fasta_file: str):
        """
        Create samtools index for FASTA file.
        
        Args:
            fasta_file: Path to FASTA file to index
            
        Raises:
            ExternalToolError: If indexing fails
        """
        try:
            cmd = ['samtools', 'faidx', fasta_file]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                error_msg = f"samtools faidx failed: {result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="samtools")
            
            logger.debug(f"Successfully indexed FASTA: {fasta_file}")
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error indexing FASTA file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="samtools") from e
    
    def _sort_and_compress_gff(self, gff_file: str, output_path: str):
        """
        Sort and compress GFF file for tabix indexing.
        
        Args:
            gff_file: Path to input GFF file
            output_path: Path for output compressed GFF
            
        Raises:
            ExternalToolError: If sorting/compression fails
        """
        try:
            # Sort GFF file by chromosome and position
            sort_cmd = ['sort', '-k1,1', '-k4,4n', gff_file]
            compress_cmd = ['bgzip', '-c']
            
            # Chain sort and bgzip
            sort_proc = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            
            with open(output_path, 'wb') as output_file:
                compress_proc = subprocess.Popen(compress_cmd, stdin=sort_proc.stdout,
                                               stdout=output_file, stderr=subprocess.PIPE)
                sort_proc.stdout.close()
                
                compress_returncode = compress_proc.wait()
                sort_returncode = sort_proc.wait()
            
            if sort_returncode != 0 or compress_returncode != 0:
                if os.path.exists(output_path):
                    os.unlink(output_path)
                error_msg = "Failed to sort and compress GFF file"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="sort/bgzip")
            
            logger.debug(f"Successfully sorted and compressed GFF: {output_path}")
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error sorting and compressing GFF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="sort/bgzip") from e
    
    def _index_gff(self, gff_file: str):
        """
        Create tabix index for GFF file.
        
        Args:
            gff_file: Path to compressed GFF file to index
            
        Raises:
            ExternalToolError: If indexing fails
        """
        try:
            cmd = ['tabix', '-p', 'gff', gff_file]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                error_msg = f"tabix indexing of GFF failed: {result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="tabix")
            
            logger.debug(f"Successfully indexed GFF: {gff_file}")
            
        except ExternalToolError:
            raise
        except Exception as e:
            error_msg = f"Error indexing GFF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="tabix") from e
    
    def _get_reference_file(self) -> str:
        """
        Get reference FASTA file for normalization.
        
        This method needs to be called in context where reference file is available.
        For now, it raises an error if called outside proper context.
        
        Returns:
            Path to reference FASTA file
            
        Raises:
            FileError: If reference file is not available in current context
        """
        # This will be set by the calling code when needed
        if hasattr(self, '_reference_file') and self._reference_file:
            return self._reference_file
        
        error_msg = "Reference FASTA file not available for normalization"
        logger.error(error_msg)
        raise FileError(error_msg)
    
    def set_reference_file(self, reference_file: str):
        """
        Set reference FASTA file for operations that require it.
        
        Args:
            reference_file: Path to reference FASTA file
            
        Raises:
            FileError: If reference file does not exist
        """
        if not os.path.exists(reference_file):
            error_msg = f"Reference FASTA file not found: {reference_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        self._reference_file = reference_file
        logger.debug(f"Set reference file: {reference_file}")
    
    def cleanup(self):
        """
        Clean up any temporary files created during preparation.
        
        Should be called when file preparation is complete to ensure
        no temporary files are left behind.
        """
        logger.debug("Cleaning up FilePreparator resources")
        
        # Clean up any tracked prepared files if needed
        if hasattr(self, 'prepared_files'):
            self.prepared_files.clear()
        
        # Note: Individual temporary files are cleaned up in their respective methods
        logger.debug("FilePreparator cleanup completed")


# Convenience function for integration with main pipeline
def prepare_pipeline_files(vcf_file: str, fasta_file: str, 
                         gff_file: Optional[str] = None) -> Dict:
    """
    Convenience function to prepare files for ddPrimer pipeline.
    
    This function provides a simple interface for file preparation that can
    be easily integrated into the main pipeline workflow.
    
    Args:
        vcf_file: Path to VCF file
        fasta_file: Path to FASTA file
        gff_file: Optional path to GFF file
        
    Returns:
        Dictionary with preparation results and final file paths
        
    Raises:
        FileError: If file preparation fails
        ExternalToolError: If external tools fail
        
    Example:
        >>> result = prepare_pipeline_files("variants.vcf", "genome.fasta", "genes.gff")
        >>> if result['success']:
        ...     vcf_path = result['vcf_file']
        ...     fasta_path = result['fasta_file']
        ...     gff_path = result['gff_file']
    """
    preparator = FilePreparator()
    
    try:
        # Set reference file for normalization operations
        preparator.set_reference_file(fasta_file)
        
        # Prepare files
        result = preparator.prepare_files(vcf_file, fasta_file, gff_file)
        
        return result
        
    finally:
        # Clean up
        preparator.cleanup()