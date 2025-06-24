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
from typing import Dict, List, Optional, Set
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
            self.temp_dir = os.path.join(Config.get_user_config_dir(), "temp")
            os.makedirs(self.temp_dir, exist_ok=True)
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
        logger.info("\nAnalyzing input files for compatibility...")
        
        try:
            # Validate input files exist
            self._validate_input_files(vcf_file, fasta_file, gff_file)
            
            # Set reference file for normalization
            self.set_reference_file(fasta_file)
            
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
            
            # CRITICAL FIX: Check if VCF preparation actually succeeded
            vcf_success = True
            if 'vcf' in prepared_files:
                vcf_path = prepared_files['vcf']
                if not vcf_path or not os.path.exists(vcf_path):
                    vcf_success = False
                    error_msg = "VCF file preparation failed"
            else:
                # Check if VCF preparation was attempted but failed
                vcf_issues = [issue for issue in analysis['issues'] 
                            if issue.get('type', '').startswith('vcf') or 
                                issue.get('type') == 'chromosome_mapping']
                if vcf_issues:
                    vcf_success = False
                    error_msg = "VCF file preparation was required but failed"
            
            if not vcf_success:
                logger.error(error_msg)
                return {
                    'success': False,
                    'reason': error_msg,
                    'issues_found': analysis['issues']
                }
            
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
            return {
                'success': False,
                'reason': error_msg,
                'issues_found': analysis.get('issues', [])
            }
    
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
            
            # Flag as needing mapping if there are actually incompatible chromosomes
            vcf_chroms = set(compat_analysis['vcf_chromosomes'].keys())
            fasta_chroms = set(compat_analysis['fasta_sequences'].keys())
            
            # Check if all VCF chromosomes are present in FASTA
            missing_in_fasta = vcf_chroms - fasta_chroms
            
            if missing_in_fasta:
                # Only create mapping issue if there are actually missing chromosomes
                issues.append({
                    'type': 'chromosome_mapping',
                    'files': [vcf_file, fasta_file],
                    'description': f'VCF chromosomes not found in FASTA: {", ".join(sorted(missing_in_fasta))}',
                    'action': 'Rename chromosomes in VCF to match FASTA',
                    'mapping': compat_analysis.get('suggested_mapping', {})
                })
            else:
                logger.debug("All VCF chromosomes found in FASTA - no mapping needed")
            
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
        
        Uses bcftools norm --check-ref to properly validate normalization.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            True if appears normalized, False otherwise
        """
        try:
            # For large VCF files, just check the first 10000 variants for normalization issues
            # Use proper subprocess chaining instead of shell pipes to handle spaces in paths
            
            # Step 1: Get first 10000 lines of VCF
            view_cmd = ['bcftools', 'view', vcf_file]
            head_cmd = ['head', '-n', '10000']
            
            # Step 2: Check normalization on the sample
            norm_cmd = ['bcftools', 'norm', '--check-ref', 'w', '-f', self._get_reference_file(), '-']
            
            logger.debug(f"Running normalization check on sample of: {vcf_file}")
            
            # Chain the commands properly using subprocess
            view_proc = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            head_proc = subprocess.Popen(head_cmd, stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            view_proc.stdout.close()  # Allow view_proc to receive a SIGPIPE if head_proc exits
            
            norm_proc = subprocess.Popen(norm_cmd, stdin=head_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            head_proc.stdout.close()  # Allow head_proc to receive a SIGPIPE if norm_proc exits
            
            # Wait for the final process and get output
            try:
                stdout, stderr = norm_proc.communicate(timeout=30)
                
                # Close any remaining processes
                view_proc.wait()
                head_proc.wait()
                
                if norm_proc.returncode == 0:
                    # Check stderr for normalization warnings
                    stderr_output = stderr.lower()
                    normalization_warnings = [
                        'not left-aligned',
                        'not normalized', 
                        'multiallelic'
                    ]
                    
                    # If no normalization warnings in the sample, assume it's normalized
                    has_warnings = any(warning in stderr_output for warning in normalization_warnings)
                    
                    if has_warnings:
                        logger.debug(f"VCF normalization warnings found in sample: {stderr}")
                        return False
                    else:
                        logger.debug("VCF sample appears to be normalized")
                        return True
                else:
                    logger.debug(f"bcftools norm check failed: {stderr}")
                    # If reference sequence errors, skip normalization 
                    if "sequence not found" in stderr.lower():
                        logger.debug("Reference sequence mismatch - assuming VCF is already normalized")
                        return True
                    return False
                    
            except subprocess.TimeoutExpired:
                logger.debug("VCF normalization check timed out - assuming normalized")
                # Kill any running processes
                try:
                    norm_proc.kill()
                    head_proc.kill() 
                    view_proc.kill()
                except:
                    pass
                return True
                
        except Exception as e:
            logger.debug(f"Error checking VCF normalization: {str(e)}")
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
        logger.info("\nTo proceed with the pipeline, corrected copies of your files need to be created.")
        
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
        Prepare corrected VCF file with immediate cleanup to minimize disk usage.
        
        Args:
            vcf_file: Original VCF file path
            issues: List of VCF-related issues to fix
            
        Returns:
            Path to prepared VCF file, or None if preparation failed
        """
        logger.debug("=== VCF PREPARATION DEBUG ===")
        logger.debug("Preparing VCF file with optimized disk usage...")
        
        try:
            # Create output path
            base_name = Path(vcf_file).stem
            if base_name.endswith('.vcf'):
                base_name = base_name[:-4]
            
            output_dir = os.path.dirname(os.path.abspath(vcf_file))
            output_path = os.path.join(output_dir, f"{base_name}_prepared.vcf.gz")
            
            logger.debug(f"Input VCF: {vcf_file}")
            logger.debug(f"Output VCF: {output_path}")
            
            # Clean up any existing broken prepared file
            if os.path.exists(output_path):
                logger.debug("Removing existing prepared file")
                os.unlink(output_path)
                for ext in ['.tbi', '.csi']:
                    idx_file = f"{output_path}{ext}"
                    if os.path.exists(idx_file):
                        os.unlink(idx_file)
            
            # Start with the original file
            current_file = vcf_file
            previous_temp_file = None  # Track only the immediately previous temp file
            
            # Helper function for immediate cleanup
            def cleanup_previous_temp():
                nonlocal previous_temp_file
                if previous_temp_file and previous_temp_file != vcf_file and os.path.exists(previous_temp_file):
                    try:
                        os.unlink(previous_temp_file)
                        # Clean up any associated index files
                        for ext in ['.tbi', '.csi']:
                            idx_file = f"{previous_temp_file}{ext}"
                            if os.path.exists(idx_file):
                                os.unlink(idx_file)
                        logger.debug(f"Cleaned up intermediate file: {previous_temp_file}")
                    except OSError as e:
                        logger.debug(f"Could not clean up {previous_temp_file}: {e}")
                    previous_temp_file = None
            
            # Check available disk space before starting
            available_space = self._get_available_disk_space(self.temp_dir)
            file_size = os.path.getsize(vcf_file)
            
            # Estimate space needed (current file + 2 temp files worth of space for safety)
            estimated_space_needed = file_size * 3
            if available_space < estimated_space_needed:
                logger.warning(f"Low disk space detected. Available: {available_space / (1024**3):.1f}GB, "
                            f"Estimated needed: {estimated_space_needed / (1024**3):.1f}GB")
                # Continue anyway but warn user
            
            # Step 1: Ensure proper compression (ALWAYS do this for .vcf files)
            needs_compression = (
                not vcf_file.endswith('.gz') or 
                any(issue['type'] == 'vcf_compression' for issue in issues)
            )
            
            if needs_compression:
                logger.debug("Step 1: Applying bgzip compression")
                new_file = self._bgzip_vcf(current_file)
                previous_temp_file = current_file if current_file != vcf_file else None
                current_file = new_file
                cleanup_previous_temp()  # Clean up immediately
                logger.debug(f"Compression complete: {current_file}")
            
            # Step 2: Create initial index
            index_issues = [issue for issue in issues if issue['type'] == 'vcf_index']
            if index_issues or needs_compression:
                logger.debug("Step 2: Creating initial index")
                try:
                    self._index_vcf(current_file)
                    logger.debug("Initial indexing complete")
                except ExternalToolError as e:
                    logger.warning(f"Initial indexing failed: {e}. Will retry after processing.")
            
            # Step 3: Add AF field if missing
            af_issues = [issue for issue in issues if issue['type'] == 'vcf_af_field']
            if af_issues:
                logger.debug("Step 3: Adding INFO/AF field")
                new_file = self._add_af_field(current_file)
                previous_temp_file = current_file if current_file != vcf_file else None
                current_file = new_file
                cleanup_previous_temp()  # Clean up immediately
                logger.debug(f"AF field addition complete: {current_file}")
                
                # Re-index after AF field addition
                try:
                    self._index_vcf(current_file)
                except ExternalToolError:
                    logger.debug("Re-indexing after AF addition failed, will continue")
            
            # Step 4: Apply chromosome mapping if needed
            mapping_issues = [issue for issue in issues if issue['type'] == 'chromosome_mapping']
            if mapping_issues:
                logger.debug("Step 4: Applying chromosome mapping")
                mapping = mapping_issues[0].get('mapping', {})
                if mapping:
                    new_file = self._apply_chromosome_mapping(current_file, mapping)
                    previous_temp_file = current_file if current_file != vcf_file else None
                    current_file = new_file
                    cleanup_previous_temp()  # Clean up immediately
                    logger.debug(f"Chromosome mapping complete: {current_file}")
                    
                    # Re-index after chromosome mapping
                    try:
                        self._index_vcf(current_file)
                    except ExternalToolError:
                        logger.debug("Re-indexing after mapping failed, will continue")
            
            # Step 5: Normalize VCF (after chromosome mapping)
            norm_issues = [issue for issue in issues if issue['type'] == 'vcf_normalization']
            if norm_issues:
                logger.debug("Step 5: Normalizing VCF")
                try:
                    new_file = self._normalize_vcf(current_file)
                    previous_temp_file = current_file if current_file != vcf_file else None
                    current_file = new_file
                    cleanup_previous_temp()  # Clean up immediately
                    logger.debug(f"Normalization complete: {current_file}")
                except ExternalToolError as e:
                    if "chromosome not found" in str(e).lower():
                        logger.warning("VCF normalization failed due to chromosome mismatch, skipping")
                    else:
                        cleanup_previous_temp()  # Clean up on error too
                        raise
            
            # Step 6: Copy to final location
            logger.debug("Step 6: Copying to final location")
            shutil.copy2(current_file, output_path)

            # Remove temporary file after successful copy
            if current_file != vcf_file:
                try:
                    os.unlink(current_file)
                    # Also clean up any associated index files
                    for ext in ['.tbi', '.csi']:
                        idx_file = f"{current_file}{ext}"
                        if os.path.exists(idx_file):
                            os.unlink(idx_file)
                    logger.debug(f"Cleaned up final temporary file: {current_file}")
                except OSError as e:
                    logger.debug(f"Could not clean up final temp file {current_file}: {e}")
            
            # Step 7: Create final index and verify
            logger.debug("Step 7: Creating final index")
            try:
                self._index_vcf(output_path)
            except ExternalToolError as e:
                error_msg = f"VCF indexing failed, file preparation incomplete: {str(e)}"
                logger.error(error_msg)
                if os.path.exists(output_path):
                    os.unlink(output_path)
                cleanup_previous_temp()
                raise ExternalToolError(error_msg, tool_name="tabix")
            
            # Step 8: Final verification
            logger.debug("Step 8: Final verification")
            if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                error_msg = f"Prepared VCF file is missing or empty: {output_path}"
                logger.error(error_msg)
                cleanup_previous_temp()
                return None
            
            # Test that bcftools can read the final file
            test_cmd = ['bcftools', 'view', '-h', output_path]
            test_result = subprocess.run(test_cmd, capture_output=True, text=True, timeout=30)
            if test_result.returncode != 0:
                error_msg = f"Final VCF file cannot be read by bcftools: {test_result.stderr}"
                logger.error(error_msg)
                if os.path.exists(output_path):
                    os.unlink(output_path)
                cleanup_previous_temp()
                return None
            
            # Final cleanup of the last temporary file
            cleanup_previous_temp()
            
            logger.info(f"Prepared VCF file: {output_path}")
            logger.debug("=== END VCF PREPARATION DEBUG ===")
            return output_path
            
        except Exception as e:
            error_msg = f"Error preparing VCF file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            # Clean up on error
            if 'output_path' in locals() and os.path.exists(output_path):
                try:
                    os.unlink(output_path)
                except OSError:
                    pass
            
            # Clean up current temp file
            if 'previous_temp_file' in locals():
                try:
                    cleanup_previous_temp()
                except:
                    pass
            
            logger.debug("=== END VCF PREPARATION DEBUG ===")
            return None

    def _get_available_disk_space(self, path: str) -> int:
        """
        Get available disk space in bytes for the given path.
        
        Args:
            path: Directory path to check
            
        Returns:
            Available space in bytes
        """
        try:
            statvfs = os.statvfs(path)
            # Available space = fragment size * available fragments
            return statvfs.f_frsize * statvfs.f_bavail
        except Exception as e:
            logger.debug(f"Could not determine disk space for {path}: {e}")
            return float('inf')  # Assume unlimited if we can't check

    def _check_disk_space_before_operation(self, input_file: str, multiplier: float = 2.0) -> bool:
        """
        Check if there's enough disk space for an operation.
        
        Args:
            input_file: Path to input file
            multiplier: Safety multiplier for space estimation
            
        Returns:
            True if enough space, False otherwise
        """
        try:
            file_size = os.path.getsize(input_file)
            available_space = self._get_available_disk_space(self.temp_dir)
            needed_space = file_size * multiplier
            
            if available_space < needed_space:
                logger.warning(f"Insufficient disk space. Available: {available_space / (1024**3):.1f}GB, "
                            f"Needed: {needed_space / (1024**3):.1f}GB")
                return False
            return True
        except Exception:
            return True  # If we can't check, assume it's fine


    
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
                elif base_name.endswith('.gff3'):
                    base_name = base_name[:-5]
                
                # Use same directory as original file for output
                output_dir = os.path.dirname(os.path.abspath(gff_file))
                output_path = os.path.join(output_dir, f"{base_name}_prepared.gff.gz")
                
                logger.debug(f"Original GFF: {gff_file}")
                logger.debug(f"Prepared GFF will be: {output_path}")
                
                # Check if the prepared file already exists
                if os.path.exists(output_path):
                    logger.debug(f"Prepared GFF file already exists: {output_path}")
                    user_input = input(f"\nPrepared GFF file already exists: {output_path}\nOverwrite? [y/n]: ").strip().lower()
                    if user_input not in ['y', 'yes']:
                        logger.info("Using existing prepared GFF file")
                        self._index_gff(output_path)  # Ensure it's indexed
                        return output_path
                    else:
                        os.unlink(output_path)  # Remove existing file
                
                logger.debug("Sorting and compressing GFF file")
                self._sort_and_compress_gff(gff_file, output_path)
                
                # Verify the file was actually created
                if not os.path.exists(output_path):
                    error_msg = f"GFF preparation failed - output file not created: {output_path}"
                    logger.error(error_msg)
                    return None
                
                logger.debug("Creating tabix index")
                self._index_gff(output_path)
                
                # Verify index was created
                index_files = [f"{output_path}.tbi", f"{output_path}.csi"]
                if not any(os.path.exists(idx) for idx in index_files):
                    logger.warning("GFF index file was not created, but continuing")
                
                logger.info(f"Prepared GFF file: {output_path}")
                logger.debug("=== END GFF PREPARATION DEBUG ===")
                return output_path
            
            logger.info("GFF file preparation completed (no changes needed)")
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
            
            logger.debug(f"Input file: {vcf_file}")
            logger.debug(f"Output file: {temp_path}")
            logger.debug(f"Input file exists: {os.path.exists(vcf_file)}")
            logger.debug(f"Input file size: {os.path.getsize(vcf_file) if os.path.exists(vcf_file) else 'N/A'}")
            
            # Check if input is already compressed
            if vcf_file.endswith('.gz'):
                logger.debug("Input is already compressed, checking if it's bgzip...")
                
                # Test if it's already properly bgzipped
                test_cmd = ['bcftools', 'view', '-h', vcf_file]
                test_result = subprocess.run(test_cmd, capture_output=True, text=True, timeout=30)
                
                if test_result.returncode == 0:
                    logger.debug("Input is already properly bgzipped, just copying")
                    shutil.copy2(vcf_file, temp_path)
                    return temp_path
                else:
                    logger.debug("Input is compressed but not bgzip, will recompress")
                    # Decompress and recompress with bgzip
                    with gzip.open(vcf_file, 'rt') as input_file:
                        compress_cmd = ['bgzip', '-c']
                        
                        with open(temp_path, 'wb') as output_file:
                            compress_proc = subprocess.run(
                                compress_cmd, 
                                input=input_file.read(), 
                                stdout=output_file,
                                stderr=subprocess.PIPE,
                                text=True
                            )
                        
                        if compress_proc.returncode != 0:
                            os.unlink(temp_path)
                            error_msg = f"Failed to recompress with bgzip: {compress_proc.stderr}"
                            logger.error(error_msg)
                            raise ExternalToolError(error_msg, tool_name="bgzip")
            else:
                logger.debug("Input is uncompressed, compressing with bgzip...")
                # Compress uncompressed VCF
                compress_cmd = ['bgzip', '-c', vcf_file]
                
                with open(temp_path, 'wb') as output_file:
                    result = subprocess.run(
                        compress_cmd, 
                        stdout=output_file,
                        stderr=subprocess.PIPE,
                        timeout=300
                    )
                
                if result.returncode != 0:
                    os.unlink(temp_path)
                    error_msg = f"Failed to compress VCF with bgzip: {result.stderr.decode()}"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="bgzip")
            
            # Verify the output is properly bgzipped
            logger.debug("Verifying bgzip compression...")
            verify_cmd = ['bcftools', 'view', '-h', temp_path]
            verify_result = subprocess.run(verify_cmd, capture_output=True, text=True, timeout=30)
            
            if verify_result.returncode != 0:
                os.unlink(temp_path)
                error_msg = f"bgzip verification failed: {verify_result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="bgzip")
            
            # Verify file is not empty
            if not os.path.exists(temp_path) or os.path.getsize(temp_path) == 0:
                error_msg = "bgzipped VCF file is empty or was not created"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="bgzip")
            
            logger.debug(f"Successfully bgzipped VCF: {temp_path}")
            logger.debug(f"Output file size: {os.path.getsize(temp_path)} bytes")
            return temp_path
            
        except subprocess.TimeoutExpired:
            if 'temp_path' in locals() and os.path.exists(temp_path):
                os.unlink(temp_path)
            error_msg = "bgzip compression timed out"
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="bgzip")
        except ExternalToolError:
            raise
        except Exception as e:
            if 'temp_path' in locals() and os.path.exists(temp_path):
                os.unlink(temp_path)
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
        
        IMPROVED: Better error handling and reference file validation.
        
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
            
            # Get reference file
            reference_file = self._get_reference_file()
            logger.debug(f"Using reference file for normalization: {reference_file}")
            
            # Verify reference file is indexed
            if not os.path.exists(f"{reference_file}.fai"):
                logger.debug("Reference file not indexed, creating index")
                self._index_fasta(reference_file)
            
            # Normalize with bcftools norm
            cmd = [
                'bcftools', 'norm', vcf_file,
                '-f', reference_file,
                '-m', '-any',  # Split multiallelic variants
                '-N',          # Do not normalize indels
                '-Oz', '-o', temp_path
            ]
            
            logger.debug(f"Normalization command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10800)
            
            if result.returncode != 0:
                os.unlink(temp_path)
                error_msg = f"bcftools norm failed: {result.stderr}"
                logger.error(error_msg)
                # Include more context about the error
                if "sequence not found" in result.stderr.lower():
                    error_msg += "\nThis may indicate chromosome name mismatch between VCF and reference FASTA"
                raise ExternalToolError(error_msg, tool_name="bcftools")
            
            # Log normalization stats if available
            if result.stderr:
                logger.debug(f"Normalization output: {result.stderr}")
            
            logger.debug(f"Successfully normalized VCF: {temp_path}")
            return temp_path
            
        except subprocess.TimeoutExpired:
            if 'temp_path' in locals() and os.path.exists(temp_path):
                os.unlink(temp_path)
            error_msg = "VCF normalization timed out"
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="bcftools")
        except ExternalToolError:
            raise
        except Exception as e:
            if 'temp_path' in locals() and os.path.exists(temp_path):
                os.unlink(temp_path)
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
            logger.debug(f"Input GFF file: {gff_file}")
            logger.debug(f"Output path: {output_path}")
            logger.debug(f"Input file exists: {os.path.exists(gff_file)}")
            logger.debug(f"Output directory exists: {os.path.exists(os.path.dirname(output_path))}")
            
            # Ensure output directory exists
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # Create a temporary intermediate file for sorted GFF
            temp_sorted_fd, temp_sorted_path = tempfile.mkstemp(
                suffix=".gff",
                prefix="ddprimer_sorted_",
                dir=self.temp_dir
            )
            os.close(temp_sorted_fd)
            
            try:
                # First, sort the GFF file with explicit temporary directory
                logger.debug("Sorting GFF file...")
                sort_cmd = [
                    'sort', 
                    '-k1,1', '-k4,4n',  # Sort by chromosome, then by start position
                    '-T', self.temp_dir,  # Specify temporary directory
                    '-o', temp_sorted_path,  # Output to temporary file
                    gff_file
                ]
                
                logger.debug(f"Sort command: {' '.join(sort_cmd)}")
                sort_result = subprocess.run(sort_cmd, capture_output=True, text=True, timeout=300)
                
                if sort_result.returncode != 0:
                    error_msg = f"Failed to sort GFF file: {sort_result.stderr}"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="sort")
                
                # Verify sorted file was created
                if not os.path.exists(temp_sorted_path) or os.path.getsize(temp_sorted_path) == 0:
                    error_msg = "Sorted GFF file was not created or is empty"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="sort")
                
                logger.debug(f"Successfully sorted GFF file: {temp_sorted_path}")
                logger.debug(f"Sorted file size: {os.path.getsize(temp_sorted_path)} bytes")
                
                # Now compress with bgzip
                logger.debug("Compressing sorted GFF file...")
                compress_cmd = ['bgzip', '-c', temp_sorted_path]
                
                logger.debug(f"Compress command: {' '.join(compress_cmd)}")
                with open(output_path, 'wb') as output_file:
                    compress_result = subprocess.run(compress_cmd, stdout=output_file,
                                                stderr=subprocess.PIPE, text=False, timeout=300)
                
                if compress_result.returncode != 0:
                    if os.path.exists(output_path):
                        os.unlink(output_path)
                    error_msg = f"Failed to compress GFF file: {compress_result.stderr.decode()}"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="bgzip")
                
                # Verify compressed file was created
                if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                    error_msg = "Compressed GFF file was not created or is empty"
                    logger.error(error_msg)
                    raise ExternalToolError(error_msg, tool_name="bgzip")
                
                logger.debug(f"Successfully compressed GFF file: {output_path}")
                logger.debug(f"Compressed file size: {os.path.getsize(output_path)} bytes")
                
            finally:
                # Clean up temporary sorted file
                if os.path.exists(temp_sorted_path):
                    os.unlink(temp_sorted_path)
            
            logger.debug(f"Successfully sorted and compressed GFF: {output_path}")
            
        except subprocess.TimeoutExpired:
            if os.path.exists(output_path):
                os.unlink(output_path)
            error_msg = "GFF sorting/compression timed out"
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="sort/bgzip")
        except ExternalToolError:
            raise
        except Exception as e:
            if os.path.exists(output_path):
                os.unlink(output_path)
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

    def _validate_prepared_vcf(self, vcf_file: str, original_issues: List[Dict]) -> List[str]:
        """
        Validate that a prepared VCF file actually addresses the original issues.
        
        Args:
            vcf_file: Path to prepared VCF file
            original_issues: List of issues that should have been fixed
            
        Returns:
            List of remaining issues (empty if all fixed)
        """
        remaining_issues = []
        
        try:
            # Check compression
            if any(issue['type'] == 'vcf_compression' for issue in original_issues):
                if not self._is_properly_bgzipped(vcf_file):
                    remaining_issues.append("VCF still not properly bgzipped")
            
            # Check indexing
            if any(issue['type'] == 'vcf_index' for issue in original_issues):
                if not self._is_vcf_indexed(vcf_file):
                    remaining_issues.append("VCF still not indexed")
            
            # Check AF field
            if any(issue['type'] == 'vcf_af_field' for issue in original_issues):
                if not self._has_af_field(vcf_file):
                    remaining_issues.append("VCF still missing INFO/AF field")
            
            # Check normalization
            if any(issue['type'] == 'vcf_normalization' for issue in original_issues):
                if not self._is_vcf_normalized(vcf_file):
                    remaining_issues.append("VCF still not normalized")
            
            # Check chromosome mapping
            mapping_issues = [issue for issue in original_issues if issue['type'] == 'chromosome_mapping']
            if mapping_issues:
                # Verify chromosomes were actually renamed
                vcf_chroms = set(self._get_vcf_chromosomes(vcf_file).keys())
                mapping = mapping_issues[0].get('mapping', {})
                expected_chroms = set(mapping.values()) if mapping else set()
                
                if expected_chroms and not expected_chroms.issubset(vcf_chroms):
                    remaining_issues.append("Chromosome mapping was not applied correctly")
            
            return remaining_issues
            
        except Exception as e:
            logger.warning(f"Error validating prepared VCF: {str(e)}")
            return ["Could not validate prepared VCF file"]


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