#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chromosome Mapper module for ddPrimer pipeline.

Contains functionality for:
1. Intelligent chromosome name mapping between VCF and FASTA files
2. Chromosome compatibility analysis and validation
3. Automatic mapping suggestions based on sequence analysis
4. Nuclear vs organellar chromosome filtering

This module provides intelligent chromosome name mapping between VCF and FASTA files
without relying on hardcoded organism-specific mappings, enabling cross-species
compatibility and flexible file format support.
"""

import subprocess
import gzip
import os
import re
import logging
from collections import defaultdict
from ..config import Config, FileError, SequenceProcessingError, ExternalToolError

# Set up module logger
logger = logging.getLogger(__name__)


class ChromosomeMapper:
    """
    Handles intelligent chromosome name mapping between VCF and FASTA files.
    
    This class provides methods for analyzing chromosome naming patterns in both
    VCF and FASTA files, detecting compatibility issues, and suggesting automatic
    mappings based on sequence characteristics and naming conventions.
    
    Attributes:
        None - stateless utility class
        
    Example:
        >>> mapper = ChromosomeMapper()
        >>> analysis = mapper.check_chromosome_compatibility("variants.vcf", "genome.fasta")
        >>> if analysis['needs_mapping']:
        ...     mapping = analysis['suggested_mapping']
        ...     print(f"Suggested mapping: {mapping}")
    """
    
    def __init__(self):
        """
        Initialize chromosome mapper.
        
        Creates a new ChromosomeMapper instance for analyzing chromosome
        compatibility between VCF and FASTA files.
        """
        logger.debug("Initialized ChromosomeMapper")
    
    def get_vcf_chromosomes(self, vcf_file):
        """
        Extract chromosome names and variant counts from VCF file.
        
        Uses bcftools to efficiently extract chromosome information from
        VCF files, supporting both compressed and uncompressed formats.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            Dictionary mapping chromosome names to variant counts
            
        Raises:
            FileError: If VCF file cannot be accessed or processed
            ExternalToolError: If bcftools execution fails
        """
        logger.debug("=== VCF CHROMOSOME EXTRACTION DEBUG ===")
        logger.debug(f"Extracting chromosome information from VCF: {vcf_file}")
        
        if not os.path.exists(vcf_file):
            error_msg = f"VCF file not found: {vcf_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        try:
            cmd = ['bcftools', 'query', '-f', '%CHROM\n', vcf_file]
            cmd_str = ' '.join(cmd)
            logger.debug(f"Running command: {cmd_str}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                chromosomes = result.stdout.strip().split('\n')
                chrom_counts = defaultdict(int)
                
                for chrom in chromosomes:
                    if chrom and chrom.strip():  # Skip empty lines
                        chrom_counts[chrom.strip()] += 1
                
                logger.debug(f"Found {len(chrom_counts)} unique chromosomes in VCF")
                if logger.isEnabledFor(logging.DEBUG):
                    for chrom, count in sorted(chrom_counts.items()):
                        logger.debug(f"  {chrom}: {count} variants")
                
                logger.debug("=== END VCF CHROMOSOME EXTRACTION DEBUG ===")
                return dict(chrom_counts)
            else:
                error_msg = f"bcftools execution failed for {vcf_file}"
                logger.error(error_msg)
                logger.debug(f"bcftools stderr: {result.stderr}")
                logger.debug("=== END VCF CHROMOSOME EXTRACTION DEBUG ===")
                raise ExternalToolError(error_msg, tool_name="bcftools") from subprocess.SubprocessError(result.stderr)
                
        except ExternalToolError:
            # Re-raise ExternalToolError without modification
            raise
        except subprocess.SubprocessError as e:
            error_msg = f"bcftools subprocess error: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END VCF CHROMOSOME EXTRACTION DEBUG ===")
            raise ExternalToolError(error_msg, tool_name="bcftools") from e
        except Exception as e:
            error_msg = f"Error reading VCF chromosomes from {vcf_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END VCF CHROMOSOME EXTRACTION DEBUG ===")
            raise FileError(error_msg) from e
    
    def get_fasta_sequences(self, fasta_file):
        """
        Extract sequence names and lengths from FASTA file.
        
        Parses FASTA files to extract sequence identifiers and calculate
        sequence lengths, supporting both compressed and uncompressed formats.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Dictionary mapping sequence names to lengths
            
        Raises:
            FileError: If FASTA file cannot be accessed or parsed
        """
        logger.debug("=== FASTA SEQUENCE EXTRACTION DEBUG ===")
        logger.debug(f"Extracting sequence information from FASTA: {fasta_file}")
        
        if not os.path.exists(fasta_file):
            error_msg = f"FASTA file not found: {fasta_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        try:
            sequences = {}
            current_seq = None
            current_length = 0
            
            # Determine file opener based on extension
            opener = gzip.open if fasta_file.endswith('.gz') else open
            mode = 'rt' if fasta_file.endswith('.gz') else 'r'
            logger.debug(f"Using opener: {opener.__name__}, mode: {mode}")
                
            with opener(fasta_file, mode) as f:
                line_count = 0
                for line in f:
                    line_count += 1
                    line = line.strip()
                    
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_seq is not None:
                            sequences[current_seq] = current_length
                            logger.debug(f"Processed sequence: {current_seq} ({current_length:,} bp)")
                        
                        # Start new sequence
                        seq_id = line[1:].split()[0]  # Take first part after >
                        current_seq = seq_id
                        current_length = 0
                        logger.debug(f"Started new sequence: {current_seq}")
                    elif line:
                        # Count sequence characters
                        current_length += len(line)
                
                # Don't forget the last sequence
                if current_seq is not None:
                    sequences[current_seq] = current_length
                    logger.debug(f"Processed final sequence: {current_seq} ({current_length:,} bp)")
            
            logger.debug(f"Processed {line_count} lines, found {len(sequences)} sequences in FASTA")
            logger.debug("=== END FASTA SEQUENCE EXTRACTION DEBUG ===")
            return sequences
            
        except (OSError, IOError) as e:
            error_msg = f"Error reading FASTA file {fasta_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END FASTA SEQUENCE EXTRACTION DEBUG ===")
            raise FileError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error reading FASTA sequences from {fasta_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END FASTA SEQUENCE EXTRACTION DEBUG ===")
            raise FileError(error_msg) from e
    
    def extract_numeric_component(self, name):
        """
        Extract numeric component from chromosome name for intelligent sorting.
        
        Analyzes chromosome names to extract numeric components for proper
        sorting, handling various naming conventions including Chr1, Chromosome1,
        accession numbers, and special chromosomes.
        
        Args:
            name: Chromosome or sequence name
            
        Returns:
            Numeric component if found, large number for special/unrecognized names
            
        Example:
            >>> mapper = ChromosomeMapper()
            >>> mapper.extract_numeric_component("Chr5")
            5
            >>> mapper.extract_numeric_component("ChrX")
            100
        """
        if not name or not isinstance(name, str):
            logger.debug(f"Invalid name provided to extract_numeric_component: {name}")
            return float('inf')
            
        name_upper = name.upper()
        logger.debug(f"Extracting numeric component from: '{name}'")
        
        # Pattern 1: Simple numbers (1, 2, 3, etc.)
        if name.isdigit():
            result = int(name)
            logger.debug(f"  Pattern 1 (simple number): {result}")
            return result
        
        # Pattern 2: Chr1, Chr2, etc.
        if name_upper.startswith('CHR') and len(name) > 3:
            chr_part = name[3:]
            if chr_part.isdigit():
                result = int(chr_part)
                logger.debug(f"  Pattern 2 (Chr prefix): {result}")
                return result
        
        # Pattern 3: Chromosome1, Chromosome2, etc.
        if name_upper.startswith('CHROMOSOME') and len(name) > 10:
            chr_part = name[10:]
            if chr_part.isdigit():
                result = int(chr_part)
                logger.debug(f"  Pattern 3 (Chromosome prefix): {result}")
                return result
        
        # Pattern 4: For accession numbers like CP002684.1, extract the main number but ignore version
        base_name = name.split('.')[0] if '.' in name else name
        numbers = re.findall(r'\d+', base_name)
        if numbers:
            # Use the last/longest number found (usually the main identifier)
            result = int(numbers[-1])
            logger.debug(f"  Pattern 4 (accession number): {result} from {numbers}")
            return result
        
        # Pattern 5: Special chromosomes (X, Y, MT, etc.)
        special_chroms = {
            'X': 100,
            'Y': 101, 
            'MT': 102,
            'MITO': 102,
            'MITOCHONDRIAL': 102,
            'CHLOROPLAST': 103,
            'PLASTID': 103,
            'PT': 103,
            'CP': 103
        }
        
        for special, value in special_chroms.items():
            if special in name_upper:
                logger.debug(f"  Pattern 5 (special chromosome): {value} for '{special}'")
                return value
        
        logger.debug(f"  No pattern matched, using infinity for sorting")
        return float('inf')  # Put unrecognized names at end
    
    def suggest_automatic_mapping(self, vcf_chroms, fasta_seqs):
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
        logger.debug("=== AUTOMATIC MAPPING SUGGESTION DEBUG ===")
        logger.debug("Analyzing files for automatic chromosome mapping")
        
        # Sort VCF chromosomes by numeric component, then alphabetically
        vcf_sorted = sorted(vcf_chroms.keys(), 
                           key=lambda x: (self.extract_numeric_component(x), x))
        
        # Filter FASTA sequences to exclude likely organellar genomes
        nuclear_fasta = self._filter_nuclear_chromosomes(fasta_seqs)
        
        # Sort nuclear FASTA sequences by numeric component, then alphabetically
        fasta_sorted = sorted(nuclear_fasta.keys(), 
                             key=lambda x: (self.extract_numeric_component(x), x))
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("VCF chromosomes (sorted by numeric order):")
            for i, chrom in enumerate(vcf_sorted, 1):
                count = vcf_chroms[chrom]
                numeric = self.extract_numeric_component(chrom)
                logger.debug(f"  {i}. '{chrom}' (numeric: {numeric}, {count:,} variants)")
            
            logger.debug("Nuclear FASTA sequences (sorted by numeric order):")
            for i, seq in enumerate(fasta_sorted, 1):
                length = nuclear_fasta[seq]
                numeric = self.extract_numeric_component(seq)
                logger.debug(f"  {i}. '{seq}' (numeric: {numeric}, {length:,} bp)")
            
            # Show excluded sequences
            excluded = {k: v for k, v in fasta_seqs.items() if k not in nuclear_fasta}
            if excluded:
                logger.debug("Excluded organellar/small sequences:")
                for seq, length in excluded.items():
                    logger.debug(f"  '{seq}' ({length:,} bp)")
        
        # Main strategy: Map chromosomes in order after sorting both by numeric component
        if len(vcf_sorted) <= len(fasta_sorted):
            # Take the first N nuclear FASTA sequences
            main_fasta = fasta_sorted[:len(vcf_sorted)]
            mapping = dict(zip(vcf_sorted, main_fasta))
            
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("Suggested mapping (sorted order matching):")
                for vcf_chr, fasta_chr in mapping.items():
                    logger.debug(f"  '{vcf_chr}' -> '{fasta_chr}'")
            
            if len(vcf_sorted) < len(fasta_sorted):
                unused_count = len(fasta_sorted) - len(vcf_sorted)
                logger.debug(f"Note: {unused_count} nuclear FASTA sequences will not be used")
            
            logger.debug("=== END AUTOMATIC MAPPING SUGGESTION DEBUG ===")
            return mapping
        
        else:
            # VCF has more chromosomes than nuclear FASTA sequences - problematic
            logger.warning(f"Problem: VCF has {len(vcf_sorted)} chromosomes but only {len(fasta_sorted)} nuclear FASTA sequences")
            logger.warning("Cannot create automatic mapping - manual intervention required")
            logger.debug("=== END AUTOMATIC MAPPING SUGGESTION DEBUG ===")
            return {}
    
    def _filter_nuclear_chromosomes(self, fasta_seqs):
        """
        Filter FASTA sequences to keep only likely nuclear chromosomes.
        
        Identifies and filters out organellar genomes (mitochondrial, chloroplast)
        and small sequences that are unlikely to be nuclear chromosomes based
        on size and naming patterns.
        
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
            
            # Also check size - organellar genomes are typically < 1 Mb
            is_too_small = length < MIN_NUCLEAR_SIZE
            
            if is_organellar:
                logger.debug(f"Excluding organellar sequence: {seq_name} ({length:,} bp)")
            elif is_too_small:
                logger.debug(f"Excluding small sequence (likely organellar): {seq_name} ({length:,} bp)")
            else:
                nuclear_seqs[seq_name] = length
                logger.debug(f"Including nuclear sequence: {seq_name} ({length:,} bp)")
        
        logger.debug(f"Filtered to {len(nuclear_seqs)} nuclear sequences from {len(fasta_seqs)} total sequences")
        return nuclear_seqs
    
    def check_chromosome_compatibility(self, vcf_file, fasta_file):
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
            
        Raises:
            FileError: If files cannot be processed
            ExternalToolError: If external tools fail
        """
        logger.debug("=== CHROMOSOME COMPATIBILITY CHECK DEBUG ===")
        logger.debug("Checking chromosome compatibility between VCF and FASTA files")
        
        try:
            # Get chromosome information from both files
            vcf_chroms = self.get_vcf_chromosomes(vcf_file)
            fasta_seqs = self.get_fasta_sequences(fasta_file)
            
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
                logger.debug(f"{len(exact_matches)} chromosome(s) match exactly:")
                for chrom in sorted(exact_matches):
                    logger.debug(f"  {chrom}")
            
            if vcf_only:
                logger.debug(f"{len(vcf_only)} chromosome(s) only in VCF:")
                for chrom in sorted(vcf_only):
                    count = vcf_chroms[chrom]
                    logger.debug(f"  {chrom} ({count:,} variants)")
            
            if fasta_only:
                logger.debug(f"{len(fasta_only)} sequence(s) only in FASTA:")
                for seq in sorted(fasta_only):
                    length = fasta_seqs[seq]
                    logger.debug(f"  {seq} ({length:,} bp)")
            
            # Provide recommendations
            if analysis['compatible'] and not analysis['needs_mapping']:
                logger.debug("Files are fully compatible - no chromosome renaming needed")
                analysis['action'] = 'none'
            elif analysis['needs_mapping']:
                logger.debug("Files need chromosome name mapping")
                analysis['action'] = 'mapping'
                # Generate automatic mapping suggestion
                analysis['suggested_mapping'] = self.suggest_automatic_mapping(vcf_chroms, fasta_seqs)
            else:
                logger.debug("No compatible chromosomes found")
                analysis['action'] = 'error'
            
            logger.debug("=== END CHROMOSOME COMPATIBILITY CHECK DEBUG ===")
            return analysis
            
        except (FileError, ExternalToolError):
            # Re-raise specific exceptions without modification
            raise
        except Exception as e:
            error_msg = f"Error in compatibility check between {vcf_file} and {fasta_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END CHROMOSOME COMPATIBILITY CHECK DEBUG ===")
            raise SequenceProcessingError(error_msg) from e