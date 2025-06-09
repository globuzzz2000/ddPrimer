#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chromosome Mapper module for ddPrimer pipeline.

This module provides intelligent chromosome name mapping between VCF and FASTA files
without relying on hardcoded organism-specific mappings.
"""

import subprocess
import gzip
import os
import re
import logging
from collections import defaultdict

# Import package modules
from ..config import Config, FileError, SequenceProcessingError

# Set up logging
logger = logging.getLogger(__name__)


class ChromosomeMapper:
    """Handles intelligent chromosome name mapping between VCF and FASTA files."""
    
    def __init__(self):
        """Initialize chromosome mapper."""
        logger.debug("Initialized ChromosomeMapper")
    
    def get_vcf_chromosomes(self, vcf_file):
        """
        Extract chromosome names and variant counts from VCF file.
        
        Args:
            vcf_file (str): Path to VCF file
            
        Returns:
            dict: Dictionary mapping chromosome names to variant counts
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        logger.debug(f"Extracting chromosome information from VCF: {vcf_file}")
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
        
        try:
            cmd = ['bcftools', 'query', '-f', '%CHROM\n', vcf_file]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                chromosomes = result.stdout.strip().split('\n')
                chrom_counts = defaultdict(int)
                for chrom in chromosomes:
                    if chrom:  # Skip empty lines
                        chrom_counts[chrom] += 1
                
                logger.debug(f"Found {len(chrom_counts)} unique chromosomes in VCF")
                return dict(chrom_counts)
            else:
                raise subprocess.SubprocessError(f"bcftools failed: {result.stderr}")
                
        except subprocess.SubprocessError:
            raise
        except Exception as e:
            logger.error(f"Error reading VCF chromosomes: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Failed to read VCF chromosomes: {e}")
    
    def get_fasta_sequences(self, fasta_file):
        """
        Extract sequence names and lengths from FASTA file.
        
        Args:
            fasta_file (str): Path to FASTA file
            
        Returns:
            dict: Dictionary mapping sequence names to lengths
            
        Raises:
            FileError: If FASTA file cannot be processed
        """
        logger.debug(f"Extracting sequence information from FASTA: {fasta_file}")
        
        if not os.path.exists(fasta_file):
            logger.error(f"FASTA file not found: {fasta_file}")
            raise FileError(f"FASTA file not found: {fasta_file}")
        
        try:
            sequences = {}
            current_seq = None
            current_length = 0
            
            opener = gzip.open if fasta_file.endswith('.gz') else open
            mode = 'rt' if fasta_file.endswith('.gz') else 'r'
                
            with opener(fasta_file, mode) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_seq is not None:
                            sequences[current_seq] = current_length
                        
                        seq_id = line[1:].split()[0]
                        current_seq = seq_id
                        current_length = 0
                    elif line:
                        current_length += len(line)
                
                # Don't forget the last sequence
                if current_seq is not None:
                    sequences[current_seq] = current_length
            
            logger.debug(f"Found {len(sequences)} sequences in FASTA")
            return sequences
            
        except Exception as e:
            logger.error(f"Error reading FASTA sequences: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Failed to read FASTA sequences: {e}")
    
    def extract_numeric_component(self, name):
        """
        Extract numeric component from chromosome name for intelligent sorting.
        
        Args:
            name (str): Chromosome or sequence name
            
        Returns:
            int: Numeric component if found, large number otherwise
        """
        # Handle common chromosome naming patterns
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
        
        # Pattern 4: For accession numbers like CP002684.1, extract the main number but ignore version
        # Remove version numbers (anything after last dot)
        base_name = name.split('.')[0] if '.' in name else name
        
        # Extract all numbers from the base name
        numbers = re.findall(r'\d+', base_name)
        if numbers:
            # Use the last/longest number found (usually the main identifier)
            return int(numbers[-1])
        
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
                return value
        
        return float('inf')  # Put unrecognized names at end
    
    def suggest_automatic_mapping(self, vcf_chroms, fasta_seqs):
        """
        Suggest intelligent mapping based on analysis of both files.
        
        Args:
            vcf_chroms (dict): VCF chromosomes and their variant counts
            fasta_seqs (dict): FASTA sequences and their lengths
            
        Returns:
            dict: Suggested mapping from VCF chromosome to FASTA sequence
        """
        logger.debug("Analyzing files for automatic chromosome mapping")
        
        # Sort VCF chromosomes by numeric component, then alphabetically
        vcf_sorted = sorted(vcf_chroms.keys(), 
                           key=lambda x: (self.extract_numeric_component(x), x))
        
        # Filter FASTA sequences to exclude likely organellar genomes
        nuclear_fasta = self._filter_nuclear_chromosomes(fasta_seqs)
        
        # Sort nuclear FASTA sequences by numeric component, then alphabetically
        fasta_sorted = sorted(nuclear_fasta.keys(), 
                             key=lambda x: (self.extract_numeric_component(x), x))
        
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
            
            logger.debug("Suggested mapping (sorted order matching):")
            for vcf_chr, fasta_chr in mapping.items():
                logger.debug(f"  '{vcf_chr}' -> '{fasta_chr}'")
            
            if len(vcf_sorted) < len(fasta_sorted):
                unused_count = len(fasta_sorted) - len(vcf_sorted)
                logger.debug(f"Note: {unused_count} nuclear FASTA sequences will not be used")
            
            return mapping
        
        else:
            # VCF has more chromosomes than nuclear FASTA sequences - problematic
            logger.warning(f"Problem: VCF has {len(vcf_sorted)} chromosomes but only {len(fasta_sorted)} nuclear FASTA sequences")
            logger.warning("Cannot create automatic mapping - manual intervention required")
            return {}
    
    
    def _filter_nuclear_chromosomes(self, fasta_seqs):
        """
        Filter FASTA sequences to keep only likely nuclear chromosomes.
        
        Args:
            fasta_seqs (dict): Dictionary mapping sequence names to lengths
            
        Returns:
            dict: Filtered dictionary with only nuclear chromosomes
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
        
        logger.debug(f"Filtered to {len(nuclear_seqs)} nuclear sequences from {len(fasta_seqs)} total sequences")
        return nuclear_seqs
    

    
    def check_chromosome_compatibility(self, vcf_file, fasta_file):
        """
        Check if VCF and FASTA files have compatible chromosome names.
        
        Args:
            vcf_file (str): Path to VCF file
            fasta_file (str): Path to FASTA file
            
        Returns:
            dict: Analysis results with compatibility info and suggested actions
            
        Raises:
            FileError: If files cannot be processed
        """
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
            
            return analysis
            
        except Exception as e:
            logger.error(f"Error in compatibility check: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Compatibility check failed: {e}")