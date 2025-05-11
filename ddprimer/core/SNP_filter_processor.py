#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP Filter Processor for ddPrimer pipeline.

This module contains the implementation of post-design SNP filtering:
1. Extract variants from VCF files
2. Filter primers that contain SNPs
3. Optionally filter amplicons with high SNP density in strict mode
"""

import os
import logging
import subprocess
import gzip
import pandas as pd
from tqdm import tqdm

from ..config import Config
from ..config.exceptions import FileError, SequenceProcessingError

class SNPFilterProcessor:
    """Handles filtering of primers based on SNP positions after primer design."""
    
    def __init__(self, strict_mode=False, amplicon_threshold=0.05):
        """
        Initialize SNP filter processor.
        
        Args:
            strict_mode (bool): Enable strict filtering including amplicon check
            amplicon_threshold (float): Maximum allowed SNP density in amplicons (0.05 = 5%)
        """
        self.strict_mode = strict_mode
        self.amplicon_threshold = amplicon_threshold
        self.logger = logging.getLogger("ddPrimer.snp_filter_processor")
    
    def get_variant_positions(self, vcf_file, chromosome=None):
        """
        Extract variant positions from the VCF file using bcftools.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        self.logger.debug(f"\nFetching variant positions from {vcf_file}")
        
        if not os.path.exists(vcf_file):
            self.logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
        
        # Base command
        command = f'bcftools query -f "%CHROM\\t%POS\\n" "{vcf_file}"'
        
        # Add chromosome filter if specified
        if chromosome:
            command += f' -r "{chromosome}"'
            self.logger.debug(f"Filtering variants for chromosome: {chromosome}")
        
        # Run command
        self.logger.debug(f"Running command: {command}")
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                self.logger.error(f"Error running bcftools: {result.stderr}")
                raise FileError(f"bcftools query failed: {result.stderr}")
            
            # Parse results
            variants = {}
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                    
                parts = line.split("\t")
                if len(parts) == 2:
                    chrom, pos = parts
                    try:
                        pos = int(pos)
                        
                        if chrom not in variants:
                            variants[chrom] = set()
                            
                        variants[chrom].add(pos)
                    except ValueError:
                        self.logger.warning(f"Invalid position value in VCF: {pos}")
            
            total_variants = sum(len(positions) for positions in variants.values())
            self.logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes")
            return variants
            
        except subprocess.SubprocessError as e:
            self.logger.error(f"Failed to run bcftools: {e}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            # Fall back to manual parsing
            self.logger.info("Falling back to manual VCF parsing")
            return self._extract_variants_manually(vcf_file, chromosome)
    
    def _extract_variants_manually(self, vcf_file, chromosome=None):
        """
        Extract variants from a VCF file manually (fallback method).
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
        """
        variants = {}
        
        try:
            # Determine if we need to handle gzip compression
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with open_func(vcf_file, mode) as vcf:
                for line in vcf:
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Parse VCF data line
                    fields = line.strip().split('\t')
                    if len(fields) < 5:  # Minimum valid VCF line
                        continue
                    
                    # Check chromosome match if filtering by chromosome
                    vcf_chrom = fields[0]
                    if chromosome and vcf_chrom != chromosome:
                        continue
                    
                    # Get position
                    try:
                        pos = int(fields[1])
                        if vcf_chrom not in variants:
                            variants[vcf_chrom] = set()
                        variants[vcf_chrom].add(pos)
                    except ValueError:
                        self.logger.warning(f"Invalid position value in VCF: {fields[1]}")
                        continue
            
            total_variants = sum(len(positions) for positions in variants.values())
            self.logger.info(f"Manually extracted {total_variants} variants from VCF file")
            
        except Exception as e:
            self.logger.error(f"Error in manual VCF parsing: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Failed to extract variants manually: {str(e)}")
        
        return variants
    
    def extract_variants_by_regions(self, vcf_file, regions):
        """
        Extract variants from VCF file only for specific regions.
        
        Args:
            vcf_file (str): Path to VCF file
            regions (dict): Dictionary mapping chromosome names to lists of regions
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
        """
        self.logger.debug(f"Extracting variants for specific regions from {vcf_file}")
        
        if not os.path.exists(vcf_file):
            self.logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
            
        if not regions:
            self.logger.warning("No regions provided")
            return {}
        
        # Dictionary to store variants by chromosome
        variants = {}
        
        # Process each chromosome and its regions
        for chrom, chrom_regions in regions.items():
            if not chrom_regions:
                continue
                
            variants[chrom] = set()
            self.logger.debug(f"Processing {len(chrom_regions)} regions for {chrom}")
            
            # Process each region
            for region in chrom_regions:
                start = region.get('start', 0)
                end = region.get('end', 0)
                
                # Skip very small regions
                if end - start < 10:
                    self.logger.debug(f"Skipping small region {chrom}:{start}-{end} (< 10 bp)")
                    continue
                    
                # Extract variants for this specific region
                region_variants = self.get_region_variants(vcf_file, chrom, start, end)
                variants[chrom].update(region_variants)
        
        # Count total variants
        total_variants = sum(len(positions) for positions in variants.values())
        self.logger.info(f"Extracted {total_variants} variants for specified regions across {len(variants)} chromosomes")
        
        return variants

    def get_region_variants(self, vcf_file, chromosome, start_pos, end_pos):
        """
        Extract variant positions from the VCF file for a specific genomic region.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str): Chromosome name
            start_pos (int): Start position of the region
            end_pos (int): End position of the region
            
        Returns:
            set: Set of variant positions within the specified region
        """
        self.logger.debug(f"Fetching variants for region {chromosome}:{start_pos}-{end_pos}")
        
        # Use bcftools to query the region
        command = f'bcftools query -f "%POS\\n" -r "{chromosome}:{start_pos}-{end_pos}" "{vcf_file}"'
        
        self.logger.debug(f"Running command: {command}")
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                # Parse results from bcftools
                variants = set()
                for line in result.stdout.strip().split("\n"):
                    if line:  # Skip empty lines
                        try:
                            pos = int(line.strip())
                            variants.add(pos)
                        except ValueError:
                            self.logger.warning(f"Invalid position value: {line}")
                            continue
                
                self.logger.debug(f"Extracted {len(variants)} variants in region {chromosome}:{start_pos}-{end_pos}")
                return variants
            else:
                self.logger.warning(f"bcftools query failed: {result.stderr}")
                # Fall back to manual extraction
        except Exception as e:
            self.logger.warning(f"Error running bcftools: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
        
        # Manual extraction fallback
        self.logger.debug("Falling back to manual region extraction")
        return self._extract_region_variants_manually(vcf_file, chromosome, start_pos, end_pos)
    
    def _extract_region_variants_manually(self, vcf_file, chrom, start, end):
        """
        Extract variants from a VCF file for a specific region using manual parsing.
        
        Args:
            vcf_file (str): Path to VCF file
            chrom (str): Chromosome name
            start (int): Start position
            end (int): End position
            
        Returns:
            set: Set of variant positions in the region
        """
        variants = set()
        
        try:
            # Determine if we need to handle gzip compression
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with open_func(vcf_file, mode) as vcf:
                for line in vcf:
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Parse VCF data line
                    fields = line.strip().split('\t')
                    if len(fields) < 5:  # Minimum valid VCF line
                        continue
                    
                    # Check chromosome match
                    vcf_chrom = fields[0]
                    if vcf_chrom != chrom:
                        continue
                    
                    # Check position is in our region
                    try:
                        pos = int(fields[1])
                        if start <= pos <= end:
                            variants.add(pos)
                    except ValueError:
                        self.logger.warning(f"Invalid position value in VCF: {fields[1]}")
                        continue
            
            self.logger.debug(f"Manually extracted {len(variants)} variants in region {chrom}:{start}-{end}")
            
        except Exception as e:
            self.logger.error(f"Error in manual VCF region parsing for {chrom}:{start}-{end}: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
        
        return variants
        
    def filter_primers_by_snp(self, primers_df, variants, original_sequences):
        """
        Filter primers based on SNP positions.
        
        Args:
            primers_df (pandas.DataFrame): DataFrame with primer information
            variants (dict): Dictionary mapping chromosomes to sets of variant positions
            original_sequences (dict): Dictionary of original (unmasked) sequences
            
        Returns:
            pandas.DataFrame: Filtered primer DataFrame
        """
        self.logger.info(f"Filtering primers based on SNP positions (strict mode: {self.strict_mode})")
        filtered_df = primers_df.copy()
        
        # Skip if no variants provided
        if not variants:
            self.logger.warning("No variants provided for SNP filtering. Returning unfiltered primers.")
            return filtered_df
        
        rows_to_drop = []
        
        # Use progress bar for filtering if enabled
        if Config.SHOW_PROGRESS:
            row_iterator = tqdm(filtered_df.iterrows(), total=len(filtered_df), desc="Filtering primers for SNPs")
        else:
            row_iterator = filtered_df.iterrows()
        
        for idx, row in row_iterator:
            chromosome = row.get("Chromosome", "")
            if not chromosome or chromosome not in variants:
                continue
                
            chrom_variants = variants[chromosome]
            if not chrom_variants:
                continue
            
            # Default fragment offset
            fragment_start = 0
            
            # Get primer positions
            f_start = row.get("Primer F Start", None)
            f_len = row.get("Primer F Len", None)
            r_start = row.get("Primer R Start", None)
            r_len = row.get("Primer R Len", None)
            
            # Check if primers overlap with SNPs
            has_snp = False
            
            # Check forward primer
            if f_start is not None and f_len is not None:
                f_end = f_start + f_len - 1
                
                for pos in range(f_start, f_end + 1):
                    # Convert to absolute genomic position if needed
                    abs_pos = pos + fragment_start
                    if abs_pos in chrom_variants:
                        self.logger.debug(f"Forward primer at idx {idx} overlaps with SNP at position {abs_pos}")
                        has_snp = True
                        break
                        
            # Check reverse primer if forward is clean
            if not has_snp and r_start is not None and r_len is not None:
                r_end = r_start
                r_begin = r_start - r_len + 1
                
                for pos in range(r_begin, r_end + 1):
                    # Convert to absolute genomic position if needed
                    abs_pos = pos + fragment_start
                    if abs_pos in chrom_variants:
                        self.logger.debug(f"Reverse primer at idx {idx} overlaps with SNP at position {abs_pos}")
                        has_snp = True
                        break
                        
            # Check probe if present and previous checks are clean
            if not has_snp and "Probe Start" in row and "Probe Len" in row:
                p_start = row.get("Probe Start", None)
                p_len = row.get("Probe Len", None)
                
                if p_start is not None and p_len is not None and p_start > 0 and p_len > 0:
                    p_end = p_start + p_len - 1
                    
                    for pos in range(p_start, p_end + 1):
                        # Convert to absolute genomic position if needed
                        abs_pos = pos + fragment_start
                        if abs_pos in chrom_variants:
                            self.logger.debug(f"Probe at idx {idx} overlaps with SNP at position {abs_pos}")
                            has_snp = True
                            break
            
            # For strict mode, check entire amplicon
            if not has_snp and self.strict_mode:
                amplicon = row.get("Amplicon", "")
                if amplicon and f_start is not None and r_start is not None:
                    # Amplicon extends from start of forward primer to end of reverse primer
                    amp_start = f_start + fragment_start
                    amp_end = r_start + fragment_start
                    
                    # Get all SNPs in amplicon region
                    amp_snps = [pos for pos in chrom_variants if amp_start <= pos <= amp_end]
                    amp_snp_count = len(amp_snps)
                    
                    # Calculate SNP density
                    amp_length = len(amplicon)
                    snp_density = amp_snp_count / amp_length if amp_length > 0 else 0
                    
                    if snp_density > self.amplicon_threshold:
                        self.logger.debug(f"Amplicon at idx {idx} has high SNP density: {snp_density:.2%} ({amp_snp_count}/{amp_length})")
                        has_snp = True
            
            # Mark row for dropping if it has SNPs
            if has_snp:
                rows_to_drop.append(idx)
        
        # Drop rows with SNPs
        filtered_df = filtered_df.drop(rows_to_drop)
        self.logger.info(f"Removed {len(rows_to_drop)} primers that contain SNPs")
        self.logger.info(f"Retained {len(filtered_df)} primers after SNP filtering")
        
        return filtered_df