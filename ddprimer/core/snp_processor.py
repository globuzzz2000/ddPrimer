#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP masking processor for ddPrimer pipeline.

Handles masking of SNPs in sequences to prepare them for primer design.
Provides functionality for:
1. VCF file processing with chromosome mapping
2. Variant extraction with quality and allele frequency filtering
3. Sequence masking with configurable parameters
4. Fixed SNP substitution - automatically substituting high-frequency variants (AF=1.0 or above threshold)
5. Region-specific variant extraction for memory optimization

This module integrates with the broader ddPrimer pipeline to provide
robust SNP masking capabilities for primer design workflows.
"""

import os
import logging
import subprocess
import gzip
from Bio import SeqIO
from tqdm import tqdm
from collections import namedtuple

# Import package modules
from ..config import Config, FileError, SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)

# Define variant structure
Variant = namedtuple('Variant', ['position', 'ref', 'alt', 'qual', 'af', 'is_fixed'])


class SNPMaskingProcessor:
    """
    Handles masking of SNPs in sequences to prepare them for primer design.
    
    This class provides methods for processing VCF files, extracting variants
    with quality filtering, and applying masking to sequences for primer design.
    Automatically substitutes high-frequency variants (AF=1.0 or above threshold) 
    and masks variable variants.
    """
    
    def __init__(self):
        """Initialize SNP masking processor."""
        logger.debug("Initialized SNPMaskingProcessor with fixed SNP support")
    
    def get_variants(self, vcf_file: str, chromosome: str = None, 
                    min_af: float = None, min_qual: float = None) -> dict:
        """
        Extract variant information from VCF file with allele frequency and quality filtering.
        
        Args:
            vcf_file: Path to VCF file
            chromosome: Specific chromosome to filter
            min_af: Minimum allele frequency threshold (0.0-1.0)
            min_qual: Minimum QUAL score threshold
            
        Returns:
            Dictionary mapping chromosomes to lists of Variant objects
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        logger.debug(f"Fetching variants from {vcf_file}")
        logger.debug(f"Filters: min_af={min_af}, min_qual={min_qual}, chromosome={chromosome}")
        
        if not os.path.exists(vcf_file):
            error_msg = f"VCF file not found: {vcf_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        # Try bcftools approach first
        try:
            return self._extract_with_bcftools(vcf_file, chromosome, min_af, min_qual)
        except Exception as e:
            logger.warning(f"bcftools extraction failed: {str(e)}")
            logger.info("Falling back to manual VCF parsing")
            return self._extract_variants_manually(vcf_file, chromosome, min_af, min_qual)
    
    def _extract_with_bcftools(self, vcf_file: str, chromosome: str = None, 
                            min_af: float = None, min_qual: float = None) -> dict:
        """
        FIXED: Extract variants using bcftools with proper filtering.
        
        Key fixes:
        - Better case handling in variant creation
        - More robust AF parsing
        """
        
        # Build filters carefully
        filters = []
        
        # Quality filter
        if min_qual is not None:
            filters.append(f"QUAL>={min_qual}")
        
        # For AF filtering, we'll include it in the filter but handle errors gracefully
        if min_af is not None:
            filters.append(f"INFO/AF>={min_af}")
        
        # Build command to get detailed information
        if filters:
            filter_string = " && ".join(filters)
            command = f'bcftools query -i "{filter_string}" -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        else:
            command = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        
        # Add chromosome filter if specified
        if chromosome:
            command += f' -r "{chromosome}"'
            logger.debug(f"Filtering variants for chromosome: {chromosome}")
        
        logger.debug(f"Running bcftools command: {command}")
        
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                # If AF filtering fails, try without AF filter but with manual AF filtering
                if min_af is not None and "INFO/AF" in result.stderr:
                    logger.warning("AF field might not exist in VCF INFO, trying alternative approach")
                    return self._extract_with_bcftools_fallback(vcf_file, chromosome, min_af, min_qual)
                else:
                    raise subprocess.SubprocessError(f"bcftools failed: {result.stderr}")
            
            # Parse results
            variants = {}
            total_processed = 0
            total_kept = 0
            fixed_count = 0
            
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                    
                parts = line.split("\t")
                if len(parts) >= 4:
                    chrom = parts[0]
                    try:
                        pos = int(parts[1])
                        ref = parts[2]  # FIXED: Keep original case from VCF
                        alt = parts[3]  # FIXED: Keep original case from VCF
                        qual = float(parts[4]) if len(parts) > 4 and parts[4] != "." else None
                        af = None
                        
                        if len(parts) > 5 and parts[5] not in [".", ""]:
                            try:
                                af_values = [float(x.replace(",", ".")) for x in parts[5].split(",")]
                                af = sum(af_values)
                            except ValueError:
                                logger.warning(f"Invalid AF value: {parts[5]}")
                        
                        total_processed += 1
                        
                        # Apply manual filters as backup
                        passes_filter = True
                        
                        if min_qual is not None and qual is not None and qual < min_qual:
                            passes_filter = False
                        
                        if min_af is not None and af is not None and af < min_af:
                            passes_filter = False
                        
                        if passes_filter:
                            # Determine if this is a fixed SNP
                            is_fixed = False
                            if af is not None:
                                if af >= 1.0:
                                    is_fixed = True
                                elif hasattr(self, 'Config') and hasattr(self.Config, 'SNP_ALLELE_FREQUENCY_THRESHOLD') and self.Config.SNP_ALLELE_FREQUENCY_THRESHOLD is not None:
                                    fixed_af_threshold = 1.0 - self.Config.SNP_ALLELE_FREQUENCY_THRESHOLD
                                    if af >= fixed_af_threshold:
                                        if not hasattr(self.Config, 'SNP_QUALITY_THRESHOLD') or self.Config.SNP_QUALITY_THRESHOLD is None or (qual is not None and qual > self.Config.SNP_QUALITY_THRESHOLD):
                                            is_fixed = True
                            
                            if is_fixed:
                                fixed_count += 1
                            
                            # Take only the first alternate allele for simplicity
                            alt_allele = alt.split(',')[0] if ',' in alt else alt
                            
                            variant = Variant(
                                position=pos,
                                ref=ref,  # FIXED: Use original case
                                alt=alt_allele,  # FIXED: Use original case
                                qual=qual,
                                af=af,
                                is_fixed=is_fixed
                            )
                            
                            if chrom not in variants:
                                variants[chrom] = []
                            variants[chrom].append(variant)
                            total_kept += 1
                            
                    except ValueError:
                        logger.warning(f"Invalid position value in VCF: {parts[1]}")
                        continue
            
            logger.debug(f"bcftools processing: {total_processed} variants examined, {total_kept} kept after filtering")
            logger.debug(f"Fixed SNPs identified: {fixed_count}")
            if min_af is not None:
                logger.debug(f"AF threshold applied: >= {min_af}")
            if min_qual is not None:
                logger.debug(f"QUAL threshold applied: >= {min_qual}")
            
            return variants
            
        except subprocess.SubprocessError:
            # Re-raise subprocess errors
            raise
        except Exception as e:
            error_msg = f"Error in bcftools processing: {str(e)}"
            logger.warning(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise
    
    def _extract_with_bcftools_fallback(self, vcf_file: str, chromosome: str = None, 
                                       min_af: float = None, min_qual: float = None) -> dict:
        """Fallback bcftools approach when AF field filtering fails."""
        # Try with only QUAL filter first, then manually filter AF
        filters = []
        if min_qual is not None:
            filters.append(f"QUAL>={min_qual}")
        
        if filters:
            filter_string = " && ".join(filters)
            command = f'bcftools query -i "{filter_string}" -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        else:
            command = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        
        if chromosome:
            command += f' -r "{chromosome}"'
        
        logger.debug(f"Fallback bcftools command: {command}")
        
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise subprocess.SubprocessError(f"Fallback bcftools failed: {result.stderr}")
        
        # Parse and manually filter
        variants = {}
        total_processed = 0
        total_kept = 0
        fixed_count = 0
        
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
                
            parts = line.split("\t")
            if len(parts) >= 4:
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                    ref = parts[2]
                    alt = parts[3]
                    qual = float(parts[4]) if len(parts) > 4 and parts[4] != "." else None
                    af_str = parts[5] if len(parts) > 5 else ""
                    
                    total_processed += 1
                    
                    # Parse AF from INFO field if needed
                    if min_af is not None:
                        af = self._parse_af_from_string(af_str)
                        if af is None or af < min_af:
                            continue
                    else:
                        af = self._parse_af_from_string(af_str)
                    
                    # QUAL should already be filtered by bcftools, but double-check
                    if min_qual is not None and qual is not None and qual < min_qual:
                        continue
                    
                    # FIXED: Determine if this is a fixed SNP
                    is_fixed = False
                    if af is not None:
                        if af >= 1.0:
                            # Always substitute perfect frequency variants
                            is_fixed = True
                        elif Config.SNP_ALLELE_FREQUENCY_THRESHOLD is not None:
                            # Use threshold if configured
                            fixed_af_threshold = 1.0 - Config.SNP_ALLELE_FREQUENCY_THRESHOLD
                            if af >= fixed_af_threshold:
                                # Also check quality if threshold is set
                                if Config.SNP_QUALITY_THRESHOLD is None or (qual is not None and qual > Config.SNP_QUALITY_THRESHOLD):
                                    is_fixed = True
                    
                    if is_fixed:
                        fixed_count += 1
                    
                    # Take only the first alternate allele for simplicity
                    alt_allele = alt.split(',')[0] if ',' in alt else alt
                    
                    variant = Variant(
                        position=pos,
                        ref=ref,
                        alt=alt_allele,
                        qual=qual,
                        af=af,
                        is_fixed=is_fixed
                    )
                    
                    if chrom not in variants:
                        variants[chrom] = []
                    variants[chrom].append(variant)
                    total_kept += 1
                    
                except ValueError:
                    logger.warning(f"Invalid position value: {parts[1]}")
                    continue
        
        logger.info(f"bcftools fallback: {total_processed} variants examined, {total_kept} kept after filtering")
        logger.info(f"Fixed SNPs identified: {fixed_count}")
        return variants
    
    def _extract_variants_manually(self, vcf_file: str, chromosome: str = None, 
                                min_af: float = None, min_qual: float = None) -> dict:
        """
        FIXED: Extract variants using manual VCF parsing with filtering applied globally.
        
        Key fixes:
        - Preserve original case from VCF
        - Better error handling
        """
        variants = {}
        total_processed = 0
        total_kept = 0
        fixed_count = 0
        
        try:
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with open_func(vcf_file, mode) as vcf:
                for line in vcf:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue
                    
                    vcf_chrom = fields[0]
                    if chromosome and vcf_chrom != chromosome:
                        continue
                    
                    try:
                        pos = int(fields[1])
                        ref = fields[3]  # FIXED: Keep original case
                        alt = fields[4]  # FIXED: Keep original case
                        total_processed += 1
                        
                        # Apply quality filter
                        qual = None
                        if min_qual is not None:
                            try:
                                qual = float(fields[5])
                                if qual < min_qual:
                                    continue
                            except (ValueError, IndexError):
                                # Skip if QUAL cannot be parsed and filter is required
                                continue
                        else:
                            try:
                                qual = float(fields[5])
                            except (ValueError, IndexError):
                                qual = None
                        
                        # Apply allele frequency filter
                        af = None
                        if min_af is not None:
                            info_field = fields[7] if len(fields) > 7 else ""
                            af_value = self._parse_af_from_info(info_field)
                            if af_value is None or af_value < min_af:
                                continue
                            af = af_value
                        else:
                            info_field = fields[7] if len(fields) > 7 else ""
                            af = self._parse_af_from_info(info_field)
                        
                        # FIXED: Determine if this is a fixed SNP
                        is_fixed = False
                        if af is not None:
                            if af >= 1.0:
                                # Always substitute perfect frequency variants
                                is_fixed = True
                            elif Config.SNP_ALLELE_FREQUENCY_THRESHOLD is not None:
                                # Use threshold if configured
                                fixed_af_threshold = 1.0 - Config.SNP_ALLELE_FREQUENCY_THRESHOLD
                                if af >= fixed_af_threshold:
                                    # Also check quality if threshold is set
                                    if Config.SNP_QUALITY_THRESHOLD is None or (qual is not None and qual > Config.SNP_QUALITY_THRESHOLD):
                                        is_fixed = True
                        
                        if is_fixed:
                            fixed_count += 1
                        
                        # Take only the first alternate allele for simplicity
                        alt_allele = alt.split(',')[0] if ',' in alt else alt
                        
                        variant = Variant(
                            position=pos,
                            ref=ref,  # FIXED: Use original case
                            alt=alt_allele,  # FIXED: Use original case
                            qual=qual,
                            af=af,
                            is_fixed=is_fixed
                        )
                        
                        if vcf_chrom not in variants:
                            variants[vcf_chrom] = []
                        variants[vcf_chrom].append(variant)
                        total_kept += 1
                        
                    except ValueError:
                        logger.warning(f"Invalid position value: {fields[1]}")
                        continue
            
            logger.info(f"Manual parsing: {total_processed} variants examined, {total_kept} kept after filtering")
            logger.info(f"Fixed SNPs identified: {fixed_count}")
            if min_af is not None:
                logger.info(f"Manual AF filtering applied: >= {min_af}")
            if min_qual is not None:
                logger.info(f"Manual QUAL filtering applied: >= {min_qual}")
            
        except Exception as e:
            error_msg = f"Error in manual VCF parsing: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        
        return variants
    
    def mask_and_substitute_variants(self, sequence: str, variants: list, 
                                flanking_size: int = 0, use_soft_masking: bool = False) -> str:
        """
        Mask or substitute variants in sequence based on their type (fixed vs. variable).
        
        Fixed SNPs (AF=1.0 or above threshold) are substituted with the alternate allele.
        Variable SNPs are masked according to the masking parameters.
        Now supports indels (insertions/deletions) for both substitution and masking.
        
        Args:
            sequence: Input DNA sequence
            variants: List of Variant objects
            flanking_size: Number of bases to mask around each variable variant
            use_soft_masking: Use lowercase letters instead of 'N' for masking
            
        Returns:
            Modified sequence with fixed SNPs substituted and variable SNPs masked
        """
        if not variants:
            logger.debug("No variants to process")
            return sequence
            
        sequence_list = list(sequence)
        masked_positions = set()
        substituted_count = 0
        masked_count = 0
        
        # Sort variants by position for consistent processing
        variants_sorted = sorted(variants, key=lambda v: v.position)
        
        logger.debug(f"Processing {len(variants_sorted)} variants with flanking_size={flanking_size}")
        
        # Track all substitutions for detailed logging
        substitutions_made = []
        
        if Config.SHOW_PROGRESS and len(variants_sorted) > 1000:
            variant_iter = tqdm(variants_sorted, desc="Processing variants")
        else:
            variant_iter = variants_sorted
            
        for variant in variant_iter:
            # Convert 1-based VCF position to 0-based sequence index
            center_idx = variant.position - 1
            
            # Check bounds
            if center_idx < 0 or center_idx >= len(sequence_list):
                logger.warning(f"Variant position {variant.position} is out of sequence bounds")
                continue
            
            # Determine variant type
            ref_len = len(variant.ref)
            alt_len = len(variant.alt)
            is_snp = (ref_len == 1 and alt_len == 1)
            is_deletion = (ref_len > alt_len)
            is_insertion = (ref_len < alt_len)
            
            if variant.is_fixed:
                # For fixed variants, try substitution
                if is_snp:
                    # Handle single nucleotide variants (SNPs)
                    if center_idx < len(sequence_list):
                        original_base = sequence_list[center_idx]
                        
                        # Capture pre-substitution context
                        context_start = max(0, center_idx - 10)
                        context_end = min(len(sequence_list), center_idx + 11)
                        pre_context = "".join(sequence_list[context_start:context_end])
                        
                        # Verify the reference matches (case-insensitive)
                        if original_base.upper() == variant.ref.upper():
                            # Make the substitution
                            sequence_list[center_idx] = variant.alt.upper()
                            substituted_count += 1
                            
                            # Capture post-substitution context
                            post_context = "".join(sequence_list[context_start:context_end])
                            
                            # Store substitution details for logging
                            substitutions_made.append({
                                'position': variant.position,
                                'ref': variant.ref,
                                'alt': variant.alt,
                                'af': variant.af,
                                'qual': variant.qual,
                                'pre_context': pre_context,
                                'post_context': post_context,
                                'context_start': context_start + 1,  # Convert back to 1-based
                                'center_offset': center_idx - context_start,
                                'variant_type': 'SNP'
                            })
                            
                            logger.debug(f"Substituted {variant.ref} -> {variant.alt} at position {variant.position} (AF={variant.af})")
                        else:
                            logger.warning(f"Reference mismatch at position {variant.position}: "
                                        f"expected {variant.ref}, found {original_base}. Masking instead.")
                            # DIAGNOSTIC: Add detailed mismatch info
                            logger.debug(f"MISMATCH DIAGNOSTIC: VCF pos {variant.position} (0-based: {center_idx})")
                            logger.debug(f"  Sequence context: ...{sequence[max(0, center_idx-5):center_idx+6]}...")
                            logger.debug(f"  Expected: {variant.ref} at position {center_idx}")
                            logger.debug(f"  Found: {original_base} at position {center_idx}")
                            
                            # Fall back to masking if reference doesn't match
                            self._apply_masking_at_position(sequence_list, center_idx, flanking_size, 
                                                        use_soft_masking, masked_positions)
                            masked_count += 1
                else:
                    # Handle indels (insertions/deletions) for fixed variants
                    # Capture pre-substitution context for indels
                    context_start = max(0, center_idx - 10)
                    context_end = min(len(sequence_list), center_idx + 10 + ref_len)
                    pre_context = "".join(sequence_list[context_start:context_end])
                    
                    if self._handle_indel_substitution(sequence_list, variant, center_idx):
                        substituted_count += 1
                        
                        # Capture post-substitution context (sequence_list has changed!)
                        # Adjust context_end based on the change in sequence length
                        length_change = alt_len - ref_len
                        post_context_end = min(len(sequence_list), context_end + length_change)
                        post_context = "".join(sequence_list[context_start:post_context_end])
                        
                        substitutions_made.append({
                            'position': variant.position,
                            'ref': variant.ref,
                            'alt': variant.alt,
                            'af': variant.af,
                            'qual': variant.qual,
                            'pre_context': pre_context,
                            'post_context': post_context,
                            'context_start': context_start + 1,
                            'center_offset': center_idx - context_start,
                            'variant_type': 'DELETION' if is_deletion else 'INSERTION',
                            'ref_len': ref_len,
                            'alt_len': alt_len
                        })
                        
                        logger.debug(f"Applied indel substitution {variant.ref} -> {variant.alt} at position {variant.position} (AF={variant.af})")
                    else:
                        logger.warning(f"Could not apply indel substitution at position {variant.position}: "
                                    f"{variant.ref} -> {variant.alt}. Masking instead.")
                        # Fall back to masking the affected region
                        self._apply_indel_masking(sequence_list, variant, center_idx, flanking_size, 
                                                use_soft_masking, masked_positions)
                        masked_count += 1
            else:
                # For variable variants, apply masking
                if is_snp:
                    # Handle SNP masking
                    self._apply_masking_at_position(sequence_list, center_idx, flanking_size, 
                                                use_soft_masking, masked_positions)
                    masked_count += 1
                else:
                    # Handle indel masking
                    self._apply_indel_masking(sequence_list, variant, center_idx, flanking_size, 
                                            use_soft_masking, masked_positions)
                    masked_count += 1
        
        modified_sequence = "".join(sequence_list)
        
        # Log detailed substitution information (permanent feature)
        if substitutions_made:
            logger.debug(f"\n=== DETAILED SUBSTITUTION REPORT ===")
            for i, sub in enumerate(substitutions_made[:5]):  # Show first 5 substitutions
                logger.debug(f"Substitution {i+1} ({sub['variant_type']}):")
                logger.debug(f"  Position: {sub['position']} (AF={sub['af']:.3f}, QUAL={sub['qual']})")
                logger.debug(f"  Change: {sub['ref']} -> {sub['alt']}")
                
                if sub['variant_type'] == 'SNP':
                    logger.debug(f"  Context positions: {sub['context_start']}-{sub['context_start'] + len(sub['pre_context']) - 1}")
                    
                    # Mark the substitution position with brackets
                    pre_marked = list(sub['pre_context'])
                    post_marked = list(sub['post_context'])
                    center_pos = sub['center_offset']
                    
                    if 0 <= center_pos < len(pre_marked):
                        pre_marked[center_pos] = f"[{pre_marked[center_pos]}]"
                    if 0 <= center_pos < len(post_marked):
                        post_marked[center_pos] = f"[{post_marked[center_pos]}]"
                    
                    logger.debug(f"  Before: {''.join(pre_marked)}")
                    logger.debug(f"  After:  {''.join(post_marked)}")
                else:
                    # FIXED: Better indel bracket display
                    logger.debug(f"  Context positions: {sub['context_start']}-{sub['context_start'] + len(sub['pre_context']) - 1}")
                    
                    # For indels, mark the exact affected region
                    pre_marked = list(sub['pre_context'])
                    post_marked = list(sub['post_context'])
                    center_pos = sub['center_offset']
                    ref_len = sub['ref_len']
                    alt_len = sub['alt_len']
                    
                    # FIXED: Mark the reference sequence in the before context
                    if center_pos >= 0 and center_pos + ref_len <= len(pre_marked):
                        if ref_len == 1:
                            # Single character: [X]
                            pre_marked[center_pos] = f"[{pre_marked[center_pos]}]"
                        else:
                            # Multiple characters: [XXX]
                            for j in range(ref_len):
                                if j == 0:
                                    pre_marked[center_pos + j] = f"[{pre_marked[center_pos + j]}"
                                elif j == ref_len - 1:
                                    pre_marked[center_pos + j] = f"{pre_marked[center_pos + j]}]"
                                # middle characters unchanged
                    
                    # FIXED: Mark the alternate sequence in the after context  
                    if center_pos >= 0 and center_pos + alt_len <= len(post_marked):
                        if alt_len == 1:
                            # Single character: [X]
                            post_marked[center_pos] = f"[{post_marked[center_pos]}]"
                        else:
                            # Multiple characters: [XXX]
                            for j in range(alt_len):
                                if j == 0:
                                    post_marked[center_pos + j] = f"[{post_marked[center_pos + j]}"
                                elif j == alt_len - 1:
                                    post_marked[center_pos + j] = f"{post_marked[center_pos + j]}]"
                                # middle characters unchanged
                    
                    logger.debug(f"  Before: {''.join(pre_marked)}")
                    logger.debug(f"  After:  {''.join(post_marked)}")
                logger.debug("")
            
            if len(substitutions_made) > 5:
                logger.debug(f"  ... and {len(substitutions_made) - 5} more substitutions")
            logger.debug(f"=== END SUBSTITUTION REPORT ===\n")
        
        # Calculate statistics
        total_modified = substituted_count + len(masked_positions)
        modification_ratio = total_modified / len(sequence) if len(sequence) > 0 else 0
        
        mask_type = "soft" if use_soft_masking else "hard"
        logger.debug(f"Sequence modification summary:")
        logger.debug(f"  Fixed SNPs substituted: {substituted_count}")
        logger.debug(f"  Variable SNPs masked ({mask_type}): {len(masked_positions)}")
        logger.debug(f"  Total modifications: {total_modified} ({modification_ratio:.2%} of sequence)")
        
        if flanking_size > 0:
            logger.debug(f"  Flanking region size: {flanking_size} bases around each variable variant")
        
        return modified_sequence

    def _handle_indel_substitution(self, sequence_list: list, variant, center_idx: int) -> bool:
        """
        FIXED: Handle indel substitution for fixed variants.
        
        Key fixes:
        - Proper bounds checking for multi-base references
        - Case-insensitive sequence matching
        - Better error diagnostics
        
        Args:
            sequence_list: Mutable list of sequence characters
            variant: Variant object with ref, alt, position
            center_idx: 0-based center index
            
        Returns:
            True if substitution was successful, False otherwise
        """
        ref_len = len(variant.ref)
        alt_len = len(variant.alt)
        
        # FIXED: Proper bounds checking for the full reference length
        if center_idx < 0 or center_idx + ref_len > len(sequence_list):
            logger.warning(f"Not enough sequence to match reference {variant.ref} at position {variant.position}")
            logger.debug(f"  center_idx: {center_idx}, ref_len: {ref_len}, sequence_length: {len(sequence_list)}")
            logger.debug(f"  Required end position: {center_idx + ref_len}, available: {len(sequence_list)}")
            return False
        
        # Extract the reference sequence from our sequence
        ref_sequence = "".join(sequence_list[center_idx:center_idx + ref_len])
        
        # FIXED: Case-insensitive reference matching
        if ref_sequence.upper() != variant.ref.upper():
            logger.warning(f"Reference sequence mismatch at position {variant.position}: "
                        f"expected '{variant.ref}', found '{ref_sequence}'")
            
            # Enhanced diagnostic info
            logger.debug(f"INDEL MISMATCH DIAGNOSTIC: VCF pos {variant.position} (0-based: {center_idx})")
            logger.debug(f"  Expected ref sequence: '{variant.ref}' (length {ref_len})")
            logger.debug(f"  Found sequence: '{ref_sequence}' (length {len(ref_sequence)})")
            logger.debug(f"  Case-insensitive comparison: '{ref_sequence.upper()}' vs '{variant.ref.upper()}'")
            context_start = max(0, center_idx - 10)
            context_end = min(len(sequence_list), center_idx + 10 + ref_len)
            context = "".join(sequence_list[context_start:context_end])
            logger.debug(f"  Sequence context: {context}")
            logger.debug(f"  Context start: {context_start + 1}")
            
            return False
        
        # Perform the substitution
        # Remove the reference sequence and insert the alternate
        del sequence_list[center_idx:center_idx + ref_len]
        
        # Insert the alternate sequence
        for i, base in enumerate(variant.alt.upper()):
            sequence_list.insert(center_idx + i, base)
        
        logger.debug(f"Successfully applied indel: {variant.ref} -> {variant.alt} at position {variant.position}")
        return True

    def _apply_indel_masking(self, sequence_list: list, variant, center_idx: int, 
                            flanking_size: int, use_soft_masking: bool, masked_positions: set):
        """
        Apply masking for indel variants.
        
        Args:
            sequence_list: Mutable list of sequence characters
            variant: Variant object
            center_idx: 0-based center index
            flanking_size: Number of flanking bases to mask
            use_soft_masking: Use lowercase instead of 'N'
            masked_positions: Set to track masked positions
        """
        ref_len = len(variant.ref)
        
        # We want to mask the reference sequence + flanking
        start_idx = max(0, center_idx - flanking_size)
        end_idx = min(len(sequence_list), center_idx + ref_len + flanking_size)
        
        # Apply masking to the range
        for idx in range(start_idx, end_idx):
            original_base = sequence_list[idx]
            
            if use_soft_masking:
                sequence_list[idx] = original_base.lower()
            else:
                sequence_list[idx] = 'N'
            
            if original_base.upper() in 'ATCG':
                masked_positions.add(idx)
        
        # FIXED: Better logging with 1-based coordinates
        logger.debug(f"Masked indel region {start_idx + 1}-{end_idx} for variant {variant.ref} -> {variant.alt} at position {variant.position}")
        logger.debug(f"  Reference length: {ref_len}, flanking: {flanking_size}, total masked: {end_idx - start_idx} bases")

    def _apply_masking_at_position(self, sequence_list: list, center_idx: int, 
                                flanking_size: int, use_soft_masking: bool, 
                                masked_positions: set):
        """Apply masking at a specific position with flanking regions."""
        # Add flanking positions
        for offset in range(-flanking_size, flanking_size + 1):
            mask_idx = center_idx + offset
            if 0 <= mask_idx < len(sequence_list):
                original_base = sequence_list[mask_idx]
                
                if use_soft_masking:
                    # Convert to lowercase for soft masking
                    sequence_list[mask_idx] = original_base.lower()
                else:
                    # Use 'N' for hard masking
                    sequence_list[mask_idx] = 'N'
                
                if original_base.upper() in 'ATCG':  # Only count actual nucleotides as masked
                    masked_positions.add(mask_idx)
    
    def mask_sequences_for_primer_design(self, sequences: dict, variants: dict, 
                                        flanking_size: int = 0, use_soft_masking: bool = False,
                                        min_af: float = None, min_qual: float = None) -> dict:
        """
        Mask and substitute variants in sequences for primer design.
        
        This method handles both masking of variable SNPs and substitution of fixed SNPs.
        Fixed SNPs are determined by AF=1.0 or above threshold with quality checks.
        
        Args:
            sequences: Dictionary of sequences
            variants: Dictionary mapping chromosomes to variant data
            flanking_size: Number of bases to mask around each variable variant
            use_soft_masking: Use lowercase letters instead of 'N'
            min_af: Minimum allele frequency threshold for filtering
            min_qual: Minimum QUAL score threshold for filtering
            
        Returns:
            Dictionary of processed sequences with fixed SNPs substituted and variable SNPs masked
        """
        mask_type = "soft" if use_soft_masking else "hard"
        logger.info(f"\nProcessing variants from VCF file...")
        
        if flanking_size > 0:
            logger.debug(f"Using flanking region size: {flanking_size} bases")
        if min_af is not None:
            logger.debug(f"AF filter: variants with AF >= {min_af}")
        if min_qual is not None:
            logger.debug(f"QUAL filter: variants with QUAL >= {min_qual}")
        
        # Log fixed SNP criteria
        logger.debug(f"Fixed SNP criteria: AF=1.0 always substituted")
        if Config.SNP_ALLELE_FREQUENCY_THRESHOLD is not None:
            fixed_af_threshold = 1.0 - Config.SNP_ALLELE_FREQUENCY_THRESHOLD
            if Config.SNP_QUALITY_THRESHOLD is not None:
                logger.debug(f"Additional fixed SNP criteria: AF >= {fixed_af_threshold:.3f} AND QUAL > {Config.SNP_QUALITY_THRESHOLD}")
            else:
                logger.debug(f"Additional fixed SNP criteria: AF >= {fixed_af_threshold:.3f}")
        
        processed_sequences = {}
        
        if Config.SHOW_PROGRESS:
            sequence_iter = tqdm(sequences.items(), desc=f"Processing sequences:")
        else:
            sequence_iter = sequences.items()
        
        total_substitutions = 0
        total_maskings = 0
        
        for seq_id, sequence in sequence_iter:
            # Try to find matching variants for this sequence
            if seq_id in variants:
                variant_data = variants[seq_id]
                
                # Use the new method for variants
                processed_sequence = self.mask_and_substitute_variants(
                    sequence, 
                    variant_data, 
                    flanking_size=flanking_size,
                    use_soft_masking=use_soft_masking
                )
                
                # Count fixed substitutions for logging
                fixed_variants = [v for v in variant_data if v.is_fixed]
                total_substitutions += len(fixed_variants)
                total_maskings += len(variant_data) - len(fixed_variants)
                
                processed_sequences[seq_id] = processed_sequence
            else:
                # No variants found for this sequence, use original
                processed_sequences[seq_id] = sequence
                logger.debug(f"No variants found for sequence {seq_id}, using original")
        
        logger.debug(f"Completed variant processing for {len(processed_sequences)} sequences:")
        logger.info(f"Fixed SNPs substituted: {total_substitutions}")
        logger.info(f"SNPs masked: {total_maskings}")
        
        return processed_sequences
    
    def process_vcf_with_chromosome_mapping(self, vcf_file: str, fasta_file: str, 
                                           chromosome: str = None, min_af: float = None, 
                                           min_qual: float = None) -> dict:
        """
        Process VCF file with automatic chromosome name mapping to match FASTA file.
        
        Args:
            vcf_file: Path to VCF file
            fasta_file: Path to FASTA file for chromosome name reference
            chromosome: Specific chromosome to filter
            min_af: Minimum allele frequency threshold
            min_qual: Minimum QUAL score threshold
            
        Returns:
            Dictionary mapping FASTA chromosome names to variant data
            
        Raises:
            SequenceProcessingError: If VCF processing fails
        """
        logger.debug("Processing VCF file with automatic chromosome mapping")
        
        try:
            # Import here to avoid circular imports
            from ..utils import ChromosomeMapper
            
            # Initialize chromosome mapper
            mapper = ChromosomeMapper()
            
            # Check compatibility and get analysis
            try:
                logger.debug("Checking chromosome compatibility")
                analysis = mapper.check_chromosome_compatibility(vcf_file, fasta_file)
                
                if not analysis['needs_mapping']:
                    # Files are compatible - use direct processing
                    logger.debug("Chromosome names are compatible - processing directly")
                    return self.get_variants(vcf_file, chromosome=chromosome, 
                                           min_af=min_af, min_qual=min_qual)
                
                # Files need mapping - get the suggested mapping
                mapping = analysis.get('suggested_mapping')
                if not mapping:
                    logger.warning("Could not generate automatic chromosome mapping")
                    logger.debug("Falling back to direct VCF processing")
                    return self.get_variants(vcf_file, chromosome=chromosome, 
                                           min_af=min_af, min_qual=min_qual)
                
                # Mapping was generated - files need chromosome name conversion
                logger.debug("Chromosome name mismatch detected - applying dynamic mapping")
                logger.info(f"Generated chromosome mapping for {len(mapping)} chromosomes")
                for vcf_chrom, fasta_chrom in mapping.items():
                    logger.debug(f"  '{vcf_chrom}' -> '{fasta_chrom}'")
                
                # Handle chromosome filtering with mapping
                target_vcf_chromosome = None
                if chromosome:
                    # User specified a FASTA chromosome name, find corresponding VCF name
                    reverse_mapping = {v: k for k, v in mapping.items()}
                    if chromosome in reverse_mapping:
                        target_vcf_chromosome = reverse_mapping[chromosome]
                        logger.debug(f"Filtering for FASTA chromosome '{chromosome}' "
                            f"(VCF chromosome '{target_vcf_chromosome}')")
                    elif chromosome in mapping:
                        # User provided VCF chromosome name
                        target_vcf_chromosome = chromosome
                        logger.debug(f"Filtering for VCF chromosome '{chromosome}'")
                    else:
                        logger.warning(f"Chromosome '{chromosome}' not found in mapping. "
                                    f"Available: {list(mapping.values())}")
                        return {}
                
                # Extract variants using original VCF chromosome names
                logger.debug("Extracting variants from VCF with original chromosome names")
                original_variants = self.get_variants(
                    vcf_file, 
                    chromosome=target_vcf_chromosome, 
                    min_af=min_af, 
                    min_qual=min_qual
                )
                
                # Apply dynamic mapping to convert chromosome names
                logger.debug("Applying chromosome name mapping to variants")
                mapped_variants = {}
                total_variants = 0
                mapped_chromosomes = 0
                
                for vcf_chrom, variant_data in original_variants.items():
                    if vcf_chrom in mapping:
                        fasta_chrom = mapping[vcf_chrom]
                        mapped_variants[fasta_chrom] = variant_data
                        total_variants += len(variant_data)
                        mapped_chromosomes += 1
                        logger.debug(f"Mapped {len(variant_data)} variants from "
                            f"'{vcf_chrom}' to '{fasta_chrom}'")
                    else:
                        logger.warning(f"No mapping found for VCF chromosome '{vcf_chrom}' "
                                    f"- skipping {len(variant_data)} variants")
                
                logger.info(f"Successfully mapped {total_variants} variants across "
                        f"{mapped_chromosomes} chromosomes")
                
                # Validate that we have variants for requested chromosome
                if chromosome and chromosome not in mapped_variants:
                    logger.warning(f"No variants found for requested chromosome '{chromosome}'")
                
                return mapped_variants
                    
            except Exception as mapping_error:
                logger.warning(f"Chromosome mapping failed: {mapping_error}")
                logger.debug("Falling back to direct VCF processing")
                
                # Fall back to direct processing (might fail if chromosome names don't match)
                return self.get_variants(vcf_file, chromosome=chromosome, 
                                       min_af=min_af, min_qual=min_qual)
                
        except Exception as e:
            error_msg = f"Error in VCF processing with chromosome mapping: {e}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
    
    def _parse_af_from_info(self, info_field: str) -> float:
        """
        Parse allele frequency from INFO field.
        
        Args:
            info_field: INFO field from VCF line
            
        Returns:
            Allele frequency if found, None otherwise
        """
        try:
            # Look for AF= in the INFO field
            for info_item in info_field.split(';'):
                if info_item.startswith('AF='):
                    af_str = info_item[3:]  # Remove 'AF=' prefix
                    # Handle multiple values (take the first one)
                    if ',' in af_str:
                        af_str = af_str.split(',')[0]
                    return float(af_str)
        except (ValueError, AttributeError):
            pass
        
        return None
    
    def _parse_af_from_string(self, af_str: str) -> float:
        """
        Parse allele frequency from a string that might be from INFO/AF field.
        
        Args:
            af_str: String containing AF value
            
        Returns:
            Allele frequency if found, None otherwise
        """
        try:
            if af_str and af_str != "." and af_str != "":
                # Handle comma-separated values (take first)
                if ',' in af_str:
                    af_str = af_str.split(',')[0]
                return float(af_str)
        except (ValueError, AttributeError):
            pass
        
        return None