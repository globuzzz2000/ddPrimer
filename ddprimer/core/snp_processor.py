#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP masking processor module for ddPrimer pipeline.

Contains functionality for:
1. VCF-based sequence masking and substitution for primer design
2. Fixed variant substitution (AF=1.0) vs variable variant masking (AF<1.0)
3. Processing of prepared VCF files (bgzipped, normalized, indexed)
4. Quality and allele frequency filtering with bcftools integration

This module works with prepared VCF files from FilePreparator to provide
streamlined variant processing for the ddPrimer pipeline. File preparation
(compression, normalization, chromosome mapping) is handled upstream.
"""

import os
import subprocess
import logging
from typing import Dict, List, Optional
from ..config import Config, FileError, ExternalToolError, SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)


class SNPMaskingProcessor:
    """
    Handles VCF-based sequence masking and substitution for primer design.
    
    This class processes prepared VCF variants to either substitute fixed variants (AF=1.0)
    into sequences or mask variable variants (AF<1.0) with 'N' or soft masking.
    Works with VCF files that have been prepared by FilePreparator (bgzipped,
    normalized, indexed, chromosome names harmonized).
    
    Attributes:
        reference_file: Path to reference FASTA file
        
    Example:
        >>> processor = SNPMaskingProcessor("genome.fasta")
        >>> processed_seq = processor.process_sequence_with_vcf(
        ...     sequence="ATCGATCG", 
        ...     vcf_path="prepared_variants.vcf.gz",
        ...     chromosome="chr1"
        ... )
    """
    
    def __init__(self, reference_file: str):
        """
        Initialize SNP masking processor.
        
        Args:
            reference_file: Path to reference FASTA file
            
        Raises:
            FileError: If reference file cannot be accessed
            ExternalToolError: If bcftools is not available
        """
        logger.debug("=== SNP PROCESSOR INITIALIZATION DEBUG ===")
        logger.debug(f"Initializing SNPMaskingProcessor with reference: {reference_file}")
        
        if not os.path.exists(reference_file):
            error_msg = f"Reference FASTA file not found: {reference_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        self.reference_file = reference_file
        
        # Validate bcftools availability
        try:
            Config.validate_vcf_dependencies()
            logger.debug("VCF dependencies validated successfully")
        except Exception as e:
            error_msg = f"VCF processing dependencies not available: {str(e)}"
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="bcftools") from e
        
        logger.debug("=== END SNP PROCESSOR INITIALIZATION DEBUG ===")
    
    def process_sequence_with_vcf(self, sequence: str, vcf_path: str, 
                                chromosome: str, **kwargs) -> str:
        """
        Process sequence with prepared VCF variants.
        
        Main entry point for sequence processing. Handles both fixed variant
        substitution (AF=1.0) and variable variant masking (AF<1.0) based on
        configuration settings. Assumes VCF is already prepared (normalized,
        bgzipped, indexed, chromosome names harmonized).
        
        Args:
            sequence: DNA sequence to process
            vcf_path: Path to prepared VCF file
            chromosome: Chromosome/sequence identifier (should match VCF)
            **kwargs: Additional processing settings from Config
            
        Returns:
            Processed sequence with variants applied
            
        Raises:
            SequenceProcessingError: If sequence processing fails
        """
        logger.debug("=== SEQUENCE PROCESSING DEBUG ===")
        logger.debug(f"Processing sequence for chromosome: {chromosome}")
        logger.debug(f"Original sequence length: {len(sequence)}")
        
        if not sequence or not isinstance(sequence, str):
            logger.warning("Empty or invalid sequence provided")
            return sequence
        
        try:
            # Get processing settings
            settings = {
                'snp_allele_frequency_threshold': kwargs.get('snp_allele_frequency_threshold', Config.VCF_ALLELE_FREQUENCY_THRESHOLD),
                'snp_quality_threshold': kwargs.get('snp_quality_threshold', Config.VCF_QUALITY_THRESHOLD),
                'snp_flanking_mask_size': kwargs.get('snp_flanking_mask_size', Config.VCF_FLANKING_MASK_SIZE),
                'snp_use_soft_masking': kwargs.get('snp_use_soft_masking', Config.VCF_USE_SOFT_MASKING),
            }
            
            logger.debug(f"Processing settings: {settings}")
            
            # Parse variants from prepared VCF file
            variants = self._parse_variants_from_prepared_vcf(vcf_path, chromosome, **settings)
            logger.debug(f"Found {len(variants)} applicable variants for chromosome {chromosome}")
            
            if not variants:
                logger.debug("No variants to process - returning original sequence")
                return sequence
            
            # Apply variants to sequence
            processed_sequence = self._apply_variants_to_sequence(sequence, variants, **settings)
            
            # Log processing statistics
            original_len = len(sequence)
            processed_len = len(processed_sequence)
            n_count = processed_sequence.count('N')
            soft_count = sum(1 for c in processed_sequence if c.islower())
            
            logger.debug(f"Processing complete: {original_len} -> {processed_len} bp")
            logger.debug(f"Hard masked (N): {n_count}, Soft masked (lowercase): {soft_count}")
            logger.debug("=== END SEQUENCE PROCESSING DEBUG ===")
            
            return processed_sequence
            
        except Exception as e:
            error_msg = f"Error processing sequence for chromosome {chromosome}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END SEQUENCE PROCESSING DEBUG ===")
            raise SequenceProcessingError(error_msg) from e
    
    def _parse_variants_from_prepared_vcf(self, vcf_path: str, chromosome: str, **kwargs) -> List[Dict]:
        """
        Parse variants from prepared VCF file for specific chromosome.
        
        Uses bcftools to extract variants from prepared VCF file (already
        bgzipped, normalized, indexed). No additional preparation needed.
        
        Args:
            vcf_path: Path to prepared VCF file
            chromosome: Chromosome to extract variants for
            **kwargs: Processing settings for filtering
            
        Returns:
            List of variant dictionaries
            
        Raises:
            ExternalToolError: If bcftools execution fails
        """
        logger.debug(f"=== PARSING VARIANTS FROM PREPARED VCF DEBUG ===")
        logger.debug(f"Extracting variants for chromosome {chromosome} from {vcf_path}")
        
        try:
            # Use bcftools query to extract variants for specific chromosome
            cmd = [
                'bcftools', 'query', 
                '-r', chromosome,  # Restrict to chromosome
                '-f', '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%AF\\n',
                vcf_path
            ]
            
            logger.debug(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                error_msg = f"bcftools query failed: {result.stderr}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="bcftools")
            
            # Parse variants from output
            variants = []
            for line_num, line in enumerate(result.stdout.strip().split('\n'), 1):
                if line.strip():
                    try:
                        variant = self._parse_variant_line(line.strip())
                        if self._should_process_variant(variant, **kwargs):
                            variants.append(variant)
                    except Exception as e:
                        logger.debug(f"Error parsing variant line {line_num}: {str(e)}")
                        continue
            
            logger.debug(f"Parsed {len(variants)} applicable variants for chromosome {chromosome}")
            logger.debug("=== END PARSING VARIANTS FROM PREPARED VCF DEBUG ===")
            return variants
            
        except ExternalToolError:
            # Re-raise ExternalToolError without modification
            raise
        except Exception as e:
            error_msg = f"Error parsing variants from prepared VCF: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END PARSING VARIANTS FROM PREPARED VCF DEBUG ===")
            raise ExternalToolError(error_msg, tool_name="bcftools") from e
    
    def _parse_variant_line(self, line: str) -> Dict:
        """
        Parse single variant line from bcftools query output.
        
        Args:
            line: Tab-separated variant line (CHROM POS REF ALT QUAL AF)
            
        Returns:
            Dictionary containing variant information
            
        Raises:
            ValueError: If line format is invalid
        """
        try:
            parts = line.split('\t')
            if len(parts) < 6:
                raise ValueError(f"Invalid variant line format: expected 6 fields, got {len(parts)}")
            
            # Parse basic fields
            chrom, pos, ref, alt, qual, af = parts[:6]
            
            # Convert types
            pos = int(pos)
            qual = float(qual) if qual != '.' else None
            
            # Parse AF field (can be comma-separated for multiallelic)
            if af == '.' or af == '':
                af_value = None
            else:
                # Take first AF value for multiallelic sites
                af_value = float(af.split(',')[0])
            
            variant = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'af': af_value
            }
            
            logger.debug(f"Parsed variant: {variant}")
            return variant
            
        except (ValueError, IndexError) as e:
            error_msg = f"Error parsing variant line '{line}': {str(e)}"
            logger.debug(error_msg)
            raise ValueError(error_msg) from e
    
    def _should_process_variant(self, variant: Dict, **kwargs) -> bool:
        """
        Determine if variant should be processed based on quality and AF thresholds.
        
        Args:
            variant: Variant dictionary
            **kwargs: Processing settings
            
        Returns:
            True if variant should be processed
        """
        # Check quality threshold
        qual_threshold = kwargs.get('snp_quality_threshold')
        if qual_threshold is not None and variant.get('qual') is not None:
            if variant['qual'] < qual_threshold:
                logger.debug(f"Skipping variant at {variant['pos']} due to low quality: {variant['qual']}")
                return False
        
        # Always process if AF is not available
        if variant.get('af') is None:
            logger.debug(f"Processing variant at {variant['pos']} (no AF information)")
            return True
        
        # Check AF threshold
        af_threshold = kwargs.get('snp_allele_frequency_threshold')
        if af_threshold is not None:
            if variant['af'] < af_threshold:
                logger.debug(f"Skipping rare variant at {variant['pos']} (AF={variant['af']:.3f})")
                return False
        
        logger.debug(f"Processing variant at {variant['pos']} (AF={variant['af']:.3f}, QUAL={variant.get('qual', 'N/A')})")
        return True
    
    def _classify_variant(self, variant: Dict, **kwargs) -> str:
        """
        Classify variant action based on allele frequency.
        
        Args:
            variant: Variant dictionary
            **kwargs: Processing settings
            
        Returns:
            'substitute' for fixed variants (AF=1.0), 'mask' for variable variants
        """
        af = variant.get('af')
        
        # Fixed variants (AF = 1.0) should be substituted
        if af is not None and abs(af - 1.0) < 0.001:  # Account for floating point precision
            logger.debug(f"Variant at {variant['pos']} classified as 'substitute' (AF={af:.3f})")
            return 'substitute'
        
        # All other variants should be masked
        logger.debug(f"Variant at {variant['pos']} classified as 'mask' (AF={af})")
        return 'mask'
    
    def _apply_variants_to_sequence(self, sequence: str, variants: List[Dict], **kwargs) -> str:
        """
        Apply all variants to sequence in coordinate order.
        
        Processes variants in reverse order (high to low position) to maintain
        coordinate integrity during insertions and deletions.
        
        Args:
            sequence: Original DNA sequence
            variants: List of variant dictionaries
            **kwargs: Processing settings
            
        Returns:
            Modified sequence with variants applied
        """
        if not variants:
            logger.debug("No variants to apply")
            return sequence
        
        # Sort variants by position in descending order
        # Process from end to beginning to maintain coordinate integrity
        sorted_variants = sorted(variants, key=lambda v: v['pos'], reverse=True)
        
        logger.debug(f"Applying {len(sorted_variants)} variants to sequence (length: {len(sequence)})")
        
        modified_sequence = sequence
        applied_count = 0
        
        for variant in sorted_variants:
            try:
                action = self._classify_variant(variant, **kwargs)
                
                # Apply variant based on classification
                if action == 'substitute':
                    modified_sequence = self._apply_substitution(modified_sequence, variant, **kwargs)
                elif action == 'mask':
                    modified_sequence = self._apply_masking(modified_sequence, variant, **kwargs)
                
                applied_count += 1
                
            except Exception as e:
                logger.warning(f"Error applying variant at position {variant['pos']}: {str(e)}")
                logger.debug(f"Variant details: {variant}")
                continue
        
        logger.debug(f"Successfully applied {applied_count}/{len(sorted_variants)} variants")
        return modified_sequence
    
    def _apply_substitution(self, sequence: str, variant: Dict, **kwargs) -> str:
        """
        Apply variant substitution to sequence (for fixed variants AF=1.0).
        
        Args:
            sequence: DNA sequence
            variant: Variant dictionary
            **kwargs: Processing settings
            
        Returns:
            Sequence with variant substituted
        """
        pos = variant['pos'] - 1  # Convert to 0-based indexing
        ref = variant['ref']
        alt = variant['alt']
        
        # Validate position and reference
        if pos < 0 or pos >= len(sequence):
            logger.debug(f"Position {variant['pos']} out of sequence bounds")
            return sequence
        
        if pos + len(ref) > len(sequence):
            logger.debug(f"Reference allele extends beyond sequence end at position {variant['pos']}")
            return sequence
        
        # Check reference match
        sequence_ref = sequence[pos:pos + len(ref)]
        if sequence_ref.upper() != ref.upper():
            logger.debug(f"Reference mismatch at {variant['pos']}: expected {ref}, found {sequence_ref}")
            return sequence
        
        # Apply substitution
        new_sequence = sequence[:pos] + alt + sequence[pos + len(ref):]
        
        logger.debug(f"Substituted {ref}->{alt} at position {variant['pos']}")
        return new_sequence
    
    def _apply_masking(self, sequence: str, variant: Dict, **kwargs) -> str:
        """
        Apply variant masking to sequence (for variable variants AF<1.0).
        
        Args:
            sequence: DNA sequence
            variant: Variant dictionary
            **kwargs: Processing settings
            
        Returns:
            Sequence with variant region masked
        """
        pos = variant['pos'] - 1  # Convert to 0-based indexing
        ref = variant['ref']
        
        # Validate position
        if pos < 0 or pos >= len(sequence):
            logger.debug(f"Position {variant['pos']} out of sequence bounds")
            return sequence
        
        if pos + len(ref) > len(sequence):
            logger.debug(f"Reference allele extends beyond sequence end at position {variant['pos']}")
            return sequence
        
        # Determine masking strategy
        use_soft_masking = kwargs.get('snp_use_soft_masking', False)
        flanking_mask_size = kwargs.get('snp_flanking_mask_size', 0)
        
        # Calculate masking region
        start_pos = max(0, pos - flanking_mask_size)
        end_pos = min(len(sequence), pos + len(ref) + flanking_mask_size)
        
        # Apply masking
        if use_soft_masking:
            # Soft masking - convert to lowercase
            masked_region = sequence[start_pos:end_pos].lower()
            logger.debug(f"Soft masked region {start_pos+1}-{end_pos} at variant {variant['pos']}")
        else:
            # Hard masking - replace with 'N'
            masked_region = 'N' * (end_pos - start_pos)
            logger.debug(f"Hard masked region {start_pos+1}-{end_pos} at variant {variant['pos']}")
        
        new_sequence = sequence[:start_pos] + masked_region + sequence[end_pos:]
        return new_sequence