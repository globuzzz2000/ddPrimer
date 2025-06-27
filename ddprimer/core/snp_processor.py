#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP masking processor module for ddPrimer pipeline.

Handles VCF-based sequence masking and substitution for primer design
including fixed variant substitution (AF=1.0) vs variable variant masking (AF<1.0),
processing of prepared VCF files, and quality/allele frequency filtering.

Contains functionality for:
1. VCF-based sequence masking and substitution for primer design
2. Fixed variant substitution (AF=1.0) vs variable variant masking (AF<1.0)
3. Processing of prepared VCF files (bgzipped, normalized, indexed)
4. Quality and allele frequency filtering with bcftools integration
5. Parallel sequence processing for improved performance

COORDINATE SYSTEM:
- INPUT: VCF coordinates (1-based, as per VCF specification)
- INTERNAL: Converts to 0-based for sequence operations
- OUTPUT: Modified sequences (coordinate-independent strings)
- CONVERSION POINT: _apply_substitution() and _apply_masking() methods
- All sequence indexing uses 0-based Python string operations

COORDINATE CONVERSION: VCF positions are converted from 1-based to 0-based
via `pos = variant['pos'] - 1` before any sequence operations.
"""

import os
import subprocess
import logging
from typing import Dict, List
from concurrent.futures import ProcessPoolExecutor, as_completed

# Import package modules
from ..config import Config, FileError, ExternalToolError, SequenceProcessingError, PipelineError, DebugLogLimiter

# Set up module logger
logger = logging.getLogger(__name__)


class SNPMaskingProcessor:
    """
    Handles VCF-based sequence masking and substitution for primer design.
    
    This class processes prepared VCF variants to either substitute fixed variants (AF=1.0)
    into sequences or mask variable variants (AF<1.0) with 'N' or soft masking.
    Works with VCF files that have been prepared by FilePreparator.
    
    COORDINATE HANDLING:
    - Input: VCF coordinates (1-based, per VCF specification)
    - Processing: Automatic conversion to 0-based for sequence operations
    - Output: Modified sequences (coordinate-independent)
    
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
    
    #############################################################################
    #                           Workflow Wrappers
    #############################################################################
    
    @classmethod
    def process_sequences_with_vcf_batch(cls, sequences: Dict[str, str], 
                                        vcf_file: str, reference_file: str) -> Dict[str, str]:
        """
        Process multiple sequences using VCF normalization approach with parallel processing.
        
        Efficiently processes large numbers of sequences by utilizing parallel
        processing while maintaining sequence integrity and error handling.
        
        Args:
            sequences: Dictionary of {seq_id: sequence_string}
            vcf_file: Path to VCF file
            reference_file: Path to reference FASTA file
            
        Returns:
            Dictionary of processed sequences with variants applied
            
        Raises:
            SequenceProcessingError: If VCF processing fails
            PipelineError: If workflow coordination fails
        """
        logger.debug("=== WORKFLOW: VCF SEQUENCE PROCESSING ===")
        logger.debug(f"Processing {len(sequences)} sequences with VCF: {vcf_file}")
        
        try:
            # Use parallel processing for large sequence sets
            if len(sequences) >= 10 and Config.NUM_PROCESSES > 1:
                processed_sequences, total_stats = cls._process_sequences_parallel(
                    sequences, vcf_file, reference_file
                )
            else:
                # Use sequential processing for small sets
                processed_sequences, total_stats = cls._process_sequences_sequential(
                    sequences, vcf_file, reference_file
                )
            
            # Report variant statistics instead of sequence count
            logger.info(f"Masked {total_stats['masked']} variable variants, substituted {total_stats['substituted']} fixed variants")
            logger.debug("=== END WORKFLOW: VCF SEQUENCE PROCESSING ===")
            
            return processed_sequences
            
        except Exception as e:
            error_msg = f"Error in VCF sequence processing workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: VCF SEQUENCE PROCESSING ===")
            
            if isinstance(e, (SequenceProcessingError, FileError, ExternalToolError)):
                raise
            else:
                raise PipelineError(error_msg) from e
    
    #############################################################################
    
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

    @classmethod
    def _process_sequences_parallel(cls, sequences: Dict[str, str], 
                                   vcf_file: str, reference_file: str) -> tuple:
        """
        Process sequences using parallel workers.
        
        Args:
            sequences: Dictionary of sequences to process
            vcf_file: Path to VCF file
            reference_file: Path to reference FASTA file
            
        Returns:
            Tuple of (processed_sequences, total_stats)
        """
        logger.debug(f"Using parallel processing with {Config.NUM_PROCESSES} workers")
        
        # Split sequences into chunks for parallel processing
        sequence_items = list(sequences.items())
        chunk_size = max(1, len(sequence_items) // (Config.NUM_PROCESSES * 2))
        chunks = [sequence_items[i:i + chunk_size] 
                  for i in range(0, len(sequence_items), chunk_size)]
        
        processed_sequences = {}
        total_stats = {'masked': 0, 'substituted': 0}
        failed_count = 0
        
        with ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            # Submit chunks for processing
            future_to_chunk = {
                executor.submit(cls._process_sequence_chunk, chunk, vcf_file, reference_file): chunk 
                for chunk in chunks
            }
            
            # Collect results
            for future in as_completed(future_to_chunk):
                try:
                    chunk_results, chunk_failed, chunk_stats = future.result()
                    processed_sequences.update(chunk_results)
                    failed_count += chunk_failed
                    total_stats['masked'] += chunk_stats['masked']
                    total_stats['substituted'] += chunk_stats['substituted']
                except Exception as e:
                    logger.error(f"Error processing sequence chunk: {e}")
                    # Add original sequences for failed chunk
                    chunk = future_to_chunk[future]
                    for seq_id, sequence in chunk:
                        processed_sequences[seq_id] = sequence
                    failed_count += len(chunk)
        
        if failed_count > 0:
            logger.debug(f"Parallel processing completed with {failed_count} failures")
        
        return processed_sequences, total_stats
    
    @classmethod
    def _process_sequence_chunk(cls, chunk: List[tuple], vcf_file: str, 
                               reference_file: str) -> tuple:
        """
        Process a chunk of sequences in a worker process.
        
        Args:
            chunk: List of (seq_id, sequence) tuples
            vcf_file: Path to VCF file
            reference_file: Path to reference FASTA file
            
        Returns:
            Tuple of (processed_sequences_dict, failed_count, chunk_stats)
        """
        processor = cls(reference_file)
        chunk_results = {}
        chunk_stats = {'masked': 0, 'substituted': 0}
        failed_count = 0
        
        for seq_id, sequence in chunk:
            try:
                modified_sequence, seq_stats = processor.process_sequence_with_vcf(
                    sequence=sequence,
                    vcf_path=vcf_file,
                    chromosome=seq_id,
                    return_stats=True
                )
                chunk_results[seq_id] = modified_sequence
                chunk_stats['masked'] += seq_stats['masked']
                chunk_stats['substituted'] += seq_stats['substituted']
                
                # Limited logging for chunk processing
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('snp_chunk_processing', interval=1000, max_initial=3)):
                    logger.debug(f"Processed sequence {seq_id}: {seq_stats['masked']} masked, {seq_stats['substituted']} substituted")
                    
            except Exception as e:
                # Keep original sequence on error
                chunk_results[seq_id] = sequence
                failed_count += 1
                
                # Limited error logging to avoid spam
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('snp_chunk_errors', interval=50, max_initial=3)):
                    logger.debug(f"Error processing sequence {seq_id}: {e}")
        
        return chunk_results, failed_count, chunk_stats
    
    @classmethod
    def _process_sequences_sequential(cls, sequences: Dict[str, str], 
                                     vcf_file: str, reference_file: str) -> tuple:
        """
        Process sequences sequentially (fallback or small datasets).
        
        Args:
            sequences: Dictionary of sequences to process
            vcf_file: Path to VCF file
            reference_file: Path to reference FASTA file
            
        Returns:
            Tuple of (processed_sequences, total_stats)
        """
        logger.debug("Using sequential processing")
        
        processor = cls(reference_file)
        processed_sequences = {}
        total_stats = {'masked': 0, 'substituted': 0}
        failed_count = 0
        processing_count = 0
        
        for seq_id, sequence in sequences.items():
            try:
                modified_sequence, seq_stats = processor.process_sequence_with_vcf(
                    sequence=sequence,
                    vcf_path=vcf_file,
                    chromosome=seq_id,
                    return_stats=True
                )
                processed_sequences[seq_id] = modified_sequence
                total_stats['masked'] += seq_stats['masked']
                total_stats['substituted'] += seq_stats['substituted']
                
                # Limited logging for sequential processing
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('snp_sequential_processing', interval=1000, max_initial=3)):
                    logger.debug(f"Processed sequence {processing_count+1}/{len(sequences)}: {seq_id}")
                    
            except Exception as e:
                failed_count += 1
                
                # Limited failure logging
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('snp_sequential_errors', interval=50, max_initial=2)):
                    logger.error(f"Error processing sequence {seq_id}: {e}")
                    
                # Keep original sequence if processing fails
                processed_sequences[seq_id] = sequence
            
            processing_count += 1
        
        if failed_count > 0:
            logger.debug(f"Sequential processing: {failed_count} failures out of {len(sequences)} sequences")
        
        return processed_sequences, total_stats
    
    def process_sequence_with_vcf(self, sequence: str, vcf_path: str, 
                                chromosome: str, return_stats: bool = False, **kwargs):
        """
        Process sequence with prepared VCF variants.
        
        Main entry point for sequence processing. Handles both fixed variant
        substitution (AF=1.0) and variable variant masking (AF<1.0) based on
        configuration settings.
        
        COORDINATE HANDLING:
        - Receives VCF with 1-based positions
        - Converts to 0-based for sequence operations
        - Returns modified sequence string (coordinate-independent)
        
        Args:
            sequence: DNA sequence to process
            vcf_path: Path to prepared VCF file
            chromosome: Chromosome/sequence identifier (should match VCF)
            return_stats: If True, return tuple of (sequence, stats)
            **kwargs: Additional processing settings from Config
            
        Returns:
            Processed sequence with variants applied, or tuple if return_stats=True
            
        Raises:
            SequenceProcessingError: If sequence processing fails
        """
        logger.debug("=== SEQUENCE PROCESSING DEBUG ===")
        logger.debug(f"Processing sequence for chromosome: {chromosome}")
        logger.debug(f"Original sequence length: {len(sequence)}")
        
        if not sequence or not isinstance(sequence, str):
            logger.debug("Empty or invalid sequence provided")
            if return_stats:
                return sequence, {'masked': 0, 'substituted': 0}
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
                if return_stats:
                    return sequence, {'masked': 0, 'substituted': 0}
                return sequence
            
            # Apply variants to sequence
            processed_sequence, stats = self._apply_variants_to_sequence(sequence, variants, **settings)
            
            # Log processing statistics
            original_len = len(sequence)
            processed_len = len(processed_sequence)
            n_count = processed_sequence.count('N')
            soft_count = sum(1 for c in processed_sequence if c.islower())
            
            logger.debug(f"Processing complete: {original_len} -> {processed_len} bp")
            logger.debug(f"Hard masked (N): {n_count}, Soft masked (lowercase): {soft_count}")
            logger.debug("=== END SEQUENCE PROCESSING DEBUG ===")
            
            if return_stats:
                return processed_sequence, stats
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
        
        Uses bcftools to extract variants from prepared VCF file.
        
        Args:
            vcf_path: Path to prepared VCF file
            chromosome: Chromosome to extract variants for
            **kwargs: Processing settings for filtering
            
        Returns:
            List of variant dictionaries with 1-based VCF coordinates
            
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
            
            # Parse variants from output with optimized logging
            variants = []
            failed_parse_count = 0
            
            lines = result.stdout.strip().split('\n') if result.stdout.strip() else []
            
            for line_num, line in enumerate(lines, 1):
                if line.strip():
                    try:
                        variant = self._parse_variant_line(line.strip())
                        if self._should_process_variant(variant, **kwargs):
                            variants.append(variant)
                            
                            # Limited logging for variant parsing
                            if (logger.isEnabledFor(logging.DEBUG) and 
                                DebugLogLimiter.should_log('snp_variant_parsing', interval=1000, max_initial=3)):
                                logger.debug(f"Parsed variant {len(variants)}: "
                                           f"{variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']} "
                                           f"(AF={variant.get('af', 'N/A')})")
                            
                    except Exception as e:
                        failed_parse_count += 1
                        
                        # Limited error logging for parsing failures
                        if (logger.isEnabledFor(logging.DEBUG) and 
                            DebugLogLimiter.should_log('snp_variant_parse_errors', interval=100, max_initial=3)):
                            logger.debug(f"Error parsing variant line {line_num}: {str(e)}")
                        continue
            
            # Summary logging
            logger.debug(f"Variant parsing summary: {len(variants)} applicable variants "
                       f"from {len(lines)} total lines, {failed_parse_count} parse failures")
            
            if failed_parse_count > 5:
                logger.debug(f"(Detailed parsing errors shown for first few failures only)")
            
            logger.debug("=== END PARSING VARIANTS FROM PREPARED VCF DEBUG ===")
            return variants
            
        except ExternalToolError:
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
            Dictionary containing variant information with 1-based VCF coordinates
            
        Raises:
            ValueError: If line format is invalid
        """
        try:
            parts = line.split('\t')
            if len(parts) < 6:
                raise ValueError(f"Invalid variant line format: expected 6 fields, got {len(parts)}")
            
            # Parse basic fields
            chrom, pos, ref, alt, qual, af = parts[:6]
            
            # Convert types - keep VCF 1-based position as-is
            pos = int(pos)  # Keep 1-based VCF coordinate
            qual = float(qual) if qual != '.' else None
            
            # Parse AF field (can be comma-separated for multiallelic)
            if af == '.' or af == '':
                af_value = None
            else:
                # Take first AF value for multiallelic sites
                af_value = float(af.split(',')[0])
            
            variant = {
                'chrom': chrom,
                'pos': pos,  # 1-based VCF position (will be converted during application)
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'af': af_value
            }
            
            return variant
            
        except (ValueError, IndexError) as e:
            error_msg = f"Error parsing variant line '{line}': {str(e)}"
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('snp_parse_line_errors', interval=50, max_initial=2)):
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
                return False
        
        # Always process if AF is not available
        if variant.get('af') is None:
            return True
        
        # Check AF threshold
        af_threshold = kwargs.get('snp_allele_frequency_threshold')
        if af_threshold is not None:
            if variant['af'] < af_threshold:
                return False
        
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
            return 'substitute'
        
        # All other variants should be masked
        return 'mask'
    
    def _apply_variants_to_sequence(self, sequence: str, variants: List[Dict], **kwargs):
        """
        Apply all variants to sequence in coordinate order.
        
        Processes variants in reverse order (high to low position) to maintain
        coordinate integrity during insertions and deletions.
        
        COORDINATE CONVERSION: VCF positions converted to 0-based here
        
        Args:
            sequence: Original DNA sequence
            variants: List of variant dictionaries with 1-based VCF positions
            **kwargs: Processing settings
            
        Returns:
            Tuple of (modified_sequence, stats)
        """
        if not variants:
            logger.debug("No variants to apply")
            return sequence, {'masked': 0, 'substituted': 0}
        
        # Sort variants by position in descending order
        # Process from end to beginning to maintain coordinate integrity
        sorted_variants = sorted(variants, key=lambda v: v['pos'], reverse=True)
        
        logger.debug(f"Applying {len(sorted_variants)} variants to sequence (length: {len(sequence)})")
        logger.debug("COORDINATE CONVERSION: VCF positions (1-based) will be converted to 0-based for sequence operations")
        
        modified_sequence = sequence
        stats = {'applied': 0, 'failed': 0, 'substituted': 0, 'masked': 0}
        
        for i, variant in enumerate(sorted_variants):
            try:
                action = self._classify_variant(variant, **kwargs)
                
                # Apply variant based on classification
                if action == 'substitute':
                    modified_sequence = self._apply_substitution(modified_sequence, variant, **kwargs)
                    stats['substituted'] += 1
                elif action == 'mask':
                    modified_sequence = self._apply_masking(modified_sequence, variant, **kwargs)
                    stats['masked'] += 1
                
                stats['applied'] += 1
                
                # Limited logging for variant application
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('snp_variant_application', interval=1000, max_initial=3)):
                    logger.debug(f"Applied variant {i+1}/{len(sorted_variants)} at VCF position "
                               f"{variant['pos']} ({action}): {variant['ref']}>{variant['alt']}")
                
            except Exception as e:
                stats['failed'] += 1
                
                # Limited error logging for application failures
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('snp_variant_apply_errors', interval=100, max_initial=3)):
                    logger.debug(f"Error applying variant at VCF position {variant['pos']}: {str(e)}")
                continue
        
        # Summary logging
        logger.debug(f"Variant application complete: {stats['applied']}/{len(sorted_variants)} applied successfully")
        logger.debug(f"  Substituted (AF=1.0): {stats['substituted']}, Masked (AF<1.0): {stats['masked']}, Failed: {stats['failed']}")
        
        return modified_sequence, stats
    
    def _apply_substitution(self, sequence: str, variant: Dict, **kwargs) -> str:
        """
        Apply variant substitution to sequence (for fixed variants AF=1.0).
        
        COORDINATE CONVERSION: VCF position (1-based) → sequence index (0-based)
        
        Args:
            sequence: DNA sequence
            variant: Variant dictionary with 1-based VCF position
            **kwargs: Processing settings
            
        Returns:
            Sequence with variant substituted
        """
        # COORDINATE CONVERSION: VCF (1-based) to Python sequence (0-based)
        vcf_pos = variant['pos']  # 1-based VCF position
        seq_index = vcf_pos - 1   # Convert to 0-based sequence index
        
        ref = variant['ref']
        alt = variant['alt']
        
        # Validate position and reference
        if seq_index < 0 or seq_index >= len(sequence):
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('snp_substitution_bounds_errors', interval=100, max_initial=2)):
                logger.debug(f"Substitution position out of bounds: VCF pos {vcf_pos} → seq index {seq_index}, seq length {len(sequence)}")
            return sequence
        
        if seq_index + len(ref) > len(sequence):
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('snp_substitution_bounds_errors', interval=100, max_initial=2)):
                logger.debug(f"Substitution extends beyond sequence: VCF pos {vcf_pos}, ref length {len(ref)}, seq length {len(sequence)}")
            return sequence
        
        # Check reference match
        sequence_ref = sequence[seq_index:seq_index + len(ref)]
        if sequence_ref.upper() != ref.upper():
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('snp_substitution_ref_mismatch', interval=50, max_initial=3)):
                logger.debug(f"Reference mismatch at VCF pos {vcf_pos}: expected '{ref}', found '{sequence_ref}'")
            return sequence
        
        # Apply substitution using 0-based sequence indexing
        new_sequence = sequence[:seq_index] + alt + sequence[seq_index + len(ref):]
        
        # Limited success logging
        if (logger.isEnabledFor(logging.DEBUG) and 
            DebugLogLimiter.should_log('snp_substitution_success', interval=1000, max_initial=2)):
            logger.debug(f"Substitution applied: VCF pos {vcf_pos} (seq[{seq_index}]) {ref}→{alt}")
        
        return new_sequence
    
    def _apply_masking(self, sequence: str, variant: Dict, **kwargs) -> str:
        """
        Apply variant masking to sequence (for variable variants AF<1.0).
        
        COORDINATE CONVERSION: VCF position (1-based) → sequence index (0-based)
        
        Args:
            sequence: DNA sequence
            variant: Variant dictionary with 1-based VCF position
            **kwargs: Processing settings
            
        Returns:
            Sequence with variant region masked
        """
        # COORDINATE CONVERSION: VCF (1-based) to Python sequence (0-based)
        vcf_pos = variant['pos']  # 1-based VCF position
        seq_index = vcf_pos - 1   # Convert to 0-based sequence index
        
        ref = variant['ref']
        
        # Validate position
        if seq_index < 0 or seq_index >= len(sequence):
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('snp_masking_bounds_errors', interval=100, max_initial=2)):
                logger.debug(f"Masking position out of bounds: VCF pos {vcf_pos} → seq index {seq_index}, seq length {len(sequence)}")
            return sequence
        
        if seq_index + len(ref) > len(sequence):
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('snp_masking_bounds_errors', interval=100, max_initial=2)):
                logger.debug(f"Masking extends beyond sequence: VCF pos {vcf_pos}, ref length {len(ref)}, seq length {len(sequence)}")
            return sequence
        
        # Determine masking strategy
        use_soft_masking = kwargs.get('snp_use_soft_masking', False)
        flanking_mask_size = kwargs.get('snp_flanking_mask_size', 0)
        
        # Calculate masking region using 0-based coordinates
        start_index = max(0, seq_index - flanking_mask_size)
        end_index = min(len(sequence), seq_index + len(ref) + flanking_mask_size)
        
        # Apply masking
        if use_soft_masking:
            # Soft masking - convert to lowercase
            masked_region = sequence[start_index:end_index].lower()
        else:
            # Hard masking - replace with 'N'
            masked_region = 'N' * (end_index - start_index)
        
        new_sequence = sequence[:start_index] + masked_region + sequence[end_index:]
        
        # Limited success logging
        if (logger.isEnabledFor(logging.DEBUG) and 
            DebugLogLimiter.should_log('snp_masking_success', interval=1000, max_initial=2)):
            mask_type = "soft" if use_soft_masking else "hard"
            logger.debug(f"{mask_type.capitalize()} masking applied: VCF pos {vcf_pos} "
                       f"(seq[{start_index}:{end_index}]) with {flanking_mask_size}bp flanking")
        
        return new_sequence