#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fragment-based SNP masking processor with chromosome-by-chromosome processing to avoid memory exhaustion.

This version processes chromosomes sequentially to handle very large VCF files:
1. Processes one chromosome at a time to minimize memory usage
2. Properly handles coordinate shifts from insertions/deletions
3. Uses tolerance for reference validation to handle edge cases
4. Implements proper variant grouping and processing order
5. Forces garbage collection between chromosomes

COORDINATE SYSTEM:
- VCF input: 1-based positions (standard VCF format)
- Fragment coordinates: 0-based half-open intervals [start, end)
- Internal processing: 0-based indexing for sequence operations
"""

import os
import subprocess
import logging
import gc
from typing import Dict, List, Tuple
from collections import defaultdict
from tqdm import tqdm

# Import package modules
from ..config import Config, FileError, ExternalToolError, SequenceProcessingError, DebugLogLimiter

# Set up module logger
logger = logging.getLogger(__name__)


class SNPMaskingProcessor:
    """
    Fragment-based SNP processor with chromosome-by-chromosome processing to avoid memory exhaustion.
    
    Key improvements:
    - Sequential chromosome processing to minimize memory usage
    - Fixed coordinate conversion between VCF and fragment positions
    - Proper handling of multiple variants at same position
    - Enhanced reference validation with tolerance for modified sequences
    - Correct processing order to avoid coordinate shift issues
    - Memory cleanup between chromosomes
    """
    
    @classmethod
    def process_fragments_with_vcf(cls, restriction_fragments: List[Dict], 
                                 vcf_file: str, reference_file: str) -> List[Dict]:
        """
        Process restriction fragments with VCF variants using chromosome-by-chromosome processing
        to avoid memory exhaustion.
        
        Args:
            restriction_fragments: List of fragment dictionaries with 0-based coordinates
            vcf_file: Path to VCF file
            reference_file: Path to reference FASTA file
            
        Returns:
            List of processed fragments with SNP variants applied
        """
        logger.info("Masking Variants with VCF file...")
        logger.debug("=== WORKFLOW: CHROMOSOME-BY-CHROMOSOME VCF PROCESSING ===")
        logger.debug(f"Processing {len(restriction_fragments)} fragments")
        logger.debug(f"VCF file: {vcf_file}")
        logger.debug(f"Reference file: {reference_file}")
        
        try:
            # Reset log limiters for a fresh run
            DebugLogLimiter.reset_all()

            if not restriction_fragments:
                logger.warning("No fragments provided for SNP processing")
                return []
            
            # Initialize processor
            processor = cls(reference_file)
            
            # Step 1: Group fragments by chromosome for efficient processing
            fragments_by_chr = cls._group_fragments_by_chromosome(restriction_fragments)
            chroms_in_fragments = set(fragments_by_chr.keys())

            # Step 1.5: Identify chromosomes to process based on size
            min_size = Config.MIN_CHROMOSOME_SIZE
            chroms_to_process = chroms_in_fragments
            if min_size and min_size > 0:
                logger.debug(f"Applying minimum chromosome size filter: {min_size} bp.")
                try:
                    chromosome_sizes = processor._get_chromosome_sizes()
                    large_enough_chroms = {
                        chrom for chrom, size in chromosome_sizes.items() if size >= min_size
                    }
                    
                    # The set of chromosomes we will actually query in the VCF
                    chroms_to_process = chroms_in_fragments.intersection(large_enough_chroms)
                    
                    # Log which chromosomes present in fragments are being skipped
                    skipped_chroms = chroms_in_fragments - chroms_to_process
                    if skipped_chroms:
                        logger.debug(
                            f"Skipping SNP processing for {len(skipped_chroms)} chromosomes "
                            f"that are smaller than {min_size} bp."
                        )
                        logger.debug(
                            f"Skipped chromosomes (sample): {', '.join(list(skipped_chroms)[:5])}"
                        )

                except FileError as e:
                    logger.error(f"Could not filter chromosomes by size: {e}")
                    logger.warning("Proceeding without chromosome size filtering.")
                    chroms_to_process = chroms_in_fragments # Reset to all chroms on error
                except Exception as e:
                    logger.error(f"An unexpected error occurred during chromosome size filtering: {e}", exc_info=True)
                    logger.warning("Proceeding without chromosome size filtering.")
                    chroms_to_process = chroms_in_fragments # Reset to all chroms on error

            # MAIN CHANGE: Process each chromosome separately to avoid memory exhaustion
            processed_fragments = []
            total_stats = {
                'masked': 0, 'substituted': 0, 'applied': 0, 'failed': 0, 
                'reference_mismatches': 0, 'coordinate_errors': 0, 'validation_failures': 0
            }
            
            # Process chromosomes one by one
            if chroms_to_process:
                for i, chromosome in enumerate(sorted(chroms_to_process), 1):
                    # Get fragments for this chromosome only
                    chr_fragments = fragments_by_chr[chromosome]
                    logger.debug(f"  Fragments in {chromosome}: {len(chr_fragments)}")
                    
                    # Load variants for this chromosome only
                    chr_variants = processor._parse_variants_single_chromosome(vcf_file, chromosome)
                    logger.debug(f"  Variants in {chromosome}: {len(chr_variants)}")
                    
                    # Process fragments for this chromosome
                    chr_processed_fragments = cls._process_chromosome_fragments(
                        chr_fragments, chr_variants, chromosome
                    )
                    
                    # Accumulate results
                    processed_fragments.extend(chr_processed_fragments)
                    
                    # Accumulate statistics
                    for fragment in chr_processed_fragments:
                        fragment_stats = fragment.get('_processing_stats', {})
                        for key in total_stats:
                            total_stats[key] += fragment_stats.get(key, 0)
                    
                    # Force garbage collection after each chromosome to free memory
                    del chr_variants
                    del chr_processed_fragments
                    gc.collect()
                    
                    logger.debug(f"  Completed {chromosome}. Memory cleanup performed.")
            else:
                logger.info("No chromosomes met the size criteria; skipping VCF parsing.")
            
            # Process fragments from chromosomes that don't meet size criteria (no variants applied)
            skipped_chroms = chroms_in_fragments - chroms_to_process
            if skipped_chroms:
                logger.debug(f"Processing {len(skipped_chroms)} skipped chromosomes without variant application")
                for chromosome in skipped_chroms:
                    chr_fragments = fragments_by_chr[chromosome]
                    for fragment in chr_fragments:
                        # No variant processing, just copy the fragment with empty stats
                        processed_fragment = fragment.copy()
                        processed_fragment['_processing_stats'] = {
                            'applied': 0, 'failed': 0, 'substituted': 0, 'masked': 0, 
                            'validation_failures': 0, 'reference_mismatches': 0, 'coordinate_errors': 0
                        }
                        processed_fragments.append(processed_fragment)
            
            # Report comprehensive statistics
            logger.info(f"Masked {total_stats['masked']} variable variants, substituted {total_stats['substituted']} fixed variants")
            
            # Enhanced error reporting
            if total_stats.get('reference_mismatches', 0) > 0:
                logger.warning(f"Reference validation failed for {total_stats['reference_mismatches']} variants")
                logger.warning("This usually indicates VCF/reference version mismatch or coordinate issues")
                
            if total_stats.get('coordinate_errors', 0) > 0:
                logger.warning(f"Coordinate validation failed for {total_stats['coordinate_errors']} variants")
                
            logger.debug("=== END WORKFLOW: CHROMOSOME-BY-CHROMOSOME VCF PROCESSING ===")
            return processed_fragments
            
        except Exception as e:
            error_msg = f"Error in chromosome-by-chromosome VCF processing: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: CHROMOSOME-BY-CHROMOSOME VCF PROCESSING ===")
            
            if isinstance(e, (SequenceProcessingError, FileError, ExternalToolError)):
                raise
            else:
                raise SequenceProcessingError(error_msg) from e
    
    @classmethod
    def _debug_fragment_info(cls, fragments: List[Dict]):
        """Debug fragment information for troubleshooting."""
        if not logger.isEnabledFor(logging.DEBUG):
            return
            
        logger.debug("=== FRAGMENT DEBUG INFO ===")
        logger.debug(f"Total fragments: {len(fragments)}")
        
        # Sample fragment analysis
        for i, fragment in enumerate(fragments[:3]):
            frag_id = fragment.get('id', 'unknown')
            chromosome = fragment.get('chromosome', 'unknown')
            start = fragment.get('start', 'unknown')
            end = fragment.get('end', 'unknown')
            seq_len = len(fragment.get('sequence', ''))
            coord_span = end - start if isinstance(start, int) and isinstance(end, int) else 'unknown'
            
            logger.debug(f"Fragment {i+1}: {frag_id}")
            logger.debug(f"  Chromosome: {chromosome}")
            logger.debug(f"  Coordinates: {start}-{end} (span: {coord_span})")
            logger.debug(f"  Sequence length: {seq_len}")
            
            if isinstance(coord_span, int) and coord_span != seq_len:
                logger.warning(f"  COORDINATE MISMATCH: span={coord_span} != seq_len={seq_len}")
        
        if len(fragments) > 3:
            logger.debug(f"  ... and {len(fragments) - 3} more fragments")
        
        logger.debug("=== END FRAGMENT DEBUG INFO ===")
    
    @classmethod 
    def _debug_variant_info(cls, all_variants: Dict[str, List[Dict]]):
        """Debug variant information for troubleshooting."""
        if not logger.isEnabledFor(logging.DEBUG):
            return
            
        logger.debug("=== VARIANT DEBUG INFO ===")
        total_variants = sum(len(variants) for variants in all_variants.values())
        logger.debug(f"Total variants loaded: {total_variants}")
        
        for chromosome, variants in all_variants.items():
            if variants:
                logger.debug(f"Chromosome {chromosome}: {len(variants)} variants")
                
                # Sample variant analysis
                for i, variant in enumerate(variants[:2]):
                    pos = variant.get('pos', 'unknown')
                    ref = variant.get('ref', 'unknown')
                    alt = variant.get('alt', 'unknown')
                    af = variant.get('af', 'unknown')
                    
                    logger.debug(f"  Variant {i+1}: pos={pos}, ref='{ref}', alt='{alt}', af={af}")
                    
                    # Check for potential issues
                    if isinstance(ref, str):
                        if ref.endswith(' ') or ref.startswith(' '):
                            logger.debug(f"  WARNING: Reference has whitespace: '{ref}'")
                        if len(ref) > 1:
                            logger.debug(f"  NOTE: Multi-character reference: '{ref}' (length={len(ref)})")
        
        logger.debug("=== END VARIANT DEBUG INFO ===")
    
    @classmethod
    def _group_fragments_by_chromosome(cls, restriction_fragments: List[Dict]) -> Dict[str, List[Dict]]:
        """Group fragments by chromosome with validation."""
        fragments_by_chr = defaultdict(list)
        
        for fragment in restriction_fragments:
            chromosome = fragment.get('chromosome')
            if chromosome:
                fragments_by_chr[chromosome].append(fragment)
            else:
                logger.warning(f"Fragment {fragment.get('id', 'unknown')} missing chromosome information")
        
        logger.debug(f"Grouped fragments across {len(fragments_by_chr)} chromosomes")
        return dict(fragments_by_chr)
    
    @classmethod
    def _process_chromosome_fragments(cls, chr_fragments: List[Dict], 
                                    chr_variants: List[Dict], 
                                    chromosome: str) -> List[Dict]:
        """Process fragments for a single chromosome with optimized variant lookup."""
        processed_fragments = []
        
        if not chr_variants:
            for fragment in tqdm(chr_fragments, desc=f"Processing {chromosome} (no variants)", unit="frag"):
                processed_fragment = cls._process_single_fragment_with_variants(fragment, [])
                processed_fragments.append(processed_fragment)
            return processed_fragments
        
        # Create spatial index for fast variant lookup
        logger.debug(f"Creating spatial index for {len(chr_variants)} variants in {chromosome}")
        
        # Simple binning approach for immediate speedup
        BIN_SIZE = 10000  # 10kb bins
        variant_bins = {}
        
        for variant in chr_variants:
            genomic_pos = variant['pos'] - 1  # Convert to 0-based
            bin_id = genomic_pos // BIN_SIZE
            if bin_id not in variant_bins:
                variant_bins[bin_id] = []
            variant_bins[bin_id].append(variant)
        
        with tqdm(chr_fragments, desc=f"Processing {chromosome}", unit="frag") as pbar:
            for fragment in pbar:
                # Fast variant lookup using bins
                fragment_start = fragment.get('start', 0)
                fragment_end = fragment.get('end', 0)
                
                start_bin = fragment_start // BIN_SIZE
                end_bin = fragment_end // BIN_SIZE
                
                # Get variants from overlapping bins only
                relevant_variants = []
                for bin_id in range(start_bin, end_bin + 1):
                    if bin_id in variant_bins:
                        relevant_variants.extend(variant_bins[bin_id])
                
                # Map only relevant variants to fragment
                mapped_variants = cls._map_variants_to_fragment(fragment, relevant_variants)
                
                # Process fragment with variants
                processed_fragment = cls._process_single_fragment_with_variants(fragment, mapped_variants)
                processed_fragments.append(processed_fragment)
        
        return processed_fragments
    
    def __init__(self, reference_file: str):
        """Initialize processor with reference file validation."""
        logger.debug(f"Initializing fragment-based SNPMaskingProcessor with reference: {reference_file}")
        
        if not os.path.exists(reference_file):
            error_msg = f"Reference FASTA file not found: {reference_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        self.reference_file = reference_file
    
    def _get_chromosome_sizes(self) -> Dict[str, int]:
        """
        Parse reference FASTA index to get chromosome sizes.

        Returns:
            Dictionary mapping chromosome names to their lengths.

        Raises:
            FileError: If the FASTA index file (.fai) is not found.
        """
        index_file = f"{self.reference_file}.fai"
        logger.debug(f"Reading chromosome sizes from FASTA index: {index_file}")
        
        if not os.path.exists(index_file):
            error_msg = (
                f"FASTA index file not found: {index_file}. "
                "Please index the reference file with 'samtools faidx' before running."
            )
            logger.error(error_msg)
            raise FileError(error_msg)
            
        sizes = {}
        try:
            with open(index_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        chrom_name, length = parts[0], parts[1]
                        sizes[chrom_name] = int(length)
        except (IOError, ValueError) as e:
            error_msg = f"Error reading or parsing FASTA index file {index_file}: {e}"
            logger.error(error_msg)
            raise FileError(error_msg) from e
            
        logger.debug(f"Found sizes for {len(sizes)} chromosomes.")
        return sizes

    def _parse_variants_single_chromosome(self, vcf_file: str, chromosome: str) -> List[Dict]:
        """Parse variants for a single chromosome to minimize memory usage."""
        logger.debug(f"Loading variants for chromosome: {chromosome}")
        
        try:
            cmd = [
                'bcftools', 'query',
                '-r', chromosome,  # Single chromosome only
                '-f', '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AF\n',
                vcf_file
            ]
            
            logger.debug(f"Running bcftools command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            
            if result.returncode != 0:
                logger.warning(f"bcftools query failed for chromosome {chromosome}: {result.stderr}")
                return []
            
            # Parse variants
            variants = []
            lines = result.stdout.strip().split('\n') if result.stdout.strip() else []
            
            parsing_errors = 0
            for line_num, line in enumerate(lines, 1):
                if not line.strip():
                    continue
                
                try:
                    variant = self._parse_variant_line(line.strip())
                    if self._validate_variant(variant, line_num):
                        if self._should_process_variant(variant):
                            variants.append(variant)
                except Exception as e:
                    parsing_errors += 1
                    if parsing_errors <= 5:  # Show first 5 parsing errors per chromosome
                        logger.debug(f"Error parsing variant line {line_num} in {chromosome}: {str(e)}")
            
            logger.debug(f"Parsed {len(variants)} variants from {len(lines)} lines for {chromosome} ({parsing_errors} parsing errors)")
            return variants
            
        except Exception as e:
            logger.error(f"Error parsing variants for chromosome {chromosome}: {str(e)}")
            return []
    
    def _parse_variant_line(self, line: str) -> Dict:
        """Parse variant line with enhanced validation."""
        try:
            parts = line.split('\t')
            if len(parts) < 6:
                raise ValueError(f"Invalid variant line format: expected 6 fields, got {len(parts)}")
            
            chrom, pos, ref, alt, qual, af = parts[:6]
            
            # Parse position
            pos = int(pos)
            
            # Parse quality
            qual = float(qual) if qual != '.' else None
            
            # Parse allele frequency
            if af == '.' or af == '':
                af_value = None
            else:
                af_value = float(af.split(',')[0])
            
            variant = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'af': af_value,
                '_original_line': line  # Keep for debugging
            }
            
            return variant
            
        except (ValueError, IndexError) as e:
            error_msg = f"Error parsing variant line '{line}': {str(e)}"
            raise ValueError(error_msg) from e
    
    def _validate_variant(self, variant: Dict, line_num: int) -> bool:
        """Validate variant with detailed error reporting."""
        ref = variant.get('ref', '')
        alt = variant.get('alt', '')
        pos = variant.get('pos')
        
        # Check for whitespace issues
        if ref != ref.strip():
            logger.debug(f"Line {line_num}: Reference has whitespace: '{ref}' -> '{ref.strip()}'")
            variant['ref'] = ref.strip()  # Clean it
            
        if alt != alt.strip():
            logger.warning(f"Line {line_num}: Alternative has whitespace: '{alt}' -> '{alt.strip()}'")
            variant['alt'] = alt.strip()  # Clean it
        
        # Check for empty sequences
        if not variant['ref'] or not variant['alt']:
            logger.warning(f"Line {line_num}: Empty ref/alt sequences: ref='{variant['ref']}', alt='{variant['alt']}'")
            return False
        
        # Check position validity
        if pos <= 0:
            logger.warning(f"Line {line_num}: Invalid position: {pos}")
            return False
        
        return True
    
    def _should_process_variant(self, variant: Dict) -> bool:
        """Determine if variant should be processed."""
        # Check quality threshold
        qual_threshold = getattr(Config, 'VCF_QUALITY_THRESHOLD', None)
        if qual_threshold is not None and variant.get('qual') is not None:
            if variant['qual'] < qual_threshold:
                return False
        
        # Check AF threshold
        af_threshold = getattr(Config, 'VCF_ALLELE_FREQUENCY_THRESHOLD', None)
        if af_threshold is not None and variant.get('af') is not None:
            if variant['af'] < af_threshold:
                return False
        
        return True
    
    @classmethod
    def _map_variants_to_fragment(cls, fragment: Dict, chr_variants: List[Dict]) -> List[Dict]:
        """
        Map variants to fragment with CORRECTED coordinate handling.
        
        This is the key fix: properly handle the coordinate conversion from
        VCF genomic positions to fragment-relative positions.
        """
        fragment_start = fragment.get('start', 0)  # Fragment start in genomic coordinates (0-based)
        fragment_end = fragment.get('end', 0)      # Fragment end in genomic coordinates (0-based)
        fragment_id = fragment.get('id', 'unknown')
        fragment_variants = []
        
        logger.debug(f"Mapping variants for fragment {fragment_id} (genomic: {fragment_start}-{fragment_end})")
        
        for variant in chr_variants:
            # COORDINATE CONVERSION: VCF position (1-based) to genomic position (0-based)
            vcf_pos_1based = variant['pos']         # VCF position (1-based)
            genomic_pos_0based = vcf_pos_1based - 1 # Convert to 0-based genomic coordinate
            
            # Check if variant falls within fragment boundaries (using 0-based coordinates)
            # For fragment [start, end), variant at genomic_pos should satisfy: start <= genomic_pos < end
            if fragment_start <= genomic_pos_0based < fragment_end:
                # Calculate position within fragment sequence (0-based relative to fragment start)
                fragment_relative_pos = genomic_pos_0based - fragment_start
                
                # Validate fragment relative position
                fragment_seq_len = len(fragment.get('sequence', ''))
                if fragment_relative_pos < 0 or fragment_relative_pos >= fragment_seq_len:
                    if DebugLogLimiter.should_log('invalid_frag_pos', max_initial=5, interval=500):
                        logger.debug(f"  Variant at genomic pos {genomic_pos_0based} maps to invalid fragment position {fragment_relative_pos} for fragment {fragment_id} (len: {fragment_seq_len}).")
                    continue
                
                # Create fragment-specific variant with comprehensive coordinate info
                fragment_variant = variant.copy()
                fragment_variant['fragment_pos_0based'] = fragment_relative_pos
                fragment_variant['original_genomic_pos_0based'] = genomic_pos_0based
                fragment_variant['original_vcf_pos_1based'] = vcf_pos_1based
                fragment_variant['fragment_id'] = fragment_id
                fragment_variant['fragment_start'] = fragment_start
                fragment_variant['fragment_end'] = fragment_end
                
                fragment_variants.append(fragment_variant)
                
                # Debug coordinate mapping
                if DebugLogLimiter.should_log('variant_mapping', max_initial=10, interval=1000):
                    logger.debug(f"  Mapped variant: VCF pos {vcf_pos_1based} -> genomic {genomic_pos_0based} -> fragment pos {fragment_relative_pos}")
        
        logger.debug(f"  Total variants mapped to fragment {fragment_id}: {len(fragment_variants)}")
        return fragment_variants
    
    @classmethod
    def _process_single_fragment_with_variants(cls, fragment: Dict, 
                                             fragment_variants: List[Dict]) -> Dict:
        """Process fragment with comprehensive error tracking and coordinate updates."""
        processed_fragment = fragment.copy()
        sequence = fragment.get('sequence', '')
        
        # Initialize stats and exit early for fragments with no variants to apply
        if not fragment_variants or not sequence:
            stats = {
                'applied': 0, 'failed': 0, 'substituted': 0, 'masked': 0, 
                'validation_failures': 0, 'reference_mismatches': 0, 'coordinate_errors': 0
            }
            if not sequence and fragment_variants:
                 logger.warning(f"Fragment {fragment.get('id', 'unknown')} has variants but no sequence")
                 stats['failed'] = len(fragment_variants)
            processed_fragment['_processing_stats'] = stats
            return processed_fragment
        
        # Map variants to their positions within the fragment sequence
        mapped_variants = cls._map_variants_to_fragment(processed_fragment, fragment_variants)

        # Apply variants using the robust rebuilding method
        modified_sequence, processing_stats = cls._apply_variants_to_fragment_sequence(
            sequence, mapped_variants, fragment.get('id', 'unknown')
        )
        
        # *** BEGIN FIX: Update coordinates to reflect length changes from indels ***
        original_length = len(sequence)
        new_length = len(modified_sequence)
        length_change = new_length - original_length
        
        processed_fragment['sequence'] = modified_sequence
        processed_fragment['length'] = new_length
        if length_change != 0:
            # Update the end coordinate to maintain consistency
            original_end = processed_fragment.get('end', original_length)
            processed_fragment['end'] = original_end + length_change
            logger.debug(f"Fragment {fragment.get('id', 'unknown')} length changed by {length_change} bp due to indel. "
                         f"End coordinate updated to {processed_fragment['end']}.")
        # *** END FIX ***
        
        processed_fragment['_processing_stats'] = processing_stats
        
        return processed_fragment
    
    @classmethod
    def _apply_variants_to_fragment_sequence(cls, sequence: str, fragment_variants: List[Dict],
                                           fragment_id: str) -> Tuple[str, Dict]:
        """
        FIXED: Apply variants by rebuilding the sequence, not modifying it in-place.
        This is a robust method that avoids coordinate shift errors.
        """
        stats = {
            'applied': 0, 'failed': 0, 'substituted': 0, 'masked': 0,
            'validation_failures': 0, 'coordinate_errors': 0, 'reference_mismatches': 0
        }
        
        if not fragment_variants:
            return sequence, stats

        # Sort variants by position in FORWARD order for sequential processing
        sorted_variants = sorted(fragment_variants, key=lambda v: v['fragment_pos_0based'])
        
        new_sequence_parts = []
        last_pos = 0

        for variant in sorted_variants:
            pos_in_fragment = variant['fragment_pos_0based']
            ref = variant['ref']
            alt = variant['alt']
            
            # Defensive check for overlapping variants (simple case)
            if pos_in_fragment < last_pos:
                if DebugLogLimiter.should_log('overlapping_variant', max_initial=5, interval=100):
                    logger.debug(f"Skipping overlapping variant at {pos_in_fragment} in fragment {fragment_id}")
                stats['failed'] += 1
                continue
            
            # Reference validation is now always against the original, pristine sequence
            fragment_ref = sequence[pos_in_fragment : pos_in_fragment + len(ref)]
            
            if not cls._reference_match(fragment_ref, ref, variant, fragment_id):
                stats['reference_mismatches'] += 1
                stats['failed'] += 1
                continue
            
            # Append the sequence slice from the end of the last variant to the start of this one
            new_sequence_parts.append(sequence[last_pos:pos_in_fragment])
            
            # Apply the variant by appending the new allele
            action = cls._classify_variant(variant)
            if action == 'substitute':
                new_sequence_parts.append(alt)
                stats['substituted'] += 1
            elif action == 'mask':
                # Mask by replacing the reference allele with 'N's of the same length
                # This prevents length changes for non-substituted heterozygous indels
                new_sequence_parts.append('N' * len(ref))
                stats['masked'] += 1
            
            stats['applied'] += 1
            
            # Update the last position to be the end of the current reference allele
            last_pos = pos_in_fragment + len(ref)

        # Append the final part of the sequence after the last variant
        new_sequence_parts.append(sequence[last_pos:])

        logger.debug(f"Fragment {fragment_id} processing complete: {stats}")
        return "".join(new_sequence_parts), stats

    @classmethod
    def _reference_match(cls, fragment_ref: str, vcf_ref: str, variant: Dict, fragment_id: str, **kwargs) -> bool:
        """
        Simple reference matching with detailed logging.
        The `has_multiple_variants` tolerance is less critical with the new robust method
        but can be kept for edge cases.
        """
        fragment_ref_clean = fragment_ref.upper().strip()
        vcf_ref_clean = vcf_ref.upper().strip()
        
        if fragment_ref_clean == vcf_ref_clean:
            return True
        
        # Enhanced debugging for mismatches
        if DebugLogLimiter.should_log('ref_mismatch', max_initial=10, interval=100):
            genomic_pos = variant.get('original_genomic_pos_0based', 'unknown')
            vcf_pos = variant.get('original_vcf_pos_1based', 'unknown')
            fragment_pos = variant.get('fragment_pos_0based', 'unknown')
            
            logger.debug(f"REFERENCE MISMATCH in fragment {fragment_id}:")
            logger.debug(f"  VCF position: {vcf_pos} (1-based)")
            logger.debug(f"  Genomic position: {genomic_pos} (0-based)")
            logger.debug(f"  Fragment position: {fragment_pos} (0-based)")
            logger.debug(f"  Expected (VCF): '{vcf_ref}' (cleaned: '{vcf_ref_clean}', length: {len(vcf_ref_clean)})")
            logger.debug(f"  Found (fragment): '{fragment_ref}' (cleaned: '{fragment_ref_clean}', length: {len(fragment_ref_clean)})")
        
        return False
    
    @classmethod
    def _classify_variant(cls, variant: Dict) -> str:
        """Classify variant based on allele frequency."""
        af = variant.get('af')
        
        # If AF indicates the variant is essentially fixed in the population, substitute it.
        if af is not None and af > 0.95:
            return 'substitute'
        
        # Otherwise, mask the site as it's variable.
        return 'mask'