#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence processing module for ddPrimer pipeline.

Handles sequence processing and filtering operations including restriction site cutting,
gene overlap filtering, sequence validation and coordinate handling for primer design workflows.

Contains functionality for:
1. Restriction site cutting with configurable patterns
2. Gene overlap filtering and fragment generation
3. Sequence validation and coordinate handling
4. Fragment length filtering and optimization

This module integrates with the broader ddPrimer pipeline to provide
robust sequence manipulation capabilities for primer design workflows.

COORDINATE SYSTEM:
- Internal coordinates: 0-based, half-open intervals [start, end)
- All fragments use consistent 0-based coordinate system
- No dual coordinate systems - single consistent format
"""

import re
import logging
from typing import Dict, List

# Import package modules
from ..config import Config, SequenceProcessingError, CoordinateValidationError, DebugLogLimiter

# Set up module logger
logger = logging.getLogger(__name__)


class SequenceProcessor:
    """
    Handles sequence processing and filtering operations.
    
    This class provides comprehensive sequence manipulation capabilities
    including restriction site analysis, fragment generation, and
    coordinate-based filtering for primer design workflows.
    
    INTERNAL COORDINATE SYSTEM: 0-based, half-open intervals [start, end)
    - start: inclusive (first nucleotide position)
    - end: exclusive (position after last nucleotide)  
    - length = end - start
    
    Features:
        - Configurable restriction site cutting
        - Fragment length validation
        - Consistent 0-based coordinate system for all operations
        - Comprehensive coordinate validation
        - Gene overlap detection and filtering
        
    Example:
        >>> fragments = SequenceProcessor.cut_at_restriction_sites(sequences)
        >>> # All fragments will have 0-based coordinates
    """
    
    #############################################################################
    #                           Workflow Wrappers
    #############################################################################
    
    @staticmethod
    def process_restriction_sites_workflow(processed_sequences: Dict[str, str]) -> List[Dict]:
        """
        Cut sequences at restriction sites for workflow integration.
        
        Processes sequences (which may already be masked/substituted) to identify 
        and cut at restriction sites, creating fragments suitable for primer design.
        All output fragments use 0-based coordinate system.
        
        Args:
            processed_sequences: Dictionary of sequences (masked and/or with SNPs substituted)
            
        Returns:
            List of restriction fragments with 0-based coordinates
            
        Raises:
            SequenceProcessingError: If there's an error in restriction site processing
        """
        logger.debug("=== WORKFLOW: RESTRICTION SITE PROCESSING ===")
        logger.debug(f"Processing {len(processed_sequences)} sequences for restriction sites")
        logger.debug("COORDINATE SYSTEM: All output fragments will use 0-based coordinates")

        try:
            restriction_fragments = SequenceProcessor.cut_at_restriction_sites(processed_sequences)
            
            # Validate all generated fragments
            valid_fragments = []
            validation_errors = 0
            
            for fragment in restriction_fragments:
                validation_result = SequenceProcessor._validate_fragment_coordinates(fragment)
                if validation_result["valid"]:
                    valid_fragments.append(fragment)
                else:
                    validation_errors += 1
                    if validation_errors <= 5:  # Log first 5 validation errors
                        logger.warning(f"Fragment {fragment.get('id', 'unknown')} failed validation: {validation_result['errors']}")
            
            if validation_errors > 0:
                logger.warning(f"Total validation errors: {validation_errors} fragments had invalid coordinates")
            
            logger.debug(f"Generated {len(valid_fragments)} valid fragments after restriction site cutting")
            logger.debug("=== END WORKFLOW: RESTRICTION SITE PROCESSING ===")
            return valid_fragments
            
        except Exception as e:
            error_msg = f"Error in restriction site processing workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: RESTRICTION SITE PROCESSING ===")
            raise SequenceProcessingError(error_msg) from e
    
    #############################################################################
    #                           Coordinate Validation
    #############################################################################
    
    @staticmethod
    def _validate_fragment_coordinates(fragment: Dict) -> Dict[str, any]:
        """
        Validate fragment coordinates for consistency with 0-based system.
        
        Args:
            fragment: Fragment dictionary with 0-based coordinates
            
        Returns:
            Dictionary with validation results
        """
        errors = []
        warnings = []
        
        # Check required fields
        start = fragment.get("start")
        end = fragment.get("end") 
        sequence = fragment.get("sequence", "")
        chromosome = fragment.get("chromosome")
        fragment_id = fragment.get("id")
        
        if start is None:
            errors.append("Missing 'start' coordinate")
        if end is None:
            errors.append("Missing 'end' coordinate")
        if not chromosome:
            errors.append("Missing 'chromosome' field")
        if not sequence:
            errors.append("Missing or empty 'sequence' field")
        if not fragment_id:
            warnings.append("Missing 'id' field")
            
        # Validate coordinate logic for 0-based system
        if start is not None and end is not None:
            if start < 0:
                errors.append(f"Start coordinate cannot be negative: {start}")
            if end <= start:
                errors.append(f"End coordinate ({end}) must be greater than start ({start}) in 0-based system")
            
            expected_length = end - start
            actual_length = len(sequence)
            
            if sequence and expected_length != actual_length:
                errors.append(f"Coordinate span ({expected_length}) doesn't match sequence length ({actual_length})")
            
            # Check for minimum segment length
            if expected_length > 0 and expected_length < Config.MIN_SEGMENT_LENGTH:
                warnings.append(f"Fragment length ({expected_length}) below minimum ({Config.MIN_SEGMENT_LENGTH})")
        
        # Check for deprecated coordinate fields that should not be present
        deprecated_fields = ["genomic_start", "genomic_end", "seq_start", "seq_end"]
        found_deprecated = [field for field in deprecated_fields if field in fragment]
        if found_deprecated:
            warnings.append(f"Fragment contains deprecated coordinate fields: {found_deprecated}")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "warnings": warnings,
            "coordinates": {
                "start": start,
                "end": end,
                "length": end - start if (start is not None and end is not None) else None,
                "sequence_length": len(sequence)
            },
            "coordinate_system": "0-based_half_open"
        }
    
    @staticmethod
    def _validate_sequence_boundaries(sequence: str, start: int, end: int) -> Dict[str, any]:
        """
        Validate that start/end coordinates are valid for the given sequence.
        
        Args:
            sequence: DNA sequence string
            start: Start position (0-based)
            end: End position (0-based, exclusive)
            
        Returns:
            Dictionary with validation results
        """
        errors = []
        
        if not isinstance(sequence, str):
            errors.append("Sequence must be a string")
            return {"valid": False, "errors": errors}
        
        seq_length = len(sequence)
        
        if start < 0:
            errors.append(f"Start position ({start}) cannot be negative")
        if end < 0:
            errors.append(f"End position ({end}) cannot be negative")
        if start >= end:
            errors.append(f"Start ({start}) must be less than end ({end})")
        if start >= seq_length:
            errors.append(f"Start position ({start}) exceeds sequence length ({seq_length})")
        if end > seq_length:
            errors.append(f"End position ({end}) exceeds sequence length ({seq_length})")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "boundaries": {
                "start": start,
                "end": end,
                "sequence_length": seq_length,
                "extracted_length": end - start if start < end else 0
            }
        }
    
    #############################################################################
    #                           Restriction Site Processing
    #############################################################################
    
    @staticmethod
    def cut_at_restriction_sites(sequences, restriction_site=None):
        """
        Cut sequences at restriction sites to generate fragments with 0-based coordinates.
        
        Processes DNA sequences by identifying restriction enzyme cut sites
        and generating fragments suitable for primer design, with consistent
        0-based coordinate system throughout.
        
        Args:
            sequences: Dictionary mapping sequence IDs to DNA sequences
            restriction_site: Restriction site pattern, defaults to Config.RESTRICTION_SITE
            
        Returns:
            List of fragment dictionaries with 0-based coordinates
            
        Raises:
            SequenceProcessingError: If restriction site pattern is invalid or processing fails
            
        Example:
            >>> sequences = {"chr1": "GAATTCATCGAATTCGCTA"}
            >>> fragments = SequenceProcessor.cut_at_restriction_sites(sequences, "GAATTC")
            >>> # All fragments will have 0-based coordinates: start, end, length
        """
        logger.info("\nFiltering sequences by restriction sites...")       
        logger.debug(f"=== RESTRICTION SITE CUTTING ===")
        logger.debug(f"Processing {len(sequences)} sequences for restriction site cutting")
        logger.debug("COORDINATE SYSTEM: All fragments will use 0-based coordinates")
        
        if restriction_site is None:
            restriction_site = Config.RESTRICTION_SITE
        
        if not restriction_site:
            logger.debug("No restriction site pattern defined - keeping sequences intact")
            fragments = []
            
            for seq_id, sequence in sequences.items():
                # Validate sequence first
                if not sequence or not isinstance(sequence, str):
                    logger.warning(f"Invalid sequence for {seq_id}: empty or not string")
                    continue
                
                fragment = {
                    "id": seq_id,
                    "sequence": sequence,
                    "chromosome": seq_id,
                    "start": 0,                    # 0-based start
                    "end": len(sequence),          # 0-based end (exclusive)
                    "length": len(sequence)
                }
                
                # Validate the created fragment
                validation_result = SequenceProcessor._validate_fragment_coordinates(fragment)
                if validation_result["valid"]:
                    fragments.append(fragment)
                else:
                    logger.warning(f"Created fragment for {seq_id} failed validation: {validation_result['errors']}")
                
            logger.debug(f"Created {len(fragments)} intact fragments with 0-based coordinates")
            return fragments
            
        try:
            restriction_pattern = re.compile(restriction_site, re.IGNORECASE)
            logger.debug(f"Using restriction site pattern: {restriction_site}")
        except re.error as e:
            error_msg = f"Invalid restriction site regex pattern: {restriction_site}"
            logger.error(error_msg)
            logger.debug(f"Regex compilation error: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
        
        fragments = []
        total_sites_found = 0
        sequences_with_sites = 0
        total_fragments_created = 0
        coordinate_validation_errors = 0
        
        for seq_num, (seq_id, sequence) in enumerate(sequences.items()):
            try:
                # Validate input sequence
                if not sequence or not isinstance(sequence, str):
                    logger.warning(f"Invalid sequence for {seq_id}: empty or not string")
                    continue
                
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('restriction_sequence_processing', interval=200, max_initial=3)):
                    logger.debug(f"Processing sequence {seq_num+1}/{len(sequences)}: {seq_id} ({len(sequence)} bp)")
                
                matches = list(restriction_pattern.finditer(sequence))
                sites_in_sequence = len(matches)
                total_sites_found += sites_in_sequence
                
                if not matches:
                    # No restriction sites - keep entire sequence
                    fragment = {
                        "id": seq_id,
                        "sequence": sequence,
                        "chromosome": seq_id,
                        "start": 0,                # 0-based start
                        "end": len(sequence),      # 0-based end (exclusive) 
                        "length": len(sequence)
                    }
                    
                    # Validate fragment coordinates
                    validation_result = SequenceProcessor._validate_fragment_coordinates(fragment)
                    if validation_result["valid"]:
                        fragments.append(fragment)
                        total_fragments_created += 1
                    else:
                        coordinate_validation_errors += 1
                        logger.warning(f"Fragment validation failed for {seq_id}: {validation_result['errors']}")
                    
                    if (logger.isEnabledFor(logging.DEBUG) and 
                        DebugLogLimiter.should_log('restriction_no_sites', interval=500, max_initial=2)):
                        logger.debug(f"No restriction sites found in {seq_id} - keeping entire sequence")
                    continue
                
                sequences_with_sites += 1
                
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('restriction_sites_found', interval=100, max_initial=3)):
                    logger.debug(f"Found {sites_in_sequence} restriction sites in {seq_id}")
                
                fragment_count = 0
                last_end = 0
                fragments_for_sequence = []
                
                # Process fragments between restriction sites
                for i, match in enumerate(matches):
                    # Calculate fragment boundaries in 0-based coordinates
                    frag_start = last_end         # 0-based start
                    frag_end = match.start()      # 0-based end (exclusive)
                    fragment_length = frag_end - frag_start
                    
                    if fragment_length >= Config.MIN_SEGMENT_LENGTH:
                        # Validate sequence boundaries before extraction
                        boundary_validation = SequenceProcessor._validate_sequence_boundaries(
                            sequence, frag_start, frag_end
                        )
                        
                        if boundary_validation["valid"]:
                            fragment = {
                                "id": f"{seq_id}_frag{fragment_count}",
                                "sequence": sequence[frag_start:frag_end],
                                "chromosome": seq_id,
                                "start": frag_start,      # 0-based start
                                "end": frag_end,          # 0-based end (exclusive)
                                "length": fragment_length
                            }
                            
                            # Final validation of created fragment
                            validation_result = SequenceProcessor._validate_fragment_coordinates(fragment)
                            if validation_result["valid"]:
                                fragments.append(fragment)
                                fragments_for_sequence.append(fragment)
                                fragment_count += 1
                                
                                if (logger.isEnabledFor(logging.DEBUG) and 
                                    DebugLogLimiter.should_log('restriction_fragment_creation', interval=500, max_initial=3)):
                                    logger.debug(f"Created fragment {fragment['id']}: {fragment_length} bp "
                                               f"(coordinates: [{frag_start}, {frag_end}))")
                            else:
                                coordinate_validation_errors += 1
                                logger.warning(f"Fragment validation failed for {fragment['id']}: {validation_result['errors']}")
                        else:
                            coordinate_validation_errors += 1
                            logger.warning(f"Boundary validation failed for fragment at {frag_start}-{frag_end}: {boundary_validation['errors']}")
                    
                    last_end = match.end()
                
                # Process final fragment after last restriction site
                final_start = last_end
                final_end = len(sequence)
                final_length = final_end - final_start
                
                if final_length >= Config.MIN_SEGMENT_LENGTH:
                    # Validate sequence boundaries
                    boundary_validation = SequenceProcessor._validate_sequence_boundaries(
                        sequence, final_start, final_end
                    )
                    
                    if boundary_validation["valid"]:
                        fragment = {
                            "id": f"{seq_id}_frag{fragment_count}",
                            "sequence": sequence[final_start:final_end],
                            "chromosome": seq_id,
                            "start": final_start,     # 0-based start
                            "end": final_end,         # 0-based end (exclusive)
                            "length": final_length
                        }
                        
                        # Final validation of created fragment
                        validation_result = SequenceProcessor._validate_fragment_coordinates(fragment)
                        if validation_result["valid"]:
                            fragments.append(fragment)
                            fragments_for_sequence.append(fragment)
                            fragment_count += 1
                            
                            if (logger.isEnabledFor(logging.DEBUG) and 
                                DebugLogLimiter.should_log('restriction_final_fragment', interval=500, max_initial=2)):
                                logger.debug(f"Created final fragment {fragment['id']}: {final_length} bp "
                                           f"(coordinates: [{final_start}, {final_end}))")
                        else:
                            coordinate_validation_errors += 1
                            logger.warning(f"Final fragment validation failed for {fragment['id']}: {validation_result['errors']}")
                    else:
                        coordinate_validation_errors += 1
                        logger.warning(f"Final fragment boundary validation failed for {final_start}-{final_end}: {boundary_validation['errors']}")
                
                total_fragments_created += len(fragments_for_sequence)
                
                if (logger.isEnabledFor(logging.DEBUG) and 
                    (fragment_count > 10 or DebugLogLimiter.should_log('restriction_sequence_summary', interval=50, max_initial=3))):
                    logger.debug(f"Generated {fragment_count} valid fragments from {seq_id}")
                    
            except Exception as e:
                error_msg = f"Error processing sequence {seq_id} for restriction sites"
                logger.error(error_msg)
                logger.debug(f"Processing error: {str(e)}", exc_info=True)
                continue

        logger.info(f"Generated {len(fragments)} fragments after restriction site cutting")
        logger.debug(f"=== RESTRICTION SITE CUTTING RESULTS ===")
        logger.debug(f"Input sequences: {len(sequences)}")
        logger.debug(f"Sequences with restriction sites: {sequences_with_sites}")
        logger.debug(f"Total restriction sites found: {total_sites_found}")
        logger.debug(f"Output fragments: {len(fragments)} (all with 0-based coordinates)")
        logger.debug(f"Coordinate validation errors: {coordinate_validation_errors}")
        logger.debug(f"Minimum fragment length: {Config.MIN_SEGMENT_LENGTH} bp")
        
        if logger.isEnabledFor(logging.DEBUG) and fragments:
            fragment_lengths = [frag['length'] for frag in fragments]
            avg_length = sum(fragment_lengths) / len(fragment_lengths)
            min_length = min(fragment_lengths)
            max_length = max(fragment_lengths)
            logger.debug(f"Fragment length statistics: avg={avg_length:.0f}, min={min_length}, max={max_length}")
            
            # Validate a sample of coordinates to ensure consistency
            sample_size = min(5, len(fragments))
            logger.debug(f"Sample coordinate validation (showing {sample_size} fragments):")
            for i, frag in enumerate(fragments[:sample_size]):
                expected_len = frag['end'] - frag['start']
                actual_len = len(frag['sequence'])
                coord_consistent = expected_len == actual_len == frag.get('length', -1)
                logger.debug(f"  {frag['id']}: [{frag['start']}, {frag['end']}) = {expected_len}bp, "
                           f"seq_len={actual_len}bp, consistent={coord_consistent}")
        
        logger.debug(f"=== END RESTRICTION SITE CUTTING ===")
        
        return fragments
    
    #############################################################################
    #                           Legacy Gene Overlap Support
    #############################################################################
    
    @staticmethod
    def filter_by_gene_overlap(fragments, genes, margin=None):
        """
        Filter fragments by gene overlap with coordinate-based truncation.
        
        NOTE: This method is kept for backward compatibility but should not be used
        in the new coordinate system. Use AnnotationProcessor.filter_by_gene_overlap instead.
        
        Args:
            fragments: List of fragment dictionaries with 0-based coordinates
            genes: List of gene dictionaries from GFF annotations
            margin: Gene overlap margin in base pairs, defaults to Config.GENE_OVERLAP_MARGIN
            
        Returns:
            List of filtered fragment dictionaries truncated to gene regions
            
        Raises:
            SequenceProcessingError: If gene overlap processing fails
        """
        logger.warning("filter_by_gene_overlap is deprecated. Use AnnotationProcessor.filter_by_gene_overlap instead.")
        logger.debug(f"=== LEGACY GENE OVERLAP FILTERING ===")
        logger.debug(f"Filtering {len(fragments)} fragments by gene overlap")
        
        if margin is None:
            margin = Config.GENE_OVERLAP_MARGIN
        
        logger.debug(f"Using gene overlap margin: {margin} bp")
        
        if not genes:
            logger.warning("No gene annotations provided for filtering - returning original fragments")
            return fragments
            
        # Import the proper implementation
        try:
            from . import AnnotationProcessor
            return AnnotationProcessor.filter_by_gene_overlap(fragments, genes)
        except ImportError:
            logger.error("Cannot import AnnotationProcessor for gene overlap filtering")
            return fragments