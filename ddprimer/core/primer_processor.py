#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer processing module for ddPrimer pipeline.

Handles primer record creation and amplicon extraction operations.

Contains functionality for:
1. Primer record creation from Primer3 output
2. Amplicon sequence extraction and validation  
3. Coordinate system handling and validation

COORDINATE SYSTEM:
- INPUT: Primer3 coordinates (0-based, from primer3-py bindings)
- INTERNAL: 0-based coordinates for all sequence operations
- OUTPUT: Primer records with coordinate metadata (mixed for compatibility)
- Amplicon extraction uses 0-based string indexing
- Fragment coordinates are passed through from AnnotationProcessor/SequenceProcessor
"""

import logging
import pandas as pd
from typing import Dict, List, Optional, Any

# Import package modules
from ..config import Config, SequenceProcessingError, PrimerDesignError, DebugLogLimiter
from .filter_processor import FilterProcessor

# Set up module logger
logger = logging.getLogger(__name__)


class PrimerProcessor:
    """
    Handles primer record creation and amplicon extraction operations ONLY.
    
    This class provides primer processing capabilities but does NOT handle
    filtering - that responsibility belongs to FilterProcessor.
    
    Single Responsibility: Create and validate primer records and amplicons.
    """
    
    #############################################################################
    #                           Workflow Methods
    #############################################################################
    
    @classmethod
    def create_primer_records_workflow(cls, primer3_results: List[Dict]) -> List[Dict]:
        """
        Create primer records from Primer3 results for workflow integration.
        
        This workflow wrapper transforms raw Primer3 output into standardized primer
        records containing all necessary information for downstream filtering and
        analysis. It handles coordinate validation, amplicon extraction, and sequence
        matching while maintaining clean separation from filtering operations.
        
        Args:
            primer3_results: List of Primer3 result dictionaries containing fragment
                           data and primer pair information from the design process
            
        Returns:
            List of primer record dictionaries (unfiltered) ready for downstream
            processing, including amplicon sequences and validation flags
            
        Raises:
            PrimerDesignError: If primer record creation fails during workflow coordination
        """
        logger.debug("=== WORKFLOW: PRIMER RECORD CREATION ===")
        logger.debug(f"Creating primer records from {len(primer3_results)} Primer3 results")
        
        try:
            primer_records = []
            processor = cls()
            
            for result in primer3_results:
                # Extract fragment and primer pairs from result
                fragment = result.get("fragment")
                primer_pairs = result.get("primer_pairs", [])
                
                if not fragment or not primer_pairs:
                    continue
                
                for pair_index, pair_data in enumerate(primer_pairs):
                    try:
                        primer_record = processor.create_primer_record(
                            fragment, pair_data, pair_index
                        )
                        if primer_record:
                            primer_records.append(primer_record)
                    except Exception as e:
                        if DebugLogLimiter.should_log('primer_record_creation_errors', interval=100, max_initial=2):
                            logger.debug(f"Error creating primer record for pair {pair_index}: {e}")
                        continue
            
            logger.debug(f"Created {len(primer_records)} primer records")
            logger.debug("=== END WORKFLOW: PRIMER RECORD CREATION ===")
            
            return primer_records
            
        except Exception as e:
            error_msg = f"Error in primer record creation workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: PRIMER RECORD CREATION ===")
            raise PrimerDesignError(error_msg) from e
    
    @classmethod
    def filter_primers_workflow(cls, primer_records: List[Dict]) -> Optional['pd.DataFrame']:
        """
        Delegate primer filtering to FilterProcessor to maintain separation of concerns.
        
        This workflow wrapper provides a clean delegation interface to the FilterProcessor
        while maintaining the architectural separation between primer record creation
        and filtering operations. It serves as the workflow coordination point for
        the filtering phase of the pipeline.
        
        Args:
            primer_records: List of primer record dictionaries containing sequences,
                          penalties, amplicons, and validation information
            
        Returns:
            Filtered DataFrame from FilterProcessor containing only primers that
            meet all filtering criteria, or None if no primers pass
            
        Raises:
            PrimerDesignError: If filtering workflow coordination fails
        """
        logger.debug("=== WORKFLOW: DELEGATING TO FILTER PROCESSOR ===")
        
        return FilterProcessor.filter_primers_workflow(primer_records)
    
    #############################################################################
    
    def __init__(self):
        """Initialize the primer processor."""
        pass
    
    def create_primer_record(self, fragment: Dict, pair_data: Dict, pair_index: int) -> Optional[Dict]:
        """
        Create a primer record from fragment and primer pair data, calculating
        correct genomic coordinates for the amplicon.

        Args:
            fragment: Fragment dictionary containing template sequence and metadata,
                        including its 0-based genomic start coordinate in the 'start' field.
            pair_data: Primer pair data from Primer3.
            pair_index: Index of the primer pair.

        Returns:
            Complete primer record dictionary or None if creation fails.
        """
        if DebugLogLimiter.should_log('primer_record_creation', interval=200, max_initial=2):
            logger.debug(f"Creating primer record for pair {pair_index}")

        try:
            forward_primer = pair_data.get("left_sequence")
            reverse_primer = pair_data.get("right_sequence")

            if not forward_primer or not reverse_primer:
                if DebugLogLimiter.should_log('primer_record_missing_sequences', interval=100, max_initial=1):
                    logger.debug(f"Missing primer sequences for pair {pair_index}")
                return None

            # Get primer coordinates relative to the fragment (0-based)
            forward_start_relative = pair_data.get("left_start", 0)
            forward_length = len(forward_primer)
            # The 'right_start' from primer3 is the 0-based index of the 3' (last) base
            reverse_end_relative = pair_data.get("right_start", 0)
            reverse_length = len(reverse_primer)

            # Get Tm values
            forward_tm = pair_data.get("left_tm", 0.0)
            reverse_tm = pair_data.get("right_tm", 0.0)

            # Process probe
            probe_seq = self._process_probe_sequence_safe(pair_data, pair_index)
            probe_tm = pair_data.get("internal_tm", 0.0) if probe_seq else 0.0
            probe_penalty = pair_data.get("internal_penalty", 0.0) if probe_seq else 0.0  # ✅ ADDED THIS LINE

            # --- CORRECTED GENOMIC COORDINATE CALCULATION ---
            # Get the TRUE 0-based genomic start position of this fragment.
            # This is now reliably stored in the 'start' field.
            fragment_genomic_start = fragment.get("start", 0)
            fragment_chromosome = fragment.get("chromosome", "")
            
            # Calculate the TRUE 0-based genomic coordinates of the amplicon.
            # Add the relative primer positions to the fragment's true genomic start.
            # We add 1 to the start to convert from 0-based to 1-based for the final output.
            amplicon_genomic_start = fragment_genomic_start + forward_start_relative + 1
            # The end coordinate is the fragment's start + the relative end of the reverse primer.
            amplicon_genomic_end = fragment_genomic_start + reverse_end_relative + 1

            # Debugging coordinate calculation
            if DebugLogLimiter.should_log('genomic_coordinate_debug', interval=500, max_initial=3):
                logger.debug(f"Coordinate calculation for {fragment.get('id', 'unknown')}:")
                logger.debug(f"  Fragment genomic start (0-based): {fragment_genomic_start}")
                logger.debug(f"  Forward relative start (0-based): {forward_start_relative}")
                logger.debug(f"  Calculated amplicon genomic start (1-based): {amplicon_genomic_start}")

            # Extract and validate amplicon
            amplicon_result = self._extract_and_validate_amplicon(
                fragment, forward_start_relative, forward_length, reverse_end_relative, reverse_length
            )
            
            # Create primer record
            primer_record = {
                "Fragment ID": fragment.get("id", "unknown"),
                "Fragment": fragment.get("sequence", ""),
                "Sequence (F)": forward_primer,
                "Sequence (R)": reverse_primer,
                "Tm (F)": forward_tm,
                "Tm (R)": reverse_tm,
                "Penalty": pair_data.get("penalty", 0.0),
                "Penalty (F)": pair_data.get("left_penalty", 0.0),
                "Penalty (R)": pair_data.get("right_penalty", 0.0),
                "Chr": fragment_chromosome,
                "Start": amplicon_genomic_start, # Correct genomic start (1-based)
                "End": amplicon_genomic_end,     # Correct genomic end (1-based)
                "Gene": fragment.get("Gene", fragment.get("gene", "")),
            }

            # Add probe if available
            if probe_seq:
                primer_record["Sequence (P)"] = probe_seq
                primer_record["Tm (P)"] = probe_tm
                primer_record["Penalty (P)"] = probe_penalty  # ✅ ADDED THIS LINE

            # Add amplicon
            primer_record["Sequence (A)"] = amplicon_result["amplicon"]
            primer_record["Length"] = amplicon_result["amplicon_length"]
            
            return primer_record

        except Exception as e:
            if DebugLogLimiter.should_log('primer_record_creation_errors', interval=50, max_initial=1):
                logger.error(f"Error creating primer record for pair {pair_index}: {str(e)}")
            return None
    
    def _validate_primer3_sequences(self, primer3_forward: str, primer3_reverse: str, 
                                   extracted_forward: str, extracted_reverse: str) -> Dict:
        """
        Validate that Primer3 sequences match what we extract from the template.
        
        Args:
            primer3_forward: Forward sequence from Primer3
            primer3_reverse: Reverse sequence from Primer3  
            extracted_forward: Forward sequence extracted from template
            extracted_reverse: Reverse sequence extracted from template
            
        Returns:
            Dictionary with validation results
        """
        errors = []
        
        if primer3_forward.upper() != extracted_forward.upper():
            errors.append(f"Forward primer mismatch")
        
        if primer3_reverse.upper() != extracted_reverse.upper():
            errors.append(f"Reverse primer mismatch")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors
        }
    
    def _process_probe_sequence_safe(self, pair_data: Dict, pair_index: int) -> Optional[str]:
        """
        Safely process probe sequence with comprehensive error handling.
        
        Args:
            pair_data: Dictionary containing primer pair information
            pair_index: Index of the current primer pair
            
        Returns:
            Processed probe sequence or None if not available
        """
        try:
            # Get probe sequence from pair data
            probe_seq = pair_data.get("internal_sequence")
            
            if not probe_seq or probe_seq == "":
                return None
            
            # Try to apply C/G filtering with comprehensive error handling
            try:
                result = FilterProcessor.ensure_more_c_than_g(probe_seq)
                
                if result is None:
                    return probe_seq
                
                if isinstance(result, tuple):
                    if len(result) == 2:
                        processed_seq, was_reversed = result
                        if was_reversed and DebugLogLimiter.should_log('probe_reversed', interval=20, max_initial=2):
                            logger.debug(f"Probe sequence was reversed for pair {pair_index}")
                        return processed_seq
                    else:
                        return probe_seq
                
                if isinstance(result, str):
                    return result
                
                return probe_seq
                
            except Exception as filter_error:
                if DebugLogLimiter.should_log('probe_filter_error', interval=100, max_initial=1):
                    logger.debug(f"FilterProcessor error for probe: {str(filter_error)}")
                return probe_seq
            
        except Exception as e:
            if DebugLogLimiter.should_log('probe_processing_error', interval=100, max_initial=1):
                logger.warning(f"Failed to process probe sequence for pair {pair_index}: {str(e)}")
            return None
        
    
    def _extract_and_validate_amplicon(self, fragment: Dict, forward_start: int, forward_length: int,
                                     reverse_start: int, reverse_length: int) -> Dict:
        """
        Extract amplicon and provide validation information.
        
        Args:
            fragment: Fragment containing template sequence
            forward_start: Forward primer start position (0-based)
            forward_length: Forward primer length
            reverse_start: Reverse primer start position (0-based)
            reverse_length: Reverse primer length
            
        Returns:
            Dictionary with amplicon and validation information
        """
        if DebugLogLimiter.should_log('amplicon_extraction', interval=500, max_initial=2):
            logger.debug(f"Amplicon extraction: forward={forward_start},{forward_length}, reverse={reverse_start},{reverse_length}")
            logger.debug(f"Template length: {len(fragment.get('sequence', ''))} bp")
        
        amplicon_result = self.extract_amplicon(
            fragment, forward_start, forward_length, reverse_start, reverse_length
        )
        
        if amplicon_result["validated"]:
            if DebugLogLimiter.should_log('amplicon_success', interval=500, max_initial=1):
                logger.debug(f"Amplicon extraction successful: {amplicon_result['amplicon_length']} bp")
        else:
            if DebugLogLimiter.should_log('amplicon_failure', interval=200, max_initial=2):
                logger.debug(f"Amplicon extraction failed: {amplicon_result.get('validation_errors', [])}")
        
        return amplicon_result
    
    def extract_amplicon(self, fragment: Dict, forward_start: int, forward_length: int,
                        reverse_start: int, reverse_length: int) -> Dict:
        """
        Extract amplicon sequence from template with comprehensive validation.

        - LEFT primer: start position is the FIRST base (0-based)
        - RIGHT primer: start position is the LAST base (0-based)
        """
        template = fragment.get("sequence", "")
        template_length = len(template)
        
        # Calculate actual reverse primer start position
        # Primer3 reports the LAST base position, we need the FIRST base position
        reverse_actual_start = reverse_start - reverse_length + 1
        
        # Initialize result dictionary
        result = {
            "amplicon": "",
            "amplicon_length": 0,
            "validated": False,
            "validation_errors": [],
            "forward_primer_extracted": "",
            "reverse_primer_extracted": "",
            "coordinates_valid": False,
            "primers_match": False
        }
        
        # Validate coordinates
        validation_errors = []
        
        # Check forward primer bounds
        if forward_start < 0:
            validation_errors.append(f"Forward primer start ({forward_start}) is negative")
        if forward_start >= template_length:
            validation_errors.append(f"Forward primer start ({forward_start}) exceeds template length ({template_length})")
        if forward_start + forward_length > template_length:
            validation_errors.append(f"Forward primer end ({forward_start + forward_length}) exceeds template length ({template_length})")
        
        # Check reverse primer bounds
        if reverse_actual_start < 0:
            validation_errors.append(f"Reverse primer actual start ({reverse_actual_start}) is negative")
        if reverse_start >= template_length:
            validation_errors.append(f"Reverse primer end ({reverse_start}) exceeds template length ({template_length})")
        if reverse_actual_start >= template_length:
            validation_errors.append(f"Reverse primer actual start ({reverse_actual_start}) exceeds template length ({template_length})")
        
        # Check primer order
        forward_end = forward_start + forward_length
        if forward_start >= reverse_actual_start:
            validation_errors.append(f"Forward primer start ({forward_start}) must be before reverse primer start ({reverse_actual_start})")
        if forward_end > reverse_actual_start:
            validation_errors.append(f"Forward primer end ({forward_end}) overlaps with reverse primer start ({reverse_actual_start})")
        
        result["validation_errors"] = validation_errors
        result["coordinates_valid"] = len(validation_errors) == 0
        
        if validation_errors:
            if DebugLogLimiter.should_log('amplicon_coordinate_errors', interval=100, max_initial=1):
                logger.debug(f"Coordinate validation failed: {validation_errors}")
            return result
        
        # Extract primers for validation using correct coordinates
        try:
            # Forward primer: standard extraction
            forward_extracted = template[forward_start:forward_start + forward_length]
            
            # Reverse primer: extract from template and reverse complement
            reverse_template_sequence = template[reverse_actual_start:reverse_actual_start + reverse_length]
            reverse_extracted = FilterProcessor.reverse_complement(reverse_template_sequence)
            
            result["forward_primer_extracted"] = forward_extracted
            result["reverse_primer_extracted"] = reverse_extracted
            
        except Exception as e:
            error_msg = f"Error extracting primers: {str(e)}"
            validation_errors.append(error_msg)
            result["validation_errors"] = validation_errors
            return result
        
        # Extract amplicon (from start of forward primer to end of reverse primer)
        try:
            amplicon_start = forward_start
            amplicon_end = reverse_start + 1  # +1 because reverse_start is the last base (inclusive)
            amplicon_sequence = template[amplicon_start:amplicon_end]
            
            result["amplicon"] = amplicon_sequence
            result["amplicon_length"] = len(amplicon_sequence)
            
            # Validate primer-amplicon consistency
            primer_validation_result = self._validate_primer_amplicon_match(
                amplicon_sequence, forward_extracted, reverse_extracted
            )
            
            result["primers_match"] = primer_validation_result["valid"]
            if not primer_validation_result["valid"]:
                validation_errors.extend(primer_validation_result["errors"])
                result["validation_errors"] = validation_errors
                if DebugLogLimiter.should_log('primer_amplicon_mismatch', interval=100, max_initial=1):
                    logger.warning(f"Primer-amplicon validation failed")
            
            # Overall validation passes if coordinates are valid AND primers match
            result["validated"] = result["coordinates_valid"] and result["primers_match"]
            
        except Exception as e:
            error_msg = f"Error extracting amplicon: {str(e)}"
            validation_errors.append(error_msg)
            result["validation_errors"] = validation_errors
        
        return result
    
    def _validate_primer_amplicon_match(self, amplicon: str, forward_primer: str, reverse_primer: str) -> Dict:
        """
        Validate that primers properly match the amplicon sequence.
        
        Args:
            amplicon: Full amplicon sequence
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            
        Returns:
            Dictionary with validation results
        """
        errors = []
        
        # Check if amplicon starts with forward primer
        if not amplicon.startswith(forward_primer):
            amplicon_start = amplicon[:len(forward_primer)] if len(amplicon) >= len(forward_primer) else amplicon
            errors.append(f"Amplicon start '{amplicon_start}' does not match forward primer '{forward_primer}'")
        
        # Check if amplicon ends with reverse complement of reverse primer
        try:
            reverse_primer_rc = FilterProcessor.reverse_complement(reverse_primer)
            
            if not amplicon.endswith(reverse_primer_rc):
                amplicon_end = amplicon[-len(reverse_primer_rc):] if len(amplicon) >= len(reverse_primer_rc) else amplicon
                errors.append(f"Amplicon end '{amplicon_end}' does not match reverse primer RC '{reverse_primer_rc}'")
                
        except Exception as e:
            errors.append(f"Failed to generate reverse complement of reverse primer: {str(e)}")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors
        }