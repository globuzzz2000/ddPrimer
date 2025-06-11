#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer3 processing module for ddPrimer pipeline.

Handles primer design through Primer3 including:
1. Running primer3 to design primers and probes using Python bindings
2. Parsing primer3 output with comprehensive error handling
3. Amplicon extraction and validation with coordinate correction
4. Parallel processing for improved performance
5. K-mer template masking integration using Primer3's masking functionality
6. Result formatting and validation
"""

import logging
import re
import os
import primer3
import multiprocessing
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from pathlib import Path

from ..config import Config, SequenceProcessingError, PrimerDesignError
from ..utils import SequenceUtils, Selector

# Set up module logger
logger = logging.getLogger(__name__)


class Primer3Processor:
    """
    Handles Primer3 operations for primer design using the primer3-py package.
    
    This class provides comprehensive primer design capabilities including
    parallel processing, result parsing, amplicon extraction, and k-mer
    template masking integration with proper coordinate handling for various pipeline modes.
    
    Attributes:
        config: Configuration object with primer3 settings
        kmer_library_path: Path to k-mer library file for template masking
        
    Example:
        >>> processor = Primer3Processor()
        >>> results = processor.run_primer3_batch_parallel(input_blocks)
        >>> primer_records = processor.parse_primer3_batch(results, fragment_info)
    """
    
    def __init__(self, config=None):
        """
        Initialize Primer3Processor with configuration and k-mer template masking setup.
        """
        self.config = config if config is not None else Config
        
        # Initialize k-mer masking attributes
        self.kmer_folder_path = None
        self.kmer_file_prefix = None
        
        # Initialize k-mer template masking if enabled
        self._setup_kmer_masking()
        
        logger.debug("Initialized Primer3Processor with k-mer template masking integration")

    def _setup_kmer_masking(self):
        """
        Set up k-mer masking parameters for Primer3.
        
        Loads folder path and file prefix from saved user settings.
        """
        logger.debug("=== K-MER MASKING SETUP DEBUG ===")
        
        try:
            # Check if k-mer masking is enabled in config
            if not getattr(self.config, 'KMER_MASKING_ENABLED', True):
                logger.debug("K-mer masking disabled in configuration")
                return
            
            # Load k-mer selection from user settings
            from ..utils import Selector
            kmer_config = Selector.load_kmer_selection()
            
            if kmer_config:
                # Use saved folder path and file prefix
                self.kmer_folder_path = kmer_config.get("folder_path")
                self.kmer_file_prefix = kmer_config.get("file_prefix")
                organism = kmer_config.get("organism", "unknown")
                
                if self.kmer_folder_path and self.kmer_file_prefix and os.path.exists(self.kmer_folder_path):
                    logger.debug(f"Using saved k-mer selection: {organism}")
                    logger.debug(f"K-mer folder: {self.kmer_folder_path}")
                    logger.debug(f"K-mer prefix: {self.kmer_file_prefix}")
                    logger.info(f"K-mer template masking enabled for: {organism}")
                    return
                else:
                    logger.debug(f"Saved k-mer folder path invalid or missing: {self.kmer_folder_path}")
            
            # Try to auto-select from available k-mer lists
            kmer_lists_by_organism = Selector.find_kmer_lists()
            
            if kmer_lists_by_organism:
                folder_path, file_prefix, organism = self._auto_select_kmer_settings(kmer_lists_by_organism)
                
                if folder_path and file_prefix:
                    self.kmer_folder_path = folder_path
                    self.kmer_file_prefix = file_prefix
                    logger.debug(f"Auto-selected k-mer folder: {folder_path}")
                    logger.debug(f"Auto-selected k-mer prefix: {file_prefix}")
                    logger.info(f"K-mer template masking enabled using: {organism} ({file_prefix})")
                    return
            
            # No k-mer libraries available
            logger.debug("No k-mer libraries available - primer design will proceed without k-mer masking")
            logger.info("K-mer masking disabled - no k-mer libraries found")
            
        except Exception as e:
            error_msg = f"Error setting up k-mer masking: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            # Don't raise exception - k-mer masking is optional
            self.kmer_folder_path = None
            self.kmer_file_prefix = None
            logger.warning("K-mer masking will be disabled due to setup error")
        
        logger.debug("=== END K-MER MASKING SETUP DEBUG ===")

    def _auto_select_kmer_settings(self, kmer_lists_by_organism: Dict) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Automatically select k-mer folder, prefix, and organism from available options.
        
        Args:
            kmer_lists_by_organism: Dictionary of k-mer lists grouped by organism
            
        Returns:
            Tuple of (folder_path, file_prefix, organism_name) or (None, None, None)
        """
        logger.debug("Auto-selecting k-mer settings")
        
        # Priority-based organism selection
        priority_organisms = [
            'arabidopsis thaliana', 'homo sapiens', 'mus musculus',
            'drosophila melanogaster', 'saccharomyces cerevisiae',
            'escherichia coli', 'caenorhabditis elegans'
        ]
        
        selected_organism = None
        for priority_org in priority_organisms:
            for organism in kmer_lists_by_organism.keys():
                if priority_org in organism.lower():
                    selected_organism = organism
                    break
            if selected_organism:
                break
        
        # If no priority organism found, use the first available
        if not selected_organism:
            selected_organism = list(kmer_lists_by_organism.keys())[0]
        
        organism_kmers = kmer_lists_by_organism[selected_organism]
        logger.debug(f"Auto-selected organism: {selected_organism}")
        
        if organism_kmers:
            # Get folder and prefix from first file
            first_file_info = organism_kmers[0]
            first_file_path = first_file_info['path']
            
            folder_path = str(Path(first_file_path).parent)
            
            # Extract prefix from filename (e.g., "thaliana_11.list" -> "thaliana")
            filename = Path(first_file_path).stem
            prefix_parts = filename.split('_')[:-1]  # Remove the size part
            file_prefix = '_'.join(prefix_parts)
            
            logger.debug(f"Selected folder: {folder_path}")
            logger.debug(f"Selected prefix: {file_prefix}")
            
            return folder_path, file_prefix, selected_organism
        
        return None, None, None

    def get_primer3_global_args_with_kmer(self) -> Dict:
        """
        Get Primer3 global arguments with k-mer template masking integration.
        
        Uses the saved folder path and file prefix for k-mer masking.
        
        Returns:
            Dictionary of Primer3 arguments including k-mer template masking settings
        """
        # Get base arguments
        args = self.config.get_primer3_global_args()
        
        # Early exit if k-mer masking is disabled
        if not getattr(self.config, 'KMER_MASKING_ENABLED', True):
            logger.info("K-mer masking is DISABLED - using base Primer3 settings only")
            return args
        
        # Add k-mer template masking if folder and prefix are available
        if (self.kmer_folder_path and self.kmer_file_prefix and 
            os.path.exists(self.kmer_folder_path)):
            
            logger.debug(f"Adding k-mer template masking to Primer3 args:")
            logger.debug(f"  Folder: {self.kmer_folder_path}")
            logger.debug(f"  Prefix: {self.kmer_file_prefix}")
            
            # Primer3 k-mer masking parameters
            template_masking_args = {
                'PRIMER_MASK_TEMPLATE': 1,
                'PRIMER_MASK_KMERLIST_PATH': self.kmer_folder_path,
                'PRIMER_MASK_KMERLIST_PREFIX': self.kmer_file_prefix,
            }
            
            # Debug: Print template masking arguments being added
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("Adding template masking arguments:")
                for key, value in template_masking_args.items():
                    logger.debug(f"  {key} = {value}")
            
            args.update(template_masking_args)
            logger.info(f"Template masking enabled: {Path(self.kmer_folder_path).name}/{self.kmer_file_prefix}_*")
        else:
            logger.debug("No k-mer folder/prefix available - using standard Primer3 arguments")
        
        return args
    
    def run_primer3_batch(self, input_blocks):
        """
        Run primer3 on a batch of input blocks using Python bindings.
        
        Processes multiple primer design requests sequentially using the
        primer3-py package, formatting results to maintain compatibility
        with existing parsing code. Now includes k-mer template masking if available.
        
        Args:
            input_blocks: List of dictionaries containing primer3 parameters
            
        Returns:
            Combined primer3 output string in primer3_core format
            
        Raises:
            PrimerDesignError: If primer3 execution fails
            
        Example:
            >>> input_blocks = [{"SEQUENCE_ID": "seq1", "SEQUENCE_TEMPLATE": "ATCG"}]
            >>> output = processor.run_primer3_batch(input_blocks)
        """
        logger.debug(f"Running Primer3 on {len(input_blocks)} input blocks")
        
        if not input_blocks:
            logger.warning("No input blocks provided for Primer3 processing")
            return ""
        
        # Log k-mer masking status
        if self._is_kmer_masking_available():
            logger.debug(f"K-mer template masking enabled with folder: {self.kmer_folder_path}")
        else:
            logger.debug("K-mer template masking disabled")
        
        results = []
        
        try:
            # Get global arguments with k-mer template masking settings
            global_args = self.get_primer3_global_args_with_kmer()
            
            # Process each input block using primer3-py
            for block in input_blocks:
                sequence_id = block.get("SEQUENCE_ID", "UNKNOWN")
                
                if logger.isEnabledFor(logging.DEBUG):
                    template_length = len(block.get("SEQUENCE_TEMPLATE", ""))
                    logger.debug(f"Processing sequence {sequence_id} with {template_length} bp template")
                
                # Run primer3 design
                try:
                    primer_result = primer3.bindings.design_primers(
                        seq_args=block,
                        global_args=global_args
                    )
                    
                    # Format result to match primer3_core output
                    formatted_result = self._format_primer3_result(sequence_id, block, primer_result)
                    results.append(formatted_result)
                    
                except Exception as e:
                    error_msg = f"Primer3 design failed for sequence {sequence_id}"
                    logger.error(error_msg)
                    logger.debug(f"Primer3 error details: {str(e)}", exc_info=True)
                    # Add empty result to maintain sequence order
                    results.append(f"SEQUENCE_ID={sequence_id}\nPRIMER_PAIR_NUM_RETURNED=0\n=")
            
            logger.debug(f"Successfully processed {len(results)} primer3 blocks")
            return "\n".join(results)
            
        except Exception as e:
            error_msg = f"Batch primer3 processing failed"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise PrimerDesignError(error_msg) from e
    
    def run_primer3_batch_parallel(self, input_blocks, max_workers=None):
        """
        Run primer3 on input blocks using parallel processing.
        
        Distributes primer design tasks across multiple processes for
        improved performance, with progress tracking and error handling.
        Now includes k-mer template masking support.
        
        Args:
            input_blocks: List of dictionaries containing primer3 parameters
            max_workers: Maximum number of worker processes, defaults to CPU count
            
        Returns:
            Combined primer3 output string
            
        Raises:
            PrimerDesignError: If parallel processing setup or execution fails
            
        Example:
            >>> output = processor.run_primer3_batch_parallel(blocks, max_workers=4)
        """
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), len(input_blocks))
        
        logger.debug(f"Running Primer3 in parallel: {max_workers} workers for {len(input_blocks)} blocks")
        
        # Log k-mer masking status - FIXED: use correct attributes
        if self._is_kmer_masking_available():
            logger.info(f"K-mer template masking enabled using: {self.kmer_file_prefix}")
        
        # Log input statistics in debug mode
        if logger.isEnabledFor(logging.DEBUG):
            self._log_input_statistics(input_blocks)
        
        try:
            # Split input blocks into chunks for parallel processing
            chunk_size = max(1, len(input_blocks) // max_workers)
            chunks = [input_blocks[i:i + chunk_size] for i in range(0, len(input_blocks), chunk_size)]
            
            logger.debug(f"Split input into {len(chunks)} chunks of approximately {chunk_size} blocks each")
            
            # Process chunks in parallel
            results = []
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(self.run_primer3_batch, chunk) for chunk in chunks]
                
                # Add progress tracking if enabled
                if self.config.SHOW_PROGRESS:
                    for future in tqdm(futures, total=len(futures), desc="Running Primer3"):
                        try:
                            result = future.result()
                            results.append(result)
                        except Exception as e:
                            logger.error(f"Primer3 chunk processing failed: {str(e)}")
                            results.append("")  # Add empty result to maintain order
                else:
                    for future in futures:
                        try:
                            result = future.result()
                            results.append(result)
                        except Exception as e:
                            logger.error(f"Primer3 chunk processing failed: {str(e)}")
                            results.append("")
            
            combined_output = "\n".join(filter(None, results))  # Filter out empty results
            
            # Analyze output in debug mode
            if logger.isEnabledFor(logging.DEBUG):
                self._log_output_statistics(combined_output, len(input_blocks))
            
            logger.debug("Completed parallel Primer3 processing")
            return combined_output
            
        except Exception as e:
            error_msg = f"Parallel Primer3 processing failed"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise PrimerDesignError(error_msg) from e
        
    def _is_kmer_masking_available(self):
        """
        Check if k-mer masking is available and properly configured.
        
        Returns:
            bool: True if k-mer masking is available
        """
        return (getattr(self.config, 'KMER_MASKING_ENABLED', True) and 
                self.kmer_folder_path and 
                self.kmer_file_prefix and 
                os.path.exists(self.kmer_folder_path))
    
    def _log_input_statistics(self, input_blocks):
        """
        Log detailed statistics about input blocks for debugging.
        
        Args:
            input_blocks: List of input blocks to analyze
        """
        seq_lengths = [len(block.get('SEQUENCE_TEMPLATE', '')) for block in input_blocks]
        if seq_lengths:
            avg_len = sum(seq_lengths) / len(seq_lengths)
            min_len = min(seq_lengths)
            max_len = max(seq_lengths)
            logger.debug(f"Input statistics: {len(input_blocks)} blocks")
            logger.debug(f"  Template lengths: avg={avg_len:.0f}, min={min_len}, max={max_len}")
            
            # Length distribution analysis
            short_count = sum(1 for length in seq_lengths if length < 150)
            medium_count = sum(1 for length in seq_lengths if 150 <= length < 500)
            long_count = sum(1 for length in seq_lengths if length >= 500)
            logger.debug(f"  Length distribution: <150bp={short_count}, 150-500bp={medium_count}, ≥500bp={long_count}")
            
            # K-mer masking status - FIXED: use helper method
            if self._is_kmer_masking_available():
                logger.debug(f"  K-mer template masking: ENABLED ({self.kmer_file_prefix})")
            else:
                logger.debug(f"  K-mer template masking: DISABLED")

    def _log_output_statistics(self, output, input_count):
        """
        Analyze and log Primer3 output statistics for debugging.
        
        Args:
            output: Combined Primer3 output string
            input_count: Number of input blocks processed
        """
        logger.debug("=== PRIMER3 OUTPUT ANALYSIS ===")
        
        # Count sequence blocks in output
        sequence_blocks = output.split("SEQUENCE_ID=")
        num_blocks = len(sequence_blocks) - 1  # First split is empty
        logger.debug(f"Primer3 output contains {num_blocks} blocks from {input_count} input blocks")
        
        if num_blocks < input_count:
            logger.warning(f"Lost {input_count - num_blocks} blocks during Primer3 processing")
        
        # Analyze success and failure rates
        successful_designs = 0
        failed_designs = 0
        error_types = {}
        
        for block in sequence_blocks[1:]:  # Skip first empty block
            if "PRIMER_PAIR_NUM_RETURNED=" in block:
                match = re.search(r'PRIMER_PAIR_NUM_RETURNED=(\d+)', block)
                if match:
                    pairs_returned = int(match.group(1))
                    if pairs_returned > 0:
                        successful_designs += 1
                    else:
                        failed_designs += 1
                        
                        # Identify failure reasons
                        if "PRIMER_ERROR=" in block:
                            error_match = re.search(r'PRIMER_ERROR=(.*?)[\n=]', block)
                            if error_match:
                                error_msg = error_match.group(1).strip()
                                error_types[error_msg] = error_types.get(error_msg, 0) + 1
        
        logger.debug(f"Design results: {successful_designs} successful, {failed_designs} failed")
        
        if error_types:
            logger.debug("Failure reasons:")
            for error, count in error_types.items():
                logger.debug(f"  {error}: {count}")
        
        # Count total primer pairs found
        total_pairs = output.count("PRIMER_PAIR_0_PENALTY=")
        logger.debug(f"Total primer pairs generated: {total_pairs}")
        
        # Log k-mer masking status - FIXED: use helper method
        if self._is_kmer_masking_available():
            logger.debug(f"  K-mer template masking: ACTIVE ({self.kmer_file_prefix})")
        else:
            logger.debug(f"  K-mer template masking: INACTIVE")
            
        logger.debug("=== END PRIMER3 OUTPUT ANALYSIS ===")
    
    def _format_primer3_result(self, sequence_id, input_block, primer_result):
        """
        Format primer3-py result to match primer3_core output format.
        
        Converts the Python dictionary output from primer3-py into the
        text format expected by the existing parsing infrastructure.
        
        Args:
            sequence_id: ID of the processed sequence
            input_block: Original input parameters
            primer_result: Results from primer3.bindings.design_primers
            
        Returns:
            Formatted output string matching primer3_core format
            
        Example:
            >>> formatted = processor._format_primer3_result("seq1", input_dict, result_dict)
        """
        lines = [f"SEQUENCE_ID={sequence_id}"]
        
        # Add template sequence
        if "SEQUENCE_TEMPLATE" in input_block:
            lines.append(f"SEQUENCE_TEMPLATE={input_block['SEQUENCE_TEMPLATE']}")
        
        # Count primer pairs found
        num_pairs = 0
        while f"PRIMER_LEFT_{num_pairs}_SEQUENCE" in primer_result:
            num_pairs += 1
        
        lines.append(f"PRIMER_PAIR_NUM_RETURNED={num_pairs}")
        
        # Add any error messages
        if "PRIMER_ERROR" in primer_result:
            lines.append(f"PRIMER_ERROR={primer_result['PRIMER_ERROR']}")
        
        # Format each primer pair
        for i in range(num_pairs):
            self._format_primer_pair(lines, primer_result, i)
        
        # Add end marker
        lines.append("=")
        
        return "\n".join(lines)
    
    def _format_primer_pair(self, lines, primer_result, pair_index):
        """
        Format a single primer pair for output.
        
        Args:
            lines: List to append formatted lines to
            primer_result: Primer3 result dictionary
            pair_index: Index of the primer pair to format
        """
        i = pair_index
        
        # Penalty and product size
        if f"PRIMER_PAIR_{i}_PENALTY" in primer_result:
            lines.append(f"PRIMER_PAIR_{i}_PENALTY={primer_result[f'PRIMER_PAIR_{i}_PENALTY']}")
        
        if f"PRIMER_PAIR_{i}_PRODUCT_SIZE" in primer_result:
            lines.append(f"PRIMER_PAIR_{i}_PRODUCT_SIZE={primer_result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}")
        
        # Format primer components
        for side in ["LEFT", "RIGHT", "INTERNAL"]:
            self._format_primer_component(lines, primer_result, i, side)
    
    def _format_primer_component(self, lines, primer_result, pair_index, side):
        """
        Format a single primer component (left, right, or internal).
        
        Args:
            lines: List to append formatted lines to
            primer_result: Primer3 result dictionary
            pair_index: Index of the primer pair
            side: Component side ("LEFT", "RIGHT", or "INTERNAL")
        """
        i = pair_index
        
        # Sequence
        if f"PRIMER_{side}_{i}_SEQUENCE" in primer_result:
            lines.append(f"PRIMER_{side}_{i}_SEQUENCE={primer_result[f'PRIMER_{side}_{i}_SEQUENCE']}")
        
        # Position and length
        if f"PRIMER_{side}_{i}" in primer_result:
            start, length = primer_result[f"PRIMER_{side}_{i}"]
            lines.append(f"PRIMER_{side}_{i}={start},{length}")
        
        # Temperature and penalty
        if f"PRIMER_{side}_{i}_TM" in primer_result:
            lines.append(f"PRIMER_{side}_{i}_TM={primer_result[f'PRIMER_{side}_{i}_TM']}")
        
        if f"PRIMER_{side}_{i}_PENALTY" in primer_result:
            lines.append(f"PRIMER_{side}_{i}_PENALTY={primer_result[f'PRIMER_{side}_{i}_PENALTY']}")
    
    def get_amplicon(self, seq, left_start, left_len, right_start, right_len):
        """
        Extract amplicon sequence with corrected coordinate handling.
        
        Returns the substring from the 5' end of the forward primer
        to the 3' end of the reverse primer, addressing systematic
        coordinate alignment issues with Primer3.
        
        Args:
            seq: Template sequence
            left_start: Start position of forward primer (1-based)
            left_len: Length of forward primer
            right_start: Start position of reverse primer (1-based)
            right_len: Length of reverse primer
            
        Returns:
            Amplicon sequence string, empty if extraction fails
            
        Raises:
            ValueError: If coordinate parameters are invalid
            
        Example:
            >>> amplicon = processor.get_amplicon(template, 10, 20, 200, 20)
            >>> print(f"Amplicon length: {len(amplicon)} bp")
        """
        # Validate input parameters
        if None in (left_start, left_len, right_start, right_len):
            logger.warning("Missing required parameters for amplicon extraction")
            return ""
        
        try:
            left_start = int(left_start)
            left_len = int(left_len)
            right_start = int(right_start)
            right_len = int(right_len)
        except (TypeError, ValueError) as e:
            error_msg = f"Invalid numerical values for amplicon extraction"
            logger.error(error_msg)
            logger.debug(f"Parameter conversion error: {str(e)}", exc_info=True)
            raise ValueError(error_msg) from e
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Amplicon extraction: left={left_start},{left_len}, right={right_start},{right_len}")
            
        # Adjust coordinates to address off-by-one issues
        adj_left_start = max(1, left_start - 1)
        
        # Handle edge cases
        if right_start > len(seq):
            logger.debug(f"Adjusting right_start from {right_start} to {len(seq)} (sequence length)")
            right_start = len(seq)
        
        # Calculate amplicon boundaries
        amp_start = adj_left_start
        amp_end = right_start
        
        # Validate coordinates
        if amp_start < 1 or amp_end > len(seq) or amp_start > amp_end:
            logger.warning(f"Invalid amplicon coordinates: start={amp_start}, end={amp_end}, seq_len={len(seq)}")
            return ""
        
        # Extract amplicon sequence (convert to 0-based for Python)
        amplicon = seq[amp_start - 1 : amp_end]
        
        # Debug validation in debug mode
        if logger.isEnabledFor(logging.DEBUG):
            self._validate_amplicon_extraction(seq, left_start, left_len, right_start, right_len, 
                                             adj_left_start, amplicon)
        
        return amplicon
    
    def _validate_amplicon_extraction(self, seq, left_start, left_len, right_start, right_len, 
                                    adj_left_start, amplicon):
        """
        Validate amplicon extraction by comparing primer sequences.
        
        Args:
            seq: Template sequence
            left_start: Original left start position
            left_len: Left primer length
            right_start: Right start position
            right_len: Right primer length
            adj_left_start: Adjusted left start position
            amplicon: Extracted amplicon sequence
        """
        # Extract primer sequences for validation
        forward_original = seq[left_start - 1 : left_start - 1 + left_len]
        forward_adjusted = seq[adj_left_start - 1 : adj_left_start - 1 + left_len]
        reverse_primer = seq[right_start - right_len : right_start]
        
        logger.debug(f"Primer validation:")
        logger.debug(f"  Forward (original pos): {forward_original}")
        logger.debug(f"  Forward (adjusted pos): {forward_adjusted}")
        logger.debug(f"  Amplicon start: {amplicon[:left_len] if len(amplicon) >= left_len else amplicon}")
        logger.debug(f"  Reverse primer: {reverse_primer}")
        logger.debug(f"  Amplicon end: {amplicon[-right_len:] if len(amplicon) >= right_len else amplicon}")
        
        # Check alignment
        if len(amplicon) >= left_len and amplicon[:left_len] != forward_adjusted:
            logger.debug("WARNING: Amplicon start does not match adjusted forward primer")
        
        if len(amplicon) >= right_len and amplicon[-right_len:] != reverse_primer:
            logger.debug("WARNING: Amplicon end does not match reverse primer")
    
    def parse_primer3_batch(self, stdout_data: str, fragment_info: Optional[Dict] = None) -> List[Dict]:
        """
        Parse primer3 output for a batch with comprehensive error handling.
        
        Processes the text output from primer3 into structured primer records
        with detailed validation and error reporting.
        
        Args:
            stdout_data: Primer3 stdout output string
            fragment_info: Dictionary mapping fragment IDs to coordinate information
            
        Returns:
            List of primer record dictionaries ready for downstream processing
            
        Raises:
            SequenceProcessingError: If parsing fails critically
            
        Example:
            >>> records = processor.parse_primer3_batch(output_string, fragment_mapping)
            >>> print(f"Parsed {len(records)} primer records")
        """
        logger.debug("=== PRIMER3 BATCH PARSING ===")
        
        fragment_info = fragment_info or {}
        debug_mode = logger.isEnabledFor(logging.DEBUG)
        
        # Log parsing overview
        if debug_mode:
            input_blocks = len(fragment_info)
            output_blocks = stdout_data.count("SEQUENCE_ID=")
            logger.debug(f"Parsing: {input_blocks} input fragments → {output_blocks} output blocks")
        
        try:
            # Parse sequence blocks from output
            sequence_blocks = self._parse_sequence_blocks(stdout_data)
            
            if debug_mode:
                logger.debug(f"Successfully parsed {len(sequence_blocks)} sequence blocks")
            
            # Process each block into primer records
            records = []
            total_pairs_found = 0
            blocks_with_pairs = 0
            blocks_without_pairs = 0
            
            for block in sequence_blocks:
                try:
                    block_records = self._process_sequence_block(block, fragment_info, debug_mode)
                    records.extend(block_records)
                    
                    if debug_mode:
                        pairs_in_block = len(block_records)
                        total_pairs_found += pairs_in_block
                        if pairs_in_block > 0:
                            blocks_with_pairs += 1
                        else:
                            blocks_without_pairs += 1
                            
                except Exception as e:
                    seq_id = block.get('sequence_id', 'unknown')
                    logger.error(f"Failed to process sequence block {seq_id}")
                    logger.debug(f"Block processing error: {str(e)}", exc_info=True)
                    continue
            
            # Log comprehensive parsing results
            if debug_mode:
                logger.debug(f"Parsing results:")
                logger.debug(f"  Blocks with primers: {blocks_with_pairs}")
                logger.debug(f"  Blocks without primers: {blocks_without_pairs}")
                logger.debug(f"  Total primer pairs: {total_pairs_found}")
                
                # Log k-mer masking impact if enabled - FIXED: use helper method
                if self._is_kmer_masking_available():
                    logger.debug(f"  K-mer template masking: ACTIVE ({self.kmer_file_prefix})")
            
            logger.debug(f"Successfully parsed {len(records)} primer records from Primer3 output")
            logger.debug("=== END PRIMER3 BATCH PARSING ===")
            
            return records
            
        except Exception as e:
            error_msg = f"Critical failure in Primer3 batch parsing"
            logger.error(error_msg)
            logger.debug(f"Parsing error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e

    def _parse_sequence_blocks(self, stdout_data: str) -> List[Dict]:
        """
        Parse stdout data into individual sequence blocks.
        
        Args:
            stdout_data: Raw Primer3 output string
            
        Returns:
            List of parsed sequence block dictionaries
        """
        blocks = []
        current_block = {
            'sequence_id': None,
            'sequence_template': '',
            'primer_pairs': []
        }
        
        lines = stdout_data.splitlines()
        for line in lines:
            line = line.strip()
            
            if line.startswith("SEQUENCE_ID="):
                # Save previous block if it exists
                if current_block['sequence_id']:
                    blocks.append(current_block)
                    
                # Start new block
                current_block = {
                    'sequence_id': line.split("=", 1)[1],
                    'sequence_template': '',
                    'primer_pairs': []
                }
                
            elif line.startswith("SEQUENCE_TEMPLATE="):
                current_block['sequence_template'] = line.split("=", 1)[1].upper()
                
            elif line == "=":
                # End of current block
                if current_block['sequence_id']:
                    blocks.append(current_block)
                current_block = {
                    'sequence_id': None,
                    'sequence_template': '',
                    'primer_pairs': []
                }
                
            else:
                # Parse primer data lines
                self._parse_primer_data_line(line, current_block)
        
        # Handle final block without closing marker
        if current_block['sequence_id']:
            blocks.append(current_block)
        
        return blocks

    def _parse_primer_data_line(self, line: str, block: Dict) -> None:
        """
        Parse a single line of primer data into the sequence block.
        
        Args:
            line: Line to parse
            block: Block dictionary to update
        """
        # Primer pair-level data
        if match := re.match(r'^PRIMER_PAIR_(\d+)_PENALTY=(.*)', line):
            idx, val = int(match.group(1)), float(match.group(2))
            pair = self._get_or_create_primer_pair(block, idx)
            pair["pair_penalty"] = val
            
        elif match := re.match(r'^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(.*)', line):
            idx, val = int(match.group(1)), int(match.group(2))
            pair = self._get_or_create_primer_pair(block, idx)
            pair["product_size"] = val
            
        # Individual primer data
        elif match := re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)_(SEQUENCE|TM|PENALTY)=(.*)', line):
            side, idx, attr, val = match.groups()
            idx = int(idx)
            pair = self._get_or_create_primer_pair(block, idx)
            
            attr_key = f"{side.lower()}_{attr.lower()}"
            if attr in ["TM", "PENALTY"]:
                pair[attr_key] = float(val)
            else:
                pair[attr_key] = val.upper()
                
        # Primer positions
        elif match := re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)=(\d+),(\d+)', line):
            side, idx, start, length = match.groups()
            idx = int(idx)
            start = int(start) + 1  # Convert to 1-based coordinates
            length = int(length)
            pair = self._get_or_create_primer_pair(block, idx)
            
            pair[f"{side.lower()}_start"] = start
            pair[f"{side.lower()}_len"] = length

    def _get_or_create_primer_pair(self, block: Dict, idx: int) -> Dict:
        """
        Get existing primer pair or create new one by index.
        
        Args:
            block: Sequence block dictionary
            idx: Primer pair index
            
        Returns:
            Primer pair dictionary
        """
        # Look for existing pair
        for pair in block['primer_pairs']:
            if pair.get('idx') == idx:
                return pair
        
        # Create new pair if not found
        new_pair = {"idx": idx}
        block['primer_pairs'].append(new_pair)
        return new_pair

    def _process_sequence_block(self, block: Dict, fragment_info: Dict, debug_mode: bool) -> List[Dict]:
        """
        Process a single sequence block into primer records.
        
        Args:
            block: Parsed sequence block
            fragment_info: Fragment information mapping
            debug_mode: Whether debug logging is enabled
            
        Returns:
            List of primer record dictionaries
        """
        if not block['sequence_id'] or not block['primer_pairs']:
            if debug_mode:
                logger.debug(f"Skipping block {block.get('sequence_id', 'unknown')}: no primers found")
            return []
        
        # Log primer pairs in debug mode
        if debug_mode:
            self._log_primer_pairs_debug(block)
        
        # Filter and limit primer pairs
        acceptable_pairs = self._filter_primer_pairs(block['primer_pairs'])
        
        if debug_mode and len(acceptable_pairs) != len(block['primer_pairs']):
            logger.debug(f"Filtered {len(block['primer_pairs'])} → {len(acceptable_pairs)} primer pairs")
        
        # Generate primer records
        records = []
        for pair in acceptable_pairs:
            try:
                record = self._create_primer_record(block, pair, fragment_info, debug_mode)
                if record:
                    records.append(record)
            except Exception as e:
                logger.error(f"Failed to create primer record for pair {pair.get('idx', 'unknown')}")
                logger.debug(f"Record creation error: {str(e)}", exc_info=True)
                continue
        
        return records

    def _log_primer_pairs_debug(self, block: Dict) -> None:
        """
        Log detailed information about primer pairs for debugging.
        
        Args:
            block: Sequence block containing primer pairs
        """
        sequence_id = block['sequence_id']
        sequence_template = block['sequence_template']
        primer_pairs = block['primer_pairs']
        
        logger.debug(f"=== PRIMER PAIRS FOR {sequence_id} ===")
        logger.debug(f"Template length: {len(sequence_template)} bp")
        logger.debug(f"Primer pairs found: {len(primer_pairs)}")
        
        # Sort by penalty for better readability
        sorted_pairs = sorted(primer_pairs, key=lambda p: p.get('pair_penalty', 999))
        
        for i, pair in enumerate(sorted_pairs):
            left_seq = pair.get('left_sequence', '')
            right_seq = pair.get('right_sequence', '')
            penalty = pair.get('pair_penalty', 'N/A')
            product_size = pair.get('product_size', 'N/A')
            
            logger.debug(f"Pair {i+1}: Penalty={penalty}, Size={product_size}")
            logger.debug(f"  Forward: {left_seq}")
            logger.debug(f"  Reverse: {right_seq}")

    def _filter_primer_pairs(self, primer_pairs: List[Dict]) -> List[Dict]:
        """
        Filter and limit primer pairs based on configuration.
        
        Args:
            primer_pairs: List of primer pairs to filter
            
        Returns:
            Filtered and limited primer pairs
        """
        acceptable = primer_pairs
        
        # Apply maximum pairs per segment limit
        if (hasattr(self.config, 'MAX_PRIMER_PAIRS_PER_SEGMENT') and 
            self.config.MAX_PRIMER_PAIRS_PER_SEGMENT > 0):
            
            original_count = len(acceptable)
            acceptable = acceptable[:self.config.MAX_PRIMER_PAIRS_PER_SEGMENT]
            
            if logger.isEnabledFor(logging.DEBUG) and original_count > self.config.MAX_PRIMER_PAIRS_PER_SEGMENT:
                logger.debug(f"Limited to {self.config.MAX_PRIMER_PAIRS_PER_SEGMENT} primer pairs "
                           f"(had {original_count} total)")
        
        return acceptable

    def _create_primer_record(self, block: Dict, pair: Dict, fragment_info: Dict, debug_mode: bool) -> Optional[Dict]:
        """
        Create a primer record dictionary from sequence block and primer pair.
        
        Args:
            block: Sequence block containing template
            pair: Primer pair data
            fragment_info: Fragment information mapping
            debug_mode: Whether debug logging is enabled
            
        Returns:
            Complete primer record dictionary or None if creation fails
            
        Raises:
            ValueError: If required primer data is missing
        """
        try:
            # Extract basic primer sequences
            left_seq = pair.get("left_sequence", "")
            right_seq = pair.get("right_sequence", "")
            
            if not left_seq or not right_seq:
                logger.warning(f"Missing primer sequences for pair {pair.get('idx', 'unknown')}")
                return None
            
            # Process probe sequence
            probe_seq, probe_reversed = self._process_probe_sequence(pair)
            
            # Get position data
            ls, ll = pair.get("left_start"), pair.get("left_len")
            rs, rl = pair.get("right_start"), pair.get("right_len")
            
            # Create amplicon
            amplicon = self._create_amplicon(block['sequence_template'], pair, debug_mode)
            
            # Get location information
            location_info = self._get_location_info(block['sequence_id'], pair, fragment_info)
            
            # Calculate product size
            product_size = pair.get("product_size")
            if product_size is None and amplicon:
                product_size = len(amplicon)
            
            # Build comprehensive primer record
            record = {
                "Gene": location_info["gene"],
                "Index": pair["idx"],
                "Template": block['sequence_template'],
                "Primer F": left_seq,
                "Tm F": pair.get("left_tm"),
                "Penalty F": pair.get("left_penalty"),
                "Primer F Start": ls,
                "Primer F Len": ll,
                "Primer R": right_seq,
                "Tm R": pair.get("right_tm"),
                "Penalty R": pair.get("right_penalty"),
                "Primer R Start": rs,
                "Primer R Len": rl,
                "Pair Penalty": pair.get("pair_penalty"),
                "Amplicon": amplicon,
                "Length": product_size,
                "Chromosome": location_info["chromosome"],
                "Location": location_info["location"]
            }
            
            # Add probe-related fields if internal oligos are enabled
            if not self.config.DISABLE_INTERNAL_OLIGO:
                internal_start = pair.get("internal_start")
                internal_len = pair.get("internal_len")
                
                record.update({
                    "Probe": probe_seq,
                    "Probe Tm": pair.get("internal_tm"),
                    "Probe Penalty": pair.get("internal_penalty"),
                    "Probe Start": internal_start,
                    "Probe Len": internal_len,
                    "Probe Reversed": probe_reversed
                })
            
            return record
            
        except Exception as e:
            error_msg = f"Failed to create primer record for pair {pair.get('idx', 'unknown')}"
            logger.error(error_msg)
            logger.debug(f"Record creation error: {str(e)}", exc_info=True)
            return None

    def _process_probe_sequence(self, pair: Dict) -> Tuple[str, bool]:
        """
        Process probe sequence with optional reverse complementation.
        
        Args:
            pair: Primer pair containing probe data
            
        Returns:
            Tuple of (probe_sequence, was_reversed)
        """
        if self.config.DISABLE_INTERNAL_OLIGO:
            return "", False
        
        probe_seq = pair.get("internal_sequence", "")
        probe_reversed = False
        
        # Apply reverse complementation if configured and probe exists
        if self.config.PREFER_PROBE_MORE_C_THAN_G and probe_seq:
            try:
                probe_seq, probe_reversed = SequenceUtils.ensure_more_c_than_g(probe_seq)
                
                if probe_reversed and logger.isEnabledFor(logging.DEBUG):
                    logger.debug(f"Reversed probe for optimal C>G ratio")
                    
            except Exception as e:
                logger.warning(f"Failed to process probe sequence: {str(e)}")
                logger.debug(f"Probe processing error: {str(e)}", exc_info=True)
        
        return probe_seq, probe_reversed

    def _create_amplicon(self, template: str, pair: Dict, debug_mode: bool) -> str:
        """
        Create amplicon sequence from template and primer positions.
        
        Args:
            template: Template sequence
            pair: Primer pair with position information
            debug_mode: Whether debug logging is enabled
            
        Returns:
            Amplicon sequence string, empty if creation fails
        """
        ls, ll = pair.get("left_start"), pair.get("left_len")
        rs, rl = pair.get("right_start"), pair.get("right_len")
        
        # Validate position data
        if not all([ls, ll, rs, rl]):
            if debug_mode:
                logger.debug(f"Missing position data for amplicon creation in pair {pair.get('idx', 'unknown')}")
                logger.debug(f"  Left: {ls},{ll}, Right: {rs},{rl}")
            return ""
        
        try:
            # Use existing amplicon extraction method
            amplicon = self.get_amplicon(template, ls, ll, rs, rl)
            
            # Detailed debug logging
            if debug_mode:
                if not amplicon:
                    logger.debug(f"WARNING: Failed to create amplicon for pair {pair.get('idx', 'unknown')}")
                    logger.debug(f"  Coordinates: left={ls},{ll}, right={rs},{rl}")
                    logger.debug(f"  Template length: {len(template)}")
                    
                    # Attempt alternative extraction
                    amplicon = self._try_reconstruct_amplicon(template, ls, rs, debug_mode)
                else:
                    logger.debug(f"Created amplicon for pair {pair.get('idx', 'unknown')}: {len(amplicon)} bp")
                    
                    # Validate primer alignment
                    self._validate_primer_alignment(pair, amplicon, debug_mode)
            
            return amplicon
            
        except Exception as e:
            error_msg = f"Error creating amplicon for pair {pair.get('idx', 'unknown')}"
            logger.error(error_msg)
            logger.debug(f"Amplicon creation error: {str(e)}", exc_info=True)
            return ""

    def _try_reconstruct_amplicon(self, template: str, ls: int, rs: int, debug_mode: bool) -> str:
        """
        Attempt to reconstruct amplicon when normal extraction fails.
        
        Args:
            template: Template sequence
            ls: Left start position
            rs: Right start position
            debug_mode: Whether debug logging is enabled
            
        Returns:
            Reconstructed amplicon or empty string
        """
        if debug_mode:
            logger.debug("Attempting amplicon reconstruction")
        
        # Try direct template extraction
        if (ls is not None and rs is not None and ls <= rs and 
            ls >= 1 and rs <= len(template)):
            
            direct_amplicon = template[ls-1:rs]
            
            if debug_mode:
                logger.debug(f"Reconstructed amplicon: {len(direct_amplicon)} bp")
                
            if len(direct_amplicon) > 0:
                return direct_amplicon
        
        if debug_mode:
            logger.debug(f"Amplicon reconstruction failed: ls={ls}, rs={rs}, template_len={len(template)}")
        
        return ""

    def _validate_primer_alignment(self, pair: Dict, amplicon: str, debug_mode: bool) -> None:
        """
        Validate that amplicon contains expected primer sequences.
        
        Args:
            pair: Primer pair data
            amplicon: Extracted amplicon sequence
            debug_mode: Whether debug logging is enabled
        """
        if not debug_mode or not amplicon:
            return
        
        left_seq = pair.get("left_sequence", "")
        right_seq = pair.get("right_sequence", "")
        
        # Check forward primer alignment
        if left_seq and len(amplicon) >= len(left_seq):
            if not amplicon.startswith(left_seq[:min(len(left_seq), 10)]):
                logger.debug("WARNING: Amplicon does not start with forward primer")
        
        # Check reverse primer alignment (reverse complement)
        if right_seq and len(amplicon) >= len(right_seq):
            try:
                rc_right = SequenceUtils.reverse_complement(right_seq)
                if rc_right[:min(len(right_seq), 10)] not in amplicon[-len(right_seq):]:
                    logger.debug("WARNING: Amplicon does not end with reverse primer complement")
            except Exception as e:
                logger.debug(f"Could not validate reverse primer alignment: {e}")

    def _get_location_info(self, sequence_id: str, pair: Dict, fragment_info: Dict) -> Dict[str, str]:
        """
        Get location information for primer record.
        
        Args:
            sequence_id: Sequence identifier
            pair: Primer pair data
            fragment_info: Fragment information mapping
            
        Returns:
            Dictionary with gene, chromosome, and location information
        """
        # Get fragment information
        frag_info = fragment_info.get(sequence_id, {})
        
        # Determine gene name
        gene = frag_info.get("gene", sequence_id)
        
        # Get chromosome
        chromosome = frag_info.get("chr", "")
        
        # Calculate absolute genomic position if possible
        location = ""
        ls = pair.get("left_start")
        if ls is not None and "start" in frag_info:
            try:
                fragment_start = frag_info.get("start", 1)
                abs_left_start = fragment_start + ls - 1
                location = str(abs_left_start)
            except (TypeError, ValueError):
                logger.debug(f"Could not calculate absolute position for {sequence_id}")
        
        return {
            "gene": gene,
            "chromosome": chromosome, 
            "location": location
        }