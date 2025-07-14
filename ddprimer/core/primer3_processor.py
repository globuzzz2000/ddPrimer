#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer3 processing module for ddPrimer pipeline.

Handles primer design through Primer3 including primer design using Python bindings,
parsing output with comprehensive error handling, amplicon extraction and validation,
parallel processing, and result formatting.

Contains functionality for:
1. Running primer3 to design primers and probes using Python bindings
2. Parsing primer3 output with comprehensive error handling
3. Parallel processing for improved performance
4. Result formatting and validation

COORDINATE SYSTEM:
- INPUT: Fragment sequences with 0-based coordinates from AnnotationProcessor
- PRIMER3: Uses 0-based coordinates internally (primer3-py bindings)
- OUTPUT: Primer records with 0-based Primer3 coordinates
- All sequence template indexing is 0-based
- No coordinate conversion needed (primer3-py handles this correctly)
"""

import logging
import re
import primer3
import multiprocessing
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Import package modules
from ..config import Config, SequenceProcessingError, PrimerDesignError, DebugLogLimiter
from .primer_processor import PrimerProcessor

# Set up module logger
logger = logging.getLogger(__name__)


class Primer3Processor:
    """
    Handles Primer3 operations for primer design using the primer3-py package.
    
    This class provides comprehensive primer design capabilities including
    parallel processing, result parsing, and amplicon extraction for various pipeline modes.
    
    Attributes:
        config: Configuration object with primer3 settings
        enable_internal_oligo: Whether to design internal oligos (probes)
    """
    
    def __init__(self, config=None, enable_internal_oligo=True):
        """
        Initialize Primer3Processor with configuration.
        
        Args:
            config: Configuration object, defaults to Config
            enable_internal_oligo: Whether to design internal oligos (probes)
        """
        self.config = config if config is not None else Config
        self.enable_internal_oligo = enable_internal_oligo
        logger.debug(f"Initialized Primer3Processor with internal oligo design: {enable_internal_oligo}")

    #############################################################################
    #                           Workflow Wrappers
    #############################################################################
    
    def design_primers_workflow(self, fragments: List[Dict]) -> List[Dict]:
        """
        Run Primer3 design on the provided fragments for workflow integration.
        
        This workflow wrapper coordinates the complete Primer3 design process including
        input preparation, parallel execution, and result parsing. It prepares fragments
        for Primer3 processing, executes primer design using optimal parallelization,
        and transforms the output into standardized primer records ready for downstream
        filtering and analysis operations.
        
        Args:
            fragments: List of sequence fragments containing template sequences,
                      genomic coordinates, and metadata from restriction site cutting
                      and gene overlap filtering
            
        Returns:
            List of primer design results containing fragment information and
            primer pair data ready for record creation and filtering
            
        Raises:
            PrimerDesignError: If there's an error in Primer3 design workflow coordination
        """
        logger.info("\nDesigning primers with Primer3...")
        logger.debug("=== WORKFLOW: PRIMER3 DESIGN ===")
        
        try:
            primer3_inputs, fragment_info = self._prepare_primer3_inputs(fragments)
            
            if not primer3_inputs:
                logger.warning("No valid fragments for primer design.")
                return []
                
            primer3_output = self.run_primer3_batch_parallel(primer3_inputs)
            primer_results = self.parse_primer3_batch(primer3_output, fragment_info)
            
            logger.debug(f"Primer design complete: {len(primer_results)} primer designs generated")
            logger.debug("=== END WORKFLOW: PRIMER3 DESIGN ===")
            
            return primer_results
            
        except Exception as e:
            error_msg = f"Error in Primer3 design workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: PRIMER3 DESIGN ===")
            raise PrimerDesignError(error_msg) from e
    
    def _prepare_primer3_inputs(self, fragments: List[Dict]) -> Tuple[List[Dict], Dict]:
        """
        Prepare Primer3-compatible input dictionaries and store fragment metadata.
        """
        primer3_inputs = []
        fragment_info = {}
        
        for fragment in fragments:
            # The fragment's 'start' field now holds the true 0-based genomic start
            fragment_info[fragment["id"]] = {
                "chr": fragment.get("chromosome", ""),
                "start": fragment.get("start", 0),  # Use the fragment's genomic start
                "end": fragment.get("end", len(fragment["sequence"])),
                "gene": fragment.get("Gene", fragment.get("gene", fragment["id"].split("_")[-1]))
            }
            
            primer3_input = {
                "SEQUENCE_ID": fragment["id"],
                "SEQUENCE_TEMPLATE": fragment["sequence"],
            }

            primer3_inputs.append(primer3_input)
        
        return primer3_inputs, fragment_info

    #############################################################################
    #                           Primer3 Execution
    #############################################################################

    def get_primer3_global_args(self) -> Dict:
        """
        Get Primer3 global arguments adjusted for internal oligo settings.
        
        Returns standard Primer3 arguments from configuration with adjustments
        for internal oligo (probe) design based on the processor configuration
        settings determined during initialization.
        
        Returns:
            Dictionary of Primer3 arguments ready for primer3.bindings.design_primers,
            with internal oligo setting configured according to processor state
        """
        args = self.config.get_primer3_global_args()
        args["PRIMER_PICK_INTERNAL_OLIGO"] = 1 if self.enable_internal_oligo else 0
        return args
    
    def run_primer3_batch(self, input_blocks):
        """
        Run primer3 on a batch of input blocks using Python bindings.
        
        Processes multiple primer design requests sequentially using the primer3-py
        package, formatting results to maintain compatibility with existing parsing
        infrastructure. Handles individual block failures gracefully while maintaining
        overall batch processing integrity.
        
        Args:
            input_blocks: List of dictionaries containing primer3 parameters including
                         sequence templates and design constraints
            
        Returns:
            Combined primer3 output string in primer3_core format ready for parsing
            by the batch parsing infrastructure
            
        Raises:
            PrimerDesignError: If primer3 execution fails at the batch level
        """
        logger.debug(f"Running Primer3 on {len(input_blocks)} input blocks")
        
        if not input_blocks:
            logger.warning("No input blocks provided for Primer3 processing")
            return ""
        
        results = []
        processing_stats = {'successful': 0, 'failed': 0}
        
        try:
            global_args = self.get_primer3_global_args()
            
            for block_num, block in enumerate(input_blocks):
                sequence_id = block.get("SEQUENCE_ID", "UNKNOWN")
                
                if DebugLogLimiter.should_log('primer3_block_processing', interval=500, max_initial=2):
                    template_length = len(block.get("SEQUENCE_TEMPLATE", ""))
                    logger.debug(f"Processing block {block_num+1}/{len(input_blocks)}: "
                               f"sequence {sequence_id} with {template_length} bp template")
                
                try:
                    primer_result = primer3.bindings.design_primers(
                        seq_args=block,
                        global_args=global_args
                    )
                    
                    formatted_result = self._format_primer3_result(sequence_id, block, primer_result)
                    results.append(formatted_result)
                    processing_stats['successful'] += 1
                    
                except Exception as e:
                    processing_stats['failed'] += 1
                    if DebugLogLimiter.should_log('primer3_processing_errors', interval=20, max_initial=3):
                        logger.debug(f"Primer3 design failed for sequence {sequence_id}")
                        logger.debug(f"Primer3 error details: {str(e)}", exc_info=True)
                    
                    results.append(f"SEQUENCE_ID={sequence_id}\nPRIMER_PAIR_NUM_RETURNED=0\n=")
            
            logger.debug(f"Primer3 batch processing: {processing_stats['successful']} successful, "
                       f"{processing_stats['failed']} failed from {len(input_blocks)} total")
            
            return "\n".join(results)
            
        except Exception as e:
            error_msg = f"Batch primer3 processing failed"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise PrimerDesignError(error_msg) from e
    
    def run_primer3_batch_parallel(self, input_blocks, max_workers=None):
        """
        Run primer3 on input blocks using parallel processing for improved performance.
        
        Distributes primer design tasks across multiple processes to optimize throughput
        while maintaining result integrity. Includes comprehensive progress tracking,
        error handling for individual chunks, and automatic worker scaling based on
        system capabilities and workload size.
        
        Args:
            input_blocks: List of dictionaries containing primer3 parameters for
                         parallel processing across multiple worker processes
            max_workers: Maximum number of worker processes to spawn, defaults to
                        CPU count with automatic scaling based on workload
            
        Returns:
            Combined primer3 output string from all workers, formatted for downstream
            parsing with proper error handling for failed chunks
            
        Raises:
            PrimerDesignError: If parallel processing setup or execution fails at
                              the coordination level
        """
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), len(input_blocks))
        
        logger.debug(f"Running Primer3 in parallel: {max_workers} workers for {len(input_blocks)} blocks")
        
        if logger.isEnabledFor(logging.DEBUG):
            self._log_input_statistics(input_blocks)
        
        try:
            chunk_size = max(1, len(input_blocks) // max_workers)
            chunks = [input_blocks[i:i + chunk_size] for i in range(0, len(input_blocks), chunk_size)]
            
            logger.debug(f"Split input into {len(chunks)} chunks of approximately {chunk_size} blocks each")
            
            results = []
            chunk_stats = {'successful': 0, 'failed': 0}
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(self.run_primer3_batch, chunk) for chunk in chunks]
                
                if self.config.SHOW_PROGRESS:
                    for future in tqdm(futures, total=len(futures), desc="Running Primer3"):
                        try:
                            result = future.result()
                            results.append(result)
                            chunk_stats['successful'] += 1
                        except Exception as e:
                            logger.error(f"Primer3 chunk processing failed: {str(e)}")
                            results.append("")
                            chunk_stats['failed'] += 1
                else:
                    for future in futures:
                        try:
                            result = future.result()
                            results.append(result)
                            chunk_stats['successful'] += 1
                        except Exception as e:
                            logger.error(f"Primer3 chunk processing failed: {str(e)}")
                            results.append("")
                            chunk_stats['failed'] += 1
            
            combined_output = "\n".join(filter(None, results))
            
            if logger.isEnabledFor(logging.DEBUG):
                self._log_output_statistics(combined_output, len(input_blocks))
            
            logger.debug(f"Parallel processing complete: {chunk_stats['successful']}/{len(chunks)} chunks successful")
            return combined_output
            
        except Exception as e:
            error_msg = f"Parallel Primer3 processing failed"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise PrimerDesignError(error_msg) from e
    
    def _log_input_statistics(self, input_blocks):
        """
        Log summary statistics about input blocks for debugging.
        """
        seq_lengths = [len(block.get('SEQUENCE_TEMPLATE', '')) for block in input_blocks]
        if seq_lengths:
            avg_len = sum(seq_lengths) / len(seq_lengths)
            min_len = min(seq_lengths)
            max_len = max(seq_lengths)
            logger.debug(f"Input statistics: {len(input_blocks)} blocks, "
                       f"lengths avg={avg_len:.0f}, min={min_len}, max={max_len}")

    def _log_output_statistics(self, output, input_count):
        """
        Analyze and log Primer3 output statistics for debugging.
        """
        sequence_blocks = output.split("SEQUENCE_ID=")
        num_blocks = len(sequence_blocks) - 1
        logger.debug(f"Primer3 output: {num_blocks} blocks from {input_count} input blocks")
        
        if num_blocks < input_count:
            logger.warning(f"Lost {input_count - num_blocks} blocks during Primer3 processing")
        
        successful_designs = 0
        failed_designs = 0
        total_pairs = 0
        
        for block in sequence_blocks[1:]:
            if "PRIMER_PAIR_NUM_RETURNED=" in block:
                match = re.search(r'PRIMER_PAIR_NUM_RETURNED=(\d+)', block)
                if match:
                    pairs_returned = int(match.group(1))
                    total_pairs += pairs_returned
                    
                    if pairs_returned > 0:
                        successful_designs += 1
                    else:
                        failed_designs += 1
        
        logger.debug(f"Design results: {successful_designs} successful, {failed_designs} failed, "
                   f"{total_pairs} total pairs")
    
    #############################################################################
    #                           Result Formatting
    #############################################################################
    
    def _format_primer3_result(self, sequence_id, input_block, primer_result):
        """
        Format primer3-py result to match primer3_core output format.
        """
        lines = [f"SEQUENCE_ID={sequence_id}"]
        
        if "SEQUENCE_TEMPLATE" in input_block:
            lines.append(f"SEQUENCE_TEMPLATE={input_block['SEQUENCE_TEMPLATE']}")
        
        num_pairs = 0
        while f"PRIMER_LEFT_{num_pairs}_SEQUENCE" in primer_result:
            num_pairs += 1
        
        lines.append(f"PRIMER_PAIR_NUM_RETURNED={num_pairs}")
        
        if "PRIMER_ERROR" in primer_result:
            lines.append(f"PRIMER_ERROR={primer_result['PRIMER_ERROR']}")
        
        for i in range(num_pairs):
            self._format_primer_pair(lines, primer_result, i)
        
        lines.append("=")
        
        return "\n".join(lines)
    
    def _format_primer_pair(self, lines, primer_result, pair_index):
        """
        Format a single primer pair for output.
        """
        i = pair_index
        
        if f"PRIMER_PAIR_{i}_PENALTY" in primer_result:
            lines.append(f"PRIMER_PAIR_{i}_PENALTY={primer_result[f'PRIMER_PAIR_{i}_PENALTY']}")
        
        if f"PRIMER_PAIR_{i}_PRODUCT_SIZE" in primer_result:
            lines.append(f"PRIMER_PAIR_{i}_PRODUCT_SIZE={primer_result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}")
        
        for side in ["LEFT", "RIGHT"]:
            self._format_primer_component(lines, primer_result, i, side)
        
        if self.enable_internal_oligo:
            self._format_primer_component(lines, primer_result, i, "INTERNAL")
    
    def _format_primer_component(self, lines, primer_result, pair_index, side):
        """
        Format a single primer component (left, right, or internal).
        """
        i = pair_index
        
        if f"PRIMER_{side}_{i}_SEQUENCE" in primer_result:
            lines.append(f"PRIMER_{side}_{i}_SEQUENCE={primer_result[f'PRIMER_{side}_{i}_SEQUENCE']}")
        
        if f"PRIMER_{side}_{i}" in primer_result:
            start, length = primer_result[f"PRIMER_{side}_{i}"]
            lines.append(f"PRIMER_{side}_{i}={start},{length}")
        
        if f"PRIMER_{side}_{i}_TM" in primer_result:
            lines.append(f"PRIMER_{side}_{i}_TM={primer_result[f'PRIMER_{side}_{i}_TM']}")
        
        if f"PRIMER_{side}_{i}_PENALTY" in primer_result:
            lines.append(f"PRIMER_{side}_{i}_PENALTY={primer_result[f'PRIMER_{side}_{i}_PENALTY']}")
    
    #############################################################################
    #                           Result Parsing
    #############################################################################

    def parse_primer3_batch(self, stdout_data: str, fragment_info: Optional[Dict] = None) -> List[Dict]:
        """
        Parse primer3 output for a batch with comprehensive error handling.
        
        Processes the text output from primer3 into structured primer records with
        detailed validation and error reporting. Handles parsing of complex primer3
        output format, coordinates fragment metadata reconstruction, and provides
        robust error handling for malformed or incomplete results while maintaining
        processing continuity for valid data.
        
        Args:
            stdout_data: Primer3 stdout output string in primer3_core format containing
                        sequence blocks with primer pair data and metadata
            fragment_info: Optional dictionary mapping fragment IDs to coordinate
                          information for genomic position reconstruction
            
        Returns:
            List of primer record dictionaries ready for downstream processing,
            with complete primer information, coordinates, and validation status
            
        Raises:
            SequenceProcessingError: If parsing fails critically at the batch level
                                   affecting overall workflow continuity
        """
        logger.debug("=== PRIMER3 BATCH PARSING ===")
        
        fragment_info = fragment_info or {}
        debug_mode = logger.isEnabledFor(logging.DEBUG)
        
        if debug_mode:
            input_blocks = len(fragment_info)
            output_blocks = stdout_data.count("SEQUENCE_ID=")
            logger.debug(f"Parsing: {input_blocks} input fragments → {output_blocks} output blocks")
        
        try:
            sequence_blocks = self._parse_sequence_blocks(stdout_data)
            
            if debug_mode:
                logger.debug(f"Successfully parsed {len(sequence_blocks)} sequence blocks")
            
            records = []
            stats = {'blocks_with_pairs': 0, 'blocks_without_pairs': 0, 'total_pairs': 0, 'failed_blocks': 0}
            
            for block_num, block in enumerate(sequence_blocks):
                try:
                    block_records = self._process_sequence_block(block, fragment_info, debug_mode)
                    records.extend(block_records)
                    
                    pairs_in_block = len(block_records)
                    stats['total_pairs'] += pairs_in_block
                    
                    if pairs_in_block > 0:
                        stats['blocks_with_pairs'] += 1
                        
                        if DebugLogLimiter.should_log('primer3_block_success', interval=200, max_initial=2):
                            seq_id = block.get('sequence_id', 'unknown')
                            logger.debug(f"Block {seq_id}: generated {pairs_in_block} primer records")
                    else:
                        stats['blocks_without_pairs'] += 1
                        
                        if DebugLogLimiter.should_log('primer3_block_no_primers', interval=500, max_initial=2):
                            seq_id = block.get('sequence_id', 'unknown')
                            logger.debug(f"Block {seq_id}: no primers generated")
                            
                except Exception as e:
                    seq_id = block.get('sequence_id', 'unknown')
                    stats['failed_blocks'] += 1
                    
                    if DebugLogLimiter.should_log('primer3_block_errors', interval=100, max_initial=2):
                        logger.error(f"Failed to process sequence block {seq_id}")
                        logger.debug(f"Block processing error: {str(e)}", exc_info=True)
                    continue
            
            if debug_mode:
                logger.debug(f"Parsing results: {stats['blocks_with_pairs']} with primers, "
                           f"{stats['blocks_without_pairs']} without, {stats['failed_blocks']} failed, "
                           f"{stats['total_pairs']} total pairs")
            
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
                if current_block['sequence_id']:
                    blocks.append(current_block)
                    
                current_block = {
                    'sequence_id': line.split("=", 1)[1],
                    'sequence_template': '',
                    'primer_pairs': []
                }
                
            elif line.startswith("SEQUENCE_TEMPLATE="):
                current_block['sequence_template'] = line.split("=", 1)[1].upper()
                
            elif line == "=":
                if current_block['sequence_id']:
                    blocks.append(current_block)
                current_block = {
                    'sequence_id': None,
                    'sequence_template': '',
                    'primer_pairs': []
                }
                
            else:
                self._parse_primer_data_line(line, current_block)
        
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
        if match := re.match(r'^PRIMER_PAIR_(\d+)_PENALTY=(.*)', line):
            idx, val = int(match.group(1)), float(match.group(2))
            pair = self._get_or_create_primer_pair(block, idx)
            pair["penalty"] = val
            
        elif match := re.match(r'^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(.*)', line):
            idx, val = int(match.group(1)), int(match.group(2))
            pair = self._get_or_create_primer_pair(block, idx)
            pair["product_size"] = val
            
        elif match := re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)_(SEQUENCE|TM|PENALTY)=(.*)', line):
            side, idx, attr, val = match.groups()
            idx = int(idx)
            pair = self._get_or_create_primer_pair(block, idx)
            
            if side == "LEFT":
                if attr == "SEQUENCE":
                    pair["left_sequence"] = val.upper()
                elif attr == "TM":
                    pair["left_tm"] = float(val)
                elif attr == "PENALTY":
                    pair["left_penalty"] = float(val)
            elif side == "RIGHT":
                if attr == "SEQUENCE":
                    pair["right_sequence"] = val.upper()
                elif attr == "TM":
                    pair["right_tm"] = float(val)
                elif attr == "PENALTY":
                    pair["right_penalty"] = float(val)
            elif side == "INTERNAL":
                if attr == "SEQUENCE":
                    pair["internal_sequence"] = val.upper()
                elif attr == "TM":
                    pair["internal_tm"] = float(val)
                elif attr == "PENALTY":
                    pair["internal_penalty"] = float(val)
                
        elif match := re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)=(\d+),(\d+)', line):
            side, idx, start, length = match.groups()
            idx = int(idx)
            start = int(start)
            length = int(length)
            pair = self._get_or_create_primer_pair(block, idx)
            
            if side == "LEFT":
                pair["left_start"] = start
                pair["left_length"] = length
            elif side == "RIGHT":
                pair["right_start"] = start
                pair["right_length"] = length
            elif side == "INTERNAL":
                pair["internal_start"] = start
                pair["internal_length"] = length

    def _get_or_create_primer_pair(self, block: Dict, idx: int) -> Dict:
        """
        Get existing primer pair or create new one by index.
        
        Args:
            block: Sequence block dictionary
            idx: Primer pair index
            
        Returns:
            Primer pair dictionary
        """
        for pair in block['primer_pairs']:
            if pair.get('idx') == idx:
                return pair
        
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
            return []
        
        records = []
        processor = PrimerProcessor()
        
        for pair_num, pair in enumerate(block['primer_pairs']):
            try:
                # Create a fragment from the block with genomic coordinate information
                fragment = {
                    "id": block['sequence_id'],
                    "sequence": block['sequence_template']
                }
                
                # Add fragment info (now includes proper genomic coordinates)
                if block['sequence_id'] in fragment_info:
                    info = fragment_info[block['sequence_id']]
                    fragment.update({
                        "chromosome": info.get("chr", ""),
                        "start": info.get("start", 0), # This is the fragment's genomic start
                        "end": info.get("end", len(block['sequence_template'])),
                        "Gene": info.get("gene", block['sequence_id'].split("_")[-1])
                    })
                else:
                    # This shouldn't happen with proper workflow, but provide fallback
                    logger.warning(f"Missing fragment info for {block['sequence_id']}")
                    fragment.update({
                        "chromosome": block['sequence_id'],
                        "start": 0, # Fallback to 0
                        "end": len(block['sequence_template']),
                        "Gene": block['sequence_id'].split("_")[-1]
                    })
                
                # Create primer record
                record = processor.create_primer_record(fragment, pair, pair_num)
                if record:
                    records.append(record)
                        
            except Exception as e:
                logger.debug(f"Failed to create primer record for pair {pair.get('idx', 'unknown')}: {e}")
                continue
        
        return records

    def _log_primer_pairs_debug(self, block: Dict) -> None:
        """
        Log primer pair details for debugging.
        """
        sequence_id = block['sequence_id']
        sequence_template = block['sequence_template']
        primer_pairs = block['primer_pairs']
        
        logger.debug(f"=== PRIMER PAIRS FOR {sequence_id} ===")
        logger.debug(f"Template: {len(sequence_template)} bp, Pairs: {len(primer_pairs)}")
        
        sorted_pairs = sorted(primer_pairs, key=lambda p: p.get('penalty', 999))
        pairs_to_log = min(3, len(sorted_pairs))
        
        for i in range(pairs_to_log):
            pair = sorted_pairs[i]
            left_seq = pair.get('left_sequence', '')
            right_seq = pair.get('right_sequence', '')
            penalty = pair.get('penalty', 'N/A')
            product_size = pair.get('product_size', 'N/A')
            
            logger.debug(f"Pair {i+1}: Penalty={penalty}, Size={product_size}")
            logger.debug(f"  Forward: {left_seq}")
            logger.debug(f"  Reverse: {right_seq}")
        
        if len(sorted_pairs) > pairs_to_log:
            logger.debug(f"... and {len(sorted_pairs) - pairs_to_log} more pairs")