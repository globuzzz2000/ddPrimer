#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer3 processing module for ddPrimer pipeline.

Handles primer design through Primer3 including:
1. Running primer3 to design primers and probes
2. Parsing primer3 output
3. Handling amplicon extraction and validation
"""

import logging
import re
import primer3
import multiprocessing
from typing import Dict, List, Optional
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Import package modules
from ..config import Config
from ..utils import SequenceUtils
from ..config import SequenceProcessingError

class Primer3Processor:
    """Handles Primer3 operations for primer design using the primer3-py package."""
    
    # Get module logger
    logger = logging.getLogger("ddPrimer.primer3_processor")
    
    def __init__(self, config=None):
        """
        Initialize Primer3Processor.
        
        Args:
            config (Config, optional): Configuration class with primer3 settings. Defaults to Config.
        """
        # Default to the shared Config if none provided
        self.config = config if config is not None else Config
        self.logger.debug("Initialized Primer3Processor")
    
    def run_primer3_batch(self, input_blocks):
        """
        Run primer3 on a batch of input blocks using Python bindings.
        
        Args:
            input_blocks (list): List of dictionaries containing primer3 parameters
            
        Returns:
            str: Combined primer3 output that mimics primer3_core output format
        """
        self.logger.debug(f"Running Primer3 on {len(input_blocks)} input blocks")
        results = []
        
        try:
            # Process each input block using primer3-py
            for block in input_blocks:
                # Extract the sequence ID and template
                sequence_id = block.get("SEQUENCE_ID", "UNKNOWN")
                
                # Run primer3 design
                primer_result = primer3.bindings.design_primers(
                    seq_args=block,
                    global_args=self.config.get_primer3_global_args()
                )
                
                # Format the result to match primer3_core output format
                formatted_result = self._format_primer3_result(sequence_id, block, primer_result)
                results.append(formatted_result)
            
            self.logger.debug(f"Successfully processed {len(results)} primer3 blocks")
            return "\n".join(results)
        except Exception as e:
            self.logger.error(f"Error running Primer3: {e}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            return ""
    
    def run_primer3_batch_parallel(self, input_blocks, max_workers=None):
        """
        Run primer3 on a batch of input blocks using parallel processing.
        
        Args:
            input_blocks (list): List of dictionaries containing primer3 parameters
            max_workers (int, optional): Maximum number of worker processes. Defaults to CPU count.
            
        Returns:
            str: Combined primer3 output
        """
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), len(input_blocks))
        
        self.logger.debug(f"Running Primer3 in parallel with {max_workers} workers for {len(input_blocks)} blocks")
        
        # Split input blocks into chunks for parallel processing
        chunk_size = max(1, len(input_blocks) // max_workers)
        chunks = [input_blocks[i:i + chunk_size] for i in range(0, len(input_blocks), chunk_size)]
        
        # Create a processor instance for each worker
        processor = self
        
        results = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(processor.run_primer3_batch, chunk) for chunk in chunks]
            
            # Use tqdm for progress indication if enabled
            if self.config.SHOW_PROGRESS:
                for future in tqdm(futures, total=len(futures), desc="Running Primer3"):
                    results.append(future.result())
            else:
                for future in futures:
                    results.append(future.result())
        
        self.logger.debug(f"Completed parallel Primer3 processing")
        return "\n".join(results)
    
    def _format_primer3_result(self, sequence_id, input_block, primer_result):
        """
        Format primer3-py result to match primer3_core output format.
        This ensures compatibility with the existing parser.
        
        Args:
            sequence_id (str): ID of the sequence
            input_block (dict): Input parameters
            primer_result (dict): Results from primer3.bindings.design_primers
            
        Returns:
            str: Formatted output string
        """
        lines = [f"SEQUENCE_ID={sequence_id}"]
        
        # Add template sequence
        if "SEQUENCE_TEMPLATE" in input_block:
            lines.append(f"SEQUENCE_TEMPLATE={input_block['SEQUENCE_TEMPLATE']}")
        
        # Get number of primers found
        num_pairs = 0
        while f"PRIMER_LEFT_{num_pairs}_SEQUENCE" in primer_result:
            num_pairs += 1
        
        lines.append(f"PRIMER_PAIR_NUM_RETURNED={num_pairs}")
        
        # Add each primer pair's info
        for i in range(num_pairs):
            # Penalty
            if f"PRIMER_PAIR_{i}_PENALTY" in primer_result:
                lines.append(f"PRIMER_PAIR_{i}_PENALTY={primer_result[f'PRIMER_PAIR_{i}_PENALTY']}")
            
            # Product size
            if f"PRIMER_PAIR_{i}_PRODUCT_SIZE" in primer_result:
                lines.append(f"PRIMER_PAIR_{i}_PRODUCT_SIZE={primer_result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}")
            
            # Left primer
            if f"PRIMER_LEFT_{i}_SEQUENCE" in primer_result:
                lines.append(f"PRIMER_LEFT_{i}_SEQUENCE={primer_result[f'PRIMER_LEFT_{i}_SEQUENCE']}")
            if f"PRIMER_LEFT_{i}" in primer_result:
                start, length = primer_result[f"PRIMER_LEFT_{i}"]
                # Convert to the expected format (0-based to 1-based will be handled in parser)
                lines.append(f"PRIMER_LEFT_{i}={start},{length}")
            if f"PRIMER_LEFT_{i}_TM" in primer_result:
                lines.append(f"PRIMER_LEFT_{i}_TM={primer_result[f'PRIMER_LEFT_{i}_TM']}")
            if f"PRIMER_LEFT_{i}_PENALTY" in primer_result:
                lines.append(f"PRIMER_LEFT_{i}_PENALTY={primer_result[f'PRIMER_LEFT_{i}_PENALTY']}")
            
            # Right primer
            if f"PRIMER_RIGHT_{i}_SEQUENCE" in primer_result:
                lines.append(f"PRIMER_RIGHT_{i}_SEQUENCE={primer_result[f'PRIMER_RIGHT_{i}_SEQUENCE']}")
            if f"PRIMER_RIGHT_{i}" in primer_result:
                start, length = primer_result[f"PRIMER_RIGHT_{i}"]
                lines.append(f"PRIMER_RIGHT_{i}={start},{length}")
            if f"PRIMER_RIGHT_{i}_TM" in primer_result:
                lines.append(f"PRIMER_RIGHT_{i}_TM={primer_result[f'PRIMER_RIGHT_{i}_TM']}")
            if f"PRIMER_RIGHT_{i}_PENALTY" in primer_result:
                lines.append(f"PRIMER_RIGHT_{i}_PENALTY={primer_result[f'PRIMER_RIGHT_{i}_PENALTY']}")
            
            # Internal oligo (probe)
            if f"PRIMER_INTERNAL_{i}_SEQUENCE" in primer_result:
                lines.append(f"PRIMER_INTERNAL_{i}_SEQUENCE={primer_result[f'PRIMER_INTERNAL_{i}_SEQUENCE']}")
            if f"PRIMER_INTERNAL_{i}" in primer_result:
                start, length = primer_result[f"PRIMER_INTERNAL_{i}"]
                lines.append(f"PRIMER_INTERNAL_{i}={start},{length}")
            if f"PRIMER_INTERNAL_{i}_TM" in primer_result:
                lines.append(f"PRIMER_INTERNAL_{i}_TM={primer_result[f'PRIMER_INTERNAL_{i}_TM']}")
            if f"PRIMER_INTERNAL_{i}_PENALTY" in primer_result:
                lines.append(f"PRIMER_INTERNAL_{i}_PENALTY={primer_result[f'PRIMER_INTERNAL_{i}_PENALTY']}")
        
        # Add end marker
        lines.append("=")
        
        return "\n".join(lines)
    
    def get_amplicon(self, seq, left_start, left_len, right_start, right_len):
        """
        Return the substring from the 5' end of the forward primer
        to the 3' end of the reverse primer, with corrected alignment.
        
        Addresses systematic off-by-one issues with Primer3 coordinates.
        
        Args:
            seq (str): Template sequence
            left_start (int): Start position of the forward primer (1-based)
            left_len (int): Length of the forward primer
            right_start (int): Start position of the reverse primer (1-based)
            right_len (int): Length of the reverse primer
            
        Returns:
            str: Amplicon sequence
        """
        # Check if any required parameter is None
        if None in (left_start, left_len, right_start, right_len):
            self.logger.warning("Missing required parameters for amplicon extraction")
            return ""
        
        # Convert to integers just in case
        try:
            left_start = int(left_start)
            left_len = int(left_len)
            right_start = int(right_start)
            right_len = int(right_len)
        except (TypeError, ValueError) as e:
            self.logger.error(f"Invalid numerical values for amplicon extraction: {e}")
            return ""
        
        # Debug info to help diagnose issues
        self.logger.debug(f"Extracting amplicon - Original: left_start={left_start}, left_len={left_len}, right_start={right_start}, right_len={right_len}")
            
        # Adjust left_start (considering possible off-by-one error for the forward primer)
        # Based on our observations, it seems we need to adjust by -1
        adj_left_start = left_start - 1
        if adj_left_start < 1:
            adj_left_start = 1  # Can't go below 1
            self.logger.debug(f"Adjusted left_start to minimum value (1)")
            
        # Handle edge case: if right_start is beyond the end of the template
        if right_start > len(seq):
            self.logger.debug(f"Adjusting right_start from {right_start} to {len(seq)} (template length)")
            right_start = len(seq)
        
        # The amplicon starts at the adjusted position and ends at the right position
        amp_start = adj_left_start
        amp_end = right_start
        
        # Validate coordinates
        if amp_start < 1 or amp_end > len(seq) or amp_start > amp_end:
            self.logger.debug(f"Invalid amplicon coordinates: start={amp_start}, end={amp_end}, seq_len={len(seq)}")
            return ""
        
        # Extract the amplicon sequence (convert to 0-based indexing for Python)
        amplicon = seq[amp_start - 1 : amp_end]
        
        # For diagnostic purposes, extract what we think the primers are:
        forward_primer = seq[left_start - 1 : left_start - 1 + left_len]
        reverse_primer = seq[right_start - right_len : right_start]
        
        # And also extract from our adjusted positions
        adj_forward_primer = seq[adj_left_start - 1 : adj_left_start - 1 + left_len]
        
        self.logger.debug(f"Forward primer from original position: {forward_primer}")
        self.logger.debug(f"Forward primer from adjusted position: {adj_forward_primer}")
        self.logger.debug(f"Start of amplicon: {amplicon[:left_len]}")
        self.logger.debug(f"Reverse primer: {reverse_primer}")
        self.logger.debug(f"End of amplicon: {amplicon[-right_len:]}")
        
        return amplicon
    
    def parse_primer3_batch(self, stdout_data: str, fragment_info: Optional[Dict] = None) -> List[Dict]:
        """
        Parse primer3 output for a batch.
        
        DROP-IN REPLACEMENT for the existing method in Primer3Processor.
        
        Args:
            stdout_data (str): Primer3 stdout output
            fragment_info (dict, optional): Dictionary mapping fragment IDs to chromosome and location info
            
        Returns:
            list: List of primer record dictionaries
        """
        self.logger.debug("Parsing Primer3 batch output")
        
        # Initialize parsing state
        records = []
        fragment_info = fragment_info or {}
        debug_mode = hasattr(self.config, 'DEBUG_MODE') and self.config.DEBUG_MODE
        
        # Parse all sequence blocks from the output
        sequence_blocks = self._parse_sequence_blocks(stdout_data)
        
        # Process each sequence block
        for block in sequence_blocks:
            block_records = self._process_sequence_block(block, fragment_info, debug_mode)
            records.extend(block_records)
        
        self.logger.debug(f"Parsed {len(records)} primer records from Primer3 output")
        return records

    def _parse_sequence_blocks(self, stdout_data: str) -> List[Dict]:
        """
        Parse the stdout data into individual sequence blocks.
        
        Args:
            stdout_data (str): Raw Primer3 output
            
        Returns:
            List[Dict]: List of parsed sequence blocks
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
                # Start of new block - save previous if it exists
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
        
        # Handle case where last block doesn't end with '='
        if current_block['sequence_id']:
            blocks.append(current_block)
        
        return blocks

    def _parse_primer_data_line(self, line: str, block: Dict) -> None:
        """
        Parse a single line of primer data and add it to the sequence block.
        
        Args:
            line (str): Line to parse
            block (Dict): Block to add data to
        """
        # Primer pair penalty
        if match := re.match(r'^PRIMER_PAIR_(\d+)_PENALTY=(.*)', line):
            idx, val = int(match.group(1)), float(match.group(2))
            pair = self._get_or_create_primer_pair(block, idx)
            pair["pair_penalty"] = val
            
        # Primer pair product size
        elif match := re.match(r'^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(.*)', line):
            idx, val = int(match.group(1)), int(match.group(2))
            pair = self._get_or_create_primer_pair(block, idx)
            pair["product_size"] = val
            
        # Primer sequences, TM, and penalties
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
        Get or create a primer pair by index within a sequence block.
        
        Args:
            block (Dict): Sequence block
            idx (int): Primer pair index
            
        Returns:
            Dict: Existing or new primer pair dictionary
        """
        # Look for existing pair with this index
        for pair in block['primer_pairs']:
            if pair.get('idx') == idx:
                return pair
        
        # Create new pair if not found
        new_pair = {"idx": idx}
        block['primer_pairs'].append(new_pair)
        return new_pair

    def _process_sequence_block(self, block: Dict, fragment_info: Dict, debug_mode: bool) -> List[Dict]:
        """
        Process a single sequence block and generate primer records.
        
        Args:
            block (Dict): Parsed sequence block
            fragment_info (Dict): Fragment information mapping
            debug_mode (bool): Whether debug logging is enabled
            
        Returns:
            List[Dict]: List of primer record dictionaries
        """
        if not block['sequence_id'] or not block['primer_pairs']:
            return []
        
        # Log all primer pairs in debug mode
        if debug_mode:
            self._log_all_primer_pairs(block)
        
        # Filter and limit primer pairs
        acceptable_pairs = self._filter_primer_pairs(block['primer_pairs'])
        
        # Generate records for each acceptable pair
        records = []
        for pair in acceptable_pairs:
            record = self._create_primer_record(block, pair, fragment_info, debug_mode)
            if record:
                records.append(record)
        
        return records

    def _filter_primer_pairs(self, primer_pairs: List[Dict]) -> List[Dict]:
        """
        Filter and limit primer pairs based on configuration.
        
        Args:
            primer_pairs (List[Dict]): List of primer pairs to filter
            
        Returns:
            List[Dict]: Filtered and limited primer pairs
        """
        # MODIFIED: Skip penalty filtering completely (keep all pairs)
        acceptable = primer_pairs
        
        # Take up to MAX_PRIMER_PAIRS_PER_SEGMENT if configured
        if (hasattr(self.config, 'MAX_PRIMER_PAIRS_PER_SEGMENT') and 
            self.config.MAX_PRIMER_PAIRS_PER_SEGMENT > 0):
            acceptable = acceptable[:self.config.MAX_PRIMER_PAIRS_PER_SEGMENT]
        
        return acceptable

    def _create_primer_record(self, block: Dict, pair: Dict, fragment_info: Dict, debug_mode: bool) -> Optional[Dict]:
        """
        Create a primer record dictionary from a sequence block and primer pair.
        
        Args:
            block (Dict): Sequence block containing template
            pair (Dict): Primer pair data
            fragment_info (Dict): Fragment information
            debug_mode (bool): Whether debug logging is enabled
            
        Returns:
            Optional[Dict]: Primer record or None if creation fails
        """
        # Get basic sequences
        left_seq = pair.get("left_sequence", "")
        right_seq = pair.get("right_sequence", "")
        
        # Handle probe sequence based on configuration
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
        
        # Build the record
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
            internal_start = pair.get("internal_start") if not self.config.DISABLE_INTERNAL_OLIGO else None
            internal_len = pair.get("internal_len") if not self.config.DISABLE_INTERNAL_OLIGO else None
            
            record.update({
                "Probe": probe_seq,
                "Probe Tm": pair.get("internal_tm"),
                "Probe Penalty": pair.get("internal_penalty"),
                "Probe Start": internal_start,
                "Probe Len": internal_len,
                "Probe Reversed": probe_reversed
            })
        
        return record

    def _process_probe_sequence(self, pair: Dict) -> tuple:
        """
        Process probe sequence with optional reverse complementation.
        
        Args:
            pair (Dict): Primer pair containing probe data
            
        Returns:
            tuple: (probe_sequence, was_reversed)
        """
        if self.config.DISABLE_INTERNAL_OLIGO:
            return "", False
        
        probe_seq = pair.get("internal_sequence", "")
        probe_reversed = False
        
        # Check if we need to reverse complement the probe
        if self.config.PREFER_PROBE_MORE_C_THAN_G and probe_seq:
            probe_seq, probe_reversed = SequenceUtils.ensure_more_c_than_g(probe_seq)
        
        return probe_seq, probe_reversed

    def _create_amplicon(self, template: str, pair: Dict, debug_mode: bool) -> str:
        """
        Create amplicon sequence from template and primer positions.
        
        Args:
            template (str): Template sequence
            pair (Dict): Primer pair with position information
            debug_mode (bool): Whether debug logging is enabled
            
        Returns:
            str: Amplicon sequence
        """
        ls, ll = pair.get("left_start"), pair.get("left_len")
        rs, rl = pair.get("right_start"), pair.get("right_len")
        
        # Check if we have the required position data
        if not all([ls, ll, rs, rl]):
            if debug_mode:
                self.logger.debug(f"Missing position data for amplicon creation in pair {pair['idx']}")
            return ""
        
        # Use the existing get_amplicon method for consistency
        amplicon = self.get_amplicon(template, ls, ll, rs, rl)
        
        # Debug logging for amplicon creation
        if debug_mode:
            current_id = pair.get('sequence_id', 'unknown')
            if not amplicon:
                self.logger.debug(f"WARNING: Could not create amplicon for pair {pair['idx']} in {current_id}")
                self.logger.debug(f"  Left start: {ls}, Left len: {ll}")
                self.logger.debug(f"  Right start: {rs}, Right len: {rl}")
                self.logger.debug(f"  Template length: {len(template)}")
                
                # Try to reconstruct amplicon
                amplicon = self._try_reconstruct_amplicon(template, ls, rs, debug_mode)
            else:
                self.logger.debug(f"Created amplicon for pair {pair['idx']}: {len(amplicon)}bp")
                self.logger.debug(f"  Left start: {ls}, Left len: {ll}")
                self.logger.debug(f"  Right start: {rs}, Right len: {rl}")
                if len(amplicon) > 40:
                    self.logger.debug(f"  Amplicon (excerpt): {amplicon[:20]}...{amplicon[-20:]}")
                else:
                    self.logger.debug(f"  Amplicon: {amplicon}")
                
                # Verify amplicon contains primers
                left_seq = pair.get("left_sequence", "")
                right_seq = pair.get("right_sequence", "")
                if left_seq and not amplicon.startswith(left_seq[:min(len(left_seq), 10)]):
                    self.logger.debug(f"  WARNING: Amplicon does not start with forward primer")
                if right_seq and SequenceUtils.reverse_complement(right_seq)[:min(len(right_seq), 10)] not in amplicon[-len(right_seq):]:
                    self.logger.debug(f"  WARNING: Amplicon does not end with reverse primer complement")
        
        return amplicon

    def _try_reconstruct_amplicon(self, template: str, ls: int, rs: int, debug_mode: bool) -> str:
        """
        Try to reconstruct amplicon when normal extraction fails.
        
        Args:
            template (str): Template sequence
            ls (int): Left start position
            rs (int): Right start position
            debug_mode (bool): Whether debug logging is enabled
            
        Returns:
            str: Reconstructed amplicon or empty string
        """
        if debug_mode:
            self.logger.debug(f"Attempting to reconstruct missing amplicon")
        
        # Try with direct template extraction
        if (ls is not None and rs is not None and ls <= rs and 
            ls >= 1 and rs <= len(template)):
            direct_amplicon = template[ls-1:rs]
            if debug_mode:
                self.logger.debug(f"Reconstructed amplicon: {len(direct_amplicon)}bp")
            if len(direct_amplicon) > 0:
                if debug_mode:
                    self.logger.debug(f"Successfully reconstructed amplicon")
                return direct_amplicon
        
        if debug_mode:
            self.logger.debug(f"Could not recover amplicon. ls={ls}, rs={rs}, seq_len={len(template)}")
        
        return ""

    def _get_location_info(self, sequence_id: str, pair: Dict, fragment_info: Dict) -> Dict[str, str]:
        """
        Get location information for the primer record.
        
        Args:
            sequence_id (str): Sequence identifier
            pair (Dict): Primer pair data
            fragment_info (Dict): Fragment information mapping
            
        Returns:
            Dict[str, str]: Location information with keys: gene, chromosome, location
        """
        # Get fragment info
        frag_info = fragment_info.get(sequence_id, {})
        
        # Determine gene name
        gene = frag_info.get("gene", sequence_id)
        
        # Get chromosome
        chromosome = frag_info.get("chr", "")
        
        # Calculate absolute position if we have fragment coordinates
        location = ""
        ls = pair.get("left_start")
        if ls is not None and "start" in frag_info:
            fragment_start = frag_info.get("start", 1)
            abs_left_start = fragment_start + ls - 1
            location = str(abs_left_start)
        
        return {
            "gene": gene,
            "chromosome": chromosome, 
            "location": location
        }

    def _log_all_primer_pairs(self, block: Dict) -> None:
        """
        Log all primer pairs for debugging purposes.
        
        Args:
            block (Dict): Sequence block to log
        """
        sequence_id = block['sequence_id']
        sequence_template = block['sequence_template']
        primer_pairs = block['primer_pairs']
        
        self.logger.debug(f"===== ALL PRIMER PAIRS FOR {sequence_id} =====")
        self.logger.debug(f"Sequence length: {len(sequence_template)} bp")
        self.logger.debug(f"Total pairs found: {len(primer_pairs)}")
        
        # Sort pairs by penalty for logging
        sorted_pairs = sorted(primer_pairs, key=lambda p: p.get('pair_penalty', 999))
        
        for i, pair in enumerate(sorted_pairs):
            left_seq = pair.get('left_sequence', '')
            right_seq = pair.get('right_sequence', '')
            penalty = pair.get('pair_penalty', 'N/A')
            product_size = pair.get('product_size', 'N/A')
            
            self.logger.debug(f"Pair #{i+1}: Penalty={penalty}")
            self.logger.debug(f"  Forward: {left_seq}")
            self.logger.debug(f"  Reverse: {right_seq}")
            self.logger.debug(f"  Product size: {product_size}")