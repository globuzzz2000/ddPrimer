#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer3 processing module for ddPrimer pipeline.
"""

import primer3
import re
import os
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from tqdm import tqdm

# Change to use the ddPrimer package config
from ..config import Config
from ..utils.sequence_utils import SequenceUtils  # Use SequenceUtils from your package

class Primer3Processor:
    """Handles Primer3 operations for primer design using the primer3-py package."""
    
    def __init__(self, config=None):
        """
        Initialize Primer3Processor.
        
        Args:
            config (Config, optional): Configuration class with primer3 settings
        """
        # Default to the shared Config if none provided
        self.config = config if config is not None else Config
    
    def run_primer3_batch(self, input_blocks):
        """
        Run primer3 on a batch of input blocks using Python bindings.
        
        Args:
            input_blocks (list): List of dictionaries containing primer3 parameters
            
        Returns:
            str: Combined primer3 output that mimics primer3_core output format
        """
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
            
            return "\n".join(results)
        except Exception as e:
            print(f"Error running Primer3: {e}")
            return ""
    
    def run_primer3_batch_parallel(self, input_blocks, max_workers=None):
        """
        Run primer3 on a batch of input blocks using parallel processing.
        
        Args:
            input_blocks (list): List of dictionaries containing primer3 parameters
            max_workers (int): Maximum number of worker processes
            
        Returns:
            str: Combined primer3 output
        """
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), len(input_blocks))
        
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
        For both forward and reverse primers, adjusts positions based on common 
        Primer3 coordinate interpretation issues.
        """
        import logging
        logger = logging.getLogger("ddPrimer")
        
        # Check if any required parameter is None
        if None in (left_start, left_len, right_start, right_len):
            return ""
        
        # Convert to integers just in case
        try:
            left_start = int(left_start)
            left_len = int(left_len)
            right_start = int(right_start)
            right_len = int(right_len)
        except (TypeError, ValueError):
            return ""
        
        # Debug info to help diagnose issues
        logger.debug(f"Extracting amplicon - Original: left_start={left_start}, left_len={left_len}, right_start={right_start}, right_len={right_len}")
            
        # Adjust left_start (considering possible off-by-one error for the forward primer)
        # Look at the logs to see if we need left_start-1 or left_start
        # Based on our observations, it seems we need to adjust by -1
        adj_left_start = left_start - 1
        if adj_left_start < 1:
            adj_left_start = 1  # Can't go below 1
            
        # Handle edge case: if right_start is beyond the end of the template
        if right_start > len(seq):
            logger.debug(f"Adjusting right_start from {right_start} to {len(seq)} (template length)")
            right_start = len(seq)
        
        # The amplicon starts at the adjusted position and ends at the right position
        amp_start = adj_left_start
        amp_end = right_start
        
        # Validate coordinates
        if amp_start < 1 or amp_end > len(seq) or amp_start > amp_end:
            logger.debug(f"Invalid amplicon coordinates: start={amp_start}, end={amp_end}, seq_len={len(seq)}")
            return ""
        
        # Extract the amplicon sequence (convert to 0-based indexing for Python)
        amplicon = seq[amp_start - 1 : amp_end]
        
        # For diagnostic purposes, extract what we think the primers are:
        forward_primer = seq[left_start - 1 : left_start - 1 + left_len]
        reverse_primer = seq[right_start - right_len : right_start]
        
        # And also extract from our adjusted positions
        adj_forward_primer = seq[adj_left_start - 1 : adj_left_start - 1 + left_len]
        
        logger.debug(f"Forward primer from original position: {forward_primer}")
        logger.debug(f"Forward primer from adjusted position: {adj_forward_primer}")
        logger.debug(f"Start of amplicon: {amplicon[:left_len]}")
        logger.debug(f"Reverse primer: {reverse_primer}")
        logger.debug(f"End of amplicon: {amplicon[-right_len:]}")
        
        return amplicon
    
    def parse_primer3_batch(self, stdout_data, fragment_info=None):
        """
        Parse primer3 output for a batch.
        Returns a list of records.
        
        Args:
            stdout_data (str): Primer3 stdout output
            fragment_info (dict, optional): Dictionary mapping fragment IDs to chromosome and location info
            
        Returns:
            list: List of primer record dictionaries
        """
        import logging
        logger = logging.getLogger("ddPrimer")
        
        records = []
        current_id = None
        sequence_template = ""
        pairs = []
        
        # Add detailed logging for direct mode
        debug_mode = hasattr(self.config, 'DEBUG_MODE') and self.config.DEBUG_MODE
        
        # Create empty fragment_info if none provided
        if fragment_info is None:
            fragment_info = {}
        
        # Helper function: log all primer pairs for a sequence
        def log_all_primer_pairs(record_id, sequence, pairs_data):
            """Helper to log all primer pairs found for a sequence"""
            if not debug_mode:
                return
                
            logger.debug(f"===== ALL PRIMER PAIRS FOR {record_id} =====")
            logger.debug(f"Sequence length: {len(sequence)} bp")
            logger.debug(f"Total pairs found: {len(pairs_data)}")
            
            # Sort pairs by penalty
            sorted_pairs = sorted(pairs_data, key=lambda p: p.get('pair_penalty', 999))
            
            # Log each pair's details
            for i, pair in enumerate(sorted_pairs):
                left_seq = pair.get('left_sequence', '')
                right_seq = pair.get('right_sequence', '')
                penalty = pair.get('pair_penalty', 'N/A')
                product_size = pair.get('product_size', 'N/A')
                
                logger.debug(f"Pair #{i+1}: Penalty={penalty}")
                logger.debug(f"  Forward: {left_seq}")
                logger.debug(f"  Reverse: {right_seq}")
                logger.debug(f"  Product size: {product_size}")

        # Helper function: finalize the current record
        def finalize_record():
            """
            Save all primer pairs regardless of penalty.
            """
            nonlocal current_id, pairs, sequence_template

            if not current_id or not pairs:
                return

            # MODIFIED: Skip penalty filtering completely
            acceptable = pairs
                
            # Take up to MAX_PRIMER_PAIRS_PER_SEGMENT, but only if we have limits
            if hasattr(self.config, 'MAX_PRIMER_PAIRS_PER_SEGMENT') and self.config.MAX_PRIMER_PAIRS_PER_SEGMENT > 0:
                acceptable = acceptable[:self.config.MAX_PRIMER_PAIRS_PER_SEGMENT]
                
            # Debug log amplicon creation
            if debug_mode:
                logger.debug(f"Creating amplicons for {current_id}...")
            
            for p in acceptable:
                left_seq = p.get("left_sequence", "")
                right_seq = p.get("right_sequence", "")
                
                # Only get probe sequence if internal oligos are enabled
                probe_seq = ""
                probe_reversed = False
                if not self.config.DISABLE_INTERNAL_OLIGO:
                    probe_seq = p.get("internal_sequence", "")
                    
                    # Check if we need to reverse complement the probe based on C/G content
                    if self.config.PREFER_PROBE_MORE_C_THAN_G and probe_seq:
                        probe_seq, probe_reversed = SequenceUtils.ensure_more_c_than_g(probe_seq)
                        
                ls, ll = p.get("left_start"), p.get("left_len")
                rs, rl = p.get("right_start"), p.get("right_len")
                
                # Only get internal positions if internal oligos are enabled
                internal_start = None
                internal_len = None
                if not self.config.DISABLE_INTERNAL_OLIGO:
                    internal_start = p.get("internal_start")
                    internal_len = p.get("internal_len")

                # Get amplicon and log details
                ampseq = self.get_amplicon(sequence_template, ls, ll, rs, rl)
                
                if debug_mode:
                    if not ampseq:
                        logger.debug(f"WARNING: Could not create amplicon for pair {p['idx']} in {current_id}")
                        logger.debug(f"  Left start: {ls}, Left len: {ll}")
                        logger.debug(f"  Right start: {rs}, Right len: {rl}")
                        logger.debug(f"  Template length: {len(sequence_template)}")
                    else:
                        logger.debug(f"Created amplicon for pair {p['idx']} in {current_id}: {len(ampseq)}bp")
                        logger.debug(f"  Left start: {ls}, Left len: {ll}")
                        logger.debug(f"  Right start: {rs}, Right len: {rl}")
                        if len(ampseq) > 40:
                            logger.debug(f"  Amplicon (excerpt): {ampseq[:20]}...{ampseq[-20:]}")
                        else:
                            logger.debug(f"  Amplicon: {ampseq}")
                        
                        # Verify that the amplicon actually contains the primers
                        if left_seq and not ampseq.startswith(left_seq[:min(len(left_seq), 10)]):
                            logger.debug(f"  WARNING: Amplicon does not start with forward primer")
                        if right_seq and SequenceUtils.reverse_complement(right_seq)[:min(len(right_seq), 10)] not in ampseq[-len(right_seq):]:
                            logger.debug(f"  WARNING: Amplicon does not end with reverse primer complement")
                
                # Get chromosome and location info from fragment_info
                frag_info = fragment_info.get(current_id, {})
                chromosome = frag_info.get("chr", "")
                
                # Calculate absolute positions based on fragment coordinates
                fragment_start = frag_info.get("start", 1)
                
                # Adjust the primer positions to absolute coordinates
                abs_left_start = None
                if ls is not None:
                    abs_left_start = fragment_start + ls - 1
                    
                # Format the location as a range
                location = ""
                if abs_left_start is not None:
                    location = f"{abs_left_start}"
                
                # For debugging, check the amplicon - if it's empty, create it again
                if debug_mode and not ampseq and sequence_template:
                    logger.debug(f"Attempting to reconstruct missing amplicon for {current_id}")
                    # Try with direct template extraction
                    if ls is not None and rs is not None and ls <= rs and ls >= 1 and rs <= len(sequence_template):
                        direct_amplicon = sequence_template[ls-1:rs]
                        logger.debug(f"Reconstructed amplicon: {len(direct_amplicon)}bp")
                        if len(direct_amplicon) > 0:
                            ampseq = direct_amplicon
                            logger.debug(f"Successfully reconstructed amplicon")
                
                # Check for valid amplicon
                if not ampseq:
                    # Try calculating product size based on primer positions
                    product_size = rs - ls + 1 if ls is not None and rs is not None else None
                    
                    # Try extracting the amplicon one more time
                    ampseq = sequence_template[ls-1:rs] if (ls is not None and rs is not None and 
                                                        1 <= ls <= len(sequence_template) and 
                                                        1 <= rs <= len(sequence_template) and 
                                                        ls <= rs) else ""
                    
                    if debug_mode:
                        if ampseq:
                            logger.debug(f"Recovered amplicon from sequence: {len(ampseq)} bp")
                        else:
                            logger.debug(f"Could not recover amplicon. ls={ls}, rs={rs}, seq_len={len(sequence_template)}")
                else:
                    # Use the reported product size if available
                    product_size = p.get("product_size", None)
                    
                    # If not available, calculate from the amplicon
                    if product_size is None and ampseq:
                        product_size = len(ampseq)
                
                # Create the record
                rec = {
                    "Gene": frag_info.get("gene", current_id),
                    "Index": p["idx"],
                    "Template": sequence_template,
                    "Primer F": left_seq,
                    "Tm F": p.get("left_tm", None),
                    "Penalty F": p.get("left_penalty", None),
                    "Primer F Start": ls,
                    "Primer F Len": ll,
                    "Primer R": right_seq,
                    "Tm R": p.get("right_tm", None),
                    "Penalty R": p.get("right_penalty", None),
                    "Primer R Start": rs,
                    "Primer R Len": rl,
                    "Pair Penalty": p.get("pair_penalty", None),
                    "Amplicon": ampseq,
                    "Length": product_size,
                    "Chromosome": chromosome,
                    "Location": location
                }
                
                # Only add probe-related fields if internal oligos are enabled
                if not self.config.DISABLE_INTERNAL_OLIGO:
                    rec.update({
                        "Probe": probe_seq,
                        "Probe Tm": p.get("internal_tm", None),
                        "Probe Penalty": p.get("internal_penalty", None),
                        "Probe Start": internal_start,
                        "Probe Len": internal_len,
                        "Probe Reversed": probe_reversed
                    })
                    
                records.append(rec)

            current_id = None
            sequence_template = ""
            pairs = []

        lines = stdout_data.splitlines()
        for line in lines:
            line = line.strip()
            if line.startswith("SEQUENCE_ID="):
                # finalize any existing record first
                finalize_record()
                current_id = line.split("=", 1)[1]
                sequence_template = ""
                pairs = []

            elif line.startswith("SEQUENCE_TEMPLATE="):
                sequence_template = line.split("=", 1)[1].upper()

            elif re.match(r'^PRIMER_PAIR_(\d+)_PENALTY=', line):
                match = re.match(r'^PRIMER_PAIR_(\d+)_PENALTY=(.*)', line)
                idx, val = int(match.group(1)), float(match.group(2))
                pair = next((p for p in pairs if p['idx'] == idx), None)
                if not pair:
                    pair = {"idx": idx}
                    pairs.append(pair)
                pair["pair_penalty"] = val

            elif re.match(r'^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=', line):
                match = re.match(r'^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(.*)', line)
                idx, val = int(match.group(1)), int(match.group(2))
                pair = next((p for p in pairs if p['idx'] == idx), None)
                if pair:
                    pair["product_size"] = val

            elif re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)_(SEQUENCE|TM|PENALTY)=', line):
                match = re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)_(SEQUENCE|TM|PENALTY)=(.*)', line)
                side, idx, attr, val = match.groups()
                idx = int(idx)
                pair = next((p for p in pairs if p['idx'] == idx), None)
                if not pair:
                    pair = {"idx": idx}
                    pairs.append(pair)
                attr_key = f"{side.lower()}_{attr.lower()}"
                if attr in ["TM", "PENALTY"]:
                    pair[attr_key] = float(val)
                else:
                    pair[attr_key] = val.upper()

            elif re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)=(\d+),(\d+)', line):
                match = re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)=(\d+),(\d+)', line)
                side, idx, start, length = match.groups()
                idx = int(idx)
                start = int(start) + 1  # Convert to 0-based to 1-based for consistency
                length = int(length)
                pair = next((p for p in pairs if p['idx'] == idx), None)
                if not pair:
                    pair = {"idx": idx}
                    pairs.append(pair)
                pair[f"{side.lower()}_start"] = start
                pair[f"{side.lower()}_len"] = length

            elif line == "=":
                # Before finalizing record, log all pairs found for this sequence
                if debug_mode and current_id and pairs:
                    log_all_primer_pairs(current_id, sequence_template, pairs)
                
                # end of current record
                finalize_record()

        # If the last record never ended with '='
        if current_id and pairs and debug_mode:
            log_all_primer_pairs(current_id, sequence_template, pairs)
        finalize_record()
        
        return records