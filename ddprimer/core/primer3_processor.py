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
        records = []
        current_id = None
        sequence_template = ""
        pairs = []
        
        # Create empty fragment_info if none provided
        if fragment_info is None:
            fragment_info = {}
        
        def get_amplicon(seq, left_start, left_len, right_start, right_len):
            """
            Return the substring from the 5' end of the forward primer
            to the 3' end of the reverse primer, 1-based coords in the template.
            """
            if None in (left_start, left_len, right_start, right_len):
                return ""
            amp_start = left_start
            amp_end = right_start
            if amp_start < 1 or amp_end > len(seq) or amp_start > amp_end:
                return ""
            return seq[amp_start - 1 : amp_end]

        # Helper function: finalize the current record
        def finalize_record():
                """
                Save up to 3 pairs that meet penalty requirements.
                """
                nonlocal current_id, pairs, sequence_template

                if not current_id or not pairs:
                    return

                acceptable = []
                for p in pairs:
                    pair_penalty = p.get("pair_penalty", 999)
                    if not self.config.DISABLE_INTERNAL_OLIGO:
                        # Only check probe penalty if internal oligos are enabled
                        probe_penalty = p.get("internal_penalty", 999)
                        if pair_penalty <= self.config.PENALTY_MAX and probe_penalty <= self.config.PENALTY_MAX:
                            acceptable.append(p)
                    else:
                        # Ignore probe penalty if internal oligos are disabled
                        if pair_penalty <= self.config.PENALTY_MAX:
                            acceptable.append(p)

                # take up to MAX_PRIMER_PAIRS_PER_SEGMENT
                acceptable = acceptable[:self.config.MAX_PRIMER_PAIRS_PER_SEGMENT]
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

                    ampseq = get_amplicon(sequence_template, ls, ll, rs, rl)
                    
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
                    
                    rec = {
                        "Sequence": current_id, 
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
                        "Length": p.get("product_size", None),
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
                start = int(start) + 1  # Convert to 1-based for consistency
                length = int(length)
                pair = next((p for p in pairs if p['idx'] == idx), None)
                if not pair:
                    pair = {"idx": idx}
                    pairs.append(pair)
                pair[f"{side.lower()}_start"] = start
                pair[f"{side.lower()}_len"] = length

            elif line == "=":
                # end of current record
                finalize_record()

        # If the last record never ended with '='
        finalize_record()
        
        return records