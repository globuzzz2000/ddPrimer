import os
import time
import tempfile
import subprocess
import concurrent.futures
import pandas as pd
from tqdm import tqdm
from ..config import Config

class BlastProcessor:
    """Handles BLAST operations for primer specificity checking."""
    
    @staticmethod
    def blast_short_seq(seq, db=None):
        """
        Runs BLASTn for short sequences and returns the two lowest e-values separately.
        
        Args:
            seq (str): Sequence to BLAST
            db (str): BLAST database path (default: from Config)
            
        Returns:
            tuple: (best_evalue, second_best_evalue)
        """
        if db is None:
            db = f'"{Config.DB_PATH}"'
            
        if not seq or not isinstance(seq, str) or not seq.strip():
            return None, None  # Ensuring consistency in output

        # Use a temporary file for the query sequence
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_query:
            tmp_query.write(f">seq\n{seq}\n")
            tmp_query.flush()
            tmp_filename = tmp_query.name

        # Run BLASTn
        try:
            result = subprocess.run(
                [
                    "blastn",
                    "-task", "blastn-short",
                    "-db", db,
                    "-query", tmp_filename,
                    "-word_size", str(Config.BLAST_WORD_SIZE),
                    "-evalue", str(Config.BLAST_EVALUE),
                    "-reward", str(Config.BLAST_REWARD),
                    "-penalty", str(Config.BLAST_PENALTY),
                    "-gapopen", str(Config.BLAST_GAPOPEN),
                    "-gapextend", str(Config.BLAST_GAPEXTEND),
                    "-max_target_seqs", str(Config.BLAST_MAX_TARGET_SEQS),
                    "-outfmt", "6 evalue"
                ],
                text=True,
                capture_output=True
            )
        finally:
            # Ensure we clean up even if an exception occurs
            try:
                os.remove(tmp_filename)  # Clean up temp file
            except OSError:
                pass

        # If BLAST fails, print the error and return None values
        if result.returncode != 0:
            print(f"BLAST Error: {result.stderr}")
            return None, None

        # Parse BLAST output correctly
        try:
            evalues = sorted([float(line.strip()) for line in result.stdout.strip().split("\n") if line.strip()])
        except ValueError:
            evalues = []

        if not evalues:
            return None, None

        # Return the two lowest e-values as separate values
        best = evalues[0] if len(evalues) > 0 else None
        second = evalues[1] if len(evalues) > 1 else None

        return best, second  # Separate columns

    @classmethod
    def process_blast_batch(cls, batch_data):
        """
        Process a batch of sequences for BLAST in parallel.
        
        Args:
            batch_data (tuple): Tuple of (batch, col_name)
            
        Returns:
            list: List of (best_evalue, second_best_evalue) tuples
        """
        batch, col_name = batch_data
        results = []
        
        for seq in batch:
            if pd.notnull(seq):
                blast1, blast2 = cls.blast_short_seq(seq)
                results.append((blast1, blast2))
            else:
                results.append((None, None))
        
        return results
    
    @staticmethod
    def passes_blast_filter(row, col):
        """
        Check if a sequence passes the BLAST specificity filter.
        
        Args:
            row (pandas.Series): DataFrame row
            col (str): Column prefix
            
        Returns:
            bool: True if passes filter, False otherwise
        """
        best = row[f"{col} BLAST1"]
        second = row[f"{col} BLAST2"]
        # If best is None, no hits => discard
        if pd.isna(best):
            return False
        # If no second best, it's effectively unique
        if pd.isna(second):
            return True
        # best must be at least FILTER_FACTOR times smaller => best * FACTOR <= second
        return best * Config.BLAST_FILTER_FACTOR <= second

    @classmethod
    def validate_sequences_batch(cls, batch, ref_genome, qry_genome, extra_genomes=None):
        """
        Process a batch of sequences in parallel with optimized validation.
        
        Args:
            batch (list): List of sequence dictionaries to validate
            ref_genome (dict): Reference genome dictionary
            qry_genome (dict): Query genome dictionary
            extra_genomes (list): List of additional genome dictionaries
            
        Returns:
            list: List of valid sequence dictionaries
        """
        if extra_genomes is None:
            extra_genomes = []
        
        # Build indices if needed (once per genome)
        if Config.VALIDATION_MODE == "TOLERANT":
            cls.build_sequence_index(ref_genome)
            cls.build_sequence_index(qry_genome)
            for genome in extra_genomes:
                cls.build_sequence_index(genome)
            
        start_time = time.time()
        Config.debug(f"Validating batch of {len(batch)} sequences...")
        
        # Use ThreadPoolExecutor for shared memory advantage with the indices
        with concurrent.futures.ThreadPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            # Process each sequence
            future_to_seg = {
                executor.submit(cls.validate_sequence, seg, ref_genome, qry_genome, extra_genomes): seg 
                for seg in batch
            }
            
            # Collect valid sequences
            valid_sequences = []
            
            # Use tqdm for progress if enabled
            if Config.SHOW_PROGRESS:
                futures_completed = tqdm(
                    concurrent.futures.as_completed(future_to_seg),
                    total=len(future_to_seg),
                    desc="Validating sequences"
                )
            else:
                futures_completed = concurrent.futures.as_completed(future_to_seg)
                
            for future in futures_completed:
                segment = future_to_seg[future]
                try:
                    is_valid = future.result()
                    if is_valid:
                        valid_sequences.append(segment)
                        
                    seq_preview = segment['sequence'][:30] + "..." if len(segment['sequence']) > 30 else segment['sequence']
                    Config.debug(f"Sequence validation result: {seq_preview} - {'Valid' if is_valid else 'Invalid'}")
                except Exception as e:
                    Config.debug(f"Error validating sequence: {e}")
        
        elapsed = time.time() - start_time
        Config.debug(f"Batch validation completed: {len(valid_sequences)}/{len(batch)} valid in {elapsed:.2f}s")
        return valid_sequences