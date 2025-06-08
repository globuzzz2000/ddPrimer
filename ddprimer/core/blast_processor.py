import os
import time
import tempfile
import subprocess
import concurrent.futures
import logging
import pandas as pd
from tqdm import tqdm

# Import package modules
from ..config import Config, SequenceProcessingError

class BlastProcessor:
    """Handles BLAST operations for primer specificity checking."""
    
    # Get module logger
    logger = logging.getLogger("ddPrimer.blast_processor")
    
    @staticmethod
    def blast_short_seq(seq, db=None):
        """
        Runs BLASTn for short sequences and returns the two lowest e-values separately.
        
        Args:
            seq (str): Sequence to BLAST
            db (str, optional): BLAST database path. Defaults to Config.DB_PATH.
            
        Returns:
            tuple: (best_evalue, second_best_evalue)
            
        Raises:
            SequenceProcessingError: When BLAST operation fails
        """
        if db is None:
            db = f'"{Config.DB_PATH}"'
            
        if not seq or not isinstance(seq, str) or not seq.strip():
            return None, None  # Ensuring consistency in output

        # Use a temporary file for the query sequence
        tmp_filename = None
        try:
            with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_query:
                tmp_query.write(f">seq\n{seq}\n")
                tmp_query.flush()
                tmp_filename = tmp_query.name

            # Run BLASTn
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
            
            # If BLAST fails, log the error and raise an exception
            if result.returncode != 0:
                BlastProcessor.logger.error(f"BLAST Error: {result.stderr}")
                raise SequenceProcessingError(f"BLAST execution failed: {result.stderr}")

            # Parse BLAST output
            try:
                evalues = sorted([float(line.strip()) for line in result.stdout.strip().split("\n") if line.strip()])
            except ValueError as e:
                BlastProcessor.logger.warning(f"Error parsing BLAST output: {e}")
                evalues = []

            if not evalues:
                return None, None

            # Return the two lowest e-values as separate values
            best = evalues[0] if len(evalues) > 0 else None
            second = evalues[1] if len(evalues) > 1 else None

            return best, second  # Separate columns
            
        except Exception as e:
            BlastProcessor.logger.error(f"Error in BLAST operation: {e}")
            BlastProcessor.logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None, None
            
        finally:
            # Ensure we clean up even if an exception occurs
            if tmp_filename:
                try:
                    os.remove(tmp_filename)  # Clean up temp file
                except OSError as e:
                    BlastProcessor.logger.debug(f"Failed to remove temp file: {e}")

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