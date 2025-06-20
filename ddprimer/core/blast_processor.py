#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST processing module for ddPrimer pipeline.

Handles BLAST operations for primer specificity checking including:
1. BLASTn execution for short sequences
2. E-value parsing and analysis
3. Batch processing capabilities
4. Specificity filtering
"""

import os
import tempfile
import subprocess
import logging
import pandas as pd

from ..config import Config, SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)


class BlastProcessor:
    """
    Handles BLAST operations for primer specificity checking.
    
    This class provides methods for running BLASTn searches on short DNA sequences
    and evaluating their specificity based on e-value distributions.
    
    Example:
        >>> blast1, blast2 = BlastProcessor.blast_short_seq("ATCGATCGATCG")
        >>> if blast1 and blast2:
        ...     specificity_ratio = blast1 / blast2
    """
    
    @staticmethod
    def blast_short_seq(seq, db=None):
        """
        Run BLASTn for short sequences and return the two best e-values separately.
        
        Executes a BLASTn search optimized for short sequences and extracts
        the best and second-best e-values for specificity assessment.
        
        Args:
            seq: DNA sequence to BLAST against database
            db: BLAST database path, defaults to Config.DB_PATH
            
        Returns:
            Tuple of (best_evalue, second_best_evalue) or (None, None) if failed
            
        Raises:
            SequenceProcessingError: If BLAST execution fails
            
        Example:
            >>> best, second = BlastProcessor.blast_short_seq("ATCGATCG")
            >>> if best and second:
            ...     print(f"Best e-value: {best}, Second: {second}")
        """
        if db is None:
            db = f'"{Config.DB_PATH}"'
            
        if not seq or not isinstance(seq, str) or not seq.strip():
            logger.debug("Empty or invalid sequence provided to BLAST")
            return None, None

        tmp_filename = None
        try:
            # Create temporary file for query sequence
            centralized_temp = os.path.join(Config.get_user_config_dir(), "temp")
            os.makedirs(centralized_temp, exist_ok=True)

            with tempfile.NamedTemporaryFile(mode="w+", delete=False, dir=centralized_temp) as tmp_query:
                tmp_query.write(f">seq\n{seq}\n")
                tmp_query.flush()
                tmp_filename = tmp_query.name

            logger.debug(f"Running BLAST for sequence: {seq[:20]}{'...' if len(seq) > 20 else ''}")

            # Execute BLASTn command
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
            
            if result.returncode != 0:
                error_msg = f"BLAST execution failed for sequence {seq[:20]}..."
                logger.error(error_msg)
                logger.debug(f"BLAST stderr: {result.stderr}", exc_info=True)
                raise SequenceProcessingError(error_msg)

            # Parse BLAST output for e-values
            try:
                evalues = sorted([
                    float(line.strip()) 
                    for line in result.stdout.strip().split("\n") 
                    if line.strip()
                ])
            except ValueError as e:
                logger.warning(f"Error parsing BLAST output for sequence {seq[:20]}...")
                logger.debug(f"BLAST parsing error: {str(e)}", exc_info=True)
                evalues = []

            if not evalues:
                logger.debug(f"No BLAST hits found for sequence {seq[:20]}...")
                return None, None

            # Return best and second-best e-values
            best = evalues[0] if len(evalues) > 0 else None
            second = evalues[1] if len(evalues) > 1 else None

            logger.debug(f"BLAST results - Best: {best}, Second: {second}")
            return best, second
            
        except SequenceProcessingError:
            # Re-raise without wrapping
            raise
        except Exception as e:
            error_msg = f"Unexpected BLAST error for sequence {seq[:20]}..."
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None, None
            
        finally:
            # Clean up temporary file
            if tmp_filename:
                try:
                    os.remove(tmp_filename)
                except OSError as e:
                    logger.debug(f"Failed to remove temp file {tmp_filename}: {e}")

    @classmethod
    def process_blast_batch(cls, batch_data):
        """
        Process a batch of sequences for BLAST in parallel.
        
        Processes multiple sequences concurrently for improved performance
        in batch BLAST operations.
        
        Args:
            batch_data: Tuple of (batch_sequences, column_name)
            
        Returns:
            List of (best_evalue, second_best_evalue) tuples
            
        Example:
            >>> sequences = ["ATCG", "GCTA", "TTTT"]
            >>> results = BlastProcessor.process_blast_batch((sequences, "primers"))
        """
        batch, col_name = batch_data
        results = []
        
        logger.debug(f"Processing BLAST batch of {len(batch)} sequences for {col_name}")
        
        for seq in batch:
            if pd.notnull(seq):
                try:
                    blast1, blast2 = cls.blast_short_seq(seq)
                    results.append((blast1, blast2))
                except Exception as e:
                    logger.debug(f"BLAST failed for sequence in batch: {str(e)}")
                    results.append((None, None))
            else:
                results.append((None, None))
        
        logger.debug(f"Completed BLAST batch processing: {len(results)} results")
        return results
    
    @staticmethod
    def passes_blast_filter(row, col):
        """
        Check if a sequence passes the BLAST specificity filter.
        
        Evaluates whether the ratio between best and second-best e-values
        meets the specificity threshold defined by the filter factor.
        
        Args:
            row: DataFrame row containing BLAST results
            col: Column prefix for BLAST results (e.g., "Primer F")
            
        Returns:
            True if sequence passes specificity filter, False otherwise
            
        Example:
            >>> passes = BlastProcessor.passes_blast_filter(row, "Primer F")
            >>> if passes:
            ...     print("Sequence has sufficient specificity")
        """
        best = row[f"{col} BLAST1"]
        second = row[f"{col} BLAST2"]
        
        # No hits found - discard for lack of specificity data
        if pd.isna(best):
            return False
            
        # Only one hit found - effectively unique
        if pd.isna(second):
            return True
            
        # Check specificity ratio: best * factor <= second
        passes = best * Config.BLAST_FILTER_FACTOR <= second
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"BLAST filter check: {best} * {Config.BLAST_FILTER_FACTOR} <= {second} = {passes}")
        
        return passes