#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST processing module for ddPrimer pipeline.

Handles BLAST operations for primer specificity checking including:
1. BLASTn execution for short sequences
2. E-value parsing and analysis
3. Identity filtering (100% identity requirement)
4. Batch processing capabilities
5. Specificity filtering
"""

import os
import tempfile
import subprocess
import logging
import pandas as pd
from tqdm import tqdm
from typing import Optional

# Import package modules
from ..config import Config, SequenceProcessingError, PrimerDesignError, DebugLogLimiter
from ..core import FilterProcessor

# Set up module logger
logger = logging.getLogger(__name__)


class BlastProcessor:
    """
    Handles BLAST operations for primer specificity checking.
    
    This class provides methods for running BLASTn searches on short DNA sequences
    and evaluating their specificity based on e-value distributions and identity filtering.
    
    Example:
        >>> blast1, blast2 = BlastProcessor.blast_short_seq("ATCGATCGATCG")
        >>> if blast1 and blast2:
        ...     specificity_ratio = blast1 / blast2
    """
    
    #############################################################################
    #                           Workflow Wrappers
    #############################################################################
    
    @classmethod
    def run_blast_specificity_workflow(cls, df: pd.DataFrame) -> Optional[pd.DataFrame]:
        """
        Run BLAST for primer specificity checking for workflow integration.
        
        Executes BLAST analysis for all primers and probes to assess specificity,
        then filters results based on BLAST e-value thresholds and identity requirements.
        
        Args:
            df: DataFrame with primer information
            
        Returns:
            DataFrame with added BLAST results and filtered for specificity,
            or None if no primers pass BLAST filtering
            
        Raises:
            PrimerDesignError: If there's an error in BLAST execution or filtering
        """
        logger.info("\nRunning BLAST for specificity checking...")
        logger.debug("=== WORKFLOW: BLAST SPECIFICITY CHECKING ===")
        logger.debug(f"DataFrame columns: {df.columns.tolist()}")
        
        try:
            # Use new column naming convention
            forward_col = "Sequence (F)"
            reverse_col = "Sequence (R)"
            probe_col = "Sequence (P)"
            
            if forward_col not in df.columns or reverse_col not in df.columns:
                error_msg = f"Required columns not found: {forward_col}, {reverse_col}"
                logger.error(error_msg)
                logger.debug(f"Available columns: {df.columns.tolist()}")
                raise PrimerDesignError(error_msg)
            
            logger.debug(f"Using new column naming convention")
            
            # Run BLAST for forward primers
            blast_results_f = []
            primers_f = df[forward_col].tolist()
            if Config.SHOW_PROGRESS:
                primers_f_iter = tqdm(primers_f, total=len(primers_f), desc="BLASTing forward primers")
            else:
                primers_f_iter = primers_f
                
            for primer_f in primers_f_iter:
                blast1, blast2 = cls.blast_short_seq(primer_f)
                blast_results_f.append((blast1, blast2))
            
            # Add both BLAST values to DataFrame for filtering
            df["BLAST (F)"] = [result[0] for result in blast_results_f]   # Best e-value
            df["BLAST (F2)"] = [result[1] for result in blast_results_f]  # Second best e-value
            
            # Run BLAST for reverse primers
            blast_results_r = []
            primers_r = df[reverse_col].tolist()
            if Config.SHOW_PROGRESS:
                primers_r_iter = tqdm(primers_r, total=len(primers_r), desc="BLASTing reverse primers")
            else:
                primers_r_iter = primers_r
                
            for primer_r in primers_r_iter:
                blast1, blast2 = cls.blast_short_seq(primer_r)
                blast_results_r.append((blast1, blast2))
            
            # Add both reverse BLAST values for filtering
            df["BLAST (R)"] = [result[0] for result in blast_results_r]   # Best e-value
            df["BLAST (R2)"] = [result[1] for result in blast_results_r]  # Second best e-value
            
            # Run BLAST for probes if present
            if probe_col in df.columns:
                blast_results_p = []
                probes = df[probe_col].tolist()
                if Config.SHOW_PROGRESS:
                    probes_iter = tqdm(probes, total=len(probes), desc="BLASTing probes")
                else:
                    probes_iter = probes
                    
                for probe in probes_iter:
                    if pd.notnull(probe) and probe:
                        blast1, blast2 = cls.blast_short_seq(probe)
                        blast_results_p.append((blast1, blast2))
                    else:
                        blast_results_p.append((None, None))
                
                # Add both probe BLAST values for filtering
                df["BLAST (P)"] = [result[0] for result in blast_results_p]   # Best e-value
                df["BLAST (P2)"] = [result[1] for result in blast_results_p]  # Second best e-value
            
            # Filter by BLAST specificity using FilterProcessor (uses both values)
            initial_count = len(df)
            df = FilterProcessor.filter_by_blast(df)
            
            if len(df) == 0:
                logger.warning("No primers passed BLAST specificity filtering.")
                logger.debug("=== END WORKFLOW: BLAST SPECIFICITY CHECKING ===")
                return None
            
            # REMOVE the second BLAST columns before returning (keep only best e-values in output)
            columns_to_remove = []
            if "BLAST (F2)" in df.columns:
                columns_to_remove.append("BLAST (F2)")
            if "BLAST (R2)" in df.columns:
                columns_to_remove.append("BLAST (R2)")
            if "BLAST (P2)" in df.columns:
                columns_to_remove.append("BLAST (P2)")
            
            if columns_to_remove:
                df = df.drop(columns=columns_to_remove)
                logger.debug(f"Removed temporary BLAST columns: {columns_to_remove}")
            
            if len(df) == 0:
                logger.warning("No primers passed BLAST specificity filtering.")
                logger.debug("=== END WORKFLOW: BLAST SPECIFICITY CHECKING ===")
                return None
            
            logger.info(f"Retained {len(df)} primer pairs after BLAST specificity filtering.")
            logger.debug("=== END WORKFLOW: BLAST SPECIFICITY CHECKING ===")
            
            return df
            
        except Exception as e:
            error_msg = f"Error in BLAST specificity workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: BLAST SPECIFICITY CHECKING ===")
            raise PrimerDesignError(error_msg) from e
    
    #############################################################################
    #                           BLAST Execution
    #############################################################################
   
    @staticmethod
    def blast_short_seq(seq, db=None):
        """
        Run BLASTn for short sequences and return the two best e-values separately.
        
        Executes a BLASTn search optimized for short sequences and extracts
        the best and second-best e-values for specificity assessment. Only considers
        hits with 100% identity to ensure specificity.
        
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
            if DebugLogLimiter.should_log('blast_invalid_seq', interval=500, max_initial=2):
                logger.debug("Empty or invalid sequence provided to BLAST")
            return None, None

        tmp_filename = None
        try:
            # Create temporary file for query sequence
            if not hasattr(BlastProcessor, "_centralized_temp_dir"):
                BlastProcessor._centralized_temp_dir = os.path.join(Config.get_user_config_dir(), "temp")
            os.makedirs(BlastProcessor._centralized_temp_dir, exist_ok=True)

            with tempfile.NamedTemporaryFile(mode="w+", delete=False, dir=BlastProcessor._centralized_temp_dir) as tmp_query:
                tmp_query.write(f">seq\n{seq}\n")
                tmp_query.flush()
                tmp_filename = tmp_query.name

            # Execute BLASTn command with extended output format including identity
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
                    "-outfmt", "6 evalue pident"  # Include percent identity
                ],
                text=True,
                capture_output=True
            )
            
            if result.returncode != 0:
                error_msg = f"BLAST execution failed for sequence {seq[:20]}... (length: {len(seq)}, db: {db})"
                logger.error(error_msg)
                logger.debug(f"BLAST stderr: {result.stderr}", exc_info=True)
                raise SequenceProcessingError(error_msg)

            # Parse BLAST output for e-values with 100% identity requirement
            try:
                valid_evalues = []
                for line in result.stdout.strip().split("\n"):
                    if line.strip():
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            evalue = float(parts[0])
                            pident = float(parts[1])
                            
                            # Only include hits with 100% identity
                            if pident == 100.0:
                                valid_evalues.append(evalue)
                            elif DebugLogLimiter.should_log('blast_identity_filter', interval=1000, max_initial=2):
                                logger.debug(f"Filtered BLAST hit with {pident}% identity (sequence {seq[:15]}...)")
                
                # Sort by e-value (best first)
                evalues = sorted(valid_evalues)
                
            except ValueError as e:
                if DebugLogLimiter.should_log('blast_parsing_errors', interval=200, max_initial=3):
                    logger.warning(f"Error parsing BLAST output for sequence {seq[:20]}... (length: {len(seq)})")
                    logger.debug(f"BLAST parsing error: {str(e)}", exc_info=True)
                evalues = []

            if not evalues:
                if DebugLogLimiter.should_log('blast_no_hits', interval=1000, max_initial=2):
                    logger.debug(f"No BLAST hits with 100% identity found for sequence {seq[:20]}... (length: {len(seq)})")
                return None, None

            # Return best and second-best e-values
            best = evalues[0] if len(evalues) > 0 else None
            second = evalues[1] if len(evalues) > 1 else None

            # Limited logging for BLAST results
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('blast_results_details', interval=1000, max_initial=2)):
                logger.debug(f"BLAST results for {seq[:20]}... -> Best: {best}, Second: {second} (100% identity only)")
            
            return best, second
            
        except SequenceProcessingError:
            # Re-raise without wrapping
            raise
        except Exception as e:
            error_msg = f"Unexpected BLAST error for sequence {seq[:20]}... (length: {len(seq)})"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None, None
            
        finally:
            # Clean up temporary file
            if tmp_filename:
                try:
                    os.remove(tmp_filename)
                except OSError as e:
                    if DebugLogLimiter.should_log('blast_cleanup_errors', interval=500, max_initial=2):
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
        
        logger.debug(f"=== BLAST BATCH PROCESSING DEBUG ===")
        logger.debug(f"Processing BLAST batch of {len(batch)} sequences for {col_name}")
        
        failed_count = 0
        processed_count = 0
        
        for seq in batch:
            if pd.notnull(seq):
                try:
                    blast1, blast2 = cls.blast_short_seq(seq)
                    results.append((blast1, blast2))
                    
                    if blast1 is None:
                        failed_count += 1
                        
                    # Limited logging using DebugLogLimiter
                    if (logger.isEnabledFor(logging.DEBUG) and 
                        DebugLogLimiter.should_log('blast_batch_processing', interval=500, max_initial=3)):
                        logger.debug(f"BLAST {processed_count}: {seq[:15]}... -> {blast1}, {blast2}")
                        
                except Exception as e:
                    if DebugLogLimiter.should_log('blast_batch_errors', interval=200, max_initial=3):
                        logger.debug(f"BLAST failed for sequence {processed_count} in batch: {str(e)}")
                    results.append((None, None))
                    failed_count += 1
            else:
                results.append((None, None))
            
            processed_count += 1
        
        logger.debug(f"Completed BLAST batch processing: {len(results)} results, {failed_count} failed")
        logger.debug(f"=== END BLAST BATCH PROCESSING DEBUG ===")
        return results