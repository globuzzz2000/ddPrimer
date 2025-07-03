#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter processing module for ddPrimer pipeline.

Contains all primer filtering functionality including:
1. Penalty thresholding and validation
2. GC content analysis and filtering
3. Repeat sequence detection and filtering
4. BLAST specificity evaluation
5. Sequence utility functions (GC content, reverse complement, etc.)

COORDINATE SYSTEM:
- This module works with primer sequences and metadata only
- Input: Primer sequences as strings + metadata DataFrames
- Output: Filtered DataFrames (no coordinate processing)
- Sequence operations are coordinate-independent
- No genomic coordinates involved in filtering logic
"""

import pandas as pd
import logging
from typing import List, Dict, Optional

# Import package modules
from ..config import Config, DebugLogLimiter, PrimerDesignError

# Set up module logger
logger = logging.getLogger(__name__)


class FilterProcessor:
    """
    Handles all primer filtering operations and sequence utilities.
    
    This class provides comprehensive filtering capabilities for primer pairs
    including penalty scoring, sequence composition analysis, specificity
    validation, and sequence utility functions.
    """
    
    #############################################################################
    #                           Workflow Wrappers
    #############################################################################
    
    @classmethod
    def filter_primers_workflow(cls, primer_results: List[Dict]) -> Optional[pd.DataFrame]:
        """
        Filter primer records using comprehensive filtering criteria for workflow integration.
        
        Applies multiple filtering steps including penalty thresholds, repeat sequences,
        and GC content validation to ensure high-quality primer selection. This workflow
        wrapper coordinates all filtering operations and provides standardized error
        handling for the pipeline orchestration layer.
        
        Args:
            primer_results: List of primer record dictionaries containing sequences,
                          penalties, and amplicon information
            
        Returns:
            Filtered DataFrame containing only primers meeting all criteria, or None
            if no primers pass filtering
            
        Raises:
            PrimerDesignError: If there's an error in primer filtering workflow coordination
        """
        logger.debug("=== WORKFLOW: PRIMER FILTERING ===")
        logger.debug(f"Filtering {len(primer_results)} primer records")
        
        try:
            if not primer_results:
                logger.warning("No primer results provided for filtering")
                logger.debug("=== END WORKFLOW: PRIMER FILTERING ===")
                return None
            
            # Convert to DataFrame for filtering
            df = pd.DataFrame(primer_results)
            logger.debug(f"Created DataFrame with {len(df)} rows and {len(df.columns)} columns")
            logger.debug(f"DataFrame columns: {df.columns.tolist()}")
            
            # Step 1: Filter by penalty scores
            logger.debug("Applying penalty filter...")
            df = cls.filter_by_penalty(df)
            if df.empty:
                logger.warning("No primers passed penalty filtering")
                logger.debug("=== END WORKFLOW: PRIMER FILTERING ===")
                return None
            
            # Step 2: Filter by repeat sequences
            logger.debug("Applying repeat sequence filter...")
            df = cls.filter_by_repeats(df)
            if df.empty:
                logger.warning("No primers passed repeat filtering")
                logger.debug("=== END WORKFLOW: PRIMER FILTERING ===")
                return None
            
            # Step 3: Filter by GC content
            logger.debug("Applying GC content filter...")
            df = cls.filter_by_gc_content(df)
            if df.empty:
                logger.warning("No primers passed GC content filtering")
                logger.debug("=== END WORKFLOW: PRIMER FILTERING ===")
                return None
            
            logger.info(f"Retained {len(df)} primer pairs after filtering")
            logger.debug("=== END WORKFLOW: PRIMER FILTERING ===")
            
            return df
            
        except Exception as e:
            error_msg = f"Error in primer filtering workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: PRIMER FILTERING ===")
            raise PrimerDesignError(error_msg) from e
    
    #############################################################################
    
    @staticmethod
    def filter_by_penalty(df, max_penalty=None):
        """
        Filter primers by penalty scores.
        
        CLEANED PENALTY HANDLING:
        - Uses only the "Penalty" column (overall primer pair penalty)
        - No more complex logic with "Pair Penalty" vs "Penalty"
        - Single source of truth for penalty filtering
        
        Args:
            df: DataFrame containing primer data with penalty information
            max_penalty: Maximum allowed penalty score, defaults to Config.PENALTY_MAX
            
        Returns:
            Filtered DataFrame containing only primers within penalty threshold
        """
        if max_penalty is None:
            max_penalty = Config.PENALTY_MAX
        
        logger.debug(f"Penalty threshold: {max_penalty}")
        
        if df.empty:
            logger.warning("Empty DataFrame provided for penalty filtering")
            return df
        
        # Check for the penalty column
        if "Penalty" not in df.columns:
            error_msg = "Required 'Penalty' column not found in DataFrame"
            logger.error(error_msg)
            logger.debug(f"Available columns: {df.columns.tolist()}")
            raise ValueError(error_msg)
        
        initial_count = len(df)
        
        # Collect statistics for logging
        if logger.isEnabledFor(logging.DEBUG):
            failed_primers = []
            invalid_penalties = 0
            
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                penalty = row.get("Penalty", "N/A")
                
                try:
                    will_pass = float(penalty) <= max_penalty
                    
                    if not will_pass:
                        failed_primers.append({
                            'gene': gene,
                            'penalty': penalty,
                        })
                    
                    if DebugLogLimiter.should_log('penalty_filter_analysis', interval=500, max_initial=2):
                        penalty_f = row.get("Penalty (F)", "N/A")
                        penalty_r = row.get("Penalty (R)", "N/A")
                        status = "PASS" if will_pass else "FAIL"
                        logger.debug(f"Gene {gene}: F={penalty_f}, R={penalty_r}, Overall={penalty} [{status}]")
                        
                except (ValueError, TypeError):
                    invalid_penalties += 1
                    if DebugLogLimiter.should_log('penalty_filter_invalid', interval=200, max_initial=1):
                        logger.debug(f"Gene {gene}: Invalid penalty value - Overall={penalty}")
            
            logger.debug(f"Penalty analysis: {len(failed_primers)} failed, {invalid_penalties} invalid from {initial_count} total")
        
        # Apply penalty filter - simple and clean
        df_filtered = df[df["Penalty"] <= max_penalty].reset_index(drop=True)
        
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"Penalty filtering: kept {filtered_count}, removed {removed_count}")
        
        return df_filtered
    
    @staticmethod
    def filter_by_repeats(primers):
        """
        Filter primers containing disallowed repeat sequences.
        
        Removes primers with problematic repeats (GGGG, CCCC) that can cause
        PCR artifacts. Only examines primer sequences, not probes.
        
        Args:
            primers: DataFrame or list of primer dictionaries
            
        Returns:
            Filtered primers in same format as input
        """
        # Convert to DataFrame if needed
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        initial_count = len(df)
        
        # Check for expected columns - handle both possible naming conventions
        forward_col = None
        reverse_col = None
        
        if "Sequence (F)" in df.columns:
            forward_col = "Sequence (F)"
        elif "Primer F" in df.columns:
            forward_col = "Primer F"
        else:
            error_msg = "No forward primer sequence column found (expected 'Sequence (F)' or 'Primer F')"
            logger.error(error_msg)
            logger.debug(f"Available columns: {df.columns.tolist()}")
            raise ValueError(error_msg)
        
        if "Sequence (R)" in df.columns:
            reverse_col = "Sequence (R)"
        elif "Primer R" in df.columns:
            reverse_col = "Primer R"
        else:
            error_msg = "No reverse primer sequence column found (expected 'Sequence (R)' or 'Primer R')"
            logger.error(error_msg)
            logger.debug(f"Available columns: {df.columns.tolist()}")
            raise ValueError(error_msg)
        
        logger.debug(f"Using columns: forward='{forward_col}', reverse='{reverse_col}'")
        
        # Collect statistics for logging
        if logger.isEnabledFor(logging.DEBUG):
            failed_primers = []
            
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                primer_f = row.get(forward_col, "")
                primer_r = row.get(reverse_col, "")
                
                has_f_repeats = FilterProcessor.has_disallowed_repeats(primer_f)
                has_r_repeats = FilterProcessor.has_disallowed_repeats(primer_r)
                will_fail = has_f_repeats or has_r_repeats
                
                if will_fail:
                    failed_primers.append({
                        'gene': gene,
                        'has_f_repeats': has_f_repeats,
                        'has_r_repeats': has_r_repeats,
                    })
                
                # Limited logging using DebugLogLimiter
                if DebugLogLimiter.should_log('repeat_filter_analysis', interval=50, max_initial=3):
                    if will_fail:
                        logger.debug(f"Gene {gene}: FILTERED due to primer repeats")
                    else:
                        logger.debug(f"Gene {gene}: PASS - no disallowed repeats")
            
            logger.debug(f"Repeat analysis: {len(failed_primers)} with disallowed repeats from {initial_count} total")
        
        # Apply repeat filtering
        df["Has_Repeats_F"] = df[forward_col].apply(FilterProcessor.has_disallowed_repeats)
        df["Has_Repeats_R"] = df[reverse_col].apply(FilterProcessor.has_disallowed_repeats)
        
        df_filtered = df[~(df["Has_Repeats_F"] | df["Has_Repeats_R"])].reset_index(drop=True)
        
        df_filtered = df_filtered.drop(columns=[col for col in df_filtered.columns if col.startswith("Has_Repeats_")])
        
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"Repeat filtering: kept {filtered_count}, removed {removed_count}")
        
        # Return in original format
        return df_filtered.to_dict('records') if not isinstance(primers, pd.DataFrame) else df_filtered

    @staticmethod
    def filter_by_gc_content(primers, min_gc=None, max_gc=None):
        """
        Filter primers by amplicon GC content.
        
        Args:
            primers: DataFrame or list of primer dictionaries
            min_gc: Minimum GC percentage, defaults to Config.SEQUENCE_MIN_GC
            max_gc: Maximum GC percentage, defaults to Config.SEQUENCE_MAX_GC
            
        Returns:
            Filtered primers in same format as input
        """
        if min_gc is None:
            min_gc = Config.SEQUENCE_MIN_GC
        if max_gc is None:
            max_gc = Config.SEQUENCE_MAX_GC
        
        logger.debug(f"GC content range: {min_gc}% - {max_gc}%")
        
        # Convert to DataFrame if needed
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        initial_count = len(df)
        
        # Check for amplicon column - handle both possible naming conventions
        amplicon_col = None
        if "Sequence (A)" in df.columns:
            amplicon_col = "Sequence (A)"
        elif "Amplicon" in df.columns:
            amplicon_col = "Amplicon"
        else:
            error_msg = "No amplicon sequence column found (expected 'Sequence (A)' or 'Amplicon')"
            logger.error(error_msg)
            logger.debug(f"Available columns: {df.columns.tolist()}")
            raise ValueError(error_msg)
        
        logger.debug(f"Using amplicon column: '{amplicon_col}'")
        
        # Check for missing amplicons
        missing_amplicons = df[amplicon_col].isna() | (df[amplicon_col] == "")
        missing_count = missing_amplicons.sum()
        
        if missing_count > 0:
            logger.warning(f"{missing_count} primers have missing amplicon sequences")
        
        # Collect statistics for logging
        if logger.isEnabledFor(logging.DEBUG):
            failed_primers = []
            gc_values = []
            
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                amplicon = row.get(amplicon_col, "")
                
                if not amplicon:
                    continue
                
                gc_content = FilterProcessor.calculate_gc(amplicon)
                gc_values.append(gc_content)
                in_range = min_gc <= gc_content <= max_gc
                
                if not in_range:
                    failed_primers.append({
                        'gene': gene,
                        'gc_content': gc_content,
                        'reason': "too low" if gc_content < min_gc else "too high"
                    })
                
                if DebugLogLimiter.should_log('gc_filter_analysis', interval=500, max_initial=2):
                    status = "PASS" if in_range else "FAIL"
                    logger.debug(f"Gene {gene}: GC={gc_content:.1f}% [{status}]")
            
            if gc_values:
                avg_gc = sum(gc_values) / len(gc_values)
                min_gc_found = min(gc_values)
                max_gc_found = max(gc_values)
                logger.debug(f"GC analysis: avg={avg_gc:.1f}%, range={min_gc_found:.1f}-{max_gc_found:.1f}%, {len(failed_primers)} failed")
        
        # Remove rows with missing amplicons and calculate GC content
        df_with_amplicons = df.dropna(subset=[amplicon_col]).copy()
        df_with_amplicons = df_with_amplicons[df_with_amplicons[amplicon_col] != ""].reset_index(drop=True)
        
        df_with_amplicons["Amplicon GC%"] = df_with_amplicons[amplicon_col].apply(FilterProcessor.calculate_gc)
        df_filtered = df_with_amplicons[
            (df_with_amplicons["Amplicon GC%"] >= min_gc) & 
            (df_with_amplicons["Amplicon GC%"] <= max_gc)
        ].reset_index(drop=True)
        
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"GC content filtering: kept {filtered_count}, removed {removed_count}")
        
        # Return in original format
        return df_filtered.to_dict('records') if not isinstance(primers, pd.DataFrame) else df_filtered

    @staticmethod
    def filter_by_blast(primers, blast_filter_factor=None):
        """
        Filter primers by BLAST specificity results.
        
        Args:
            primers: DataFrame or list of primer dictionaries
            blast_filter_factor: BLAST specificity threshold, defaults to Config.BLAST_FILTER_FACTOR
            
        Returns:
            Filtered primers passing BLAST specificity criteria
        """
        if blast_filter_factor is None:
            blast_filter_factor = Config.BLAST_FILTER_FACTOR
        
        logger.debug(f"BLAST filter factor: {blast_filter_factor}")
        
        # Convert to DataFrame if needed
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        initial_count = len(df)
        
        # Check for required BLAST columns with new naming
        required_cols = ["BLAST (F)", "BLAST (F2)", "BLAST (R)", "BLAST (R2)"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logger.warning(f"Missing required BLAST columns: {missing_cols}")
            logger.warning("Cannot apply BLAST filtering - skipping")
            return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
        
        logger.debug("Using BLAST columns: BLAST (F), BLAST (F2), BLAST (R), BLAST (R2)")
        if "BLAST (P)" in df.columns and "BLAST (P2)" in df.columns:
            logger.debug("Probe BLAST columns also available: BLAST (P), BLAST (P2)")
        
        keep_indices = []
        failed_primers = []
        
        for i, row in df.iterrows():
            gene = row.get("Gene", f"Unknown_{i}")
            
            # Check forward primer specificity using original logic
            f_blast1 = row.get("BLAST (F)")
            f_blast2 = row.get("BLAST (F2)")
            keep_f = FilterProcessor._passes_blast_filter(f_blast1, f_blast2, blast_filter_factor)
            
            # Check reverse primer specificity using original logic
            r_blast1 = row.get("BLAST (R)")
            r_blast2 = row.get("BLAST (R2)")
            keep_r = FilterProcessor._passes_blast_filter(r_blast1, r_blast2, blast_filter_factor)
            
            # Check probe if columns exist
            keep_p = True
            if "BLAST (P)" in df.columns and "BLAST (P2)" in df.columns:
                probe_seq = row.get("Sequence (P)")
                if pd.notnull(probe_seq) and probe_seq:
                    p_blast1 = row.get("BLAST (P)")
                    p_blast2 = row.get("BLAST (P2)")
                    keep_p = FilterProcessor._passes_blast_filter(p_blast1, p_blast2, blast_filter_factor)
            
            overall_pass = keep_f and keep_r and keep_p
            
            if overall_pass:
                keep_indices.append(i)
            else:
                failure_reasons = []
                if not keep_f:
                    failure_reasons.append("forward primer specificity")
                if not keep_r:
                    failure_reasons.append("reverse primer specificity")
                if not keep_p:
                    failure_reasons.append("probe specificity")
                    
                failed_primers.append({
                    'gene': gene,
                    'reasons': failure_reasons,
                })
                
            if DebugLogLimiter.should_log('blast_filter_analysis', interval=500, max_initial=2):
                status = "PASS" if overall_pass else "FAIL"
                logger.debug(f"Gene {gene}: F={keep_f}, R={keep_r}, P={keep_p} -> {status}")
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"BLAST analysis: {len(failed_primers)} failed from {initial_count} total")
        
        # Apply filtering
        df_filtered = df.loc[keep_indices].reset_index(drop=True)
        
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"BLAST filtering: kept {filtered_count}, removed {removed_count}")
        
        # Return in original format
        return df_filtered.to_dict('records') if not isinstance(primers, pd.DataFrame) else df_filtered

    @staticmethod
    def _passes_blast_filter(blast1, blast2, filter_factor):
        """
        Original BLAST filtering logic using two e-values.
        
        Args:
            blast1: Best e-value (lowest)
            blast2: Second best e-value
            filter_factor: Specificity factor from Config.BLAST_FILTER_FACTOR
            
        Returns:
            True if primer passes specificity filter, False otherwise
        """
        # No hits found - insufficient specificity data
        if blast1 is None:
            return False
            
        # Only one hit found - effectively unique (specific)
        if blast2 is None:
            return True
            
        # Check specificity ratio: blast1 * filter_factor <= blast2
        # This means the second-best hit should be significantly worse than the best hit
        return blast1 * filter_factor <= blast2
    

    @staticmethod
    def has_disallowed_repeats(seq):
        """
        Check for disallowed repeats in a DNA sequence.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            True if disallowed repeats found, False otherwise
        """
        if not isinstance(seq, str):
            if DebugLogLimiter.should_log('disallowed_repeats_invalid', interval=200, max_initial=1):
                logger.debug("Invalid sequence type provided to has_disallowed_repeats")
            return True
            
        if not seq:
            return False
            
        seq_upper = seq.upper()
        has_repeats = "CCCC" in seq_upper or "GGGG" in seq_upper
        
        if (has_repeats and DebugLogLimiter.should_log('disallowed_repeats', interval=200, max_initial=3)):
            repeat_types = []
            if "CCCC" in seq_upper:
                repeat_types.append("CCCC")
            if "GGGG" in seq_upper:
                repeat_types.append("GGGG")
            logger.debug(f"Disallowed repeats found ({', '.join(repeat_types)}): {seq[:20]}...")
            
        return has_repeats
    
    @staticmethod
    def calculate_gc(seq):
        """
        Calculate GC content of a DNA sequence.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            GC content as a percentage (0-100)
        """
        if not seq or not isinstance(seq, str):
            if DebugLogLimiter.should_log('gc_calc_invalid', interval=200, max_initial=1):
                logger.debug("Invalid sequence provided to calculate_gc")
            return 0.0
            
        if not seq.strip():
            return 0.0
            
        seq_upper = seq.upper()
        gc_count = sum(1 for base in seq_upper if base in "GC")
        total_bases = len(seq_upper)
        
        if total_bases == 0:
            return 0.0
            
        gc_percentage = (gc_count / total_bases) * 100
        
        return gc_percentage
    
    @staticmethod
    def reverse_complement(seq):
        """
        Generate the reverse complement of a DNA sequence.
        
        Args:
            seq: DNA sequence to reverse complement
            
        Returns:
            Reverse complement sequence
        """
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to reverse_complement")
            return ""
        
        if not seq.strip():
            return ""
        
        seq_upper = seq.upper()
        
        # Define complement mapping including IUPAC ambiguous codes
        complement = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
            'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
            'D': 'H', 'H': 'D', 'V': 'B'
        }
        
        # Check for invalid characters
        invalid_chars = set(seq_upper) - set(complement.keys())
        if invalid_chars:
            error_msg = f"Invalid nucleotide characters found: {', '.join(invalid_chars)}"
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        # Generate reverse complement
        try:
            rev_comp = ''.join(complement.get(base, base) for base in reversed(seq_upper))
            return rev_comp
        except Exception as e:
            error_msg = f"Error generating reverse complement for sequence: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ValueError(error_msg) from e
    
    @staticmethod
    def ensure_more_c_than_g(seq):
        """
        Ensure a sequence has more Cs than Gs, reversing if necessary.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            Tuple of (possibly_reversed_sequence, was_reversed)
        """
        if not seq or not isinstance(seq, str):
            if DebugLogLimiter.should_log('ensure_c_g_invalid', interval=20, max_initial=2):
                logger.debug("Invalid sequence provided to ensure_more_c_than_g")
            return seq, False
        
        if not seq.strip():
            return seq, False
        
        seq_upper = seq.upper()
        c_count = seq_upper.count('C')
        g_count = seq_upper.count('G')
        
        if c_count >= g_count:
            return seq, False
        
        try:
            rev_comp = FilterProcessor.reverse_complement(seq)
            if DebugLogLimiter.should_log('c_g_reversed', interval=200, max_initial=1):
                logger.debug("Sequence reversed to ensure more Cs than Gs")
            return rev_comp, True
        except ValueError as e:
            error_msg = f"Failed to reverse complement sequence: {str(e)}"
            logger.error(error_msg)
            return seq, False
        except Exception as e:
            error_msg = f"Unexpected error in ensure_more_c_than_g: {str(e)}"
            logger.error(error_msg)
            return seq, False