#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer processing module for ddPrimer pipeline.

Contains functionality for:
1. Primer filtering based on various criteria
2. Penalty thresholding and validation
3. GC content analysis and filtering
4. Repeat sequence detection and filtering
5. BLAST specificity evaluation
6. Internal oligo processing and optimization
"""

import pandas as pd
import logging

from ..config import Config
from ..utils import SequenceUtils

# Set up module logger
logger = logging.getLogger(__name__)


class PrimerProcessor:
    """
    Handles primer processing and filtering operations.
    
    This class provides comprehensive filtering capabilities for primer pairs
    including penalty scoring, sequence composition analysis, and specificity
    validation through multiple filtering stages.
    
    Example:
        >>> df = PrimerProcessor.filter_by_penalty(primer_df)
        >>> df = PrimerProcessor.filter_by_gc_content(df)
        >>> df = PrimerProcessor.filter_by_blast(df)
    """
    
    @staticmethod
    def filter_by_penalty(df, max_penalty=None):
        """
        Filter primers by penalty scores.
        
        Removes primer pairs with penalty scores exceeding the specified
        threshold, keeping only the highest-quality designs.
        
        Args:
            df: DataFrame containing primer data with penalty information
            max_penalty: Maximum allowed penalty score, defaults to Config.PENALTY_MAX
            
        Returns:
            Filtered DataFrame containing only primers within penalty threshold
            
        Raises:
            ValueError: If penalty columns are missing or invalid
            
        Example:
            >>> filtered_df = PrimerProcessor.filter_by_penalty(primers_df, max_penalty=2.0)
            >>> print(f"Kept {len(filtered_df)} primers after penalty filtering")
        """
        if max_penalty is None:
            max_penalty = Config.PENALTY_MAX
        
        logger.debug("=== PENALTY FILTER DEBUG ===")
        logger.debug(f"Penalty threshold: {max_penalty}")
        
        if df.empty:
            logger.warning("Empty DataFrame provided for penalty filtering")
            return df
        
        # Ensure Pair Penalty column exists
        if "Pair Penalty" not in df.columns:
            logger.warning("'Pair Penalty' column not found - attempting to create from individual penalties")
            
            if "Penalty F" in df.columns and "Penalty R" in df.columns:
                logger.debug("Creating 'Pair Penalty' from individual F and R penalties")
                df["Pair Penalty"] = df["Penalty F"] + df["Penalty R"]
            else:
                error_msg = "Cannot create 'Pair Penalty' - missing individual penalty columns"
                logger.error(error_msg)
                logger.debug(f"Available columns: {df.columns.tolist()}")
                raise ValueError(error_msg)
        
        initial_count = len(df)
        
        # Log detailed penalty information in debug mode
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Penalty analysis before filtering:")
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                penalty_f = row.get("Penalty F", "N/A")
                penalty_r = row.get("Penalty R", "N/A") 
                penalty_pair = row.get("Pair Penalty", "N/A")
                
                try:
                    will_pass = float(penalty_pair) <= max_penalty
                    status = "PASS" if will_pass else "FAIL"
                    logger.debug(f"Gene {gene}: F={penalty_f}, R={penalty_r}, Pair={penalty_pair} [{status}]")
                    
                    if not will_pass:
                        primer_f = row.get("Primer F", "N/A")
                        primer_r = row.get("Primer R", "N/A")
                        logger.debug(f"  FILTERED: {gene} penalty {penalty_pair} > {max_penalty}")
                        logger.debug(f"  Forward: {primer_f}")
                        logger.debug(f"  Reverse: {primer_r}")
                        
                except (ValueError, TypeError):
                    logger.debug(f"Gene {gene}: Invalid penalty values - F={penalty_f}, R={penalty_r}, Pair={penalty_pair}")
        
        # Apply penalty filter
        df_filtered = df[df["Pair Penalty"] <= max_penalty].reset_index(drop=True)
        
        # Log filtering results
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"Penalty filtering results:")
        logger.debug(f"  Initial primers: {initial_count}")
        logger.debug(f"  Primers kept: {filtered_count}")
        logger.debug(f"  Primers removed: {removed_count}")
        
        if logger.isEnabledFor(logging.DEBUG):
            kept_genes = set(df_filtered["Gene"]) if not df_filtered.empty else set()
            removed_genes = set(df["Gene"]) - kept_genes
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(sorted(kept_genes))}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(sorted(removed_genes))}")
        
        logger.debug("=== END PENALTY FILTER DEBUG ===")
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
            
        Example:
            >>> clean_primers = PrimerProcessor.filter_by_repeats(primer_data)
            >>> print("Removed primers with GGGG or CCCC repeats")
        """
        logger.debug("=== REPEAT FILTER DEBUG ===")
        
        # Convert to DataFrame if needed
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        initial_count = len(df)
        
        # Debug logging for repeat analysis
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Analyzing primers for repeat sequences:")
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                primer_f = row.get("Primer F", "")
                primer_r = row.get("Primer R", "")
                probe = row.get("Probe", "") if "Probe" in df.columns else ""
                
                has_f_repeats = SequenceUtils.has_disallowed_repeats(primer_f)
                has_r_repeats = SequenceUtils.has_disallowed_repeats(primer_r)
                has_probe_repeats = SequenceUtils.has_disallowed_repeats(probe) if probe else False
                
                will_fail = has_f_repeats or has_r_repeats
                
                if will_fail:
                    logger.debug(f"Gene {gene}: WILL BE FILTERED due to primer repeats")
                    if has_f_repeats:
                        repeat_positions = []
                        if "GGGG" in primer_f:
                            repeat_positions.append(f"GGGG at {primer_f.find('GGGG')}")
                        if "CCCC" in primer_f:
                            repeat_positions.append(f"CCCC at {primer_f.find('CCCC')}")
                        logger.debug(f"  Forward primer: {primer_f} ({'; '.join(repeat_positions)})")
                    
                    if has_r_repeats:
                        repeat_positions = []
                        if "GGGG" in primer_r:
                            repeat_positions.append(f"GGGG at {primer_r.find('GGGG')}")
                        if "CCCC" in primer_r:
                            repeat_positions.append(f"CCCC at {primer_r.find('CCCC')}")
                        logger.debug(f"  Reverse primer: {primer_r} ({'; '.join(repeat_positions)})")
                    
                    # Note probe repeats but don't filter based on them
                    if has_probe_repeats:
                        logger.debug(f"  NOTE: Probe has repeats but will NOT cause filtering: {probe}")
                else:
                    logger.debug(f"Gene {gene}: PASS - no disallowed repeats in primers")
        
        # Apply repeat filtering (primers only, not probes)
        df["Has_Repeats_F"] = df["Primer F"].apply(SequenceUtils.has_disallowed_repeats)
        df["Has_Repeats_R"] = df["Primer R"].apply(SequenceUtils.has_disallowed_repeats)
        
        df_filtered = df[~(df["Has_Repeats_F"] | df["Has_Repeats_R"])].reset_index(drop=True)
        
        # Clean up temporary columns
        df_filtered = df_filtered.drop(columns=[col for col in df_filtered.columns if col.startswith("Has_Repeats_")])
        
        # Log filtering results
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"Repeat filtering results:")
        logger.debug(f"  Initial primers: {initial_count}")
        logger.debug(f"  Primers kept: {filtered_count}")
        logger.debug(f"  Primers removed: {removed_count}")
        
        if logger.isEnabledFor(logging.DEBUG):
            kept_genes = set(df_filtered["Gene"]) if not df_filtered.empty else set()
            removed_genes = set(df["Gene"]) - kept_genes
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(sorted(kept_genes))}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(sorted(removed_genes))}")
        
        logger.debug("=== END REPEAT FILTER DEBUG ===")
        
        # Return in original format
        return df_filtered.to_dict('records') if not isinstance(primers, pd.DataFrame) else df_filtered
    
    @staticmethod
    def filter_by_gc_content(primers, min_gc=None, max_gc=None):
        """
        Filter primers by amplicon GC content.
        
        Removes primers with amplicons outside the acceptable GC content
        range to ensure optimal PCR performance.
        
        Args:
            primers: DataFrame or list of primer dictionaries
            min_gc: Minimum GC percentage, defaults to Config.SEQUENCE_MIN_GC
            max_gc: Maximum GC percentage, defaults to Config.SEQUENCE_MAX_GC
            
        Returns:
            Filtered primers in same format as input
            
        Raises:
            ValueError: If amplicon sequences are missing
            
        Example:
            >>> filtered = PrimerProcessor.filter_by_gc_content(primers, min_gc=40, max_gc=60)
            >>> print(f"Kept primers with 40-60% GC content")
        """
        if min_gc is None:
            min_gc = Config.SEQUENCE_MIN_GC
        if max_gc is None:
            max_gc = Config.SEQUENCE_MAX_GC
        
        logger.debug("=== GC CONTENT FILTER DEBUG ===")
        logger.debug(f"GC content range: {min_gc}% - {max_gc}%")
        
        # Convert to DataFrame if needed
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        initial_count = len(df)
        
        # Check for missing amplicon sequences
        missing_amplicons = df["Amplicon"].isna() | (df["Amplicon"] == "")
        if missing_amplicons.any():
            missing_count = missing_amplicons.sum()
            logger.warning(f"{missing_count} primers have missing amplicon sequences")
            
            if logger.isEnabledFor(logging.DEBUG):
                for _, row in df[missing_amplicons].iterrows():
                    gene = row.get("Gene", "Unknown")
                    logger.debug(f"Missing amplicon: Gene={gene}, F={row.get('Primer F', 'N/A')}, R={row.get('Primer R', 'N/A')}")
        
        # Debug analysis of GC content
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("GC content analysis:")
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                amplicon = row.get("Amplicon", "")
                
                if not amplicon:
                    logger.debug(f"Gene {gene}: WILL BE FILTERED - no amplicon sequence")
                    continue
                
                gc_content = SequenceUtils.calculate_gc(amplicon)
                in_range = min_gc <= gc_content <= max_gc
                
                # Detailed sequence analysis
                g_count = amplicon.upper().count('G')
                c_count = amplicon.upper().count('C')
                total_gc = g_count + c_count
                length = len(amplicon)
                
                status = "PASS" if in_range else "FAIL"
                logger.debug(f"Gene {gene}: GC={gc_content:.1f}% ({total_gc}/{length}) [{status}]")
                
                if not in_range:
                    reason = "too low" if gc_content < min_gc else "too high"
                    logger.debug(f"  FILTERED: GC content {reason} ({gc_content:.1f}% not in {min_gc}-{max_gc}%)")
                    
                    # Log amplicon sequence for detailed inspection
                    if len(amplicon) > 40:
                        logger.debug(f"  Amplicon: {amplicon[:20]}...{amplicon[-20:]} ({length} bp)")
                    else:
                        logger.debug(f"  Amplicon: {amplicon} ({length} bp)")
        
        # Remove rows with missing amplicons
        df_with_amplicons = df.dropna(subset=["Amplicon"]).copy()
        df_with_amplicons = df_with_amplicons[df_with_amplicons["Amplicon"] != ""].reset_index(drop=True)
        
        # Calculate and filter by GC content
        df_with_amplicons["Amplicon GC%"] = df_with_amplicons["Amplicon"].apply(SequenceUtils.calculate_gc)
        df_filtered = df_with_amplicons[
            (df_with_amplicons["Amplicon GC%"] >= min_gc) & 
            (df_with_amplicons["Amplicon GC%"] <= max_gc)
        ].reset_index(drop=True)
        
        # Log filtering results
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"GC content filtering results:")
        logger.debug(f"  Initial primers: {initial_count}")
        logger.debug(f"  Primers kept: {filtered_count}")
        logger.debug(f"  Primers removed: {removed_count}")
        
        if logger.isEnabledFor(logging.DEBUG):
            kept_genes = set(df_filtered["Gene"]) if not df_filtered.empty else set()
            removed_genes = set(df["Gene"]) - kept_genes
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(sorted(kept_genes))}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(sorted(removed_genes))}")
        
        logger.debug("=== END GC CONTENT FILTER DEBUG ===")
        
        # Return in original format
        return df_filtered.to_dict('records') if not isinstance(primers, pd.DataFrame) else df_filtered
    
    @staticmethod
    def process_internal_oligos(primers):
        """
        Process internal oligos with reverse complementation optimization.
        
        Reverse complements internal oligos (probes) that have more G than C
        to optimize binding characteristics and probe performance.
        
        Args:
            primers: DataFrame or list of primer dictionaries
            
        Returns:
            Processed primers with optimized probe orientations
            
        Example:
            >>> processed = PrimerProcessor.process_internal_oligos(primer_data)
            >>> # Check which probes were reversed
            >>> reversed_probes = processed[processed['Probe Reversed'] == True]
        """
        # Convert to DataFrame if needed
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # Skip processing if no probes present
        if "Probe" not in df.columns:
            logger.debug("No internal oligos/probes found - skipping processing")
            return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
        
        logger.debug("Processing internal oligos for optimal C>G ratio")
        processed_count = 0
        reversed_count = 0
        
        # Process each probe
        for i, row in df.iterrows():
            probe_seq = row.get("Probe")
            
            if pd.isnull(probe_seq) or not probe_seq:
                df.at[i, "Probe Reversed"] = False
                continue
            
            # Apply reverse complementation if needed
            optimized_probe, was_reversed = SequenceUtils.ensure_more_c_than_g(probe_seq)
            
            df.at[i, "Probe"] = optimized_probe
            df.at[i, "Probe Reversed"] = was_reversed
            
            processed_count += 1
            if was_reversed:
                reversed_count += 1
                logger.debug(f"Reversed probe for better C>G ratio: {probe_seq} -> {optimized_probe}")
        
        logger.debug(f"Internal oligo processing completed: {processed_count} probes, {reversed_count} reversed")
        
        # Return in original format
        return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
    
    @staticmethod
    def filter_by_blast(primers, blast_filter_factor=None):
        """
        Filter primers by BLAST specificity results.
        
        Removes primers that fail specificity criteria based on the ratio
        between best and second-best BLAST e-values.
        
        Args:
            primers: DataFrame or list of primer dictionaries
            blast_filter_factor: BLAST specificity threshold, defaults to Config.BLAST_FILTER_FACTOR
            
        Returns:
            Filtered primers passing BLAST specificity criteria
            
        Example:
            >>> specific_primers = PrimerProcessor.filter_by_blast(blast_results)
            >>> print(f"Kept {len(specific_primers)} specific primers")
        """
        if blast_filter_factor is None:
            blast_filter_factor = Config.BLAST_FILTER_FACTOR
        
        logger.debug("=== BLAST FILTER DEBUG ===")
        logger.debug(f"BLAST filter factor: {blast_filter_factor}")
        
        # Convert to DataFrame if needed
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        initial_count = len(df)
        
        # Check for required BLAST columns
        required_cols = ["Primer F BLAST1", "Primer F BLAST2", "Primer R BLAST1", "Primer R BLAST2"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logger.warning(f"Missing BLAST columns: {missing_cols} - skipping BLAST filtering")
            return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
        
        # Debug analysis of BLAST results
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("BLAST specificity analysis:")
            
        keep_indices = []
        
        for i, row in df.iterrows():
            gene = row.get("Gene", f"Unknown_{i}")
            
            # Check forward primer specificity
            keep_f = PrimerProcessor._passes_blast_filter(row, "Primer F", blast_filter_factor)
            
            # Check reverse primer specificity  
            keep_r = PrimerProcessor._passes_blast_filter(row, "Primer R", blast_filter_factor)
            
            # Check probe specificity if present
            if "Probe" in df.columns and "Probe BLAST1" in df.columns:
                probe_seq = row.get("Probe")
                keep_p = (pd.isna(probe_seq) or probe_seq == "" or 
                         PrimerProcessor._passes_blast_filter(row, "Probe", blast_filter_factor))
            else:
                keep_p = True
            
            # Log detailed BLAST analysis
            if logger.isEnabledFor(logging.DEBUG):
                f_blast1 = row.get("Primer F BLAST1", "N/A")
                f_blast2 = row.get("Primer F BLAST2", "N/A") 
                r_blast1 = row.get("Primer R BLAST1", "N/A")
                r_blast2 = row.get("Primer R BLAST2", "N/A")
                
                logger.debug(f"Gene {gene}:")
                logger.debug(f"  Forward: Best={f_blast1}, Second={f_blast2}, Pass={keep_f}")
                logger.debug(f"  Reverse: Best={r_blast1}, Second={r_blast2}, Pass={keep_r}")
                
                if "Probe" in df.columns and "Probe BLAST1" in df.columns and not pd.isna(row.get("Probe")):
                    p_blast1 = row.get("Probe BLAST1", "N/A")
                    p_blast2 = row.get("Probe BLAST2", "N/A")
                    logger.debug(f"  Probe: Best={p_blast1}, Second={p_blast2}, Pass={keep_p}")
            
            # Keep primer if all components pass
            overall_pass = keep_f and keep_r and keep_p
            
            if overall_pass:
                keep_indices.append(i)
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(f"  Result: PASS - primer kept")
            else:
                if logger.isEnabledFor(logging.DEBUG):
                    reasons = []
                    if not keep_f:
                        reasons.append("forward primer specificity")
                    if not keep_r:
                        reasons.append("reverse primer specificity")
                    if not keep_p:
                        reasons.append("probe specificity")
                    logger.debug(f"  Result: FAIL - filtered due to {', '.join(reasons)}")
        
        # Apply filtering
        df_filtered = df.loc[keep_indices].reset_index(drop=True)
        
        # Log filtering results
        filtered_count = len(df_filtered)
        removed_count = initial_count - filtered_count
        
        logger.debug(f"BLAST filtering results:")
        logger.debug(f"  Initial primers: {initial_count}")
        logger.debug(f"  Primers kept: {filtered_count}")
        logger.debug(f"  Primers removed: {removed_count}")
        
        if logger.isEnabledFor(logging.DEBUG):
            kept_genes = set(df_filtered["Gene"]) if not df_filtered.empty else set()
            removed_genes = set(df["Gene"]) - kept_genes
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(sorted(kept_genes))}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(sorted(removed_genes))}")
        
        logger.debug("=== END BLAST FILTER DEBUG ===")
        
        # Return in original format
        return df_filtered.to_dict('records') if not isinstance(primers, pd.DataFrame) else df_filtered
    
    @staticmethod
    def _passes_blast_filter(row, col_prefix, filter_factor):
        """
        Check if a sequence passes the BLAST specificity filter.
        
        Evaluates the ratio between best and second-best e-values against
        the specificity threshold to determine if a sequence is sufficiently unique.
        
        Args:
            row: DataFrame row containing BLAST results
            col_prefix: Column prefix for BLAST data (e.g., "Primer F")
            filter_factor: Minimum specificity ratio required
            
        Returns:
            True if sequence passes specificity filter, False otherwise
            
        Example:
            >>> passes = PrimerProcessor._passes_blast_filter(row, "Primer F", 100.0)
        """
        best = row.get(f"{col_prefix} BLAST1")
        second = row.get(f"{col_prefix} BLAST2")
        
        # No hits found - insufficient specificity data
        if pd.isna(best):
            return False
            
        # Only one hit found - effectively unique
        if pd.isna(second):
            return True
            
        # Check specificity ratio: best * filter_factor <= second
        return best * filter_factor <= second