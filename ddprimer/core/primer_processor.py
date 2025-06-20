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
7. Generic primer record creation and validation
8. Amplicon extraction and coordinate handling
9. Sequence utility functions (GC content, reverse complement, etc.)

This module now includes all generic primer processing operations and sequence utilities,
providing a comprehensive toolkit for primer analysis and manipulation.
"""

import pandas as pd
import logging
from typing import Dict, List, Optional, Tuple

from ..config import Config

# Set up module logger
logger = logging.getLogger(__name__)


class PrimerProcessor:
    """
    Handles primer processing and filtering operations.
    
    This class provides comprehensive filtering capabilities for primer pairs
    including penalty scoring, sequence composition analysis, specificity
    validation, amplicon extraction, and sequence utility functions.
    
    Example:
        >>> df = PrimerProcessor.filter_by_penalty(primer_df)
        >>> df = PrimerProcessor.filter_by_gc_content(df)
        >>> df = PrimerProcessor.filter_by_blast(df)
        >>> amplicon = PrimerProcessor.get_amplicon(seq, start, len, start2, len2)
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
        
        # Optimized logging - collect statistics instead of logging each primer
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Analyzing penalty distribution...")
            failed_primers = []
            invalid_penalties = 0
            
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                penalty_pair = row.get("Pair Penalty", "N/A")
                
                try:
                    will_pass = float(penalty_pair) <= max_penalty
                    
                    if not will_pass:
                        failed_primers.append({
                            'gene': gene,
                            'penalty': penalty_pair,
                            'penalty_f': row.get("Penalty F", "N/A"),
                            'penalty_r': row.get("Penalty R", "N/A"),
                            'primer_f': row.get("Primer F", "N/A"),
                            'primer_r': row.get("Primer R", "N/A")
                        })
                    
                    # Sample logging - only log every 100th primer OR first 5 failures
                    if (index % 100 == 0 or 
                        (not will_pass and len(failed_primers) <= 5)):
                        penalty_f = row.get("Penalty F", "N/A")
                        penalty_r = row.get("Penalty R", "N/A")
                        status = "PASS" if will_pass else "FAIL"
                        logger.debug(f"Gene {gene}: F={penalty_f}, R={penalty_r}, Pair={penalty_pair} [{status}]")
                        
                        if not will_pass and len(failed_primers) <= 5:
                            logger.debug(f"  FILTERED: {gene} penalty {penalty_pair} > {max_penalty}")
                            logger.debug(f"  Forward: {row.get('Primer F', 'N/A')}")
                            logger.debug(f"  Reverse: {row.get('Primer R', 'N/A')}")
                        
                except (ValueError, TypeError):
                    invalid_penalties += 1
                    if invalid_penalties <= 5:  # Only log first 5 invalid penalties
                        logger.debug(f"Gene {gene}: Invalid penalty values - Pair={penalty_pair}")
            
            # Summary logging
            logger.debug(f"Penalty analysis summary:")
            logger.debug(f"  Total primers analyzed: {initial_count}")
            logger.debug(f"  Primers failing penalty filter: {len(failed_primers)}")
            logger.debug(f"  Invalid penalty values: {invalid_penalties}")
            
            if len(failed_primers) > 5:
                logger.debug(f"  (Detailed logging shown for first 5 failures only)")
        
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
            # Limit gene list display
            kept_sample = sorted(list(kept_genes)[:10])
            removed_sample = sorted(list(removed_genes)[:10])
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(kept_sample)}{'...' if len(kept_genes) > 10 else ''}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(removed_sample)}{'...' if len(removed_genes) > 10 else ''}")
        
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
        
        # Optimized logging - collect statistics instead of logging each primer
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Analyzing primers for repeat sequences...")
            failed_primers = []
            probe_repeats_count = 0
            
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                primer_f = row.get("Primer F", "")
                primer_r = row.get("Primer R", "")
                probe = row.get("Probe", "") if "Probe" in df.columns else ""
                
                has_f_repeats = PrimerProcessor.has_disallowed_repeats(primer_f)
                has_r_repeats = PrimerProcessor.has_disallowed_repeats(primer_r)
                has_probe_repeats = PrimerProcessor.has_disallowed_repeats(probe) if probe else False
                
                will_fail = has_f_repeats or has_r_repeats
                
                if will_fail:
                    repeat_info = {
                        'gene': gene,
                        'has_f_repeats': has_f_repeats,
                        'has_r_repeats': has_r_repeats,
                        'primer_f': primer_f,
                        'primer_r': primer_r
                    }
                    failed_primers.append(repeat_info)
                
                if has_probe_repeats:
                    probe_repeats_count += 1
                
                # Sample logging - only log every 100th primer OR first 5 failures
                if (index % 100 == 0 or 
                    (will_fail and len(failed_primers) <= 5)):
                    
                    if will_fail and len(failed_primers) <= 5:
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
                        
                        if has_probe_repeats:
                            logger.debug(f"  NOTE: Probe has repeats but will NOT cause filtering: {probe}")
                    
                    elif index % 100 == 0 and not will_fail:
                        logger.debug(f"Gene {gene}: PASS - no disallowed repeats in primers")
            
            # Summary logging
            logger.debug(f"Repeat analysis summary:")
            logger.debug(f"  Total primers analyzed: {initial_count}")
            logger.debug(f"  Primers with disallowed repeats: {len(failed_primers)}")
            logger.debug(f"  Probes with repeats (not filtered): {probe_repeats_count}")
            
            if len(failed_primers) > 5:
                logger.debug(f"  (Detailed logging shown for first 5 failures only)")
        
        # Apply repeat filtering (primers only, not probes)
        df["Has_Repeats_F"] = df["Primer F"].apply(PrimerProcessor.has_disallowed_repeats)
        df["Has_Repeats_R"] = df["Primer R"].apply(PrimerProcessor.has_disallowed_repeats)
        
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
            # Limit gene list display
            kept_sample = sorted(list(kept_genes)[:10])
            removed_sample = sorted(list(removed_genes)[:10])
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(kept_sample)}{'...' if len(kept_genes) > 10 else ''}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(removed_sample)}{'...' if len(removed_genes) > 10 else ''}")
        
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
        
        # Check for missing amplicons
        missing_amplicons = df["Amplicon"].isna() | (df["Amplicon"] == "")
        missing_count = missing_amplicons.sum()
        
        if missing_count > 0:
            logger.warning(f"{missing_count} primers have missing amplicon sequences")
            
            # Only show details for first few missing amplicons
            if logger.isEnabledFor(logging.DEBUG):
                missing_samples = []
                for _, row in df[missing_amplicons].iterrows():
                    if len(missing_samples) < 5:
                        gene = row.get("Gene", "Unknown")
                        missing_samples.append({
                            'gene': gene,
                            'primer_f': row.get("Primer F", "N/A"),
                            'primer_r': row.get("Primer R", "N/A")
                        })
                
                logger.debug("Sample missing amplicons:")
                for sample in missing_samples:
                    logger.debug(f"  Gene={sample['gene']}, F={sample['primer_f']}, R={sample['primer_r']}")
                
                if missing_count > 5:
                    logger.debug(f"  ... and {missing_count - 5} more primers with missing amplicons")
        
        # Optimized logging for GC content analysis
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Analyzing GC content distribution...")
            failed_primers = []
            gc_values = []
            
            for index, row in df.iterrows():
                gene = row.get("Gene", f"Unknown_{index}")
                amplicon = row.get("Amplicon", "")
                
                if not amplicon:
                    continue
                
                gc_content = PrimerProcessor.calculate_gc(amplicon)
                gc_values.append(gc_content)
                in_range = min_gc <= gc_content <= max_gc
                
                if not in_range:
                    failed_primers.append({
                        'gene': gene,
                        'gc_content': gc_content,
                        'amplicon': amplicon,
                        'reason': "too low" if gc_content < min_gc else "too high"
                    })
                
                # Sample logging - only log every 100th primer OR first 5 failures
                if (index % 100 == 0 or 
                    (not in_range and len(failed_primers) <= 5)):
                    
                    g_count = amplicon.upper().count('G')
                    c_count = amplicon.upper().count('C')
                    total_gc = g_count + c_count
                    length = len(amplicon)
                    
                    status = "PASS" if in_range else "FAIL"
                    logger.debug(f"Gene {gene}: GC={gc_content:.1f}% ({total_gc}/{length}) [{status}]")
                    
                    if not in_range and len(failed_primers) <= 5:
                        reason = "too low" if gc_content < min_gc else "too high"
                        logger.debug(f"  FILTERED: GC content {reason} ({gc_content:.1f}% not in {min_gc}-{max_gc}%)")
                        
                        # Log amplicon sequence for detailed inspection
                        if len(amplicon) > 40:
                            logger.debug(f"  Amplicon: {amplicon[:20]}...{amplicon[-20:]} ({length} bp)")
                        else:
                            logger.debug(f"  Amplicon: {amplicon} ({length} bp)")
            
            # Summary statistics
            if gc_values:
                avg_gc = sum(gc_values) / len(gc_values)
                min_gc_found = min(gc_values)
                max_gc_found = max(gc_values)
                logger.debug(f"GC content analysis summary:")
                logger.debug(f"  Amplicons analyzed: {len(gc_values)}")
                logger.debug(f"  GC content range found: {min_gc_found:.1f}% - {max_gc_found:.1f}%")
                logger.debug(f"  Average GC content: {avg_gc:.1f}%")
                logger.debug(f"  Primers failing GC filter: {len(failed_primers)}")
                
                if len(failed_primers) > 5:
                    logger.debug(f"  (Detailed logging shown for first 5 failures only)")
        
        # Remove rows with missing amplicons
        df_with_amplicons = df.dropna(subset=["Amplicon"]).copy()
        df_with_amplicons = df_with_amplicons[df_with_amplicons["Amplicon"] != ""].reset_index(drop=True)
        
        # Calculate and filter by GC content
        df_with_amplicons["Amplicon GC%"] = df_with_amplicons["Amplicon"].apply(PrimerProcessor.calculate_gc)
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
            # Limit gene list display
            kept_sample = sorted(list(kept_genes)[:10])
            removed_sample = sorted(list(removed_genes)[:10])
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(kept_sample)}{'...' if len(kept_genes) > 10 else ''}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(removed_sample)}{'...' if len(removed_genes) > 10 else ''}")
        
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
            optimized_probe, was_reversed = PrimerProcessor.ensure_more_c_than_g(probe_seq)
            
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
        
        # Optimized logging for BLAST analysis
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Analyzing BLAST specificity results...")
            failed_primers = []
            no_hits_count = 0
            unique_hits_count = 0
        
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
            
            # Overall pass/fail
            overall_pass = keep_f and keep_r and keep_p
            
            # Collect statistics for summary
            if logger.isEnabledFor(logging.DEBUG):
                f_blast1 = row.get("Primer F BLAST1")
                f_blast2 = row.get("Primer F BLAST2")
                r_blast1 = row.get("Primer R BLAST1")
                r_blast2 = row.get("Primer R BLAST2")
                
                # Count no hits and unique hits
                if pd.isna(f_blast1) and pd.isna(r_blast1):
                    no_hits_count += 1
                elif (pd.isna(f_blast2) and pd.notna(f_blast1)) or (pd.isna(r_blast2) and pd.notna(r_blast1)):
                    unique_hits_count += 1
                
                if not overall_pass:
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
                        'f_blast1': f_blast1,
                        'f_blast2': f_blast2,
                        'r_blast1': r_blast1,
                        'r_blast2': r_blast2
                    })
                
                # Sample logging - only log every 100th primer OR first 5 failures
                if (i % 100 == 0 or 
                    (not overall_pass and len(failed_primers) <= 5)):
                    
                    logger.debug(f"Gene {gene}:")
                    logger.debug(f"  Forward: Best={f_blast1}, Second={f_blast2}, Pass={keep_f}")
                    logger.debug(f"  Reverse: Best={r_blast1}, Second={r_blast2}, Pass={keep_r}")
                    
                    if "Probe" in df.columns and "Probe BLAST1" in df.columns and not pd.isna(row.get("Probe")):
                        p_blast1 = row.get("Probe BLAST1", "N/A")
                        p_blast2 = row.get("Probe BLAST2", "N/A")
                        logger.debug(f"  Probe: Best={p_blast1}, Second={p_blast2}, Pass={keep_p}")
                    
                    if overall_pass:
                        logger.debug(f"  Result: PASS - primer kept")
                    elif len(failed_primers) <= 5:
                        logger.debug(f"  Result: FAIL - filtered due to {', '.join(failed_primers[-1]['reasons'])}")
            
            # Keep primer if all components pass
            if overall_pass:
                keep_indices.append(i)
        
        # Summary logging
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"BLAST analysis summary:")
            logger.debug(f"  Total primers analyzed: {initial_count}")
            logger.debug(f"  Primers with no BLAST hits: {no_hits_count}")
            logger.debug(f"  Primers with unique hits: {unique_hits_count}")
            logger.debug(f"  Primers failing BLAST filter: {len(failed_primers)}")
            
            if len(failed_primers) > 5:
                logger.debug(f"  (Detailed logging shown for first 5 failures only)")
        
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
            # Limit gene list display
            kept_sample = sorted(list(kept_genes)[:10])
            removed_sample = sorted(list(removed_genes)[:10])
            logger.debug(f"  Genes kept ({len(kept_genes)}): {', '.join(kept_sample)}{'...' if len(kept_genes) > 10 else ''}")
            logger.debug(f"  Genes removed ({len(removed_genes)}): {', '.join(removed_sample)}{'...' if len(removed_genes) > 10 else ''}")
        
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

    @staticmethod
    def get_amplicon(seq, left_start, left_len, right_start, right_len):
        """
        Extract amplicon sequence with corrected coordinate handling.
        
        Returns the substring from the 5' end of the forward primer
        to the 3' end of the reverse primer, addressing systematic
        coordinate alignment issues with Primer3.
        
        Args:
            seq: Template sequence
            left_start: Start position of forward primer (1-based)
            left_len: Length of forward primer
            right_start: Start position of reverse primer (1-based)
            right_len: Length of reverse primer
            
        Returns:
            Amplicon sequence string, empty if extraction fails
            
        Raises:
            ValueError: If coordinate parameters are invalid
            
        Example:
            >>> amplicon = PrimerProcessor.get_amplicon(template, 10, 20, 200, 20)
            >>> print(f"Amplicon length: {len(amplicon)} bp")
        """
        # Validate input parameters
        if None in (left_start, left_len, right_start, right_len):
            logger.warning("Missing required parameters for amplicon extraction")
            return ""
        
        try:
            left_start = int(left_start)
            left_len = int(left_len)
            right_start = int(right_start)
            right_len = int(right_len)
        except (TypeError, ValueError) as e:
            error_msg = f"Invalid numerical values for amplicon extraction"
            logger.error(error_msg)
            logger.debug(f"Parameter conversion error: {str(e)}", exc_info=True)
            raise ValueError(error_msg) from e
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Amplicon extraction: left={left_start},{left_len}, right={right_start},{right_len}")
            
        # Adjust coordinates to address off-by-one issues
        adj_left_start = max(1, left_start - 1)
        
        # Handle edge cases
        if right_start > len(seq):
            logger.debug(f"Adjusting right_start from {right_start} to {len(seq)} (sequence length)")
            right_start = len(seq)
        
        # Calculate amplicon boundaries
        amp_start = adj_left_start
        amp_end = right_start
        
        # Validate coordinates
        if amp_start < 1 or amp_end > len(seq) or amp_start > amp_end:
            logger.warning(f"Invalid amplicon coordinates: start={amp_start}, end={amp_end}, seq_len={len(seq)}")
            return ""
        
        # Extract amplicon sequence (convert to 0-based for Python)
        amplicon = seq[amp_start - 1 : amp_end]
        
        # Debug validation in debug mode
        if logger.isEnabledFor(logging.DEBUG):
            PrimerProcessor._validate_amplicon_extraction(seq, left_start, left_len, right_start, right_len, 
                                             adj_left_start, amplicon)
        
        return amplicon

    @staticmethod
    def _validate_amplicon_extraction(seq, left_start, left_len, right_start, right_len, 
                                    adj_left_start, amplicon):
        """
        Validate amplicon extraction by comparing primer sequences.
        
        Args:
            seq: Template sequence
            left_start: Original left start position
            left_len: Left primer length
            right_start: Right start position
            right_len: Right primer length
            adj_left_start: Adjusted left start position
            amplicon: Extracted amplicon sequence
        """
        # Extract primer sequences for validation
        forward_original = seq[left_start - 1 : left_start - 1 + left_len]
        forward_adjusted = seq[adj_left_start - 1 : adj_left_start - 1 + left_len]
        reverse_primer = seq[right_start - right_len : right_start]
        
        logger.debug(f"Primer validation:")
        logger.debug(f"  Forward (original pos): {forward_original}")
        logger.debug(f"  Forward (adjusted pos): {forward_adjusted}")
        logger.debug(f"  Amplicon start: {amplicon[:left_len] if len(amplicon) >= left_len else amplicon}")
        logger.debug(f"  Reverse primer: {reverse_primer}")
        logger.debug(f"  Amplicon end: {amplicon[-right_len:] if len(amplicon) >= right_len else amplicon}")
        
        # Check alignment
        if len(amplicon) >= left_len and amplicon[:left_len] != forward_adjusted:
            logger.debug("WARNING: Amplicon start does not match adjusted forward primer")
        
        if len(amplicon) >= right_len and amplicon[-right_len:] != reverse_primer:
            logger.debug("WARNING: Amplicon end does not match reverse primer")

    @staticmethod
    def filter_primer_pairs(primer_pairs: List[Dict], config) -> List[Dict]:
        """
        Filter and limit primer pairs based on configuration.
        
        Args:
            primer_pairs: List of primer pairs to filter
            config: Configuration object with filtering parameters
            
        Returns:
            Filtered and limited primer pairs
        """
        acceptable = primer_pairs
        
        # Apply maximum pairs per segment limit
        if (hasattr(config, 'MAX_PRIMER_PAIRS_PER_SEGMENT') and 
            config.MAX_PRIMER_PAIRS_PER_SEGMENT > 0):
            
            original_count = len(acceptable)
            acceptable = acceptable[:config.MAX_PRIMER_PAIRS_PER_SEGMENT]
            
            if logger.isEnabledFor(logging.DEBUG) and original_count > config.MAX_PRIMER_PAIRS_PER_SEGMENT:
                logger.debug(f"Limited to {config.MAX_PRIMER_PAIRS_PER_SEGMENT} primer pairs "
                           f"(had {original_count} total)")
        
        return acceptable

    @staticmethod
    def create_primer_record(block: Dict, pair: Dict, fragment_info: Dict, config, debug_mode: bool = False) -> Optional[Dict]:
        """
        Create a primer record dictionary from sequence block and primer pair.
        
        Args:
            block: Sequence block containing template
            pair: Primer pair data
            fragment_info: Fragment information mapping
            config: Configuration object
            debug_mode: Whether debug logging is enabled
            
        Returns:
            Complete primer record dictionary or None if creation fails
            
        Raises:
            ValueError: If required primer data is missing
        """
        try:
            # Extract basic primer sequences
            left_seq = pair.get("left_sequence", "")
            right_seq = pair.get("right_sequence", "")
            
            if not left_seq or not right_seq:
                logger.warning(f"Missing primer sequences for pair {pair.get('idx', 'unknown')}")
                return None
            
            # Process probe sequence
            probe_seq, probe_reversed = PrimerProcessor._process_probe_sequence(pair, config)
            
            # Get position data
            ls, ll = pair.get("left_start"), pair.get("left_len")
            rs, rl = pair.get("right_start"), pair.get("right_len")
            
            # Create amplicon
            amplicon = PrimerProcessor._create_amplicon(block['sequence_template'], pair, debug_mode)
            
            # Get location information
            location_info = PrimerProcessor._get_location_info(block['sequence_id'], pair, fragment_info)
            
            # Calculate product size
            product_size = pair.get("product_size")
            if product_size is None and amplicon:
                product_size = len(amplicon)
            
            # Build comprehensive primer record
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
            if not config.DISABLE_INTERNAL_OLIGO:
                internal_start = pair.get("internal_start")
                internal_len = pair.get("internal_len")
                
                record.update({
                    "Probe": probe_seq,
                    "Probe Tm": pair.get("internal_tm"),
                    "Probe Penalty": pair.get("internal_penalty"),
                    "Probe Start": internal_start,
                    "Probe Len": internal_len,
                    "Probe Reversed": probe_reversed
                })
            
            return record
            
        except Exception as e:
            error_msg = f"Failed to create primer record for pair {pair.get('idx', 'unknown')}"
            logger.error(error_msg)
            logger.debug(f"Record creation error: {str(e)}", exc_info=True)
            return None

    @staticmethod
    def _process_probe_sequence(pair: Dict, config) -> Tuple[str, bool]:
        """
        Process probe sequence with optional reverse complementation.
        
        Args:
            pair: Primer pair containing probe data
            config: Configuration object
            
        Returns:
            Tuple of (probe_sequence, was_reversed)
        """
        if config.DISABLE_INTERNAL_OLIGO:
            return "", False
        
        probe_seq = pair.get("internal_sequence", "")
        probe_reversed = False
        
        # Apply reverse complementation if configured and probe exists
        if config.PREFER_PROBE_MORE_C_THAN_G and probe_seq:
            try:
                probe_seq, probe_reversed = PrimerProcessor.ensure_more_c_than_g(probe_seq)
                
                if probe_reversed and logger.isEnabledFor(logging.DEBUG):
                    logger.debug(f"Reversed probe for optimal C>G ratio")
                    
            except Exception as e:
                logger.warning(f"Failed to process probe sequence: {str(e)}")
                logger.debug(f"Probe processing error: {str(e)}", exc_info=True)
        
        return probe_seq, probe_reversed

    @staticmethod
    def _create_amplicon(template: str, pair: Dict, debug_mode: bool) -> str:
        """
        Create amplicon sequence from template and primer positions.
        
        Args:
            template: Template sequence
            pair: Primer pair with position information
            debug_mode: Whether debug logging is enabled
            
        Returns:
            Amplicon sequence string, empty if creation fails
        """
        ls, ll = pair.get("left_start"), pair.get("left_len")
        rs, rl = pair.get("right_start"), pair.get("right_len")
        
        # Validate position data
        if not all([ls, ll, rs, rl]):
            if debug_mode:
                logger.debug(f"Missing position data for amplicon creation in pair {pair.get('idx', 'unknown')}")
                logger.debug(f"  Left: {ls},{ll}, Right: {rs},{rl}")
            return ""
        
        try:
            # Use existing amplicon extraction method
            amplicon = PrimerProcessor.get_amplicon(template, ls, ll, rs, rl)
            
            # Detailed debug logging
            if debug_mode:
                if not amplicon:
                    logger.debug(f"WARNING: Failed to create amplicon for pair {pair.get('idx', 'unknown')}")
                    logger.debug(f"  Coordinates: left={ls},{ll}, right={rs},{rl}")
                    logger.debug(f"  Template length: {len(template)}")
                    
                    # Attempt alternative extraction
                    amplicon = PrimerProcessor._try_reconstruct_amplicon(template, ls, rs, debug_mode)
                else:
                    logger.debug(f"Created amplicon for pair {pair.get('idx', 'unknown')}: {len(amplicon)} bp")
                    
                    # Validate primer alignment
                    PrimerProcessor._validate_primer_alignment(pair, amplicon, debug_mode)
            
            return amplicon
            
        except Exception as e:
            error_msg = f"Error creating amplicon for pair {pair.get('idx', 'unknown')}"
            logger.error(error_msg)
            logger.debug(f"Amplicon creation error: {str(e)}", exc_info=True)
            return ""

    @staticmethod
    def _try_reconstruct_amplicon(template: str, ls: int, rs: int, debug_mode: bool) -> str:
        """
        Attempt to reconstruct amplicon when normal extraction fails.
        
        Args:
            template: Template sequence
            ls: Left start position
            rs: Right start position
            debug_mode: Whether debug logging is enabled
            
        Returns:
            Reconstructed amplicon or empty string
        """
        if debug_mode:
            logger.debug("Attempting amplicon reconstruction")
        
        # Try direct template extraction
        if (ls is not None and rs is not None and ls <= rs and 
            ls >= 1 and rs <= len(template)):
            
            direct_amplicon = template[ls-1:rs]
            
            if debug_mode:
                logger.debug(f"Reconstructed amplicon: {len(direct_amplicon)} bp")
                
            if len(direct_amplicon) > 0:
                return direct_amplicon
        
        if debug_mode:
            logger.debug(f"Amplicon reconstruction failed: ls={ls}, rs={rs}, template_len={len(template)}")
        
        return ""

    @staticmethod
    def _validate_primer_alignment(pair: Dict, amplicon: str, debug_mode: bool) -> None:
        """
        Validate that amplicon contains expected primer sequences.
        
        Args:
            pair: Primer pair data
            amplicon: Extracted amplicon sequence
            debug_mode: Whether debug logging is enabled
        """
        if not debug_mode or not amplicon:
            return
        
        left_seq = pair.get("left_sequence", "")
        right_seq = pair.get("right_sequence", "")
        
        # Check forward primer alignment
        if left_seq and len(amplicon) >= len(left_seq):
            if not amplicon.startswith(left_seq[:min(len(left_seq), 10)]):
                logger.debug("WARNING: Amplicon does not start with forward primer")
        
        # Check reverse primer alignment (reverse complement)
        if right_seq and len(amplicon) >= len(right_seq):
            try:
                rc_right = PrimerProcessor.reverse_complement(right_seq)
                if rc_right[:min(len(right_seq), 10)] not in amplicon[-len(right_seq):]:
                    logger.debug("WARNING: Amplicon does not end with reverse primer complement")
            except Exception as e:
                logger.debug(f"Could not validate reverse primer alignment: {e}")

    @staticmethod
    def _get_location_info(sequence_id: str, pair: Dict, fragment_info: Dict) -> Dict[str, str]:
        """
        Get location information for primer record.
        
        Args:
            sequence_id: Sequence identifier
            pair: Primer pair data
            fragment_info: Fragment information mapping
            
        Returns:
            Dictionary with gene, chromosome, and location information
        """
        # Get fragment information
        frag_info = fragment_info.get(sequence_id, {})
        
        # Determine gene name
        gene = frag_info.get("gene", sequence_id)
        
        # Get chromosome
        chromosome = frag_info.get("chr", "")
        
        # Calculate absolute genomic position if possible
        location = ""
        ls = pair.get("left_start")
        if ls is not None and "start" in frag_info:
            try:
                fragment_start = frag_info.get("start", 1)
                abs_left_start = fragment_start + ls - 1
                location = str(abs_left_start)
            except (TypeError, ValueError):
                logger.debug(f"Could not calculate absolute position for {sequence_id}")
        
        return {
            "gene": gene,
            "chromosome": chromosome, 
            "location": location
        }

    # =============================================================================
    # SEQUENCE UTILITY METHODS (moved from sequence_utils.py)
    # =============================================================================

    @staticmethod
    def has_disallowed_repeats(seq):
        """
        Check for disallowed repeats in a DNA sequence.
        
        Identifies sequences containing runs of four or more consecutive
        G or C nucleotides, which can cause issues in PCR amplification
        and primer binding.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            True if disallowed repeats found, False otherwise
            
        Example:
            >>> PrimerProcessor.has_disallowed_repeats("ATCGGGGATC")
            True
            >>> PrimerProcessor.has_disallowed_repeats("ATCGATCG")
            False
        """
        if not isinstance(seq, str):
            logger.debug("Invalid sequence type provided to has_disallowed_repeats")
            return True
            
        if not seq:
            logger.debug("Empty sequence provided to has_disallowed_repeats")
            return False
            
        seq_upper = seq.upper()
        has_repeats = "CCCC" in seq_upper or "GGGG" in seq_upper
        
        if logger.isEnabledFor(logging.DEBUG) and has_repeats:
            logger.debug(f"Disallowed repeats found in sequence: {seq[:20]}...")
            
        return has_repeats
    
    @staticmethod
    def calculate_gc(seq):
        """
        Calculate GC content of a DNA sequence.
        
        Computes the percentage of guanine and cytosine nucleotides
        in the provided DNA sequence, which is important for primer
        design and PCR optimization.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            GC content as a percentage (0-100)
            
        Example:
            >>> PrimerProcessor.calculate_gc("ATCGGCTA")
            37.5
        """
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to calculate_gc")
            return 0.0
            
        if not seq.strip():
            logger.debug("Empty sequence provided to calculate_gc")
            return 0.0
            
        seq_upper = seq.upper()
        gc_count = sum(1 for base in seq_upper if base in "GC")
        total_bases = len(seq_upper)
        
        if total_bases == 0:
            return 0.0
            
        gc_percentage = (gc_count / total_bases) * 100
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"GC content calculation: {gc_count}/{total_bases} = {gc_percentage:.1f}%")
            
        return gc_percentage
    
    @staticmethod
    def reverse_complement(seq):
        """
        Generate the reverse complement of a DNA sequence.
        
        Creates the reverse complement by reversing the sequence and
        replacing each nucleotide with its complement (A↔T, G↔C).
        Supports IUPAC ambiguous nucleotide codes.
        
        Args:
            seq: DNA sequence to reverse complement
            
        Returns:
            Reverse complement sequence
            
        Raises:
            ValueError: If sequence contains invalid characters
            
        Example:
            >>> PrimerProcessor.reverse_complement("ATCG")
            'CGAT'
        """
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to reverse_complement")
            return ""
        
        if not seq.strip():
            logger.debug("Empty sequence provided to reverse_complement")
            return ""
        
        seq_upper = seq.upper()
        
        # Define the complement mapping including IUPAC codes
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
            
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(f"Reverse complement: {seq[:20]}... -> {rev_comp[:20]}...")
                
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
        
        Checks the C/G ratio in a sequence and returns the reverse
        complement if there are more Gs than Cs. This is useful for
        probe design where C-rich sequences are preferred.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            Tuple of (possibly_reversed_sequence, was_reversed)
            
        Example:
            >>> seq, reversed = PrimerProcessor.ensure_more_c_than_g("GGGATC")
            >>> print(f"Sequence: {seq}, Was reversed: {reversed}")
        """
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to ensure_more_c_than_g")
            return seq, False
        
        if not seq.strip():
            logger.debug("Empty sequence provided to ensure_more_c_than_g")
            return seq, False
        
        seq_upper = seq.upper()
        c_count = seq_upper.count('C')
        g_count = seq_upper.count('G')
        
        logger.debug(f"C/G analysis: C={c_count}, G={g_count}")
        
        if c_count >= g_count:
            logger.debug("C count >= G count, no reversal needed")
            return seq, False
        
        # Need to reverse complement
        try:
            rev_comp = PrimerProcessor.reverse_complement(seq)
            logger.debug("Sequence reversed to ensure more Cs than Gs")
            return rev_comp, True
        except ValueError as e:
            error_msg = f"Failed to reverse complement sequence: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            # Return original sequence if reverse complement fails
            return seq, False
        except Exception as e:
            error_msg = f"Unexpected error in ensure_more_c_than_g: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return seq, False