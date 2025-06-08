#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer processing module for ddPrimer pipeline.

Contains functionality for:
1. Primer filtering based on various criteria
2. Primer penalty calculation and thresholding
3. GC content analysis
4. Repeat sequence detection
5. BLAST specificity filtering
6. Internal oligo processing
"""

import pandas as pd
import logging

# Import package modules
from ..config import Config
from ..utils import SequenceUtils


class PrimerProcessor:
    """Handles primer processing and filtering operations."""
    
    @staticmethod
    def filter_by_penalty(df, max_penalty=None):
        """
        Filter primers by penalty scores.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer data
            max_penalty (float): Maximum allowed penalty (default: from Config)
            
        Returns:
            pandas.DataFrame: Filtered DataFrame
        """
        logger = logging.getLogger()
        
        if max_penalty is None:
            max_penalty = Config.PENALTY_MAX
        
        logger.debug(f"=== PENALTY FILTER DEBUG ===")
        logger.debug(f"Penalty threshold: {max_penalty}")
        
        # Track passed and failed sequences
        passed_genes = []
        failed_genes = []
        
        # Handle empty DataFrame
        if df.empty:
            logging.warning("Empty DataFrame - no primers to filter")
            return df
        
        # Check if required column exists
        if "Pair Penalty" not in df.columns:
            logging.warning("'Pair Penalty' column not found in filter_by_penalty method")
            logging.debug(f"Available columns: {df.columns.tolist()}")
            
            # Try to create it from individual penalties
            if "Penalty F" in df.columns and "Penalty R" in df.columns:
                logging.debug("Creating 'Pair Penalty' from individual penalties")
                df["Pair Penalty"] = df["Penalty F"] + df["Penalty R"]
            else:
                logging.warning("Cannot create 'Pair Penalty' - missing individual penalties")
                # Add a default high value
                df["Pair Penalty"] = max_penalty * 2
        
        # Debug log all primers before filtering
        logger.debug("PRIMERS BEFORE PENALTY FILTERING:")
        for index, row in df.iterrows():
            gene = row.get("Gene", f"Unknown_{index}")
            p_f = row.get("Penalty F", "N/A")
            p_r = row.get("Penalty R", "N/A")
            p_pair = row.get("Pair Penalty", "N/A")
            
            try:
                will_pass = float(p_pair) <= max_penalty
                status = "PASS" if will_pass else "FAIL"
                
                if will_pass:
                    passed_genes.append(gene)
                else:
                    failed_genes.append(gene)
                    # Get primer sequences for deeper debugging
                    primer_f = row.get("Primer F", "N/A")
                    primer_r = row.get("Primer R", "N/A")
                    logger.debug(f"Gene: {gene} will be FILTERED OUT due to high penalty")
                    logger.debug(f"  Penalties - F: {p_f}, R: {p_r}, Pair: {p_pair} > threshold {max_penalty}")
                    logger.debug(f"  Forward primer: {primer_f}")
                    logger.debug(f"  Reverse primer: {primer_r}")
                
                logger.debug(f"Gene: {gene} - F: {p_f}, R: {p_r}, Pair: {p_pair} [{status}]")
            except (ValueError, TypeError):
                logger.debug(f"Gene: {gene} - F: {p_f}, R: {p_r}, Pair: {p_pair} [ERROR - invalid penalty values]")
                failed_genes.append(gene)
        
        # Store original DataFrame for comparison
        df_before = df.copy()
        
        # Now filter
        df = df[df["Pair Penalty"] <= max_penalty].reset_index(drop=True)
        
        # Log filtering results
        filtered_out = set(df_before["Gene"]) - set(df["Gene"])
        logger.debug(f"\nPENALTY FILTER RESULTS:")
        logger.debug(f"Total primers before: {len(df_before)}")
        logger.debug(f"Total primers after: {len(df)}")
        logger.debug(f"Primers filtered out: {len(filtered_out)}")
        
        # Log summary of passed and failed genes
        logger.debug(f"PASSED ({len(set(df['Gene']))}): {', '.join(set(df['Gene']))}")
        logger.debug(f"FAILED ({len(filtered_out)}): {', '.join(filtered_out)}")
        
        logger.debug(f"=== END PENALTY FILTER DEBUG ===")
        
        return df
    
    @staticmethod
    def filter_by_repeats(primers):
        """
        Filter primers by disallowed repeats (GGGG, CCCC) - only in PRIMERS, NOT in PROBES.
        
        Args:
            primers (list): List of primer dictionaries or DataFrame
            
        Returns:
            list: Filtered primer dictionaries
        """
        logger = logging.getLogger()
        
        logger.debug(f"=== REPEAT FILTER DEBUG ===")
        
        # Track passed and failed sequences
        passed_genes = []
        failed_genes = []
        
        # Convert to DataFrame if not already
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # Debug log all primers before filtering
        logger.debug("PRIMERS BEFORE REPEAT FILTERING:")
        for index, row in df.iterrows():
            gene = row.get("Gene", f"Unknown_{index}")
            primer_f = row.get("Primer F", "")
            primer_r = row.get("Primer R", "")
            probe = row.get("Probe", "") if "Probe" in df.columns else ""
            
            # Check only the primers for repeats - NOT the probe
            has_f_repeats = SequenceUtils.has_disallowed_repeats(primer_f)
            has_r_repeats = SequenceUtils.has_disallowed_repeats(primer_r)
            
            will_fail = has_f_repeats or has_r_repeats
            
            if will_fail:
                failed_genes.append(gene)
                logger.debug(f"Gene: {gene} will be FILTERED OUT due to repeats")
                
                # Identify the specific repeats in each sequence
                if has_f_repeats:
                    logger.debug(f"  Forward primer: {primer_f}")
                    if "GGGG" in primer_f:
                        logger.debug(f"    Contains GGGG at position {primer_f.find('GGGG')}")
                    if "CCCC" in primer_f:
                        logger.debug(f"    Contains CCCC at position {primer_f.find('CCCC')}")
                
                if has_r_repeats:
                    logger.debug(f"  Reverse primer: {primer_r}")
                    if "GGGG" in primer_r:
                        logger.debug(f"    Contains GGGG at position {primer_r.find('GGGG')}")
                    if "CCCC" in primer_r:
                        logger.debug(f"    Contains CCCC at position {primer_r.find('CCCC')}")
                
                # Just log probes with repeats for information, but don't filter based on them
                if probe and SequenceUtils.has_disallowed_repeats(probe):
                    logger.debug(f"  NOTE: Probe has repeats but will NOT be filtered out for this reason:")
                    logger.debug(f"  Probe: {probe}")
                    if "GGGG" in probe:
                        logger.debug(f"    Contains GGGG at position {probe.find('GGGG')}")
                    if "CCCC" in probe:
                        logger.debug(f"    Contains CCCC at position {probe.find('CCCC')}")
            else:
                passed_genes.append(gene)
                logger.debug(f"Gene: {gene} - No repeats found in primers, will PASS")
        
        # Store original DataFrame for comparison
        df_before = df.copy()
        
        # Check forward and reverse primers for repeats (but NOT the probe)
        df["Has_Repeats_F"] = df["Primer F"].apply(SequenceUtils.has_disallowed_repeats)
        df["Has_Repeats_R"] = df["Primer R"].apply(SequenceUtils.has_disallowed_repeats)
        
        # Filter by primer repeats only
        df = df[~(df["Has_Repeats_F"] | df["Has_Repeats_R"])].reset_index(drop=True)
        
        # Remove temporary columns
        df = df.drop(columns=[col for col in df.columns if col.startswith("Has_Repeats_")])
        
        # Log filtering results
        filtered_out = set(df_before["Gene"]) - set(df["Gene"])
        logger.debug(f"\nREPEAT FILTER RESULTS:")
        logger.debug(f"Total primers before: {len(df_before)}")
        logger.debug(f"Total primers after: {len(df)}")
        logger.debug(f"Primers filtered out: {len(filtered_out)}")
        
        # Log summary of passed and failed genes
        logger.debug(f"PASSED ({len(set(df['Gene']))}): {', '.join(set(df['Gene']))}")
        logger.debug(f"FAILED ({len(filtered_out)}): {', '.join(filtered_out)}")
        
        logger.debug(f"=== END REPEAT FILTER DEBUG ===")
        
        # Return as list of dictionaries or DataFrame
        return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
    
    @staticmethod
    def filter_by_gc_content(primers, min_gc=None, max_gc=None):
        """
        Filter primers by amplicon GC content.
        
        Args:
            primers (list): List of primer dictionaries or DataFrame
            min_gc (float): Minimum GC percentage (default: from Config)
            max_gc (float): Maximum GC percentage (default: from Config)
            
        Returns:
            list: Filtered primer dictionaries
        """
        import logging
        logger = logging.getLogger()
        
        logger.debug(f"=== GC CONTENT FILTER DEBUG ===")
        
        if min_gc is None:
            min_gc = Config.SEQUENCE_MIN_GC
        if max_gc is None:
            max_gc = Config.SEQUENCE_MAX_GC
        
        logger.debug(f"GC content limits: Min={min_gc}%, Max={max_gc}%")
        
        # Track passed and failed sequences
        passed_genes = []
        failed_genes = []
        
        # Convert to DataFrame if not already
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # First check if we're missing amplicons and log this issue
        missing_amplicons = df["Amplicon"].isna() | (df["Amplicon"] == "")
        if missing_amplicons.any():
            logger.debug(f"WARNING: {missing_amplicons.sum()} rows have missing amplicons")
            for _, row in df[missing_amplicons].iterrows():
                gene = row.get("Gene", "Unknown")
                logger.debug(f"Row with missing amplicon: Gene={gene}, Primer F={row.get('Primer F', 'N/A')}, Primer R={row.get('Primer R', 'N/A')}")
        
        # Debug log all amplicons before filtering
        logger.debug("AMPLICONS BEFORE GC CONTENT FILTERING:")
        for index, row in df.iterrows():
            gene = row.get("Gene", f"Unknown_{index}")
            amplicon = row.get("Amplicon", "")
            
            if not amplicon:
                logger.debug(f"Gene: {gene} - No amplicon, will be FILTERED OUT")
                logger.debug(f"  Forward primer: {row.get('Primer F', 'N/A')}")
                logger.debug(f"  Reverse primer: {row.get('Primer R', 'N/A')}")
                failed_genes.append(gene)
                continue
                
            # Calculate GC content
            gc_content = SequenceUtils.calculate_gc(amplicon)
            in_range = min_gc <= gc_content <= max_gc
            
            # Count G and C bases for verification
            g_count = amplicon.upper().count('G')
            c_count = amplicon.upper().count('C')
            total_gc = g_count + c_count
            total_len = len(amplicon)
            manual_gc_pct = (total_gc / total_len) * 100 if total_len > 0 else 0
            
            # Check for calculation discrepancies
            if abs(manual_gc_pct - gc_content) > 0.1:  # Allow for minor float differences
                logger.debug(f"WARNING for {gene}: GC% calculation discrepancy: Manual {manual_gc_pct:.1f}% vs SequenceUtils {gc_content:.1f}%")
            
            if in_range:
                passed_genes.append(gene)
                logger.debug(f"Gene: {gene} - GC: {gc_content:.1f}% [PASS: within {min_gc}-{max_gc}%]")
                logger.debug(f"  G: {g_count}, C: {c_count}, Total GC: {total_gc}/{total_len} ({manual_gc_pct:.1f}%)")
                # Print part of the amplicon
                if len(amplicon) > 40:
                    logger.debug(f"  Amplicon (excerpt): {amplicon[:20]}...{amplicon[-20:]}")
                else:
                    logger.debug(f"  Amplicon: {amplicon}")
            else:
                failed_genes.append(gene)
                logger.debug(f"Gene: {gene} - GC: {gc_content:.1f}% [FAIL: outside {min_gc}-{max_gc}%]")
                logger.debug(f"  G: {g_count}, C: {c_count}, Total GC: {total_gc}/{total_len} ({manual_gc_pct:.1f}%)")
                
                # Print full amplicon for detailed inspection
                logger.debug(f"  Full amplicon sequence [{len(amplicon)} bp]:")
                logger.debug(f"  {amplicon}")
                
                # Explain why it failed
                if gc_content < min_gc:
                    logger.debug(f"  REASON: GC content too low ({gc_content:.1f}% < {min_gc}%)")
                else:
                    logger.debug(f"  REASON: GC content too high ({gc_content:.1f}% > {max_gc}%)")
        
        # Store original DataFrame for comparison
        df_before = df.copy()
        
        # Filter out rows with missing amplicons first
        df = df.dropna(subset=["Amplicon"]).copy()
        df = df[df["Amplicon"] != ""].reset_index(drop=True)
        
        # Calculate GC content 
        df["Amplicon GC%"] = df["Amplicon"].apply(SequenceUtils.calculate_gc)
        
        # Filter by GC content
        df = df[(df["Amplicon GC%"] >= min_gc) & 
                (df["Amplicon GC%"] <= max_gc)].reset_index(drop=True)
        
        # Log filtering results
        filtered_out = set(df_before["Gene"]) - set(df["Gene"])
        logger.debug(f"\nGC CONTENT FILTER RESULTS:")
        logger.debug(f"Total primers before: {len(df_before)}")
        logger.debug(f"Total primers after: {len(df)}")
        logger.debug(f"Primers filtered out: {len(filtered_out)}")
        
        # Log summary of passed and failed genes
        logger.debug(f"PASSED ({len(set(df['Gene']))}): {', '.join(set(df['Gene']))}")
        logger.debug(f"FAILED ({len(filtered_out)}): {', '.join(filtered_out)}")
        
        logger.debug(f"=== END GC CONTENT FILTER DEBUG ===")
        
        # Return as list of dictionaries or DataFrame
        return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
    
    @staticmethod
    def process_internal_oligos(primers):
        """
        Reverse complement internal oligos that have more G than C.
        
        Args:
            primers (list): List of primer dictionaries or DataFrame
            
        Returns:
            list: Processed primer dictionaries
        """
        # Convert to DataFrame if not already
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # Skip if no internal oligos/probes
        if "Probe" not in df.columns:
            return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
        
        # Apply reverse complementation if needed
        for i, row in df.iterrows():
            if pd.isnull(row["Probe"]) or not row["Probe"]:
                continue
                
            probe_seq = row["Probe"]
            reversed_probe, was_reversed = SequenceUtils.ensure_more_c_than_g(probe_seq)
            
            df.at[i, "Probe"] = reversed_probe
            df.at[i, "Probe Reversed"] = was_reversed
        
        # Return as list of dictionaries or DataFrame
        return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
    
    @staticmethod
    def filter_by_blast(primers, blast_filter_factor=None):
        """
        Filter primers by BLAST specificity.
        
        Args:
            primers (list): List of primer dictionaries or DataFrame
            blast_filter_factor (float): BLAST filter factor (default: from Config)
            
        Returns:
            list: Filtered primer dictionaries
        """
        logger = logging.getLogger()
        
        logger.debug(f"=== BLAST FILTER DEBUG ===")
        
        if blast_filter_factor is None:
            blast_filter_factor = Config.BLAST_FILTER_FACTOR
        
        logger.debug(f"BLAST filter factor: {blast_filter_factor}")
        
        # Convert to DataFrame if not already
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # Skip if no BLAST results
        required_columns = ["Primer F BLAST1", "Primer F BLAST2", "Primer R BLAST1", "Primer R BLAST2"]
        if not all(col in df.columns for col in required_columns):
            logger.debug(f"Missing required BLAST columns, skipping filter")
            return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
        
        # Track passed and failed sequences
        passed_genes = []
        failed_genes = []
        
        # Debug log all primers before filtering
        logger.debug("PRIMERS BEFORE BLAST FILTERING:")
        
        # Filter by BLAST specificity
        keep_indices = []
        
        for i, row in df.iterrows():
            gene = row.get("Gene", f"Unknown_{i}")
            
            # Check forward primer
            keep_f = PrimerProcessor._passes_blast_filter(row, "Primer F", blast_filter_factor)
            
            # Check reverse primer
            keep_r = PrimerProcessor._passes_blast_filter(row, "Primer R", blast_filter_factor)
            
            # Check probe if present
            if "Probe" in df.columns and "Probe BLAST1" in df.columns:
                keep_p = pd.isna(row["Probe"]) or row["Probe"] == "" or \
                         PrimerProcessor._passes_blast_filter(row, "Probe", blast_filter_factor)
            else:
                keep_p = True
            
            # Log BLAST scores
            f_blast1 = row.get("Primer F BLAST1", "N/A")
            f_blast2 = row.get("Primer F BLAST2", "N/A")
            r_blast1 = row.get("Primer R BLAST1", "N/A")
            r_blast2 = row.get("Primer R BLAST2", "N/A")
            
            # Debug log BLAST results
            logger.debug(f"Gene: {gene}")
            logger.debug(f"  Forward BLAST: Best={f_blast1}, Second={f_blast2}, Pass={keep_f}")
            logger.debug(f"  Reverse BLAST: Best={r_blast1}, Second={r_blast2}, Pass={keep_r}")
            
            if "Probe" in df.columns and "Probe BLAST1" in df.columns and not pd.isna(row["Probe"]):
                p_blast1 = row.get("Probe BLAST1", "N/A")
                p_blast2 = row.get("Probe BLAST2", "N/A")
                logger.debug(f"  Probe BLAST: Best={p_blast1}, Second={p_blast2}, Pass={keep_p}")
            
            # Keep only if all pass
            will_pass = keep_f and keep_r and keep_p
            if will_pass:
                keep_indices.append(i)
                passed_genes.append(gene)
                logger.debug(f"  OVERALL: PASS - will be kept")
            else:
                failed_genes.append(gene)
                
                # Log detailed reason for failure
                logger.debug(f"  OVERALL: FAIL - will be filtered out")
                if not keep_f:
                    logger.debug(f"    REASON: Forward primer BLAST specificity (best={f_blast1}, second={f_blast2})")
                    if pd.notna(f_blast1) and pd.notna(f_blast2):
                        logger.debug(f"    Ratio check: {f_blast1} * {blast_filter_factor} <= {f_blast2} is {f_blast1 * blast_filter_factor <= f_blast2}")
                
                if not keep_r:
                    logger.debug(f"    REASON: Reverse primer BLAST specificity (best={r_blast1}, second={r_blast2})")
                    if pd.notna(r_blast1) and pd.notna(r_blast2):
                        logger.debug(f"    Ratio check: {r_blast1} * {blast_filter_factor} <= {r_blast2} is {r_blast1 * blast_filter_factor <= r_blast2}")
                
                if not keep_p and "Probe" in df.columns and "Probe BLAST1" in df.columns and not pd.isna(row["Probe"]):
                    logger.debug(f"    REASON: Probe BLAST specificity (best={p_blast1}, second={p_blast2})")
                    if pd.notna(p_blast1) and pd.notna(p_blast2):
                        logger.debug(f"    Ratio check: {p_blast1} * {blast_filter_factor} <= {p_blast2} is {p_blast1 * blast_filter_factor <= p_blast2}")
        
        # Store original DataFrame for comparison
        df_before = df.copy()
        
        # Apply filter
        df = df.loc[keep_indices].reset_index(drop=True)
        
        # Log filtering results
        filtered_out = set(df_before["Gene"]) - set(df["Gene"])
        logger.debug(f"\nBLAST FILTER RESULTS:")
        logger.debug(f"Total primers before: {len(df_before)}")
        logger.debug(f"Total primers after: {len(df)}")
        logger.debug(f"Primers filtered out: {len(filtered_out)}")
        
        # Log summary of passed and failed genes
        logger.debug(f"PASSED ({len(set(df['Gene']))}): {', '.join(set(df['Gene']))}")
        logger.debug(f"FAILED ({len(filtered_out)}): {', '.join(filtered_out)}")
        
        logger.debug(f"=== END BLAST FILTER DEBUG ===")
        
        # Return as list of dictionaries or DataFrame
        return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
    
    @staticmethod
    def _passes_blast_filter(row, col_prefix, filter_factor):
        """
        Helper method to check if a sequence passes the BLAST filter.
        
        Args:
            row (pandas.Series): DataFrame row
            col_prefix (str): Column prefix (e.g., "Primer F")
            filter_factor (float): BLAST filter factor
            
        Returns:
            bool: True if passes filter, False otherwise
        """
        best = row.get(f"{col_prefix} BLAST1")
        second = row.get(f"{col_prefix} BLAST2")
        
        # If best is None or NaN, no hits => discard
        if pd.isna(best):
            return False
            
        # If no second best, it's effectively unique
        if pd.isna(second):
            return True
            
        # best must be at least filter_factor times smaller => best * filter_factor <= second
        return best * filter_factor <= second