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
from ..config import Config
from ..utils.sequence_utils import SequenceUtils


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
        from ..config import Config
        import logging
        
        if max_penalty is None:
            max_penalty = Config.PENALTY_MAX
        
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
        
        # Now filter
        df = df[df["Pair Penalty"] <= max_penalty].reset_index(drop=True)
        return df
    
    @staticmethod
    def filter_by_repeats(primers):
        """
        Filter primers by disallowed repeats (GGGG, CCCC).
        
        Args:
            primers (list): List of primer dictionaries or DataFrame
            
        Returns:
            list: Filtered primer dictionaries
        """
        # Convert to DataFrame if not already
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # Check forward and reverse primers for repeats
        df["Has_Repeats_F"] = df["Primer F"].apply(SequenceUtils.has_disallowed_repeats)
        df["Has_Repeats_R"] = df["Primer R"].apply(SequenceUtils.has_disallowed_repeats)
        
        # Check internal oligo (probe) if present
        if "Probe" in df.columns:
            df["Has_Repeats_P"] = df["Probe"].apply(
                lambda x: SequenceUtils.has_disallowed_repeats(x) if pd.notnull(x) and x else False
            )
            df = df[~(df["Has_Repeats_F"] | df["Has_Repeats_R"] | df["Has_Repeats_P"])].reset_index(drop=True)
        else:
            df = df[~(df["Has_Repeats_F"] | df["Has_Repeats_R"])].reset_index(drop=True)
        
        # Remove temporary columns
        df = df.drop(columns=[col for col in df.columns if col.startswith("Has_Repeats_")])
        
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
        if min_gc is None:
            min_gc = Config.PRIMER_MIN_GC
        if max_gc is None:
            max_gc = Config.PRIMER_MAX_GC
        
        # Convert to DataFrame if not already
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # Calculate GC content if not already present
        if "Amplicon GC%" not in df.columns:
            df["Amplicon GC%"] = df["Amplicon"].apply(SequenceUtils.calculate_gc)
        
        # Filter by GC content
        df = df[(df["Amplicon GC%"] >= min_gc) & 
                (df["Amplicon GC%"] <= max_gc)].reset_index(drop=True)
        
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
        if blast_filter_factor is None:
            blast_filter_factor = Config.BLAST_FILTER_FACTOR
        
        # Convert to DataFrame if not already
        df = pd.DataFrame(primers) if not isinstance(primers, pd.DataFrame) else primers
        
        # Skip if no BLAST results
        required_columns = ["Primer F BLAST1", "Primer F BLAST2", "Primer R BLAST1", "Primer R BLAST2"]
        if not all(col in df.columns for col in required_columns):
            return df.to_dict('records') if not isinstance(primers, pd.DataFrame) else df
        
        # Filter by BLAST specificity
        keep_indices = []
        
        for i, row in df.iterrows():
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
            
            # Keep only if all pass
            if keep_f and keep_r and keep_p:
                keep_indices.append(i)
        
        # Apply filter
        df = df.loc[keep_indices].reset_index(drop=True)
        
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
