#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to test primer filtering on real data from Direct.xlsx.
This script uses your actual SequenceUtils implementation and assigns penalties
since they are not present in the input file.
"""

import pandas as pd
import os
import sys
import logging
import random
from pprint import pprint
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

# Now try importing with the exact module name (check case sensitivity)
from ddprimer.config import Config
from ddprimer.utils.sequence_utils import SequenceUtils
from ddprimer.core.primer_processor import PrimerProcessor

# Configure logging
logging.basicConfig(level=logging.DEBUG, 
                   format='%(levelname)s - %(message)s')
logger = logging.getLogger()

# Add the parent directory to sys.path to import your modules
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

def load_primer_data(excel_file):
    """Load primer data from the Excel file and assign penalties."""
    try:
        df = pd.read_excel(excel_file)
        logger.info(f"Loaded {len(df)} primers from {excel_file}")
        
        # Rename columns to match expected format
        column_mapping = {
            'sequence name': 'Gene',
            'sequence': 'Amplicon',
            'F': 'Primer F',
            'R': 'Primer R'
        }
        
        # Rename only the columns that exist
        existing_columns = set(df.columns).intersection(column_mapping.keys())
        rename_dict = {col: column_mapping[col] for col in existing_columns}
        df = df.rename(columns=rename_dict)
        
        # Add the missing penalty columns
        if 'Penalty F' not in df.columns:
            # Generate realistic penalties (typically between 0 and 4, with most being under 1)
            df['Penalty F'] = [round(random.random() * 1.5, 2) for _ in range(len(df))]
            logger.info("Added 'Penalty F' column with generated values")
        
        if 'Penalty R' not in df.columns:
            df['Penalty R'] = [round(random.random() * 1.5, 2) for _ in range(len(df))]
            logger.info("Added 'Penalty R' column with generated values")
        
        if 'Pair Penalty' not in df.columns:
            # Pair penalty is typically the sum of forward and reverse penalties
            df['Pair Penalty'] = df['Penalty F'] + df['Penalty R']
            logger.info("Added 'Pair Penalty' column with generated values")
        
        # Add Probe column if missing (set to N/A)
        if 'Probe' not in df.columns:
            df['Probe'] = 'N/A'
            logger.info("Added 'Probe' column with N/A values")
        
        # Log the updated DataFrame structure
        logger.info(f"Updated DataFrame columns: {', '.join(df.columns)}")
        
        return df
    except Exception as e:
        logger.error(f"Error loading Excel file: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error

def test_filter_by_penalty(df):
    """Test the filter_by_penalty method."""
    logger.info("\n" + "="*40)
    logger.info("TESTING PENALTY FILTER")
    logger.info("="*40)
    
    # Display the penalty threshold
    logger.info(f"Using PENALTY_MAX: {Config.PENALTY_MAX}")
    
    # Log all primer penalties before filtering
    logger.info("PENALTIES BEFORE FILTERING:")
    for _, row in df.iterrows():
        gene = row.get("Gene", "Unknown")
        p_f = row.get("Penalty F", "N/A")
        p_r = row.get("Penalty R", "N/A")
        p_pair = row.get("Pair Penalty", "N/A")
        
        logger.info(f"Gene: {gene} - F: {p_f}, R: {p_r}, Pair: {p_pair}")
    
    # Apply the filter
    filtered_df = PrimerProcessor.filter_by_penalty(df)
    
    # Log the primers that passed the filter
    logger.info("\nPENALTIES AFTER FILTERING:")
    for _, row in filtered_df.iterrows():
        gene = row.get("Gene", "Unknown")
        p_f = row.get("Penalty F", "N/A")
        p_r = row.get("Penalty R", "N/A")
        p_pair = row.get("Pair Penalty", "N/A")
        
        logger.info(f"Gene: {gene} - F: {p_f}, R: {p_r}, Pair: {p_pair}")
    
    # Find primers that were filtered out
    filtered_out = set(df["Gene"]) - set(filtered_df["Gene"])
    if filtered_out:
        logger.info("\nPRIMERS FILTERED OUT BY PENALTY:")
        for gene in filtered_out:
            row = df[df["Gene"] == gene].iloc[0]
            p_f = row.get("Penalty F", "N/A")
            p_r = row.get("Penalty R", "N/A")
            p_pair = row.get("Pair Penalty", "N/A")
            logger.info(f"Gene: {gene} - F: {p_f}, R: {p_r}, Pair: {p_pair}")
    else:
        logger.info("\nNo primers were filtered out by penalty")
    
    logger.info(f"\nPenalty filter: {len(df)} -> {len(filtered_df)} primers")
    return filtered_df

def test_filter_by_repeats(df):
    """Test the filter_by_repeats method."""
    logger.info("\n" + "="*40)
    logger.info("TESTING REPEAT FILTER")
    logger.info("="*40)
    
    # Log sequences before filtering
    logger.info("SEQUENCES BEFORE FILTERING:")
    for _, row in df.iterrows():
        gene = row.get("Gene", "Unknown")
        primer_f = row.get("Primer F", "N/A")
        primer_r = row.get("Primer R", "N/A")
        probe = row.get("Probe", "N/A")
        
        # Check if this sequence has disallowed repeats using the actual SequenceUtils
        has_f_repeats = "REPEAT!" if SequenceUtils.has_disallowed_repeats(primer_f) else "OK"
        has_r_repeats = "REPEAT!" if SequenceUtils.has_disallowed_repeats(primer_r) else "OK"
        has_p_repeats = "REPEAT!" if probe != "N/A" and SequenceUtils.has_disallowed_repeats(probe) else "OK"
        
        logger.info(f"Gene: {gene}")
        logger.info(f"  F: {primer_f} [{has_f_repeats}]")
        logger.info(f"  R: {primer_r} [{has_r_repeats}]")
        if probe != "N/A":
            logger.info(f"  P: {probe} [{has_p_repeats}]")
    
    # Apply the filter
    filtered_df = PrimerProcessor.filter_by_repeats(df)
    
    # Find primers that were filtered out
    filtered_out = set(df["Gene"]) - set(filtered_df["Gene"])
    if filtered_out:
        logger.info("\nPRIMERS FILTERED OUT BY REPEATS:")
        for gene in filtered_out:
            row = df[df["Gene"] == gene].iloc[0]
            primer_f = row.get("Primer F", "N/A")
            primer_r = row.get("Primer R", "N/A")
            probe = row.get("Probe", "N/A")
            
            has_f_repeats = SequenceUtils.has_disallowed_repeats(primer_f)
            has_r_repeats = SequenceUtils.has_disallowed_repeats(primer_r)
            has_p_repeats = probe != "N/A" and SequenceUtils.has_disallowed_repeats(probe)
            
            logger.info(f"Gene: {gene}")
            if has_f_repeats:
                logger.info(f"  F: {primer_f} [FILTERED - has repeats]")
            if has_r_repeats:
                logger.info(f"  R: {primer_r} [FILTERED - has repeats]")
            if has_p_repeats:
                logger.info(f"  P: {probe} [FILTERED - has repeats]")
    else:
        logger.info("\nNo primers were filtered out by repeats")
    
    logger.info(f"\nRepeat filter: {len(df)} -> {len(filtered_df)} primers")
    return filtered_df

def test_filter_by_gc_content(df):
    """Test the filter_by_gc_content method."""
    logger.info("\n" + "="*40)
    logger.info("TESTING GC CONTENT FILTER")
    logger.info("="*40)
    
    # Display GC content limits
    logger.info(f"Using GC content limits: Min={Config.SEQUENCE_MIN_GC}%, Max={Config.SEQUENCE_MAX_GC}%")
    
    # Calculate and log GC content before filtering
    logger.info("GC CONTENT BEFORE FILTERING:")
    for _, row in df.iterrows():
        gene = row.get("Gene", "Unknown")
        amplicon = row.get("Amplicon", "N/A")
        
        if amplicon != "N/A":
            gc_content = SequenceUtils.calculate_gc(amplicon)
            status = "OK" if Config.SEQUENCE_MIN_GC <= gc_content <= Config.SEQUENCE_MAX_GC else "OUT OF RANGE"
            logger.info(f"Gene: {gene} - GC: {gc_content:.1f}% [{status}]")
            # Print part of the amplicon for validation
            amplicon_summary = f"{amplicon[:20]}...{amplicon[-20:]}" if len(amplicon) > 40 else amplicon
            logger.info(f"  Amplicon: {amplicon_summary}")
        else:
            logger.info(f"Gene: {gene} - No amplicon")
    
    # Apply the filter
    filtered_df = PrimerProcessor.filter_by_gc_content(df)
    
    # Find primers that were filtered out
    filtered_out = set(df["Gene"]) - set(filtered_df["Gene"])
    if filtered_out:
        logger.info("\nPRIMERS FILTERED OUT BY GC CONTENT:")
        for gene in filtered_out:
            row = df[df["Gene"] == gene].iloc[0]
            amplicon = row.get("Amplicon", "N/A")
            
            if amplicon != "N/A":
                gc_content = SequenceUtils.calculate_gc(amplicon)
                logger.info(f"Gene: {gene} - GC: {gc_content:.1f}% [FILTERED - outside {Config.SEQUENCE_MIN_GC}-{Config.SEQUENCE_MAX_GC}%]")
            else:
                logger.info(f"Gene: {gene} - No amplicon [FILTERED]")
    else:
        logger.info("\nNo primers were filtered out by GC content")
    
    logger.info(f"\nGC content filter: {len(df)} -> {len(filtered_df)} primers")
    return filtered_df

def deliberately_fail_some_penalties(df):
    """
    Deliberately adjust some penalties to exceed the threshold
    to demonstrate the penalty filter working.
    """
    logger.info("\n" + "="*40)
    logger.info("ADJUSTING SOME PENALTIES FOR TESTING")
    logger.info("="*40)
    
    if len(df) >= 2:
        # Set high penalties for ~20% of the primers
        num_to_fail = max(1, int(len(df) * 0.2))
        
        # Get random indices to fail
        indices_to_fail = random.sample(range(len(df)), num_to_fail)
        
        for idx in indices_to_fail:
            # Set penalties above the threshold
            high_penalty = Config.PENALTY_MAX + random.uniform(0.5, 2.0)
            df.loc[idx, 'Penalty F'] = round(high_penalty, 2)
            df.loc[idx, 'Penalty R'] = round(high_penalty * 0.8, 2)
            df.loc[idx, 'Pair Penalty'] = round(high_penalty * 1.8, 2)
            
            gene = df.loc[idx, 'Gene']
            logger.info(f"Adjusted penalties for {gene} to exceed threshold:")
            logger.info(f"  F: {df.loc[idx, 'Penalty F']}")
            logger.info(f"  R: {df.loc[idx, 'Penalty R']}")
            logger.info(f"  Pair: {df.loc[idx, 'Pair Penalty']}")
    
    return df

def main():
    """Main function to test the primer filters on real data."""
    # Path to your Excel file - using the file provided instead of hardcoding the path
    excel_file = "/Users/jakob/ddPrimer/test_files/Primers/Direct.xlsx"
    
    logger.info("Starting primer filter tests on real data...")
    logger.info(f"Testing file: {excel_file}")
    
    # Load the data
    df = load_primer_data(excel_file)
    if df.empty:
        logger.error("No data to process. Exiting.")
        return
    
    # Log some basic info about the data
    logger.info(f"Columns in the dataset: {', '.join(df.columns)}")
    logger.info(f"Number of primers: {len(df)}")
    
    # Make a copy for sequential filtering tests
    working_df = df.copy()
    
    # Deliberately adjust some penalties to demonstrate filtering
    working_df = deliberately_fail_some_penalties(working_df)
    
    # Test each filter one by one
    after_penalty = test_filter_by_penalty(working_df)
    after_repeats = test_filter_by_repeats(after_penalty)
    after_gc = test_filter_by_gc_content(after_repeats)
    
    # Final results
    logger.info("\n" + "="*40)
    logger.info("FINAL RESULTS")
    logger.info("="*40)
    logger.info(f"Starting primers: {len(df)}")
    logger.info(f"After penalty filter: {len(after_penalty)}")
    logger.info(f"After repeat filter: {len(after_repeats)}")
    logger.info(f"After GC content filter: {len(after_gc)}")
    
    # List of all genes that passed all filters
    logger.info("\nGenes that passed all filters:")
    for gene in sorted(after_gc["Gene"].unique()):
        logger.info(f"- {gene}")
    
    # List of all genes that failed any filter
    failed_genes = set(df["Gene"]) - set(after_gc["Gene"])
    if failed_genes:
        logger.info("\nGenes that failed filtering:")
        for gene in sorted(failed_genes):
            # Determine which filter caused the failure
            in_penalty = gene in set(after_penalty["Gene"])
            in_repeats = gene in set(after_repeats["Gene"])
            
            if not in_penalty:
                filter_failed = "penalty"
            elif not in_repeats:
                filter_failed = "repeats"
            else:
                filter_failed = "GC content"
            
            logger.info(f"- {gene} (failed {filter_failed} filter)")

if __name__ == "__main__":
    main()