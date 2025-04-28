#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common utilities for all primer design modes in the ddPrimer pipeline.

This module contains functions shared between different pipeline modes:
- Standard mode
- Direct mode
- Alignment mode
"""

import os
import logging
import pandas as pd
from tqdm import tqdm

from ..config import Config
from ..utils import FileUtils
from ..core import (
    PrimerProcessor,
    NupackProcessor,
    BlastProcessor,
    SequenceProcessor,
    Primer3Processor
)

# Set up logging
logger = logging.getLogger("ddPrimer")


def run_primer_design_workflow(masked_sequences, output_dir, reference_file, mode='standard', 
                              genes=None, coordinate_map=None, gff_file=None, temp_dir=None,
                              second_fasta=None, skip_annotation_filtering=False, matching_status=None,
                              all_sequences=None, add_rows_function=None):
    """
    Unified primer design workflow for all modes.
    
    Args:
        masked_sequences: Dictionary of masked sequences
        output_dir: Output directory path
        reference_file: Path to reference file (FASTA, CSV, MAF, etc.)
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        genes: Gene annotations (optional, not used in direct mode)
        coordinate_map: Coordinate mapping for alignment mode (optional)
        gff_file: Path to GFF file (optional)
        temp_dir: Temporary directory (optional)
        second_fasta: Path to second FASTA file for alignment mode (optional)
        skip_annotation_filtering: Skip gene annotation filtering (optional)
        matching_status: Dictionary with reference matching status for sequences (optional)
        all_sequences: Dictionary with all original sequences (optional, for direct mode)
        add_rows_function: Function to add rows for sequences without primers (optional, for direct mode)
        
    Returns:
        bool: Success or failure
    """
    logger = logging.getLogger("ddPrimer")
    
    try:
        # Step 1: Cut sequences at restriction sites
        restriction_fragments = cut_restriction_sites(masked_sequences)
        
        if not restriction_fragments:
            logger.warning("No valid fragments after restriction site filtering. Exiting.")
            return False
        
        # Step 2: Filter fragments based on mode
        if mode == 'direct' or skip_annotation_filtering:
            # Process fragments for direct mode or when skipping annotation filtering
            filtered_fragments = process_direct_mode_fragments(restriction_fragments)
            if skip_annotation_filtering and mode != 'direct':
                logger.debug("Skipping gene annotation filtering as requested")
            logger.debug(f"Prepared {len(filtered_fragments)} fragments for processing")
        else:
            # Filter by gene overlap for standard/alignment modes
            if not genes:
                logger.warning("Gene annotations not provided for standard/alignment mode. Exiting.")
                return False
                    
            logger.info("\nFiltering sequences by gene overlap...")
            try:
                logger.debug(f"Using gene overlap margin: {Config.GENE_OVERLAP_MARGIN}")
                filtered_fragments = SequenceProcessor.filter_by_gene_overlap(restriction_fragments, genes)
                logger.debug("Gene overlap filtering completed successfully")
            except Exception as e:
                logger.error(f"Error in gene overlap filtering: {e}")
                logger.debug(e, exc_info=True)
                return False
            logger.info(f"Retained {len(filtered_fragments)} fragments after gene overlap filtering")
        
        # Step 3: Design primers with Primer3
        primer_results = run_primer3_design(filtered_fragments, mode)
        
        if not primer_results:
            logger.warning("No primers were designed by Primer3. Exiting.")
            return False
        
        # Step 4: Filter primers
        df = filter_primers(primer_results)
        
        if df is None or len(df) == 0:
            logger.warning("No primers passed filtering criteria. Exiting.")
            return False
        
        # Step 5: Calculate thermodynamic properties
        df = calculate_thermodynamics(df)
        
        # Step 6: Run BLAST for specificity
        df = run_blast(df)
        
        if df is None or len(df) == 0:
            logger.warning("No primers passed BLAST filtering. Exiting.")
            return False
        
        # Step 7: Process alignment coordinates if applicable
        if mode == 'alignment' and coordinate_map:
            df = map_alignment_coordinates(df, coordinate_map)
        
        # Step 8: Add rows for sequences without primers (only in direct mode)
        if mode == 'direct' and add_rows_function and all_sequences:
            logger.debug("Adding rows for sequences without primers and reference match failures")
            df = add_rows_function(df, all_sequences, matching_status)
        
        # Step 9: Save results to Excel file
        from ..utils import FileUtils
        output_path = FileUtils.save_results(
            df, 
            output_dir, 
            reference_file, 
            mode=mode,
            second_fasta=second_fasta,
            logger=logger
        )
        
        if output_path:
            logger.info(f"Results saved to: {output_path}")
            return True
        else:
            logger.error("Failed to save results.")
            return False
            
    except Exception as e:
        logger.error(f"Error in primer design workflow: {e}")
        logger.debug(e, exc_info=True)
        return False


def process_direct_mode_fragments(restriction_fragments):
    """
    Process fragments specifically for direct mode.
    
    Args:
        restriction_fragments: List of restriction fragments
        
    Returns:
        list: Simplified fragments for direct mode
    """
    filtered_fragments = []
    for fragment in restriction_fragments:
        # Extract the base sequence ID (without any fragment numbers)
        base_id = fragment["id"].split("_frag")[0]
        
        # Create a simplified fragment with ONLY what we need
        # Explicitly exclude location data for direct mode
        simplified_fragment = {
            "id": fragment["id"],
            "sequence": fragment["sequence"],
            "Gene": base_id  # Explicitly set Gene to base_id for direct mode
        }
        filtered_fragments.append(simplified_fragment)
    
    return filtered_fragments


def run_primer3_design(fragments, mode='standard'):
    """
    Run Primer3 design on the provided fragments.
    
    Args:
        fragments: List of sequence fragments
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        
    Returns:
        list: Primer design results
    """
    logger.info("\nDesigning primers with Primer3...")
    
    # Initialize Primer3Processor
    primer3_processor = Primer3Processor(Config)
    
    # Prepare input blocks for Primer3
    primer3_inputs, fragment_info = prepare_primer3_inputs(fragments, mode=mode)
    
    if not primer3_inputs:
        logger.warning("No valid fragments for primer design. Exiting.")
        return None
        
    # Run Primer3
    logger.debug(f"Running Primer3 on {len(primer3_inputs)} fragments...")
    try:
        # Use parallel processing with progress bar
        primer3_output = primer3_processor.run_primer3_batch_parallel(primer3_inputs)
        logger.debug(f"Primer3 execution completed successfully")
        
        # Parse the results
        primer_results = primer3_processor.parse_primer3_batch(primer3_output, fragment_info)
        logger.debug(f"Parsed {len(primer_results)} primer pairs from Primer3 output")
        
        return primer_results
        
    except Exception as e:
        logger.error(f"Error running Primer3: {e}")
        logger.debug(e, exc_info=True)
        return None


def map_alignment_coordinates(df, coordinate_map):
    """
    Map primer coordinates to the second genome for alignment mode.
    
    Args:
        df: DataFrame with primer results
        coordinate_map: Coordinate mapping dictionary
        
    Returns:
        pandas.DataFrame: Updated DataFrame with mapped coordinates
    """
    logger.debug("\nMapping primer coordinates to second genome...")
    
    # Add columns for second genome data
    df["Qry Chromosome"] = None
    df["Qry Location"] = None
    mapped_count = 0
    
    for idx, row in df.iterrows():
        chrom = row["Chromosome"]
        
        # Handle the case where Location might not contain a dash
        if "-" in str(row["Location"]):
            parts = row["Location"].split("-")
            start_pos = int(parts[0])
        else:
            try:
                start_pos = int(row["Location"])
            except (ValueError, TypeError):
                logger.warning(f"Invalid location format for index {idx}: {row['Location']}")
                continue
        
        if chrom in coordinate_map:
            # Try to map forward primer position
            try:
                for pos in range(start_pos, start_pos + len(row["Primer F"])):
                    if pos in coordinate_map[chrom]:
                        mapping = coordinate_map[chrom][pos]
                        df.at[idx, "Qry Chromosome"] = mapping["qry_src"]
                        df.at[idx, "Qry Location"] = mapping["qry_pos"]
                        mapped_count += 1
                        break
            except Exception as e:
                logger.warning(f"Error mapping primer at index {idx}: {e}")
                continue
    
    logger.debug(f"Successfully mapped {mapped_count}/{len(df)} primers to second genome")
    return df


def cut_restriction_sites(masked_sequences):
    """
    Cut sequences at restriction sites.
    
    Args:
        masked_sequences: Dictionary of sequences
        
    Returns:
        list: List of restriction fragments
    """
    logger.info("\nFiltering sequences by restriction sites...")
    try:
        logger.debug(f"Using restriction site pattern: {Config.RESTRICTION_SITE}")
        
        # Log sequence stats before restriction site cutting
        logger.debug(f"Sequences before restriction site cutting:")
        for seq_id, seq in masked_sequences.items():
            logger.debug(f"  {seq_id}: {len(seq)} bp")
        
        # Use the standard restriction site method for all modes
        restriction_fragments = SequenceProcessor.cut_at_restriction_sites(masked_sequences)
        
        # Log detailed information about restriction fragments
        logger.debug("Restriction fragments after cutting:")
        for fragment in restriction_fragments:
            logger.debug(f"  {fragment['id']}: {len(fragment['sequence'])} bp, chr={fragment.get('chr', 'NA')}, "
                        f"start={fragment.get('start', 'NA')}, end={fragment.get('end', 'NA')}")
        
        logger.debug("Restriction site filtering completed successfully")
        logger.info(f"Generated {len(restriction_fragments)} fragments after restriction site cutting")
        
        return restriction_fragments
    except Exception as e:
        logger.error(f"Error in restriction site filtering: {e}")
        logger.debug(e, exc_info=True)
        raise


def prepare_primer3_inputs(fragments, mode='standard'):
    """
    Prepare input blocks for Primer3 from sequence fragments.
    
    Args:
        fragments: List of sequence fragments
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        
    Returns:
        tuple: (primer3_inputs, fragment_info) - Lists of Primer3 input dicts and fragment info dict
    """
    primer3_inputs = []
    fragment_info = {}  # Dictionary to store fragment information
    
    for fragment in fragments:
        # Store fragment information differently based on mode
        if mode == 'direct':
            # For direct mode, only store minimal information - explicitly exclude location
            fragment_info[fragment["id"]] = {
                "gene": fragment.get("Gene", fragment["id"])
            }
        else:
            # For standard/alignment modes, store location information
            fragment_info[fragment["id"]] = {
                "chr": fragment.get("chr", ""),
                "start": fragment.get("start", 1),
                "end": fragment.get("end", len(fragment["sequence"])),
                "gene": fragment.get("Gene", fragment["id"].split("_")[-1])
            }
        
        # Create primer3 input
        primer3_input = {
            "SEQUENCE_ID": fragment["id"],
            "SEQUENCE_TEMPLATE": fragment["sequence"],
        }
        
        # Add target region if sequence is long
        if len(fragment["sequence"]) > 200:
            target_start = len(fragment["sequence"]) // 4
            target_len = len(fragment["sequence"]) // 2
            primer3_input["SEQUENCE_TARGET"] = [target_start, target_len]
            if Config.DEBUG_MODE:
                logger.debug(f"Added target region [{target_start}, {target_len}] for {fragment['id']}")
        
        primer3_inputs.append(primer3_input)
    
    return primer3_inputs, fragment_info


def filter_primers(primer_results):
    """
    Filter primers based on various criteria.
    
    Args:
        primer_results: List of primer results from Primer3
        
    Returns:
        pandas.DataFrame: Filtered primer DataFrame
    """
    logger.info("\nFiltering primers...")

    # Convert to DataFrame
    df = pd.DataFrame(primer_results)
    initial_count = len(df)
    logger.debug(f"Initial primer count: {initial_count}")

    # Filter by penalty
    try:
        logger.debug(f"Filtering by penalty with threshold: {Config.PENALTY_MAX}")
        df_before = df.copy() if Config.DEBUG_MODE else None
        df = PrimerProcessor.filter_by_penalty(df)
        
        # Log what was filtered in debug mode
        if Config.DEBUG_MODE and df_before is not None:
            filtered_ids = set(df_before['Gene']) - set(df['Gene'])
            if filtered_ids:
                logger.debug(f"Penalty filtering removed primers for: {', '.join(filtered_ids)}")
        
        logger.debug("Penalty filtering completed successfully")
    except Exception as e:
        logger.error(f"Error in penalty filtering: {e}")
        logger.debug(e, exc_info=True)
        raise
    logger.debug(f"After penalty filtering: {len(df)}/{initial_count} primers")
    
    # Filter by repeats
    try:
        logger.debug("Filtering primers by repeat sequences")
        df = PrimerProcessor.filter_by_repeats(df)
        logger.debug("Repeat filtering completed successfully")
    except Exception as e:
        logger.error(f"Error in repeat filtering: {e}")
        logger.debug(e, exc_info=True)
        raise
    logger.debug(f"After repeat filtering: {len(df)}/{initial_count} primers")
    
    # Filter by GC content
    try:
        logger.debug(f"Filtering by GC content: Min={Config.SEQUENCE_MIN_GC}, Max={Config.SEQUENCE_MAX_GC}")
        df = PrimerProcessor.filter_by_gc_content(df)
        logger.debug("GC content filtering completed successfully")
    except Exception as e:
        logger.error(f"Error in GC content filtering: {e}")
        logger.debug(e, exc_info=True)
        raise
    logger.info(f"After filtering: {len(df)}/{initial_count} primers")
    
    # Process internal oligos (reverse complement if needed)
    try:
        logger.debug("Processing internal oligos")
        df = PrimerProcessor.process_internal_oligos(df)
        logger.debug("Internal oligo processing completed successfully")
    except Exception as e:
        logger.error(f"Error in internal oligo processing: {e}")
        logger.debug(e, exc_info=True)
        raise
    
    if len(df) == 0:
        logger.warning("No primers passed filtering. Exiting.")
        return None
        
    return df


def calculate_thermodynamics(df):
    """
    Calculate thermodynamic properties using NUPACK.
    
    Args:
        df: DataFrame with primer information
        
    Returns:
        pandas.DataFrame: DataFrame with added thermodynamic properties
    """
    logger.info("\nCalculating thermodynamic properties with NUPACK...")
    
    # Calculate deltaG for forward primers
    logger.debug("Calculating ΔG for forward primers...")
    try:
        logger.debug(f"NUPACK settings: Temp={Config.NUPACK_TEMPERATURE}°C, Na={Config.NUPACK_SODIUM}M, Mg={Config.NUPACK_MAGNESIUM}M")
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing forward primers")
            df["Primer F dG"] = df["Primer F"].progress_apply(NupackProcessor.calc_deltaG)
        else:
            df["Primer F dG"] = df["Primer F"].apply(NupackProcessor.calc_deltaG)
        logger.debug("Forward primer deltaG calculation completed successfully")
    except Exception as e:
        logger.error(f"Error calculating forward primer deltaG: {e}")
        logger.debug(e, exc_info=True)
        raise
    
    # Calculate deltaG for reverse primers
    logger.debug("Calculating ΔG for reverse primers...")
    try:
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing reverse primers")
            df["Primer R dG"] = df["Primer R"].progress_apply(NupackProcessor.calc_deltaG)
        else:
            df["Primer R dG"] = df["Primer R"].apply(NupackProcessor.calc_deltaG)
        logger.debug("Reverse primer deltaG calculation completed successfully")
    except Exception as e:
        logger.error(f"Error calculating reverse primer deltaG: {e}")
        logger.debug(e, exc_info=True)
        raise
    
    # Calculate deltaG for probes if present
    if "Probe" in df.columns:
        logger.debug("Calculating ΔG for probes...")
        try:
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc="Processing probes")
                df["Probe dG"] = df["Probe"].progress_apply(lambda x: 
                                              NupackProcessor.calc_deltaG(x) 
                                              if pd.notnull(x) and x else None)
            else:
                df["Probe dG"] = df["Probe"].apply(lambda x: 
                                                 NupackProcessor.calc_deltaG(x) 
                                                 if pd.notnull(x) and x else None)
            logger.debug("Probe deltaG calculation completed successfully")
        except Exception as e:
            logger.error(f"Error calculating probe deltaG: {e}")
            logger.debug(e, exc_info=True)
            raise
    
    # Calculate deltaG for amplicons
    logger.debug("Calculating ΔG for amplicons...")
    try:
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing amplicons")
            df["Amplicon dG"] = df["Amplicon"].progress_apply(NupackProcessor.calc_deltaG)
        else:
            df["Amplicon dG"] = df["Amplicon"].apply(NupackProcessor.calc_deltaG)
        logger.debug("Amplicon deltaG calculation completed successfully")
    except Exception as e:
        logger.error(f"Error calculating amplicon deltaG: {e}")
        logger.debug(e, exc_info=True)
        raise
        
    return df


def run_blast(df):
    """
    Run BLAST for primer specificity checking.
    
    Args:
        df: DataFrame with primer information
        
    Returns:
        pandas.DataFrame: DataFrame with added BLAST results
    """
    logger.info("\nRunning BLAST for specificity checking...")
    
    # Run BLAST for forward primers
    try:
        logger.debug(f"BLAST settings: Word Size={Config.BLAST_WORD_SIZE}, E-value={Config.BLAST_EVALUE}")
        blast_results_f = []
        primers_f = df["Primer F"].tolist()
        if Config.SHOW_PROGRESS:
            primers_f_iter = tqdm(primers_f, total=len(primers_f), desc="BLASTing forward primers")
        else:
            primers_f_iter = primers_f
            
        for primer_f in primers_f_iter:
            blast1, blast2 = BlastProcessor.blast_short_seq(primer_f)
            blast_results_f.append((blast1, blast2))
        
        df["Primer F BLAST1"], df["Primer F BLAST2"] = zip(*blast_results_f)
        logger.debug("Forward primer BLAST completed successfully")
    except Exception as e:
        logger.error(f"Error in forward primer BLAST: {e}")
        logger.debug(e, exc_info=True)
        raise
    
    # Run BLAST for reverse primers
    try:
        blast_results_r = []
        primers_r = df["Primer R"].tolist()
        if Config.SHOW_PROGRESS:
            primers_r_iter = tqdm(primers_r, total=len(primers_r), desc="BLASTing reverse primers")
        else:
            primers_r_iter = primers_r
            
        for primer_r in primers_r_iter:
            blast1, blast2 = BlastProcessor.blast_short_seq(primer_r)
            blast_results_r.append((blast1, blast2))
        
        df["Primer R BLAST1"], df["Primer R BLAST2"] = zip(*blast_results_r)
        logger.debug("Reverse primer BLAST completed successfully")
    except Exception as e:
        logger.error(f"Error in reverse primer BLAST: {e}")
        logger.debug(e, exc_info=True)
        raise
    
    # Run BLAST for probes if present
    if "Probe" in df.columns:
        try:
            blast_results_p = []
            probes = df["Probe"].tolist()
            if Config.SHOW_PROGRESS:
                probes_iter = tqdm(probes, total=len(probes), desc="BLASTing probes")
            else:
                probes_iter = probes
                
            for probe in probes_iter:
                if pd.notnull(probe) and probe:
                    blast1, blast2 = BlastProcessor.blast_short_seq(probe)
                else:
                    blast1, blast2 = None, None
                blast_results_p.append((blast1, blast2))
            
            df["Probe BLAST1"], df["Probe BLAST2"] = zip(*blast_results_p)
            logger.debug("Probe BLAST completed successfully")
        except Exception as e:
            logger.error(f"Error in probe BLAST: {e}")
            logger.debug(e, exc_info=True)
            raise
    
    # Filter by BLAST specificity
    try:
        logger.debug(f"Filtering by BLAST specificity, filter factor: {Config.BLAST_FILTER_FACTOR}")
        initial_count = len(df)
        df = PrimerProcessor.filter_by_blast(df)
        logger.debug("BLAST filtering completed successfully")
    except Exception as e:
        logger.error(f"Error in BLAST filtering: {e}")
        logger.debug(e, exc_info=True)
        raise
    logger.info(f"After BLAST filtering: {len(df)}/{initial_count} primers")
    
    if len(df) == 0:
        logger.warning("No primers passed BLAST filtering. Exiting.")
        return None
        
    return df


def save_results(df, output_dir, input_file, mode='standard'):
    """
    Save results to an Excel file.
    
    Args:
        df: DataFrame with primer results
        output_dir: Output directory
        input_file: Path to the input file (FASTA, CSV, etc.)
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        
    Returns:
        str: Path to the output file
    """
    from ..utils import FileUtils
    return FileUtils.save_results(df, output_dir, input_file, mode, logger=logger)