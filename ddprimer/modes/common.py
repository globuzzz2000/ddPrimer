#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common utilities for all primer design modes in the ddPrimer pipeline.

This module contains functions shared between different pipeline modes:
- Standard mode
- Direct mode  
- Alignment mode

Contains functionality for:
1. Unified primer design workflow execution
2. Fragment filtering based on pipeline mode
3. Primer3 design and result processing
4. Thermodynamic property calculations
5. BLAST specificity checking
6. Coordinate mapping for alignment mode
7. Enhanced SNP processing with fixed SNP substitution support

This module integrates with the broader ddPrimer pipeline to provide
shared workflow components across all primer design modes.
"""

import logging
import pandas as pd
from tqdm import tqdm

# Import package modules
from ..config import Config, SequenceProcessingError, PrimerDesignError
from ..utils import FileIO
from ..core import PrimerProcessor, BlastProcessor, SequenceProcessor, Primer3Processor, ViennaRNAProcessor 

# Set up module logger
logger = logging.getLogger(__name__)


def run_primer_design_workflow(masked_sequences, output_dir, reference_file, mode='standard', 
                              genes=None, coordinate_map=None, gff_file=None, temp_dir=None,
                              second_fasta=None, skip_annotation_filtering=False, matching_status=None,
                              all_sequences=None, add_rows_function=None):
    """
    Unified primer design workflow for all modes.
    
    Executes the complete primer design pipeline including restriction site processing,
    fragment filtering, primer design, thermodynamic calculations, and BLAST analysis.
    Now supports enhanced SNP processing with fixed SNP substitution.
    
    Args:
        masked_sequences: Dictionary of processed sequences (masked/substituted)
        output_dir: Output directory path
        reference_file: Path to reference file (FASTA, CSV, MAF, etc.)
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        genes: Gene annotations (not used in direct mode)
        coordinate_map: Coordinate mapping for alignment mode
        gff_file: Path to GFF file
        temp_dir: Temporary directory
        second_fasta: Path to second FASTA file for alignment mode
        skip_annotation_filtering: Skip gene annotation filtering
        matching_status: Dictionary with reference matching status for sequences
        all_sequences: Dictionary with all original sequences (for direct mode)
        add_rows_function: Function to add rows for sequences without primers
        
    Returns:
        True if workflow completed successfully, False otherwise
        
    Raises:
        SequenceProcessingError: If sequence processing fails
        PrimerDesignError: If primer design fails
    """
    logger.debug(f"Starting primer design workflow in {mode} mode")
    logger.debug("Enhanced SNP processing: Fixed SNPs substituted, variable SNPs masked")
    
    try:
        # Step 1: Cut sequences at restriction sites
        restriction_fragments = process_restriction_sites(masked_sequences)
        logger.debug(f"CHECKPOINT 1: {len(restriction_fragments)} restriction fragments")
        
        if not restriction_fragments:
            logger.warning("No valid fragments after restriction site filtering. Exiting.")
            return False
        
        # Step 2: Filter fragments based on mode
        filtered_fragments = filter_fragments_by_mode(
            restriction_fragments, 
            mode, 
            genes, 
            skip_annotation_filtering
        )
        
        # Check if filtering returned empty list or None
        if not filtered_fragments:  # This handles both None and empty list
            logger.warning("No fragments passed filtering. Exiting.")
            if mode in ['standard', 'alignment'] and not skip_annotation_filtering:
                logger.info("Suggestion: Try using --noannotation flag to skip gene filtering")
            return False
        
        logger.debug(f"CHECKPOINT 2: {len(filtered_fragments)} filtered fragments")
        
        # Step 3: Design primers with Primer3
        primer_results = design_primers_with_primer3(filtered_fragments, mode)
        logger.debug(f"CHECKPOINT 3: {len(primer_results)} primer results from Primer3")
        
        if not primer_results:
            logger.warning("No primers were designed by Primer3. Exiting.")
            return False
        
        # Create a temporary DataFrame to check what we have before filtering
        temp_df = pd.DataFrame(primer_results)
        logger.debug(f"CHECKPOINT 3b: {len(temp_df)} primer records, {len(temp_df['Gene'].unique())} unique genes")
        
        # Log penalty distribution
        if 'Pair Penalty' in temp_df.columns:
            penalty_stats = temp_df['Pair Penalty'].describe()
            logger.debug(f"Penalty distribution: min={penalty_stats['min']:.2f}, mean={penalty_stats['mean']:.2f}, max={penalty_stats['max']:.2f}")
            below_threshold = (temp_df['Pair Penalty'] <= Config.PENALTY_MAX).sum()
            logger.debug(f"Primers with penalty <= {Config.PENALTY_MAX}: {below_threshold}/{len(temp_df)}")
        
        # Step 4: Filter primers
        df = filter_primers(primer_results)
        logger.debug(f"CHECKPOINT 4: {len(df)} primers after all filtering")
        
        if df is None or len(df) == 0:
            logger.warning("No primers passed filtering criteria. Exiting.")
            return False
        
        # Step 5: Calculate thermodynamic properties
        df = calculate_thermodynamics(df)
        
        # Step 6: Run BLAST for specificity
        df = run_blast_specificity(df)
        
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
        output_path = FileIO.save_results(
            df, 
            output_dir, 
            reference_file, 
            mode=mode,
            second_fasta=second_fasta
        )
        
        if output_path:
            logger.info(f"\nResults saved to: {output_path}")
            return True
        else:
            logger.error("Failed to save results.")
            return False
            
    except SequenceProcessingError as e:
        error_msg = f"Sequence processing error: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False
    except PrimerDesignError as e:
        error_msg = f"Primer design error: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False
    except Exception as e:
        error_msg = f"Error in primer design workflow: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False


def filter_fragments_by_mode(restriction_fragments, mode, genes, skip_annotation_filtering=False):
    """
    Filter fragments based on pipeline mode.
    
    Applies mode-specific filtering to restriction fragments, including
    gene overlap extraction for standard/alignment modes or direct
    processing for direct mode.
    
    Args:
        restriction_fragments: List of restriction fragments
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        genes: Gene annotations
        skip_annotation_filtering: Skip gene annotation filtering
        
    Returns:
        List of filtered fragments
        
    Raises:
        SequenceProcessingError: If there's an error in fragment filtering
    """
    # Process fragments based on mode and annotation settings
    if mode == 'direct' or skip_annotation_filtering:
        # Process fragments for direct mode or when skipping annotation filtering
        filtered_fragments = process_direct_mode_fragments(restriction_fragments)
        if skip_annotation_filtering and mode != 'direct':
            logger.debug("Skipping gene annotation filtering as requested")
        logger.debug(f"Prepared {len(filtered_fragments)} fragments for processing")
        return filtered_fragments
    else:
        # Extract ALL gene overlaps for standard/alignment modes
        if not genes:
            logger.error("Gene annotations not provided for standard/alignment mode.")
            logger.info("Hint: Use --noannotation flag to skip gene filtering, or provide a GFF file")
            return []  # Return empty list instead of None
                
        logger.info("Extracting gene-overlapping regions...")
        try:
            # Import AnnotationProcessor here to avoid circular imports
            from ..core.annotation_processor import AnnotationProcessor
            
            logger.debug(f"Using gene overlap margin: {Config.GENE_OVERLAP_MARGIN}")
            # Use the enhanced function that extracts ALL gene overlaps
            filtered_fragments = AnnotationProcessor.filter_by_gene_overlap_enhanced(restriction_fragments, genes)
            logger.debug("Gene overlap extraction completed successfully")
            logger.info(f"Extracted {len(filtered_fragments)} gene fragments from {len(restriction_fragments)} restriction fragments")
            return filtered_fragments
        except Exception as e:
            error_msg = f"Error in gene overlap extraction: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e


def process_direct_mode_fragments(restriction_fragments):
    """
    Process fragments specifically for direct mode.
    
    Creates simplified fragments for direct mode with minimal metadata,
    explicitly excluding location data that's not relevant for this mode.
    
    Args:
        restriction_fragments: List of restriction fragments
        
    Returns:
        List of simplified fragments for direct mode
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


def design_primers_with_primer3(fragments, mode='standard'):
    """
    Run Primer3 design on the provided fragments.
    
    Prepares input for Primer3, executes primer design in parallel,
    and parses results into a standardized format.
    
    Args:
        fragments: List of sequence fragments
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        
    Returns:
        List of primer design results
        
    Raises:
        PrimerDesignError: If there's an error in primer design
    """
    logger.info("\nDesigning primers with Primer3...")
    
    # Initialize Primer3Processor
    primer3_processor = Primer3Processor(Config)
    
    # Prepare input blocks for Primer3
    primer3_inputs, fragment_info = prepare_primer3_inputs(fragments, mode=mode)
    
    if not primer3_inputs:
        logger.warning("No valid fragments for primer design.")
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
        error_msg = f"Error running Primer3: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e


def map_alignment_coordinates(df, coordinate_map):
    """
    Map primer coordinates to the second genome for alignment mode.
    
    Translates primer coordinates from the reference genome to the
    corresponding positions in the second genome using the provided
    coordinate mapping.
    
    Args:
        df: DataFrame with primer results
        coordinate_map: Coordinate mapping dictionary
        
    Returns:
        Updated DataFrame with mapped coordinates
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
                logger.warning(f"Error mapping primer at index {idx}: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                continue
    
    logger.debug(f"Successfully mapped {mapped_count}/{len(df)} primers to second genome")
    return df


def process_restriction_sites(processed_sequences):
    """
    Cut sequences at restriction sites.
    
    Processes sequences (which may already be masked/substituted) to identify 
    and cut at restriction sites, creating fragments suitable for primer design.
    
    Args:
        processed_sequences: Dictionary of sequences (masked and/or with fixed SNPs substituted)
        
    Returns:
        List of restriction fragments
        
    Raises:
        SequenceProcessingError: If there's an error in restriction site processing
    """
    logger.info("\nFiltering sequences by restriction sites...")
    logger.debug("Processing sequences that may contain masked variants and fixed SNP substitutions")
    
    try:
        logger.debug(f"Using restriction site pattern: {Config.RESTRICTION_SITE}")
        
        # Log sequence stats before restriction site cutting
        if Config.DEBUG_MODE:
            logger.debug(f"Sequences before restriction site cutting:")
            for seq_id, seq in processed_sequences.items():
                n_count = seq.count('N')
                lowercase_count = sum(1 for c in seq if c.islower())
                logger.debug(f"  {seq_id}: {len(seq)} bp, {n_count} hard masked, {lowercase_count} soft masked")
        
        # Use the standard restriction site method for all modes
        restriction_fragments = SequenceProcessor.cut_at_restriction_sites(processed_sequences)
        
        # Log detailed information about restriction fragments
        if Config.DEBUG_MODE:
            logger.debug("Restriction fragments after cutting:")
            for fragment in restriction_fragments:
                frag_seq = fragment['sequence']
                n_count = frag_seq.count('N')
                lowercase_count = sum(1 for c in frag_seq if c.islower())
                logger.debug(f"  {fragment['id']}: {len(frag_seq)} bp, chr={fragment.get('chr', 'NA')}, "
                            f"start={fragment.get('start', 'NA')}, end={fragment.get('end', 'NA')}, "
                            f"{n_count} hard masked, {lowercase_count} soft masked")
        
        logger.debug("Restriction site filtering completed successfully")
        logger.info(f"Generated {len(restriction_fragments)} fragments after restriction site cutting")
        
        return restriction_fragments
    except Exception as e:
        error_msg = f"Error in restriction site filtering: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise SequenceProcessingError(error_msg) from e


def prepare_primer3_inputs(fragments, mode='standard'):
    """
    Prepare input blocks for Primer3 from sequence fragments.
    
    Creates Primer3-compatible input dictionaries from fragment data,
    storing fragment metadata for later result processing.
    No target regions - let Primer3 freely choose amplicon locations based on product size range.
    
    Args:
        fragments: List of sequence fragments
        mode: Pipeline mode ('standard', 'direct', or 'alignment')
        
    Returns:
        Tuple of (primer3_inputs, fragment_info) - Lists of Primer3 input dicts and fragment info dict
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

        primer3_inputs.append(primer3_input)
    
    return primer3_inputs, fragment_info


def filter_primers(primer_results):
    """
    Filter primers based on various criteria.
    
    Applies sequential filtering based on penalty scores, repeat sequences,
    and GC content to ensure high-quality primer selection.
    
    Args:
        primer_results: List of primer results from Primer3
        
    Returns:
        Filtered primer DataFrame, or None if no primers pass filtering
        
    Raises:
        PrimerDesignError: If there's an error in primer filtering
    """
    logger.debug("Filtering primers...")

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
        error_msg = f"Error in penalty filtering: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    logger.debug(f"After penalty filtering: {len(df)}/{initial_count} primers")
    
    # Filter by repeats
    try:
        logger.debug("Filtering primers by repeat sequences")
        df = PrimerProcessor.filter_by_repeats(df)
        logger.debug("Repeat filtering completed successfully")
    except Exception as e:
        error_msg = f"Error in repeat filtering: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    logger.debug(f"After repeat filtering: {len(df)}/{initial_count} primers")
    
    # Filter by GC content
    try:
        logger.debug(f"Filtering by GC content: Min={Config.SEQUENCE_MIN_GC}, Max={Config.SEQUENCE_MAX_GC}")
        df = PrimerProcessor.filter_by_gc_content(df)
        logger.debug("GC content filtering completed successfully")
    except Exception as e:
        error_msg = f"Error in GC content filtering: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    logger.info(f"After filtering: {len(df)}/{initial_count} primers")
    
    # Process internal oligos (reverse complement if needed)
    try:
        logger.debug("Processing internal oligos")
        df = PrimerProcessor.process_internal_oligos(df)
        logger.debug("Internal oligo processing completed successfully")
    except Exception as e:
        error_msg = f"Error in internal oligo processing: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    
    if len(df) == 0:
        logger.warning("No primers passed filtering.")
        return None
        
    return df


def calculate_thermodynamics(df):
    """
    Calculate thermodynamic properties using ViennaRNA.
    
    Computes minimum free energy (ΔG) for primers, probes, and amplicons
    using ViennaRNA with DNA-specific parameters.
    
    Args:
        df: DataFrame with primer information
        
    Returns:
        DataFrame with added thermodynamic properties
        
    Raises:
        PrimerDesignError: If there's an error in thermodynamic calculations
    """
        
    # Add empty columns to maintain DataFrame structure
    df["Primer F dG"] = None
    df["Primer R dG"] = None
    if "Probe" in df.columns:
        df["Probe dG"] = None
    df["Amplicon dG"] = None
        
    logger.info("\nCalculating thermodynamic properties with ViennaRNA...")
    
    # Calculate deltaG for forward primers
    logger.debug("Calculating ΔG for forward primers...")
    try:
        logger.debug(f"Thermodynamic settings: Temp={Config.THERMO_TEMPERATURE}°C, Na={Config.THERMO_SODIUM}M, Mg={Config.THERMO_MAGNESIUM}M")
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing forward primers with ViennaRNA")
            df["Primer F dG"] = df["Primer F"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        else:
            df["Primer F dG"] = df["Primer F"].apply(ViennaRNAProcessor.calc_deltaG)
        logger.debug("Forward primer deltaG calculation completed successfully")
    except Exception as e:
        error_msg = f"Error calculating forward primer deltaG: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    
    # Calculate deltaG for reverse primers
    logger.debug("Calculating ΔG for reverse primers...")
    try:
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing reverse primers with ViennaRNA")
            df["Primer R dG"] = df["Primer R"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        else:
            df["Primer R dG"] = df["Primer R"].apply(ViennaRNAProcessor.calc_deltaG)
        logger.debug("Reverse primer deltaG calculation completed successfully")
    except Exception as e:
        error_msg = f"Error calculating reverse primer deltaG: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    
    # Calculate deltaG for probes if present
    if "Probe" in df.columns:
        logger.debug("Calculating ΔG for probes...")
        try:
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc="Processing probes with ViennaRNA")
                df["Probe dG"] = df["Probe"].progress_apply(lambda x: 
                                            ViennaRNAProcessor.calc_deltaG(x) 
                                            if pd.notnull(x) and x else None)
            else:
                df["Probe dG"] = df["Probe"].apply(lambda x: 
                                               ViennaRNAProcessor.calc_deltaG(x) 
                                               if pd.notnull(x) and x else None)
            logger.debug("Probe deltaG calculation completed successfully")
        except Exception as e:
            error_msg = f"Error calculating probe deltaG: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise PrimerDesignError(error_msg) from e
    
    # Calculate deltaG for amplicons
    logger.debug("Calculating ΔG for amplicons...")
    try:
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing amplicons with ViennaRNA")
            df["Amplicon dG"] = df["Amplicon"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        else:
            df["Amplicon dG"] = df["Amplicon"].apply(ViennaRNAProcessor.calc_deltaG)
        logger.debug("Amplicon deltaG calculation completed successfully")
    except Exception as e:
        error_msg = f"Error calculating amplicon deltaG: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
        
    return df


def run_blast_specificity(df):
    """
    Run BLAST for primer specificity checking.
    
    Executes BLAST analysis for all primers and probes to assess specificity,
    then filters results based on BLAST e-value thresholds.
    
    Args:
        df: DataFrame with primer information
        
    Returns:
        DataFrame with added BLAST results and filtered for specificity,
        or None if no primers pass BLAST filtering
        
    Raises:
        PrimerDesignError: If there's an error in BLAST execution or filtering
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
        error_msg = f"Error in forward primer BLAST: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    
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
        error_msg = f"Error in reverse primer BLAST: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    
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
            error_msg = f"Error in probe BLAST: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise PrimerDesignError(error_msg) from e
    
    # Filter by BLAST specificity
    try:
        logger.debug(f"Filtering by BLAST specificity, filter factor: {Config.BLAST_FILTER_FACTOR}")
        initial_count = len(df)
        df = PrimerProcessor.filter_by_blast(df)
        logger.debug("BLAST filtering completed successfully")
        logger.info(f"After BLAST filtering: {len(df)}/{initial_count} primers")
    except Exception as e:
        error_msg = f"Error in BLAST filtering: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise PrimerDesignError(error_msg) from e
    
    if len(df) == 0:
        logger.warning("No primers passed BLAST filtering.")
        return None
        
    return df