#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct mode for ddPrimer pipeline with fixed VCF prompting.

This module contains the implementation of the direct mode workflow:
1. Load sequences directly from CSV or Excel files
2. Run primer design on these sequences using the common pipeline
3. Optionally filter primers based on SNP positions if --snp is enabled
"""

import os
import logging
import pandas as pd

# Import package modules
from ..config import Config
from ..config.exceptions import FileSelectionError, SequenceProcessingError
from ..utils.file_io import FileIO
from ..core import SNPFilterProcessor
from . import common
from ..helpers import DirectMasking
from ..helpers.sequence_analyzer import SequenceAnalyzer  # Updated import for SequenceAnalyzer

# Set up logger
logger = logging.getLogger("ddPrimer")


def run(args):
    """
    Run the direct mode primer design workflow with improved file handling.
    
    Args:
        args (argparse.Namespace): Command line arguments
        
    Returns:
        bool: True if the workflow completed successfully, False otherwise
    """
    logger.info("=== Direct Mode Workflow ===")
    
    try:
        # Get the input file - handle both when --direct is specified alone or with a file
        sequence_file = None
        if args.direct is True:  # --direct was used without a specified file
            logger.info("\nSelecting sequence file (CSV or Excel)")
            try:
                sequence_file = FileIO.select_sequences_file()
                logger.debug(f"Selected sequence file: {sequence_file}")
            except FileSelectionError as e:
                logger.error(f"Sequence file selection failed: {str(e)}")
                return False
        else:  # --direct was used with a specified file
            sequence_file = args.direct
            logger.debug(f"Using specified sequence file: {sequence_file}")
        
        # Check if SNP filtering is enabled
        snp_filtering_enabled = hasattr(args, 'snp_enabled') and args.snp_enabled
        
        # Get reference FASTA and VCF files if SNP filtering is enabled
        ref_fasta = args.fasta
        ref_vcf = args.vcf
        variants = None
        
        # FIXED: Explicitly check for SNP filtering before checking if VCF exists
        if snp_filtering_enabled:
            logger.debug("\n>>> SNP filtering is enabled <<<")
            
            # If not provided, prompt for reference files
            if not ref_fasta:
                logger.info("\n>>> Please select FASTA reference file <<<")
                try:
                    ref_fasta = FileIO.select_fasta_file("Select FASTA reference file")
                    logger.debug(f"Selected reference FASTA file: {ref_fasta}")
                except FileSelectionError as e:
                    logger.error(f"Reference FASTA file selection failed: {str(e)}")
                    # Disable SNP filtering if FASTA selection fails
                    args.snp_enabled = False
                    logger.warning("SNP filtering has been disabled due to missing reference FASTA file")
            
            if args.snp_enabled and not ref_vcf:
                logger.info("\n>>> Please select VCF variant file <<<")
                try:
                    ref_vcf = FileIO.select_file(
                        "Select VCF variant file", 
                        [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                    )
                    logger.debug(f"Selected VCF file: {ref_vcf}")
                except FileSelectionError as e:
                    logger.error(f"VCF file selection failed: {str(e)}")
                    # Disable SNP filtering if VCF selection fails
                    args.snp_enabled = False
                    logger.warning("SNP filtering has been disabled due to missing VCF file")
            
            # Extract variants from VCF if all required files are available
            if args.snp_enabled and ref_fasta and ref_vcf:
                try:
                    snp_processor = SNPFilterProcessor()
                    variants = snp_processor.get_variant_positions(ref_vcf)
                    logger.debug(f"Variants extracted successfully from {ref_vcf}")
                    
                    total_variants = sum(len(positions) for positions in variants.values())
                    logger.debug(f"Extracted {total_variants} variants from {len(variants)} chromosomes")
                except Exception as e:
                    logger.error(f"Error extracting variants: {str(e)}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    # Disable SNP filtering if variant extraction fails
                    args.snp_enabled = False
                    logger.warning("SNP filtering has been disabled due to error extracting variants")
        else:
            logger.debug("\n>>> SNP filtering is disabled <<<")
            
        # Signal that all file selections are complete
        FileIO.mark_selection_complete()
        
        # Set up output directory
        output_dir = _setup_output_directory(args, sequence_file)
        
        # Before loading, analyze the file structure to provide useful info to the user
        logger.info(f"\nAnalyzing file structure: {sequence_file}")
        analysis = SequenceAnalyzer.analyze_file(sequence_file)
        SequenceAnalyzer.print_analysis(analysis)
        
        # Dictionary to track reference matching status
        matching_status = {}
        sequences = {}
        all_sequences = {}  # Store all original sequences
        
        # Load sequences directly from the provided file
        logger.info("\nLoading sequences from input file...")
        try:
            sequences = FileIO.load_sequences_from_table(sequence_file)
            
            if sequences:
                logger.debug(f"Loaded {len(sequences)} sequences from {sequence_file}")
                
                # Log sample of loaded sequences
                sample_size = min(3, len(sequences))
                sample_items = list(sequences.items())[:sample_size]
                logger.info(f"Successfully loaded {len(sequences)} sequences")
                logger.debug("Sample of loaded sequences:")
                for seq_id, sequence in sample_items:
                    # Truncate sequence for display
                    display_seq = f"{sequence[:30]}..." if len(sequence) > 30 else sequence
                    logger.debug(f"  {seq_id}: {display_seq}")
                
                all_sequences = sequences.copy()  # Keep a copy of all sequences
                
                # Initialize matching status for all sequences
                for seq_id in sequences:
                    matching_status[seq_id] = "Not attempted"
            else:
                logger.error("No valid sequences found in the input file")
                return False
        except Exception as e:
            logger.error(f"Error loading sequences from file: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        # Check if we have any sequences to process
        if not sequences:
            logger.warning("No sequences found in the input file. Exiting.")
            return False
        
        # For SNP filtering in direct mode, we need to map sequences to the reference
        # to get proper coordinates for variant filtering if SNP filtering is enabled
        if args.snp_enabled and ref_fasta and variants:
            sequence_mapping = {}
            
            # Process each sequence to find its location in the reference genome
            for seq_id, sequence in sequences.items():
                # Find sequence location in reference genome using BLAST
                source_chrom, start_pos, end_pos, identity = DirectMasking.find_location(
                    sequence, ref_fasta
                )
                
                # If sequence was successfully matched to reference
                if source_chrom:
                    matching_status[seq_id] = "Success"
                    logger.debug(f"Matched sequence {seq_id} to {source_chrom}:{start_pos}-{end_pos} "
                               f"with {identity}% identity")
                    
                    # Store mapping information
                    sequence_mapping[seq_id] = {
                        "chr": source_chrom,
                        "start": start_pos,
                        "end": end_pos,
                        "identity": identity
                    }
                else:
                    # Sequence could not be matched to reference
                    matching_status[seq_id] = "Failure"
                    logger.debug(f"Could not match sequence {seq_id} to reference genome")
            
            # Log mapping statistics
            success_count = sum(1 for s in matching_status.values() if s == "Success")
            logger.info(f"Successfully matched {success_count}/{len(sequences)} sequences to reference")
            logger.debug(f"Failed to match {len(sequences) - success_count}/{len(sequences)} sequences to reference")
            
            # Adjust variant positions based on sequence mapping
            if sequence_mapping:
                # Create a new variants dictionary with adjusted positions
                adjusted_variants = {}
                
                for seq_id, mapping in sequence_mapping.items():
                    chr_name = mapping["chr"]
                    if chr_name in variants:
                        # Get variants for this chromosome
                        chr_variants = variants[chr_name]
                        
                        # Filter to only include variants in this sequence's region
                        region_start = mapping["start"]
                        region_end = mapping["end"]
                        
                        # Adjust variant positions to be relative to sequence start
                        seq_variants = set()
                        for var_pos in chr_variants:
                            if region_start <= var_pos <= region_end:
                                # Convert to sequence-relative position (1-based)
                                seq_pos = var_pos - region_start + 1
                                seq_variants.add(seq_pos)
                        
                        if seq_variants:
                            adjusted_variants[seq_id] = seq_variants
                
                # Use the adjusted variants for filtering
                variants = adjusted_variants
                
        # Debug logging for sequences
        _log_sequence_info(sequences, matching_status)
        
        # Use the common workflow function to handle the primer design
        success = common.run_primer_design_workflow(
            masked_sequences=sequences,
            output_dir=output_dir,
            reference_file=sequence_file,
            mode='direct',
            genes=None,
            coordinate_map=None,
            gff_file=None,
            skip_annotation_filtering=args.noannotation,
            matching_status=matching_status,
            all_sequences=all_sequences,
            add_rows_function=DirectMasking.add_missing_sequences,
            # Pass SNP filtering parameters
            variants=variants,
            strict_mode=args.snp_strict if hasattr(args, 'snp_strict') else False,
            amplicon_threshold=args.snp_threshold if hasattr(args, 'snp_threshold') else 0.05
        )
        
        return success
            
    except SequenceProcessingError as e:
        logger.error(f"Sequence processing error: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False
    except Exception as e:
        logger.error(f"Error in direct mode workflow: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False


def _setup_output_directory(args, sequence_file):
    """
    Set up the output directory for direct mode.
    
    Args:
        args (argparse.Namespace): Command line arguments
        sequence_file (str): Path to the sequence file
        
    Returns:
        str: Path to the output directory
    """
    if args.output:
        output_dir = args.output
    else:
        # Use the directory of the reference file
        input_dir = os.path.dirname(os.path.abspath(sequence_file))
        output_dir = os.path.join(input_dir, "Primers")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logger.debug(f"Created output directory: {output_dir}")
    
    return output_dir


def _log_sequence_info(sequences, matching_status):
    """
    Log information about sequences for debugging.
    
    Args:
        sequences (dict): Dictionary of sequences
        matching_status (dict): Dictionary of matching statuses
    """
    if Config.DEBUG_MODE:
        logger.debug("\n=== DIRECT MODE SEQUENCE INFO ===")
        for i, (seq_id, seq) in enumerate(sequences.items()):
            if i < 5 or i >= len(sequences) - 5:  # Log first and last 5 sequences
                logger.debug(f"  {seq_id}: length={len(seq)} bp, status={matching_status[seq_id]}")
            elif i == 5 and len(sequences) > 10:
                logger.debug(f"  ... ({len(sequences) - 10} more sequences) ...")
        logger.debug("=== END SEQUENCE INFO ===\n")