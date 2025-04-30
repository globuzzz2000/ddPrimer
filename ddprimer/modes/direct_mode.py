#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct mode for ddPrimer pipeline.

This module contains the implementation of the direct mode workflow:
1. Load sequences directly from CSV or Excel files
2. Run primer design on these sequences using the common pipeline
3. Optionally mask SNPs in sequences if --snp is enabled
"""

import os
import logging
import pandas as pd

# Import package modules
from ..config import Config
from ..config.exceptions import FileSelectionError, SequenceProcessingError
from ..utils.file_io import FileIO
from ..core import SNPMaskingProcessor
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
        
        # Set up output directory
        output_dir = _setup_output_directory(args, sequence_file)
        
        # Before loading, analyze the file structure to provide useful info to the user
        logger.info(f"\nAnalyzing file structure: {sequence_file}")
        analysis = SequenceAnalyzer.analyze_file(sequence_file)
        SequenceAnalyzer.print_analysis(analysis)
        
        # Dictionary to track reference matching status
        matching_status = {}
        masked_sequences = {}
        all_sequences = {}  # Store all original sequences, including those that fail matching
        
        # Handle SNP masking if enabled
        if args.snp:
            # Process SNP masking workflow
            masked_sequences, matching_status, all_sequences = _process_snp_masking(args, sequence_file)
            if not masked_sequences:
                return False
        else:
            # If SNP masking is disabled, use original sequences
            logger.debug("\n>>> SNP masking is disabled <<<")
            # Signal that all file selections are complete (if we didn't already do it in SNP mode)
            FileIO.mark_selection_complete()
            
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
                    masked_sequences = sequences.copy()  # No masking, so use original
                    
                    # Update matching status for all sequences
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
        if not masked_sequences:
            logger.warning("No sequences found in the input file. Exiting.")
            return False
        
        # Debug logging for sequences
        _log_sequence_info(masked_sequences, matching_status)
        
        # Use the common workflow function to handle the rest
        success = common.run_primer_design_workflow(
            masked_sequences=masked_sequences,
            output_dir=output_dir,
            reference_file=sequence_file,
            mode='direct',
            genes=None,
            coordinate_map=None,
            gff_file=None,
            skip_annotation_filtering=args.noannotation,
            matching_status=matching_status,  # Pass matching status
            all_sequences=all_sequences,  # Pass all sequences including those that failed matching
            add_rows_function=DirectMasking.add_missing_sequences  # Pass the static method
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
        
        # Check if we have any sequences to process
        if not masked_sequences:
            logger.warning("No sequences found in the input file. Exiting.")
            return False
        
        # Debug logging for sequences
        _log_sequence_info(masked_sequences, matching_status)
        
        # Use the common workflow function to handle the rest
        success = common.run_primer_design_workflow(
            masked_sequences=masked_sequences,
            output_dir=output_dir,
            reference_file=sequence_file,
            mode='direct',
            genes=None,
            coordinate_map=None,
            gff_file=None,
            skip_annotation_filtering=args.noannotation,
            matching_status=matching_status,  # Pass matching status
            all_sequences=all_sequences,  # Pass all sequences including those that failed matching
            add_rows_function=DirectMasking.add_missing_sequences  # Pass the static method
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


def _process_snp_masking(args, sequence_file):
    """
    Process SNP masking for direct mode sequences.
    
    Args:
        args (argparse.Namespace): Command line arguments
        sequence_file (str): Path to the sequence file
        
    Returns:
        tuple: (masked_sequences, matching_status, all_sequences)
    """
    logger.debug("\n>>> SNP masking is enabled <<<")
    
    # Get reference FASTA and VCF files
    ref_fasta = args.fasta
    ref_vcf = args.vcf
    
    # If not provided, prompt for them
    if not ref_fasta:
        logger.info("\n>>> Please select FASTA sequence file <<<")
        try:
            ref_fasta = FileIO.select_fasta_file("Select FASTA file")
            logger.debug(f"Selected reference FASTA file: {ref_fasta}")
        except FileSelectionError as e:
            logger.error(f"Reference FASTA file selection failed: {str(e)}")
            args.snp = False
            return {}, {}, {}
    
    if args.snp and not ref_vcf:
        logger.info("\n>>> Please select VCF variant file <<<")
        try:
            ref_vcf = FileIO.select_file(
                "Select VCF variant file", 
                [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
            )
            logger.debug(f"Selected VCF file: {ref_vcf}")
        except FileSelectionError as e:
            logger.error(f"VCF file selection failed: {str(e)}")
            args.snp = False
            return {}, {}, {}
    
    # Signal that all file selections are complete
    FileIO.mark_selection_complete()
    
    # Load sequences directly from the provided file
    logger.info("\nLoading sequences from input file...")
    
    try:
        # First analyze the file to provide better feedback
        from ..helpers.sequence_analyzer import SequenceAnalyzer
        logger.debug("Analyzing sequence file structure...")
        analysis = SequenceAnalyzer.analyze_file(sequence_file)
        SequenceAnalyzer.print_analysis(analysis)
        
        # Get recommended columns
        name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
        if name_col and seq_col:
            logger.debug(f"Using column '{name_col}' for sequence names and '{seq_col}' for sequences")
        
        # Load the sequences
        sequences = FileIO.load_sequences_from_table(sequence_file)
        
        if sequences:
            # Log a sample of the loaded sequences for verification
            sample_size = min(3, len(sequences))
            sample_items = list(sequences.items())[:sample_size]
            
            logger.info(f"Successfully loaded {len(sequences)} sequences")
            logger.debug("Sample of sequences loaded:")
            for seq_id, sequence in sample_items:
                # Truncate sequence for display
                display_seq = f"{sequence[:50]}..." if len(sequence) > 50 else sequence
                logger.debug(f"  {seq_id}: {display_seq}")
            
            all_sequences = sequences.copy()  # Keep a copy of all sequences
        else:
            logger.warning("No valid sequences found in the input file.")
            return {}, {}, {}
            
    except Exception as e:
        logger.error(f"Error loading sequences from file: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return {}, {}, {}
    
    if not sequences:
        logger.warning("No sequences found in the input file. Exiting.")
        return {}, {}, {}
        
    # Apply SNP masking if enabled and required files are present
    if args.snp and ref_fasta and ref_vcf:
        return _mask_sequences_with_snps(sequences, ref_fasta, ref_vcf, all_sequences)
    else:
        # If SNP masking failed to initialize, use original sequences
        masked_sequences = sequences
        matching_status = {seq_id: "Not attempted" for seq_id in sequences}
        return masked_sequences, matching_status, all_sequences

def _mask_sequences_with_snps(sequences, ref_fasta, ref_vcf, all_sequences):
    """
    Apply SNP masking to sequences.
    
    Args:
        sequences (dict): Dictionary of sequences
        ref_fasta (str): Path to reference FASTA file
        ref_vcf (str): Path to VCF file
        all_sequences (dict): Dictionary of all sequences
        
    Returns:
        tuple: (masked_sequences, matching_status, all_sequences)
    """
    logger.info(f"\nMasking SNPs using reference FASTA and VCF...")
    logger.debug(f"Reference FASTA: {ref_fasta}")
    logger.debug(f"Reference VCF: {ref_vcf}")
    
    # Initialize processors
    snp_processor = SNPMaskingProcessor()
    matching_status = {}
    masked_sequences = {}
    
    try:
        # Process each sequence individually
        for seq_id, sequence in sequences.items():
            logger.debug(f"Processing sequence {seq_id} ({len(sequence)} bp)")
            
            # Find sequence location in reference genome using BLAST
            source_chrom, start_pos, end_pos, identity = DirectMasking.find_location(
                sequence, ref_fasta
            )
            
            # If sequence was successfully matched to reference
            if source_chrom:
                matching_status[seq_id] = "Success"
                logger.debug(f"Matched sequence {seq_id} to {source_chrom}:{start_pos}-{end_pos} "
                           f"with {identity}% identity")
                
                # Extract only the variants for this specific region using SNPMaskingProcessor
                region_variants = snp_processor.get_region_variants(ref_vcf, source_chrom, start_pos, end_pos)
                
                # Adjust variant positions to sequence coordinates
                adjusted_variants = set()
                for var_pos in region_variants:
                    # Convert 1‑based genome coordinate to 1‑based sequence coordinate
                    seq_pos = (var_pos - start_pos) + 1
                    adjusted_variants.add(seq_pos)
                
                # Mask the sequence with adjusted variants
                if adjusted_variants:
                    logger.debug(f"Masking sequence {seq_id} with {len(adjusted_variants)} variants")
                    masked_sequence = snp_processor.mask_variants(sequence, adjusted_variants)
                    masked_sequences[seq_id] = masked_sequence
                else:
                    logger.debug(f"No variants found in region, using original sequence")
                    masked_sequences[seq_id] = sequence
            else:
                # Sequence could not be matched to reference
                matching_status[seq_id] = "Failure"
                logger.debug(f"Could not match sequence {seq_id} to reference genome")
                # We'll exclude this sequence from masked_sequences
        
        logger.debug(f"SNP masking completed")
        logger.info(f"Successfully matched {sum(1 for s in matching_status.values() if s == 'Success')}/{len(sequences)} sequences to reference")
        logger.debug(f"Failed to match {sum(1 for s in matching_status.values() if s == 'Failure')}/{len(sequences)} sequences to reference")
        
        return masked_sequences, matching_status, all_sequences
        
    except Exception as e:
        logger.error(f"Error during SNP masking: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        logger.warning("Using original sequences without SNP masking")
        
        # If SNP masking fails, use original sequences
        masked_sequences = sequences
        matching_status = {seq_id: "Not attempted" for seq_id in sequences}
        
        return masked_sequences, matching_status, all_sequences


def _log_sequence_info(masked_sequences, matching_status):
    """
    Log information about sequences for debugging.
    
    Args:
        masked_sequences (dict): Dictionary of masked sequences
        matching_status (dict): Dictionary of matching statuses
    """
    if Config.DEBUG_MODE:
        logger.debug("\n=== DIRECT MODE SEQUENCE INFO ===")
        for i, (seq_id, seq) in enumerate(masked_sequences.items()):
            if i < 5 or i >= len(masked_sequences) - 5:  # Log first and last 5 sequences
                logger.debug(f"  {seq_id}: length={len(seq)} bp, status={matching_status[seq_id]}")
            elif i == 5 and len(masked_sequences) > 10:
                logger.debug(f"  ... ({len(masked_sequences) - 10} more sequences) ...")
        logger.debug("=== END SEQUENCE INFO ===\n")