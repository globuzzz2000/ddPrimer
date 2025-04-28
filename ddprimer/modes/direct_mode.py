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
import tempfile
import subprocess
import io
import re

# Import package modules
from ..config import Config
from ..utils import FileUtils
from ..core import SNPMaskingProcessor
from . import common
from ..helpers import DirectMasking

# Set up logging
logger = logging.getLogger("ddPrimer")

def run(args):
    """
    Run the direct mode primer design workflow.
    
    Args:
        args: Command line arguments
        
    Returns:
        bool: Success or failure
    """
    logger.info("=== Direct Mode Workflow ===")
    
    try:
        # Get the input file - handle both when --direct is specified alone or with a file
        if args.direct is True:  # --direct was used without a specified file
            sequence_file = FileUtils.get_sequences_file()
        else:  # --direct was used with a specified file
            sequence_file = args.direct
            
        logger.debug(f"Direct sequence file: {sequence_file}")
        
        # Set up output directory
        if args.output:
            output_dir = args.output
        else:
            # Use the directory of the reference file
            input_dir = os.path.dirname(os.path.abspath(sequence_file))
            output_dir = os.path.join(input_dir, "Primers")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logger.debug(f"Created output directory: {output_dir}")
        
        # Dictionary to track reference matching status
        matching_status = {}
        masked_sequences = {}
        all_sequences = {}  # Store all original sequences, including those that fail matching
        
        # Check if SNP masking is enabled
        if args.snp:
            logger.info("\n>>> SNP masking is enabled <<<")
            
            # Get reference FASTA and VCF files
            ref_fasta = args.fasta
            ref_vcf = args.vcf
            
            # If not provided, prompt for them
            if not ref_fasta:
                logger.info("\n>>> Please select reference FASTA file <<<")
                try:
                    ref_fasta = FileUtils.get_file(
                        "Select reference FASTA file", 
                        [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting reference FASTA file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    args.snp = False
            
            if args.snp and not ref_vcf:
                logger.info("\n>>> Please select VCF file with variants <<<")
                try:
                    ref_vcf = FileUtils.get_file(
                        "Select VCF file with variants", 
                        [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting VCF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    args.snp = False
            
            # Signal that all file selections are complete
            FileUtils.mark_selection_complete()
        else:
            logger.info("\n>>> SNP masking is disabled <<<")
            # Signal that all file selections are complete
            FileUtils.mark_selection_complete()
        
        # Load sequences directly from the provided file
        logger.info("\nLoading sequences from input file...")
        try:
            sequences = FileUtils.load_sequences_from_table(sequence_file)
            logger.debug(f"Loaded {len(sequences)} sequences from {sequence_file}")
            all_sequences = sequences.copy()  # Keep a copy of all sequences
        except Exception as e:
            logger.error(f"Error loading sequences from file: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        if not sequences:
            logger.warning("No sequences found in the input file. Exiting.")
            return False
        
        # Apply SNP masking if enabled
        if args.snp and ref_fasta and ref_vcf:
            logger.info(f"\nMasking SNPs using reference FASTA and VCF...")
            logger.debug(f"Reference FASTA: {ref_fasta}")
            logger.debug(f"Reference VCF: {ref_vcf}")
            
            # Initialize processors
            snp_processor = SNPMaskingProcessor()
            
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
                
                logger.info(f"SNP masking completed")
                logger.info(f"Successfully matched {sum(1 for s in matching_status.values() if s == 'Success')}/{len(sequences)} sequences to reference")
                logger.debug(f"Failed to match {sum(1 for s in matching_status.values() if s == 'Failure')}/{len(sequences)} sequences to reference")
            except Exception as e:
                logger.error(f"Error during SNP masking: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                logger.warning("Using original sequences without SNP masking")
                masked_sequences = sequences
                for seq_id in sequences:
                    matching_status[seq_id] = "Not attempted"
        else:
            # If SNP masking is disabled, use original sequences
            masked_sequences = sequences
            for seq_id in sequences:
                matching_status[seq_id] = "Not attempted"
        
        # Debug logging for sequences
        if Config.DEBUG_MODE:
            logger.debug("\n=== DIRECT MODE SEQUENCE INFO ===")
            for i, (seq_id, seq) in enumerate(masked_sequences.items()):
                if i < 5 or i >= len(masked_sequences) - 5:  # Log first and last 5 sequences
                    logger.debug(f"  {seq_id}: length={len(seq)} bp, status={matching_status[seq_id]}")
                elif i == 5 and len(masked_sequences) > 10:
                    logger.debug(f"  ... ({len(masked_sequences) - 10} more sequences) ...")
            logger.debug("=== END SEQUENCE INFO ===\n")
        
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
            
    except Exception as e:
        logger.error(f"Error in direct mode workflow: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False