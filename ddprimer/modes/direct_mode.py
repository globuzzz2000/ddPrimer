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
from ..utils import FileUtils
from ..core import SNPMaskingProcessor
from . import common

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
        else:
            logger.info("\n>>> SNP masking is disabled <<<")
        
        # Load sequences directly from the provided file
        logger.info("\nLoading sequences from input file...")
        try:
            sequences = FileUtils.load_sequences_from_table(sequence_file)
            logger.debug(f"Loaded {len(sequences)} sequences from {sequence_file}")
        except Exception as e:
            logger.error(f"Error loading sequences from file: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        if not sequences:
            logger.warning("No sequences found in the input file. Exiting.")
            return False
        
        # Apply SNP masking if enabled
        masked_sequences = {}
        if args.snp and ref_fasta and ref_vcf:
            logger.info(f"\nMasking SNPs using reference FASTA and VCF...")
            logger.info(f"Reference FASTA: {ref_fasta}")
            logger.info(f"Reference VCF: {ref_vcf}")
            
            # Initialize SNP masking processor
            snp_processor = SNPMaskingProcessor()
            
            try:
                # Extract variants from VCF
                variants = snp_processor.get_variant_positions(ref_vcf)
                logger.debug(f"Extracted variants from VCF file")
                
                # Apply masking to each sequence
                for seq_id, sequence in sequences.items():
                    # For direct mode, we don't know the chromosome context,
                    # so we'll use a simplified approach just for masking
                    base_variants = set()
                    for chrom_variants in variants.values():
                        # For each position in any chromosome, check if it's within sequence length
                        for pos in chrom_variants:
                            if pos <= len(sequence):
                                base_variants.add(pos)
                    
                    # Mask the sequence
                    if base_variants:
                        logger.debug(f"Masking variants in sequence {seq_id}")
                        masked_sequence = snp_processor.mask_variants(sequence, base_variants)
                        masked_sequences[seq_id] = masked_sequence
                    else:
                        logger.debug(f"No variants to mask in sequence {seq_id}")
                        masked_sequences[seq_id] = sequence
                
                logger.info(f"SNP masking completed for {len(masked_sequences)} sequences")
            except Exception as e:
                logger.error(f"Error during SNP masking: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                logger.warning("Using original sequences without SNP masking")
                masked_sequences = sequences
        else:
            # If SNP masking is disabled, use original sequences
            masked_sequences = sequences
        
        # Debug logging for sequences
        if Config.DEBUG_MODE:
            logger.debug("\n=== DIRECT MODE SEQUENCE INFO ===")
            for i, (seq_id, seq) in enumerate(masked_sequences.items()):
                if i < 5 or i >= len(masked_sequences) - 5:  # Log first and last 5 sequences
                    logger.debug(f"  {seq_id}: length={len(seq)} bp")
                elif i == 5 and len(masked_sequences) > 10:
                    logger.debug(f"  ... ({len(masked_sequences) - 10} more sequences) ...")
            logger.debug("=== END SEQUENCE INFO ===\n")
        
        # Use the common workflow function to handle the rest
        # Note: passing None for genes as gene filtering is not needed in direct mode
        # and None for coordinate_map as it's not needed for direct mode
        success = common.run_primer_design_workflow(
            masked_sequences=masked_sequences,
            output_dir=output_dir,
            reference_file=sequence_file,
            mode='direct',
            genes=None,
            coordinate_map=None,
            gff_file=None,
            skip_annotation_filtering=args.noannotation
        )
        
        return success
            
    except Exception as e:
        logger.error(f"Error in direct mode workflow: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False