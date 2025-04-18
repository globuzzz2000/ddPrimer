#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct mode for ddPrimer pipeline.

This module contains the implementation of the direct mode workflow:
1. Load sequences directly from CSV or Excel files
2. Run primer design on these sequences using the common pipeline
3. Optionally check primers/probes for SNP overlaps with a reference genome
"""

import os
import logging
import pandas as pd

# Import package modules
from ..config import Config
from ..utils import (FileUtils, SNPChecker)
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
        
        # Check if SNP checking is enabled
        snp_checker = None
        if args.check_snps:
            logger.info("\nSNP checking is enabled")
            
            # Get reference FASTA and VCF files
            ref_fasta = args.fasta
            ref_vcf = args.vcf
            
            # If not provided, prompt for them
            if not ref_fasta:
                ref_fasta = FileUtils.get_file("Select reference FASTA file", 
                                              [("FASTA files", "*.fa *.fasta *.fna")])
                if not ref_fasta:
                    logger.warning("Reference FASTA file not provided. Disabling SNP checking.")
                    args.check_snps = False
            
            if args.check_snps and not ref_vcf:
                ref_vcf = FileUtils.get_file("Select VCF file with variants", 
                                           [("VCF files", "*.vcf *.vcf.gz")])
                if not ref_vcf:
                    logger.warning("VCF file not provided. Disabling SNP checking.")
                    args.check_snps = False
            
            if args.check_snps:
                logger.info(f"Reference FASTA: {ref_fasta}")
                logger.info(f"Reference VCF: {ref_vcf}")
                
                # Initialize SNP checker
                snp_checker = SNPChecker(ref_fasta, ref_vcf)
        
        # Load sequences directly from the provided file
        logger.info("\nLoading sequences from input file...")
        try:
            masked_sequences = FileUtils.load_sequences_from_table(sequence_file)
            logger.debug(f"Loaded {len(masked_sequences)} sequences from {sequence_file}")
        except Exception as e:
            logger.error(f"Error loading sequences from file: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        if not masked_sequences:
            logger.warning("No sequences found in the input file. Exiting.")
            return False
        
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
            snp_checker=snp_checker
        )
        
        return success
            
    except Exception as e:
        logger.error(f"Error in direct mode workflow: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False