#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standard mode for ddPrimer pipeline.

This module contains the implementation of the standard mode workflow:
1. Load FASTA, VCF and GFF files
2. Extract variants from VCF and mask the FASTA file
3. Pass masked sequences to the common primer design workflow
"""

import os
import logging
from tqdm import tqdm

# Import package modules
from ..config import Config
from ..utils import FileUtils
from ..core import (
    SNPMaskingProcessor,
    AnnotationProcessor
)
from . import common  # Import common module functions

# Set up logging
logger = logging.getLogger("ddPrimer")


def run(args):
    """
    Run the standard mode primer design workflow.
    
    Args:
        args: Command line arguments
        
    Returns:
        bool: Success or failure
    """
    logger.info("=== Standard Mode Workflow ===")
    
    try:
        # Get input files if not provided in args
        if not args.fasta:
            logger.info("\n>>> Please select reference FASTA file <<<")
            try:
                args.fasta = FileUtils.get_file(
                    "Select reference FASTA file", 
                    [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                )
                logger.debug(f"FASTA file selection successful: {args.fasta}")
            except Exception as e:
                logger.error(f"Error selecting FASTA file: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
        
        # VCF file selection
        if not args.vcf:
            logger.info("\n>>> Please select VCF file with variants <<<")
            try:
                args.vcf = FileUtils.get_file(
                    "Select VCF file with variants", 
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                )
                logger.debug(f"VCF file selection successful: {args.vcf}")
            except Exception as e:
                logger.error(f"Error selecting VCF file: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
        
        # GFF file selection
        if not args.gff:
            logger.info("\n>>> Please select GFF annotation file <<<")
            try:
                args.gff = FileUtils.get_file(
                    "Select GFF annotation file", 
                    [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
                )
                logger.debug(f"GFF file selection successful: {args.gff}")
            except Exception as e:
                logger.error(f"Error selecting GFF file: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
        
        # Set up output directory
        if args.output:
            output_dir = args.output
        else:
            # Use the directory of the reference file
            input_dir = os.path.dirname(os.path.abspath(args.fasta))
            output_dir = os.path.join(input_dir, "Primers")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logger.debug(f"Created output directory: {output_dir}")
        
        # Extract variants and mask sequences
        logger.info("\nExtracting variants from VCF file...")
        snp_processor = SNPMaskingProcessor()
        try:
            variants = snp_processor.get_variant_positions(args.vcf)
            logger.debug(f"Variants extracted successfully from {args.vcf}")
        except Exception as e:
            logger.error(f"Error extracting variants: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        total_variants = sum(len(positions) for positions in variants.values())
        logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes")
        
        logger.info("\nLoading sequences from FASTA file...")
        try:
            sequences = FileUtils.load_fasta(args.fasta)
            logger.debug(f"Sequences loaded successfully from {args.fasta}")
        except Exception as e:
            logger.error(f"Error loading FASTA: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        logger.debug(f"Loaded {len(sequences)} sequences")
        
        logger.info("\nMasking variants in sequences...")
        
        masked_sequences = {}
        seq_items = list(sequences.items())
        if Config.SHOW_PROGRESS:
            seq_iter = tqdm(seq_items, total=len(seq_items), desc="Masking sequences")
        else:
            seq_iter = seq_items
            
        for seq_id, sequence in seq_iter:
            # Get variants for this sequence/chromosome
            seq_variants = variants.get(seq_id, set())
            
            if seq_variants:
                logger.debug(f"Masking {len(seq_variants)} variants in {seq_id}...")
                try:
                    masked_seq = snp_processor.mask_variants(sequence, seq_variants)
                    masked_sequences[seq_id] = masked_seq
                except Exception as e:
                    logger.error(f"Error masking variants in {seq_id}: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            else:
                logger.debug(f"No variants to mask in {seq_id}")
                masked_sequences[seq_id] = sequence
        
        # Load gene annotations
        logger.info("\nLoading gene annotations from GFF file...")
        try:
            genes = AnnotationProcessor.load_genes_from_gff(args.gff)
            logger.debug(f"Gene annotations loaded successfully from {args.gff}")
            logger.info(f"Loaded {len(genes)} gene annotations")
        except Exception as e:
            logger.error(f"Error loading gene annotations: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
            
        # Use the common workflow function to handle the rest of the pipeline
        success = common.run_primer_design_workflow(
            masked_sequences=masked_sequences,
            output_dir=output_dir,
            reference_file=args.fasta,
            mode='standard',
            genes=genes,
            coordinate_map=None,
            gff_file=args.gff
        )
        
        return success
            
    except Exception as e:
        logger.error(f"Error in standard mode workflow: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False