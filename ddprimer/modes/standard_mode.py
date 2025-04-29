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
from ..config.exceptions import FileSelectionError, SequenceProcessingError
from ..utils.file_io import FileIO
from ..core import (
    SNPMaskingProcessor,
    AnnotationProcessor
)
from . import common  # Import common module functions

# Set up logger
logger = logging.getLogger("ddPrimer")


def run(args):
    """
    Run the standard mode primer design workflow.
    
    Args:
        args (argparse.Namespace): Command line arguments
        
    Returns:
        bool: True if the workflow completed successfully, False otherwise
    """
    logger.info("=== Standard Mode Workflow ===")
    
    try:
        # Get input files if not provided in args
        if not args.fasta:
            logger.info("\n>>> Please select FASTA sequence file <<<")
            try:
                args.fasta = FileIO.select_fasta_file("Select FASTA sequence file")
                logger.debug(f"Selected FASTA file: {args.fasta}")
            except FileSelectionError as e:
                logger.error(f"FASTA file selection failed: {str(e)}")
                return False
        
        # VCF file selection
        if not args.vcf:
            logger.info("\n>>> Please select VCF variant file <<<")
            try:
                args.vcf = FileIO.select_file(
                    "Select VCF variant file", 
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                )
                logger.debug(f"Selected VCF file: {args.vcf}")
            except FileSelectionError as e:
                logger.error(f"VCF file selection failed: {str(e)}")
                return False
        
        # GFF file selection
        if not args.noannotation and not args.gff:
            logger.info("\n>>> Please select GFF annotation file <<<")
            try:
                args.gff = FileIO.select_file(
                    "Select GFF annotation file", 
                    [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
                )
                logger.debug(f"Selected GFF file: {args.gff}")
            except FileSelectionError as e:
                logger.error(f"GFF file selection failed: {str(e)}")
                return False
        elif args.noannotation:
            logger.info("\nSkipping GFF annotation file selection")
            # Set args.gff to None so it's consistent
            args.gff = None
        
        # Signal that all file selections are complete
        FileIO.mark_selection_complete()
        
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
        
        # Extract variants from VCF file
        logger.info("\nExtracting variants from VCF file...")
        try:
            snp_processor = SNPMaskingProcessor()
            variants = snp_processor.get_variant_positions(args.vcf)
            logger.debug(f"Variants extracted successfully from {args.vcf}")
            
            total_variants = sum(len(positions) for positions in variants.values())
            logger.debug(f"Extracted {total_variants} variants from {len(variants)} chromosomes")
        except Exception as e:
            logger.error(f"Error extracting variants: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(f"Failed to extract variants: {str(e)}")
        
        # Load sequences from FASTA file
        logger.info("Loading sequences from FASTA file...")
        try:
            sequences = FileIO.load_fasta(args.fasta)
            logger.debug(f"Sequences loaded successfully from {args.fasta}")
            logger.debug(f"Loaded {len(sequences)} sequences")
        except Exception as e:
            logger.error(f"Error loading FASTA: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        # Mask variants in sequences
        logger.debug("\nMasking variants in sequences...")
        
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
                    logger.error(f"Error masking variants in {seq_id}: {str(e)}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise SequenceProcessingError(f"Failed to mask variants in {seq_id}: {str(e)}")
            else:
                logger.debug(f"No variants to mask in {seq_id}")
                masked_sequences[seq_id] = sequence
        
        # Load gene annotations if needed
        if not args.noannotation:
            logger.info("\nLoading gene annotations from GFF file...")
            try:
                genes = AnnotationProcessor.load_genes_from_gff(args.gff)
                logger.debug(f"Gene annotations loaded successfully from {args.gff}")
                logger.debug(f"Loaded {len(genes)} gene annotations")
            except Exception as e:
                logger.error(f"Error loading gene annotations: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
        else:
            logger.info("\nSkipping gene annotation loading (--noannotation specified)")
            genes = None  # Set to None when annotation filtering is disabled
            
        # Use the common workflow function to handle the rest of the pipeline
        success = common.run_primer_design_workflow(
            masked_sequences=masked_sequences,
            output_dir=output_dir,
            reference_file=args.fasta,
            mode='standard',
            genes=genes,
            coordinate_map=None,
            gff_file=args.gff,
            skip_annotation_filtering=args.noannotation
        )
        
        return success
            
    except SequenceProcessingError as e:
        logger.error(f"Sequence processing error: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False
    except Exception as e:
        logger.error(f"Error in standard mode workflow: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False