#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standard mode for ddPrimer pipeline.

Contains functionality for:
1. File selection and validation for standard mode inputs
2. VCF variant extraction with chromosome mapping
3. Sequence processing with SNP masking and fixed SNP substitution
4. Gene annotation loading and processing
5. Integration with common primer design workflow

This module integrates with the broader ddPrimer pipeline to provide
the standard primer design workflow for genomic sequences.
"""

import os
import logging

# Import package modules
from ..config import Config, FileSelectionError, SequenceProcessingError
from ..utils import FileIO
from ..core import SNPMaskingProcessor, AnnotationProcessor
from . import common

# Set up module logger
logger = logging.getLogger(__name__)


def run(args):
    """
    Run the standard mode primer design workflow.
    
    Executes the complete standard mode pipeline including file selection,
    variant extraction, sequence processing, gene annotation loading, and
    primer design workflow execution. Supports fixed SNP substitution.
    
    Args:
        args: Command line arguments containing file paths and options
        
    Returns:
        True if the workflow completed successfully, False otherwise
        
    Raises:
        FileSelectionError: If required input files cannot be selected
        SequenceProcessingError: If sequence processing fails
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
                error_msg = f"FASTA file selection failed: {str(e)}"
                logger.error(error_msg)
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
                error_msg = f"VCF file selection failed: {str(e)}"
                logger.error(error_msg)
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
                error_msg = f"GFF file selection failed: {str(e)}"
                logger.error(error_msg)
                return False
        elif args.noannotation:
            logger.info("\nSkipping GFF annotation file selection")
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
            
            # Extract variants with chromosome mapping
            variants = snp_processor.process_vcf_with_chromosome_mapping(
                args.vcf,
                args.fasta,
                min_af=Config.SNP_ALLELE_FREQUENCY_THRESHOLD,
                min_qual=Config.SNP_QUALITY_THRESHOLD
            )
            
            logger.debug(f"Variants extracted successfully from {args.vcf}")
            
            # Count variants and fixed SNPs
            total_variants = sum(len(variant_list) for variant_list in variants.values())
            fixed_variants = sum(sum(1 for v in variant_list if v.is_fixed) 
                               for variant_list in variants.values())
            logger.debug(f"Extracted {total_variants} variants from {len(variants)} chromosomes")
            logger.debug(f"Fixed SNPs identified: {fixed_variants}")
                
        except Exception as e:
            error_msg = f"Error extracting variants: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
        
        # Load sequences from FASTA file
        logger.info("Loading sequences from FASTA file...")
        try:
            sequences = FileIO.load_fasta(args.fasta)
            logger.debug(f"Sequences loaded successfully from {args.fasta}")
            logger.debug(f"Loaded {len(sequences)} sequences")
        except Exception as e:
            error_msg = f"Error loading FASTA: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        # Process variants in sequences (masking + substitution)
        logger.debug("\nProcessing variants in sequences (masking + fixed SNP substitution)...")
        try:
            masked_sequences = snp_processor.mask_sequences_for_primer_design(
                sequences=sequences,
                variants=variants,
                flanking_size=Config.SNP_FLANKING_MASK_SIZE,
                use_soft_masking=Config.SNP_USE_SOFT_MASKING,
                min_af=Config.SNP_ALLELE_FREQUENCY_THRESHOLD,
                min_qual=Config.SNP_QUALITY_THRESHOLD
            )
            logger.debug("Variant processing completed successfully")
        except Exception as e:
            error_msg = f"Error processing variants: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
        
        # Load gene annotations if needed
        if not args.noannotation:
            logger.info("\nLoading gene annotations from GFF file...")
            try:
                genes = AnnotationProcessor.load_genes_from_gff(args.gff)
                logger.debug(f"Gene annotations loaded successfully from {args.gff}")
                logger.debug(f"Loaded {len(genes)} gene annotations")
            except Exception as e:
                error_msg = f"Error loading gene annotations: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
        else:
            logger.info("\nSkipping gene annotation loading (--noannotation specified)")
            genes = None
            
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
        error_msg = f"Sequence processing error: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False
    except Exception as e:
        error_msg = f"Error in standard mode workflow: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False