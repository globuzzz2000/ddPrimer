#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Alignment mode for ddPrimer pipeline.

This module contains the implementation of the alignment mode workflow:
1. Load FASTA, VCF, and GFF files for both species
2. Run LastZ alignment or use pre-computed MAF
3. Process alignment and identify conserved regions
4. Hand over to common pipeline for primer design
"""

import os
import logging
import shutil
import tempfile

# Import package modules
from ..config import Config
from ..utils import FileUtils
from ..core import AnnotationProcessor
from ..alignment import AlignmentWorkflow
from . import common  # Import common module functions

# Set up logging
logger = logging.getLogger("ddPrimer")


def run(args):
    """
    Run the alignment mode primer design workflow.
    
    Args:
        args: Command line arguments
        
    Returns:
        bool: Success or failure
    """
    logger.info("=== Alignment Mode Primer Design ===")
    
    try:
        # Check and get required input files
        if args.maf_file:
            # Using pre-computed MAF, only need VCF and GFF files
            logger.info("\n>>> Using pre-computed MAF file <<<")
            reference_file = args.maf_file
            
            if not args.vcf:
                logger.info("\n>>> Please select REFERENCE species VCF file <<<")
                try:
                    args.vcf = FileUtils.get_file(
                        "Select REFERENCE species VCF file", 
                        [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting reference VCF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            
            if not args.second_vcf:
                logger.info("\n>>> Please select SECOND species VCF file <<<")
                try:
                    args.second_vcf = FileUtils.get_file(
                        "Select SECOND species VCF file", 
                        [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting second species VCF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            
            if not args.gff:
                logger.info("\n>>> Please select GFF annotation file <<<")
                try:
                    args.gff = FileUtils.get_file(
                        "Select GFF annotation file", 
                        [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting GFF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
        else:
            # Need all files for alignment workflow
            if not args.fasta:
                logger.info("\n>>> Please select REFERENCE species FASTA file <<<")
                try:
                    args.fasta = FileUtils.get_file(
                        "Select REFERENCE species FASTA file", 
                        [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                    )
                    reference_file = args.fasta
                except Exception as e:
                    logger.error(f"Error selecting reference FASTA file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            else:
                reference_file = args.fasta
            
            if not args.vcf:
                logger.info("\n>>> Please select REFERENCE species VCF file <<<")
                try:
                    args.vcf = FileUtils.get_file(
                        "Select REFERENCE species VCF file", 
                        [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting reference VCF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            
            if not args.second_fasta:
                logger.info("\n>>> Please select SECOND species FASTA file <<<")
                try:
                    args.second_fasta = FileUtils.get_file(
                        "Select SECOND species FASTA file", 
                        [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting second species FASTA file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            
            if not args.second_vcf:
                logger.info("\n>>> Please select SECOND species VCF file <<<")
                try:
                    args.second_vcf = FileUtils.get_file(
                        "Select SECOND species VCF file", 
                        [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting second species VCF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            
            if not args.gff:
                logger.info("\n>>> Please select GFF annotation file <<<")
                try:
                    args.gff = FileUtils.get_file(
                        "Select GFF annotation file", 
                        [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting GFF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
        
        # Set up output directory
        if args.output:
            output_dir = args.output
        elif reference_file:
            if args.maf_file:
                # For MAF files, use the parent directory of the "Alignments" folder
                maf_dir = os.path.dirname(os.path.abspath(args.maf_file))
                if os.path.basename(maf_dir) == "Alignments":
                    # If it's in an Alignments folder, go up one level
                    output_dir = os.path.join(os.path.dirname(maf_dir), "Primers")
                else:
                    # Otherwise use the parent directory of the MAF file
                    output_dir = os.path.join(os.path.dirname(maf_dir), "Primers")
            else:
                # Use the directory of the reference file for non-MAF files
                input_dir = os.path.dirname(os.path.abspath(reference_file))
                output_dir = os.path.join(input_dir, "Primers")
        else:
            # Fallback to current directory if no reference file
            output_dir = os.path.join(os.getcwd(), "Primers")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logger.debug(f"Created output directory: {output_dir}")
        
        # Create temporary directory for intermediate files if needed
        temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
        logger.debug(f"Created temporary directory: {temp_dir}")
        
        try:
            # Call the AlignmentWorkflow function to handle the alignment and masking
            logger.info("\n>>> Running Alignment Primer Design Workflow <<<")
            try:
                masked_sequences, coordinate_map = AlignmentWorkflow(args, output_dir, logger)
                logger.info(f"Alignment workflow completed with {len(masked_sequences)} masked sequences")
            except Exception as e:
                logger.error(f"Error in alignment workflow: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
            
            if not masked_sequences:
                logger.warning("No masked sequences were generated. Exiting.")
                return False
            
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
                reference_file=reference_file,
                mode='alignment',
                genes=genes,
                coordinate_map=coordinate_map,
                gff_file=args.gff,
                temp_dir=temp_dir
            )
            
            return success
            
        finally:
            # Clean up temporary directory
            try:
                if temp_dir and os.path.exists(temp_dir):
                    logger.debug(f"Cleaning up temporary directory: {temp_dir}")
                    shutil.rmtree(temp_dir)
            except Exception as e:
                logger.warning(f"Error cleaning up temporary files: {e}")
            
    except Exception as e:
        logger.error(f"Error in alignment mode workflow: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False