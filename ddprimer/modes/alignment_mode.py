#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Alignment mode for ddPrimer pipeline.

This module contains the implementation of the alignment mode workflow:
1. Load FASTA, VCF, and GFF files for both genomes
2. Run LastZ alignment or use pre-computed MAF
3. Process alignment and identify conserved regions
4. Hand over to common pipeline for primer design
"""

import os
import logging

# Import package modules
from ..config import Config, FileSelectionError, AlignmentError
from ..utils import FileIO, TempDirectoryManager
from ..core import AnnotationProcessor
from ..helpers import run_alignment_workflow, LastZRunner
from . import common

# Set up logger
logger = logging.getLogger(__name__)


def run(args):
    """
    Run the alignment mode primer design workflow.
    
    Args:
        args (argparse.Namespace): Command line arguments
        
    Returns:
        bool: True if the workflow completed successfully, False otherwise
    """
    logger.info("=== Alignment Mode Primer Design ===")
    
    try:
        # For LastZ-only mode, we just need reference and second FASTA files
        if args.lastzonly:
            return _run_lastz_only_mode(args)
            
        # Standard alignment mode
        # Check if SNP masking is enabled
        if args.snp:
            logger.debug("\n>>> SNP masking is enabled <<<")
        else:
            logger.debug("\n>>> SNP masking is disabled <<<")
        
        # Check and get required input files
        # Modified check for maf to handle the case when the flag is used without a value
        if hasattr(args, 'maf') and args.maf is not None:
            # Using pre-computed MAF workflow
            return _run_maf_workflow(args)
        else:
            # Using direct alignment workflow
            return _run_direct_alignment_workflow(args)
            
    except AlignmentError as e:
        logger.error(f"Alignment error: {e}")
        return False
    except Exception as e:
        logger.error(f"Error in alignment mode workflow: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False


def _run_lastz_only_mode(args):
    """
    Run the LastZ-only mode (alignment without primer design).
    
    Args:
        args (argparse.Namespace): Command line arguments
        
    Returns:
        bool: True if the alignment completed successfully, False otherwise
    """
    logger.info("=== LastZ-Only Mode ===")
    
    # Only get required input files
    if not args.fasta:
        logger.info("\n>>> Please select REFERENCE FASTA sequence file <<<")
        try:
            args.fasta = FileIO.select_fasta_file("Select REFERENCE FASTA file")
            logger.debug(f"Selected reference FASTA file: {args.fasta}")
        except FileSelectionError as e:
            logger.error(f"Reference FASTA file selection failed: {e}")
            return False
    
    if not args.second_fasta:
        logger.info("\n>>> Please select SECOND FASTA sequence file <<<")
        try:
            args.second_fasta = FileIO.select_fasta_file("Select SECOND FASTA file")
            logger.debug(f"Selected second FASTA file: {args.second_fasta}")
        except FileSelectionError as e:
            logger.error(f"Second FASTA file selection failed: {e}")
            return False
    
    # Signal that all file selections are complete
    FileIO.mark_selection_complete()
    
    # Set up output directory
    if args.output:
        output_dir = args.output
    else:
        # Use the directory of the reference FASTA
        input_dir = os.path.dirname(os.path.abspath(args.fasta))
        output_dir = os.path.join(input_dir, "Alignments")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logger.debug(f"Created output directory: {output_dir}")
    
    # Run only the LastZ portion of the alignment workflow
    logger.info("\n>>> Running LastZ alignment only <<<")
    
    try:
        # Create LastZ runner instance
        runner = LastZRunner()
        
        # Run LastZ alignment
        output_file = runner.run_parallel_alignment(
            args.fasta,
            args.second_fasta,
            output_dir,
            Config.LASTZ_OPTIONS,
            processes=Config.NUM_PROCESSES,
            keep_temp=args.debug  # Keep temp files in debug mode
        )
        
        logger.info(f"LastZ alignment completed successfully!")
        logger.info(f"Alignment file saved to: {output_file}")
        return True
        
    except Exception as e:
        logger.error(f"Error running LastZ alignment: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise AlignmentError(f"LastZ alignment failed: {e}")


def _run_maf_workflow(args):
    """
    Run the MAF-based alignment workflow (using pre-computed MAF file).
    
    Args:
        args (argparse.Namespace): Command line arguments
        
    Returns:
        bool: True if the workflow completed successfully, False otherwise
    """
    logger.debug("\n>>> Using pre-computed MAF file workflow <<<")
    
    # If --maf was used but no value provided, prompt for file selection
    if args.maf is True or args.maf == '':
        logger.info("\n>>> Please select MAF alignment file <<<")
        try:
            args.maf = FileIO.select_file(
                "Select MAF alignment file", 
                [("MAF Files", "*.maf"), ("All Files", "*")]
            )
            logger.debug(f"Selected MAF file: {args.maf}")
        except FileSelectionError as e:
            logger.error(f"MAF file selection failed: {e}")
            return False
    
    # Using pre-computed MAF, only need VCF and GFF files
    logger.debug("\n>>> Using pre-computed MAF file <<<")
    reference_file = args.maf
    
    if args.snp:
        # Only require VCF files if SNP masking is enabled
        if not args.vcf:
            logger.info("\n>>> Please select REFERENCE VCF variant file <<<")
            try:
                args.vcf = FileIO.select_file(
                    "Select REFERENCE VCF file", 
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                )
                logger.debug(f"Selected reference VCF variant file: {args.vcf}")
            except FileSelectionError as e:
                logger.error(f"Reference VCF file selection failed: {e}")
                return False
        
        if not args.second_vcf:
            logger.info("\n>>> Please select SECOND VCF variant file <<<")
            try:
                args.second_vcf = FileIO.select_file(
                    "Select SECOND VCF variant file", 
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                )
                logger.debug(f"Selected second VCF file: {args.second_vcf}")
            except FileSelectionError as e:
                logger.error(f"Second VCF file selection failed: {e}")
                return False
    
    # Only prompt for GFF if annotation filtering is not disabled
    if not args.noannotation and not args.gff:
        logger.info("\n>>> Please select GFF annotation file <<<")
        try:
            args.gff = FileIO.select_file(
                "Select GFF annotation file", 
                [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
            )
            logger.debug(f"Selected GFF file: {args.gff}")
        except FileSelectionError as e:
            logger.error(f"GFF file selection failed: {e}")
            return False
    elif args.noannotation:
        logger.info("\nSkipping GFF annotation file selection")
    
    # Signal that all file selections are complete
    FileIO.mark_selection_complete()
    
    # Set up output directory and run the workflow
    return _setup_and_run_alignment_workflow(args, reference_file)


def _run_direct_alignment_workflow(args):
    """
    Run the direct alignment workflow (calculating alignment from two FASTA files).
    
    Args:
        args (argparse.Namespace): Command line arguments
        
    Returns:
        bool: True if the workflow completed successfully, False otherwise
    """
    logger.debug("\n>>> Using direct alignment workflow <<<")
    
    # Need all files for alignment workflow
    if not args.fasta:
        logger.info("\n>>> Please select REFERENCE FASTA sequence file <<<")
        try:
            args.fasta = FileIO.select_fasta_file("Select REFERENCE FASTA sequence file")
            reference_file = args.fasta
            logger.debug(f"Selected reference FASTA file: {args.fasta}")
        except FileSelectionError as e:
            logger.error(f"Reference FASTA file selection failed: {e}")
            return False
    else:
        reference_file = args.fasta
    
    if args.snp:
        # Only require VCF files if SNP masking is enabled
        if not args.vcf:
            logger.info("\n>>> Please select REFERENCE VCF variant file <<<")
            try:
                args.vcf = FileIO.select_file(
                    "Select REFERENCE VCF variant file", 
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                )
                logger.debug(f"Selected reference VCF file: {args.vcf}")
            except FileSelectionError as e:
                logger.error(f"Reference VCF file selection failed: {e}")
                return False
    
    if not args.second_fasta:
        logger.info("\n>>> Please select SECOND FASTA sequence file <<<")
        try:
            args.second_fasta = FileIO.select_fasta_file("Select SECOND FASTA sequence file")
            logger.debug(f"Selected second FASTA file: {args.second_fasta}")
        except FileSelectionError as e:
            logger.error(f"Second FASTA file selection failed: {e}")
            return False
    
    if args.snp:
        # Only require VCF files if SNP masking is enabled
        if not args.second_vcf:
            logger.info("\n>>> Please select SECOND VCF variant file <<<")
            try:
                args.second_vcf = FileIO.select_file(
                    "Select SECOND VCF variant file", 
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                )
                logger.debug(f"Selected second VCF file: {args.second_vcf}")
            except FileSelectionError as e:
                logger.error(f"Second VCF file selection failed: {e}")
                return False
    
    # Only prompt for GFF if annotation filtering is not disabled
    if not args.noannotation and not args.gff:
        logger.info("\n>>> Please select GFF annotation file <<<")
        try:
            args.gff = FileIO.select_file(
                "Select GFF annotation file", 
                [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
            )
            logger.debug(f"Selected GFF file: {args.gff}")
        except FileSelectionError as e:
            logger.error(f"GFF file selection failed: {e}")
            return False
    elif args.noannotation:
        logger.info("\nSkipping GFF annotation file selection (--noannotation specified)")
    
    # Signal that all file selections are complete
    FileIO.mark_selection_complete()
    
    # Set up output directory and run the workflow
    return _setup_and_run_alignment_workflow(args, reference_file)


def _setup_and_run_alignment_workflow(args, reference_file):
    """
    Set up output directory and run the alignment workflow.
    
    Args:
        args (argparse.Namespace): Command line arguments
        reference_file (str): Path to the reference file
        
    Returns:
        bool: True if the workflow completed successfully, False otherwise
    """
    # Set up output directory
    if args.output:
        output_dir = args.output
    elif reference_file:
        # Get the directory of the original input file
        if hasattr(args, 'maf') and args.maf:
            # For MAF files, use the parent directory of the MAF file's directory
            maf_dir = os.path.dirname(os.path.abspath(args.maf))
            output_dir = os.path.join(os.path.dirname(maf_dir), "Primers")
        else:
            # For FASTA files, use the directory of the reference FASTA
            input_dir = os.path.dirname(os.path.abspath(reference_file))
            output_dir = os.path.join(input_dir, "Primers")
        logger.debug(f"Setting output directory to: {output_dir}")
    else:
        # Fallback to current directory if no reference file
        output_dir = os.path.join(os.getcwd(), "Primers")
        logger.debug(f"Using current directory for output: {output_dir}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logger.debug(f"Created output directory: {output_dir}")
    
    # Use context manager for temporary directory
    with TempDirectoryManager(output_dir) as temp_dir:
        # Call the AlignmentWorkflow function to handle the alignment and masking
        logger.info("\n>>> Running Alignment Primer Design Workflow <<<")
        try:
            # Pass only needed arguments to run_alignment_workflow
            masked_sequences, coordinate_map = run_alignment_workflow(args, output_dir)
            logger.debug(f"Alignment workflow completed with {len(masked_sequences)} masked sequences")
            
            if not masked_sequences:
                logger.warning("No masked sequences were generated. Exiting.")
                return False
            
            # Load gene annotations if needed
            if not args.noannotation:
                logger.info("\nLoading gene annotations from GFF file...")
                try:
                    genes = AnnotationProcessor.load_genes_from_gff(args.gff)
                    logger.debug(f"Loaded {len(genes)} gene annotations")
                except Exception as e:
                    logger.error(f"Error loading gene annotations: {e}")
                    return False
            else:
                logger.info("\nSkipping gene annotation loading (--noannotation specified)")
                genes = None  # Set to None when annotation filtering is disabled
            
            # Before running primer design workflow, make sure we have the correct reference paths
            if hasattr(args, 'maf') and args.maf:
                # For MAF-based alignments, use the MAF file for naming
                reference_file = args.maf
                second_fasta_path = None
                logger.debug(f"Using MAF file for reference: {reference_file}")
            elif args.fasta and args.second_fasta:
                # For direct FASTA alignments, use both FASTA files for naming
                reference_file = args.fasta
                second_fasta_path = args.second_fasta
                logger.debug(f"Using reference FASTA: {reference_file}")
                logger.debug(f"Using second FASTA: {second_fasta_path}")
            else:
                # Fallback - use what we have
                reference_file = reference_file
                second_fasta_path = None
                logger.debug(f"Using fallback reference file: {reference_file}")
            
            # Use the common workflow function to handle the rest of the pipeline
            success = common.run_primer_design_workflow(
                masked_sequences=masked_sequences,
                output_dir=output_dir,
                reference_file=reference_file,
                mode='alignment',
                genes=genes,
                coordinate_map=coordinate_map,
                gff_file=args.gff,
                temp_dir=temp_dir,
                second_fasta=second_fasta_path,
                skip_annotation_filtering=args.noannotation
            )
            
            return success
            
        except Exception as e:
            logger.error(f"Error in alignment workflow: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Alignment workflow failed: {e}")