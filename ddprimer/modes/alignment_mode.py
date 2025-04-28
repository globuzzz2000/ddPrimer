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
from ..helpers import AlignmentWorkflow
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
        # For LastZ-only mode, we just need reference and second FASTA files
        if args.lastzonly:
            logger.info("=== LastZ-Only Mode ===")
            
            # Only get required input files
            if not args.fasta:
                logger.info("\n>>> Please select REFERENCE species FASTA file <<<")
                try:
                    args.fasta = FileUtils.get_file(
                        "Select REFERENCE species FASTA file", 
                        [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting reference FASTA file: {e}")
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
            
            # Import lastz_runner directly
            from ..alignment.lastz_runner import LastZRunner
            
            try:
                # Create LastZ runner instance
                runner = LastZRunner()
                
                # Run LastZ alignment
                output_file = runner.run_parallel_alignment(
                    args.fasta,
                    args.second_fasta,
                    output_dir,
                    args.lastz_options,
                    processes=Config.NUM_PROCESSES,
                    keep_temp=args.debug  # Keep temp files in debug mode
                )
                
                logger.info(f"LastZ alignment completed successfully!")
                logger.info(f"Alignment file saved to: {output_file}")
                return True
                
            except Exception as e:
                logger.error(f"Error running LastZ alignment: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
                
        # Standard alignment mode - continue with existing code
        # Check if SNP masking is enabled
        if args.snp:
            logger.info("\n>>> SNP masking is enabled <<<")
        else:
            logger.info("\n>>> SNP masking is disabled <<<")
        
        # Check and get required input files
        # Modified check for maf_file to handle the case when the flag is used without a value
        if hasattr(args, 'maf_file') and args.maf_file is not None:
            # If --maf-file was used but no value provided, prompt for file selection
            if args.maf_file is True or args.maf_file == '':
                logger.info("\n>>> Please select MAF alignment file <<<")
                try:
                    args.maf_file = FileUtils.get_file(
                        "Select MAF alignment file", 
                        [("MAF Files", "*.maf"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting MAF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            
            # Using pre-computed MAF, only need VCF and GFF files
            logger.info("\n>>> Using pre-computed MAF file <<<")
            reference_file = args.maf_file
            
            if args.snp:
                # Only require VCF files if SNP masking is enabled
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
            
            # Only prompt for GFF if annotation filtering is not disabled
            if not args.noannotation and not args.gff:
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
            elif args.noannotation:
                logger.info("\n>>> Skipping GFF annotation file selection (--noannotation specified) <<<")
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
            
            if args.snp:
                # Only require VCF files if SNP masking is enabled
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
            
            if args.snp:
                # Only require VCF files if SNP masking is enabled
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
            
            # Only prompt for GFF if annotation filtering is not disabled
            if not args.noannotation and not args.gff:
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
            elif args.noannotation:
                logger.info("\n>>> Skipping GFF annotation file selection (--noannotation specified) <<<")
        
        # Set up output directory
        if args.output:
            output_dir = args.output
        elif reference_file:
            # Get the directory of the original input file (test_data directory)
            if args.maf_file:
                # For MAF files, use the directory the MAF file is in
                input_dir = os.path.dirname(os.path.abspath(args.maf_file))
            else:
                # For FASTA files, use the directory of the reference FASTA
                input_dir = os.path.dirname(os.path.abspath(reference_file))
            
            # Create Primers directory in the same location as input files
            output_dir = os.path.join(input_dir, "Primers")
            logger.debug(f"Setting output directory to: {output_dir}")
        else:
            # Fallback to current directory if no reference file
            output_dir = os.path.join(os.getcwd(), "Primers")
            logger.debug(f"Using current directory for output: {output_dir}")

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
                # Pass snp flag to AlignmentWorkflow
                masked_sequences, coordinate_map = AlignmentWorkflow(args, output_dir, logger)
                logger.debug(f"Alignment workflow completed with {len(masked_sequences)} masked sequences")
            except Exception as e:
                logger.error(f"Error in alignment workflow: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return False
            
            if not masked_sequences:
                logger.warning("No masked sequences were generated. Exiting.")
                return False
            
            # Load gene annotations if needed
            if not args.noannotation:
                logger.info("\nLoading gene annotations from GFF file...")
                try:
                    genes = AnnotationProcessor.load_genes_from_gff(args.gff)
                    logger.debug(f"Gene annotations loaded successfully from {args.gff}")
                    logger.info(f"Loaded {len(genes)} gene annotations")
                except Exception as e:
                    logger.error(f"Error loading gene annotations: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return False
            else:
                logger.info("\nSkipping gene annotation loading (--noannotation specified)")
                genes = None  # Set to None when annotation filtering is disabled
            
            # Before running primer design workflow, make sure we have the correct reference paths
            if args.maf_file:
                # For MAF-based alignments, use the MAF file for naming
                reference_file = args.maf_file
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