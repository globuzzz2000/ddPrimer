#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Primer Design Pipeline

A streamlined script that:
1. Prompts the user to provide a FASTA and a VCF file
2. Extracts variants from VCF and masks the FASTA file
3. Filters sequences based on restriction sites and gene overlap
4. Runs Primer3 on masked and filtered sequences
5. Filters primers and oligos 
6. Runs NUPACK for thermodynamic properties
7. BLASTs primers and oligos for specificity
8. Saves results to an Excel file
"""

import os
import sys
import argparse
import logging
import traceback
from datetime import datetime
import tempfile
import shutil

# Import package modules
from .config import Config, setup_logging, display_config, display_primer3_settings
from .config.exceptions import DDPrimerError, BlastError, FileSelectionError
from .utils.file_io import FileIO
from .utils.blast_db_creator import BlastDBCreator
from .utils.blast_verification import BlastVerification
from .modes import run_alignment_mode, run_direct_mode, run_standard_mode

# Define logger
logger = logging.getLogger("ddPrimer")


def parse_arguments():
    """
    Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
        
    Raises:
        argparse.ArgumentError: If there are conflicting or invalid arguments
    """
    parser = argparse.ArgumentParser(
        description='ddPrimer: A pipeline for primer design and filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='ddprimer [--direct [.csv, .xlsx]] [--alignment] [-h] [--debug] [--config [.json]] \n                [--cli] [--nooligo] [--snp] [--noannotation] [--lastzonly] [--createdb [.fasta, .fna, .fa [DB_NAME]]] \n                [--fasta [.fasta, .fna, .fa]] [--second-fasta [.fasta, .fna, .fa]] [--vcf [.vcf, .vcf.gz]] [--second-vcf [.vcf, .vcf.gz]] [--gff [.gff, .gff3]] [--maf [.maf]] [--output <output_dir>]'
    )

    # Create argument groups for better organization
    mode_group = parser.add_argument_group('modes')
    option_group = parser.add_argument_group('options')
    input_group = parser.add_argument_group('inputs (optional)')

    # Modes
    mode_group.add_argument('--direct', metavar='[.csv, .xlsx]', nargs='?', const=True, help='Enable target-sequence based primer design workflow')
    mode_group.add_argument("--alignment", action="store_true", help="Enable alignment-based primer design workflow")

    # Options - DO NOT add help here, argparse adds it automatically
    option_group.add_argument('--debug', action='store_true', help='Enable debug mode')
    option_group.add_argument('--output', metavar='<output_dir>', help='Output directory (for results and config templates)')
    option_group.add_argument('--config', metavar='[.json]', nargs='?', const='DISPLAY', 
                      help='Configuration file path or special mode ("all", "basic", or "template")')

    option_group.add_argument('--cli', action='store_true', help='Force CLI mode')
    option_group.add_argument('--nooligo', action='store_true', help='Disable internal oligo (probe) design')
    option_group.add_argument('--snp', action='store_true', help='For direct and alignment mode: Enable SNP masking in sequences (requires fasta and vcf file)')
    option_group.add_argument('--noannotation', action='store_true', help='For standard and alignment mode: Disable gene annotation filtering in all modes')
    option_group.add_argument('--lastzonly', action='store_true', help='Run only LastZ alignments (no primer design)')
    option_group.add_argument('--createdb', nargs='+', metavar=('[.fasta, .fna, .fa]', '[DB_NAME]'),
                      help='Create a BLAST database. Optional arguments: [.fasta, .fna, .fa] [DB_NAME]')
    
    # Input files
    input_group.add_argument('--fasta', metavar='[.fasta, .fna, .fa]', help='Reference genome FASTA file')
    input_group.add_argument('--second-fasta', metavar='[.fasta, .fna, .fa]', help='Second genome FASTA file')
    input_group.add_argument('--vcf', metavar='[.vcf, .vcf.gz]', help='Variant Call Format (VCF) file with variants')
    input_group.add_argument('--second-vcf', metavar='[.vcf, .vcf.gz]', help='VCF file for the second genome')
    input_group.add_argument('--gff', metavar='[.gff, .gff3]', help='GFF annotation file')
    input_group.add_argument('--maf', metavar='[.maf]', nargs='?', const=True, help="For alignment mode: Pre-computed MAF alignment file (skips LastZ alignment)")
    
    args = parser.parse_args()

    # Check for conflicting options
    if args.direct and args.alignment:
        parser.error("--direct cannot be used with --alignment")
    
    # If --maf is used, automatically activate alignment mode
    if args.maf:
        args.alignment = True
    
    # Handle --snp requirements
    if args.snp and args.cli:
        if not args.fasta or not args.vcf:
            parser.error("SNP masking (--snp) requires --fasta and --vcf")
    
    # Extra validation for lastzonly mode
    if args.lastzonly:
        # Force alignment mode when lastzonly is specified
        args.alignment = True
        
        # Only need reference FASTA and second FASTA for lastzonly mode
        if args.cli and not args.fasta:
            parser.error("--lastzonly requires --fasta (reference genome)")
        if args.cli and not args.second_fasta:
            parser.error("--lastzonly requires --second-fasta (second genome)")
        # These flags are incompatible with lastzonly
        if args.snp:
            parser.error("--lastzonly cannot be used with --snp")
        if args.direct:
            parser.error("--lastzonly cannot be used with --direct")
        if args.noannotation:
            parser.error("--lastzonly cannot be used with --noannotation")
    
    # Alignment mode validation (only in CLI mode)
    if args.cli and args.alignment and not args.lastzonly:
        if not (args.maf or args.second_fasta):
            parser.error("Alignment mode requires either --maf or --second-fasta")
        if not args.maf and not args.fasta:
            parser.error("Reference genome FASTA (--fasta) is required for alignment")
        if args.snp:
            # Only require VCF files if SNP masking is enabled
            if args.second_fasta and not args.second_vcf:
                parser.error("Second species VCF file (--second-vcf) is required for SNP masking")
    
    # Process createdb arguments
    if args.createdb is not None:
        # First argument is the FASTA file (if provided)
        if len(args.createdb) >= 1:
            args.createdb_fasta = args.createdb[0]
        else:
            args.createdb_fasta = True  # Will prompt for file selection
        
        # Second argument is the database name (if provided)
        if len(args.createdb) >= 2:
            args.dbname = args.createdb[1]
        else:
            args.dbname = None
            
        # Use output directory for BLAST database if provided
        args.dboutdir = args.output if args.output else None
    
    return args


class WorkflowFactory:
    """Factory for creating the appropriate workflow based on command line arguments."""
    
    @staticmethod
    def create_workflow(args):
        """
        Create and return the appropriate workflow based on command line arguments.
        
        Args:
            args (argparse.Namespace): Command line arguments
            
        Returns:
            callable: Workflow function to execute
        """
        # For both lastzonly and alignment, use alignment mode
        if args.alignment or args.lastzonly:
            # Alignment mode (also handles lastzonly)
            return run_alignment_mode
        elif args.direct:
            # Direct mode
            return run_direct_mode
        else:
            # Standard mode
            return run_standard_mode


class TempDirectoryManager:
    """Context manager for temporary directory creation and cleanup."""
    
    def __init__(self, base_dir=None):
        """
        Initialize the temporary directory manager.
        
        Args:
            base_dir (str, optional): Base directory to create the temp directory in
        """
        self.temp_dir = None
        self.base_dir = base_dir
        self.logger = logging.getLogger("ddPrimer")
        
    def __enter__(self):
        """
        Create and return the temporary directory path.
        
        Returns:
            str: Path to the temporary directory
        """
        self.temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=self.base_dir)
        self.logger.debug(f"Created temporary directory: {self.temp_dir}")
        return self.temp_dir
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Clean up the temporary directory.
        
        Args:
            exc_type: Exception type if an exception was raised
            exc_val: Exception value if an exception was raised
            exc_tb: Exception traceback if an exception was raised
        """
        try:
            if self.temp_dir and os.path.exists(self.temp_dir):
                self.logger.debug(f"Cleaning up temporary directory: {self.temp_dir}")
                shutil.rmtree(self.temp_dir)
        except Exception as e:
            self.logger.warning(f"Error cleaning up temporary files: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)


def run_pipeline():
    """
    Run the primer design pipeline.
    
    Returns:
        bool: True if the pipeline completed successfully, False otherwise
    """
    # Initialize logger first to handle any early errors
    logger = None
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Setup logging
        log_file = setup_logging(debug=args.debug if args is not None else False)
        # Hide the NUPACK shutdown banner:
        logging.getLogger("nupack.rebind.render").setLevel(logging.WARNING)
        
        logger = logging.getLogger("ddPrimer")
        logger.debug("Starting pipeline execution")
        logger.debug(f"Arguments: {args}")
        logger.debug(f"Config settings: NUM_PROCESSES={Config.NUM_PROCESSES}, BATCH_SIZE={Config.BATCH_SIZE}")
        
        # Display configuration if --config is provided with special values
        if args.config in ['DISPLAY', 'all', 'basic', 'template']:
            if args.config == 'template':
                # Generate template configuration file
                from .config.template_generator import generate_config_template
                # Use the output directory if provided
                output_dir = args.output if hasattr(args, 'output') and args.output else None
                generate_config_template(Config, output_dir=output_dir)
            else:
                # Display configuration
                display_config(Config)
                # Show Primer3 settings only for 'all' mode
                if args.config == 'all':
                    display_primer3_settings(Config)
            return True
            
        # Force CLI mode if specified
        if args.cli:
            FileIO.use_cli = True
            logger.debug("CLI mode enforced via command line argument")

        # Load custom configuration if provided
        if args.config and args.config not in ['DISPLAY', 'all', 'basic', 'template']:
            logger.debug(f"Loading custom configuration from {args.config}")
            Config.load_from_file(args.config)
            
        # Apply nooligo setting if specified
        if args.nooligo:
            logger.info("Internal oligo (probe) design is disabled")
            # Modify settings
            Config.PRIMER3_SETTINGS["PRIMER_PICK_INTERNAL_OLIGO"] = 0
            Config.DISABLE_INTERNAL_OLIGO = True
        
        # Process BLAST database arguments
        if args.createdb is not None:
            logger.info("BLAST database creation requested")
            try:
                blast_db_creator = BlastDBCreator()
                
                # Get the FASTA file path from args.createdb_fasta
                fasta_file = args.createdb_fasta if isinstance(args.createdb_fasta, str) else None
                
                # If no FASTA file was provided, prompt for one
                if not fasta_file:
                    logger.info("No FASTA file provided, prompting for file selection")
                    try:
                        fasta_file = FileIO.select_fasta_file("Select FASTA file for BLAST database creation")
                    except FileSelectionError as e:
                        logger.error(f"Failed to select FASTA file: {str(e)}")
                        return False
                        
                logger.info(f"Creating BLAST database from {fasta_file}")
                
                # Use the output directory if provided, otherwise use default
                output_dir = args.output if args.output else None
                
                db_path = blast_db_creator.create_database(
                    fasta_file,
                    args.dbname,  # This comes from the second argument to --createdb
                    output_dir
                )
                Config.DB_PATH = db_path
                Config.USE_CUSTOM_DB = True
                logger.info(f"BLAST database created: {db_path}")
            except Exception as e:
                logger.error(f"Error creating BLAST database: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                # Mark file selection as complete
                FileIO.mark_selection_complete()
                return False
        
        # Skip BLAST verification for lastzonly mode
        if not args.lastzonly:
            # Verify BLAST database
            logger.debug("Verifying BLAST database...")
            
            if not BlastVerification.verify_blast_database(logger):
                logger.error("\n======================================")
                logger.error("ERROR: BLAST database verification failed!")
                logger.error("\nPossible solutions:")
                logger.error("1. Rebuild the database with the makeblastdb command:")
                logger.error("   makeblastdb -in your_genome.fasta -dbtype nucl -out your_db_name")
                logger.error("\n2. Create a BLAST database directly with ddprimer:")
                logger.error("   ddprimer --createdb path/to/your_genome.fasta [optional_db_name] --output output_dir")
                logger.error("\n3. Set a different database path in your configuration file:")
                logger.error("   ddprimer --config your_config.json")
                logger.error("\n4. For memory map errors, try closing any other BLAST processes and restart.")
                logger.error("   Some memory map errors are recoverable - use --debug for detailed information.")
                logger.error("======================================")
                # Mark file selection as complete
                FileIO.mark_selection_complete()
                return False
        
        logger.debug("=== Primer Design Pipeline ===")
        
        # Use factory pattern to get the appropriate workflow
        workflow = WorkflowFactory.create_workflow(args)
        
        # Execute the workflow
        success = workflow(args)
        
        if success:
            logger.info("\n=== Pipeline execution completed successfully! ===")
            logger.info("")
            return True
        else:
            logger.error("Pipeline execution failed")
            return False
            
    except DDPrimerError as e:
        # Handle application-specific exceptions
        if logger:
            logger.error(f"Pipeline error: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        else:
            print(f"Pipeline error: {str(e)}")
            print(f"Run with --debug for more detailed error information")
        return False
    except Exception as e:
        # Handle unexpected exceptions
        if logger:
            logger.error(f"Unhandled exception during pipeline execution: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        else:
            print(f"Unhandled exception during pipeline execution: {str(e)}")
            print(f"Run with --debug for more detailed error information")
            print(traceback.format_exc())
        return False


def main():
    """
    Entry point when running the script directly.
    
    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    success = run_pipeline()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()