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
6. Runs ViennaRNA for thermodynamic properties
7. BLASTs primers and oligos for specificity
8. Saves results to an Excel file

Now with support for automatic model organism database creation!
"""

import sys
import argparse
import logging
import traceback
from pathlib import Path
from typing import Optional, Union

# Import package modules
from .config import Config, setup_logging, display_config, display_primer3_settings, DDPrimerError, FileError, ExternalToolError
from .utils import FileIO, BlastDBCreator, BlastVerification
from .modes import run_alignment_mode, run_direct_mode, run_standard_mode

# Type alias for path inputs
PathLike = Union[str, Path]

# Define logger
logger = logging.getLogger("ddPrimer.pipeline")


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
        usage='ddprimer [--direct [.csv, .xlsx]] [--alignment] [-h] [--debug] [--config [.json]] \n'
            '                [--cli] [--nooligo] [--snp] [--noannotation] [--lastzonly] \n'
            '                [--db [.fasta, .fna, .fa] [DB_NAME]] \n'
            '                [--fasta [.fasta, .fna, .fa]] [--second-fasta [.fasta, .fna, .fa]] [--vcf [.vcf, .vcf.gz]] \n'
            '                [--second-vcf [.vcf, .vcf.gz]] [--gff [.gff, .gff3]] [--maf [.maf]] [--output <output_dir>]'
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
    option_group.add_argument('--db', nargs='*', metavar=('[.fasta, .fna, .fa]', '[DB_NAME]'),
                    help='Create or select a BLAST database. With no arguments, shows model organism selection menu. '
                        'With one argument, creates database from the specified FASTA file. '
                        'With two arguments, creates database from the first argument and uses the second as the database name.')
        
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
                parser.error("Second VCF file (--second-vcf) is required for SNP masking")
    
    # Process --db arguments
    if args.db is not None:
        # Set the db output directory if provided
        args.dboutdir = args.output if args.output else None
        
        if len(args.db) == 0:
            # User ran "--db" without any arguments, show the selection menu
            args.db_action = 'select'
            args.db_fasta = None
            args.db_name = None
        elif len(args.db) >= 1:
            # User provided a specific file path: "--db path/to/file.fasta"
            args.db_action = 'file'
            args.db_fasta = args.db[0]
            
            # Check for optional database name
            if len(args.db) >= 2:
                args.db_name = args.db[1]
            else:
                args.db_name = None
    else:
        args.db_action = None
        args.db_fasta = None
        args.db_name = None
    
    return args


class WorkflowFactory:
    """
    Factory for creating the appropriate workflow based on command line arguments.
    
    This class provides a factory method to determine and return the correct
    workflow function based on the command line arguments provided by the user.
    
    Example:
        >>> args = parse_arguments()
        >>> workflow = WorkflowFactory.create_workflow(args)
        >>> success = workflow(args)
    """
    
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
    """
    Context manager for temporary directory creation and cleanup.
    
    Provides a safe way to create and automatically clean up temporary
    directories used during pipeline execution.
    
    Attributes:
        temp_dir: Path to the created temporary directory
        base_dir: Base directory for temporary directory creation
        logger: Logger instance for this manager
        
    Example:
        >>> with TempDirectoryManager() as temp_dir:
        ...     # Use temp_dir for operations
        ...     pass
        # temp_dir is automatically cleaned up
    """
    
    def __init__(self, base_dir: Optional[PathLike] = None):
        """
        Initialize the temporary directory manager.
        
        Args:
            base_dir: Base directory to create the temp directory in
        """
        self.temp_dir = None
        self.base_dir = str(base_dir) if base_dir else None
        self.logger = logging.getLogger("ddPrimer.temp_manager")
        
    def __enter__(self):
        """
        Create and return the temporary directory path.
        
        Returns:
            str: Path to the temporary directory
            
        Raises:
            FileError: If temporary directory creation fails
        """
        try:
            import tempfile
            self.temp_dir = tempfile.mkdtemp(
                prefix=Config.BLAST_TEMP_FILE_PREFIX,
                dir=self.base_dir
            )
            self.logger.debug(f"Created temporary directory: {self.temp_dir}")
            return self.temp_dir
        except OSError as e:
            error_msg = f"Failed to create temporary directory: {str(e)}"
            self.logger.error(error_msg)
            raise FileError(error_msg) from e
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Clean up the temporary directory.
        
        Args:
            exc_type: Exception type if an exception was raised
            exc_val: Exception value if an exception was raised
            exc_tb: Exception traceback if an exception was raised
        """
        if self.temp_dir and Path(self.temp_dir).exists():
            try:
                import shutil
                self.logger.debug(f"Cleaning up temporary directory: {self.temp_dir}")
                shutil.rmtree(self.temp_dir)
            except OSError as e:
                self.logger.warning(f"Error cleaning up temporary files: {str(e)}")
                self.logger.debug(f"Error details: {str(e)}", exc_info=True)


def run_pipeline():
    """
    Run the primer design pipeline.
    
    Main entry point for pipeline execution. Handles argument parsing,
    configuration loading, database verification, and workflow execution.
    
    Returns:
        bool: True if the pipeline completed successfully, False otherwise
        
    Raises:
        DDPrimerError: For application-specific errors
        FileError: For file-related operations
        ExternalToolError: For external tool failures
    """
    # Initialize logger first to handle any early errors
    logger = None

    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Setup logging
        log_file = setup_logging(debug=args.debug if args is not None else False)
        
        logger = logging.getLogger("ddPrimer.pipeline")
        logger.debug("Initializing Config settings")
        Config.get_instance()
        logger.debug(f"Initial DB_PATH: {Config.DB_PATH}")
        logger.debug("Starting pipeline execution")
        logger.debug(f"Arguments: {args}")
        logger.debug(f"Config settings: NUM_PROCESSES={Config.NUM_PROCESSES}, BATCH_SIZE={Config.BATCH_SIZE}")
        
        # Display configuration if --config is provided with special values
        if args.config in ['DISPLAY', 'all', 'basic', 'template']:
            if args.config == 'template':
                # Generate template configuration file
                from .config import generate_config_template
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
            try:
                Config.load_from_file(args.config)
            except FileNotFoundError as e:
                error_msg = f"Configuration file not found: {args.config}"
                logger.error(error_msg)
                raise FileError(error_msg) from e
            except Exception as e:
                error_msg = f"Error loading configuration file {args.config}: {str(e)}"
                logger.error(error_msg)
                raise FileError(error_msg) from e
            
        # Apply nooligo setting if specified
        if args.nooligo:
            logger.info("Internal oligo (probe) design is disabled")
            # Modify settings
            Config.PRIMER3_SETTINGS["PRIMER_PICK_INTERNAL_OLIGO"] = 0
            Config.DISABLE_INTERNAL_OLIGO = True
        
        # Process BLAST database arguments - simplified handling
        if args.db is not None:
            logger.info("BLAST database operation requested")
            try:
                blast_db_creator = BlastDBCreator()
                
                # Handle model organism / existing db selection
                if args.db_action == 'select':
                    logger.info("Database selection requested")
                    from .utils import ModelOrganismManager as MOManager
                    organism_key, organism_name, fasta_file = MOManager.select_model_organism(logger)
                    
                    if organism_key is None and fasta_file is None:
                        logger.info("Database selection canceled. Exiting...")
                        return False
                    
                    # Handle the case of selecting an existing database
                    if organism_key == 'existing_db':
                        # fasta_file actually contains the database path in this case
                        selected_db_path = fasta_file
                        Config.DB_PATH = selected_db_path
                        Config.save_database_config(selected_db_path)
                        Config.USE_CUSTOM_DB = True
                        logger.info(f"Now using BLAST database: {selected_db_path}")
                        
                        # If only running database operations, exit successfully
                        if not args.fasta:
                            return True
                        else:
                            # Continue with other operations using the selected database
                            pass
                    
                    # For model organism or custom file cases, continue with database creation
                    if organism_key is not None and args.db_name is None:
                        # For model organism, use scientific name as the database name
                        organism_name = MOManager.MODEL_ORGANISMS[organism_key]["name"]
                        scientific_name = organism_name.split(' (')[0] if ' (' in organism_name else organism_name
                        db_name = scientific_name.replace(' ', '_')
                        logger.info(f"Using database name: {db_name}")
                    else:
                        db_name = args.db_name
                
                # Handle direct file path
                elif args.db_action == 'file':
                    fasta_file = args.db_fasta
                    if not Path(fasta_file).exists():
                        error_msg = f"FASTA file not found: {fasta_file}"
                        logger.error(error_msg)
                        raise FileError(error_msg)
                    organism_key = None  # Not a model organism
                    db_name = args.db_name
                    
                # If we have a file to create a database from
                if fasta_file and organism_key != 'existing_db':
                    
                    # Use the output directory if provided, otherwise use default
                    output_dir = args.output if args.output else None
                    
                    # Create the database
                    db_path = blast_db_creator.create_database(
                        fasta_file,
                        db_name,  # Custom database name
                        output_dir
                    )
                    
                    # Clean up genome file if it was from a model organism
                    if organism_key is not None and organism_key != 'existing_db':
                        from .utils import ModelOrganismManager
                        ModelOrganismManager.cleanup_genome_file(fasta_file, logger)

                    # Check if there's already a database path set
                    if Config.DB_PATH and Config.DB_PATH != "/Library/Application Support/Blast_DBs/Tair DB/TAIR10":
                        # If there's already a non-default database, ask if we should use the new one
                        logger.info(f"Current BLAST database path: {Config.DB_PATH}")
                        use_new_db = input("Use the newly created database instead? [Y/n]: ").strip().lower()
                        if use_new_db == "" or use_new_db.startswith("y"):
                            Config.DB_PATH = db_path
                            Config.save_database_config(db_path)
                            Config.USE_CUSTOM_DB = True
                            logger.info(f"Now using new BLAST database: {db_path}")
                        else:
                            logger.info(f"Keeping current BLAST database: {Config.DB_PATH}")
                    else:
                        # If no database or default database, automatically use the new one
                        Config.DB_PATH = db_path
                        Config.save_database_config(db_path)
                        Config.USE_CUSTOM_DB = True
                        logger.info(f"BLAST database created and set as active: {db_path}")
                        
                    # If only running database operations, exit successfully
                    if not args.fasta:
                        return True
                        
            except FileError:
                # Re-raise FileError without modification
                raise
            except ExternalToolError:
                # Re-raise ExternalToolError without modification
                raise
            except Exception as e:
                error_msg = f"Unexpected error during database operation: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                # Mark file selection as complete
                FileIO.mark_selection_complete()
                raise DDPrimerError(error_msg) from e
                
        # Skip BLAST verification for lastzonly mode
        if not args.lastzonly:
            # Verify BLAST database
            logger.debug("Verifying BLAST database...")
            
            try:
                if not BlastVerification.verify_blast_database(logger):
                    error_msg = "BLAST database verification failed, and no new database created."
                    logger.error(error_msg)
                    # Mark file selection as complete
                    FileIO.mark_selection_complete()
                    raise ExternalToolError(error_msg, tool_name="blastn")
            except Exception as e:
                error_msg = f"Error during BLAST verification: {str(e)}"
                logger.error(error_msg)
                # Mark file selection as complete
                FileIO.mark_selection_complete()
                raise ExternalToolError(error_msg, tool_name="blastn") from e
        
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