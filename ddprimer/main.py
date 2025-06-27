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

Contains functionality for:
1. Command line argument parsing and validation
2. Configuration management and display
3. Database creation and verification
4. Complete primer design workflow execution
5. Temporary directory management
6. Unified pipeline orchestration and error handling

This module serves as the main entry point for the ddPrimer pipeline,
coordinating all components to provide a complete primer design solution.
"""

import os
import sys
import argparse
import logging
import traceback
import pandas as pd
from pathlib import Path
from typing import Optional, Union
from tqdm import tqdm

# Import package modules
from .config import Config, setup_logging, display_config, display_primer3_settings, DDPrimerError, FileError, ExternalToolError, SequenceProcessingError, PrimerDesignError, FileSelectionError
from .utils import run_direct_mode

# Type alias for path inputs
PathLike = Union[str, Path]

# Set up module logger
logger = logging.getLogger(__name__)


#############################################################################
#                          Command Line Argument Parsing
#############################################################################

def parse_arguments():
    """
    Parse command line arguments with comprehensive help.
    
    Validates argument combinations and provides comprehensive help
    for all available options and modes.
    
    Returns:
        Parsed arguments namespace
        
    Raises:
        SystemExit: If arguments are invalid or conflicting
    """
    parser = argparse.ArgumentParser(
        description='ddPrimer: A pipeline for primer design and filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='ddprimer [-h] [--debug [MODULE...]] [--config [.json]] [--db [.fasta, .fna, .fa] [DB_NAME]]\n'
            '                [--cli]  [--nooligo] [--noannotation]\n'
            '                [--direct [.csv, .xlsx]] [--snp]\n'        
            '                [--fasta [.fasta, .fna, .fa]] [--gff [.gff, .gff3]] [--vcf [.vcf, .vcf.gz]]\\n'
            '                [--output <output_dir>]'
    )

    # Create argument groups for better organization
    option_group = parser.add_argument_group('options')
    input_group = parser.add_argument_group('inputs (optional)')

    # Pipeline Help
    option_group.add_argument('--debug', nargs='*', metavar='MODULE', 
                            help=   'Enable debug mode. Use without arguments for universal debug,'
                                    'or specify module names (e.g. "--debug blast_processor").')
    option_group.add_argument('--config', metavar='[.json]', nargs='?', const='DISPLAY', 
                            help=   'Configuration file path. With no arguments, shows config help mode. ')
    option_group.add_argument('--db', nargs='*', metavar=('[.fasta, .fna, .fa]', '[DB_NAME]'),
                            help=   'Create or select a BLAST database. With no arguments, shows database selection menu.'
                                    'Optionally use FASTA file path argument to create database,'
                                    'optional second argument to determine database name.')
    
    # Pipeline Options
    option_group.add_argument('--cli', action='store_true', 
                            help=   'Force CLI mode.')
    option_group.add_argument('--nooligo', action='store_true', 
                            help=   'Disable internal oligo (probe) design.')
    option_group.add_argument('--noannotation', action='store_true', 
                            help=   'Disable gene annotation filtering.')
    
    # Direct Mode
    option_group.add_argument('--direct', metavar='[.csv, .xlsx]', nargs='?', const=True, 
                            help=   'Enable target-sequence based primer design workflow using CSV/Excel input.')
    option_group.add_argument('--snp', action='store_true', 
                            help=   'Enable SNP masking in direct mode. Requires VCF and FASTA files.')

    # Input files
    input_group.add_argument('--fasta', metavar='[.fasta, .fna, .fa]', 
                            help=   'Reference genome FASTA file')
    input_group.add_argument('--vcf', metavar='[.vcf, .vcf.gz]', 
                            help=   'Variant Call Format (VCF) file with variants')
    input_group.add_argument('--gff', metavar='[.gff, .gff3]', 
                            help=   'GFF annotation file')
    input_group.add_argument('--output', metavar='<output_dir>', 
                            help=   'Output directory')
    
    args = parser.parse_args()

    # Process debug argument
    if args.debug is not None:
        if len(args.debug) == 0:
            # --debug with no arguments: universal debug
            args.debug = True
        else:
            # --debug with module names: specific debug
            args.debug = args.debug  # Keep as list of module names
    else:
        # --debug not specified
        args.debug = False

    # Validate direct mode arguments
    if args.direct:
        # Direct mode cannot use GFF files
        if args.gff:
            logger.error("Direct mode cannot use GFF files (--gff). Gene filtering is automatically disabled.")
            sys.exit(1)
        
    # SNP mode enables filtering in direct mode
    if args.snp:
        if not args.direct:
            logger.error("--snp flag can only be used with --direct mode.")
            sys.exit(1)
            
    # Process --db arguments (existing code remains the same)
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


#############################################################################
#                           Primer Design Workflow
#############################################################################

def run_primer_design_workflow(processed_sequences, output_dir, reference_file,
                              genes=None, gff_file=None, skip_annotation_filtering=False,
                              enable_internal_oligo=True, chromosome_map=None):
    """
    Unified primer design workflow using clean separation of concerns.
    
    Each processor has a single responsibility:
    - SequenceProcessor: Handle restriction sites and fragments
    - AnnotationProcessor: Handle gene filtering
    - Primer3Processor: Design primers with Primer3
    - PrimerProcessor: Create primer records (no filtering)
    - FilterProcessor: Filter primer records (no creation)
    - ThermoProcessor: Calculate thermodynamics
    - BlastProcessor: Check specificity
    
    Args:
        processed_sequences: Dictionary of processed sequences (masked/substituted)
        output_dir: Output directory path
        reference_file: Path to reference file (FASTA)
        genes: Gene annotations
        gff_file: Path to GFF file
        skip_annotation_filtering: Skip gene annotation filtering
        enable_internal_oligo: Whether to design internal oligos (probes)
        
    Returns:
        True if workflow completed successfully, False otherwise
        
    Raises:
        SequenceProcessingError: If sequence processing fails
        PrimerDesignError: If primer design fails
    """
    logger.debug("=== MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
    logger.debug(f"Processing {len(processed_sequences)} sequences")
    logger.debug(f"Skip annotation filtering: {skip_annotation_filtering}")
    logger.debug(f"Internal oligo design: {enable_internal_oligo}")
    
    from .utils import FileIO
    from .core import PrimerProcessor, BlastProcessor, SequenceProcessor, Primer3Processor, ThermoProcessor, AnnotationProcessor

    try:
        # Step 1: Cut sequences at restriction sites using workflow wrapper
        logger.debug("MAIN: Delegating restriction site processing")
        restriction_fragments = SequenceProcessor.process_restriction_sites_workflow(processed_sequences)
        
        if not restriction_fragments:
            logger.warning("No valid fragments after restriction site filtering. Exiting.")
            logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
            return False
        
        # Step 2: Filter fragments based on gene overlap using workflow wrapper
        logger.debug("MAIN: Delegating gene overlap filtering")
        filtered_fragments = AnnotationProcessor.filter_fragments_by_gene_overlap_workflow(
            restriction_fragments, 
            genes, 
            skip_annotation_filtering
        )
        
        # Check if filtering returned empty list or None
        if not filtered_fragments:  # This handles both None and empty list
            logger.warning("No fragments passed filtering. Exiting.")
            if not skip_annotation_filtering:
                logger.info("Suggestion: Try using --noannotation flag to skip gene filtering")
            logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
            return False
        
        # Step 3: Design primers with Primer3 using workflow wrapper
        # NOTE: This already creates primer records internally!
        logger.debug("MAIN: Delegating primer design and record creation")
        primer3_processor = Primer3Processor(Config, enable_internal_oligo=enable_internal_oligo)
        primer_records = primer3_processor.design_primers_workflow(filtered_fragments)
        
        if not primer_records:
            logger.warning("No primers were designed by Primer3. Exiting.")
            logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
            return False
        
        # REMOVED: Step 4 - Create primer records (already done in step 3)
        # The primer3_processor.design_primers_workflow() already returns primer records,
        # not raw Primer3 results, so we don't need to create them again.
        
        # Step 5: Filter primers using workflow wrapper
        logger.debug("MAIN: Delegating primer filtering")
        df = PrimerProcessor.filter_primers_workflow(primer_records)
        
        if df is None or len(df) == 0:
            logger.warning("No primers passed filtering criteria. Exiting.")
            logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
            return False
        
        # Step 6: Calculate thermodynamic properties using workflow wrapper
        logger.debug("MAIN: Delegating thermodynamic calculations")
        df = ThermoProcessor.calculate_thermodynamics_workflow(df)
        
        # Step 7: Run BLAST for specificity using workflow wrapper
        logger.debug("MAIN: Delegating BLAST specificity checking")
        df = BlastProcessor.run_blast_specificity_workflow(df)
        
        if df is None or len(df) == 0:
            logger.warning("No primers passed BLAST filtering. Exiting.")
            logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
            return False
        
        # Step 8: Save results to Excel file, passing the chromosome map for formatting
        logger.debug("MAIN: Delegating results saving")
        output_path = FileIO.save_results(
            df,
            output_dir,
            reference_file,
            mode='standard',
            chromosome_map=chromosome_map
        )
        
        if output_path:
            logger.info(f"\nResults saved to: {output_path}")
            logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
            return True
        else:
            logger.error("Failed to save results.")
            logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
            return False
            
    except SequenceProcessingError as e:
        error_msg = f"Sequence processing error: {str(e)}"
        logger.error(error_msg)
        logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
        return False
    except PrimerDesignError as e:
        error_msg = f"Primer design error: {str(e)}"
        logger.error(error_msg)
        logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
        return False
    except Exception as e:
        error_msg = f"Error in primer design workflow: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        logger.debug("=== END MAIN WORKFLOW: PRIMER DESIGN PIPELINE ===")
        return False


#############################################################################
#                              Standard Mode Execution
#############################################################################

def run_standard_mode(args):
    """
    Run the standard mode primer design workflow with file preparation.
    
    Modified version that includes automatic file preparation before
    proceeding with the primer design pipeline.
    """
    logger.info("=== Standard Mode Workflow ===")
    logger.debug("=== MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
    
    from .utils import FileIO
    from .core import SNPMaskingProcessor, AnnotationProcessor

    try:
        # Get input files if not provided in args
        if not args.fasta:
            logger.info("\n>>> Please select FASTA sequence file <<<")
            try:
                args.fasta = FileIO.select_file(
                    "Select FASTA sequence file",
                    [
                        ("FASTA Files", "*.fasta"), 
                        ("FASTA Files", "*.fa"), 
                        ("FASTA Files", "*.fna"), 
                        ("All Files", "*")]
                )
            except FileSelectionError as e:
                error_msg = f"FASTA file selection failed: {str(e)}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
                return False
        
        # VCF file selection
        if not args.vcf:
            logger.info("\n>>> Please select VCF variant file <<<")
            try:
                args.vcf = FileIO.select_file(
                    "Select VCF variant file", 
                    [
                        ("VCF Files", "*.vcf"),
                        ("Compressed VCF Files", "*.vcf.gz"),
                        ("All Files", "*.*")
                    ]
                )
            except FileSelectionError as e:
                error_msg = f"VCF file selection failed: {str(e)}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
                return False

            # GFF file selection
            if not args.noannotation and not args.gff:
                logger.info("\n>>> Please select GFF annotation file <<<")
                try:
                    args.gff = FileIO.select_file(
                        "Select GFF annotation file", 
                        [
                            ("GFF Files", "*.gff"),
                            ("GFF3 Files", "*.gff3"), 
                            ("Compressed GFF Files", "*.gff.gz"),
                            ("Compressed GFF3 Files", "*.gff3.gz"),
                            ("All Files", "*.*")
                        ]
                    )
                except FileSelectionError as e:
                    error_msg = f"GFF file selection failed: {str(e)}"
                    logger.error(error_msg)
                    logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
                    return False
        elif args.noannotation:
            logger.info("\nSkipping GFF annotation file selection")
            args.gff = None

        # Signal that all file selections are complete
        FileIO.mark_selection_complete()
        chromosome_map = None # Initialize map
        # File preparation step
        try:
            from .utils import FilePreparator
            
            logger.debug("MAIN: Delegating file preparation")
            prep_result = FilePreparator.prepare_pipeline_files_workflow(
                vcf_file=args.vcf,
                fasta_file=args.fasta,
                gff_file=args.gff
            )
            
            if not prep_result['success']:
                if prep_result.get('reason'):
                    logger.error(f"File preparation failed: {prep_result['reason']}")
                logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
                return False
            
            # Store the chromosome map for the final output
            chromosome_map = prep_result.get('chromosome_map')

            # Update file paths to use prepared files
            if prep_result['changes_made']:
                args.vcf = prep_result.get('vcf_file', args.vcf)
                args.fasta = prep_result.get('fasta_file', args.fasta)
                args.gff = prep_result.get('gff_file', args.gff)
            else:
                logger.debug("Original files are compatible and will be used as-is")
            
        except Exception as e:
            error_msg = f"File preparation failed: {str(e)}"
            logger.error(error_msg)
            logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
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
        
        # Extract variants from VCF file and process sequences using workflow wrapper
        logger.info("\nProcessing sequences with VCF variants...")
        logger.debug("MAIN: Delegating VCF sequence processing")
        try:
            # Load sequences from FASTA file first
            sequences = FileIO.load_fasta(args.fasta)
            
            # Process sequences with VCF using the workflow wrapper
            processed_sequences = SNPMaskingProcessor.process_sequences_with_vcf_batch(
                sequences=sequences,
                vcf_file=args.vcf,  # This is now the prepared VCF
                reference_file=args.fasta  # This might be the original or prepared FASTA
            )
            
        except Exception as e:
            error_msg = f"Error processing sequences with VCF: {str(e)}"
            logger.error(error_msg)
            logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
            raise SequenceProcessingError(error_msg) from e
        
        # Load gene annotations if needed
        if not args.noannotation:
            logger.info("\nLoading gene annotations from GFF file...")
            logger.debug("MAIN: Delegating gene annotation loading")
            try:
                genes = AnnotationProcessor.load_genes_from_gff(args.gff)
            except Exception as e:
                error_msg = f"Error loading gene annotations: {str(e)}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
                return False
        else:
            logger.info("\nSkipping gene annotation loading (--noannotation specified)")
            genes = None
            
        # Determine internal oligo setting from args
        enable_internal_oligo = not args.nooligo
        
        # Run the primer design workflow
        logger.debug("MAIN: Delegating to primer design workflow")
        success = run_primer_design_workflow(
            processed_sequences=processed_sequences,
            output_dir=output_dir,
            reference_file=args.fasta,
            genes=genes,
            gff_file=args.gff,
            skip_annotation_filtering=args.noannotation,
            enable_internal_oligo=not args.nooligo,
            chromosome_map=chromosome_map
        )
        
        logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
        return success
            
            
    except SequenceProcessingError as e:
        error_msg = f"Sequence processing error: {str(e)}"
        logger.error(error_msg)
        logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
        return False
    except Exception as e:
        error_msg = f"Error in standard mode workflow: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        logger.debug("=== END MAIN WORKFLOW: STANDARD MODE EXECUTION ===")
        return False


#############################################################################
#                          Main Pipeline Execution
#############################################################################

def run_pipeline():
    """
    Run the primer design pipeline.
    
    Main entry point for pipeline execution. Handles argument parsing,
    configuration loading, database verification, and workflow execution.
    
    Returns:
        True if the pipeline completed successfully, False otherwise
        
    Raises:
        DDPrimerError: For application-specific errors
        FileError: For file-related operations
        ExternalToolError: For external tool failures
    """
    logger.debug("=== MAIN WORKFLOW: PIPELINE EXECUTION ===")
    
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Setup logging
        log_file = setup_logging(debug=args.debug if args is not None else False)

        from .utils import FileIO, BlastDatabaseManager

        # Get the config
        Config.get_instance()
        
        # Display configuration if --config is provided with special values
        if args.config in ['DISPLAY', 'primer3', 'template']:
            if args.config == 'template':
                # Generate template configuration file
                from .config import generate_config_template
                # Use the output directory if provided
                output_dir = args.output if hasattr(args, 'output') and args.output else None
                generate_config_template(Config, output_dir=output_dir)
            else:
                if args.config == 'primer3':
                    display_primer3_settings(Config)
                else:
                    # Display configuration
                    display_config(Config)
            logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
            return True

        # Force CLI mode if specified
        if args.cli:
            FileIO.use_cli = True

        # Load custom configuration if provided
        if args.config and args.config not in ['DISPLAY', 'primer3', 'template']:
            try:
                Config.load_from_file(args.config)
            except FileNotFoundError as e:
                error_msg = f"Configuration file not found: {args.config}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                raise FileError(error_msg) from e
            except Exception as e:
                error_msg = f"Error loading configuration file {args.config}: {str(e)}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                raise FileError(error_msg) from e
            
        # Apply nooligo setting if specified
        if args.nooligo:
            logger.info("Internal oligo (probe) design is disabled")
            # Modify settings
            Config.PRIMER3_SETTINGS["PRIMER_PICK_INTERNAL_OLIGO"] = 0
        
        # Process BLAST database arguments
        if args.db is not None:
            try:
                logger.debug("MAIN: Delegating BLAST database setup")
                blast_db_manager = BlastDatabaseManager()
                
                # Initialize variables that will be used in database creation
                organism_key = None
                fasta_file = None
                db_name = None
                
                # Handle model organism / existing db selection
                if args.db_action == 'select':
                    logger.info("=== BLAST Database selection menu ===")
                    organism_key, organism_name, fasta_file = blast_db_manager.select_model_organism()
                    
                    if organism_key is None and fasta_file is None:
                        logger.info("Database selection canceled. Exiting...")
                        logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                        return False
                    
                    # Handle the case of selecting an existing database
                    if organism_key == 'existing_db':
                        # fasta_file actually contains the database path in this case
                        selected_db_path = fasta_file
                        if blast_db_manager.set_active_database(selected_db_path):
                            logger.info(f"\n=== Successfully selected BLAST database ===\n")
                            
                            # If only running database operations, exit successfully
                            if not args.fasta:
                                logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                                return True
                        else:
                            logger.error("Failed to set selected database as active")
                            logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                            return False
                    
                    # For model organism cases, set up database name
                    if organism_key is not None and organism_key != 'existing_db':
                        # For model organism, use scientific name as the database name
                        organism_data = blast_db_manager.MODEL_ORGANISMS[organism_key]
                        organism_name = organism_data["name"]
                        scientific_name = organism_name.split(' (')[0] if ' (' in organism_name else organism_name
                        db_name = scientific_name.replace(' ', '_')
                    
                    # For custom file cases, derive database name from filename
                    elif organism_key is None and fasta_file is not None:
                        # This is a custom file selection
                        db_name = os.path.splitext(os.path.basename(fasta_file))[0]
                
                # Handle direct file path
                elif args.db_action == 'file':
                    fasta_file = args.db_fasta
                    if not Path(fasta_file).exists():
                        error_msg = f"FASTA file not found: {fasta_file}"
                        logger.error(error_msg)
                        logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                        raise FileError(error_msg)
                    organism_key = None  # Not a model organism
                    db_name = args.db_name
                    
                # If we have a file to create a database from
                if fasta_file and organism_key != 'existing_db':
                    # Only show creation message for command line usage, not menu selection
                    if args.db_action == 'file':
                        logger.info("=== BLAST Database creation ===")
                    
                    # Use the output directory if provided, otherwise use default
                    output_dir = args.output if args.output else None
                    
                    # Create the database
                    db_path = blast_db_manager.create_database(
                        fasta_file,
                        db_name,  # Custom database name
                        output_dir
                    )
                    
                    # Clean up genome file if it was from a model organism
                    if organism_key is not None and organism_key != 'existing_db':
                        blast_db_manager._cleanup_genome_file(fasta_file)

                    # Check if there's already a database path set
                    if Config.DB_PATH:
                        # If there's already a non-default database, ask if we should use the new one
                        logger.info(f"Current BLAST database path: {Config.DB_PATH}")
                        use_new_db = input("\n>>> Use the newly created database instead? <<<\n[y/n]: ").strip().lower()
                        if use_new_db == "" or use_new_db.startswith("y"):
                            blast_db_manager.set_active_database(db_path)
                            logger.info(f"\n=== Successfully updated BLAST database ===\n")
                        else:
                            logger.info(f"\n=== Successfully created BLAST database ===\n")
                    else:
                        # If no database or default database, automatically use the new one
                        blast_db_manager.set_active_database(db_path)
                        logger.info(f"\n=== Successfully created and set BLAST database ===\n")
                        
                    # If only running database operations, exit successfully
                    if not args.fasta:
                        logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
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
                # Mark file selection as complete
                FileIO.mark_selection_complete()
                logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                raise DDPrimerError(error_msg) from e
                
        # Verify BLAST database
        try:
            logger.debug("MAIN: Delegating BLAST database verification")
            blast_db_manager = BlastDatabaseManager()
            if not blast_db_manager.verify_database():
                logger.warning("BLAST database verification failed. Attempting interactive setup...")
                if not blast_db_manager.setup_database_interactive():
                    error_msg = "BLAST database verification failed, and no new database created."
                    logger.error(error_msg)
                    # Mark file selection as complete
                    FileIO.mark_selection_complete()
                    logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
                    raise ExternalToolError(error_msg, tool_name="blastn")
        except Exception as e:
            error_msg = f"Error during BLAST verification: {str(e)}"
            logger.error(error_msg)
            # Mark file selection as complete
            FileIO.mark_selection_complete()
            logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
            raise ExternalToolError(error_msg, tool_name="blastn") from e
        
        # Validate direct mode arguments
        if args.direct:
            # Execute direct mode workflow
            logger.debug("MAIN: Delegating to direct mode workflow")
            from .utils.direct_mode import run_direct_mode
            success = run_direct_mode(args)
        else:
            # Execute standard mode workflow
            logger.debug("MAIN: Delegating to standard mode workflow")
            success = run_standard_mode(args)
            
        if success:
            logger.info("\n=== Pipeline execution completed successfully! ===\n")
            logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
            return True
        else:
            logger.error("\n=== Pipeline execution failed ===")
            logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
            return False
            
    except DDPrimerError as e:
        # Handle application-specific exceptions
        error_msg = f"Pipeline error: {str(e)}"
        logger.error(error_msg)
        logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
        return False
    except Exception as e:
        # Handle unexpected exceptions
        error_msg = f"Unhandled exception during pipeline execution: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        print(traceback.format_exc())
        logger.debug("=== END MAIN WORKFLOW: PIPELINE EXECUTION ===")
        return False


#############################################################################
#                              Entry Point
#############################################################################

def main():
    """
    Entry point when running the script directly.
    
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    success = run_pipeline()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()