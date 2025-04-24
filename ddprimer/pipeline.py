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

# Import package modules with corrected import paths
from .config import (Config, setup_logging, PrimerDesignError)
from .utils.file_utils import FileUtils
from .utils.blast_db_creator import BlastDBCreator
from .modes import run_alignment_mode, run_direct_mode, run_standard_mode

# Define logger
logger = logging.getLogger("ddPrimer")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='ddPrimer: A pipeline for primer design and filtering')

    parser.add_argument('--fasta', help='Input FASTA file (reference genome)')
    parser.add_argument('--vcf', help='VCF file with variants')
    parser.add_argument('--gff', help='GFF annotation file')
    parser.add_argument('--direct', nargs='?', const=True, help='CSV or Excel file with sequence name and sequence columns (shortcut mode)')
    parser.add_argument('--output', help='Output directory')
    parser.add_argument('--config', help='Configuration file')
    parser.add_argument('--cli', action='store_true', help='Force CLI mode')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('--nooligo', action='store_true', help='Disable internal oligo (probe) design')
    parser.add_argument('--snp', action='store_true', help='Enable SNP masking in sequences (requires --fasta and --vcf)')
    # BLAST database creation
    parser.add_argument('--dbfasta', help='Create a BLAST database from this FASTA file (overrides config)')
    parser.add_argument('--dbname', help='Custom name for the BLAST database (optional)')
    parser.add_argument('--dboutdir', help='Custom output directory for the BLAST database (optional)')
    
    # Alignment mode options
    alignment_group = parser.add_argument_group("Alignment Options")
    alignment_group.add_argument("--alignment", action="store_true", 
                            help="Enable alignment mode primer design workflow")
    alignment_group.add_argument("--maf-file", nargs='?', const=True, 
                            help="Pre-computed MAF alignment file (skips LastZ alignment)")
    alignment_group.add_argument("--second-fasta", 
                            help="Second species genome FASTA file")
    alignment_group.add_argument("--second-vcf", 
                            help="Variant Call Format (VCF) file for the second genome")
    alignment_group.add_argument("--min-identity", type=float, default=80, 
                            help="Minimum sequence identity for primer regions (default: 80)")
    alignment_group.add_argument("--min-length", type=int, default=20,
                            help="Minimum length of conserved regions (default: 20)")
    alignment_group.add_argument("--lastz-options", default="--format=maf",
                            help="Additional options for LastZ alignment")
    
    args = parser.parse_args()

    # Check for conflicting options
    if args.direct and args.alignment:
        parser.error("--direct cannot be used with --alignment")
    
    # Handle --snp requirements
    if args.snp and args.cli:
        if not args.fasta or not args.vcf:
            parser.error("SNP masking (--snp) requires --fasta and --vcf")
    
    # Alignment mode validation (only in CLI mode)
    if args.cli and args.alignment:
        if not (args.maf_file or args.second_fasta):
            parser.error("Alignment mode requires either --maf-file or --second-fasta")
        if not args.maf_file and not args.fasta:
            parser.error("Reference genome FASTA (--fasta) is required for alignment")
        if args.snp:
            # Only require VCF files if SNP masking is enabled
            if args.second_fasta and not args.second_vcf:
                parser.error("Second species VCF file (--second-vcf) is required for SNP masking")
    
    return args


class WorkflowFactory:
    """Factory for creating the appropriate workflow based on command line arguments."""
    
    @staticmethod
    def create_workflow(args):
        """
        Create and return the appropriate workflow based on command line arguments.
        
        Args:
            args: Command line arguments
            
        Returns:
            callable: Workflow function to execute
        """
        
        if args.alignment:
            # Alignment mode
            return run_alignment_mode
        elif args.direct:
            # Direct mode
            return run_direct_mode
        else:
            # Standard mode
            return run_standard_mode

class TempDirectory:
    """Context manager for temporary directory creation and cleanup."""
    
    def __init__(self, base_dir=None):
        self.temp_dir = None
        self.base_dir = base_dir
        
    def __enter__(self):
        """Create and return the temporary directory path."""
        import tempfile
        self.temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=self.base_dir)
        logger.debug(f"Created temporary directory: {self.temp_dir}")
        return self.temp_dir
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clean up the temporary directory."""
        import shutil
        try:
            if self.temp_dir and os.path.exists(self.temp_dir):
                logger.debug(f"Cleaning up temporary directory: {self.temp_dir}")
                shutil.rmtree(self.temp_dir)
        except Exception as e:
            logger.warning(f"Error cleaning up temporary files: {e}")


def run_pipeline():
    """Run the primer design pipeline."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Setup logging
        log_file = setup_logging(debug=args.debug if args is not None else False)
        
        logger = logging.getLogger("ddPrimer")
        logger.debug("Starting pipeline execution")
        logger.debug(f"Arguments: {args}")
        logger.debug(f"Config settings: NUM_PROCESSES={Config.NUM_PROCESSES}, BATCH_SIZE={Config.BATCH_SIZE}")
        
        # Force CLI mode if specified
        if args.cli:
            FileUtils.use_cli = True
            logger.debug("CLI mode enforced via command line argument")

        # Load custom configuration if provided
        if args.config:
            logger.debug(f"Loading custom configuration from {args.config}")
            Config.load_from_file(args.config)

        # Apply nooligo setting if specified
        if args.nooligo:
            logger.info("Internal oligo (probe) design is disabled")
            # Modify settings
            Config.PRIMER3_SETTINGS["PRIMER_PICK_INTERNAL_OLIGO"] = 0
            Config.DISABLE_INTERNAL_OLIGO = True
        
        # Process BLAST database arguments
        if args.dbfasta:
            logger.info(f"Creating BLAST database from {args.dbfasta}")
            try:
                blast_db_creator = BlastDBCreator()
                db_path = blast_db_creator.create_database(
                    args.dbfasta,
                    args.dbname,
                    args.dboutdir
                )
                Config.DB_PATH = db_path
                Config.USE_CUSTOM_DB = True
                logger.info(f"BLAST database created: {db_path}")
            except Exception as e:
                logger.error(f"Error creating BLAST database: {e}")
                logger.debug(traceback.format_exc())
                return False
        
        # Only verify BLAST database if we're not in direct mode or alignment mode with --snp flag
        # Direct mode doesn't require LastZ, alignment mode only needs it if --snp is used
        if not args.direct and not (args.alignment and not args.snp):
            # Verify BLAST database is accessible
            logger.info("Verifying BLAST database...")
            from .utils import verify_blast_database
            
            if not verify_blast_database(logger):
                logger.error("\n======================================")
                logger.error("ERROR: BLAST database verification failed!")
                logger.error("\nPossible solutions:")
                logger.error("1. Rebuild the database with the makeblastdb command:")
                logger.error("   makeblastdb -in your_genome.fasta -dbtype nucl -out your_db_name")
                logger.error("\n2. Create a BLAST database directly with ddprimer:")
                logger.error("   python -m ddprimer --dbfasta your_genome.fasta [--dbname custom_name] [--dboutdir output_dir]")
                logger.error("\n3. Set a different database path in your configuration file:")
                logger.error("   python -m ddprimer --config your_config.json")
                logger.error("\n4. For memory map errors, try closing any other BLAST processes and restart.")
                logger.error("   Some memory map errors are recoverable - use --debug for detailed information.")
                logger.error("======================================")
                return False
        
        logger.info("=== Primer Design Pipeline ===")
        
        # Use factory pattern to get the appropriate workflow
        workflow = WorkflowFactory.create_workflow(args)
        
        # Execute the workflow
        success = workflow(args)
        
        if success:
            logger.info("Pipeline execution completed successfully")
            return True
        else:
            logger.error("Pipeline execution failed")
            return False
            
    except PrimerDesignError as e:
        # Handle application-specific exceptions
        logger.error(f"Pipeline error: {e}")
        logger.debug(traceback.format_exc())
        return False
    except Exception as e:
        # Handle unexpected exceptions
        logger.error(f"Unhandled exception during pipeline execution: {e}")
        logger.debug(traceback.format_exc())
        return False


def main():
    """Entry point when running the script directly."""
    success = run_pipeline()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()