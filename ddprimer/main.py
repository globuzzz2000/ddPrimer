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
from .utils import FileIO, BlastDatabaseManager
from .core import PrimerProcessor, BlastProcessor, SequenceProcessor, Primer3Processor, ViennaRNAProcessor, SNPMaskingProcessor, AnnotationProcessor

# Type alias for path inputs
PathLike = Union[str, Path]

# Set up module logger
logger = logging.getLogger(__name__)


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
            '                [--fasta [.fasta, .fna, .fa]] [--vcf [.vcf, .vcf.gz]] \n'
            '                [--gff [.gff, .gff3]] [--output <output_dir>]'
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
    
    # Direct Mode (placeholder for future implementation)
    option_group.add_argument('--direct', metavar='[.csv, .xlsx]', nargs='?', const=True, 
                            help=   'Enable target-sequence based primer design workflow (not yet implemented).')
    option_group.add_argument('--snp', action='store_true', 
                            help=   'For direct mode: Enable SNP masking in sequences (not yet implemented).')

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

    # Check for direct mode and provide informative message
    if args.direct or args.snp:
        logger.error("Direct mode (--direct) and SNP masking (--snp) are not yet implemented.")
        logger.error("Currently only standard mode is available.")
        sys.exit(1)
    
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


class TempDirectoryManager:
    """
    Context manager for temporary directory creation and cleanup.
    
    Provides a safe way to create and automatically clean up temporary
    directories used during pipeline execution.
    
    Attributes:
        temp_dir: Path to the created temporary directory
        base_dir: Base directory for temporary directory creation
        
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
        
    def __enter__(self):
        """
        Create and return the temporary directory path.
        
        Returns:
            Path to the temporary directory
            
        Raises:
            FileError: If temporary directory creation fails
        """
        try:
            import tempfile
            self.temp_dir = tempfile.mkdtemp(
                prefix="ddprimer_temp_",
                dir=self.base_dir
            )
            return self.temp_dir
        except OSError as e:
            error_msg = f"Failed to create temporary directory: {str(e)}"
            logger.error(error_msg)
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
                shutil.rmtree(self.temp_dir)
            except OSError as e:
                logger.warning(f"Error cleaning up temporary files: {str(e)}")


def prepare_files_for_pipeline(vcf_file, reference_file, gff_file=None, interactive=True):
    """
    Standalone file preparation function.
    
    Returns:
        Dictionary with prepared file paths: {'vcf': path, 'fasta': path, 'gff': path}
    """
    from .utils import FilePreparator
    
    logger.info("Preparing files for ddPrimer pipeline...")
    
    with FilePreparator() as preparator:
        result = preparator.prepare_files(
            vcf_file=vcf_file,
            fasta_file=reference_file,
            gff_file=gff_file,
            interactive=interactive
        )

        return result['prepared_files']


def process_sequences_with_vcf(sequences, vcf_file, reference_file):
    """
    Process sequences using VCF normalization approach with intelligent chromosome mapping.
    
    Args:
        sequences: Dictionary of {seq_id: sequence_string}
        vcf_file: Path to VCF file
        reference_file: Path to reference FASTA file
        
    Returns:
        Dictionary of processed sequences with variants applied
        
    Raises:
        SequenceProcessingError: If VCF processing fails
    """
    
    try:
        # Validate VCF dependencies
        Config.validate_vcf_dependencies()
        
        # Initialize SNP processor
        snp_processor = SNPMaskingProcessor(reference_file)
        
        # Get processing settings
        vcf_settings = Config.get_vcf_processing_settings()
        
        processed_sequences = {}
        
        for seq_id, sequence in sequences.items():
            try:
                # Apply VCF variants to sequence
                modified_sequence = snp_processor.process_sequence_with_vcf(
                    sequence=sequence,
                    vcf_path=vcf_file,
                    chromosome=seq_id,
                    **vcf_settings
                )
                
                processed_sequences[seq_id] = modified_sequence
                
            except Exception as e:
                logger.error(f"Error processing sequence {seq_id}: {e}")
                # Keep original sequence if processing fails
                processed_sequences[seq_id] = sequence
                
        logger.info(f"Successfully processed {len(processed_sequences)} sequences with VCF")
        return processed_sequences
        
    except Exception as e:
        error_msg = f"VCF processing failed: {str(e)}"
        logger.error(error_msg)
        raise SequenceProcessingError(error_msg) from e


def process_restriction_sites(processed_sequences):
    """
    Cut sequences at restriction sites.
    
    Processes sequences (which may already be masked/substituted) to identify 
    and cut at restriction sites, creating fragments suitable for primer design.
    
    Args:
        processed_sequences: Dictionary of sequences (masked and/or with fixed SNPs substituted)
        
    Returns:
        List of restriction fragments
        
    Raises:
        SequenceProcessingError: If there's an error in restriction site processing
    """
    logger.info("\nFiltering sequences by restriction sites...")
    
    try:
        # Use the standard restriction site method
        restriction_fragments = SequenceProcessor.cut_at_restriction_sites(processed_sequences)
        
        logger.info(f"Generated {len(restriction_fragments)} fragments after restriction site cutting")
        
        return restriction_fragments
    except Exception as e:
        error_msg = f"Error in restriction site filtering: {str(e)}"
        logger.error(error_msg)
        raise SequenceProcessingError(error_msg) from e


def filter_fragments_by_gene_overlap(restriction_fragments, genes, skip_annotation_filtering=False):
    """
    Filter fragments based on gene overlap.
    
    Applies gene overlap filtering to restriction fragments, extracting
    gene-overlapping regions for primer design.
    
    Args:
        restriction_fragments: List of restriction fragments
        genes: Gene annotations
        skip_annotation_filtering: Skip gene annotation filtering
        
    Returns:
        List of filtered fragments
        
    Raises:
        SequenceProcessingError: If there's an error in fragment filtering
    """
    if skip_annotation_filtering:
        # Create simplified fragments without location data
        filtered_fragments = []
        for fragment in restriction_fragments:
            simplified_fragment = {
                "id": fragment["id"],
                "sequence": fragment["sequence"],
                "chr": fragment.get("chr", ""),
                "start": fragment.get("start", 1),
                "end": fragment.get("end", len(fragment["sequence"])),
                "Gene": fragment["id"].split("_")[-1]
            }
            filtered_fragments.append(simplified_fragment)
        return filtered_fragments
    else:
        # Extract ALL gene overlaps for standard mode
        if not genes:
            logger.error("Gene annotations not provided for standard mode.")
            logger.info("Hint: Use --noannotation flag to skip gene filtering, or provide a GFF file")
            return []
                
        logger.info("Extracting gene-overlapping regions...")
        try:
            # Use the function that extracts ALL gene overlaps
            filtered_fragments = AnnotationProcessor.filter_by_gene_overlap(restriction_fragments, genes)
            return filtered_fragments
        except Exception as e:
            error_msg = f"Error in gene overlap extraction: {str(e)}"
            logger.error(error_msg)
            raise SequenceProcessingError(error_msg) from e


def prepare_primer3_inputs(fragments):
    """
    Prepare input blocks for Primer3 from sequence fragments.
    
    Creates Primer3-compatible input dictionaries from fragment data,
    storing fragment metadata for later result processing.
    No target regions - let Primer3 freely choose amplicon locations based on product size range.
    
    Args:
        fragments: List of sequence fragments
        
    Returns:
        Tuple of (primer3_inputs, fragment_info) - Lists of Primer3 input dicts and fragment info dict
    """
    primer3_inputs = []
    fragment_info = {}  # Dictionary to store fragment information
    
    for fragment in fragments:
        # Store location information
        fragment_info[fragment["id"]] = {
            "chr": fragment.get("chr", ""),
            "start": fragment.get("start", 1),
            "end": fragment.get("end", len(fragment["sequence"])),
            "gene": fragment.get("Gene", fragment["id"].split("_")[-1])
        }
        
        # Create primer3 input
        primer3_input = {
            "SEQUENCE_ID": fragment["id"],
            "SEQUENCE_TEMPLATE": fragment["sequence"],
        }

        primer3_inputs.append(primer3_input)
    
    return primer3_inputs, fragment_info


def design_primers_with_primer3(fragments):
    """
    Run Primer3 design on the provided fragments.
    
    Prepares input for Primer3, executes primer design in parallel,
    and parses results into a standardized format.
    
    Args:
        fragments: List of sequence fragments
        
    Returns:
        List of primer design results
        
    Raises:
        PrimerDesignError: If there's an error in primer design
    """
    logger.info("\nDesigning primers with Primer3...")
    
    # Initialize Primer3Processor
    primer3_processor = Primer3Processor(Config)
    
    # Prepare input blocks for Primer3
    primer3_inputs, fragment_info = prepare_primer3_inputs(fragments)
    
    if not primer3_inputs:
        logger.warning("No valid fragments for primer design.")
        return None
        
    # Run Primer3
    try:
        # Use parallel processing with progress bar
        primer3_output = primer3_processor.run_primer3_batch_parallel(primer3_inputs)
        
        # Parse the results
        primer_results = primer3_processor.parse_primer3_batch(primer3_output, fragment_info)
        
        return primer_results
        
    except Exception as e:
        error_msg = f"Error running Primer3: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e


def filter_primers(primer_results):
    """
    Filter primers based on various criteria.
    
    Applies sequential filtering based on penalty scores, repeat sequences,
    and GC content to ensure high-quality primer selection.
    
    Args:
        primer_results: List of primer results from Primer3
        
    Returns:
        Filtered primer DataFrame, or None if no primers pass filtering
        
    Raises:
        PrimerDesignError: If there's an error in primer filtering
    """
    # Convert to DataFrame
    df = pd.DataFrame(primer_results)
    initial_count = len(df)

    # Filter by penalty
    try:
        df = PrimerProcessor.filter_by_penalty(df)
    except Exception as e:
        error_msg = f"Error in penalty filtering: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    # Filter by repeats
    try:
        df = PrimerProcessor.filter_by_repeats(df)
    except Exception as e:
        error_msg = f"Error in repeat filtering: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    # Filter by GC content
    try:
        df = PrimerProcessor.filter_by_gc_content(df)
    except Exception as e:
        error_msg = f"Error in GC content filtering: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    logger.info(f"After filtering: {len(df)}/{initial_count} primers")
    
    # Process internal oligos (reverse complement if needed)
    try:
        df = PrimerProcessor.process_internal_oligos(df)
    except Exception as e:
        error_msg = f"Error in internal oligo processing: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    if len(df) == 0:
        logger.warning("No primers passed filtering.")
        return None
        
    return df


def calculate_thermodynamics(df):
    """
    Calculate thermodynamic properties using ViennaRNA.
    
    Computes minimum free energy (ΔG) for primers, probes, and amplicons
    using ViennaRNA with DNA-specific parameters.
    
    Args:
        df: DataFrame with primer information
        
    Returns:
        DataFrame with added thermodynamic properties
        
    Raises:
        PrimerDesignError: If there's an error in thermodynamic calculations
    """
        
    # Add empty columns to maintain DataFrame structure
    df["Primer F dG"] = None
    df["Primer R dG"] = None
    if "Probe" in df.columns:
        df["Probe dG"] = None
    df["Amplicon dG"] = None
        
    logger.info("\nCalculating thermodynamic properties with ViennaRNA...")
    
    # Calculate deltaG for forward primers
    try:
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing forward primers with ViennaRNA")
            df["Primer F dG"] = df["Primer F"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        else:
            df["Primer F dG"] = df["Primer F"].apply(ViennaRNAProcessor.calc_deltaG)
    except Exception as e:
        error_msg = f"Error calculating forward primer deltaG: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    # Calculate deltaG for reverse primers
    try:
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing reverse primers with ViennaRNA")
            df["Primer R dG"] = df["Primer R"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        else:
            df["Primer R dG"] = df["Primer R"].apply(ViennaRNAProcessor.calc_deltaG)
    except Exception as e:
        error_msg = f"Error calculating reverse primer deltaG: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    # Calculate deltaG for probes if present
    if "Probe" in df.columns:
        try:
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc="Processing probes with ViennaRNA")
                df["Probe dG"] = df["Probe"].progress_apply(lambda x: 
                                            ViennaRNAProcessor.calc_deltaG(x) 
                                            if pd.notnull(x) and x else None)
            else:
                df["Probe dG"] = df["Probe"].apply(lambda x: 
                                               ViennaRNAProcessor.calc_deltaG(x) 
                                               if pd.notnull(x) and x else None)
        except Exception as e:
            error_msg = f"Error calculating probe deltaG: {str(e)}"
            logger.error(error_msg)
            raise PrimerDesignError(error_msg) from e
    
    # Calculate deltaG for amplicons
    try:
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc="Processing amplicons with ViennaRNA")
            df["Amplicon dG"] = df["Amplicon"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        else:
            df["Amplicon dG"] = df["Amplicon"].apply(ViennaRNAProcessor.calc_deltaG)
    except Exception as e:
        error_msg = f"Error calculating amplicon deltaG: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
        
    return df


def run_blast_specificity(df):
    """
    Run BLAST for primer specificity checking.
    
    Executes BLAST analysis for all primers and probes to assess specificity,
    then filters results based on BLAST e-value thresholds.
    
    Args:
        df: DataFrame with primer information
        
    Returns:
        DataFrame with added BLAST results and filtered for specificity,
        or None if no primers pass BLAST filtering
        
    Raises:
        PrimerDesignError: If there's an error in BLAST execution or filtering
    """
    logger.info("\nRunning BLAST for specificity checking...")
    
    # Run BLAST for forward primers
    try:
        blast_results_f = []
        primers_f = df["Primer F"].tolist()
        if Config.SHOW_PROGRESS:
            primers_f_iter = tqdm(primers_f, total=len(primers_f), desc="BLASTing forward primers")
        else:
            primers_f_iter = primers_f
            
        for primer_f in primers_f_iter:
            blast1, blast2 = BlastProcessor.blast_short_seq(primer_f)
            blast_results_f.append((blast1, blast2))
        
        df["Primer F BLAST1"], df["Primer F BLAST2"] = zip(*blast_results_f)
    except Exception as e:
        error_msg = f"Error in forward primer BLAST: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    # Run BLAST for reverse primers
    try:
        blast_results_r = []
        primers_r = df["Primer R"].tolist()
        if Config.SHOW_PROGRESS:
            primers_r_iter = tqdm(primers_r, total=len(primers_r), desc="BLASTing reverse primers")
        else:
            primers_r_iter = primers_r
            
        for primer_r in primers_r_iter:
            blast1, blast2 = BlastProcessor.blast_short_seq(primer_r)
            blast_results_r.append((blast1, blast2))
        
        df["Primer R BLAST1"], df["Primer R BLAST2"] = zip(*blast_results_r)
    except Exception as e:
        error_msg = f"Error in reverse primer BLAST: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    # Run BLAST for probes if present
    if "Probe" in df.columns:
        try:
            blast_results_p = []
            probes = df["Probe"].tolist()
            if Config.SHOW_PROGRESS:
                probes_iter = tqdm(probes, total=len(probes), desc="BLASTing probes")
            else:
                probes_iter = probes
                
            for probe in probes_iter:
                if pd.notnull(probe) and probe:
                    blast1, blast2 = BlastProcessor.blast_short_seq(probe)
                else:
                    blast1, blast2 = None, None
                blast_results_p.append((blast1, blast2))
            
            df["Probe BLAST1"], df["Probe BLAST2"] = zip(*blast_results_p)
        except Exception as e:
            error_msg = f"Error in probe BLAST: {str(e)}"
            logger.error(error_msg)
            raise PrimerDesignError(error_msg) from e
    
    # Filter by BLAST specificity
    try:
        initial_count = len(df)
        df = PrimerProcessor.filter_by_blast(df)
        logger.info(f"After BLAST filtering: {len(df)}/{initial_count} primers")
    except Exception as e:
        error_msg = f"Error in BLAST filtering: {str(e)}"
        logger.error(error_msg)
        raise PrimerDesignError(error_msg) from e
    
    if len(df) == 0:
        logger.warning("No primers passed BLAST filtering.")
        return None
        
    return df


def run_primer_design_workflow(processed_sequences, output_dir, reference_file, 
                              genes=None, gff_file=None, skip_annotation_filtering=False):
    """
    Unified primer design workflow.
    
    Executes the complete primer design pipeline including restriction site processing,
    fragment filtering, primer design, thermodynamic calculations, and BLAST analysis.
    
    Args:
        processed_sequences: Dictionary of processed sequences (masked/substituted)
        output_dir: Output directory path
        reference_file: Path to reference file (FASTA)
        genes: Gene annotations
        gff_file: Path to GFF file
        skip_annotation_filtering: Skip gene annotation filtering
        
    Returns:
        True if workflow completed successfully, False otherwise
        
    Raises:
        SequenceProcessingError: If sequence processing fails
        PrimerDesignError: If primer design fails
    """
    try:
        # Step 1: Cut sequences at restriction sites
        restriction_fragments = process_restriction_sites(processed_sequences)
        logger.debug(f"CHECKPOINT 1: {len(restriction_fragments)} restriction fragments")
        
        if not restriction_fragments:
            logger.warning("No valid fragments after restriction site filtering. Exiting.")
            return False
        
        # Step 2: Filter fragments based on gene overlap
        filtered_fragments = filter_fragments_by_gene_overlap(
            restriction_fragments, 
            genes, 
            skip_annotation_filtering
        )
        
        # Check if filtering returned empty list or None
        if not filtered_fragments:  # This handles both None and empty list
            logger.warning("No fragments passed filtering. Exiting.")
            if not skip_annotation_filtering:
                logger.info("Suggestion: Try using --noannotation flag to skip gene filtering")
            return False
        
        logger.debug(f"CHECKPOINT 2: {len(filtered_fragments)} filtered fragments")
        
        # Step 3: Design primers with Primer3
        primer_results = design_primers_with_primer3(filtered_fragments)
        logger.debug(f"CHECKPOINT 3: {len(primer_results)} primer results from Primer3")
        
        if not primer_results:
            logger.warning("No primers were designed by Primer3. Exiting.")
            return False
        
        # Step 4: Filter primers
        df = filter_primers(primer_results)
        logger.debug(f"CHECKPOINT 4: {len(df)} primers after all filtering")
        
        if df is None or len(df) == 0:
            logger.warning("No primers passed filtering criteria. Exiting.")
            return False
        
        # Step 5: Calculate thermodynamic properties
        df = calculate_thermodynamics(df)
        
        # Step 6: Run BLAST for specificity
        df = run_blast_specificity(df)
        
        if df is None or len(df) == 0:
            logger.warning("No primers passed BLAST filtering. Exiting.")
            return False
        
        # Step 7: Save results to Excel file
        output_path = FileIO.save_results(
            df, 
            output_dir, 
            reference_file, 
            mode='standard'
        )
        
        if output_path:
            logger.info(f"\nResults saved to: {output_path}")
            return True
        else:
            logger.error("Failed to save results.")
            return False
            
    except SequenceProcessingError as e:
        error_msg = f"Sequence processing error: {str(e)}"
        logger.error(error_msg)
        return False
    except PrimerDesignError as e:
        error_msg = f"Primer design error: {str(e)}"
        logger.error(error_msg)
        return False
    except Exception as e:
        error_msg = f"Error in primer design workflow: {str(e)}"
        logger.error(error_msg)
        return False


def run_standard_mode(args):
    """
    Run the standard mode primer design workflow with file preparation.
    
    Modified version that includes automatic file preparation before
    proceeding with the primer design pipeline.
    """
    logger.info("=== Standard Mode Workflow ===")
    
    try:
        # Get input files if not provided in args
        if not args.fasta:
            logger.info("\n>>> Please select FASTA sequence file <<<")
            try:
                args.fasta = FileIO.select_fasta_file("Select FASTA sequence file")
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
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*.*")]
                )
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
                    [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("Compressed GFF Files", "*.gff.gz"), ("All Files", "*.*")]
                )
            except FileSelectionError as e:
                error_msg = f"GFF file selection failed: {str(e)}"
                logger.error(error_msg)
                return False
        elif args.noannotation:
            logger.info("\nSkipping GFF annotation file selection")
            args.gff = None

        # Signal that all file selections are complete
        FileIO.mark_selection_complete()
        
        # File preparation step
        try:
            from .utils import prepare_pipeline_files
            
            prep_result = prepare_pipeline_files(
                vcf_file=args.vcf,
                fasta_file=args.fasta,
                gff_file=args.gff
            )
            
            if not prep_result['success']:
                if prep_result.get('reason'):
                    logger.error(f"File preparation failed: {prep_result['reason']}")
                return False
            
            # Update file paths to use prepared files
            if prep_result['changes_made']:
                if 'vcf' in prep_result.get('prepared_files', {}):
                    args.vcf = prep_result['vcf_file']
                if 'fasta' in prep_result.get('prepared_files', {}):
                    args.fasta = prep_result['fasta_file']
                if 'gff' in prep_result.get('prepared_files', {}):
                    args.gff = prep_result['gff_file']
            else:
                logger.info("Original files are compatible and will be used as-is")
            
        except Exception as e:
            error_msg = f"File preparation failed: {str(e)}"
            logger.error(error_msg)
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
        
        # Extract variants from VCF file and process sequences
        logger.info("\nProcessing sequences with VCF variants...")
        try:
            # Load sequences from FASTA file first
            sequences = FileIO.load_fasta(args.fasta)
            
            # Process sequences with VCF using the new normalization approach
            processed_sequences = process_sequences_with_vcf(
                sequences=sequences,
                vcf_file=args.vcf,  # This is now the prepared VCF
                reference_file=args.fasta  # This might be the original or prepared FASTA
            )
            
        except Exception as e:
            error_msg = f"Error processing sequences with VCF: {str(e)}"
            logger.error(error_msg)
            raise SequenceProcessingError(error_msg) from e
        
        # Load gene annotations if needed
        if not args.noannotation:
            logger.info("\nLoading gene annotations from GFF file...")
            try:
                genes = AnnotationProcessor.load_genes_from_gff(args.gff)
            except Exception as e:
                error_msg = f"Error loading gene annotations: {str(e)}"
                logger.error(error_msg)
                return False
        else:
            logger.info("\nSkipping gene annotation loading (--noannotation specified)")
            genes = None
            
        # Run the primer design workflow
        success = run_primer_design_workflow(
            processed_sequences=processed_sequences,
            output_dir=output_dir,
            reference_file=args.fasta,
            genes=genes,
            gff_file=args.gff,
            skip_annotation_filtering=args.noannotation
        )
        
        return success
            
    except SequenceProcessingError as e:
        error_msg = f"Sequence processing error: {str(e)}"
        logger.error(error_msg)
        return False
    except Exception as e:
        error_msg = f"Error in standard mode workflow: {str(e)}"
        logger.error(error_msg)
        return False


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
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Setup logging
        log_file = setup_logging(debug=args.debug if args is not None else False)
        
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
        
        # Process BLAST database arguments
        if args.db is not None:
            try:
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
                        return False
                    
                    # Handle the case of selecting an existing database
                    if organism_key == 'existing_db':
                        # fasta_file actually contains the database path in this case
                        selected_db_path = fasta_file
                        if blast_db_manager.set_active_database(selected_db_path):
                            logger.info(f"\n=== Successfully selected BLAST database ===\n")
                            
                            # If only running database operations, exit successfully
                            if not args.fasta:
                                return True
                        else:
                            logger.error("Failed to set selected database as active")
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
                raise DDPrimerError(error_msg) from e
                
        # Verify BLAST database
        try:
            blast_db_manager = BlastDatabaseManager()
            if not blast_db_manager.verify_database():
                logger.warning("BLAST database verification failed. Attempting interactive setup...")
                if not blast_db_manager.setup_database_interactive():
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
        
        # Execute the standard mode workflow (only mode currently supported)
        success = run_standard_mode(args)
        
        if success:
            logger.info("\n=== Pipeline execution completed successfully! ===\n")
            return True
        else:
            logger.error("\n=== Pipeline execution failed ===")
            return False
            
    except DDPrimerError as e:
        # Handle application-specific exceptions
        error_msg = f"Pipeline error: {str(e)}"
        logger.error(error_msg)
        return False
    except Exception as e:
        # Handle unexpected exceptions
        error_msg = f"Unhandled exception during pipeline execution: {str(e)}"
        logger.error(error_msg)
        print(traceback.format_exc())
        return False


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