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
import pandas as pd
import argparse
import logging
import traceback
from datetime import datetime
from tqdm import tqdm
from Bio import SeqIO

# Import package modules using the __init__.py structure
from .config import Config
from .utils import (
    FileUtils, 
    BlastDBCreator,
    SequenceUtils
)
from .core import (
    SNPMaskingProcessor, 
    Primer3Processor, 
    AnnotationProcessor, 
    BlastProcessor, 
    NupackProcessor, 
    SequenceProcessor, 
    PrimerProcessor
)
from .cross_species import ( 
    CrossSpeciesWorkflow,
    LastZRunner, 
    MAFParser
)

# Set up logging
logger = logging.getLogger("ddPrimer")

def setup_logging(debug=False):
    """Configure logging for the application."""
    # Use the Config.DEBUG_MODE if not explicitly set via command line
    debug_enabled = debug or Config.DEBUG_MODE
    
    log_level = logging.DEBUG if debug_enabled else logging.INFO
    
    # Different log formats based on debug mode
    if debug_enabled:
        log_format = '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
    else:
        log_format = '%(message)s'  # Simpler format for regular use
    
    # Set up logging to file
    log_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer", "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"ddPrimer_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    # Create file handler (always uses detailed format)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG if debug_enabled else logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
    ))
    
    # Create console handler with appropriate format
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(logging.Formatter(log_format))
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if debug_enabled else logging.INFO)
    
    # Clear existing handlers to avoid duplicates
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Add the handlers
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)
    
    # Configure our specific logger
    logger = logging.getLogger("ddPrimer")
    
    if debug_enabled:
        logger.debug(f"Debug mode: {debug_enabled}")
        logger.debug(f"Log file: {log_file}")
        logger.debug(f"Python version: {sys.version}")
        logger.debug(f"Platform: {sys.platform}")
    else:
        logger.info(f"Log file: {log_file}")
    
    # Set Config.DEBUG_MODE to ensure consistency
    Config.DEBUG_MODE = debug_enabled
    
    return log_file

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='ddPrimer: A pipeline for primer design and filtering')

    parser.add_argument('--fasta', help='Input FASTA file')
    parser.add_argument('--vcf', help='VCF file with variants')
    parser.add_argument('--gff', help='GFF annotation file')
    parser.add_argument('--direct', nargs='?', const=True, help='CSV or Excel file with sequence name and sequence columns (shortcut mode)')
    parser.add_argument('--output', help='Output directory')
    parser.add_argument('--config', help='Configuration file')
    parser.add_argument('--cli', action='store_true', help='Force CLI mode')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('--nooligo', action='store_true', help='Disable internal oligo (probe) design')
    # BLAST database creation
    parser.add_argument('--dbfasta', help='Create a BLAST database from this FASTA file (overrides config)')
    parser.add_argument('--dbname', help='Custom name for the BLAST database (optional)')
    parser.add_argument('--dboutdir', help='Custom output directory for the BLAST database (optional)')

    # New arguments for cross-species primer design
    alignment_group = parser.add_argument_group("Cross-Species Alignment Options")
    alignment_group.add_argument("--alignment", action="store_true", 
                            help="Enable cross-species primer design workflow")
    alignment_group.add_argument("--maf-file", 
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
    if args.direct and (args.fasta or args.vcf or args.gff or args.alignment):
        parser.error("--direct cannot be used with --fasta, --vcf, --gff, or --alignment options")
    
    if args.alignment and not (args.maf_file or args.second_fasta):
        pass
    elif args.alignment:
        if not args.maf_file and not args.second_fasta:
            parser.error("Cross-species alignment requires either --maf-file or --second-fasta")
        if not args.maf_file and not args.fasta:
            parser.error("Reference genome FASTA (--fasta) is required for alignment")
        if args.second_fasta and not args.second_vcf:
            parser.error("Second species VCF file (--second-vcf) is required for variant filtering")
    
    return args


def run_pipeline():
    """Run the primer design pipeline."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logging
    log_file = setup_logging(debug=args.debug if args is not None else False)
    
    logger.debug("Starting pipeline execution")
    logger.debug(f"Arguments: {args}")
    logger.debug(f"Config settings: NUM_PROCESSES={Config.NUM_PROCESSES}, BATCH_SIZE={Config.BATCH_SIZE}")
    
    # Force CLI mode if specified
    if args is not None and args.cli:
        FileUtils.use_cli = True
        logger.debug("CLI mode enforced via command line argument")

    # Load custom configuration if provided
    if args is not None and args.config:
        logger.debug(f"Loading custom configuration from {args.config}")
        Config.load_from_file(args.config)

    # Apply nooligo setting if specified
    if args is not None and args.nooligo:
        logger.info("Internal oligo (probe) design is disabled")
        # Modify both settings
        Config.PRIMER3_SETTINGS["PRIMER_PICK_INTERNAL_OLIGO"] = 0
        Config.DISABLE_INTERNAL_OLIGO = True
    
    # Initialize the Primer3Processor early to avoid reference issues
    primer3_processor = Primer3Processor(Config)
    
    logger.info("=== Primer Design Pipeline ===")
    
    try:
        # Check if cross-species alignment is enabled
        cross_species_mode = args is not None and args.alignment
        # Check if direct mode is enabled
        direct_sequence_mode = args is not None and args.direct
        
        # Track reference file for determining output directory
        reference_file = None
        
        # Create temporary directory for intermediate files
        import tempfile
        import shutil
        temp_dir = None  # Will be created after output_dir is set
        
        masked_sequences = {}
        coordinate_map = None  # Will store mapping between genomes if cross-species workflow is used
        
        # ---------- CROSS-SPECIES MODE ----------
        if cross_species_mode:
            logger.info("=== Cross-Species Primer Design Mode Enabled ===")
            
            # Get required input files based on workflow
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
                        logger.debug(traceback.format_exc())
                        raise
                
                if not args.second_vcf:
                    logger.info("\n>>> Please select SECOND species VCF file <<<")
                    try:
                        args.second_vcf = FileUtils.get_file(
                            "Select SECOND species VCF file", 
                            [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                        )
                    except Exception as e:
                        logger.error(f"Error selecting second species VCF file: {e}")
                        logger.debug(traceback.format_exc())
                        raise
                
                if not args.gff:
                    logger.info("\n>>> Please select GFF annotation file <<<")
                    try:
                        args.gff = FileUtils.get_file(
                            "Select GFF annotation file", 
                            [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
                        )
                    except Exception as e:
                        logger.error(f"Error selecting GFF file: {e}")
                        logger.debug(traceback.format_exc())
                        raise
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
                        logger.debug(traceback.format_exc())
                        raise
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
                        logger.debug(traceback.format_exc())
                        raise
                
                if not args.second_fasta:
                    logger.info("\n>>> Please select SECOND species FASTA file <<<")
                    try:
                        args.second_fasta = FileUtils.get_file(
                            "Select SECOND species FASTA file", 
                            [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                        )
                    except Exception as e:
                        logger.error(f"Error selecting second species FASTA file: {e}")
                        logger.debug(traceback.format_exc())
                        raise
                
                if not args.second_vcf:
                    logger.info("\n>>> Please select SECOND species VCF file <<<")
                    try:
                        args.second_vcf = FileUtils.get_file(
                            "Select SECOND species VCF file", 
                            [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                        )
                    except Exception as e:
                        logger.error(f"Error selecting second species VCF file: {e}")
                        logger.debug(traceback.format_exc())
                        raise
                
                if not args.gff:
                    logger.info("\n>>> Please select GFF annotation file <<<")
                    try:
                        args.gff = FileUtils.get_file(
                            "Select GFF annotation file", 
                            [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
                        )
                    except Exception as e:
                        logger.error(f"Error selecting GFF file: {e}")
                        logger.debug(traceback.format_exc())
                        raise
            
            # Now determine output directory based on reference file
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
            
            # Call the CrossSpeciesWorkflow function to handle the cross-species primer design
            logger.info("\n>>> Running Cross-Species Primer Design Workflow <<<")
            try:
                # Create temporary directory for intermediate files if needed by the workflow
                if not temp_dir:
                    temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
                    logger.debug(f"Created temporary directory: {temp_dir}")
                
                masked_sequences, coordinate_map = CrossSpeciesWorkflow(args, output_dir, logger)
                logger.info(f"Cross-species workflow completed with {len(masked_sequences)} masked sequences")
            except Exception as e:
                logger.error(f"Error in cross-species workflow: {e}")
                logger.debug(traceback.format_exc())
                raise
                
            # Set GFF file for downstream processing
            gff_file = args.gff
            
        # ---------- DIRECT SEQUENCE MODE ----------
        elif direct_sequence_mode:
            logger.info("=== Direct Sequence Mode Enabled ===")
            
            # Get the input file - handle both when --direct is specified alone or with a file
            if args.direct is True:  # --direct was used without a specified file
                sequence_file = FileUtils.get_sequences_file()
            else:  # --direct was used with a specified file
                sequence_file = args.direct
                
            reference_file = sequence_file
            logger.debug(f"Direct sequence file selection successful: {sequence_file}")
            
            # Now determine output directory based on reference file
            if args.output:
                output_dir = args.output
            elif reference_file:
                # Use the directory of the reference file
                input_dir = os.path.dirname(os.path.abspath(reference_file))
                output_dir = os.path.join(input_dir, "Primers")
            else:
                # Fallback to current directory if no reference file
                output_dir = os.path.join(os.getcwd(), "Primers")
            
            # Create output directory
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created output directory: {output_dir}")
            
            # Create temporary directory for intermediate files
            temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
            logger.debug(f"Created temporary directory: {temp_dir}")
            
            # Load sequences directly from the provided file
            logger.info("\nLoading sequences from input file...")
            try:
                masked_sequences = FileUtils.load_sequences_from_table(sequence_file)
                logger.debug(f"Loaded {len(masked_sequences)} sequences from {sequence_file}")
            except Exception as e:
                logger.error(f"Error loading sequences from file: {e}")
                logger.debug(traceback.format_exc())
                raise
            
            # Create an empty genes dictionary since we don't have GFF
            genes = {}
            gff_file = None  # No GFF file in this mode
            
        # ---------- SINGLE SPECIES MODE ----------
        else:
            # Standard single-species mode
            logger.info(">>> Please select reference FASTA file <<<")
            try:
                fasta_file = args.fasta if args is not None and args.fasta else FileUtils.get_file(
                    "Select reference FASTA file", 
                    [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                )
                reference_file = fasta_file
                logger.debug(f"FASTA file selection successful: {fasta_file}")
            except Exception as e:
                logger.error(f"Error selecting FASTA file: {e}")
                logger.debug(traceback.format_exc())
                raise
            
            # Now determine output directory based on reference file
            if args.output:
                output_dir = args.output
            elif reference_file:
                # Use the directory of the reference file
                input_dir = os.path.dirname(os.path.abspath(reference_file))
                output_dir = os.path.join(input_dir, "Primers")
            else:
                # Fallback to current directory if no reference file
                output_dir = os.path.join(os.getcwd(), "Primers")
            
            # Create output directory
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created output directory: {output_dir}")
            
            # Create temporary directory for intermediate files
            temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
            logger.debug(f"Created temporary directory: {temp_dir}")
            
            # VCF file selection
            logger.info("\n>>> Please select VCF file with variants <<<")
            try:
                vcf_file = args.vcf if args is not None and args.vcf else FileUtils.get_file(
                    "Select VCF file with variants", 
                    [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                )
                logger.debug(f"VCF file selection successful: {vcf_file}")
            except Exception as e:
                logger.error(f"Error selecting VCF file: {e}")
                logger.debug(traceback.format_exc())
                raise
            
            # GFF file selection
            logger.info("\n>>> Please select GFF annotation file <<<")
            try:
                gff_file = args.gff if args is not None and args.gff else FileUtils.get_file(
                    "Select GFF annotation file", 
                    [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3"), ("All Files", "*")]
                )
                logger.debug(f"GFF file selection successful: {gff_file}")
            except Exception as e:
                logger.error(f"Error selecting GFF file: {e}")
                logger.debug(traceback.format_exc())
                raise
            
            # Extract variants and mask sequences
            logger.info("\nExtracting variants from VCF file...")
            snp_processor = SNPMaskingProcessor()
            try:
                variants = snp_processor.get_variant_positions(vcf_file)
                logger.debug(f"Variants extracted successfully from {vcf_file}")
            except Exception as e:
                logger.error(f"Error extracting variants: {e}")
                logger.debug(traceback.format_exc())
                raise
            
            total_variants = sum(len(positions) for positions in variants.values())
            logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes")
            
            logger.debug("\nLoading sequences from FASTA file...")
            try:
                sequences = FileUtils.load_fasta(fasta_file)
                logger.debug(f"Sequences loaded successfully from {fasta_file}")
            except Exception as e:
                logger.error(f"Error loading FASTA: {e}")
                logger.debug(traceback.format_exc())
                raise
            logger.debug(f"Loaded {len(sequences)} sequences")
            
            logger.debug("\nMasking variants in sequences...")
            
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
                        logger.error(f"Error masking variants in {seq_id}: {e}")
                        logger.debug(traceback.format_exc())
                        raise
                else:
                    logger.debug(f"No variants to mask in {seq_id}")
                    masked_sequences[seq_id] = sequence
        
        # ----- COMMON PROCESSING FOR BOTH MODES AFTER MASKED SEQUENCES ARE GENERATED -----
        
        # Step 1: Load gene annotations (skip in direct mode)
        if not direct_sequence_mode:
            logger.info("\nLoading gene annotations from GFF file...")
            try:
                genes = AnnotationProcessor.load_genes_from_gff(gff_file)
                logger.debug(f"Gene annotations loaded successfully from {gff_file}")
                logger.info(f"Loaded {len(genes)} gene annotations")
            except Exception as e:
                logger.error(f"Error loading gene annotations: {e}")
                logger.debug(traceback.format_exc())
                raise
                
        # Step 2: Filter sequences by restriction sites
        logger.info("\nFiltering sequences by restriction sites...")
        try:
            logger.debug(f"Using restriction site pattern: {Config.RESTRICTION_SITE}")
            
            # Log sequence stats before restriction site cutting
            logger.debug(f"Sequences before restriction site cutting:")
            for seq_id, seq in masked_sequences.items():
                logger.debug(f"  {seq_id}: {len(seq)} bp")
            
            # Use the standard restriction site method for all modes
            restriction_fragments = SequenceProcessor.cut_at_restriction_sites(masked_sequences)
            
            # Log detailed information about restriction fragments
            logger.debug("Restriction fragments after cutting:")
            for fragment in restriction_fragments:
                logger.debug(f"  {fragment['id']}: {len(fragment['sequence'])} bp, chr={fragment.get('chr', 'NA')}, "
                            f"start={fragment.get('start', 'NA')}, end={fragment.get('end', 'NA')}")
            
            logger.debug("Restriction site filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in restriction site filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"Generated {len(restriction_fragments)} fragments after restriction site cutting")

        # Step 3: Filter by gene overlap (skip in direct mode)
        if direct_sequence_mode:
            logger.debug("Skipping gene overlap filtering in direct sequence mode")
            
            # Process fragments for direct mode, simplifying structure
            simplified_fragments = []
            for fragment in restriction_fragments:
                # Extract the base sequence ID (without any fragment numbers)
                base_id = fragment["id"].split("_frag")[0]
                
                # Create a simplified fragment with ONLY what we need
                # Explicitly exclude location data for direct mode
                simplified_fragment = {
                    "id": fragment["id"],
                    "sequence": fragment["sequence"],
                    "Gene": base_id  # Explicitly set Gene to base_id for direct mode
                }
                simplified_fragments.append(simplified_fragment)
            
            filtered_fragments = simplified_fragments
            logger.debug(f"Prepared {len(filtered_fragments)} fragments for direct mode")
            
            # Log fragment contents
            logger.debug("Direct mode fragments after processing:")
            for i, frag in enumerate(filtered_fragments):
                if i < 5 or i >= len(filtered_fragments) - 5:  # Log first and last 5 fragments
                    logger.debug(f"  {frag['id']}: Gene={frag['Gene']}, length={len(frag['sequence'])} bp")
                elif i == 5 and len(filtered_fragments) > 10:
                    logger.debug(f"  ... ({len(filtered_fragments) - 10} more fragments) ...")
        else:
            logger.info("\nFiltering sequences by gene overlap...")
            try:
                logger.debug(f"Using gene overlap margin: {Config.GENE_OVERLAP_MARGIN}")
                filtered_fragments = SequenceProcessor.filter_by_gene_overlap(restriction_fragments, genes)
                logger.debug("Gene overlap filtering completed successfully")
            except Exception as e:
                logger.error(f"Error in gene overlap filtering: {e}")
                logger.debug(traceback.format_exc())
                raise
            logger.info(f"Retained {len(filtered_fragments)} fragments after gene overlap filtering")
        
        if not filtered_fragments:
            logger.warning("No valid fragments for primer design. Exiting.")
            return None
        
        if direct_sequence_mode and Config.DEBUG_MODE:
            logger.debug("\n=== DIRECT MODE DIAGNOSTIC INFO ===")
            logger.debug(f"Number of input sequences: {len(masked_sequences)}")
            logger.debug(f"Number of fragments after processing: {len(filtered_fragments)}")
            
            # Check which sequences changed during processing
            changed_sequences = []
            unchanged_sequences = []
            
            for fragment in filtered_fragments:
                seq_id = fragment["id"]
                if seq_id in masked_sequences:
                    original_seq = masked_sequences[seq_id]
                    current_seq = fragment["sequence"]
                    
                    if original_seq == current_seq:
                        unchanged_sequences.append(seq_id)
                    else:
                        changed_sequences.append({
                            "id": seq_id,
                            "original_length": len(original_seq),
                            "current_length": len(current_seq),
                            "difference": len(original_seq) - len(current_seq)
                        })
                else:
                    # This is a fragment derived from an original sequence
                    logger.debug(f"Fragment {seq_id} is derived from an original sequence")
            
            if changed_sequences:
                logger.debug("Sequences changed during processing:")
                for info in changed_sequences:
                    logger.debug(f"  {info['id']}: {info['original_length']} → {info['current_length']} (diff: {info['difference']})")
            else:
                logger.debug("All sequences preserved their original content")
            
            # Check specific Primer3 parameters to ensure they match between modes
            logger.debug("\nPrimer3 Parameter Summary:")
            key_params = ["PRIMER_PRODUCT_SIZE_RANGE", "PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE", 
                        "PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM", "PRIMER_MIN_GC", "PRIMER_MAX_GC"]
            
            settings = primer3_processor.config.PRIMER3_SETTINGS
            for param in key_params:
                value = settings.get(param, "Not set")
                logger.debug(f"  {param}: {value}")
            
            logger.debug("=== END DIAGNOSTIC INFO ===\n")

        # Step 4: Design primers with Primer3
        logger.info("\nDesigning primers with Primer3...")

        # Prepare input blocks for Primer3
        primer3_inputs = []
        fragment_info = {}  # Create a dictionary to store fragment information

        # Log fragment information before Primer3 processing
        logger.debug(f"Preparing {len(filtered_fragments)} fragments for Primer3...")

        for fragment in filtered_fragments:
            # Store fragment information differently based on mode
            if direct_sequence_mode:
                # For direct mode, only store minimal information - explicitly exclude location
                fragment_info[fragment["id"]] = {
                    "gene": fragment.get("Gene", fragment["id"])
                }
            else:
                # For regular mode, store location information
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
            
            # Add target region if sequence is long
            if len(fragment["sequence"]) > 200:
                target_start = len(fragment["sequence"]) // 4
                target_len = len(fragment["sequence"]) // 2
                primer3_input["SEQUENCE_TARGET"] = [target_start, target_len]
                logger.debug(f"Added target region [{target_start}, {target_len}] for {fragment['id']}")
            
            # Log detailed information for diagnostic purposes
            if len(primer3_inputs) < 5:  # Only log the first 5 sequences to avoid overwhelming logs
                logger.debug(f"Primer3 input for {fragment['id']}:")
                logger.debug(f"  Sequence length: {len(fragment['sequence'])} bp")
                logger.debug(f"  Sequence: {fragment['sequence'][:50]}...{fragment['sequence'][-50:] if len(fragment['sequence']) > 100 else ''}")
                if direct_sequence_mode:
                    logger.debug(f"  Gene name: {fragment.get('Gene', 'Not set')}")
                else:
                    logger.debug(f"  Location: chr={fragment.get('chr', 'NA')}, start={fragment.get('start', 'NA')}, end={fragment.get('end', 'NA')}")
            
            primer3_inputs.append(primer3_input)

        # Log a sample of the fragment_info dictionary
        logger.debug("Fragment info sample (first 5 entries):")
        sample_entries = list(fragment_info.items())[:5]
        for frag_id, info in sample_entries:
            logger.debug(f"  {frag_id}: {info}")

        # Run primer3
        if not primer3_inputs:
            logger.warning("No valid fragments for primer design after filtering. Exiting.")
            return None

        logger.debug(f"Running Primer3 on {len(primer3_inputs)} fragments...")
        try:
            logger.debug(f"Primer3 settings: Min Size={Config.PRIMER_MIN_SIZE}, Opt Size={Config.PRIMER_OPT_SIZE}, Max Size={Config.PRIMER_MAX_SIZE}")
            logger.debug(f"Primer3 Tm settings: Min={Config.PRIMER_MIN_TM}, Opt={Config.PRIMER_OPT_TM}, Max={Config.PRIMER_MAX_TM}")
            
            # Use parallel processing with progress bar
            primer3_output = primer3_processor.run_primer3_batch_parallel(primer3_inputs)
            logger.debug(f"Primer3 execution completed, processing results...")
            
            # Pass fragment_info to the parse method
            primer_results = primer3_processor.parse_primer3_batch(primer3_output, fragment_info)
            logger.debug(f"Primer3 execution completed successfully, found {len(primer_results)} primer pairs")
            
            # Log a sample of the primer results
            if primer_results and len(primer_results) > 0:
                logger.debug("Primer results sample (first result):")
                sample_result = primer_results[0]
                for key, value in sample_result.items():
                    if key in ["Primer F", "Primer R", "Amplicon", "Gene"]:
                        logger.debug(f"  {key}: {value}")
                
        except Exception as e:
            logger.error(f"Error running Primer3: {e}")
            logger.debug(traceback.format_exc())
            raise

        if not primer_results or len(primer_results) == 0:
            logger.warning("No primers were designed by Primer3. Exiting.")
            return None
        
        # Step 5: Filter primers
        logger.info("\nFiltering primers...")
        
        # Convert to DataFrame
        df = pd.DataFrame(primer_results)
        initial_count = len(df)
        logger.debug(f"Initial primer count: {initial_count}")
        
        # Filter by penalty
        try:
            logger.debug(f"Filtering by penalty with threshold: {Config.PENALTY_MAX}")
            df = PrimerProcessor.filter_by_penalty(df)
            logger.debug("Penalty filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in penalty filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.debug(f"After penalty filtering: {len(df)}/{initial_count} primers")
        
        # Filter by repeats
        try:
            logger.debug("Filtering primers by repeat sequences")
            df = PrimerProcessor.filter_by_repeats(df)
            logger.debug("Repeat filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in repeat filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.debug(f"After repeat filtering: {len(df)}/{initial_count} primers")
        
        # Filter by GC content
        try:
            logger.debug(f"Filtering by GC content: Min={Config.SEQUENCE_MIN_GC}, Max={Config.SEQUENCE_MAX_GC}")
            df = PrimerProcessor.filter_by_gc_content(df)
            logger.debug("GC content filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in GC content filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"After filtering: {len(df)}/{initial_count} primers")
        
        # Process internal oligos (reverse complement if needed)
        try:
            logger.debug("Processing internal oligos")
            df = PrimerProcessor.process_internal_oligos(df)
            logger.debug("Internal oligo processing completed successfully")
        except Exception as e:
            logger.error(f"Error in internal oligo processing: {e}")
            logger.debug(traceback.format_exc())
            raise
        
        if len(df) == 0:
            logger.warning("No primers passed filtering. Exiting.")
            return None
        
        # Step 6: Run NUPACK for thermodynamics
        logger.info("\nCalculating thermodynamic properties with NUPACK...")
        
        # Calculate deltaG for forward primers
        logger.debug("Calculating ΔG for forward primers...")
        try:
            logger.debug(f"NUPACK settings: Temp={Config.NUPACK_TEMPERATURE}°C, Na={Config.NUPACK_SODIUM}M, Mg={Config.NUPACK_MAGNESIUM}M")
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc="Processing forward primers")
                df["Primer F dG"] = df["Primer F"].progress_apply(NupackProcessor.calc_deltaG)
            else:
                df["Primer F dG"] = df["Primer F"].apply(NupackProcessor.calc_deltaG)
            logger.debug("Forward primer deltaG calculation completed successfully")
        except Exception as e:
            logger.error(f"Error calculating forward primer deltaG: {e}")
            logger.debug(traceback.format_exc())
            raise
        
        # Calculate deltaG for reverse primers
        logger.debug("Calculating ΔG for reverse primers...")
        try:
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc="Processing reverse primers")
                df["Primer R dG"] = df["Primer R"].progress_apply(NupackProcessor.calc_deltaG)
            else:
                df["Primer R dG"] = df["Primer R"].apply(NupackProcessor.calc_deltaG)
            logger.debug("Reverse primer deltaG calculation completed successfully")
        except Exception as e:
            logger.error(f"Error calculating reverse primer deltaG: {e}")
            logger.debug(traceback.format_exc())
            raise
        
        # Calculate deltaG for probes if present
        if "Probe" in df.columns:
            logger.debug("Calculating ΔG for probes...")
            try:
                if Config.SHOW_PROGRESS:
                    tqdm.pandas(desc="Processing probes")
                    df["Probe dG"] = df["Probe"].progress_apply(lambda x: 
                                                  NupackProcessor.calc_deltaG(x) 
                                                  if pd.notnull(x) and x else None)
                else:
                    df["Probe dG"] = df["Probe"].apply(lambda x: 
                                                     NupackProcessor.calc_deltaG(x) 
                                                     if pd.notnull(x) and x else None)
                logger.debug("Probe deltaG calculation completed successfully")
            except Exception as e:
                logger.error(f"Error calculating probe deltaG: {e}")
                logger.debug(traceback.format_exc())
                raise
        
        # Calculate deltaG for amplicons
        logger.debug("Calculating ΔG for amplicons...")
        try:
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc="Processing amplicons")
                df["Amplicon dG"] = df["Amplicon"].progress_apply(NupackProcessor.calc_deltaG)
            else:
                df["Amplicon dG"] = df["Amplicon"].apply(NupackProcessor.calc_deltaG)
            logger.debug("Amplicon deltaG calculation completed successfully")
        except Exception as e:
            logger.error(f"Error calculating amplicon deltaG: {e}")
            logger.debug(traceback.format_exc())
            raise
        
        # Step 7: Run BLAST for specificity
        logger.info("\nRunning BLAST for specificity checking...")
        
        # Run BLAST for forward primers
        try:
            logger.debug(f"BLAST settings: Word Size={Config.BLAST_WORD_SIZE}, E-value={Config.BLAST_EVALUE}")
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
            logger.debug("Forward primer BLAST completed successfully")
        except Exception as e:
            logger.error(f"Error in forward primer BLAST: {e}")
            logger.debug(traceback.format_exc())
            raise
        
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
            logger.debug("Reverse primer BLAST completed successfully")
        except Exception as e:
            logger.error(f"Error in reverse primer BLAST: {e}")
            logger.debug(traceback.format_exc())
            raise
        
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
                logger.debug("Probe BLAST completed successfully")
            except Exception as e:
                logger.error(f"Error in probe BLAST: {e}")
                logger.debug(traceback.format_exc())
                raise
        
        # Filter by BLAST specificity
        try:
            logger.debug(f"Filtering by BLAST specificity, filter factor: {Config.BLAST_FILTER_FACTOR}")
            df = PrimerProcessor.filter_by_blast(df)
            logger.debug("BLAST filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in BLAST filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"After BLAST filtering: {len(df)}/{initial_count} primers")
        
        if len(df) == 0:
            logger.warning("No primers passed BLAST filtering. Exiting.")
            return None
        
        # Step 8: Map primer coordinates to second genome if cross-species alignment was performed
        if cross_species_mode and coordinate_map:
            logger.debug("\nMapping primer coordinates to second genome...")
            
            # Add columns for second genome data
            df["Qry Chromosome"] = None
            df["Qry Location"] = None
            mapped_count = 0
            
            for idx, row in df.iterrows():
                chrom = row["Chromosome"]
                
                # Handle the case where Location might not contain a dash
                if "-" in str(row["Location"]):
                    parts = row["Location"].split("-")
                    start_pos = int(parts[0])
                else:
                    try:
                        start_pos = int(row["Location"])
                    except (ValueError, TypeError):
                        logger.warning(f"Invalid location format for index {idx}: {row['Location']}")
                        continue
                
                if chrom in coordinate_map:
                    # Try to map forward primer position
                    try:
                        for pos in range(start_pos, start_pos + len(row["Primer F"])):
                            if pos in coordinate_map[chrom]:
                                mapping = coordinate_map[chrom][pos]
                                df.at[idx, "Qry Chromosome"] = mapping["qry_src"]
                                df.at[idx, "Qry Location"] = mapping["qry_pos"]
                                mapped_count += 1
                                break
                    except Exception as e:
                        logger.warning(f"Error mapping primer at index {idx}: {e}")
                        continue
            
            logger.debug(f"Successfully mapped {mapped_count}/{len(df)} primers to second genome")
        
        # For direct mode, check if any sequences didn't get primers and add them to the results
        if direct_sequence_mode:
            logger.debug("Checking for sequences without primers...")
            
            # Get all sequence IDs that have primers
            sequences_with_primers = set()
            if "Gene" in df.columns:
                # Get all sequence names from the Gene column that have valid primers
                valid_primers_mask = df["Primer F"] != "No suitable primers found"
                valid_genes = df.loc[valid_primers_mask, "Gene"].astype(str).unique()
                sequences_with_primers.update(valid_genes)
            
            # Get all input sequence IDs
            all_input_sequences = set(masked_sequences.keys())
            
            # Find sequences without primers
            sequences_without_primers = all_input_sequences - sequences_with_primers
            
            logger.debug(f"All input sequences: {len(all_input_sequences)}")
            logger.debug(f"Sequences with primers: {len(sequences_with_primers)}")
            logger.debug(f"Sequences without primers: {len(sequences_without_primers)}")
            
            if sequences_without_primers:
                logger.debug(f"Adding rows for {len(sequences_without_primers)} sequences without primers")
                
                # Create rows for sequences without primers
                no_primer_rows = []
                for seq_id in sequences_without_primers:
                    # Create a row with only the necessary columns to avoid dtype issues
                    row = {
                        "Gene": seq_id,
                        "Primer F": "No suitable primers found"
                    }
                    
                    # Add the row to our list
                    no_primer_rows.append(row)
                
                # Add these rows to the DataFrame
                if no_primer_rows:
                    # Create a new DataFrame with just the minimal set of columns
                    no_primer_df = pd.DataFrame(no_primer_rows)
                    
                    # Concatenate with the original DataFrame
                    # This avoids the FutureWarning by not including empty columns
                    df = pd.concat([df, no_primer_df], ignore_index=True, sort=False)
                    logger.debug(f"Added {len(no_primer_rows)} rows for sequences without primers")
            
        
        # Step 9: Save results
        logger.debug("\nSaving results...")

        # Create filename based on mode
        if cross_species_mode:
            if args.fasta:
                ref_fasta_name = os.path.splitext(os.path.basename(args.fasta))[0]
            else:
                ref_fasta_name = "reference"
            
            if args.second_fasta:
                second_fasta_name = os.path.splitext(os.path.basename(args.second_fasta))[0]
            else:
                second_fasta_name = "query"
            
            # Create filename with both species names
            output_file = os.path.join(output_dir, f"Primers_{ref_fasta_name}_vs_{second_fasta_name}.xlsx")
        elif direct_sequence_mode:
            # For direct sequence mode
            input_file_name = os.path.splitext(os.path.basename(reference_file))[0]
            output_file = os.path.join(output_dir, f"Primers_{input_file_name}.xlsx")
        else:
            # Standard mode
            if reference_file:
                fasta_name = os.path.splitext(os.path.basename(reference_file))[0]
                output_file = os.path.join(output_dir, f"Primers_{fasta_name}.xlsx")
            else:
                # Use a generic name if reference_file is None
                output_file = os.path.join(output_dir, "Primers_output.xlsx")

        logger.debug(f"Output file will be: {output_file}")

        # Define a base set of columns that should appear in all modes
        base_columns = [
            "Gene", "Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST1", "Primer F BLAST2",
            "Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST1", "Primer R BLAST2",
            "Pair Penalty", "Amplicon", "Length", "Amplicon GC%", "Amplicon dG"
        ]

        # Add location columns only if not in direct mode
        if not direct_sequence_mode:
            location_columns = ["Chromosome", "Location"]
            for col in location_columns:
                if col in df.columns:
                    base_columns.append(col)

            # Add cross-species columns if available
            cross_species_cols = ["Qry Chromosome", "Qry Location"]
            for col in cross_species_cols:
                if col in df.columns and not df[col].isna().all():
                    base_columns.append(col)

        # Add probe columns if present
        if "Probe" in df.columns:
            probe_cols = [
                "Probe", "Probe Tm", "Probe Penalty", "Probe dG", 
                "Probe BLAST1", "Probe BLAST2"
            ]
            # Insert probe columns after primer columns
            idx = base_columns.index("Pair Penalty") + 1
            for col in reversed(probe_cols):
                if col in df.columns:
                    base_columns.insert(idx, col)
            
            logger.debug(f"Added probe columns to output: {', '.join([c for c in probe_cols if c in df.columns])}")

        # For direct mode, explicitly remove any location columns that might have been added
        if direct_sequence_mode:
            columns = [col for col in base_columns if col not in ['Chromosome', 'Location', 'Qry Chromosome', 'Qry Location']]
        else:
            columns = base_columns

        # Ensure all columns in the list exist in the DataFrame
        columns = [col for col in columns if col in df.columns]

        # Reorder the columns
        df = df[columns]

        # Debug the final column list
        logger.debug(f"Final output columns: {', '.join(columns)}")

        # Save with formatting using the utility function
        try:
            output_path = FileUtils.save_formatted_excel(df, output_file, logger=logger)
            logger.info(f"Results saved to: {output_path}")
            success = True
        except Exception as e:
            # Fallback if importing FileUtils fails for some reason
            logger.error(f"Error saving Excel file: {e}")
            logger.warning("Falling back to basic Excel export")
            
            try:
                df.to_excel(output_file, index=False)
                logger.info(f"Results saved to: {output_file} (without formatting)")
                success = True
            except Exception as ex:
                logger.error(f"Failed to save results: {ex}")
                success = False
        
        # Cleanup temporary files
        try:
            if temp_dir and os.path.exists(temp_dir):
                logger.debug(f"Cleaning up temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir)
        except Exception as e:
            logger.warning(f"Error cleaning up temporary files: {e}")
    
    except Exception as e:
        logger.error(f"Unhandled exception during pipeline execution: {e}")
        logger.debug(traceback.format_exc())
        raise
    
    return success