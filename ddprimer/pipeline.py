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
    parser.add_argument('--output', help='Output directory')
    parser.add_argument('--config', help='Configuration file')
    parser.add_argument('--cli', action='store_true', help='Force CLI mode')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
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

    if args.alignment and not (args.maf_file or args.second_fasta):
        pass
    elif args.alignment:
        if not args.maf_file and not args.second_fasta:
            parser.error("Cross-species alignment requires either --maf-file or --second-fasta")
        if not args.maf_file and not args.fasta:
            parser.error("Reference genome FASTA file (--fasta) is required")
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
    
    # Handle BLAST database creation (code remains the same)...
    
    logger.info("=== Primer Design Pipeline ===")
    
    try:
        # Check if cross-species alignment is enabled
        cross_species_mode = args is not None and args.alignment
        
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
                    # This assumes your MAF file is in an "Alignments" directory
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
            
            # Create temporary directory for intermediate files
            temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
            logger.debug(f"Created temporary directory: {temp_dir}")
            
            # ---------- IMPROVED CROSS-SPECIES WORKFLOW ----------
            
            # Step 1: Handle MAF file (pre-computed or generate new one)
            maf_file = args.maf_file
            if not maf_file:
                logger.info("\n>>> Running LastZ alignment between genomes <<<")

                lastz_options = args.lastz_options
                if "--format=maf" not in lastz_options and "-f maf" not in lastz_options:
                    lastz_options += " --format=maf"
                
                alignment_dir = os.path.dirname(output_dir)
                
                # Run LastZ alignment
                try:
                    lastz_runner = LastZRunner()
                    maf_file = lastz_runner.run_parallel_alignment(
                        args.fasta,
                        args.second_fasta,
                        alignment_dir,  # Use parent directory for alignments
                        lastz_options
                    )
                    logger.info(f"LastZ alignment completed: {maf_file}")
                except Exception as e:
                    logger.error(f"Error running LastZ alignment: {e}")
                    logger.debug(traceback.format_exc())
                    raise
            
            # Step 2: Parse MAF file and identify conserved regions
            logger.info("\n>>> Parsing alignment and identifying conserved regions <<<")
            try:
                # Initialize MAF parser
                maf_parser = MAFParser()
                
                # First analyze the MAF file to understand its structure
                maf_analysis = maf_parser.analyze_maf_file(maf_file)
                logger.debug(f"MAF file contains {maf_analysis['alignment_count']} alignment blocks")
                logger.debug(f"Reference sequences: {', '.join(maf_analysis['ref_seq_ids'][:5])}...")
                logger.debug(f"Query sequences: {', '.join(maf_analysis['query_seq_ids'][:5])}...")
                
                # Parse MAF file
                alignments = maf_parser.parse_maf_file(maf_file)
                
                # Identify conserved regions
                conserved_regions = maf_parser.identify_conserved_regions(
                    args.min_identity,
                    args.min_length
                )
                
                # Generate coordinate mapping between reference and query genomes
                coordinate_map = maf_parser.generate_coordinate_map(conserved_regions)
                logger.debug("Generated coordinate mapping between reference and second genome")
                
                total_regions = sum(len(regions) for regions in conserved_regions.values())
                logger.info(f"Identified {total_regions} conserved regions across {len(conserved_regions)} chromosomes")
            except Exception as e:
                logger.error(f"Error processing MAF file: {e}")
                logger.debug(traceback.format_exc())
                raise
            
            # Step 3: Prepare for masking by getting variants from both species
            ref_variants = {}
            second_variants = {}
            
            # Get reference genome variants if VCF provided
            if args.vcf:
                logger.info("\n>>> Extracting variants from reference genome VCF <<<")
                try:
                    snp_processor = SNPMaskingProcessor()
                    ref_variants = snp_processor.get_variant_positions(args.vcf)
                    total_ref_variants = sum(len(positions) for positions in ref_variants.values())
                    logger.info(f"Extracted {total_ref_variants} variants from reference genome")
                except Exception as e:
                    logger.error(f"Error extracting reference variants: {e}")
                    logger.debug(traceback.format_exc())
                    raise
            
            # Get second species variants if VCF provided
            if args.second_vcf:
                logger.info("\n>>> Extracting variants from second species VCF <<<")
                try:
                    snp_processor = SNPMaskingProcessor() if not 'snp_processor' in locals() else snp_processor
                    second_variants = snp_processor.get_variant_positions(args.second_vcf)
                    total_second_variants = sum(len(positions) for positions in second_variants.values())
                    logger.info(f"Extracted {total_second_variants} variants from second species genome")
                except Exception as e:
                    logger.error(f"Error extracting second species variants: {e}")
                    logger.debug(traceback.format_exc())
                    raise
            
            # Step 4: Load or generate reference sequences
            reference_sequences = {}
            
            # If we have a pre-computed MAF and don't need to run alignment
            if not args.fasta:
                # Generate reference sequences from the MAF file
                logger.info("\n>>> Generating reference sequences from MAF file <<<")
                try:
                    reference_sequences = maf_parser.extract_reference_sequences_from_maf()
                    logger.debug(f"Generated {len(reference_sequences)} reference sequences from MAF file")
                    
                    # Write to a temporary FASTA file for later use if needed
                    temp_ref_fasta = os.path.join(temp_dir, "ref_from_maf.fasta")
                    with open(temp_ref_fasta, 'w') as f:
                        for seq_id, seq in reference_sequences.items():
                            f.write(f">{seq_id}\n{seq}\n")
                    
                    logger.debug(f"Wrote reference sequences to temporary file: {temp_ref_fasta}")
                    
                except Exception as e:
                    logger.error(f"Error generating reference sequences from MAF: {e}")
                    logger.debug(traceback.format_exc())
                    raise
            else:
                # Load reference sequences from provided FASTA
                logger.info("\n>>> Loading reference sequences from FASTA file <<<")
                try:
                    reference_sequences = FileUtils.load_fasta(args.fasta)
                    logger.debug(f"Loaded {len(reference_sequences)} sequences from reference FASTA")
                except Exception as e:
                    logger.error(f"Error loading reference FASTA: {e}")
                    logger.debug(traceback.format_exc())
                    raise
            
            # Step 5: Apply masking to the reference sequences
            logger.info("\n>>> Creating masked genome with conserved regions and masked variants <<<")
            
            # First create a masked reference with only conserved regions - use temporary file
            masked_fasta_path = os.path.join(temp_dir, "masked_reference.fasta")
            try:
                maf_parser.mask_non_conserved_regions(
                    reference_sequences,
                    masked_fasta_path,
                    conserved_regions,
                    args.min_identity
                )
                logger.debug(f"Created alignment-masked reference genome: {masked_fasta_path}")
            except Exception as e:
                logger.error(f"Error masking non-conserved regions: {e}")
                logger.debug(traceback.format_exc())
                raise
            
            # Load the alignment-masked sequences
            alignment_masked_sequences = FileUtils.load_fasta(masked_fasta_path)
            logger.debug(f"Loaded {len(alignment_masked_sequences)} alignment-masked sequences")
            
            # Map second species variants to reference coordinates
            if second_variants:
                logger.info("\n>>> Mapping second species variants to reference coordinates <<<")
                try:
                    mapped_second_variants = maf_parser.map_second_variants_to_reference(second_variants, coordinate_map)
                    logger.debug(f"Mapped variants from second species to {len(mapped_second_variants)} reference chromosomes")
                except Exception as e:
                    logger.error(f"Error mapping second species variants: {e}")
                    logger.debug(traceback.format_exc())
                    mapped_second_variants = {}
            else:
                mapped_second_variants = {}
            
            # Now mask variants from both species
            final_masked_sequences = {}
            
            logger.info("\n>>> Masking variants in conserved regions <<<")
            # Process each sequence
            for seq_id, sequence in alignment_masked_sequences.items():
                # Get variants for reference genome
                ref_seq_variants = ref_variants.get(seq_id, set())
                
                # Get mapped variants from second species
                second_seq_variants = mapped_second_variants.get(seq_id, set())
                
                # Combine variants from both genomes
                all_variants = ref_seq_variants.union(second_seq_variants)
                logger.debug(f"Sequence {seq_id}: {len(ref_seq_variants)} reference variants, "
                           f"{len(second_seq_variants)} mapped second species variants, "
                           f"{len(all_variants)} total variants")
                
                # Mask all variants in the sequence
                if all_variants:
                    logger.debug(f"Masking {len(all_variants)} variants in {seq_id}...")
                    try:
                        variant_masked_seq = snp_processor.mask_variants(sequence, all_variants)
                        final_masked_sequences[seq_id] = variant_masked_seq
                    except Exception as e:
                        logger.error(f"Error masking variants in {seq_id}: {e}")
                        logger.debug(traceback.format_exc())
                        raise
                else:
                    logger.debug(f"No variants to mask in {seq_id}")
                    final_masked_sequences[seq_id] = sequence
            
            # Write final masked sequences to temporary file for debugging/verification
            final_masked_path = os.path.join(temp_dir, "final_masked.fasta")
            with open(final_masked_path, 'w') as f:
                for seq_id, seq in final_masked_sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")
            
            logger.info(f"Created fully masked reference genome: {final_masked_path}")
            
            # Use the final masked sequences for primer design
            masked_sequences = final_masked_sequences
            
            # Set GFF file for downstream processing
            gff_file = args.gff
            
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
        
        # Step 4: Load gene annotations
        logger.info("\nLoading gene annotations from GFF file...")
        try:
            genes = AnnotationProcessor.load_genes_from_gff(gff_file)
            logger.debug(f"Gene annotations loaded successfully from {gff_file}")
        except Exception as e:
            logger.error(f"Error loading gene annotations: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"Loaded {len(genes)} gene annotations")
        
        # Step 5: Filter sequences by restriction sites and gene overlap
        logger.info("\nFiltering sequences by restriction sites...")
        try:
            logger.debug(f"Using restriction site pattern: {Config.RESTRICTION_SITE}")
            restriction_fragments = SequenceProcessor.cut_at_restriction_sites(masked_sequences)
            logger.debug("Restriction site filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in restriction site filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"Generated {len(restriction_fragments)} fragments after restriction site cutting")
        
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
        
        # Step 6: Design primers with Primer3
        logger.info("\nDesigning primers with Primer3...")
        primer3_processor = Primer3Processor(Config)

        # Prepare input blocks for Primer3
        primer3_inputs = []
        fragment_info = {}  # Create a dictionary to store fragment information

        for fragment in filtered_fragments:
            # Store fragment information for later use
            fragment_info[fragment["id"]] = {
                "chr": fragment.get("chr", ""),
                "start": fragment.get("start", 1),
                "end": fragment.get("end", len(fragment["sequence"]))
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
            
            primer3_inputs.append(primer3_input)

        if not primer3_inputs:
            logger.warning("No valid fragments for primer design after N filtering. Exiting.")
            return None

        logger.debug(f"Running Primer3 on {len(primer3_inputs)} fragments...")
        try:
            logger.debug(f"Primer3 settings: Min Size={Config.PRIMER_MIN_SIZE}, Opt Size={Config.PRIMER_OPT_SIZE}, Max Size={Config.PRIMER_MAX_SIZE}")
            logger.debug(f"Primer3 Tm settings: Min={Config.PRIMER_MIN_TM}, Opt={Config.PRIMER_OPT_TM}, Max={Config.PRIMER_MAX_TM}")
            
            # Use parallel processing with progress bar
            primer3_output = primer3_processor.run_primer3_batch_parallel(primer3_inputs)
            
            # Pass fragment_info to the parse method
            primer_results = primer3_processor.parse_primer3_batch(primer3_output, fragment_info)
            logger.debug("Primer3 execution completed successfully")
        except Exception as e:
            logger.error(f"Error running Primer3: {e}")
            logger.debug(traceback.format_exc())
            raise
        
        # Step 7: Filter primers
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
            logger.debug(f"Filtering by GC content: Min={Config.PRIMER_MIN_GC}, Max={Config.PRIMER_MAX_GC}")
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
        
        # Step 8: Run NUPACK for thermodynamics
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
        
        # Step 9: Run BLAST for specificity
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
        
        # Step 10: Map primer coordinates to second genome if cross-species alignment was performed
        if args is not None and args.alignment and coordinate_map:
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
                    start_pos = parts[0]
                else:
                    start_pos = str(row["Location"])
                
                if chrom in coordinate_map:
                    # Try to map forward primer position
                    try:
                        for pos in range(int(start_pos), int(start_pos) + len(row["Primer F"])):
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
        
        # Step 11: Save results
        logger.debug("\nSaving results...")

        # Create filename
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
        else:
            fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
            output_file = os.path.join(output_dir, f"Primers_{fasta_name}.xlsx")
            
        logger.debug(f"Output file will be: {output_file}")

        # Extract gene names from sequence IDs
        logger.debug("Extracting gene names from sequence IDs")
        if "Sequence" in df.columns:
            # Create a new "Gene" column with just the gene name using the function from AnnotationProcessor
            df["Gene"] = df["Sequence"].apply(AnnotationProcessor.extract_gene_name)
            # Remove the original "Sequence" column
            df = df.drop("Sequence", axis=1)

        # Define column order
        columns = [
            "Gene", "Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST",
            "Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST",
            "Pair Penalty", "Amplicon", "Length", "Amplicon GC%", "Amplicon dG", "Chromosome", "Location",
        ]

        # Add second genome columns if available
        if args is not None and args.alignment and coordinate_map:
            second_genome_cols = [
                "Qry Chromosome", "Qry Location"
            ]
            columns.extend(second_genome_cols)
            logger.debug(f"Added second genome columns to output: {', '.join(second_genome_cols)}")

        # Add probe columns if present
        if "Probe" in df.columns:
            probe_cols = [
                "Probe", "Probe Tm", "Probe Penalty", "Probe dG", 
                "Probe BLAST"
            ]
            # Insert probe columns after primer columns
            idx = columns.index("Pair Penalty") + 1
            for col in reversed(probe_cols):
                if col in df.columns:
                    columns.insert(idx, col)
            
            logger.debug(f"Added probe columns to output: {', '.join(probe_cols)}")

        # Ensure all columns in the list exist in the DataFrame
        columns = [col for col in columns if col in df.columns]

        # Reorder the columns
        df = df[columns]

        # Save with formatting using the utility function
        try:
            output_path = FileUtils.save_formatted_excel(df, output_file, logger=logger)
            logger.info(f"Results saved to: {output_path}")
        except Exception as e:
            # Fallback if importing FileUtils fails for some reason
            logger.error(f"Error using FileUtils: {e}")
            logger.warning("Falling back to basic Excel export")
            df.to_excel(output_file, index=False)
            logger.info(f"Results saved to: {output_file} (without formatting)")

        try:
            if temp_dir and os.path.exists(temp_dir):
                logger.debug(f"Cleaning up temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir)
        except Exception as e:
            logger.warning(f"Error cleaning up temporary files: {e}")

        return output_path

    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        logger.debug(traceback.format_exc())
        raise



if __name__ == "__main__":
    try:
        output_file = run_pipeline()
        if output_file:
            logger.info("\n=== Pipeline completed successfully ===")
        else:
            logger.info("\n=== Pipeline completed with no results ===")
    except Exception as e:
        logger.error(f"\n!!! Pipeline failed: {str(e)}")
        traceback.print_exc()