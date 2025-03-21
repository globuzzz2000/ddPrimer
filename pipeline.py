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

# Import package modules using the __init__.py structure
from .config import Config
from .utils import FileUtils, UIUtils
from .core import (
    SNPMaskingProcessor, 
    Primer3Processor, 
    AnnotationProcessor, 
    BlastProcessor, 
    NupackProcessor, 
    SequenceProcessor, 
    PrimerProcessor
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
    return parser.parse_args()


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
    
    logger.info("=== Primer Design Pipeline ===")
    
    try:
        # Step 1: Get input files
        logger.info("Selecting input files...")
        
        # FASTA file selection
        logger.info(">>> Please select reference FASTA file (file dialog will open) <<<")
        fasta_file = None
        try:
            fasta_file = args.fasta if args is not None and args.fasta else FileUtils.get_file(
                "Select reference FASTA file", 
                [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
            )
            logger.debug(f"FASTA file selection successful: {fasta_file}")
        except Exception as e:
            logger.error(f"Error selecting FASTA file: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.debug(f"Selected FASTA: {fasta_file}")

        # VCF file selection
        logger.info("\n>>> Please select VCF file with variants (file dialog will open) <<<")
        vcf_file = None
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
        logger.debug(f"Selected VCF: {vcf_file}")
        
        # GFF file selection
        logger.info("\n>>> Please select GFF annotation file (file dialog will open) <<<")
        gff_file = None
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
        logger.debug(f"Selected GFF: {gff_file}")
        
        # Step 2: Extract variants and mask sequences
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
        
        logger.info("\nLoading sequences from FASTA file...")
        try:
            sequences = FileUtils.load_fasta(fasta_file)
            logger.debug(f"Sequences loaded successfully from {fasta_file}")
        except Exception as e:
            logger.error(f"Error loading FASTA: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"Loaded {len(sequences)} sequences")
        
        logger.info("\nMasking variants in sequences...")
        masked_sequences = {}
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
        
        # Step 3: Load gene annotations
        logger.info("\nLoading gene annotations from GFF file...")
        try:
            genes = AnnotationProcessor.load_genes_from_gff(gff_file)
            logger.debug(f"Gene annotations loaded successfully from {gff_file}")
        except Exception as e:
            logger.error(f"Error loading gene annotations: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"Loaded {len(genes)} gene annotations")
        
        # Step 4: Filter sequences by restriction sites and gene overlap
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
        
        # Step 5: Design primers with Primer3
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

        logger.info(f"Running Primer3 on {len(primer3_inputs)} fragments...")
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
        
        # Step 6: Filter primers
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
        logger.info(f"After penalty filtering: {len(df)}/{initial_count} primers")
        
        # Filter by repeats
        try:
            logger.debug("Filtering primers by repeat sequences")
            df = PrimerProcessor.filter_by_repeats(df)
            logger.debug("Repeat filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in repeat filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"After repeat filtering: {len(df)}/{initial_count} primers")
        
        # Filter by GC content
        try:
            logger.debug(f"Filtering by GC content: Min={Config.PRIMER_MIN_GC}, Max={Config.PRIMER_MAX_GC}")
            df = PrimerProcessor.filter_by_gc_content(df)
            logger.debug("GC content filtering completed successfully")
        except Exception as e:
            logger.error(f"Error in GC content filtering: {e}")
            logger.debug(traceback.format_exc())
            raise
        logger.info(f"After GC% filtering: {len(df)}/{initial_count} primers")
        
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
        
        # Step 7: Run NUPACK for thermodynamics
        logger.info("\nCalculating thermodynamic properties with NUPACK...")
        
        # Calculate deltaG for forward primers
        logger.info("Calculating ΔG for forward primers...")
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
        logger.info("Calculating ΔG for reverse primers...")
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
            logger.info("Calculating ΔG for probes...")
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
        logger.info("Calculating ΔG for amplicons...")
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
        
        # Step 8: Run BLAST for specificity
        logger.info("\nRunning BLAST for specificity checking...")
        
        # Run BLAST for forward primers
        logger.info("BLASTing forward primers...")
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
        logger.info("BLASTing reverse primers...")
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
            logger.info("BLASTing probes...")
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
        
        # Step 9: Save results
        logger.info("\nSaving results...")
        
        # Create output directory
        output_dir = args.output if args is not None and args.output else os.path.join(os.path.dirname(fasta_file), "primers")
        os.makedirs(output_dir, exist_ok=True)
        logger.debug(f"Created output directory: {output_dir}")
        
        # Create filename with timestamp
        fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
        vcf_name = os.path.splitext(os.path.basename(vcf_file))[0]
        
        output_file = os.path.join(output_dir, f"primers_{fasta_name}_{vcf_name}.xlsx")
        logger.debug(f"Output file will be: {output_file}")
        
        # Define column order
        columns = [
            "Sequence", "Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST",
            "Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST",
            "Pair Penalty", "Amplicon", "Length", "Amplicon GC%", "Amplicon dG", "Chromosome", "Location",
        ]
        
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
        
        # Select columns that exist in the DataFrame
        columns = [col for col in columns if col in df.columns]
        
        # Reorder and save
        df = df[columns]
        try:
            df.to_excel(output_file, index=False)
            logger.debug(f"Successfully wrote results to Excel file: {output_file}")
        except Exception as e:
            logger.error(f"Error writing to Excel file: {e}")
            logger.debug(traceback.format_exc())
            raise
        
        logger.info(f"\nResults saved to: {output_file}")
        logger.info(f"Total primer pairs: {len(df)}")
        
        return output_file
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