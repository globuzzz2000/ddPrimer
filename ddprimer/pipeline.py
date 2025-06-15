#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified ddPrimer Pipeline

A streamlined primer design pipeline that supports two workflows:
1. Standard Mode: FASTA + VCF + GFF files for genome-based primer design
2. Direct Mode: CSV/Excel files with target sequences

Features:
- SNP masking and substitution for both modes
- Restriction site cutting and gene annotation filtering
- Primer3 design with thermodynamic analysis
- BLAST specificity checking
- Comprehensive Excel output

This unified module replaces the previous multi-mode architecture
for improved maintainability and reduced complexity.
"""

import sys
import os
import argparse
import logging
import traceback
import pandas as pd
from pathlib import Path
from typing import Optional, Union
from tqdm import tqdm

# Import package modules
from .config import (Config, setup_logging, display_config, display_primer3_settings, 
                    DDPrimerError, FileError, ExternalToolError)
from .utils import FileIO, BlastDBCreator, BlastVerification, DirectModeUtils
from .core import (SNPMaskingProcessor, AnnotationProcessor, PrimerProcessor, 
                   BlastProcessor, SequenceProcessor, Primer3Processor, ViennaRNAProcessor)

# Type alias for path inputs
PathLike = Union[str, Path]

# Set up module logger
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse command line arguments for unified pipeline."""
    parser = argparse.ArgumentParser(
        description='ddPrimer: Unified pipeline for primer design and filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='ddprimer [--direct [.csv, .xlsx]] [--debug [MODULE...]] [--config [.json]] \n'
              '        [--cli] [--nooligo] [--snp] [--noannotation] [--db ...] \n'
              '        [--fasta ...] [--vcf ...] [--gff ...] [--output ...]'
    )

    # Modes
    parser.add_argument('--direct', metavar='[.csv, .xlsx]', nargs='?', const=True, 
                       help='Direct mode: design primers from target sequences in CSV/Excel')

    # Options
    parser.add_argument('--debug', nargs='*', metavar='MODULE', 
                       help='Enable debug mode (universal or specific modules)')
    parser.add_argument('--output', metavar='<dir>', help='Output directory')
    parser.add_argument('--config', metavar='[.json]', nargs='?', const='DISPLAY', 
                       help='Configuration file or display mode')
    parser.add_argument('--cli', action='store_true', help='Force CLI mode')
    parser.add_argument('--nooligo', action='store_true', help='Disable probe design')
    parser.add_argument('--snp', action='store_true', help='Enable SNP masking')
    parser.add_argument('--noannotation', action='store_true', help='Skip gene filtering')
    parser.add_argument('--db', nargs='*', metavar=('[.fasta]', '[name]'),
                       help='Create/select BLAST database')

    # Input files
    parser.add_argument('--fasta', metavar='[.fasta]', help='Reference genome FASTA')
    parser.add_argument('--vcf', metavar='[.vcf]', help='Variant file')
    parser.add_argument('--gff', metavar='[.gff]', help='Annotation file')
    
    args = parser.parse_args()

    # Process debug argument
    if args.debug is not None:
        args.debug = True if len(args.debug) == 0 else args.debug
    else:
        args.debug = False

    # Validate SNP requirements
    if args.snp and args.cli and (not args.fasta or not args.vcf):
        parser.error("SNP masking requires --fasta and --vcf")

    # Process database arguments
    if args.db is not None:
        args.dboutdir = args.output
        if len(args.db) == 0:
            args.db_action, args.db_fasta, args.db_name = 'select', None, None
        else:
            args.db_action = 'file'
            args.db_fasta = args.db[0]
            args.db_name = args.db[1] if len(args.db) > 1 else None
    else:
        args.db_action = args.db_fasta = args.db_name = None

    return args


def setup_pipeline(args):
    """Set up logging, configuration, and database."""
    # Setup logging
    setup_logging(debug=args.debug)
    
    logger.debug("Initializing unified ddPrimer pipeline")
    Config.get_instance()
    
    # Handle configuration display/loading
    if args.config in ['DISPLAY', 'all', 'basic', 'template']:
        if args.config == 'template':
            from .config import generate_config_template
            generate_config_template(Config, output_dir=args.output)
        else:
            display_config(Config)
            if args.config == 'all':
                display_primer3_settings(Config)
        return False  # Exit after display
    
    # Load custom configuration
    if args.config and args.config not in ['DISPLAY', 'all', 'basic', 'template']:
        try:
            Config.load_from_file(args.config)
            logger.debug(f"Loaded configuration from {args.config}")
        except Exception as e:
            raise FileError(f"Error loading configuration: {e}")
    
    # Apply CLI and probe settings
    if args.cli:
        FileIO.use_cli = True
        logger.debug("CLI mode enabled")
    
    if args.nooligo:
        Config.PRIMER3_SETTINGS["PRIMER_PICK_INTERNAL_OLIGO"] = 0
        Config.DISABLE_INTERNAL_OLIGO = True
        logger.info("Probe design disabled")
    
    # Handle database operations
    if args.db is not None:
        handle_database_operations(args)
        if not args.fasta and not args.direct:
            return False  # Exit after database-only operations
    
    # Verify BLAST database
    if not BlastVerification.verify_blast_database():
        raise ExternalToolError("BLAST database verification failed", tool_name="blastn")
    
    return True


def handle_database_operations(args):
    """Handle database creation and selection."""
    logger.debug("Processing BLAST database operations")
    
    blast_db_creator = BlastDBCreator()
    
    if args.db_action == 'select':
        # Model organism selection
        from .utils import ModelOrganismManager as MOManager
        organism_key, organism_name, fasta_file = MOManager.select_model_organism()
        
        if organism_key is None and fasta_file is None:
            logger.info("Database selection canceled")
            return
        
        if organism_key == 'existing_db':
            # Select existing database
            Config.DB_PATH = fasta_file
            Config.save_database_config(fasta_file)
            Config.USE_CUSTOM_DB = True
            logger.info(f"Selected database: {fasta_file}")
            return
        
        # Create from model organism
        if organism_key and not args.db_name:
            scientific_name = organism_name.split(' (')[0] if ' (' in organism_name else organism_name
            db_name = scientific_name.replace(' ', '_')
        else:
            db_name = args.db_name
    
    elif args.db_action == 'file':
        # Create from specified file
        fasta_file = args.db_fasta
        if not Path(fasta_file).exists():
            raise FileError(f"FASTA file not found: {fasta_file}")
        organism_key = None
        db_name = args.db_name
    
    # Create database
    if 'fasta_file' in locals() and organism_key != 'existing_db':
        db_path = blast_db_creator.create_db(fasta_file, db_name, args.output)
        
        # Clean up model organism file
        if organism_key:
            from .utils import ModelOrganismManager
            ModelOrganismManager.cleanup_genome_file(fasta_file, logger)
        
        # Update configuration
        if (not Config.DB_PATH or 
            Config.DB_PATH == "/Library/Application Support/Blast_DBs/Tair DB/TAIR10"):
            Config.DB_PATH = db_path
            Config.save_database_config(db_path)
            Config.USE_CUSTOM_DB = True
            logger.info(f"BLAST database created: {db_path}")
        else:
            # Ask user if they want to use new database
            use_new = input(f"Current database: {Config.DB_PATH}\n"
                          f"Use new database instead? [Y/n]: ").strip().lower()
            if use_new == "" or use_new.startswith("y"):
                Config.DB_PATH = db_path
                Config.save_database_config(db_path)
                Config.USE_CUSTOM_DB = True
                logger.info(f"Now using: {db_path}")


def load_sequences_and_variants(args):
    """Load sequences and variants based on mode."""
    if args.direct:
        return load_direct_mode_data(args)
    else:
        return load_standard_mode_data(args)


def load_direct_mode_data(args):
    """Load data for direct mode."""
    logger.info("=== Direct Mode: Target Sequence Design ===")
    
    # Get sequence file
    if args.direct is True:
        sequence_file = FileIO.select_sequences_file()
    else:
        sequence_file = args.direct
    
    # Set up output directory
    output_dir = setup_output_directory(args, sequence_file)
    
    # Analyze file structure
    analysis = DirectModeUtils.analyze_sequence_file(sequence_file)
    DirectModeUtils.print_analysis(analysis)
    
    # Initialize tracking variables
    matching_status = {}
    all_sequences = {}
    
    if args.snp:
        # SNP masking workflow for direct mode
        logger.debug("SNP masking enabled for direct mode")
        
        # Get reference files
        ref_fasta = args.fasta or FileIO.select_fasta_file("Select reference FASTA")
        ref_vcf = args.vcf or FileIO.select_file("Select VCF file", 
                                               [("VCF Files", "*.vcf"), ("VCF.gz Files", "*.vcf.gz")])
        
        FileIO.mark_selection_complete()
        
        # Load target sequences
        sequences = FileIO.load_sequences_from_table(sequence_file)
        all_sequences = sequences.copy()
        
        # Apply SNP masking with location finding
        masked_sequences = apply_direct_snp_masking(sequences, ref_fasta, ref_vcf, 
                                                  matching_status, all_sequences)
    else:
        # No SNP masking
        logger.debug("SNP masking disabled")
        FileIO.mark_selection_complete()
        
        sequences = FileIO.load_sequences_from_table(sequence_file)
        all_sequences = sequences.copy()
        masked_sequences = sequences.copy()
        
        # Mark all as not attempted
        for seq_id in sequences:
            matching_status[seq_id] = "Not attempted"
    
    return {
        'sequences': masked_sequences,
        'genes': None,  # No gene filtering in direct mode
        'output_dir': output_dir,
        'reference_file': sequence_file,
        'mode': 'direct',
        'matching_status': matching_status,
        'all_sequences': all_sequences
    }


def load_standard_mode_data(args):
    """Load data for standard mode."""
    logger.info("=== Standard Mode: Genome-Based Design ===")
    
    # Get input files
    fasta_file = args.fasta or FileIO.select_fasta_file("Select reference FASTA")
    vcf_file = args.vcf or FileIO.select_file("Select VCF file", 
                                            [("VCF Files", "*.vcf"), ("VCF.gz Files", "*.vcf.gz")])
    
    gff_file = None
    if not args.noannotation:
        gff_file = args.gff or FileIO.select_file("Select GFF annotation file",
                                                [("GFF Files", "*.gff"), ("GFF3 Files", "*.gff3")])
    
    FileIO.mark_selection_complete()
    
    # Set up output directory
    output_dir = setup_output_directory(args, fasta_file)
    
    # Process variants and sequences
    logger.info("Extracting variants from VCF...")
    snp_processor = SNPMaskingProcessor()
    variants = snp_processor.process_vcf_with_chromosome_mapping(
        vcf_file, fasta_file,
        min_af=Config.SNP_ALLELE_FREQUENCY_THRESHOLD,
        min_qual=Config.SNP_QUALITY_THRESHOLD
    )
    
    # Count variants
    total_variants = sum(len(v_list) for v_list in variants.values())
    fixed_variants = sum(sum(1 for v in v_list if v.is_fixed) for v_list in variants.values())
    logger.info(f"Extracted {total_variants} variants ({fixed_variants} fixed)")
    
    # Load sequences
    logger.info("Loading sequences from FASTA...")
    sequences = FileIO.load_fasta(fasta_file)
    
    # Apply SNP processing
    logger.info("Processing variants in sequences...")
    masked_sequences = snp_processor.mask_sequences_for_primer_design(
        sequences=sequences,
        variants=variants,
        flanking_size=Config.SNP_FLANKING_MASK_SIZE,
        use_soft_masking=Config.SNP_USE_SOFT_MASKING,
        min_af=Config.SNP_ALLELE_FREQUENCY_THRESHOLD,
        min_qual=Config.SNP_QUALITY_THRESHOLD
    )
    
    # Load gene annotations
    genes = None
    if not args.noannotation:
        logger.info("Loading gene annotations...")
        genes = AnnotationProcessor.load_genes_from_gff(gff_file)
        logger.info(f"Loaded {len(genes)} gene annotations")
    
    return {
        'sequences': masked_sequences,
        'genes': genes,
        'output_dir': output_dir,
        'reference_file': fasta_file,
        'mode': 'standard',
        'matching_status': None,
        'all_sequences': None
    }


def apply_direct_snp_masking(sequences, ref_fasta, ref_vcf, matching_status, all_sequences):
    """Apply SNP masking for direct mode with location finding."""
    logger.info("Applying SNP masking with reference genome location finding...")
    
    snp_processor = SNPMaskingProcessor()
    
    # Find genomic locations
    sequence_locations = {}
    for seq_id, sequence in sequences.items():
        logger.debug(f"Finding location for sequence {seq_id}")
        
        source_chrom, start_pos, end_pos, identity = DirectModeUtils.find_sequence_location(
            sequence, ref_fasta
        )
        
        if source_chrom:
            matching_status[seq_id] = "Success"
            sequence_locations[seq_id] = {
                'chrom': source_chrom,
                'start': start_pos,
                'end': end_pos,
                'identity': identity
            }
            logger.debug(f"Located {seq_id} at {source_chrom}:{start_pos}-{end_pos}")
        else:
            matching_status[seq_id] = "Failure"
            logger.debug(f"Could not locate {seq_id} in reference genome")
    
    # Get successfully located sequences
    located_sequences = {seq_id: seq for seq_id, seq in sequences.items() 
                        if seq_id in sequence_locations}
    
    if not located_sequences:
        logger.warning("No sequences could be located in reference genome")
        return sequences
    
    # Extract all variants from VCF
    all_variants = snp_processor.process_vcf_with_chromosome_mapping(
        ref_vcf, ref_fasta,
        min_af=Config.SNP_ALLELE_FREQUENCY_THRESHOLD,
        min_qual=Config.SNP_QUALITY_THRESHOLD
    )
    
    # Apply variants to located sequences
    from .core.snp_processor import Variant
    
    masked_sequences = {}
    for seq_id, sequence in sequences.items():
        if seq_id not in sequence_locations:
            # Keep original sequence if not located
            masked_sequences[seq_id] = sequence
            continue
        
        location = sequence_locations[seq_id]
        chrom = location['chrom']
        start = location['start']
        
        # Get variants for this chromosome and adjust positions
        seq_variants = []
        if chrom in all_variants:
            for variant in all_variants[chrom]:
                if start <= variant.position <= location['end']:
                    adjusted_position = variant.position - start + 1
                    adjusted_variant = Variant(
                        position=adjusted_position,
                        ref=variant.ref,
                        alt=variant.alt,
                        qual=variant.qual,
                        af=variant.af,
                        is_fixed=variant.is_fixed
                    )
                    seq_variants.append(adjusted_variant)
        
        # Apply masking
        if seq_variants:
            masked_sequence = snp_processor.mask_and_substitute_variants(
                sequence, seq_variants,
                flanking_size=Config.SNP_FLANKING_MASK_SIZE,
                use_soft_masking=Config.SNP_USE_SOFT_MASKING
            )
            masked_sequences[seq_id] = masked_sequence
        else:
            masked_sequences[seq_id] = sequence
    
    logger.info(f"Applied SNP masking to {len(located_sequences)} located sequences")
    return masked_sequences


def setup_output_directory(args, reference_file):
    """Set up output directory based on arguments and reference file."""
    if args.output:
        output_dir = args.output
    else:
        input_dir = os.path.dirname(os.path.abspath(reference_file))
        output_dir = os.path.join(input_dir, "Primers")
    
    os.makedirs(output_dir, exist_ok=True)
    logger.debug(f"Output directory: {output_dir}")
    return output_dir


def process_restriction_sites(sequences):
    """Cut sequences at restriction sites."""
    logger.info("Processing restriction sites...")
    logger.debug(f"Using restriction site pattern: {Config.RESTRICTION_SITE}")
    
    fragments = SequenceProcessor.cut_at_restriction_sites(sequences)
    logger.info(f"Generated {len(fragments)} restriction fragments")
    return fragments


def filter_fragments(fragments, genes, mode, skip_annotation_filtering):
    """Filter fragments based on mode and gene overlap."""
    if mode == 'direct' or skip_annotation_filtering:
        # Direct mode or no annotation filtering
        logger.info("Preparing fragments for direct processing...")
        filtered = []
        for fragment in fragments:
            base_id = fragment["id"].split("_frag")[0]
            simplified = {
                "id": fragment["id"],
                "sequence": fragment["sequence"],
                "Gene": base_id
            }
            filtered.append(simplified)
        logger.info(f"Prepared {len(filtered)} fragments")
        return filtered
    else:
        # Standard mode with gene filtering
        if not genes:
            logger.error("Gene annotations required for standard mode")
            logger.info("Use --noannotation to skip gene filtering")
            return []
        
        logger.info("Extracting gene-overlapping regions...")
        filtered = AnnotationProcessor.filter_by_gene_overlap_enhanced(fragments, genes)
        logger.info(f"Extracted {len(filtered)} gene fragments")
        return filtered


def design_primers(fragments):
    """Design primers using Primer3."""
    logger.info("Designing primers with Primer3...")
    
    primer3_processor = Primer3Processor(Config)
    
    # Prepare Primer3 inputs
    primer3_inputs = []
    fragment_info = {}
    
    for fragment in fragments:
        primer3_inputs.append({
            "SEQUENCE_ID": fragment["id"],
            "SEQUENCE_TEMPLATE": fragment["sequence"]
        })
        
        # Store fragment metadata
        if "chr" in fragment and "start" in fragment:
            # Standard mode with coordinates
            fragment_info[fragment["id"]] = {
                "chr": fragment.get("chr", ""),
                "start": fragment.get("start", 1),
                "end": fragment.get("end", len(fragment["sequence"])),
                "gene": fragment.get("Gene", fragment["id"])
            }
        else:
            # Direct mode without coordinates
            fragment_info[fragment["id"]] = {
                "gene": fragment.get("Gene", fragment["id"])
            }
    
    # Run Primer3
    primer3_output = primer3_processor.run_primer3_batch_parallel(primer3_inputs)
    
    # Parse results
    primer_results = primer3_processor.parse_primer3_batch(primer3_output, fragment_info)
    logger.info(f"Designed primers for {len(primer_results)} fragments")
    
    return primer_results


def filter_and_analyze_primers(primer_results):
    """Filter primers and calculate properties."""
    logger.info("Filtering and analyzing primers...")
    
    # Convert to DataFrame
    df = pd.DataFrame(primer_results)
    initial_count = len(df)
    
    if initial_count == 0:
        logger.warning("No primer results to process")
        return df
    
    # Apply filters
    logger.debug("Applying penalty filter...")
    df = PrimerProcessor.filter_by_penalty(df)
    
    logger.debug("Applying repeat filter...")
    df = PrimerProcessor.filter_by_repeats(df)
    
    logger.debug("Applying GC content filter...")
    df = PrimerProcessor.filter_by_gc_content(df)
    
    if len(df) == 0:
        logger.warning("No primers passed filtering")
        return df
    
    logger.info(f"Filtered primers: {len(df)}/{initial_count} passed")
    
    # Process internal oligos
    df = PrimerProcessor.process_internal_oligos(df)
    
    # Calculate thermodynamic properties
    logger.info("Calculating thermodynamic properties...")
    df = calculate_thermodynamics(df)
    
    # Run BLAST analysis
    logger.info("Running BLAST specificity analysis...")
    df = run_blast_analysis(df)
    
    return df


def calculate_thermodynamics(df):
    """Calculate thermodynamic properties using ViennaRNA."""
    # Add empty columns
    df["Primer F dG"] = None
    df["Primer R dG"] = None
    if "Probe" in df.columns:
        df["Probe dG"] = None
    df["Amplicon dG"] = None
    
    # Calculate ΔG values
    if Config.SHOW_PROGRESS:
        tqdm.pandas(desc="Forward primers")
        df["Primer F dG"] = df["Primer F"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        
        tqdm.pandas(desc="Reverse primers")
        df["Primer R dG"] = df["Primer R"].progress_apply(ViennaRNAProcessor.calc_deltaG)
        
        if "Probe" in df.columns:
            tqdm.pandas(desc="Probes")
            df["Probe dG"] = df["Probe"].progress_apply(
                lambda x: ViennaRNAProcessor.calc_deltaG(x) if pd.notnull(x) and x else None
            )
        
        tqdm.pandas(desc="Amplicons")
        df["Amplicon dG"] = df["Amplicon"].progress_apply(ViennaRNAProcessor.calc_deltaG)
    else:
        df["Primer F dG"] = df["Primer F"].apply(ViennaRNAProcessor.calc_deltaG)
        df["Primer R dG"] = df["Primer R"].apply(ViennaRNAProcessor.calc_deltaG)
        
        if "Probe" in df.columns:
            df["Probe dG"] = df["Probe"].apply(
                lambda x: ViennaRNAProcessor.calc_deltaG(x) if pd.notnull(x) and x else None
            )
        
        df["Amplicon dG"] = df["Amplicon"].apply(ViennaRNAProcessor.calc_deltaG)
    
    return df


def run_blast_analysis(df):
    """Run BLAST analysis for primer specificity."""
    # BLAST forward primers
    blast_results_f = []
    primers_f = df["Primer F"].tolist()
    if Config.SHOW_PROGRESS:
        primers_f = tqdm(primers_f, desc="BLAST forward primers")
    
    for primer_f in primers_f:
        blast1, blast2 = BlastProcessor.blast_short_seq(primer_f)
        blast_results_f.append((blast1, blast2))
    
    df["Primer F BLAST1"], df["Primer F BLAST2"] = zip(*blast_results_f)
    
    # BLAST reverse primers
    blast_results_r = []
    primers_r = df["Primer R"].tolist()
    if Config.SHOW_PROGRESS:
        primers_r = tqdm(primers_r, desc="BLAST reverse primers")
    
    for primer_r in primers_r:
        blast1, blast2 = BlastProcessor.blast_short_seq(primer_r)
        blast_results_r.append((blast1, blast2))
    
    df["Primer R BLAST1"], df["Primer R BLAST2"] = zip(*blast_results_r)
    
    # BLAST probes if present
    if "Probe" in df.columns:
        blast_results_p = []
        probes = df["Probe"].tolist()
        if Config.SHOW_PROGRESS:
            probes = tqdm(probes, desc="BLAST probes")
        
        for probe in probes:
            if pd.notnull(probe) and probe:
                blast1, blast2 = BlastProcessor.blast_short_seq(probe)
            else:
                blast1, blast2 = None, None
            blast_results_p.append((blast1, blast2))
        
        df["Probe BLAST1"], df["Probe BLAST2"] = zip(*blast_results_p)
    
    # Filter by BLAST specificity
    initial_count = len(df)
    df = PrimerProcessor.filter_by_blast(df)
    logger.info(f"BLAST filtering: {len(df)}/{initial_count} primers passed")
    
    return df


def save_results(df, data_info):
    """Save results to Excel file."""
    if len(df) == 0:
        logger.warning("No results to save")
        return None
    
    # Add missing sequences for direct mode
    if data_info['mode'] == 'direct' and data_info['matching_status'] and data_info['all_sequences']:
        df = DirectModeUtils.add_missing_sequences(
            df, data_info['all_sequences'], data_info['matching_status']
        )
    
    # Save to Excel
    output_path = FileIO.save_results(
        df, 
        data_info['output_dir'], 
        data_info['reference_file'], 
        mode=data_info['mode']
    )
    
    if output_path:
        logger.info(f"Results saved to: {output_path}")
        return output_path
    else:
        logger.error("Failed to save results")
        return None


def run_unified_workflow(args):
    """Main unified workflow for both standard and direct modes."""
    try:
        # Load sequences and variants
        data_info = load_sequences_and_variants(args)
        
        if not data_info['sequences']:
            logger.error("No sequences loaded")
            return False
        
        # Process restriction sites
        fragments = process_restriction_sites(data_info['sequences'])
        if not fragments:
            logger.error("No valid restriction fragments generated")
            return False
        
        # Filter fragments
        filtered_fragments = filter_fragments(
            fragments, 
            data_info['genes'], 
            data_info['mode'], 
            args.noannotation
        )
        
        if not filtered_fragments:
            logger.error("No fragments passed filtering")
            if data_info['mode'] == 'standard' and not args.noannotation:
                logger.info("Try using --noannotation to skip gene filtering")
            return False
        
        # Design primers
        primer_results = design_primers(filtered_fragments)
        if not primer_results:
            logger.error("No primers designed")
            return False
        
        # Filter and analyze primers
        df = filter_and_analyze_primers(primer_results)
        if len(df) == 0:
            logger.error("No primers passed analysis")
            return False
        
        # Save results
        output_path = save_results(df, data_info)
        return output_path is not None
        
    except Exception as e:
        logger.error(f"Workflow error: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False


def run_pipeline():
    """Main pipeline entry point."""
    logger_instance = None
    
    try:
        # Parse arguments and setup
        args = parse_arguments()
        
        # Setup pipeline (returns False if should exit early)
        if not setup_pipeline(args):
            return True
        
        logger_instance = logging.getLogger(__name__)
        logger_instance.info("=== ddPrimer Unified Pipeline ===")
        
        # Run workflow
        success = run_unified_workflow(args)
        
        if success:
            logger_instance.info("=== Pipeline completed successfully! ===")
            return True
        else:
            logger_instance.error("Pipeline execution failed")
            return False
            
    except DDPrimerError as e:
        if logger_instance:
            logger_instance.error(f"Pipeline error: {e}")
        else:
            print(f"Pipeline error: {e}")
        return False
    except Exception as e:
        if logger_instance:
            logger_instance.error(f"Unexpected error: {e}")
            logger_instance.debug(f"Error details: {str(e)}", exc_info=True)
        else:
            print(f"Unexpected error: {e}")
            print(traceback.format_exc())
        return False


def main():
    """Entry point for direct execution."""
    success = run_pipeline()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()