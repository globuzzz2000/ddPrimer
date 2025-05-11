#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Alignment workflow for ddPrimer pipeline.

This module implements the cross-species alignment workflow:
1. Run alignment first (or use pre-computed MAF)
2. Identify conserved regions
3. Design primers on alignment-masked sequences
4. Optionally filter primers based on SNP/indel positions from both genomes (if --snp is enabled)
"""

import os
import logging
import tempfile
import shutil

from ..config import Config
from ..config.exceptions import AlignmentError, FileFormatError
from ..utils.file_io import FileIO
from ..core import SNPFilterProcessor
from .lastz_runner import LastZRunner
from .maf_parser import MAFParser

# Set up logger
logger = logging.getLogger("ddPrimer.helpers")


def run_alignment_workflow(args, output_dir):
    """
    Run the alignment primer design workflow with improved approach:
    1. Run alignment first (or use pre-computed MAF)
    2. Identify conserved regions
    3. Create alignment-masked sequences (without SNP masking)
    4. Extract variants for later post-design filtering (if --snp is enabled)
    
    Args:
        args (argparse.Namespace): Command line arguments
        output_dir (str): Output directory path
    
    Returns:
        tuple: (alignment_masked_sequences, coordinate_map, variants) - Masked sequences, 
               coordinate mapping, and variants for post-design filtering
        
    Raises:
        AlignmentError: If there's an error in the alignment workflow
    """
    # Initialize processors
    maf_parser = MAFParser()
    snp_processor = SNPFilterProcessor()
    
    try:
        # Step 1: Handle MAF file (pre-computed or generate new one)
        maf_file, ref_fasta_required, second_fasta_required = _handle_maf_file(args, output_dir)
        
        # Step 2: Parse MAF file and identify conserved regions
        conserved_regions, coordinate_map = _process_maf_file(maf_file, maf_parser, args)
        
        # Step 3: Load reference sequences (either from FASTA or MAF)
        reference_sequences, temp_dir = _load_reference_sequences(
            args, maf_parser, ref_fasta_required, output_dir
        )
        
        # Step 4: Create alignment-masked sequences (mask non-conserved regions only)
        alignment_masked_sequences = _create_alignment_masked_sequences(
            args, maf_parser, reference_sequences, conserved_regions, output_dir
        )
        
        # Step 5: Extract variant positions for post-design filtering if SNP filtering is enabled
        variants = None
        if args.snp:
            # Extract variants from both genomes
            ref_variants, second_variants = _extract_variant_positions(
                args, snp_processor, conserved_regions, coordinate_map
            )
            
            # Combine variants from both genomes into a single variant map for filtering
            variants = ref_variants
            
            # Map variants from second genome to reference coordinates if available
            if second_variants:
                second_mapped_variants = _map_second_variants_to_reference(
                    second_variants, coordinate_map
                )
                
                # Merge with reference variants
                for chrom, positions in second_mapped_variants.items():
                    if chrom in variants:
                        variants[chrom].update(positions)
                    else:
                        variants[chrom] = positions
        
        # Clean up intermediate files if not in debug mode
        if not Config.DEBUG_MODE:
            _clean_up_intermediate_files(args, output_dir, temp_dir)
        
        return alignment_masked_sequences, coordinate_map, variants
        
    except Exception as e:
        logger.error(f"Error in alignment workflow: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise AlignmentError(f"Alignment workflow failed: {str(e)}")


def _handle_maf_file(args, output_dir):
    """
    Handle the MAF file, either using a pre-computed one or generating a new one.
    
    Args:
        args (argparse.Namespace): Command line arguments
        output_dir (str): Output directory path
        
    Returns:
        tuple: (maf_file, ref_fasta_required, second_fasta_required)
        
    Raises:
        AlignmentError: If handling the MAF file fails
    """
    if hasattr(args, 'maf') and args.maf is not None:
        # If the maf argument is True (flag used without value), it will be handled
        # by the alignment_mode.py module which prompts for file selection
        logger.debug("\n>>> Using pre-computed MAF file <<<")
        maf_file = args.maf
        logger.debug(f"Using MAF file: {maf_file}")
        
        # With pre-computed MAF, we don't need the FASTA files
        ref_fasta_required = False
        second_fasta_required = False
        return maf_file, ref_fasta_required, second_fasta_required
    else:
        # Need FASTA files to run alignment
        ref_fasta_required = True
        second_fasta_required = True
        
        # Check if required files are provided
        if not args.fasta:
            raise AlignmentError("Reference genome FASTA (--fasta) is required for alignment")
        if not args.second_fasta:
            raise AlignmentError("Second FASTA (--second-fasta) is required for alignment")
        
        return _run_lastz_alignment(args, output_dir), ref_fasta_required, second_fasta_required


def _run_lastz_alignment(args, output_dir):
    """
    Run LastZ alignment between reference and query genomes.
    
    Args:
        args (argparse.Namespace): Command line arguments
        output_dir (str): Output directory path
        
    Returns:
        str: Path to the output MAF file
        
    Raises:
        AlignmentError: If LastZ alignment fails
    """
    logger.debug("Running LastZ alignment between genomes...")

    # Create alignment directory
    input_dir = os.path.dirname(os.path.abspath(args.fasta))
    alignments_dir = os.path.join(input_dir, "Alignments")
    os.makedirs(alignments_dir, exist_ok=True)
    
    # Ensure lastz_options include MAF format
    lastz_options = Config.LASTZ_OPTIONS
    if "--format=maf" not in lastz_options and "-f maf" not in lastz_options:
        lastz_options += " --format=maf"
    
    # Run LastZ alignment
    try:
        lastz_runner = LastZRunner()
        maf_file = lastz_runner.run_parallel_alignment(
            args.fasta,
            args.second_fasta,
            alignments_dir,
            lastz_options
        )
        logger.info(f"LastZ alignment completed: {maf_file}")
        return maf_file
    except Exception as e:
        logger.error(f"Error running LastZ alignment: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise AlignmentError(f"LastZ alignment failed: {str(e)}")


def _process_maf_file(maf_file, maf_parser, args):
    """
    Parse MAF file and identify conserved regions.
    
    Args:
        maf_file (str): Path to MAF file
        maf_parser (MAFParser): MAF parser instance
        args (argparse.Namespace): Command line arguments
        
    Returns:
        tuple: (conserved_regions, coordinate_map)
        
    Raises:
        AlignmentError: If processing the MAF file fails
    """
    logger.debug("\nParsing alignment and identifying conserved regions...")
    try:
        # First analyze the MAF file to understand its structure
        maf_analysis = maf_parser.analyze_maf_file(maf_file)
        logger.debug(f"MAF file contains {maf_analysis['alignment_count']} alignment blocks")
        logger.debug(f"Reference sequences: {', '.join(maf_analysis['ref_seq_ids'][:5])}...")
        logger.debug(f"Query sequences: {', '.join(maf_analysis['query_seq_ids'][:5])}...")
        
        # Parse MAF file
        alignments = maf_parser.parse_maf_file(maf_file)
        
        # Identify conserved regions - use Config instead of args
        conserved_regions = maf_parser.identify_conserved_regions(
            Config.MIN_IDENTITY,
            Config.MIN_SEGMENT_LENGTH
        )
        
        # Generate coordinate mapping between reference and query genomes
        coordinate_map = maf_parser.generate_coordinate_map(conserved_regions)
        logger.debug("Generated coordinate mapping between reference and second genome")
        
        total_regions = sum(len(regions) for regions in conserved_regions.values())
        logger.info(f"Identified {total_regions} conserved regions across {len(conserved_regions)} chromosomes")
        
        return conserved_regions, coordinate_map
    except Exception as e:
        logger.error(f"Error processing MAF file: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise AlignmentError(f"MAF file processing failed: {str(e)}")


def _extract_variant_positions(args, snp_processor, conserved_regions, coordinate_map):
    """
    Extract variant positions from both genomes.
    
    Args:
        args (argparse.Namespace): Command line arguments
        snp_processor (SNPMaskingProcessor): SNP masking processor instance
        conserved_regions (dict): Conserved regions with mapping information
        coordinate_map (dict): Coordinate mapping from reference to query
        
    Returns:
        tuple: (ref_variants, second_variants)
    """
    ref_variants = {}
    second_variants = {}
    
    # Get reference genome variants if VCF provided and SNP filtering is enabled
    if args.snp and args.vcf:
        logger.info("Extracting variants from reference genome VCF...")
        try:
            # Use region-specific variant extraction
            ref_variants = snp_processor.extract_variants_by_regions(
                args.vcf, conserved_regions
            )
            total_ref_variants = sum(len(positions) for positions in ref_variants.values())
            logger.debug(f"Extracted {total_ref_variants} variants from reference genome for conserved regions")
        except Exception as e:
            logger.error(f"Error extracting reference variants: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.warning("Continuing without reference genome variants")

    # Get second genome variants if VCF provided and SNP filtering is enabled
    if args.snp and args.second_vcf:
        logger.info("Extracting variants from second genome VCF...")
        try:
            # Extract variants from the second genome
            # Create a mapping of conserved regions in the second genome's coordinates
            second_genome_regions = _map_regions_to_second_genome(conserved_regions)
            
            # Extract variants for these regions in the second genome
            second_variants = snp_processor.extract_variants_by_regions(
                args.second_vcf, second_genome_regions
            )
            total_second_variants = sum(len(positions) for positions in second_variants.values())
            logger.debug(f"Extracted {total_second_variants} variants from second genome for conserved regions")
        except Exception as e:
            logger.error(f"Error extracting second genome variants: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.warning("Continuing without second genome variants")
    
    return ref_variants, second_variants


def _map_regions_to_second_genome(conserved_regions):
    """
    Create a mapping of conserved regions in the second genome's coordinates.
    
    Args:
        conserved_regions (dict): Conserved regions with mapping information
        
    Returns:
        dict: Conserved regions mapped to second genome coordinates
    """
    second_genome_regions = {}
    
    for chrom, regions in conserved_regions.items():
        for region in regions:
            qry_src = region['qry_src']
            if qry_src not in second_genome_regions:
                second_genome_regions[qry_src] = []
            
            # Add the region in second genome coordinates
            second_genome_regions[qry_src].append({
                'start': region['qry_start'],
                'end': region['qry_end']
            })
    
    return second_genome_regions


def _load_reference_sequences(args, maf_parser, ref_fasta_required, output_dir):
    """
    Load reference sequences either from FASTA file or from MAF file.
    
    Args:
        args (argparse.Namespace): Command line arguments
        maf_parser (MAFParser): MAF parser instance
        ref_fasta_required (bool): Whether reference FASTA is required
        output_dir (str): Output directory path
        
    Returns:
        tuple: (reference_sequences, temp_dir)
        
    Raises:
        FileFormatError: If loading the reference sequences fails
    """
    reference_sequences = {}
    temp_dir = None
    
    # If we have a pre-computed MAF and don't need to run alignment
    if not ref_fasta_required:
        # We need to generate a minimal reference FASTA from the MAF file
        # containing only the conserved regions
        logger.info("Generating reference sequences from MAF file...")
        try:
            reference_sequences = maf_parser.extract_reference_sequences_from_maf()
            logger.debug(f"Generated {len(reference_sequences)} reference sequences from MAF file")
            
            # Create a temporary file in the temp_dir inside output_dir
            temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
            temp_ref_fasta = os.path.join(temp_dir, "ref_from_maf.fasta")
            
            with open(temp_ref_fasta, 'w') as f:
                for seq_id, seq in reference_sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")
            
            logger.debug(f"Wrote reference sequences to temporary file: {temp_ref_fasta}")
            
        except Exception as e:
            logger.error(f"Error generating reference sequences from MAF: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(f"Failed to generate reference sequences from MAF: {str(e)}")
    else:
        # Load reference sequences from provided FASTA
        logger.debug("\n>>> Loading reference sequences from FASTA file <<<")
        try:
            reference_sequences = FileIO.load_fasta(args.fasta)
            logger.debug(f"Loaded {len(reference_sequences)} sequences from reference FASTA")
        except Exception as e:
            logger.error(f"Error loading reference FASTA: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(f"Failed to load reference FASTA: {str(e)}")
    
    return reference_sequences, temp_dir


def _create_alignment_masked_sequences(args, maf_parser, reference_sequences, conserved_regions, output_dir):
    """
    Create alignment-masked sequences (mask non-conserved regions only).
    
    Args:
        args (argparse.Namespace): Command line arguments
        maf_parser (MAFParser): MAF parser instance
        reference_sequences (dict): Dictionary of reference sequences
        conserved_regions (dict): Conserved regions with mapping information
        output_dir (str): Output directory path
        
    Returns:
        dict: Alignment-masked sequences
        
    Raises:
        AlignmentError: If creating alignment-masked sequences fails
    """
    logger.debug("\n>>> Creating alignment-masked genome (conserved regions only) <<<")

    # Create a masked reference with only conserved regions
    masked_fasta_path = os.path.join(output_dir, "alignment_masked_reference.fasta")
    try:
        maf_parser.mask_non_conserved_regions(
            reference_sequences,
            masked_fasta_path,
            conserved_regions,
            Config.MIN_IDENTITY
        )
        logger.debug(f"Created alignment-masked reference genome: {masked_fasta_path}")
    except Exception as e:
        logger.error(f"Error masking non-conserved regions: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise AlignmentError(f"Failed to mask non-conserved regions: {str(e)}")

    # Load the alignment-masked sequences
    try:
        alignment_masked_sequences = FileIO.load_fasta(masked_fasta_path)
        logger.debug(f"Loaded {len(alignment_masked_sequences)} alignment-masked sequences")
        logger.debug(f"Ready for primer design on {len(alignment_masked_sequences)} alignment-masked sequences")
        return alignment_masked_sequences
    except Exception as e:
        logger.error(f"Error loading alignment-masked sequences: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise FileFormatError(f"Failed to load alignment-masked sequences: {str(e)}")


def _map_second_variants_to_reference(second_variants, coordinate_map):
    """
    Map variants from second genome to reference coordinates.
    
    Args:
        second_variants (dict): Dictionary of second genome variants
        coordinate_map (dict): Coordinate mapping between genomes
        
    Returns:
        dict: Mapped variants in reference genome coordinates
    """
    mapped_variants = {}
    
    # Process each sequence
    for seq_id, sequence_mappings in coordinate_map.items():
        mapped_variants[seq_id] = set()
        
        # Skip if no coordinate map for this sequence
        if not sequence_mappings:
            logger.debug(f"No coordinate mapping found for {seq_id}")
            continue
            
        # Process each variant from second genome
        for second_chrom, positions in second_variants.items():
            # Find mappings from second genome to reference
            for pos in positions:
                for ref_pos, mapping in sequence_mappings.items():
                    if mapping["qry_src"] == second_chrom and mapping["qry_pos"] == pos:
                        # Found a mapping from second genome variant to reference
                        mapped_variants[seq_id].add(ref_pos)
    
    # Count mapped variants
    total_mapped = sum(len(positions) for positions in mapped_variants.values())
    logger.debug(f"Mapped {total_mapped} variants from second genome to reference coordinates")
    
    return mapped_variants


def _clean_up_intermediate_files(args, output_dir, temp_dir):
    """
    Clean up intermediate files created during the workflow.
    
    Args:
        args (argparse.Namespace): Command line arguments
        output_dir (str): Output directory path
        temp_dir (str): Temporary directory path
    """
    # Check debug mode from both args and Config (for backward compatibility)
    is_debug_mode = args.debug if hasattr(args, 'debug') else False
    is_debug_mode = is_debug_mode or Config.DEBUG_MODE
    
    # Skip cleanup if in debug mode
    if is_debug_mode:
        logger.debug("Debug mode enabled - skipping cleanup of intermediate files")
        return
        
    try:
        # List of files to clean up
        files_to_clean = [
            os.path.join(output_dir, "alignment_masked_reference.fasta")
        ]
        
        # Clean up temp directory if it was created for ref_from_maf.fasta
        if temp_dir and os.path.exists(temp_dir):
            logger.debug(f"Cleaning up temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir)
        
        for file_path in files_to_clean:
            if os.path.exists(file_path):
                logger.debug(f"Cleaning up intermediate file: {file_path}")
                os.remove(file_path)
    except Exception as e:
        logger.warning(f"Error during cleanup of masked files: {str(e)}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)