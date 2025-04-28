#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct mode for ddPrimer pipeline.

This module contains the implementation of the direct mode workflow:
1. Load sequences directly from CSV or Excel files
2. Run primer design on these sequences using the common pipeline
3. Optionally mask SNPs in sequences if --snp is enabled
"""

import os
import logging
import pandas as pd
import tempfile
import subprocess
import io
import re

# Import package modules
from ..config import Config
from ..utils import FileUtils
from ..core import SNPMaskingProcessor, BlastProcessor
from . import common

# Set up logging
logger = logging.getLogger("ddPrimer")


def find_sequence_location_with_blast(sequence, ref_fasta, min_identity=90, min_coverage=90):
    """
    Find the location of a sequence in the reference genome using BLAST.
    
    Args:
        sequence: Query sequence to locate
        ref_fasta: Path to reference FASTA file
        min_identity: Minimum percent identity (default: 90)
        min_coverage: Minimum query coverage (default: 90)
        
    Returns:
        tuple: (chromosome, start_position, end_position, percent_identity) or (None, None, None, None) if not found
    """
    # Create temporary file for query sequence
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_file:
        temp_file_path = temp_file.name
        temp_file.write(">query\n")
        temp_file.write(sequence)
    
    try:
        # Run BLAST command
        cmd = [
            'blastn',
            '-query', temp_file_path,
            '-subject', ref_fasta,
            '-outfmt', '6 qseqid sseqid pident qcovs qstart qend sstart send length',
            '-max_target_seqs', '1',
            '-evalue', '1e-10'
        ]
        
        logger.debug(f"Running BLAST command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        blast_output = result.stdout.strip()
        
        # Process BLAST output
        if blast_output:
            fields = blast_output.split('\t')
            if len(fields) >= 9:
                _, chrom, percent_identity, query_coverage, _, _, start, end, _ = fields
                
                # Convert to appropriate types
                percent_identity = float(percent_identity)
                query_coverage = float(query_coverage)
                start = int(start)
                end = int(end)
                
                # Ensure start is always less than end (BLAST might report them reversed)
                if start > end:
                    start, end = end, start
                
                # Check if the alignment meets our criteria
                if percent_identity >= min_identity and query_coverage >= min_coverage:
                    logger.debug(f"Found sequence match on {chrom} at positions {start}-{end}, "
                                f"identity: {percent_identity}%, coverage: {query_coverage}%")
                    return chrom, start, end, percent_identity
        
        logger.debug("No significant matches found in BLAST search")
        return None, None, None, None
    
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST error: {e}")
        logger.debug(f"BLAST stderr: {e.stderr}")
        return None, None, None, None
    
    except Exception as e:
        logger.error(f"Error in BLAST search: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return None, None, None, None
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_file_path):
            os.unlink(temp_file_path)


def add_rows_for_sequences_without_primers(df, masked_sequences, matching_status=None):
    """
    Check for sequences that didn't get primers and add them to the results.
    Also add rows for sequences that failed to match to reference genome.
    """
    logger.debug("Checking for sequences without primers and reference match failures...")
    
    # 1) Which sequences actually got primers?
    sequences_with_primers = set()
    if "Gene" in df.columns:
        valid = df["Primer F"] != "No suitable primers found"
        sequences_with_primers.update(df.loc[valid, "Gene"].astype(str).unique())
    
    # 2) All the input sequence IDs
    all_input = set(masked_sequences.keys())
    
    # 3) Which ones were matched but yielded no primers?
    no_primer = all_input - sequences_with_primers
    
    # 4) Which ones actually *failed* to match?
    failures = set()
    if matching_status:
        failures = {seq for seq, st in matching_status.items() if st == "Failure"}
    
    # 5) Take the overlap out of the no_primer list
    no_primer -= failures
    
    # 6) Build rows
    rows = []
    # a) sequences that matched but got no primers
    for seq in sorted(no_primer):
        rows.append({
            "Gene": seq,
            "Primer F": "No suitable primers found",
            "Reference Match": matching_status.get(seq, "Not attempted") if matching_status else "Not attempted"
        })
    # b) sequences that never matched
    for seq in sorted(failures):
        rows.append({
            "Gene": seq,
            "Primer F": "Sequence could not be matched against reference",
            "Reference Match": "Failure"
        })
    
    # 7) Append them (and ensure DF has a Reference Match column)
    if rows:
        extra = pd.DataFrame(rows)
        if matching_status and "Reference Match" not in df.columns:
            df["Reference Match"] = "Success"
        df = pd.concat([df, extra], ignore_index=True, sort=False)
        logger.debug(f"Added {len(rows)} rows for no-primer / match failures")
    
    return df


def run(args):
    """
    Run the direct mode primer design workflow.
    
    Args:
        args: Command line arguments
        
    Returns:
        bool: Success or failure
    """
    logger.info("=== Direct Mode Workflow ===")
    
    try:
        # Get the input file - handle both when --direct is specified alone or with a file
        if args.direct is True:  # --direct was used without a specified file
            sequence_file = FileUtils.get_sequences_file()
        else:  # --direct was used with a specified file
            sequence_file = args.direct
            
        logger.debug(f"Direct sequence file: {sequence_file}")
        
        # Set up output directory
        if args.output:
            output_dir = args.output
        else:
            # Use the directory of the reference file
            input_dir = os.path.dirname(os.path.abspath(sequence_file))
            output_dir = os.path.join(input_dir, "Primers")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logger.debug(f"Created output directory: {output_dir}")
        
        # Dictionary to track reference matching status
        matching_status = {}
        masked_sequences = {}
        all_sequences = {}  # Store all original sequences, including those that fail matching
        
        # Check if SNP masking is enabled
        if args.snp:
            logger.info("\n>>> SNP masking is enabled <<<")
            
            # Get reference FASTA and VCF files
            ref_fasta = args.fasta
            ref_vcf = args.vcf
            
            # If not provided, prompt for them
            if not ref_fasta:
                logger.info("\n>>> Please select reference FASTA file <<<")
                try:
                    ref_fasta = FileUtils.get_file(
                        "Select reference FASTA file", 
                        [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting reference FASTA file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    args.snp = False
            
            if args.snp and not ref_vcf:
                logger.info("\n>>> Please select VCF file with variants <<<")
                try:
                    ref_vcf = FileUtils.get_file(
                        "Select VCF file with variants", 
                        [("VCF Files", "*.vcf"), ("Compressed VCF Files", "*.vcf.gz"), ("All Files", "*")]
                    )
                except Exception as e:
                    logger.error(f"Error selecting VCF file: {e}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    args.snp = False
        else:
            logger.info("\n>>> SNP masking is disabled <<<")
        
        # Load sequences directly from the provided file
        logger.info("\nLoading sequences from input file...")
        try:
            sequences = FileUtils.load_sequences_from_table(sequence_file)
            logger.debug(f"Loaded {len(sequences)} sequences from {sequence_file}")
            all_sequences = sequences.copy()  # Keep a copy of all sequences
        except Exception as e:
            logger.error(f"Error loading sequences from file: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        if not sequences:
            logger.warning("No sequences found in the input file. Exiting.")
            return False
        
        # Apply SNP masking if enabled
        if args.snp and ref_fasta and ref_vcf:
            logger.info(f"\nMasking SNPs using reference FASTA and VCF...")
            logger.debug(f"Reference FASTA: {ref_fasta}")
            logger.debug(f"Reference VCF: {ref_vcf}")
            
            # Initialize SNP masking processor
            snp_processor = SNPMaskingProcessor()
            
            try:
                # Extract variants from VCF
                variants = snp_processor.get_variant_positions(ref_vcf)
                logger.debug(f"Extracted variants from VCF file")
                
                # Process each sequence
                for seq_id, sequence in sequences.items():
                    logger.debug(f"Processing sequence {seq_id} ({len(sequence)} bp)")
                    
                    # Find sequence location in reference genome using BLAST
                    source_chrom, start_pos, end_pos, identity = find_sequence_location_with_blast(
                        sequence, ref_fasta
                    )
                    
                    # If sequence was successfully matched to reference
                    if source_chrom:
                        matching_status[seq_id] = "Success"
                        logger.debug(f"Matched sequence {seq_id} to {source_chrom}:{start_pos}-{end_pos} "
                                    f"with {identity}% identity")
                        
                        # Apply variants only from the matching chromosome
                        seq_variants = variants.get(source_chrom, set())
                        
                        # Adjust variant positions to sequence coordinates
                        adjusted_variants = set()
                        for var_pos in seq_variants:
                            if start_pos <= var_pos <= end_pos:
                                # Convert genome coordinate to sequence coordinate
                                seq_pos = var_pos - start_pos
                                adjusted_variants.add(seq_pos)
                        
                        # Mask the sequence with adjusted variants
                        if adjusted_variants:
                            logger.debug(f"Masking sequence {seq_id} with {len(adjusted_variants)} variants")
                            masked_sequence = snp_processor.mask_variants(sequence, adjusted_variants)
                            masked_sequences[seq_id] = masked_sequence
                        else:
                            logger.debug(f"No variants found in region, using original sequence")
                            masked_sequences[seq_id] = sequence
                    else:
                        # Sequence could not be matched to reference
                        matching_status[seq_id] = "Failure"
                        logger.debug(f"Could not match sequence {seq_id} to reference genome")
                        # We'll exclude this sequence from masked_sequences
                
                logger.info(f"SNP masking completed")
                logger.info(f"Successfully matched {sum(1 for s in matching_status.values() if s == 'Success')}/{len(sequences)} sequences to reference")
                logger.debug(f"Failed to match {sum(1 for s in matching_status.values() if s == 'Failure')}/{len(sequences)} sequences to reference")
            except Exception as e:
                logger.error(f"Error during SNP masking: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                logger.warning("Using original sequences without SNP masking")
                masked_sequences = sequences
                for seq_id in sequences:
                    matching_status[seq_id] = "Not attempted"
        else:
            # If SNP masking is disabled, use original sequences
            masked_sequences = sequences
            for seq_id in sequences:
                matching_status[seq_id] = "Not attempted"
        
        # Debug logging for sequences
        if Config.DEBUG_MODE:
            logger.debug("\n=== DIRECT MODE SEQUENCE INFO ===")
            for i, (seq_id, seq) in enumerate(masked_sequences.items()):
                if i < 5 or i >= len(masked_sequences) - 5:  # Log first and last 5 sequences
                    logger.debug(f"  {seq_id}: length={len(seq)} bp, status={matching_status[seq_id]}")
                elif i == 5 and len(masked_sequences) > 10:
                    logger.debug(f"  ... ({len(masked_sequences) - 10} more sequences) ...")
            logger.debug("=== END SEQUENCE INFO ===\n")
        
        # Use the common workflow function to handle the rest
        success = common.run_primer_design_workflow(
            masked_sequences=masked_sequences,
            output_dir=output_dir,
            reference_file=sequence_file,
            mode='direct',
            genes=None,
            coordinate_map=None,
            gff_file=None,
            skip_annotation_filtering=args.noannotation,
            matching_status=matching_status,  # Pass matching status
            all_sequences=all_sequences,  # Pass all sequences including those that failed matching
            add_rows_function=add_rows_for_sequences_without_primers  # Pass the function itself
        )
        
        return success
            
    except Exception as e:
        logger.error(f"Error in direct mode workflow: {e}")
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        return False