#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct masking utilities for ddPrimer pipeline.

This module provides utilities for the direct mode workflow, including:
- Finding sequence locations in the reference genome
- Handling sequences without primers or that failed to match
"""

import os
import logging
import pandas as pd
import tempfile
import subprocess

from ..config.exceptions import SequenceProcessingError

# Set up logger
logger = logging.getLogger("ddPrimer.helpers")


class DirectMasking:
    """
    Handles sequence location finding and result processing for direct mode.
    Uses SNPMaskingProcessor for variant handling.
    """
    
    @staticmethod
    def find_location(sequence, ref_fasta, min_identity=90, min_coverage=90):
        """
        Find the location of a sequence in the reference genome using BLAST.
        
        Args:
            sequence (str): Query sequence to locate
            ref_fasta (str): Path to reference FASTA file
            min_identity (float): Minimum percent identity (default: 90)
            min_coverage (float): Minimum query coverage (default: 90)
            
        Returns:
            tuple: (chromosome, start_position, end_position, percent_identity) or (None, None, None, None) if not found
            
        Raises:
            SequenceProcessingError: If there is an error during BLAST search
        """
        temp_file_path = None
        
        try:
            # Create temporary file for query sequence
            with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_file:
                temp_file_path = temp_file.name
                temp_file.write(">query\n")
                temp_file.write(sequence)
            
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
            logger.error(f"BLAST error: {str(e)}")
            logger.debug(f"BLAST stderr: {e.stderr}")
            raise SequenceProcessingError(f"BLAST error: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error in BLAST search: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(f"BLAST search failed: {str(e)}")
        
        finally:
            # Clean up temporary file
            if temp_file_path and os.path.exists(temp_file_path):
                try:
                    os.unlink(temp_file_path)
                except Exception as e:
                    logger.warning(f"Could not remove temporary file {temp_file_path}: {str(e)}")

    @staticmethod
    def add_missing_sequences(df, masked_sequences, matching_status=None):
        """
        Check for sequences that didn't get primers and add them to the results.
        Also add rows for sequences that failed to match to reference genome.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer results
            masked_sequences (dict): Dictionary of masked sequences
            matching_status (dict, optional): Dictionary with reference matching status for sequences
                
        Returns:
            pandas.DataFrame: Updated DataFrame with rows for sequences without primers
        """
        logger.debug("Checking for sequences without primers and reference match failures...")
        
        # If there's no Gene column, we can't properly identify which sequences have primers.
        # This shouldn't happen in normal operation of direct mode, as the column is explicitly
        # added in process_direct_mode_fragments, but we'll handle it defensively.
        if "Gene" not in df.columns:
            logger.debug("No 'Gene' column found in DataFrame, skipping missing sequences check")
            return df
        
        # 1) Which sequences actually got primers?
        sequences_with_primers = set()
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