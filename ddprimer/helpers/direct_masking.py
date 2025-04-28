#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct masking utilities for ddPrimer pipeline.

This module provides utilities for the direct mode workflow, including:
- Finding sequence locations in the reference genome
- Handling sequences without primers or that failed to match
- Processing VCF files and extracting region-specific variants
"""

import os
import logging
import pandas as pd
import tempfile
import subprocess
import gzip

# Set up logging
logger = logging.getLogger("ddPrimer.helpers")

class DirectMasking:
    """
    Handles all aspects of sequence masking for direct mode:
    - Sequence location finding
    - VCF preparation and region-specific variant extraction
    - Result processing for sequences without primers
    """
    
    def __init__(self):
        """Initialize DirectMasking processor."""
        pass
    
    @staticmethod
    def find_location(sequence, ref_fasta, min_identity=90, min_coverage=90):
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

    @staticmethod
    def add_missing_sequences(df, masked_sequences, matching_status=None):
        """
        Check for sequences that didn't get primers and add them to the results.
        Also add rows for sequences that failed to match to reference genome.
        
        Args:
            df: DataFrame with primer results
            masked_sequences: Dictionary of masked sequences
            matching_status: Dictionary with reference matching status for sequences
            
        Returns:
            pandas.DataFrame: Updated DataFrame with rows for sequences without primers
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
    
    def prepare_vcf_file(self, vcf_file):
        """
        Prepare a VCF file for region queries by ensuring it's compressed with bgzip and indexed with tabix.
        
        Args:
            vcf_file (str): Path to VCF file (.vcf or .vcf.gz)
            
        Returns:
            str: Path to the prepared VCF file
        """
        logger.debug(f"Preparing VCF file: {vcf_file}")
        
        # Check if bgzip and tabix are available
        try:
            subprocess.run(["bgzip", "--version"], capture_output=True, check=False)
            subprocess.run(["tabix", "--version"], capture_output=True, check=False)
        except FileNotFoundError:
            logger.warning("bgzip or tabix not found. Please install them for optimal performance.")
            logger.warning("Falling back to slower VCF parsing method.")
            return vcf_file
        
        # Determine if we need to compress the file
        if vcf_file.endswith('.vcf'):
            compressed_path = vcf_file + '.gz'
            
            # Check if the compressed file already exists and is newer than the original
            if os.path.exists(compressed_path) and os.path.getmtime(compressed_path) > os.path.getmtime(vcf_file):
                logger.debug(f"Using existing compressed VCF: {compressed_path}")
            else:
                logger.info(f"Compressing VCF file with bgzip...")
                try:
                    result = subprocess.run(
                        ["bgzip", "-c", vcf_file], 
                        stdout=open(compressed_path, 'wb'),
                        stderr=subprocess.PIPE,
                        check=True
                    )
                    logger.debug("VCF compression completed successfully")
                except subprocess.CalledProcessError as e:
                    logger.warning(f"Failed to compress VCF file: {e.stderr.decode('utf-8')}")
                    return vcf_file
            
            # Use the compressed file for subsequent operations
            vcf_file = compressed_path
        
        # Ensure the compressed file is indexed
        index_path = vcf_file + '.tbi'
        if not os.path.exists(index_path) or os.path.getmtime(index_path) < os.path.getmtime(vcf_file):
            logger.info(f"Indexing VCF file with tabix...")
            try:
                result = subprocess.run(
                    ["tabix", "-p", "vcf", vcf_file],
                    stderr=subprocess.PIPE,
                    check=True
                )
                logger.debug("VCF indexing completed successfully")
            except subprocess.CalledProcessError as e:
                logger.warning(f"Failed to index VCF file: {e.stderr.decode('utf-8')}")
                # Even if indexing fails, return the compressed file as we can still use it
        
        return vcf_file
    
    def get_region_variants(self, vcf_file, chromosome, start_pos, end_pos):
        """
        Extract variant positions from the VCF file for a specific genomic region.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str): Chromosome name
            start_pos (int): Start position of the region
            end_pos (int): End position of the region
            
        Returns:
            set: Set of variant positions within the specified region
        """
        logger.debug(f"Fetching variants for region {chromosome}:{start_pos}-{end_pos}")
        
        # Prepare the VCF file (compress and index if needed)
        prepared_vcf = self.prepare_vcf_file(vcf_file)
        
        # Use bcftools to query the region
        command = f'bcftools query -f "%POS\\n" -r "{chromosome}:{start_pos}-{end_pos}" "{prepared_vcf}"'
        
        try:
            logger.debug(f"Running command: {command}")
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                # Parse results from bcftools
                variants = set()
                for line in result.stdout.strip().split("\n"):
                    if line:  # Skip empty lines
                        try:
                            pos = int(line.strip())
                            variants.add(pos)
                        except ValueError:
                            logger.warning(f"Invalid position value: {line}")
                            continue
                
                logger.debug(f"Extracted {len(variants)} variants in region {chromosome}:{start_pos}-{end_pos}")
                return variants
            else:
                logger.warning(f"bcftools query failed: {result.stderr}")
                # Fall back to manual parsing
        except Exception as e:
            logger.warning(f"Error running bcftools: {str(e)}")
            # Fall back to manual parsing
        
        # Manual parsing approach for when bcftools fails
        logger.debug("Falling back to manual VCF parsing")
        variants = set()
        
        try:
            # Determine if we need to handle gzip compression
            open_func = gzip.open if prepared_vcf.endswith('.gz') else open
            mode = 'rt' if prepared_vcf.endswith('.gz') else 'r'
            
            with open_func(prepared_vcf, mode) as vcf:
                for line in vcf:
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Parse VCF data line
                    fields = line.strip().split('\t')
                    if len(fields) < 5:  # Minimum valid VCF line
                        continue
                    
                    # Check chromosome match
                    chrom = fields[0]
                    if chrom != chromosome:
                        continue
                    
                    # Check position is in our region
                    try:
                        pos = int(fields[1])
                        if start_pos <= pos <= end_pos:
                            variants.add(pos)
                    except ValueError:
                        logger.warning(f"Invalid position value in VCF: {fields[1]}")
                        continue
            
            logger.debug(f"Extracted {len(variants)} variants in region {chromosome}:{start_pos}-{end_pos} using manual parsing")
            return variants
            
        except Exception as e:
            logger.error(f"Error parsing VCF file {prepared_vcf}: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return set()