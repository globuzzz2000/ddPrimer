#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP checker module for ddPrimer pipeline.

This module contains functionality to:
1. Match primer sequences against a reference genome
2. Check if primers/probes overlap with SNPs
3. Filter primers based on SNP overlap
"""

import os
import logging
import pandas as pd
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import subprocess
import re
from tqdm import tqdm

from ..config import Config
from .sequence_utils import SequenceUtils


class SNPChecker:
    """Handles checking of primers/probes against SNPs in a reference genome."""
    
    def __init__(self, ref_fasta=None, ref_vcf=None):
        """
        Initialize SNP checker.
        
        Args:
            ref_fasta (str): Path to reference FASTA file
            ref_vcf (str): Path to VCF file with variants
        """
        self.ref_fasta = ref_fasta
        self.ref_vcf = ref_vcf
        self.logger = logging.getLogger("ddPrimer")
        self.variant_positions = {}
        
        # Load variant positions if VCF is provided
        if ref_vcf:
            self._load_variant_positions()
    
    def _load_variant_positions(self):
        """Load variant positions from VCF file."""
        self.logger.info(f"Loading variant positions from {self.ref_vcf}")
        try:
            from ..core.SNP_masking_processor import SNPMaskingProcessor
            
            snp_processor = SNPMaskingProcessor()
            self.variant_positions = snp_processor.get_variant_positions(self.ref_vcf)
            
            # Log variant counts per chromosome
            for chrom, positions in self.variant_positions.items():
                self.logger.debug(f"Loaded {len(positions)} variants for {chrom}")
            
            total_variants = sum(len(positions) for positions in self.variant_positions.values())
            self.logger.info(f"Loaded {total_variants} variants from {len(self.variant_positions)} chromosomes")
            
        except Exception as e:
            self.logger.error(f"Error loading variant positions: {e}")
            self.logger.debug(e, exc_info=True)
            raise
    
    def _create_temporary_fasta(self, sequences):
        """
        Create a temporary FASTA file with primer/probe sequences.
        
        Args:
            sequences (dict): Dictionary mapping sequence names to sequences
            
        Returns:
            str: Path to temporary FASTA file
        """
        try:
            # Create temporary file
            fd, temp_path = tempfile.mkstemp(suffix='.fasta')
            os.close(fd)
            
            # Write sequences to temporary file
            with open(temp_path, 'w') as f:
                for name, seq in sequences.items():
                    # Sanitize name for FASTA
                    sanitized_name = re.sub(r'\s+', '_', name)
                    f.write(f">{sanitized_name}\n{seq}\n")
            
            return temp_path
            
        except Exception as e:
            self.logger.error(f"Error creating temporary FASTA: {e}")
            self.logger.debug(e, exc_info=True)
            raise
    
    def _run_blast_alignment(self, query_file):
        """
        Run BLAST to align primers/probes against reference genome.
        
        Args:
            query_file (str): Path to FASTA file with primer/probe sequences
            
        Returns:
            list: List of BLAST hits
        """
        self.logger.info("Running BLAST to align primers/probes against reference genome")
        
        try:
            # Create temporary file for BLAST output
            fd, blast_out = tempfile.mkstemp(suffix='.blast')
            os.close(fd)
            
            # Set up BLAST command
            # Using blastn with parameters tuned for short sequences
            cmd = [
                'blastn',
                '-query', query_file,
                '-subject', self.ref_fasta,
                '-out', blast_out,
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
                '-word_size', '7',  # Smaller word size for short sequences
                '-evalue', '1000',  # More permissive e-value for short sequences
                '-dust', 'no',      # Disable dust filtering
                '-perc_identity', '95'  # High identity threshold
            ]
            
            # Run BLAST
            self.logger.debug(f"Running BLAST command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            
            # Parse BLAST output
            hits = []
            with open(blast_out, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:  # Ensure all fields are present
                        hit = {
                            'qseqid': parts[0],    # Query sequence ID
                            'sseqid': parts[1],    # Subject (reference) sequence ID
                            'pident': float(parts[2]),  # Percentage identity
                            'length': int(parts[3]),    # Alignment length
                            'qstart': int(parts[6]),    # Query start
                            'qend': int(parts[7]),      # Query end
                            'sstart': int(parts[8]),    # Subject start
                            'send': int(parts[9]),      # Subject end
                            'evalue': float(parts[10]), # E-value
                            'bitscore': float(parts[11]) # Bit score
                        }
                        hits.append(hit)
            
            # Clean up temporary files
            os.remove(blast_out)
            
            self.logger.info(f"BLAST alignment complete - found {len(hits)} hits")
            return hits
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error running BLAST: {e}")
            self.logger.debug(e, exc_info=True)
            # Since BLAST failed, try with a fallback method
            self.logger.info("Falling back to direct sequence search in reference genome")
            return self._fallback_sequence_search(query_file)
        except Exception as e:
            self.logger.error(f"Error in BLAST alignment: {e}")
            self.logger.debug(e, exc_info=True)
            raise
    
    def _fallback_sequence_search(self, query_file):
        """
        Fallback method to search for sequences in reference genome without BLAST.
        Used if BLAST fails or is not available.
        
        Args:
            query_file (str): Path to FASTA file with primer/probe sequences
            
        Returns:
            list: List of hits in BLAST-like format
        """
        self.logger.info("Using fallback sequence search method")
        
        try:
            # Load query sequences
            queries = {}
            for record in SeqIO.parse(query_file, "fasta"):
                queries[record.id] = str(record.seq).upper()
            
            # Load reference sequences
            ref_seqs = {}
            for record in SeqIO.parse(self.ref_fasta, "fasta"):
                ref_seqs[record.id] = str(record.seq).upper()
            
            # Search for sequences
            hits = []
            
            for q_id, q_seq in tqdm(queries.items(), desc="Searching sequences"):
                for ref_id, ref_seq in ref_seqs.items():
                    # Find all occurrences of query in reference
                    matches = self._find_all_occurrences(ref_seq, q_seq)
                    
                    for pos in matches:
                        hit = {
                            'qseqid': q_id,
                            'sseqid': ref_id,
                            'pident': 100.0,  # Exact match
                            'length': len(q_seq),
                            'qstart': 1,
                            'qend': len(q_seq),
                            'sstart': pos + 1,  # 1-based position
                            'send': pos + len(q_seq),
                            'evalue': 0.0,     # Placeholder
                            'bitscore': 100.0  # Placeholder
                        }
                        hits.append(hit)
            
            self.logger.info(f"Fallback search complete - found {len(hits)} hits")
            return hits
            
        except Exception as e:
            self.logger.error(f"Error in fallback sequence search: {e}")
            self.logger.debug(e, exc_info=True)
            raise
    
    def _find_all_occurrences(self, text, pattern):
        """
        Find all occurrences of a pattern in text.
        
        Args:
            text (str): Text to search in
            pattern (str): Pattern to search for
            
        Returns:
            list: List of starting positions (0-based)
        """
        positions = []
        start = 0
        
        while True:
            start = text.find(pattern, start)
            if start == -1:
                break
            positions.append(start)
            start += 1
            
        return positions
    
    def _check_snp_overlap(self, hit):
        """
        Check if a hit overlaps with any SNPs.
        
        Args:
            hit (dict): BLAST hit information
            
        Returns:
            bool, list: True if overlaps, False otherwise; and list of overlapping SNP positions
        """
        chrom = hit['sseqid']
        start = min(hit['sstart'], hit['send'])
        end = max(hit['sstart'], hit['send'])
        
        # Check if chromosome exists in variant positions
        if chrom not in self.variant_positions:
            return False, []
        
        # Check for overlaps
        overlapping_snps = []
        for pos in range(start, end + 1):
            if pos in self.variant_positions[chrom]:
                overlapping_snps.append(pos)
        
        return len(overlapping_snps) > 0, overlapping_snps
    
    def check_primers_for_snps(self, df):
        """
        Check primers and probes in DataFrame for SNP overlaps.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer information
            
        Returns:
            pandas.DataFrame: Updated DataFrame with SNP overlap information
        """
        self.logger.info("\nChecking primers and probes for SNP overlaps...")
        
        if self.ref_fasta is None or self.ref_vcf is None:
            self.logger.warning("Reference FASTA or VCF file not provided, skipping SNP check")
            return df
        
        try:
            # Extract primer and probe sequences
            primers = {}
            
            # Add forward primers
            for idx, row in df.iterrows():
                primer_f = row.get('Primer F')
                if pd.notna(primer_f) and primer_f and primer_f != "No suitable primers found":
                    primers[f"F_{idx}_{row.get('Gene', 'Unknown')}"] = primer_f
            
            # Add reverse primers
            for idx, row in df.iterrows():
                primer_r = row.get('Primer R')
                if pd.notna(primer_r) and primer_r:
                    primers[f"R_{idx}_{row.get('Gene', 'Unknown')}"] = primer_r
            
            # Add probes if present
            if 'Probe' in df.columns:
                for idx, row in df.iterrows():
                    probe = row.get('Probe')
                    if pd.notna(probe) and probe:
                        primers[f"P_{idx}_{row.get('Gene', 'Unknown')}"] = probe
            
            self.logger.debug(f"Checking {len(primers)} sequences for SNP overlaps")
            
            # Create temporary FASTA file
            if not primers:
                self.logger.warning("No valid sequences to check for SNPs")
                return df
                
            temp_fasta = self._create_temporary_fasta(primers)
            
            # Run BLAST alignment
            blast_hits = self._run_blast_alignment(temp_fasta)
            
            # Clean up temporary file
            os.remove(temp_fasta)
            
            # Process hits and check for SNP overlaps
            overlaps = {}
            overlapping_snps = {}
            
            for hit in tqdm(blast_hits, desc="Checking SNP overlaps"):
                qseqid = hit['qseqid']
                has_overlap, snps = self._check_snp_overlap(hit)
                
                if has_overlap:
                    overlaps[qseqid] = True
                    overlapping_snps[qseqid] = snps
                elif qseqid not in overlaps:
                    overlaps[qseqid] = False
                    overlapping_snps[qseqid] = []
            
            # Add SNP overlap information to DataFrame
            df['Primer F SNP Overlap'] = False
            df['Primer R SNP Overlap'] = False
            df['SNP Overlap Positions'] = ''
            
            if 'Probe' in df.columns:
                df['Probe SNP Overlap'] = False
            
            for idx, row in df.iterrows():
                # Check forward primer
                f_id = f"F_{idx}_{row.get('Gene', 'Unknown')}"
                if f_id in overlaps:
                    df.at[idx, 'Primer F SNP Overlap'] = overlaps[f_id]
                    if overlaps[f_id]:
                        df.at[idx, 'SNP Overlap Positions'] = f"F:{','.join(map(str, overlapping_snps[f_id]))}"
                
                # Check reverse primer
                r_id = f"R_{idx}_{row.get('Gene', 'Unknown')}"
                if r_id in overlaps:
                    df.at[idx, 'Primer R SNP Overlap'] = overlaps[r_id]
                    if overlaps[r_id]:
                        snp_pos = df.at[idx, 'SNP Overlap Positions']
                        df.at[idx, 'SNP Overlap Positions'] = f"{snp_pos};R:{','.join(map(str, overlapping_snps[r_id]))}"
                
                # Check probe if present
                if 'Probe' in df.columns:
                    p_id = f"P_{idx}_{row.get('Gene', 'Unknown')}"
                    if p_id in overlaps:
                        df.at[idx, 'Probe SNP Overlap'] = overlaps[p_id]
                        if overlaps[p_id]:
                            snp_pos = df.at[idx, 'SNP Overlap Positions']
                            df.at[idx, 'SNP Overlap Positions'] = f"{snp_pos};P:{','.join(map(str, overlapping_snps[p_id]))}"
            
            # Log summary of SNP overlaps
            total_f_overlaps = df['Primer F SNP Overlap'].sum()
            total_r_overlaps = df['Primer R SNP Overlap'].sum()
            
            self.logger.info(f"Found {total_f_overlaps} forward primers with SNP overlaps")
            self.logger.info(f"Found {total_r_overlaps} reverse primers with SNP overlaps")
            
            if 'Probe' in df.columns:
                total_p_overlaps = df['Probe SNP Overlap'].sum()
                self.logger.info(f"Found {total_p_overlaps} probes with SNP overlaps")
            
            return df
            
        except Exception as e:
            self.logger.error(f"Error checking primers for SNP overlaps: {e}")
            self.logger.debug(e, exc_info=True)
            # Don't raise the exception - just return the original DataFrame
            return df
    
    def filter_by_snp_overlap(self, df):
        """
        Filter primers by SNP overlap.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer information
            
        Returns:
            pandas.DataFrame: Filtered DataFrame
        """
        self.logger.info("\nFiltering primers by SNP overlap...")
        
        # Check if SNP overlap columns exist
        if 'Primer F SNP Overlap' not in df.columns or 'Primer R SNP Overlap' not in df.columns:
            self.logger.warning("SNP overlap information not found, skipping SNP filtering")
            return df
        
        # Store original DataFrame for comparison
        df_before = df.copy()
        
        # Filter out primers that overlap with SNPs
        df = df[~(df['Primer F SNP Overlap'] | df['Primer R SNP Overlap'])].copy()
        
        # If probes are present, also filter by probe SNP overlap
        if 'Probe' in df.columns and 'Probe SNP Overlap' in df.columns:
            df = df[~df['Probe SNP Overlap']].copy()
        
        # Reset index
        df = df.reset_index(drop=True)
        
        # Log filtering results
        filtered_out = set(df_before["Gene"]) - set(df["Gene"])
        self.logger.info(f"After SNP filtering: {len(df)}/{len(df_before)} primers retained")
        self.logger.debug(f"Primers filtered out due to SNP overlap: {len(filtered_out)}")
        
        if filtered_out:
            self.logger.debug(f"Filtered genes: {', '.join(filtered_out)}")
        
        return df