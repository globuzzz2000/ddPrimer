import os
import logging
import subprocess
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from tqdm import tqdm
from ..config import Config

# Set up logging
logger = logging.getLogger("ddPrimer")

class SNPMaskingProcessor:
    """Handles masking of SNPs in sequences to prepare them for primer design."""
    
    def __init__(self):
        """Initialize SNP masking processor."""
        pass
    
    def get_variant_positions(self, vcf_file, chromosome=None):
        """
        Extract variant positions from the VCF file using bcftools.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
        """
        logger.info(f"Fetching variant positions from {vcf_file}...")
        
        # Base command
        command = f'bcftools query -f "%CHROM\\t%POS\\n" "{vcf_file}"'
        
        # Add chromosome filter if specified
        if chromosome:
            command += f' -r "{chromosome}"'
        
        # Run command
        logger.debug(f"Running command: {command}")
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Error running bcftools: {result.stderr}")
            return {}
        
        # Parse results
        variants = {}
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
                
            parts = line.split("\t")
            if len(parts) == 2:
                chrom, pos = parts
                pos = int(pos)
                
                if chrom not in variants:
                    variants[chrom] = set()
                    
                variants[chrom].add(pos)
        
        logger.info(f"Extracted {sum(len(positions) for positions in variants.values())} variants from {len(variants)} chromosomes")
        return variants
    
    def extract_reference_sequences(self, fasta_file):
        """
        Retrieve all reference sequences from the FASTA file.
        
        Args:
            fasta_file (str): Path to FASTA file
            
        Returns:
            dict: Dictionary mapping sequence IDs to sequences
        """
        logger.info(f"Extracting sequences from {fasta_file}...")
        sequences = {}
        
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = str(record.seq)
                
        logger.debug(f"Extracted {len(sequences)} sequences from FASTA file")
        return sequences
    
    def mask_variants(self, sequence, variant_positions):
        """
        Simple variant masking - replace only the exact SNP positions with 'N'.
        No padding - only the exact variant position is masked.
        
        Args:
            sequence (str): Input DNA sequence
            variant_positions (set): Set of variant positions (1-based)
            
        Returns:
            str: Masked sequence with SNPs replaced by 'N'
        """
        if not variant_positions:
            logger.debug("No variant positions to mask")
            return sequence
            
        sequence_list = list(sequence)
        masked_positions = 0
        
        # Mask only the exact variant positions
        variant_pos_list = list(variant_positions)
        logger.debug(f"Masking {len(variant_pos_list)} variant positions")
        
        if Config.SHOW_PROGRESS and len(variant_pos_list) > 1000:  # Only show progress for large sets
            variant_iter = tqdm(variant_pos_list, desc="Masking variants")
        else:
            variant_iter = variant_pos_list
            
        for pos in variant_iter:
            # Convert 1-based VCF position to 0-based sequence index
            idx = pos - 1
            
            if idx < 0 or idx >= len(sequence_list):
                continue
                
            # Mask only the exact position
            if sequence_list[idx] != 'N':  # Only count newly masked positions
                sequence_list[idx] = 'N'
                masked_positions += 1
        
        masked_sequence = "".join(sequence_list)
        
        # Calculate masking statistics
        masked_ratio = masked_positions / len(sequence) if len(sequence) > 0 else 0
        logger.debug(f"Masked {masked_positions} positions ({masked_ratio:.2%} of sequence)")
        
        return masked_sequence
    
    def mask_sequences_for_primer_design(self, sequences, variant_positions):
        """
        Mask variants in sequences to prepare them for primer design.
        Only the exact SNP positions are masked.
        
        Args:
            sequences (dict): Dictionary of sequence ID to sequence
            variant_positions (dict): Dictionary of sequence ID to variant positions
            
        Returns:
            dict: Dictionary of masked sequences
        """
        masked_sequences = {}
        
        # Process each sequence
        seq_ids = list(sequences.keys())
        logger.debug(f"Masking {len(seq_ids)} sequences")
        
        if Config.SHOW_PROGRESS:
            seq_iter = tqdm(seq_ids, desc="Masking sequences")
        else:
            seq_iter = seq_ids
            
        for seq_id in seq_iter:
            sequence = sequences[seq_id]
            # Get variants for this sequence
            seq_variants = variant_positions.get(seq_id, set())
            
            if seq_variants:
                logger.debug(f"Masking {len(seq_variants)} variants in sequence {seq_id} ({len(sequence)} bp)")
                
                # Mask variants
                masked_sequence = self.mask_variants(sequence, seq_variants)
                masked_sequences[seq_id] = masked_sequence
                
                # Calculate masking statistics
                masked_count = masked_sequence.count('N')
                masked_ratio = masked_count / len(sequence) if len(sequence) > 0 else 0
                logger.debug(f"Sequence {seq_id}: {masked_count}/{len(sequence)} bases masked ({masked_ratio:.2%})")
            else:
                logger.debug(f"No variants to mask in sequence {seq_id}")
                masked_sequences[seq_id] = sequence
        
        return masked_sequences
    
    def combine_masked_sequences(self, masked_sequences1, masked_sequences2):
        """
        Combine two sets of masked sequences, masking a position if it's masked in either set.
        
        Args:
            masked_sequences1 (dict): First set of masked sequences
            masked_sequences2 (dict): Second set of masked sequences
            
        Returns:
            dict: Combined masked sequences
        """
        combined_sequences = {}
        
        # Get all sequence IDs from both sets
        all_seq_ids = set(masked_sequences1.keys()) | set(masked_sequences2.keys())
        logger.debug(f"Combining masking from two sets ({len(masked_sequences1)} and {len(masked_sequences2)} sequences)")
        
        for seq_id in all_seq_ids:
            # If sequence only exists in one set, use that version
            if seq_id not in masked_sequences1:
                logger.debug(f"Sequence {seq_id} only exists in second set")
                combined_sequences[seq_id] = masked_sequences2[seq_id]
                continue
            elif seq_id not in masked_sequences2:
                logger.debug(f"Sequence {seq_id} only exists in first set")
                combined_sequences[seq_id] = masked_sequences1[seq_id]
                continue
            
            # If sequence exists in both sets, combine the masking
            seq1 = masked_sequences1[seq_id]
            seq2 = masked_sequences2[seq_id]
            
            # Ensure sequences are the same length
            if len(seq1) != len(seq2):
                logger.warning(f"Sequence length mismatch for {seq_id}. Using first sequence.")
                combined_sequences[seq_id] = seq1
                continue
            
            # Combine masking - if either has an N, use N
            combined_seq = ""
            for i in range(len(seq1)):
                if seq1[i] == 'N' or seq2[i] == 'N':
                    combined_seq += 'N'
                else:
                    combined_seq += seq1[i]
            
            # Calculate masking statistics
            n_count = combined_seq.count('N')
            n_ratio = n_count / len(combined_seq) if len(combined_seq) > 0 else 0
            logger.debug(f"Combined sequence {seq_id}: {n_count}/{len(combined_seq)} bases masked ({n_ratio:.2%})")
            
            combined_sequences[seq_id] = combined_seq
        
        return combined_sequences