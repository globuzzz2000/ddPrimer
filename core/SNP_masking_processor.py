import os
import subprocess
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from tqdm import tqdm
from ..config import Config

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
        print(f"Fetching variant positions from {vcf_file}...")
        
        # Base command
        command = f'bcftools query -f "%CHROM\\t%POS\\n" "{vcf_file}"'
        
        # Add chromosome filter if specified
        if chromosome:
            command += f' -r "{chromosome}"'
        
        # Run command
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print("Error running bcftools:", result.stderr)
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
        
        return variants
    
    def extract_reference_sequences(self, fasta_file):
        """
        Retrieve all reference sequences from the FASTA file.
        
        Args:
            fasta_file (str): Path to FASTA file
            
        Returns:
            dict: Dictionary mapping sequence IDs to sequences
        """
        print(f"Extracting sequences from {fasta_file}...")
        sequences = {}
        
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = str(record.seq)
                
        return sequences
    
    def mask_variants(self, sequence, variant_positions, max_masked_ratio=0.05, 
                      mask_padding=3, window_sliding_step=5, max_window_slide=100):
        """
        Mask regions with variants by replacing them with 'N'.
        Also mask surrounding bases to avoid primers that end right at a SNP.
        
        Args:
            sequence (str): Input DNA sequence
            variant_positions (set): Set of variant positions (1-based)
            max_masked_ratio (float): Maximum allowed ratio of masked bases
            mask_padding (int): Number of bases to mask on each side of a variant
            window_sliding_step (int): Step size for window sliding
            max_window_slide (int): Maximum window slide distance
            
        Returns:
            str: Masked sequence
        """
        if not variant_positions:
            return sequence
            
        print(f"Masking {len(variant_positions)} variants in sequence of length {len(sequence)}.")
        sequence_list = list(sequence)
        
        # Mask variant positions with padding
        for pos in variant_positions:
            if pos < 1 or pos > len(sequence_list):
                continue
                
            # Apply padding - mask positions before and after the variant
            start = max(0, pos - 1 - mask_padding)
            end = min(len(sequence_list), pos + mask_padding)
            
            for i in range(start, end):
                sequence_list[i] = 'N'
        
        masked_sequence = "".join(sequence_list)
        
        # Check if too much of the sequence is masked
        masked_count = masked_sequence.count("N")
        total_length = len(sequence)
        masked_ratio = masked_count / total_length
        
        if masked_ratio > max_masked_ratio:
            print(f"Too many masked bases ({masked_ratio:.2%}), attempting to slide search window...")
            
            # Identify regions with fewer masked bases by sliding windows
            best_ratio = masked_ratio
            best_sequence = masked_sequence
            
            for shift in range(window_sliding_step, max_window_slide, window_sliding_step):
                # Try shifting the sequence to avoid excessive masking
                shifted_sequence = masked_sequence[shift:] + masked_sequence[:shift]
                shifted_ratio = shifted_sequence.count("N") / total_length
                
                if shifted_ratio < best_ratio:
                    best_ratio = shifted_ratio
                    best_sequence = shifted_sequence
                    
                    if best_ratio < max_masked_ratio:
                        print(f"Sliding window successful with shift of {shift} bases. New mask ratio: {best_ratio:.2%}")
                        break
            
            if best_ratio >= max_masked_ratio:
                print(f"Sliding window reduced masking from {masked_ratio:.2%} to {best_ratio:.2%}, still above threshold.")
            
            return best_sequence
        
        return masked_sequence
    
    def prepare_sequences_for_primer_design(self, sequences, variant_positions, min_length=100):
        """
        Prepare sequences for primer design by masking variants and identifying viable regions.
        
        Args:
            sequences (dict): Dictionary of sequence ID to sequence
            variant_positions (dict): Dictionary of sequence ID to variant positions
            min_length (int): Minimum length for a viable region
            
        Returns:
            dict: Dictionary mapping sequence IDs to lists of viable regions
                  Each region is a tuple of (start, end, region_sequence)
        """
        viable_regions_by_sequence = {}
        
        # Process each sequence
        seq_ids = list(sequences.keys())
        if Config.SHOW_PROGRESS:
            seq_iter = tqdm(seq_ids, desc="Finding viable regions")
        else:
            seq_iter = seq_ids
            
        for seq_id in seq_iter:
            sequence = sequences[seq_id]
            # Get variants for this sequence
            seq_variants = variant_positions.get(seq_id, set())
            
            # Mask variants
            masked_sequence = self.mask_variants(sequence, seq_variants)
            
            # Find viable regions
            viable_regions = self._find_viable_regions(masked_sequence, min_length)
            
            # Extract the actual sequences for the viable regions
            regions_with_sequences = []
            for start, end in viable_regions:
                region_seq = masked_sequence[start:end]
                # Only include regions without Ns (should not happen due to find_viable_regions)
                if 'N' not in region_seq:
                    regions_with_sequences.append((start, end, region_seq))
            
            if regions_with_sequences:
                viable_regions_by_sequence[seq_id] = regions_with_sequences
                print(f"Found {len(regions_with_sequences)} viable regions for {seq_id}")
            else:
                print(f"No viable regions found for {seq_id}")
        
        return viable_regions_by_sequence
    
    def mask_variants(self, sequence, variant_positions, max_masked_ratio=0.05, 
                      mask_padding=3, window_sliding_step=5, max_window_slide=100):
        """
        Mask regions with variants by replacing them with 'N'.
        Also mask surrounding bases to avoid primers that end right at a SNP.
        
        Args:
            sequence (str): Input DNA sequence
            variant_positions (set): Set of variant positions (1-based)
            max_masked_ratio (float): Maximum allowed ratio of masked bases
            mask_padding (int): Number of bases to mask on each side of a variant
            window_sliding_step (int): Step size for window sliding
            max_window_slide (int): Maximum window slide distance
            
        Returns:
            str: Masked sequence
        """
        if not variant_positions:
            return sequence
            
        print(f"Masking {len(variant_positions)} variants in sequence of length {len(sequence)}.")
        sequence_list = list(sequence)
        
        # Mask variant positions with padding
        variant_pos_list = list(variant_positions)
        if Config.SHOW_PROGRESS and len(variant_pos_list) > 1000:  # Only show progress for large sets
            variant_iter = tqdm(variant_pos_list, desc="Masking variants")
        else:
            variant_iter = variant_pos_list
            
        for pos in variant_iter:
            if pos < 1 or pos > len(sequence_list):
                continue
                
            # Apply padding - mask positions before and after the variant
            start = max(0, pos - 1 - mask_padding)
            end = min(len(sequence_list), pos + mask_padding)
            
            for i in range(start, end):
                sequence_list[i] = 'N'
        
        masked_sequence = "".join(sequence_list)
        
        # Check if too much of the sequence is masked
        masked_count = masked_sequence.count("N")
        total_length = len(sequence)
        masked_ratio = masked_count / total_length
        
        if masked_ratio > max_masked_ratio:
            print(f"Too many masked bases ({masked_ratio:.2%}), attempting to slide search window...")
            
            # Identify regions with fewer masked bases by sliding windows
            best_ratio = masked_ratio
            best_sequence = masked_sequence
            
            shift_range = range(window_sliding_step, max_window_slide, window_sliding_step)
            if Config.SHOW_PROGRESS:
                shift_iter = tqdm(shift_range, desc="Sliding window search")
            else:
                shift_iter = shift_range
                
            for shift in shift_iter:
                # Try shifting the sequence to avoid excessive masking
                shifted_sequence = masked_sequence[shift:] + masked_sequence[:shift]
                shifted_ratio = shifted_sequence.count("N") / total_length
                
                if shifted_ratio < best_ratio:
                    best_ratio = shifted_ratio
                    best_sequence = shifted_sequence
                    
                    if best_ratio < max_masked_ratio:
                        print(f"Sliding window successful with shift of {shift} bases. New mask ratio: {best_ratio:.2%}")
                        break
            
            if best_ratio >= max_masked_ratio:
                print(f"Sliding window reduced masking from {masked_ratio:.2%} to {best_ratio:.2%}, still above threshold.")
            
            return best_sequence
        
        return masked_sequence