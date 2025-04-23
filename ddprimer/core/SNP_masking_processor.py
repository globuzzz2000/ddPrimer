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
        
        print(f"Extracted {sum(len(v) for v in variants.values())} variant positions from {len(variants)} chromosomes")
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
    
    def mask_variants(self, sequence, variant_positions, mask_padding=3, min_region_length=100):
        """
        Mask variants in the sequence and identify unmasked regions suitable for primer design.
        
        Args:
            sequence (str): Input DNA sequence
            variant_positions (set): Set of variant positions (1-based)
            mask_padding (int): Number of bases to mask on each side of a variant
            min_region_length (int): Minimum length required for an unmasked region
            
        Returns:
            str: Masked sequence with the best unmasked regions preserved
        """
        if not variant_positions:
            print("No variants to mask, returning original sequence")
            return sequence
            
        # Step 1: Create the masked sequence
        sequence_list = list(sequence)
        total_length = len(sequence)
        
        # Track which positions are masked for variant visualization
        masked_positions = set()
        
        # Mask variant positions with padding
        variant_pos_list = sorted(list(variant_positions))
        variant_density = len(variant_pos_list) / total_length if total_length > 0 else 0
        print(f"Processing {len(variant_pos_list)} variants (density: {variant_density:.4%})")
        
        if Config.SHOW_PROGRESS and len(variant_pos_list) > 1000:  # Only show progress for large sets
            variant_iter = tqdm(variant_pos_list, desc="Masking variants")
        else:
            variant_iter = variant_pos_list
            
        for pos in variant_iter:
            # Convert to 0-based for sequence indexing (variants are 1-based)
            pos_0based = pos - 1
            
            if pos_0based < 0 or pos_0based >= len(sequence_list):
                continue
                
            # Apply padding - mask positions before and after the variant
            start = max(0, pos_0based - mask_padding)
            end = min(len(sequence_list), pos_0based + mask_padding + 1)
            
            for i in range(start, end):
                sequence_list[i] = 'N'
                masked_positions.add(i)
        
        masked_sequence = "".join(sequence_list)
        masked_count = len(masked_positions)
        masked_ratio = masked_count / total_length if total_length > 0 else 1.0
        
        print(f"Initial masking: {masked_count}/{total_length} bases ({masked_ratio:.2%})")
        
        # Step 2: Identify contiguous unmasked regions
        unmasked_regions = self.find_unmasked_regions(masked_sequence, min_length=min_region_length)
        
        if not unmasked_regions:
            print(f"No unmasked regions of length >= {min_region_length} found.")
            # Try with smaller regions if none found with default size
            smaller_regions = self.find_unmasked_regions(masked_sequence, min_length=50)
            if smaller_regions:
                print(f"Found {len(smaller_regions)} smaller unmasked regions (min length: 50)")
                unmasked_regions = smaller_regions
            else:
                print("No usable unmasked regions found even with reduced length criteria.")
                # Return sequence anyway for downstream handling
                return masked_sequence
        
        # Sort regions by length (descending)
        unmasked_regions.sort(key=lambda x: x[1] - x[0], reverse=True)
        
        print(f"Found {len(unmasked_regions)} unmasked regions, largest is {unmasked_regions[0][1] - unmasked_regions[0][0]} bp")
        
        # Step 3: Create an optimized sequence with only the best unmasked regions
        # If we have at least one good region, create a new sequence with only that region unmasked
        if unmasked_regions:
            # Choose the best region (the longest one)
            best_region = unmasked_regions[0]
            
            # Create a new sequence with everything masked except the best region
            optimized_sequence_list = ['N'] * total_length
            start, end = best_region
            
            # Copy the unmasked region
            for i in range(start, end):
                optimized_sequence_list[i] = sequence[i]
            
            # If we have multiple good regions, include the top few
            max_regions_to_include = 3  # Include up to 3 good regions
            for region in unmasked_regions[1:max_regions_to_include]:
                if (region[1] - region[0]) >= min_region_length * 0.75:  # Only include reasonably sized regions
                    region_start, region_end = region
                    for i in range(region_start, region_end):
                        optimized_sequence_list[i] = sequence[i]
            
            optimized_sequence = "".join(optimized_sequence_list)
            
            # Calculate what percentage of the sequence is now usable
            usable_bases = sum(1 for c in optimized_sequence if c != 'N')
            usable_ratio = usable_bases / total_length if total_length > 0 else 0
            
            print(f"Optimized sequence has {usable_bases}/{total_length} usable bases ({usable_ratio:.2%})")
            
            return optimized_sequence
        
        # If no good regions, return the original masked sequence
        return masked_sequence
    
    def find_unmasked_regions(self, sequence, min_length=100):
        """
        Find contiguous regions in the sequence that don't contain 'N's.
        
        Args:
            sequence (str): Input sequence with masked positions (N's)
            min_length (int): Minimum length for a region to be considered
            
        Returns:
            list: List of tuples (start, end) for each unmasked region
        """
        regions = []
        region_start = None
        
        for i, base in enumerate(sequence):
            if base != 'N':
                # Start a new region if we're not in one
                if region_start is None:
                    region_start = i
            else:
                # End the current region if we were in one
                if region_start is not None:
                    region_end = i
                    region_length = region_end - region_start
                    
                    if region_length >= min_length:
                        regions.append((region_start, region_end))
                    
                    region_start = None
        
        # Handle the case where the sequence ends with a valid region
        if region_start is not None:
            region_end = len(sequence)
            region_length = region_end - region_start
            
            if region_length >= min_length:
                regions.append((region_start, region_end))
        
        return regions
    
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
            
            # Mask variants and find optimal regions
            masked_sequence = self.mask_variants(sequence, seq_variants, min_region_length=min_length)
            
            # Find viable (unmasked) regions
            viable_regions = self.find_unmasked_regions(masked_sequence, min_length)
            
            # Extract the actual sequences for the viable regions
            regions_with_sequences = []
            for start, end in viable_regions:
                region_seq = masked_sequence[start:end]
                # Verify no N's in the region (should be guaranteed by find_unmasked_regions)
                if 'N' not in region_seq:
                    regions_with_sequences.append((start, end, region_seq))
            
            if regions_with_sequences:
                viable_regions_by_sequence[seq_id] = regions_with_sequences
                print(f"Found {len(regions_with_sequences)} viable regions for {seq_id}")
            else:
                print(f"No viable regions found for {seq_id}")
        
        return viable_regions_by_sequence