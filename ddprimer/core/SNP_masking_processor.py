import os
import logging
import subprocess
import gzip
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from tqdm import tqdm
from ..config import Config
from ..config.exceptions import FileError, SequenceProcessingError

# Set up logging
logger = logging.getLogger("ddPrimer.snp_masking_processor")

class SNPMaskingProcessor:
    """Handles masking of SNPs in sequences to prepare them for primer design."""
    
    def __init__(self):
        """Initialize SNP masking processor."""
        logger.debug("Initialized SNPMaskingProcessor")
    
    def get_variant_positions(self, vcf_file, chromosome=None):
        """
        Extract variant positions from the VCF file using bcftools.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        logger.debug(f"\nFetching variant positions from {vcf_file}")
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
        
        # Base command
        command = f'bcftools query -f "%CHROM\\t%POS\\n" "{vcf_file}"'
        
        # Add chromosome filter if specified
        if chromosome:
            command += f' -r "{chromosome}"'
            logger.debug(f"Filtering variants for chromosome: {chromosome}")
        
        # Run command
        logger.debug(f"Running command: {command}")
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"Error running bcftools: {result.stderr}")
                raise FileError(f"bcftools query failed: {result.stderr}")
            
            # Parse results
            variants = {}
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                    
                parts = line.split("\t")
                if len(parts) == 2:
                    chrom, pos = parts
                    try:
                        pos = int(pos)
                        
                        if chrom not in variants:
                            variants[chrom] = set()
                            
                        variants[chrom].add(pos)
                    except ValueError:
                        logger.warning(f"Invalid position value in VCF: {pos}")
            
            total_variants = sum(len(positions) for positions in variants.values())
            logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes")
            return variants
            
        except subprocess.SubprocessError as e:
            logger.error(f"Failed to run bcftools: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Failed to run bcftools: {e}")
    
    def prepare_vcf_file(self, vcf_file):
        """
        Prepare a VCF file for region queries by ensuring it's compressed with bgzip and indexed with tabix.
        
        Args:
            vcf_file (str): Path to VCF file (.vcf or .vcf.gz)
            
        Returns:
            str: Path to the prepared VCF file
            
        Raises:
            FileError: If VCF file cannot be prepared
        """
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
        
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
            
        Raises:
            SequenceProcessingError: If region variants cannot be extracted
        """
        logger.debug(f"Fetching variants for region {chromosome}:{start_pos}-{end_pos}")
        
        # Prepare the VCF file (compress and index if needed)
        try:
            prepared_vcf = self.prepare_vcf_file(vcf_file)
            
            # Use bcftools to query the region
            command = f'bcftools query -f "%POS\\n" -r "{chromosome}:{start_pos}-{end_pos}" "{prepared_vcf}"'
            
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
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            # Fall back to manual parsing
        
        # Manual parsing approach for when bcftools fails
        logger.debug("Falling back to manual VCF parsing")
        return self._extract_variants_manually(prepared_vcf, chromosome, start_pos, end_pos)
    
    def _extract_variants_manually(self, vcf_file, chrom, start, end):
        """
        Extract variants from a VCF file for a specific region using manual parsing.
        This is a fallback method when bcftools fails.
        
        Args:
            vcf_file (str): Path to VCF file
            chrom (str): Chromosome name
            start (int): Start position
            end (int): End position
            
        Returns:
            set: Set of variant positions in the region
        """
        variants = set()
        
        try:
            # Determine if we need to handle gzip compression
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with open_func(vcf_file, mode) as vcf:
                for line in vcf:
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Parse VCF data line
                    fields = line.strip().split('\t')
                    if len(fields) < 5:  # Minimum valid VCF line
                        continue
                    
                    # Check chromosome match
                    vcf_chrom = fields[0]
                    if vcf_chrom != chrom:
                        continue
                    
                    # Check position is in our region
                    try:
                        pos = int(fields[1])
                        if start <= pos <= end:
                            variants.add(pos)
                    except ValueError:
                        logger.warning(f"Invalid position value in VCF: {fields[1]}")
                        continue
            
            logger.debug(f"Manually extracted {len(variants)} variants in region {chrom}:{start}-{end}")
            
        except Exception as e:
            logger.error(f"Error in manual VCF parsing for region {chrom}:{start}-{end}: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        
        return variants
    
    def extract_reference_sequences(self, fasta_file):
        """
        Retrieve all reference sequences from the FASTA file.
        
        Args:
            fasta_file (str): Path to FASTA file
            
        Returns:
            dict: Dictionary mapping sequence IDs to sequences
            
        Raises:
            FileError: If FASTA file cannot be read
        """
        logger.info(f"Extracting sequences from {fasta_file}")
        
        if not os.path.exists(fasta_file):
            logger.error(f"FASTA file not found: {fasta_file}")
            raise FileError(f"FASTA file not found: {fasta_file}")
            
        sequences = {}
        
        try:
            with open(fasta_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequences[record.id] = str(record.seq)
                    
            logger.info(f"Extracted {len(sequences)} sequences from FASTA file")
            return sequences
            
        except Exception as e:
            logger.error(f"Failed to read FASTA file: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Failed to read FASTA file: {e}")
    
    def mask_variants(self, sequence, variant_positions):
        """
        Simple variant masking - replace only the exact SNP positions with 'N'.
        
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
        logger.debug(f"Processing {len(variant_pos_list)} variant positions")
        
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
    
    def mask_sequences_for_primer_design(self, sequences, variants):
        """
        Mask variants in sequences for primer design.
        
        Args:
            sequences (dict): Dictionary of sequences
            variants (dict): Dictionary mapping chromosomes to sets of variant positions
            
        Returns:
            dict: Dictionary of masked sequences
        """
        logger.info("Masking variants in sequences...")
        masked_sequences = {}
        
        if Config.SHOW_PROGRESS:
            sequence_iter = tqdm(sequences.items(), desc="Masking sequences")
        else:
            sequence_iter = sequences.items()
        
        for seq_id, sequence in sequence_iter:
            # Try to find matching variants for this sequence
            # Assume sequence ID corresponds to chromosome name
            if seq_id in variants:
                variant_positions = variants[seq_id]
                masked_sequence = self.mask_variants(sequence, variant_positions)
                masked_sequences[seq_id] = masked_sequence
            else:
                # No variants found for this sequence, use original
                masked_sequences[seq_id] = sequence
                logger.debug(f"No variants found for sequence {seq_id}, using original")
        
        logger.info(f"Completed variant masking for {len(masked_sequences)} sequences")
        return masked_sequences
    
    def extract_variants_by_regions(self, vcf_file, conserved_regions):
        """
        Extract variants from VCF file only for specific conserved regions.
        This optimizes memory usage by only loading variants in regions of interest.
        
        Args:
            vcf_file (str): Path to VCF file
            conserved_regions (dict): Dictionary mapping chromosome names to lists of conserved regions
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        logger.debug(f"Extracting variants for specific regions from {vcf_file}")
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
            
        if not conserved_regions:
            logger.warning("No conserved regions provided")
            return {}
        
        # Prepare the VCF file (compress and index if needed)
        try:
            prepared_vcf = self.prepare_vcf_file(vcf_file)
            
            # Dictionary to store variants by chromosome
            variants = {}
            
            # Process each chromosome and its conserved regions
            for chrom, regions in conserved_regions.items():
                if not regions:
                    continue
                    
                variants[chrom] = set()
                logger.debug(f"Processing {len(regions)} conserved regions for {chrom}")
                
                # Process each conserved region
                for region in regions:
                    start = region['start']
                    end = region['end']
                    
                    # Skip very small regions
                    if end - start < 10:
                        logger.debug(f"Skipping small region {chrom}:{start}-{end} (< 10 bp)")
                        continue
                        
                    # Extract variants for this specific region
                    region_variants = self.get_region_variants(prepared_vcf, chrom, start, end)
                    variants[chrom].update(region_variants)
            
            # Count total variants
            total_variants = sum(len(positions) for positions in variants.values())
            logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes")
            
            return variants
            
        except Exception as e:
            logger.error(f"Failed to extract variants by regions: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(f"Failed to extract variants by regions: {e}")
    
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
        logger.info(f"Combining masking from two sets ({len(masked_sequences1)} and {len(masked_sequences2)} sequences)")
        
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