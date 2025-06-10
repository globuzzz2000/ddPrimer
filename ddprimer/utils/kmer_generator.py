#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K-mer generator module for ddPrimer pipeline.

Contains functionality for:
1. Memory-efficient k-mer counting from large genomes
2. Frequency analysis and filtering with config defaults
3. Output file generation with comprehensive headers
4. Batch processing with progress tracking

This module provides k-mer frequency analysis capabilities that can be used
for primer design optimization and genome characterization.
"""

import os
import time
import logging
from pathlib import Path
from typing import Dict, Iterator, Tuple, Optional
import mmap
from tqdm import tqdm

from ..config import Config, FileError, SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)


class KmerGenerator:
    """
    Memory-efficient k-mer generator for large genomes.
    
    This class provides methods for counting k-mers in DNA sequences using
    memory mapping for efficient processing of large genome files. Uses
    configuration defaults for k-mer sizes and frequency thresholds.
    
    Attributes:
        kmer_size: Length of k-mers to generate
        min_frequency: Minimum frequency threshold for output
        
    Example:
        >>> generator = KmerGenerator(kmer_size=11)
        >>> kmer_counts = generator.count_kmers_from_file("genome.fasta")
        >>> generator.write_kmer_file(kmer_counts, "output.list", "genome.fasta")
    """
    
    def __init__(self, kmer_size: int):
        """
        Initialize k-mer generator with specified k-mer size.
        
        Args:
            kmer_size: Length of k-mers to generate
            
        Raises:
            ValueError: If kmer_size is invalid
        """
        if kmer_size <= 0:
            raise ValueError(f"K-mer size must be positive, got {kmer_size}")
            
        self.kmer_size = kmer_size
        self.min_frequency = Config.KMER_MIN_FREQUENCY
        
        # Pre-compile translation tables for performance
        self.complement_trans = str.maketrans('ATGC', 'TACG')
        self.iupac_to_n_trans = str.maketrans('RYSWKMBDHV', 'NNNNNNNNNN')
        
        # Valid bases for k-mer counting
        self.valid_bases = set('ATGC')
    
    def count_kmers_from_file(self, fasta_file: str) -> Dict[str, int]:
        """
        Count k-mers directly from FASTA file using memory mapping.
        
        Efficiently processes large genome files by using memory mapping
        and canonical k-mer representation for reduced memory usage.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Dictionary mapping canonical k-mers to their frequencies
            
        Raises:
            FileError: If the FASTA file cannot be accessed
            SequenceProcessingError: If there's an error processing sequences
            
        Example:
            >>> counts = generator.count_kmers_from_file("genome.fasta")
            >>> print(f"Found {len(counts)} unique k-mers")
        """
        logger.info(f"Processing {Path(fasta_file).name} for {self.kmer_size}-mers")
        
        # Validate file exists
        if not os.path.exists(fasta_file):
            error_msg = f"FASTA file not found: {fasta_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        # Dictionary to store k-mer counts
        kmer_counts = {}
        total_kmers_processed = 0
        skipped_kmers = 0
        
        start_time = time.time()
        
        try:
            with open(fasta_file, 'rb') as f:
                # Memory map the file for efficient reading
                with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mmapped_file:
                    sequence_data = mmapped_file.read().decode('utf-8', errors='ignore')
                    
                    # Process sequences
                    sequences = list(self._parse_fasta_string(sequence_data))
                    
                    if Config.SHOW_PROGRESS:
                        sequence_iter = tqdm(sequences, desc="Processing sequences")
                    else:
                        sequence_iter = sequences
                    
                    for seq_id, sequence in sequence_iter:
                        logger.debug(f"Processing sequence: {seq_id}")
                        
                        # Clean and convert sequence
                        clean_seq = self._clean_sequence(sequence)
                        
                        # Count k-mers in this sequence
                        seq_kmers_processed = 0
                        seq_kmers_skipped = 0
                        
                        # Create progress bar for k-mers if showing progress and sequence is large
                        seq_length = len(clean_seq)
                        if Config.SHOW_PROGRESS and seq_length > 1000000:  # Only for large sequences
                            kmer_iter = tqdm(
                                self._generate_kmers(clean_seq),
                                total=max(0, seq_length - self.kmer_size + 1),
                                desc=f"K-mers in {seq_id}",
                                leave=False
                            )
                        else:
                            kmer_iter = self._generate_kmers(clean_seq)
                        
                        for kmer, canonical_kmer in kmer_iter:
                            if canonical_kmer:  # Valid k-mer (no N's)
                                kmer_counts[canonical_kmer] = kmer_counts.get(canonical_kmer, 0) + 1
                                seq_kmers_processed += 1
                            else:
                                seq_kmers_skipped += 1
                        
                        total_kmers_processed += seq_kmers_processed
                        skipped_kmers += seq_kmers_skipped
                        
                        logger.debug(f"Sequence {seq_id}: {seq_kmers_processed:,} k-mers processed, {seq_kmers_skipped:,} skipped")
        
        except (OSError, IOError) as e:
            error_msg = f"Error reading FASTA file {fasta_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        except Exception as e:
            error_msg = f"Error processing sequences from {fasta_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
        
        # Filter by frequency and log results
        filtered_counts = {kmer: count for kmer, count in kmer_counts.items() 
                          if count >= self.min_frequency}
        
        total_time = time.time() - start_time
        logger.info(f"K-mer counting completed in {total_time:.1f}s")
        logger.info(f"Found {len(filtered_counts):,} unique {self.kmer_size}-mers ≥ frequency {self.min_frequency}")
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Total k-mers processed: {total_kmers_processed:,}")
            logger.debug(f"K-mers skipped (containing N): {skipped_kmers:,}")
            logger.debug(f"Unique k-mers before filtering: {len(kmer_counts):,}")
            logger.debug(f"Processing rate: {total_kmers_processed/total_time:.0f} k-mers/second")
        
        return filtered_counts
    
    def _parse_fasta_string(self, fasta_data: str) -> Iterator[Tuple[str, str]]:
        """
        Parse FASTA data string into sequences.
        
        Efficiently parses FASTA format data and yields sequence ID and
        sequence pairs for processing.
        
        Args:
            fasta_data: Complete FASTA file content as string
            
        Yields:
            Tuples of (sequence_id, sequence)
        """
        lines = fasta_data.split('\n')
        current_id = None
        current_seq = []
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Yield previous sequence
                if current_id and current_seq:
                    yield current_id, ''.join(current_seq)
                
                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Yield last sequence
        if current_id and current_seq:
            yield current_id, ''.join(current_seq)
    
    def _clean_sequence(self, sequence: str) -> str:
        """
        Clean sequence by converting IUPAC codes to N and uppercase.
        
        Prepares sequences for k-mer analysis by standardizing format
        and handling ambiguous nucleotides.
        
        Args:
            sequence: Raw DNA sequence string
            
        Returns:
            Cleaned and standardized sequence
        """
        # Convert to uppercase and handle IUPAC codes
        clean_seq = sequence.upper().translate(self.iupac_to_n_trans)
        return clean_seq
    
    def _generate_kmers(self, sequence: str) -> Iterator[Tuple[str, Optional[str]]]:
        """
        Generate k-mers from sequence with canonical representation.
        
        Yields k-mers and their canonical forms (lexicographically smaller
        of k-mer and its reverse complement) for reduced memory usage.
        
        Args:
            sequence: DNA sequence
            
        Yields:
            Tuples of (original_kmer, canonical_kmer or None if invalid)
        """
        for i in range(len(sequence) - self.kmer_size + 1):
            kmer = sequence[i:i + self.kmer_size]
            
            # Check if k-mer is valid (no N's)
            if set(kmer).issubset(self.valid_bases):
                # Get canonical k-mer (lexicographically smaller of kmer and reverse complement)
                rc_kmer = kmer.translate(self.complement_trans)[::-1]
                canonical_kmer = min(kmer, rc_kmer)
                yield kmer, canonical_kmer
            else:
                yield kmer, None
    
    def write_kmer_file(self, kmer_counts: Dict[str, int], output_file: str, 
                       genome_file: str, species_name: str):
        """
        Write k-mer counts to output file with comprehensive headers.
        
        Creates a formatted k-mer frequency file with metadata headers
        and sorted k-mer entries for analysis and visualization.
        
        Args:
            kmer_counts: Dictionary of k-mer frequencies
            output_file: Path to output file
            genome_file: Source genome file path
            species_name: Species name for headers
            
        Raises:
            FileError: If output file cannot be written
            
        Example:
            >>> generator.write_kmer_file(counts, "output.list", "genome.fasta", "A.thaliana")
        """
        logger.debug(f"Writing {len(kmer_counts):,} k-mers to {output_file}")
        
        total_count = sum(kmer_counts.values())
        
        try:
            with open(output_file, 'w') as f:
                # Write comprehensive header
                f.write(f"# K-mer frequency list generated by ddPrimer\n")
                f.write(f"# Source: {Path(genome_file).name}\n")
                f.write(f"# Species: {species_name}\n")
                f.write(f"# K-mer size: {self.kmer_size}\n")
                f.write(f"# Minimum frequency: {self.min_frequency}\n")
                f.write(f"# Unique k-mers: {len(kmer_counts):,}\n")
                f.write(f"# Total occurrences: {total_count:,}\n")
                f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"# Format: kmer\\tcount\\tfrequency\n")
                f.write(f"#\n")
                
                # Write k-mers sorted by frequency (descending)
                if Config.SHOW_PROGRESS:
                    sorted_items = tqdm(
                        sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True),
                        desc="Writing k-mers"
                    )
                else:
                    sorted_items = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
                
                for kmer, count in sorted_items:
                    frequency = count / total_count if total_count > 0 else 0
                    f.write(f"{kmer}\t{count}\t{frequency:.6e}\n")
                    
        except (OSError, IOError) as e:
            error_msg = f"Error writing k-mer file {output_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error writing k-mer file {output_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
        logger.debug(f"Successfully wrote k-mer file: {output_file}")


def run_kmer_generation(fasta_file: str, output_dir: Optional[str] = None) -> bool:
    """
    Run k-mer generation workflow with config defaults.
    
    Orchestrates the complete k-mer generation process using configuration
    defaults for k-mer sizes and frequency thresholds.
    
    Args:
        fasta_file: Path to input FASTA file
        output_dir: Directory for output files
        
    Returns:
        True if generation completed successfully, False otherwise
        
    Raises:
        FileError: If files cannot be accessed
        SequenceProcessingError: If k-mer processing fails
        
    Example:
        >>> success = run_kmer_generation("genome.fasta", "./output")
    """
    logger.debug("=== K-MER GENERATION WORKFLOW DEBUG ===")
    
    # Validate input file
    if not os.path.exists(fasta_file):
        error_msg = f"FASTA file not found: {fasta_file}"
        logger.error(error_msg)
        raise FileError(error_msg)
    
    # Create output directory - default to ~/.ddprimer/kmer_lists
    if not output_dir:
        output_dir = os.path.join(os.path.expanduser("~"), ".ddprimer", "kmer_lists")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        logger.debug(f"Output directory: {output_dir}")
    except OSError as e:
        error_msg = f"Failed to create output directory {output_dir}: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise FileError(error_msg) from e
    
    # Determine species name from file
    species_name = Path(fasta_file).stem.replace("_", " ")
    
    logger.info(f"Starting k-mer generation for {species_name}")
    logger.info(f"K-mer sizes: {Config.KMER_SIZES}")
    logger.info(f"Minimum frequency: {Config.KMER_MIN_FREQUENCY}")
    
    total_start_time = time.time()
    
    try:
        # Process each k-mer size from config
        for kmer_size in Config.KMER_SIZES:
            logger.info(f"\nProcessing {kmer_size}-mers...")
            
            # Initialize generator
            generator = KmerGenerator(kmer_size)
            
            # Count k-mers
            kmer_counts = generator.count_kmers_from_file(fasta_file)
            
            if not kmer_counts:
                logger.warning(f"No k-mers found for size {kmer_size}")
                continue
            
            # Generate output filename
            output_file = os.path.join(output_dir, f"{Path(fasta_file).stem}_{kmer_size}.list")
            
            # Write output
            generator.write_kmer_file(kmer_counts, output_file, fasta_file, species_name)
            
            logger.info(f"Wrote {len(kmer_counts):,} {kmer_size}-mers to {Path(output_file).name}")
        
        total_time = time.time() - total_start_time
        logger.info(f"\nK-mer generation completed successfully in {total_time:.1f}s")
        logger.info(f"Output files in: {output_dir}")
        
        logger.debug("=== END K-MER GENERATION WORKFLOW DEBUG ===")
        return True
        
    except (FileError, SequenceProcessingError):
        # Re-raise specific exceptions
        logger.debug("=== END K-MER GENERATION WORKFLOW DEBUG ===")
        raise
    except Exception as e:
        error_msg = f"Unexpected error in k-mer generation: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        logger.debug("=== END K-MER GENERATION WORKFLOW DEBUG ===")
        raise SequenceProcessingError(error_msg) from e