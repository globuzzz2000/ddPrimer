#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct Mode Implementation for ddPrimer Pipeline

Handles target-sequence based primer design workflow using CSV/Excel input.
This mode bypasses SNP masking and annotation processing, directly proceeding
from sequence input to restriction cutting and primer design.

Contains functionality for:
1. Flexible sequence table loading with automatic column detection
2. Sequence validation and preprocessing
3. Integration with existing pipeline components
4. Workflow orchestration for direct mode execution

This module provides a streamlined workflow for users who have pre-processed
target sequences and want to skip genome-based variant processing.
"""

import os
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Import package modules
from ..config import FileError, FileFormatError, SequenceProcessingError, PrimerDesignError
from ..utils import FileIO

# Set up module logger
logger = logging.getLogger(__name__)


class DirectModeProcessor:
    """
    Handles direct mode workflow execution with sequence table input.
    
    This processor manages the simplified workflow for direct sequence input,
    bypassing genome-based processing while maintaining compatibility with
    the existing pipeline architecture.
    
    Example:
        >>> processor = DirectModeProcessor()
        >>> sequences = processor.load_sequences_from_table("sequences.csv")
        >>> success = processor.run_direct_workflow(sequences, output_dir)
    """

    # Common column name patterns for flexible detection
    SEQUENCE_ID_PATTERNS = [
        'id', 'name', 'gene', 'target', 'sequence_id', 'seq_id', 
        'gene_name', 'target_name', 'identifier', 'label'
    ]
    
    SEQUENCE_PATTERNS = [
        'sequence', 'seq', 'dna', 'nucleotide', 'target_sequence',
        'dna_sequence', 'genomic_sequence', 'template'
    ]

    def __init__(self):
        """Initialize the DirectModeProcessor."""
        logger.debug("DirectModeProcessor initialized")

    @classmethod
    def detect_columns(cls, df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
        """
        Automatically detect sequence ID and sequence columns in DataFrame.
        
        Uses flexible pattern matching to identify the most likely columns
        containing sequence identifiers and DNA sequences.
        
        Args:
            df: Input DataFrame from CSV/Excel file
            
        Returns:
            Tuple of (id_column_name, sequence_column_name) or (None, None) if not found
        """
        logger.debug("Detecting sequence columns in DataFrame")
        
        if df.empty:
            logger.warning("DataFrame is empty")
            return None, None
        
        columns = [col.lower().strip() for col in df.columns]
        original_columns = list(df.columns)
        
        # Create mapping from lowercase to original column names
        col_mapping = {col.lower().strip(): original for col, original in zip(columns, original_columns)}
        
        # Detect sequence ID column
        id_column = None
        id_column_lower = None
        for pattern in cls.SEQUENCE_ID_PATTERNS:
            for col in columns:
                if pattern in col:
                    id_column = col_mapping[col]
                    id_column_lower = col
                    logger.debug(f"Detected ID column: '{id_column}' (matched pattern: '{pattern}')")
                    break
            if id_column:
                break
        
        # Detect sequence column - make sure it's different from ID column
        seq_column = None
        for pattern in cls.SEQUENCE_PATTERNS:
            for col in columns:
                # Skip if this is the same column we already identified as ID
                if col == id_column_lower:
                    logger.debug(f"Skipping column '{col_mapping[col]}' for sequence detection (already used as ID)")
                    continue
                    
                if pattern in col:
                    seq_column = col_mapping[col]
                    logger.debug(f"Detected sequence column: '{seq_column}' (matched pattern: '{pattern}')")
                    break
            if seq_column:
                break
        
        # Fallback: use first two columns if detection fails and they're different
        if not id_column and len(original_columns) >= 1:
            id_column = original_columns[0]
            logger.debug(f"Fallback: Using first column as ID: '{id_column}'")
            
        if not seq_column and len(original_columns) >= 2:
            # Make sure sequence column is different from ID column
            candidate_seq_col = original_columns[1]
            if candidate_seq_col != id_column:
                seq_column = candidate_seq_col
                logger.debug(f"Fallback: Using second column as sequence: '{seq_column}'")
            elif len(original_columns) >= 3:
                # Try third column if second is same as first
                seq_column = original_columns[2]
                logger.debug(f"Fallback: Using third column as sequence: '{seq_column}'")
        elif not seq_column and len(original_columns) == 1:
            # If only one column, assume it contains sequences and generate IDs
            seq_column = original_columns[0]
            id_column = None  # Will generate automatic IDs
            logger.debug(f"Single column detected, treating as sequences: '{seq_column}'")
        
        # Final validation - make sure we don't have the same column for both
        if id_column and seq_column and id_column == seq_column:
            logger.warning(f"ID and sequence columns are the same: '{id_column}'. Using fallback logic.")
            if len(original_columns) >= 2:
                id_column = original_columns[0]
                seq_column = original_columns[1]
                logger.debug(f"Fallback correction: ID='{id_column}', Sequence='{seq_column}'")
            else:
                # Only one column available
                seq_column = original_columns[0]
                id_column = None
                logger.debug(f"Fallback correction: Auto-generate IDs, Sequence='{seq_column}'")
        
        return id_column, seq_column

    @classmethod
    def load_sequences_from_table(cls, file_path: str) -> Dict[str, str]:
        """
        Load sequences from CSV or Excel file with flexible column detection.
        
        Automatically detects sequence ID and sequence columns, validates
        sequences, and returns a dictionary compatible with the pipeline.
        
        Args:
            file_path: Path to CSV or Excel file containing sequences
            
        Returns:
            Dictionary mapping sequence IDs to sequences
            
        Raises:
            FileError: If file cannot be accessed
            FileFormatError: If file format is invalid or columns cannot be detected
            SequenceProcessingError: If sequences are invalid
        """
        logger.debug("=== DIRECT MODE: SEQUENCE TABLE LOADING ===")
        logger.debug(f"Loading sequences from table: {file_path}")
        
        if not os.path.exists(file_path):
            error_msg = f"Sequence table file not found: {file_path}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        try:
            # Load DataFrame based on file extension
            file_ext = Path(file_path).suffix.lower()
            
            if file_ext == '.csv':
                # Try different encodings and separators for CSV
                for encoding in ['utf-8', 'latin-1', 'cp1252']:
                    for sep in [',', ';', '\t']:
                        try:
                            df = pd.read_csv(file_path, encoding=encoding, sep=sep)
                            if len(df.columns) > 1 or (len(df.columns) == 1 and len(df) > 0):
                                logger.debug(f"Successfully loaded CSV with encoding={encoding}, sep='{sep}'")
                                break
                        except Exception as e:
                            logger.debug(f"Failed to load CSV with encoding={encoding}, sep='{sep}': {str(e)}")
                            continue
                    else:
                        continue
                    break
                else:
                    error_msg = f"Could not parse CSV file with any encoding/separator combination: {file_path}"
                    logger.error(error_msg)
                    raise FileFormatError(error_msg)
                    
            elif file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path)
                logger.debug(f"Successfully loaded Excel file")
            else:
                error_msg = f"Unsupported file format: {file_ext}. Supported formats: .csv, .xlsx, .xls"
                logger.error(error_msg)
                raise FileFormatError(error_msg)
            
            if df.empty:
                error_msg = f"No data found in file: {file_path}"
                logger.error(error_msg)
                raise FileFormatError(error_msg)
            
            logger.debug(f"Loaded DataFrame with {len(df)} rows and {len(df.columns)} columns")
            logger.debug(f"Columns: {list(df.columns)}")
            
            # Detect sequence columns
            id_column, seq_column = cls.detect_columns(df)
            
            if not seq_column:
                error_msg = (f"Could not detect sequence column in file: {file_path}. "
                           f"Available columns: {list(df.columns)}")
                logger.error(error_msg)
                raise FileFormatError(error_msg)
            
            logger.debug(f"Using columns - ID: '{id_column}', Sequence: '{seq_column}'")
            
            # Extract sequences
            sequences = {}
            
            for idx, row in df.iterrows():
                # Get sequence ID
                if id_column:
                    seq_id = str(row[id_column]).strip()
                    if not seq_id or seq_id.lower() in ['nan', 'none', '']:
                        seq_id = f"Sequence_{idx + 1}"
                else:
                    seq_id = f"Sequence_{idx + 1}"
                
                # Get sequence
                sequence = str(row[seq_column]).strip().upper()
                
                # Skip empty sequences
                if not sequence or sequence.lower() in ['nan', 'none', '']:
                    logger.warning(f"Skipping empty sequence for ID: {seq_id}")
                    continue
                
                # Basic sequence validation
                if not cls.validate_sequence(sequence):
                    logger.warning(f"Skipping invalid sequence for ID: {seq_id} (contains non-DNA characters)")
                    continue
                
                # Handle duplicate IDs
                original_id = seq_id
                counter = 1
                while seq_id in sequences:
                    seq_id = f"{original_id}_{counter}"
                    counter += 1
                
                sequences[seq_id] = sequence
                logger.debug(f"Added sequence: {seq_id} ({len(sequence)} bp)")
            
            if not sequences:
                error_msg = f"No valid sequences found in file: {file_path}"
                logger.error(error_msg)
                raise SequenceProcessingError(error_msg)
            
            logger.debug(f"Successfully loaded {len(sequences)} sequences from table")
            logger.debug("=== END DIRECT MODE: SEQUENCE TABLE LOADING ===")
            
            return sequences
            
        except (FileError, FileFormatError, SequenceProcessingError):
            raise
        except Exception as e:
            error_msg = f"Error loading sequences from table {file_path}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END DIRECT MODE: SEQUENCE TABLE LOADING ===")
            raise FileFormatError(error_msg) from e

    @staticmethod
    def validate_sequence(sequence: str) -> bool:
        """
        Validate that a sequence contains only valid DNA characters.
        
        Uses a simple but robust approach that handles various input formats
        and provides clear validation.
        
        Args:
            sequence: DNA sequence string to validate
            
        Returns:
            True if sequence is valid, False otherwise
        """
        if not sequence:
            return False
        
        try:
            # Convert to string, strip whitespace, and convert to uppercase
            sequence_str = str(sequence).strip().upper()
            
            if not sequence_str:
                return False
            
            # Remove any internal whitespace
            sequence_clean = ''.join(sequence_str.split())
            
            if not sequence_clean:
                return False
            
            # Check minimum length (at least 10 bp for meaningful primers)
            if len(sequence_clean) < 10:
                logger.debug(f"Sequence too short: {len(sequence_clean)} bp")
                return False
            
            # Simple character validation - only allow standard DNA bases
            # Start with strict validation, can be relaxed if needed
            allowed_chars = set('ATCG')
            sequence_chars = set(sequence_clean)
            
            # Check if all characters are valid
            if not sequence_chars.issubset(allowed_chars):
                invalid_chars = sequence_chars - allowed_chars
                logger.debug(f"Invalid DNA characters found: {sorted(invalid_chars)}")
                
                # Allow some common ambiguous bases if present
                extended_chars = set('ATCGRYSWKMBDHVN')
                if sequence_chars.issubset(extended_chars):
                    logger.debug("Sequence contains ambiguous bases but is acceptable")
                    return True
                else:
                    return False
            
            return True
            
        except Exception as e:
            logger.debug(f"Error in sequence validation: {str(e)}")
            return False

    @classmethod
    def run_direct_workflow(cls, sequences: Dict[str, str], output_dir: str, 
                           input_file: str, enable_internal_oligo: bool = True,
                           enable_snp_masking: bool = False, vcf_file: str = None, 
                           fasta_file: str = None) -> bool:
        """
        Execute the direct mode primer design workflow.
        
        Runs the streamlined workflow that bypasses SNP masking and annotation
        processing, proceeding directly from input sequences to primer design.
        Optionally supports SNP masking if enabled.
        
        Args:
            sequences: Dictionary of sequence IDs to sequences
            output_dir: Output directory path
            input_file: Original input file path (for naming output)
            enable_internal_oligo: Whether to design internal oligos (probes)
            enable_snp_masking: Whether to apply SNP masking to sequences
            vcf_file: Path to VCF file (required if enable_snp_masking=True)
            fasta_file: Path to FASTA file (required if enable_snp_masking=True)
            
        Returns:
            True if workflow completed successfully, False otherwise
            
        Raises:
            SequenceProcessingError: If sequence processing fails
            PrimerDesignError: If primer design fails
        """
        logger.debug("=== DIRECT MODE WORKFLOW: PRIMER DESIGN PIPELINE ===")
        logger.debug(f"Processing {len(sequences)} input sequences")
        logger.debug(f"Internal oligo design: {enable_internal_oligo}")
        logger.debug(f"SNP masking enabled: {enable_snp_masking}")
        
        # Import required processors
        from ..core import (SequenceProcessor, Primer3Processor, 
                           PrimerProcessor, ThermoProcessor, BlastProcessor)
        from ..config import Config
        
        try:
            processed_sequences = sequences
            
            # Step 1: Apply SNP masking if enabled
            if enable_snp_masking:
                if not vcf_file or not fasta_file:
                    error_msg = "SNP masking requires both VCF and FASTA files"
                    logger.error(error_msg)
                    logger.debug("=== END DIRECT MODE WORKFLOW ===")
                    return False
                
                logger.debug("DIRECT: Applying SNP masking to sequences")
                processed_sequences = cls._apply_snp_masking_to_sequences(
                    sequences, vcf_file, fasta_file
                )
                
                if not processed_sequences:
                    logger.warning("No sequences remained after SNP masking. Exiting.")
                    logger.debug("=== END DIRECT MODE WORKFLOW ===")
                    return False
            
            # Step 2: Process restriction sites (reuse existing functionality)
            logger.debug("DIRECT: Processing restriction sites")
            restriction_fragments = SequenceProcessor.process_restriction_sites_workflow(processed_sequences)
            
            if not restriction_fragments:
                logger.warning("No valid fragments after restriction site processing. Exiting.")
                logger.debug("=== END DIRECT MODE WORKFLOW ===")
                return False
            
            # Step 3: Skip gene overlap filtering (direct mode)
            logger.debug("DIRECT: Skipping gene overlap filtering (direct mode)")
            filtered_fragments = restriction_fragments
            
            # Step 4: Design primers with Primer3
            logger.debug("DIRECT: Designing primers with Primer3")
            primer3_processor = Primer3Processor(Config, enable_internal_oligo=enable_internal_oligo)
            primer_records = primer3_processor.design_primers_workflow(filtered_fragments)
            
            if not primer_records:
                logger.warning("No primers were designed by Primer3. Exiting.")
                logger.debug("=== END DIRECT MODE WORKFLOW ===")
                return False
            
            # Step 5: Filter primers
            logger.debug("DIRECT: Filtering primers")
            df = PrimerProcessor.filter_primers_workflow(primer_records)
            
            if df is None or len(df) == 0:
                logger.warning("No primers passed filtering criteria. Exiting.")
                logger.debug("=== END DIRECT MODE WORKFLOW ===")
                return False
            
            # Step 6: Calculate thermodynamic properties
            logger.debug("DIRECT: Calculating thermodynamic properties")
            df = ThermoProcessor.calculate_thermodynamics_workflow(df)
            
            # Step 7: Run BLAST for specificity
            logger.debug("DIRECT: Running BLAST specificity checks")
            df = BlastProcessor.run_blast_specificity_workflow(df)
            
            if df is None or len(df) == 0:
                logger.warning("No primers passed BLAST filtering. Exiting.")
                logger.debug("=== END DIRECT MODE WORKFLOW ===")
                return False
            
            # Step 8: Save results to Excel file (using direct mode)
            logger.debug("DIRECT: Saving results")
            output_path = FileIO.save_results(
                df, 
                output_dir, 
                input_file, 
                mode='direct'  # Use direct mode for output formatting
            )
            
            if output_path:
                logger.info(f"\nResults saved to: {output_path}")
                logger.debug("=== END DIRECT MODE WORKFLOW ===")
                return True
            else:
                logger.error("Failed to save results.")
                logger.debug("=== END DIRECT MODE WORKFLOW ===")
                return False
                
        except SequenceProcessingError as e:
            error_msg = f"Sequence processing error in direct mode: {str(e)}"
            logger.error(error_msg)
            logger.debug("=== END DIRECT MODE WORKFLOW ===")
            return False
        except PrimerDesignError as e:
            error_msg = f"Primer design error in direct mode: {str(e)}"
            logger.error(error_msg)
            logger.debug("=== END DIRECT MODE WORKFLOW ===")
            return False
        except Exception as e:
            error_msg = f"Error in direct mode workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END DIRECT MODE WORKFLOW ===")
            return False

    @classmethod
    def _apply_snp_masking_to_sequences(cls, sequences: Dict[str, str], 
                                       vcf_file: str, fasta_file: str) -> Dict[str, str]:
        """
        Apply SNP masking to direct mode sequences by mapping them to reference.
        
        This method attempts to identify where each input sequence maps to the
        reference genome, then applies VCF-based variant masking to those regions.
        
        Args:
            sequences: Dictionary of sequence IDs to sequences
            vcf_file: Path to VCF file with variants
            fasta_file: Path to reference FASTA file
            
        Returns:
            Dictionary of masked sequences
            
        Raises:
            SequenceProcessingError: If SNP masking fails
        """
        logger.debug("=== DIRECT MODE: SNP MASKING ===")
        logger.debug(f"Applying SNP masking to {len(sequences)} sequences")
        
        try:
            from ..core import SNPMaskingProcessor
            from ..utils import FileIO
            
            # Load reference sequences
            logger.debug("Loading reference genome")
            reference_sequences = FileIO.load_fasta(fasta_file)
            
            masked_sequences = {}
            mapping_stats = {"mapped": 0, "unmapped": 0, "total_masked": 0, "total_substituted": 0}
            
            for seq_id, sequence in sequences.items():
                try:
                    # Find the best match for this sequence in the reference
                    match_info = cls._find_sequence_in_reference(
                        sequence, reference_sequences, seq_id
                    )
                    
                    if not match_info:
                        logger.debug(f"Could not map sequence {seq_id} to reference - discarding")
                        mapping_stats["unmapped"] += 1
                        continue  # Skip unmapped sequences
                    
                    ref_chrom = match_info['chromosome']
                    ref_start = match_info['start']
                    ref_end = match_info['end']
                    
                    logger.debug(f"Mapped {seq_id} to {ref_chrom}:{ref_start}-{ref_end}")
                    mapping_stats["mapped"] += 1
                    
                    # Extract the reference region
                    ref_sequence = reference_sequences[ref_chrom][ref_start:ref_end]
                    
                    # Apply SNP masking to the reference region
                    processor = SNPMaskingProcessor(fasta_file)
                    masked_ref_sequence, seq_stats = processor.process_sequence_with_vcf(
                        sequence=ref_sequence,
                        vcf_path=vcf_file,
                        chromosome=ref_chrom,
                        return_stats=True
                    )
                    
                    masked_sequences[seq_id] = masked_ref_sequence
                    mapping_stats["total_masked"] += seq_stats.get("masked", 0)
                    mapping_stats["total_substituted"] += seq_stats.get("substituted", 0)
                    logger.debug(f"Applied SNP masking to {seq_id}")
                    
                except Exception as e:
                    logger.warning(f"Error masking sequence {seq_id}: {str(e)}")
                    mapping_stats["unmapped"] += 1
                    continue  # Skip sequences that error during processing
            
            # SNP masking summary logs
            if mapping_stats["mapped"] > 0:
                logger.info(f"Successfully mapped and processed {mapping_stats['mapped']} sequences")
            if mapping_stats["unmapped"] > 0:
                logger.info(f"{mapping_stats['unmapped']} sequences could not be mapped to reference and were discarded")
            
            if mapping_stats["total_masked"] > 0 or mapping_stats["total_substituted"] > 0:
                logger.info(f"Masked {mapping_stats['total_masked']} variable variants, substituted {mapping_stats['total_substituted']} fixed variants")
            
            logger.debug(f"SNP masking complete: {len(masked_sequences)} sequences processed")
            logger.debug("=== END DIRECT MODE: SNP MASKING ===")

            return masked_sequences
            
        except Exception as e:
            error_msg = f"Error in SNP masking for direct mode: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END DIRECT MODE: SNP MASKING ===")
            raise SequenceProcessingError(error_msg) from e

    @classmethod
    def _find_sequence_in_reference(cls, query_sequence: str, 
                                   reference_sequences: Dict[str, str], 
                                   seq_id: str) -> Optional[Dict]:
        """
        Find where a query sequence maps to in the reference genome.
        
        Uses exact string matching to locate the query sequence within the
        reference genome. This is suitable for sequences that are extracted
        directly from the reference.
        
        Args:
            query_sequence: The sequence to find
            reference_sequences: Dictionary of reference chromosome sequences
            seq_id: Sequence identifier for logging
            
        Returns:
            Dictionary with mapping info: {'chromosome': str, 'start': int, 'end': int}
            or None if no match found
        """
        logger.debug(f"Searching for {seq_id} ({len(query_sequence)} bp) in reference")
        
        # Convert query to uppercase for matching
        query_upper = query_sequence.upper()
        
        # Search each chromosome
        for chrom_name, chrom_sequence in reference_sequences.items():
            chrom_upper = chrom_sequence.upper()
            
            # Look for exact match
            pos = chrom_upper.find(query_upper)
            if pos != -1:
                logger.debug(f"Found exact match for {seq_id} on {chrom_name} at position {pos}")
                return {
                    'chromosome': chrom_name,
                    'start': pos,
                    'end': pos + len(query_sequence)
                }
            
            # Also try reverse complement
            reverse_complement = cls._reverse_complement(query_upper)
            pos = chrom_upper.find(reverse_complement)
            if pos != -1:
                logger.debug(f"Found reverse complement match for {seq_id} on {chrom_name} at position {pos}")
                return {
                    'chromosome': chrom_name,
                    'start': pos,
                    'end': pos + len(reverse_complement)
                }
        
        logger.debug(f"No exact match found for {seq_id} in reference genome")
        
        # Try partial matching (at least 80% of sequence length)
        min_match_length = max(50, int(len(query_sequence) * 0.8))
        
        for chrom_name, chrom_sequence in reference_sequences.items():
            chrom_upper = chrom_sequence.upper()
            
            # Sliding window search for partial matches
            for i in range(len(query_upper) - min_match_length + 1):
                query_fragment = query_upper[i:i + min_match_length]
                pos = chrom_upper.find(query_fragment)
                
                if pos != -1:
                    # Found a partial match, extend it
                    logger.debug(f"Found partial match for {seq_id} on {chrom_name} at position {pos}")
                    
                    # Use the original query sequence coordinates but map to reference position
                    ref_start = max(0, pos - i)
                    ref_end = min(len(chrom_sequence), ref_start + len(query_sequence))
                    
                    return {
                        'chromosome': chrom_name,
                        'start': ref_start,
                        'end': ref_end
                    }
        
        return None

    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """
        Generate reverse complement of a DNA sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Reverse complement sequence
        """
        complement_map = {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D', 'N': 'N', '-': '-'
        }
        
        return ''.join(complement_map.get(base, base) for base in reversed(sequence))


#############################################################################
#                         Direct Mode Main Execution
#############################################################################

def run_direct_mode(args):
    """
    Execute the direct mode workflow with input file handling.
    
    Handles file selection, sequence loading, and workflow execution
    for the direct mode pipeline. Supports optional SNP masking.
    
    Args:
        args: Parsed command line arguments
        
    Returns:
        True if workflow completed successfully, False otherwise
    """
    logger.info("=== Direct Mode Workflow ===")
    logger.debug("=== MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
    
    try:
        # Handle input file selection
        input_file = None
        
        # Check if direct mode was called with a file argument
        if isinstance(args.direct, str):
            # Direct mode was called with a file path: --direct file.csv
            input_file = args.direct
            if not os.path.exists(input_file):
                error_msg = f"Direct mode input file not found: {input_file}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
                return False
        else:
            # Direct mode was called without file: --direct
            # Prompt for file selection
            logger.info("\n>>> Please select sequence table file (CSV/Excel) <<<")
            try:
                input_file = FileIO.select_file(
                    "Select sequence table file",
                    [
                        ("CSV Files", "*.csv"),
                        ("Excel Files", "*.xlsx"),
                        ("Excel Files", "*.xls"),
                        ("All Files", "*.*")
                    ]
                )
            except Exception as e:
                error_msg = f"File selection failed: {str(e)}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
                return False
        
        # Handle SNP masking file selection if enabled
        vcf_file = None
        fasta_file = None
        enable_snp_masking = hasattr(args, 'snp') and args.snp
        
        if enable_snp_masking:
            # Get VCF file if not provided in args
            if not hasattr(args, 'vcf') or not args.vcf:
                logger.info("\n>>> Please select VCF variant file <<<")
                try:
                    vcf_file = FileIO.select_file(
                        "Select VCF variant file", 
                        [
                            ("VCF Files", "*.vcf"),
                            ("Compressed VCF Files", "*.vcf.gz"),
                            ("All Files", "*.*")
                        ]
                    )
                except Exception as e:
                    error_msg = f"VCF file selection failed: {str(e)}"
                    logger.error(error_msg)
                    logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
                    return False
            else:
                vcf_file = args.vcf
            
            # Get FASTA file if not provided in args
            if not hasattr(args, 'fasta') or not args.fasta:
                logger.info("\n>>> Please select reference FASTA file <<<")
                try:
                    fasta_file = FileIO.select_file(
                        "Select reference FASTA file",
                        [
                            ("FASTA Files", "*.fasta"), 
                            ("FASTA Files", "*.fa"), 
                            ("FASTA Files", "*.fna"), 
                            ("All Files", "*")]
                    )
                except Exception as e:
                    error_msg = f"FASTA file selection failed: {str(e)}"
                    logger.error(error_msg)
                    logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
                    return False
            else:
                fasta_file = args.fasta
            
            # Prepare VCF and FASTA files for compatibility
            logger.info("\nPreparing files for SNP masking...")
            try:
                from ..utils import FilePreparator
                
                logger.debug("DIRECT: Delegating file preparation for SNP masking")
                prep_result = FilePreparator.prepare_pipeline_files_workflow(
                    vcf_file=vcf_file,
                    fasta_file=fasta_file,
                    gff_file=None  # No GFF needed for direct mode
                )
                
                if not prep_result['success']:
                    error_msg = f"File preparation failed: {prep_result.get('reason', 'Unknown error')}"
                    logger.error(error_msg)
                    logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
                    return False
                
                # Update file paths to use prepared files
                if prep_result.get('changes_made', False):
                    vcf_file = prep_result['vcf_file']
                    fasta_file = prep_result['fasta_file']
                    logger.debug("Files successfully prepared for SNP masking")
                else:
                    logger.debug("Files are compatible and ready for SNP masking")
                
            except Exception as e:
                error_msg = f"File preparation failed: {str(e)}"
                logger.error(error_msg)
                logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
                return False
        
        # Signal that file selection is complete
        FileIO.mark_selection_complete()
        
        # Load sequences from the table
        logger.debug(f"\nLoading sequences from: {input_file}")
        try:
            processor = DirectModeProcessor()
            sequences = processor.load_sequences_from_table(input_file)
            logger.info(f"\nLoaded {len(sequences)} sequences from table")
        except Exception as e:
            error_msg = f"Failed to load sequences from table: {str(e)}"
            logger.error(error_msg)
            logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
            return False
        
        # Set up output directory
        if hasattr(args, 'output') and args.output:
            output_dir = args.output
        else:
            # Use the directory of the input file
            input_dir = os.path.dirname(os.path.abspath(input_file))
            output_dir = os.path.join(input_dir, "Primers")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Determine internal oligo setting from args
        enable_internal_oligo = not (hasattr(args, 'nooligo') and args.nooligo)
        
        # Run the direct workflow with optional SNP masking
        logger.debug("DIRECT: Delegating to direct mode workflow")
        processor = DirectModeProcessor()
        success = processor.run_direct_workflow(
            sequences=sequences,
            output_dir=output_dir,
            input_file=input_file,
            enable_internal_oligo=enable_internal_oligo,
            enable_snp_masking=enable_snp_masking,
            vcf_file=vcf_file,
            fasta_file=fasta_file
        )
        
        logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
        return success
        
    except Exception as e:
        error_msg = f"Error in direct mode execution: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        logger.debug("=== END MAIN WORKFLOW: DIRECT MODE EXECUTION ===")
        return False