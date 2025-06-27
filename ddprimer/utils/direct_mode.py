#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct mode processor for ddPrimer pipeline.

Contains functionality for:
1. Target sequence-based primer design workflow
2. CSV/Excel input file processing with sequence analysis
3. Optional VCF variant masking for direct sequences
4. Simplified fragment preparation bypassing restriction sites
5. Integration with existing core processors

This module provides an alternative workflow for primer design when users
have specific target sequences rather than genomic regions, enabling
direct primer design with optional variant consideration.
"""

import os
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union

from ..config import Config, SequenceProcessingError, PrimerDesignError, FileError
from .file_io import FileIO
from ..core import SNPMaskingProcessor, Primer3Processor, PrimerProcessor, ThermoProcessor, BlastProcessor

# Set up module logger
logger = logging.getLogger(__name__)

# Type alias for path inputs
PathLike = Union[str, Path]


class SequenceAnalyzer:
    """
    Analyzes sequence table files to identify sequence and name columns.
    
    This class provides intelligent analysis of CSV and Excel files containing
    sequence data, automatically detecting the most likely sequence and identifier
    columns based on content patterns and column names.
    
    Example:
        >>> analysis = SequenceAnalyzer.analyze_file("sequences.csv")
        >>> name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
    """
    
    @staticmethod
    def analyze_file(file_path: PathLike) -> Dict[str, Any]:
        """
        Analyze a sequence table file to identify column structure.
        
        Args:
            file_path: Path to CSV or Excel file
            
        Returns:
            Dictionary containing analysis results with column recommendations
            
        Raises:
            FileError: If file cannot be read
        """
        logger.debug(f"=== SEQUENCE FILE ANALYSIS ===")
        logger.debug(f"Analyzing file: {file_path}")
        
        try:
            # Load file based on extension
            if str(file_path).endswith('.csv'):
                df = pd.read_csv(file_path)
            elif str(file_path).endswith(('.xlsx', '.xls')):
                df = pd.read_excel(file_path)
            else:
                error_msg = f"Unsupported file format: {file_path}"
                logger.error(error_msg)
                return {"error": error_msg}
            
            # Remove completely empty columns and rows
            df = df.dropna(axis=1, how='all')
            df = df.dropna(axis=0, how='all')
            
            if df.empty:
                error_msg = "File contains no data"
                logger.error(error_msg)
                return {"error": error_msg}
            
            logger.debug(f"File contains {len(df)} rows and {len(df.columns)} columns")
            
            analysis = {
                "file_path": str(file_path),
                "num_rows": len(df),
                "num_columns": len(df.columns),
                "columns": list(df.columns),
                "column_analysis": {}
            }
            
            # Analyze each column
            for col in df.columns:
                col_analysis = SequenceAnalyzer._analyze_column(df[col], col)
                analysis["column_analysis"][col] = col_analysis
                logger.debug(f"Column '{col}': {col_analysis['type']} (confidence: {col_analysis['confidence']:.2f})")
            
            logger.debug("=== END SEQUENCE FILE ANALYSIS ===")
            return analysis
            
        except Exception as e:
            error_msg = f"Error analyzing file {file_path}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return {"error": error_msg}
    
    @staticmethod
    def _analyze_column(series: pd.Series, column_name: str) -> Dict[str, Any]:
        """
        Analyze a single column to determine its likely content type.
        
        Args:
            series: Pandas series to analyze
            column_name: Name of the column
            
        Returns:
            Dictionary with column analysis results
        """
        # Remove null values for analysis
        non_null_values = series.dropna()
        
        if len(non_null_values) == 0:
            return {
                "type": "empty",
                "confidence": 0.0,
                "sample_values": [],
                "avg_length": 0,
                "nucleotide_ratio": 0.0
            }
        
        # Convert to strings and clean
        str_values = [str(val).strip().upper() for val in non_null_values]
        sample_values = str_values[:3]  # First 3 values for display
        
        # Calculate basic statistics
        avg_length = sum(len(val) for val in str_values) / len(str_values)
        
        # Check for DNA sequence patterns
        nucleotide_count = 0
        total_chars = 0
        
        for val in str_values:
            for char in val:
                if char in 'ATCGN':
                    nucleotide_count += 1
                total_chars += 1
        
        nucleotide_ratio = nucleotide_count / total_chars if total_chars > 0 else 0.0
        
        # Determine column type and confidence
        type_scores = {
            "sequence": SequenceAnalyzer._score_sequence_column(column_name, str_values, nucleotide_ratio, avg_length),
            "name": SequenceAnalyzer._score_name_column(column_name, str_values, avg_length),
            "position": SequenceAnalyzer._score_position_column(column_name, str_values),
            "other": 0.1  # Base score for unrecognized columns
        }
        
        # Get the type with highest score
        best_type = max(type_scores.keys(), key=lambda k: type_scores[k])
        confidence = type_scores[best_type]
        
        return {
            "type": best_type,
            "confidence": confidence,
            "sample_values": sample_values,
            "avg_length": avg_length,
            "nucleotide_ratio": nucleotide_ratio,
            "type_scores": type_scores
        }
    
    @staticmethod
    def _score_sequence_column(name: str, values: List[str], nucleotide_ratio: float, avg_length: float) -> float:
        """Score how likely a column is to contain DNA sequences."""
        score = 0.0
        
        # Name-based scoring
        name_lower = name.lower()
        if 'seq' in name_lower or 'dna' in name_lower or 'sequence' in name_lower:
            score += 0.4
        elif 'target' in name_lower or 'template' in name_lower:
            score += 0.3
        
        # Content-based scoring
        if nucleotide_ratio > 0.8:  # >80% nucleotides
            score += 0.4
        elif nucleotide_ratio > 0.6:  # >60% nucleotides
            score += 0.2
        
        # Length-based scoring (typical sequences are 50-500 bp)
        if 50 <= avg_length <= 500:
            score += 0.2
        elif 20 <= avg_length <= 1000:
            score += 0.1
        
        return min(score, 1.0)
    
    @staticmethod
    def _score_name_column(name: str, values: List[str], avg_length: float) -> float:
        """Score how likely a column is to contain sequence identifiers."""
        score = 0.0
        
        # Name-based scoring
        name_lower = name.lower()
        if name_lower in ['id', 'name', 'identifier', 'gene', 'target']:
            score += 0.4
        elif 'id' in name_lower or 'name' in name_lower:
            score += 0.3
        
        # Content-based scoring (short, alphanumeric identifiers)
        if 3 <= avg_length <= 30:
            score += 0.3
        
        # Check if values look like identifiers
        identifier_patterns = 0
        for val in values[:10]:  # Check first 10 values
            if len(val) < 50 and ('_' in val or val.isalnum()):
                identifier_patterns += 1
        
        if identifier_patterns > len(values[:10]) * 0.7:  # >70% look like IDs
            score += 0.3
        
        return min(score, 1.0)
    
    @staticmethod
    def _score_position_column(name: str, values: List[str]) -> float:
        """Score position columns (not used in direct mode but kept for completeness)."""
        # Direct mode doesn't use position information, always return low score
        return 0.0
    
    @staticmethod
    def get_recommended_columns(analysis: Dict[str, Any]) -> Tuple[Optional[str], Optional[str]]:
        """
        Get recommended name and sequence columns from analysis.
        
        Args:
            analysis: Analysis results from analyze_file()
            
        Returns:
            Tuple of (name_column, sequence_column) - either may be None
        """
        if "error" in analysis:
            logger.error(f"Cannot get recommendations due to analysis error: {analysis['error']}")
            return None, None
        
        column_analysis = analysis.get("column_analysis", {})
        
        # Find best sequence column
        sequence_candidates = []
        name_candidates = []
        
        for col, col_data in column_analysis.items():
            if col_data["type"] == "sequence":
                sequence_candidates.append((col, col_data["confidence"]))
            elif col_data["type"] == "name":
                name_candidates.append((col, col_data["confidence"]))
        
        # Sort by confidence
        sequence_candidates.sort(key=lambda x: x[1], reverse=True)
        name_candidates.sort(key=lambda x: x[1], reverse=True)
        
        # Get best candidates
        best_sequence = sequence_candidates[0][0] if sequence_candidates else None
        best_name = name_candidates[0][0] if name_candidates else None
        
        # Log recommendations
        if best_sequence:
            seq_confidence = sequence_candidates[0][1]
            logger.debug(f"Recommended sequence column: '{best_sequence}' (confidence: {seq_confidence:.2f})")
        else:
            logger.warning("No sequence column identified")
        
        if best_name:
            name_confidence = name_candidates[0][1]
            logger.debug(f"Recommended name column: '{best_name}' (confidence: {name_confidence:.2f})")
        else:
            logger.debug("No name column identified - will generate sequence IDs")
        
        return best_name, best_sequence


class DirectMode:
    """
    Processes primer design workflow for direct sequence input mode.
    
    This class handles the complete direct mode workflow including sequence
    loading from tables, optional VCF processing, simplified fragment
    preparation, and integration with existing primer design components.
    
    Attributes:
        vcf_file: Optional VCF file for variant masking
        reference_file: Reference genome for VCF processing
        
    Example:
        >>> processor = DirectMode(vcf_file="variants.vcf", reference_file="genome.fasta")
        >>> success = processor.run_workflow("sequences.csv", "output_dir")
    """
    
    def __init__(self, vcf_file: Optional[PathLike] = None, reference_file: Optional[PathLike] = None):
        """
        Initialize DirectMode.
        
        Args:
            vcf_file: Optional VCF file for variant masking
            reference_file: Reference genome file (required if vcf_file provided)
        """
        self.vcf_file = vcf_file
        self.reference_file = reference_file
        
        # Validate VCF requirements
        if self.vcf_file and not self.reference_file:
            error_msg = "Reference genome file required when VCF file is provided"
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        logger.debug(f"DirectMode initialized with VCF={bool(vcf_file)}, Reference={bool(reference_file)}")
    
    def load_sequences_from_table(self, table_file: PathLike) -> Dict[str, str]:
        """
        Load sequences from CSV or Excel table using intelligent column detection.
        
        Args:
            table_file: Path to CSV or Excel file
            
        Returns:
            Dictionary mapping sequence IDs to sequences
            
        Raises:
            FileError: If file loading fails
            SequenceProcessingError: If sequence processing fails
        """
        logger.debug("=== SEQUENCE TABLE LOADING ===")
        logger.debug(f"Loading sequences from: {table_file}")
        
        try:
            sequences = FileIO.load_sequences_from_table(table_file)
            
            if not sequences:
                error_msg = "No valid sequences found in table file"
                logger.error(error_msg)
                raise SequenceProcessingError(error_msg)
            
            logger.debug(f"Successfully loaded {len(sequences)} sequences")
            logger.debug("=== END SEQUENCE TABLE LOADING ===")
            return sequences
            
        except Exception as e:
            error_msg = f"Failed to load sequences from table: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
    
    def process_sequences_with_vcf(self, sequences: Dict[str, str]) -> Dict[str, str]:
        """
        Apply VCF variants to sequences using sequence similarity matching.
        
        For each input sequence, finds the best match in the reference genome
        and applies variants from that region.
        
        Args:
            sequences: Dictionary of sequence ID to sequence
            
        Returns:
            Dictionary of processed sequences with variants applied
            
        Raises:
            SequenceProcessingError: If VCF processing fails
        """
        if not self.vcf_file:
            logger.debug("No VCF file provided, returning original sequences")
            return sequences
        
        logger.debug("=== VCF SEQUENCE PROCESSING ===")
        logger.debug(f"Processing {len(sequences)} sequences with VCF variants")
        
        try:
            # Initialize SNP processor
            snp_processor = SNPMaskingProcessor(self.reference_file)
            
            # Load reference genome for sequence matching
            reference_sequences = FileIO.load_fasta(self.reference_file)
            
            processed_sequences = {}
            
            for seq_id, sequence in sequences.items():
                try:
                    # Find best matching chromosome/region
                    best_match = self._find_best_sequence_match(sequence, reference_sequences)
                    
                    if not best_match:
                        logger.warning(f"No suitable reference match found for sequence {seq_id}, keeping original")
                        processed_sequences[seq_id] = sequence
                        continue
                    
                    chr_name, match_start, match_end = best_match
                    logger.debug(f"Sequence {seq_id} matched to {chr_name}:{match_start}-{match_end}")
                    
                    # Extract the matching reference region
                    ref_sequence = reference_sequences[chr_name][match_start:match_end]
                    
                    # Apply VCF variants to the reference region
                    modified_sequence = snp_processor.process_sequence_with_vcf(
                        sequence=ref_sequence,
                        vcf_path=self.vcf_file,
                        chromosome=chr_name,
                        start_offset=match_start
                    )
                    
                    processed_sequences[seq_id] = modified_sequence
                    
                except Exception as e:
                    logger.warning(f"Error processing sequence {seq_id} with VCF: {e}")
                    # Keep original sequence if processing fails
                    processed_sequences[seq_id] = sequence
            
            logger.debug(f"Successfully processed {len(processed_sequences)} sequences with VCF")
            logger.debug("=== END VCF SEQUENCE PROCESSING ===")
            return processed_sequences
            
        except Exception as e:
            error_msg = f"VCF processing failed: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
    
    def _find_best_sequence_match(self, query_sequence: str, reference_sequences: Dict[str, str]) -> Optional[Tuple[str, int, int]]:
        """
        Find the best matching region in reference genome for a query sequence.
        
        Uses sliding window approach to find the best match based on sequence similarity.
        
        Args:
            query_sequence: Query sequence to match
            reference_sequences: Dictionary of reference sequences
            
        Returns:
            Tuple of (chromosome, start, end) for best match, or None if no good match
        """
        logger.debug(f"=== SEQUENCE MATCHING ===")
        logger.debug(f"Finding reference match for sequence of length {len(query_sequence)}")
        
        query_seq = query_sequence.upper().strip()
        query_len = len(query_seq)
        
        if query_len < 20:
            logger.warning("Query sequence too short for reliable matching")
            return None
        
        best_match = None
        best_score = 0.0
        min_similarity = 0.8  # Minimum 80% similarity required
        
        # Search each reference chromosome
        for chr_name, ref_seq in reference_sequences.items():
            ref_seq = ref_seq.upper()
            ref_len = len(ref_seq)
            
            # Skip if reference is shorter than query
            if ref_len < query_len:
                continue
            
            # Use step size for efficiency (every 10 bp for sequences > 100bp)
            step_size = max(1, min(10, query_len // 10))
            
            # Sliding window search
            for start in range(0, ref_len - query_len + 1, step_size):
                ref_window = ref_seq[start:start + query_len]
                
                # Calculate similarity score
                similarity = self._calculate_sequence_similarity(query_seq, ref_window)
                
                if similarity > best_score:
                    best_score = similarity
                    best_match = (chr_name, start, start + query_len)
                    
                    # Early exit if we find a perfect or near-perfect match
                    if similarity >= 0.98:
                        break
            
            # Early exit if we found a very good match
            if best_score >= 0.98:
                break
        
        if best_match and best_score >= min_similarity:
            chr_name, start, end = best_match
            logger.debug(f"Best match: {chr_name}:{start}-{end} (similarity: {best_score:.3f})")
            logger.debug("=== END SEQUENCE MATCHING ===")
            return best_match
        else:
            logger.debug(f"No suitable match found (best score: {best_score:.3f})")
            logger.debug("=== END SEQUENCE MATCHING ===")
            return None
    
    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence similarity as fraction of matching nucleotides.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Similarity score between 0.0 and 1.0
        """
        if len(seq1) != len(seq2):
            return 0.0
        
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)
    
    def prepare_fragments_for_primer_design(self, sequences: Dict[str, str]) -> List[Dict[str, Any]]:
        """
        Prepare sequence fragments for primer design using restriction site cutting.
        
        Args:
            sequences: Dictionary of processed sequences
            
        Returns:
            List of restriction fragments ready for primer design
            
        Raises:
            SequenceProcessingError: If fragment preparation fails
        """
        logger.debug("=== FRAGMENT PREPARATION ===")
        logger.debug(f"Preparing {len(sequences)} sequences for primer design")
        
        try:
            # Import here to avoid circular imports
            from ..core import SequenceProcessor
            
            # Use restriction site cutting from standard pipeline
            restriction_fragments = SequenceProcessor.cut_at_restriction_sites(sequences)
            
            # Convert to simplified fragments (no gene overlap filtering)
            simplified_fragments = []
            for fragment in restriction_fragments:
                simplified_fragment = {
                    "id": fragment["id"],
                    "sequence": fragment["sequence"],
                    "chr": fragment.get("chr", ""),
                    "start": fragment.get("start", 1),
                    "end": fragment.get("end", len(fragment["sequence"])),
                    "Gene": fragment["id"].split("_")[0]  # Use the original sequence ID as gene name
                }
                simplified_fragments.append(simplified_fragment)
            
            logger.debug(f"Generated {len(simplified_fragments)} fragments after restriction site cutting")
            logger.debug("=== END FRAGMENT PREPARATION ===")
            return simplified_fragments
            
        except Exception as e:
            error_msg = f"Fragment preparation failed: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
    
    def run_workflow(self, table_file: PathLike, output_dir: PathLike) -> bool:
        """
        Execute the complete direct mode primer design workflow.
        
        Args:
            table_file: Path to CSV/Excel file with target sequences
            output_dir: Output directory for results
            
        Returns:
            True if workflow completed successfully, False otherwise
            
        Raises:
            SequenceProcessingError: If sequence processing fails
            PrimerDesignError: If primer design fails
        """
        logger.info("=== Direct Mode Workflow ===")
        
        try:
            # Step 1: Load sequences from table
            logger.info("Loading target sequences from table...")
            sequences = self.load_sequences_from_table(table_file)
            
            if not sequences:
                logger.warning("No valid sequences found in input table")
                return False
            
            # Step 2: Process with VCF if provided
            if self.vcf_file:
                logger.info("Processing sequences with VCF variants...")
                sequences = self.process_sequences_with_vcf(sequences)
            else:
                logger.info("Skipping VCF processing (no VCF file provided)")
            
            # Step 3: Prepare fragments using restriction site cutting
            logger.info("Preparing fragments with restriction site cutting...")
            fragments = self.prepare_fragments_for_primer_design(sequences)
            
            if not fragments:
                logger.warning("No valid fragments after restriction site cutting")
                return False
            
            # Step 4: Design primers using existing workflow components
            logger.info("Designing primers with Primer3...")
            primer_results = self._design_primers_with_primer3(fragments)
            
            if not primer_results:
                logger.warning("No primers were designed by Primer3")
                return False
            
            # Step 5: Filter primers
            logger.info("Filtering primers...")
            df = self._filter_primers(primer_results)
            
            if df is None or len(df) == 0:
                logger.warning("No primers passed filtering criteria")
                return False
            
            # Step 6: Calculate thermodynamic properties
            logger.info("Calculating thermodynamic properties...")
            df = self._calculate_thermodynamics(df)
            
            # Step 7: Run BLAST for specificity
            logger.info("Running BLAST for specificity checking...")
            df = self._run_blast_specificity(df)
            
            if df is None or len(df) == 0:
                logger.warning("No primers passed BLAST filtering")
                return False
            
            # Step 8: Save results
            output_path = FileIO.save_results(
                df, 
                output_dir, 
                table_file, 
                mode='direct'
            )
            
            if output_path:
                logger.info(f"Results saved to: {output_path}")
                return True
            else:
                logger.error("Failed to save results")
                return False
                
        except SequenceProcessingError as e:
            logger.error(f"Sequence processing error: {str(e)}")
            return False
        except PrimerDesignError as e:
            logger.error(f"Primer design error: {str(e)}")
            return False
        except Exception as e:
            logger.error(f"Error in direct mode workflow: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
    
    def _design_primers_with_primer3(self, fragments: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Design primers using Primer3 (reuses existing functionality).
        
        Args:
            fragments: List of sequence fragments
            
        Returns:
            List of primer design results
        """
        # Import here to avoid circular imports
        from ..main import design_primers_with_primer3
        
        return design_primers_with_primer3(fragments)
    
    def _filter_primers(self, primer_results: List[Dict[str, Any]]) -> Optional[pd.DataFrame]:
        """
        Filter primers using existing filtering logic.
        
        Args:
            primer_results: List of primer results from Primer3
            
        Returns:
            Filtered primer DataFrame, or None if no primers pass filtering
        """
        # Import here to avoid circular imports
        from ..main import filter_primers
        
        return filter_primers(primer_results)
    
    def _calculate_thermodynamics(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate thermodynamic properties using existing functionality.
        
        Args:
            df: DataFrame with primer information
            
        Returns:
            DataFrame with added thermodynamic properties
        """
        # Import here to avoid circular imports
        from ..main import calculate_thermodynamics
        
        return calculate_thermodynamics(df)
    
    def _run_blast_specificity(self, df: pd.DataFrame) -> Optional[pd.DataFrame]:
        """
        Run BLAST specificity checking using existing functionality.
        
        Args:
            df: DataFrame with primer information
            
        Returns:
            DataFrame with BLAST results, or None if no primers pass filtering
        """
        # Import here to avoid circular imports
        from ..main import run_blast_specificity
        
        return run_blast_specificity(df)