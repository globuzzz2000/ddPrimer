#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Direct mode utilities for ddPrimer pipeline.

Contains functionality for:
1. Sequence location finding using BLAST
2. CSV/Excel file analysis and column detection
3. Missing sequence handling and result completion

This module consolidates utilities specifically needed for direct mode
sequence processing and file analysis.
"""

import os
import logging
import pandas as pd
import tempfile
import subprocess

# Import package modules
from ..config import SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)


class DirectModeUtils:
    """
    Utilities for direct mode sequence processing.
    
    This class consolidates functions needed for direct mode operation,
    including sequence location finding and file analysis capabilities.
    
    Example:
        >>> chrom, start, end, identity = DirectModeUtils.find_sequence_location(seq, ref_fasta)
        >>> analysis = DirectModeUtils.analyze_sequence_file("sequences.csv")
    """
    
    @staticmethod
    def find_sequence_location(sequence, ref_fasta, min_identity=90, min_coverage=90):
        """
        Find the location of a sequence in the reference genome using BLAST.
        
        Args:
            sequence (str): Query sequence to locate
            ref_fasta (str): Path to reference FASTA file
            min_identity (float): Minimum percent identity (default: 90)
            min_coverage (float): Minimum query coverage (default: 90)
            
        Returns:
            tuple: (chromosome, start_position, end_position, percent_identity) 
                   or (None, None, None, None) if not found
            
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
                # Split by lines first, then process each line
                lines = blast_output.strip().split('\n')
                for line in lines:
                    if not line.strip():
                        continue
                        
                    fields = line.strip().split('\t')
                    if len(fields) >= 9:
                        try:
                            # Extract fields by index to handle variable number of fields
                            qseqid = fields[0]
                            chrom = fields[1]
                            percent_identity = float(fields[2])
                            query_coverage = float(fields[3])
                            qstart = int(fields[4])
                            qend = int(fields[5])
                            start = int(fields[6])
                            end = int(fields[7])
                            length = int(fields[8])
                            
                            # Ensure start is always less than end (BLAST might report them reversed)
                            if start > end:
                                start, end = end, start
                            
                            # Check if the alignment meets our criteria
                            if percent_identity >= min_identity and query_coverage >= min_coverage:
                                logger.debug(f"Found sequence match on {chrom} at positions {start}-{end}, "
                                            f"identity: {percent_identity}%, coverage: {query_coverage}%")
                                return chrom, start, end, percent_identity
                                
                        except ValueError as parse_error:
                            logger.debug(f"Error parsing BLAST output line '{line}': {str(parse_error)}")
                            continue  # Try next line if this one is malformed
            
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
    def add_missing_sequences(df, all_sequences, matching_status=None):
        """
        Check for sequences that didn't get primers and add them to the results.
        Also add rows for sequences that failed to match to reference genome.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer results
            all_sequences (dict): Dictionary of all original sequences
            matching_status (dict, optional): Dictionary with reference matching status
                
        Returns:
            pandas.DataFrame: Updated DataFrame with rows for sequences without primers
        """
        logger.debug("Checking for sequences without primers and reference match failures...")
        
        # If there's no Gene column, we can't properly identify which sequences have primers
        if "Gene" not in df.columns:
            logger.debug("No 'Gene' column found in DataFrame, skipping missing sequences check")
            return df
        
        # 1) Which sequences actually got primers?
        sequences_with_primers = set()
        valid = df["Primer F"] != "No suitable primers found"
        sequences_with_primers.update(df.loc[valid, "Gene"].astype(str).unique())
        
        # 2) All the input sequence IDs
        all_input = set(all_sequences.keys())
        
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

    @staticmethod
    def analyze_sequence_file(file_path):
        """
        Analyze a CSV or Excel file to identify sequence columns and structure.
        
        Args:
            file_path (str): Path to the file to analyze
            
        Returns:
            dict: Analysis results including column information and recommendations
            
        Raises:
            SequenceProcessingError: If file analysis fails
        """
        logger.debug(f"Analyzing sequence file: {file_path}")
        
        try:
            # Load file based on extension
            if file_path.lower().endswith('.csv'):
                df = pd.read_csv(file_path)
            elif file_path.lower().endswith(('.xlsx', '.xls')):
                df = pd.read_excel(file_path)
            else:
                raise SequenceProcessingError(f"Unsupported file format: {file_path}")
            
            # Remove empty columns and rows
            df = df.dropna(axis=1, how='all')
            df = df.dropna(axis=0, how='all')
            
            if df.empty:
                return {"error": "File contains no data after removing empty rows/columns"}
            
            analysis = {
                "file_path": file_path,
                "shape": df.shape,
                "columns": list(df.columns),
                "column_analysis": {},
                "recommendations": {}
            }
            
            # Analyze each column
            for col in df.columns:
                col_data = df[col].dropna()
                if col_data.empty:
                    continue
                
                # Basic column statistics
                col_analysis = {
                    "non_null_count": len(col_data),
                    "sample_values": col_data.head(3).tolist(),
                    "likely_sequence": False,
                    "likely_name": False,
                    "avg_length": 0,
                    "has_dna_chars": False
                }
                
                # Check if this looks like a sequence column
                if col_data.dtype == 'object':  # String data
                    # Calculate average length
                    lengths = [len(str(val)) for val in col_data.head(10)]
                    col_analysis["avg_length"] = sum(lengths) / len(lengths) if lengths else 0
                    
                    # Check for DNA characters
                    sample_text = ' '.join(str(val) for val in col_data.head(5))
                    dna_chars = set('ATCGN')
                    text_chars = set(sample_text.upper())
                    
                    # If mostly DNA characters and reasonable length, likely a sequence
                    if (len(text_chars & dna_chars) >= 3 and 
                        col_analysis["avg_length"] > 15 and
                        len(text_chars - dna_chars - {' ', '\t', '\n'}) < 5):
                        col_analysis["likely_sequence"] = True
                        col_analysis["has_dna_chars"] = True
                    
                    # Check if this looks like names/IDs
                    elif col_analysis["avg_length"] < 50:
                        col_analysis["likely_name"] = True
                
                analysis["column_analysis"][col] = col_analysis
            
            # Generate recommendations
            sequence_candidates = [col for col, data in analysis["column_analysis"].items() 
                                 if data["likely_sequence"]]
            name_candidates = [col for col, data in analysis["column_analysis"].items() 
                             if data["likely_name"]]
            
            # Recommend best sequence column
            if sequence_candidates:
                # Pick the one with longest average length (most likely to be sequences)
                best_seq = max(sequence_candidates, 
                             key=lambda x: analysis["column_analysis"][x]["avg_length"])
                analysis["recommendations"]["sequence_column"] = best_seq
            
            # Recommend best name column
            if name_candidates:
                # Prefer first column if it looks like names, otherwise pick any name candidate
                if name_candidates[0] == df.columns[0]:
                    analysis["recommendations"]["name_column"] = name_candidates[0]
                else:
                    analysis["recommendations"]["name_column"] = name_candidates[0]
            
            # Add summary
            analysis["summary"] = {
                "total_rows": len(df),
                "total_columns": len(df.columns),
                "sequence_columns_found": len(sequence_candidates),
                "name_columns_found": len(name_candidates),
                "analysis_success": len(sequence_candidates) > 0
            }
            
            logger.debug(f"File analysis completed: {analysis['summary']}")
            return analysis
            
        except Exception as e:
            error_msg = f"Error analyzing file {file_path}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return {"error": error_msg}

    @staticmethod
    def get_recommended_columns(analysis):
        """
        Get recommended name and sequence columns from file analysis.
        
        Args:
            analysis (dict): File analysis results from analyze_sequence_file()
            
        Returns:
            tuple: (name_column, sequence_column) or (None, None) if not found
        """
        if "error" in analysis:
            logger.warning(f"Analysis contains error: {analysis['error']}")
            return None, None
        
        recommendations = analysis.get("recommendations", {})
        name_col = recommendations.get("name_column")
        seq_col = recommendations.get("sequence_column")
        
        logger.debug(f"Recommended columns: name='{name_col}', sequence='{seq_col}'")
        return name_col, seq_col

    @staticmethod
    def print_analysis(analysis):
        """
        Print a user-friendly analysis report.
        
        Args:
            analysis (dict): File analysis results from analyze_sequence_file()
        """
        if "error" in analysis:
            logger.error(f"File analysis failed: {analysis['error']}")
            return
        
        summary = analysis.get("summary", {})
        recommendations = analysis.get("recommendations", {})
        
        logger.info("=== File Analysis Results ===")
        logger.info(f"File: {os.path.basename(analysis['file_path'])}")
        logger.info(f"Dimensions: {summary['total_rows']} rows × {summary['total_columns']} columns")
        
        if recommendations.get("sequence_column"):
            logger.info(f"Sequence column detected: '{recommendations['sequence_column']}'")
        else:
            logger.warning("No sequence column detected")
        
        if recommendations.get("name_column"):
            logger.info(f"Name column detected: '{recommendations['name_column']}'")
        else:
            logger.info("No name column detected - will generate sequential IDs")
        
        # Show column details in debug mode
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Column analysis:")
            for col, data in analysis.get("column_analysis", {}).items():
                logger.debug(f"  {col}: {data['non_null_count']} values, "
                           f"avg_len={data['avg_length']:.1f}, "
                           f"sequence={data['likely_sequence']}, "
                           f"name={data['likely_name']}")


# Backward compatibility aliases for easy migration
class DirectMasking:
    """Backward compatibility wrapper for DirectModeUtils."""
    
    @staticmethod
    def find_location(*args, **kwargs):
        """Alias for find_sequence_location."""
        return DirectModeUtils.find_sequence_location(*args, **kwargs)
    
    @staticmethod
    def add_missing_sequences(*args, **kwargs):
        """Alias for add_missing_sequences."""
        return DirectModeUtils.add_missing_sequences(*args, **kwargs)


class SequenceAnalyzer:
    """Backward compatibility wrapper for DirectModeUtils."""
    
    @staticmethod
    def analyze_file(*args, **kwargs):
        """Alias for analyze_sequence_file."""
        return DirectModeUtils.analyze_sequence_file(*args, **kwargs)
    
    @staticmethod
    def get_recommended_columns(*args, **kwargs):
        """Alias for get_recommended_columns."""
        return DirectModeUtils.get_recommended_columns(*args, **kwargs)
    
    @staticmethod
    def print_analysis(*args, **kwargs):
        """Alias for print_analysis."""
        return DirectModeUtils.print_analysis(*args, **kwargs)