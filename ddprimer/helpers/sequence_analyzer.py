#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence file analysis module for ddPrimer pipeline.

This module provides utilities for analyzing sequence files and detecting:
1. Likely DNA sequence columns in CSV/Excel files
2. Likely sequence name/ID columns
3. Statistics about column content for improved file parsing
"""

import logging
import pandas as pd
from ..config import Config


class SequenceAnalyzer:
    """Helper class to analyze sequence files and detect DNA sequence columns."""
    
    @staticmethod
    def is_dna_sequence(text, min_length=15):
        """
        Check if a string is likely a DNA sequence.
        
        Args:
            text (str): The string to check
            min_length (int): Minimum length to qualify as a sequence
            
        Returns:
            bool: True if the text appears to be a DNA sequence, False otherwise
        """
        logger = logging.getLogger("ddPrimer")
        
        if not isinstance(text, str):
            logger.debug("Invalid sequence type provided to is_dna_sequence")
            return False
            
        # Minimum length check
        if len(text) < min_length:
            return False
            
        # Check if it has a high percentage of valid nucleotides
        valid_chars = set('ATGCNRYSWKMBDHV.-')
        text_upper = text.upper()
        valid_count = sum(1 for c in text_upper if c in valid_chars)
        valid_percentage = valid_count / len(text) if len(text) > 0 else 0
        
        return valid_percentage >= 0.9  # At least 90% valid nucleotides
    
    @staticmethod
    def analyze_file(file_path):
        """
        Analyze a CSV or Excel file to find potential sequence columns.
        
        Args:
            file_path (str): Path to the file to analyze
            
        Returns:
            dict: Analysis results including potential sequence and name columns
        """
        logger = logging.getLogger("ddPrimer")
        
        try:
            # Load the file with pandas
            if file_path.endswith('.csv'):
                df = pd.read_csv(file_path)
            elif file_path.endswith(('.xlsx', '.xls')):
                df = pd.read_excel(file_path)
            else:
                logger.error(f"Unsupported file format: {file_path}")
                return {"error": "Unsupported file format"}
            
            # Remove completely empty columns and rows
            df = df.dropna(axis=1, how='all')
            df = df.dropna(axis=0, how='all')
            
            # For each column, check:
            # 1. Average string length
            # 2. Percentage of values that look like DNA sequences
            # 3. Number of unique values
            column_stats = {}
            
            for col in df.columns:
                values = df[col].dropna()
                if len(values) == 0:
                    continue
                
                # Convert all to strings
                values = values.astype(str)
                
                # Calculate stats
                avg_length = sum(len(str(v)) for v in values) / len(values)
                unique_count = len(values.unique())
                unique_percentage = unique_count / len(values)
                
                # Check how many look like DNA sequences
                dna_count = sum(1 for v in values if SequenceAnalyzer.is_dna_sequence(str(v)))
                dna_percentage = dna_count / len(values) if len(values) > 0 else 0
                
                column_stats[col] = {
                    "avg_length": avg_length,
                    "unique_count": unique_count,
                    "unique_percentage": unique_percentage,
                    "dna_count": dna_count,
                    "dna_percentage": dna_percentage,
                    "total_values": len(values)
                }
            
            # Identify likely sequence columns (high DNA percentage, longer strings)
            likely_seq_cols = []
            for col, stats in column_stats.items():
                if stats["dna_percentage"] >= 0.8 and stats["avg_length"] >= 20:
                    likely_seq_cols.append((col, stats))
            
            # Sort by DNA percentage then average length
            likely_seq_cols.sort(key=lambda x: (x[1]["dna_percentage"], x[1]["avg_length"]), reverse=True)
            
            # Identify likely name columns (high uniqueness, shorter strings)
            likely_name_cols = []
            for col, stats in column_stats.items():
                if stats["unique_percentage"] >= 0.8 and stats["avg_length"] < 50 and col not in [seq_col for seq_col, _ in likely_seq_cols]:
                    likely_name_cols.append((col, stats))
            
            # Sort by uniqueness then alphabetically
            likely_name_cols.sort(key=lambda x: (x[1]["unique_percentage"], str(x[0])), reverse=True)
            
            return {
                "file_path": file_path,
                "total_columns": len(df.columns),
                "total_rows": len(df),
                "likely_sequence_columns": likely_seq_cols,
                "likely_name_columns": likely_name_cols,
                "column_stats": column_stats
            }
            
        except Exception as e:
            logger.error(f"Error analyzing file: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return {"error": str(e)}
    
    @staticmethod
    def print_analysis(analysis, logger=None):
        """
        Print analysis results in a readable format.
        
        Args:
            analysis (dict): Analysis results from analyze_file
            logger (logging.Logger, optional): Logger to use for output. If None, uses default ddPrimer logger.
        """
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
        if "error" in analysis:
            logger.error(f"Analysis error: {analysis['error']}")
            return
        
        logger.info(f"\nAnalysis of file: {analysis['file_path']}")
        logger.info(f"Total columns: {analysis['total_columns']}, Total rows: {analysis['total_rows']}")
        
        if analysis['likely_sequence_columns']:
            logger.info("\nLikely sequence columns:")
            for i, (col, stats) in enumerate(analysis['likely_sequence_columns'], 1):
                logger.info(f"  {i}. '{col}': {stats['avg_length']:.1f} chars avg, {stats['dna_percentage']*100:.1f}% DNA-like")
        else:
            logger.info("\nNo likely sequence columns detected.")
        
        if analysis['likely_name_columns']:
            logger.info("\nLikely name/ID columns:")
            for i, (col, stats) in enumerate(analysis['likely_name_columns'], 1):
                logger.info(f"  {i}. '{col}': {stats['unique_percentage']*100:.1f}% unique values, {stats['avg_length']:.1f} chars avg")
        else:
            logger.info("\nNo likely name/ID columns detected.")
    
    @staticmethod
    def get_recommended_columns(analysis):
        """
        Get recommended columns for sequence name and data based on analysis results.
        
        Args:
            analysis (dict): Analysis results from analyze_file
            
        Returns:
            tuple: (recommended_name_column, recommended_sequence_column)
        """
        logger = logging.getLogger("ddPrimer")
        
        if "error" in analysis:
            logger.error(f"Cannot get recommended columns: {analysis['error']}")
            return None, None
            
        # Default to None for both
        recommended_name_col = None
        recommended_seq_col = None
        
        # Get the best sequence column if available
        if analysis['likely_sequence_columns']:
            recommended_seq_col = analysis['likely_sequence_columns'][0][0]
            logger.debug(f"Recommended sequence column: '{recommended_seq_col}'")
        
        # Get the best name column if available
        if analysis['likely_name_columns']:
            recommended_name_col = analysis['likely_name_columns'][0][0]
            logger.debug(f"Recommended name column: '{recommended_name_col}'")
            
        return recommended_name_col, recommended_seq_col