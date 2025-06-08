#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the sequence_analyzer module.

These tests verify that the SequenceAnalyzer class functions correctly,
including DNA sequence detection, file analysis, and column recommendation.
"""

import os
import pytest
import pandas as pd
import tempfile
from unittest.mock import patch, MagicMock

# Import package modules
from ...helpers import SequenceAnalyzer


# Test is_dna_sequence
def test_is_dna_sequence_valid():
    """Test identification of valid DNA sequences."""
    # Valid DNA sequences
    assert SequenceAnalyzer.is_dna_sequence("ATGCATGCATGCATGCATGC")
    assert SequenceAnalyzer.is_dna_sequence("ATGCNNNNNATGCATGCATGC")
    assert SequenceAnalyzer.is_dna_sequence("GTACGTACGTACGTACGTAC")
    
    # Valid DNA with ambiguous bases
    assert SequenceAnalyzer.is_dna_sequence("ATGCRYSWTGCATGCATGC")
    
    # Valid DNA with gaps
    assert SequenceAnalyzer.is_dna_sequence("ATGC--ATGCATGC--ATGC")


def test_is_dna_sequence_invalid():
    """Test rejection of invalid DNA sequences."""
    # Too short
    assert not SequenceAnalyzer.is_dna_sequence("ATGC")
    
    # Not a string
    assert not SequenceAnalyzer.is_dna_sequence(12345)
    assert not SequenceAnalyzer.is_dna_sequence(None)
    
    # Invalid characters
    assert not SequenceAnalyzer.is_dna_sequence("ATGCATGCATGCATGCATGCXXXX")
    assert not SequenceAnalyzer.is_dna_sequence("This is not a DNA sequence")
    
    # Below threshold of valid nucleotides
    assert not SequenceAnalyzer.is_dna_sequence("ATGC123456789ATGC123456789")


def test_is_dna_sequence_custom_length():
    """Test DNA sequence detection with custom minimum length."""
    # Short sequences with different minimum lengths
    assert not SequenceAnalyzer.is_dna_sequence("ATGCATGC", min_length=10)
    assert SequenceAnalyzer.is_dna_sequence("ATGCATGC", min_length=5)


# Test analyze_file
@patch('ddprimer.helpers.sequence_analyzer.pd.read_csv')
def test_analyze_csv_file(mock_read_csv):
    """Test analysis of a CSV file."""
    # Setup mock DataFrame
    mock_df = pd.DataFrame({
        'ID': ['seq1', 'seq2', 'seq3'],
        'Name': ['Sequence 1', 'Sequence 2', 'Sequence 3'],
        'Sequence': ['ATGCATGCATGCATGCATGC', 'GTACGTACGTACGTACGTAC', 'ATGCATGCATGCATGCATGC'],
        'Notes': ['note1', 'note2', 'note3']
    })
    mock_read_csv.return_value = mock_df
    
    # Call the function
    result = SequenceAnalyzer.analyze_file("test.csv")
    
    # Verify read_csv was called
    mock_read_csv.assert_called_once_with("test.csv")
    
    # Check the result
    assert isinstance(result, dict)
    assert "file_path" in result
    assert "total_columns" in result
    assert "total_rows" in result
    assert "likely_sequence_columns" in result
    assert "likely_name_columns" in result
    assert "column_stats" in result
    
    # Check column stats
    assert len(result["column_stats"]) == 4  # Four columns
    assert "ID" in result["column_stats"]
    assert "Sequence" in result["column_stats"]
    
    # Check likely sequence columns
    for col, stats in result["likely_sequence_columns"]:
        if col == "Sequence":
            assert stats["dna_percentage"] == 1.0
            assert stats["avg_length"] == 20.0


@patch('ddprimer.helpers.sequence_analyzer.pd.read_excel')
def test_analyze_excel_file(mock_read_excel):
    """Test analysis of an Excel file."""
    # Setup mock DataFrame
    mock_df = pd.DataFrame({
        'ID': ['seq1', 'seq2', 'seq3'],
        'Name': ['Sequence 1', 'Sequence 2', 'Sequence 3'],
        'Sequence': ['ATGCATGCATGCATGCATGC', 'GTACGTACGTACGTACGTAC', 'ATGCATGCATGCATGCATGC'],
        'Notes': ['note1', 'note2', 'note3']
    })
    mock_read_excel.return_value = mock_df
    
    # Call the function
    result = SequenceAnalyzer.analyze_file("test.xlsx")
    
    # Verify read_excel was called
    mock_read_excel.assert_called_once_with("test.xlsx")
    
    # Check the result (similar to CSV test)
    assert isinstance(result, dict)
    assert "likely_sequence_columns" in result
    assert "likely_name_columns" in result


def test_analyze_file_unsupported_format():
    """Test handling of unsupported file formats."""
    # Call the function with an unsupported format
    result = SequenceAnalyzer.analyze_file("test.txt")
    
    # Check the result
    assert "error" in result
    assert result["error"] == "Unsupported file format"


@patch('ddprimer.helpers.sequence_analyzer.pd.read_csv')
def test_analyze_file_with_error(mock_read_csv):
    """Test error handling during file analysis."""
    # Setup mock to raise an exception
    mock_read_csv.side_effect = Exception("File reading error")
    
    # Call the function
    result = SequenceAnalyzer.analyze_file("test.csv")
    
    # Check the result
    assert "error" in result
    assert "File reading error" in result["error"]


# Test print_analysis
@patch('ddprimer.helpers.sequence_analyzer.logging.getLogger')
def test_print_analysis_success(mock_get_logger):
    """Test printing analysis results."""
    # Setup mock logger
    mock_logger = MagicMock()
    mock_get_logger.return_value = mock_logger
    
    # Setup test analysis
    analysis = {
        "file_path": "test.csv",
        "total_columns": 4,
        "total_rows": 3,
        "likely_sequence_columns": [
            ("Sequence", {"avg_length": 20.0, "dna_percentage": 1.0})
        ],
        "likely_name_columns": [
            ("ID", {"unique_percentage": 1.0, "avg_length": 4.0})
        ],
        "column_stats": {
            "Sequence": {"avg_length": 20.0, "dna_percentage": 1.0},
            "ID": {"unique_percentage": 1.0, "avg_length": 4.0}
        }
    }
    
    # Call the function
    SequenceAnalyzer.print_analysis(analysis)
    
    # Verify logger was called multiple times
    assert mock_logger.debug.call_count > 0


@patch('ddprimer.helpers.sequence_analyzer.logging.getLogger')
def test_print_analysis_with_error(mock_get_logger):
    """Test printing analysis results with error."""
    # Setup mock logger
    mock_logger = MagicMock()
    mock_get_logger.return_value = mock_logger
    
    # Setup test analysis with error
    analysis = {"error": "Analysis failed"}
    
    # Call the function
    SequenceAnalyzer.print_analysis(analysis)
    
    # Verify error was logged
    mock_logger.error.assert_called_once()
    assert "Analysis failed" in mock_logger.error.call_args[0][0]


@patch('ddprimer.helpers.sequence_analyzer.logging.getLogger')
def test_print_analysis_no_sequences(mock_get_logger):
    """Test printing analysis results with no sequence columns."""
    # Setup mock logger
    mock_logger = MagicMock()
    mock_get_logger.return_value = mock_logger
    
    # Setup test analysis without sequence columns
    analysis = {
        "file_path": "test.csv",
        "total_columns": 4,
        "total_rows": 3,
        "likely_sequence_columns": [],
        "likely_name_columns": [
            ("ID", {"unique_percentage": 1.0, "avg_length": 4.0})
        ],
        "column_stats": {
            "ID": {"unique_percentage": 1.0, "avg_length": 4.0}
        }
    }
    
    # Call the function
    SequenceAnalyzer.print_analysis(analysis)
    
    # Verify warning was logged
    mock_logger.warning.assert_called_once()
    assert "No likely sequence columns" in mock_logger.warning.call_args[0][0]


# Test get_recommended_columns
@patch('ddprimer.helpers.sequence_analyzer.logging.getLogger')
def test_get_recommended_columns_success(mock_get_logger):
    """Test getting recommended columns from analysis."""
    # Setup mock logger
    mock_logger = MagicMock()
    mock_get_logger.return_value = mock_logger
    
    # Setup test analysis
    analysis = {
        "likely_sequence_columns": [
            ("Sequence", {"dna_percentage": 1.0}),
            ("AltSequence", {"dna_percentage": 0.9})
        ],
        "likely_name_columns": [
            ("ID", {"unique_percentage": 1.0}),
            ("Name", {"unique_percentage": 0.9})
        ]
    }
    
    # Call the function
    name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
    
    # Verify the results
    assert name_col == "ID"  # First name column
    assert seq_col == "Sequence"  # First sequence column


@patch('ddprimer.helpers.sequence_analyzer.logging.getLogger')
def test_get_recommended_columns_partial(mock_get_logger):
    """Test getting recommended columns with partial data."""
    # Setup mock logger
    mock_logger = MagicMock()
    mock_get_logger.return_value = mock_logger
    
    # Setup test analysis with only sequence columns
    analysis = {
        "likely_sequence_columns": [
            ("Sequence", {"dna_percentage": 1.0})
        ],
        "likely_name_columns": []
    }
    
    # Call the function
    name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
    
    # Verify the results
    assert name_col is None  # No name columns
    assert seq_col == "Sequence"
    
    # Setup test analysis with only name columns
    analysis = {
        "likely_sequence_columns": [],
        "likely_name_columns": [
            ("ID", {"unique_percentage": 1.0})
        ]
    }
    
    # Call the function
    name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
    
    # Verify the results
    assert name_col == "ID"
    assert seq_col is None  # No sequence columns


@patch('ddprimer.helpers.sequence_analyzer.logging.getLogger')
def test_get_recommended_columns_with_error(mock_get_logger):
    """Test getting recommended columns with analysis error."""
    # Setup mock logger
    mock_logger = MagicMock()
    mock_get_logger.return_value = mock_logger
    
    # Setup test analysis with error
    analysis = {"error": "Analysis failed"}
    
    # Call the function
    name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
    
    # Verify the results
    assert name_col is None
    assert seq_col is None
    
    # Verify error was logged
    mock_logger.error.assert_called_once()
    assert "Cannot get recommended columns" in mock_logger.error.call_args[0][0]


# Integration test with real-like data
def test_analyze_with_real_data():
    """Test analysis with realistic data using a temporary file."""
    # Create test DataFrame
    df = pd.DataFrame({
        'ID': ['seq1', 'seq2', 'seq3', 'seq4'],
        'Description': ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4'],
        'DNASequence': [
            'ATGCATGCATGCATGCATGCATGCATGC',
            'GTACGTACGTACGTACGTACGTACGTAC',
            'ATATATATATATATATATATATAT',
            'GCGCGCGCGCGCGCGCGCGCGCGC'
        ],
        'Protein': [
            'MHACKSILVPRTQED',  # Made longer with amino acids not in DNA alphabet
            'VTPGSAVLNMKRGED',
            'IYIPATATDNEGEKL',
            'ARARDEDEQWSTNVP'
        ],
        'Numeric': [1, 2, 3, 4]
    })
    
    # Save to temporary CSV file
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as temp_file:
        df.to_csv(temp_file.name, index=False)
        temp_file_path = temp_file.name
    
    try:
        # Analyze the file
        result = SequenceAnalyzer.analyze_file(temp_file_path)
        
        # Verify results
        assert "DNASequence" in [col for col, _ in result["likely_sequence_columns"]]
        assert "Protein" not in [col for col, _ in result["likely_sequence_columns"]]
        assert "ID" in [col for col, _ in result["likely_name_columns"]]
        
        # Check sorting of likely_name_columns - "ID" should be first due to highest uniqueness
        name_columns = [col for col, _ in result["likely_name_columns"]]
        assert name_columns.index("ID") == 0, f"ID should be the first name column, got {name_columns}"
        
        # Get recommended columns
        name_col, seq_col = SequenceAnalyzer.get_recommended_columns(result)
        assert name_col == "ID", f"Expected name_col to be 'ID', got '{name_col}'"
        assert seq_col == "DNASequence", f"Expected seq_col to be 'DNASequence', got '{seq_col}'"
        
    finally:
        # Clean up
        if os.path.exists(temp_file_path):
            os.unlink(temp_file_path)


if __name__ == "__main__":
    pytest.main(["-v"])