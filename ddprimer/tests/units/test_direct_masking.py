#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the direct_masking module.

These tests verify that the DirectMasking class functions correctly,
including finding sequence locations and handling sequences without primers.
"""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

# Import package modules
from ...helpers import DirectMasking
from ...config import SequenceProcessingError


# Test DirectMasking.find_location
@patch('ddprimer.helpers.direct_masking.tempfile.NamedTemporaryFile')
@patch('ddprimer.helpers.direct_masking.subprocess.run')
def test_find_location_success(mock_subprocess_run, mock_temp_file):
    """Test successful sequence location finding with BLAST."""
    # Setup mock temporary file
    mock_temp_file_instance = MagicMock()
    mock_temp_file_instance.name = "/tmp/test_seq.fasta"
    mock_temp_file.return_value.__enter__.return_value = mock_temp_file_instance
    
    # Setup mock subprocess result
    mock_blast_result = MagicMock()
    mock_blast_result.stdout = "query\tchr1\t95.5\t98.2\t1\t100\t500\t600\t101"
    mock_blast_result.returncode = 0
    mock_subprocess_run.return_value = mock_blast_result
    
    # Call the function
    chrom, start, end, identity = DirectMasking.find_location(
        "ATGCATGCATGC", "test_genome.fasta", min_identity=90, min_coverage=90
    )
    
    # Verify the function wrote the query sequence to a temporary file
    mock_temp_file_instance.write.assert_any_call(">query\n")
    mock_temp_file_instance.write.assert_any_call("ATGCATGCATGC")
    
    # Verify BLAST was called with correct parameters
    mock_subprocess_run.assert_called_once()
    blastn_cmd = mock_subprocess_run.call_args[0][0]
    assert "blastn" == blastn_cmd[0]
    assert "-query" == blastn_cmd[1]
    assert "/tmp/test_seq.fasta" == blastn_cmd[2]
    assert "-subject" == blastn_cmd[3]
    assert "test_genome.fasta" == blastn_cmd[4]
    
    # Check the results
    assert chrom == "chr1"
    assert start == 500
    assert end == 600
    assert identity == 95.5


@patch('ddprimer.helpers.direct_masking.tempfile.NamedTemporaryFile')
@patch('ddprimer.helpers.direct_masking.subprocess.run')
def test_find_location_reversed_coordinates(mock_subprocess_run, mock_temp_file):
    """Test handling reversed coordinates in BLAST results."""
    # Setup mock temporary file
    mock_temp_file_instance = MagicMock()
    mock_temp_file_instance.name = "/tmp/test_seq.fasta"
    mock_temp_file.return_value.__enter__.return_value = mock_temp_file_instance
    
    # Setup mock subprocess result with reversed coordinates (end < start)
    mock_blast_result = MagicMock()
    mock_blast_result.stdout = "query\tchr1\t95.5\t98.2\t1\t100\t600\t500\t101"
    mock_blast_result.returncode = 0
    mock_subprocess_run.return_value = mock_blast_result
    
    # Call the function
    chrom, start, end, identity = DirectMasking.find_location(
        "ATGCATGCATGC", "test_genome.fasta", min_identity=90, min_coverage=90
    )
    
    # Check the results - coordinates should be swapped
    assert chrom == "chr1"
    assert start == 500  # The smaller value
    assert end == 600    # The larger value
    assert identity == 95.5


@patch('ddprimer.helpers.direct_masking.tempfile.NamedTemporaryFile')
@patch('ddprimer.helpers.direct_masking.subprocess.run')
def test_find_location_no_match(mock_subprocess_run, mock_temp_file):
    """Test behavior when no sequence match is found."""
    # Setup mock temporary file
    mock_temp_file_instance = MagicMock()
    mock_temp_file_instance.name = "/tmp/test_seq.fasta"
    mock_temp_file.return_value.__enter__.return_value = mock_temp_file_instance
    
    # Setup mock subprocess result with empty output
    mock_blast_result = MagicMock()
    mock_blast_result.stdout = ""
    mock_blast_result.returncode = 0
    mock_subprocess_run.return_value = mock_blast_result
    
    # Call the function
    chrom, start, end, identity = DirectMasking.find_location(
        "ATGCATGCATGC", "test_genome.fasta", min_identity=90, min_coverage=90
    )
    
    # Check the results - should all be None
    assert chrom is None
    assert start is None
    assert end is None
    assert identity is None


@patch('ddprimer.helpers.direct_masking.tempfile.NamedTemporaryFile')
@patch('ddprimer.helpers.direct_masking.subprocess.run')
def test_find_location_low_quality_match(mock_subprocess_run, mock_temp_file):
    """Test behavior when match quality is below thresholds."""
    # Setup mock temporary file
    mock_temp_file_instance = MagicMock()
    mock_temp_file_instance.name = "/tmp/test_seq.fasta"
    mock_temp_file.return_value.__enter__.return_value = mock_temp_file_instance
    
    # Setup mock subprocess result with low quality match
    mock_blast_result = MagicMock()
    mock_blast_result.stdout = "query\tchr1\t85.5\t80.2\t1\t100\t500\t600\t101"
    mock_blast_result.returncode = 0
    mock_subprocess_run.return_value = mock_blast_result
    
    # Call the function with higher thresholds
    chrom, start, end, identity = DirectMasking.find_location(
        "ATGCATGCATGC", "test_genome.fasta", min_identity=90, min_coverage=90
    )
    
    # Check the results - should all be None
    assert chrom is None
    assert start is None
    assert end is None
    assert identity is None


@patch('ddprimer.helpers.direct_masking.tempfile.NamedTemporaryFile')
@patch('ddprimer.helpers.direct_masking.subprocess.run')
def test_find_location_blast_error(mock_subprocess_run, mock_temp_file):
    """Test handling of BLAST command errors."""
    # Setup mock temporary file
    mock_temp_file_instance = MagicMock()
    mock_temp_file_instance.name = "/tmp/test_seq.fasta"
    mock_temp_file.return_value.__enter__.return_value = mock_temp_file_instance
    
    # Setup mock subprocess to raise an error
    mock_subprocess_run.side_effect = Exception("BLAST command failed")
    
    # Call the function and verify it raises the expected exception
    with pytest.raises(SequenceProcessingError, match="BLAST search failed"):
        DirectMasking.find_location("ATGCATGCATGC", "test_genome.fasta")


@patch('ddprimer.helpers.direct_masking.tempfile.NamedTemporaryFile')
@patch('ddprimer.helpers.direct_masking.subprocess.run')
@patch('ddprimer.helpers.direct_masking.os.unlink')
def test_find_location_temp_file_cleanup(mock_unlink, mock_subprocess_run, mock_temp_file):
    """Test cleanup of temporary files even when errors occur."""
    # Setup mock temporary file
    mock_temp_file_instance = MagicMock()
    mock_temp_file_instance.name = "/tmp/test_seq.fasta"
    mock_temp_file.return_value.__enter__.return_value = mock_temp_file_instance
    
    # Setup mock for os.path.exists to return True
    with patch('ddprimer.helpers.direct_masking.os.path.exists', return_value=True):
        # Setup mock subprocess to raise an error
        mock_subprocess_run.side_effect = Exception("BLAST command failed")
        
        # Call the function and catch the exception
        with pytest.raises(SequenceProcessingError):
            DirectMasking.find_location("ATGCATGCATGC", "test_genome.fasta")
    
    # Verify the temporary file was cleaned up
    mock_unlink.assert_called_once_with("/tmp/test_seq.fasta")


# Test DirectMasking.add_missing_sequences
def test_add_missing_sequences_no_missing():
    """Test when all sequences got primers assigned."""
    # Create test data with all sequences having primers
    df = pd.DataFrame({
        "Gene": ["seq1", "seq2", "seq3"],
        "Primer F": ["ATGC", "CGTA", "GCTA"]
    })
    
    masked_sequences = {
        "seq1": "ATGCATGC",
        "seq2": "CGTATAGC",
        "seq3": "GCTATAGC"
    }
    
    # Call the function
    result_df = DirectMasking.add_missing_sequences(df, masked_sequences)
    
    # Verify the dataframe is unchanged
    assert len(result_df) == len(df)
    assert all(result_df["Gene"].values == df["Gene"].values)
    assert all(result_df["Primer F"].values == df["Primer F"].values)


def test_add_missing_sequences_with_no_primers():
    """Test when some sequences didn't get primers assigned."""
    # Create test data with some sequences having primers
    df = pd.DataFrame({
        "Gene": ["seq1", "seq2"],
        "Primer F": ["ATGC", "CGTA"]
    })
    
    masked_sequences = {
        "seq1": "ATGCATGC",
        "seq2": "CGTATAGC",
        "seq3": "GCTATAGC",  # This sequence doesn't have primers
        "seq4": "TATAGCTC"   # This sequence doesn't have primers
    }
    
    # Call the function
    result_df = DirectMasking.add_missing_sequences(df, masked_sequences)
    
    # Verify the result
    assert len(result_df) == 4  # Original 2 + 2 missing sequences
    
    # Check the added rows
    missing_rows = result_df[result_df["Gene"].isin(["seq3", "seq4"])]
    assert len(missing_rows) == 2
    assert all(missing_rows["Primer F"] == "No suitable primers found")
    assert all(missing_rows["Reference Match"] == "Not attempted")


def test_add_missing_sequences_with_match_failures():
    """Test when some sequences failed to match to the reference."""
    # Create test data
    df = pd.DataFrame({
        "Gene": ["seq1", "seq2"],
        "Primer F": ["ATGC", "CGTA"]
    })
    
    masked_sequences = {
        "seq1": "ATGCATGC",
        "seq2": "CGTATAGC",
        "seq3": "GCTATAGC",   # No primers but successful match
        "seq4": "TATAGCTC"    # Match failure
    }
    
    # Match status dictionary
    matching_status = {
        "seq1": "Success",
        "seq2": "Success", 
        "seq3": "Success",
        "seq4": "Failure"     # This sequence failed to match
    }
    
    # Call the function
    result_df = DirectMasking.add_missing_sequences(df, masked_sequences, matching_status)
    
    # Verify the result
    assert len(result_df) == 4  # Original 2 + 1 no primer + 1 match failure
    
    # Check the successful sequences with the Reference Match column
    assert all(result_df.loc[result_df["Gene"].isin(["seq1", "seq2"]), "Reference Match"] == "Success")
    
    # Check the no-primer row
    no_primer_row = result_df[result_df["Gene"] == "seq3"]
    assert len(no_primer_row) == 1
    assert no_primer_row["Primer F"].values[0] == "No suitable primers found"
    assert no_primer_row["Reference Match"].values[0] == "Success"
    
    # Check the match failure row
    failure_row = result_df[result_df["Gene"] == "seq4"]
    assert len(failure_row) == 1
    assert failure_row["Primer F"].values[0] == "Sequence could not be matched against reference"
    assert failure_row["Reference Match"].values[0] == "Failure"


def test_add_missing_sequences_no_gene_column():
    """Test behavior when the dataframe has no Gene column."""
    # Create test data with no Gene column
    df = pd.DataFrame({
        "ID": ["seq1", "seq2"],
        "Primer F": ["ATGC", "CGTA"]
    })
    
    masked_sequences = {
        "seq1": "ATGCATGC",
        "seq2": "CGTATAGC",
        "seq3": "GCTATAGC"  # This sequence doesn't have primers
    }
    
    # Call the function
    result_df = DirectMasking.add_missing_sequences(df, masked_sequences)
    
    # Verify the dataframe is unchanged since there's no Gene column to match against
    assert len(result_df) == len(df)
    assert "Gene" not in result_df.columns


if __name__ == "__main__":
    pytest.main(["-v"])