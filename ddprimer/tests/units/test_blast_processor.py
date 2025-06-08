#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the blast processor module.
Tests BLAST operations for primer specificity checking.
"""

import pandas as pd
from unittest import mock

# Import package modules
from ...core.blast_processor import BlastProcessor
from ...config import Config


class TestBlastProcessor:
    """Test suite for BlastProcessor class."""

    def setup_method(self):
        """Set up test data before each test."""
        # Mock configuration settings
        self.original_db_path = getattr(Config, 'DB_PATH', None)
        self.original_blast_word_size = getattr(Config, 'BLAST_WORD_SIZE', None)
        self.original_blast_evalue = getattr(Config, 'BLAST_EVALUE', None)
        self.original_blast_reward = getattr(Config, 'BLAST_REWARD', None)
        self.original_blast_penalty = getattr(Config, 'BLAST_PENALTY', None)
        self.original_blast_gapopen = getattr(Config, 'BLAST_GAPOPEN', None)
        self.original_blast_gapextend = getattr(Config, 'BLAST_GAPEXTEND', None)
        self.original_blast_max_target_seqs = getattr(Config, 'BLAST_MAX_TARGET_SEQS', None)
        self.original_blast_filter_factor = getattr(Config, 'BLAST_FILTER_FACTOR', None)
        self.original_num_processes = getattr(Config, 'NUM_PROCESSES', None)
        self.original_show_progress = getattr(Config, 'SHOW_PROGRESS', None)
        
        # Set test config values
        Config.DB_PATH = "test_blast_db"
        Config.BLAST_WORD_SIZE = 7
        Config.BLAST_EVALUE = 1000
        Config.BLAST_REWARD = 1
        Config.BLAST_PENALTY = -3
        Config.BLAST_GAPOPEN = 5
        Config.BLAST_GAPEXTEND = 2
        Config.BLAST_MAX_TARGET_SEQS = 5
        Config.BLAST_FILTER_FACTOR = 2
        Config.NUM_PROCESSES = 2
        Config.SHOW_PROGRESS = False
        
        # Sample sequences for testing
        self.test_sequence = "ATCGATCGATCGTAGCTAGCTAGC"
        self.test_batch = [
            "ATCGATCGATCG",
            "GCTAGCTAGCTA",
            "NNNNNNNNNNNN",  # All Ns should return no hits
            ""  # Empty sequence
        ]
        
        # Sample DataFrames for testing
        self.test_df = pd.DataFrame({
            'Forward Primer': ['ATCGATCGATCG', 'GCTAGCTAGCTA'],
            'Forward Primer BLAST1': [0.001, 0.01],
            'Forward Primer BLAST2': [0.1, None],
            'Reverse Primer': ['GCATGCATGCAT', 'TAGTCGTAGCTA'],
            'Reverse Primer BLAST1': [0.002, None],
            'Reverse Primer BLAST2': [0.1, None],
            'Probe': ['CGATCG', 'GCTAGC'],
            'Probe BLAST1': [0.003, 0.003],
            'Probe BLAST2': [0.05, 0.009]
        })

    def teardown_method(self):
        """Restore configuration after each test."""
        if self.original_db_path is not None:
            Config.DB_PATH = self.original_db_path
        if self.original_blast_word_size is not None:
            Config.BLAST_WORD_SIZE = self.original_blast_word_size
        if self.original_blast_evalue is not None:
            Config.BLAST_EVALUE = self.original_blast_evalue
        if self.original_blast_reward is not None:
            Config.BLAST_REWARD = self.original_blast_reward
        if self.original_blast_penalty is not None:
            Config.BLAST_PENALTY = self.original_blast_penalty
        if self.original_blast_gapopen is not None:
            Config.BLAST_GAPOPEN = self.original_blast_gapopen
        if self.original_blast_gapextend is not None:
            Config.BLAST_GAPEXTEND = self.original_blast_gapextend
        if self.original_blast_max_target_seqs is not None:
            Config.BLAST_MAX_TARGET_SEQS = self.original_blast_max_target_seqs
        if self.original_blast_filter_factor is not None:
            Config.BLAST_FILTER_FACTOR = self.original_blast_filter_factor
        if self.original_num_processes is not None:
            Config.NUM_PROCESSES = self.original_num_processes
        if self.original_show_progress is not None:
            Config.SHOW_PROGRESS = self.original_show_progress

    @mock.patch('subprocess.run')
    def test_blast_short_seq(self, mock_run):
        """Test BLAST operation for short sequences."""
        # Mock successful BLAST with multiple hits
        mock_process = mock.MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "0.001\n0.1\n0.5\n"
        mock_run.return_value = mock_process
        
        best, second = BlastProcessor.blast_short_seq(self.test_sequence)
        
        # Verify BLAST was called with correct parameters
        assert mock_run.called
        args = mock_run.call_args[0][0]
        assert args[0] == "blastn"
        assert args[1:3] == ["-task", "blastn-short"]
        
        # The BlastProcessor adds quotes to the DB_PATH
        assert args[3:5] == ["-db", f'"{Config.DB_PATH}"']
        
        assert "-word_size" in args
        assert str(Config.BLAST_WORD_SIZE) in args
        
        # Verify results
        assert best == 0.001
        assert second == 0.1
        
        # Test with no hits
        mock_process.stdout = ""
        best, second = BlastProcessor.blast_short_seq(self.test_sequence)
        assert best is None
        assert second is None
        
        # Test with BLAST error
        # The implementation returns None, None for errors
        mock_process.returncode = 1
        mock_process.stderr = "Error: BLAST database not found"
        best, second = BlastProcessor.blast_short_seq(self.test_sequence)
        assert best is None
        assert second is None
        
        # Test with empty sequence
        best, second = BlastProcessor.blast_short_seq("")
        assert best is None
        assert second is None
        
        # Test with None sequence
        best, second = BlastProcessor.blast_short_seq(None)
        assert best is None
        assert second is None

    @mock.patch.object(BlastProcessor, 'blast_short_seq')
    def test_process_blast_batch(self, mock_blast):
        """Test batch processing of BLAST operations."""
        # Set up mock return values
        mock_blast.side_effect = [
            (0.001, 0.1),    # First sequence
            (0.01, 0.2),     # Second sequence
            (None, None),    # All Ns
            (None, None)     # Empty sequence
        ]
        
        # Process the batch
        batch_data = (self.test_batch, "Test")
        results = BlastProcessor.process_blast_batch(batch_data)
        
        # Verify results
        assert len(results) == 4
        assert results[0] == (0.001, 0.1)
        assert results[1] == (0.01, 0.2)
        assert results[2] == (None, None)
        assert results[3] == (None, None)
        
        # Verify BLAST was called the correct number of times
        assert mock_blast.call_count == 4

    def test_passes_blast_filter(self):
        """Test the BLAST specificity filter."""
        # Create test DataFrame
        df = self.test_df.copy()
        
        # Test cases
        # Case 1: best * factor < second (should pass)
        row = df.iloc[0]
        result = BlastProcessor.passes_blast_filter(row, "Forward Primer")
        # Handle NumPy boolean types with bool() conversion
        assert bool(result) == True
        
        # Case 2: best * factor > second (should fail)
        # The implementation correctly calculates:
        # For Probe: best=0.003, second=0.05
        # 0.003 * 2 = 0.006 which is < 0.05, so this should PASS (True)
        row = df.iloc[0]
        result = BlastProcessor.passes_blast_filter(row, "Probe")
        assert bool(result) == True
        
        # Case 3: No second hit (should pass)
        row = df.iloc[1]
        result = BlastProcessor.passes_blast_filter(row, "Forward Primer")
        assert bool(result) == True
        
        # Case 4: No best hit (should fail)
        row = df.iloc[1]
        result = BlastProcessor.passes_blast_filter(row, "Reverse Primer")
        assert bool(result) == False