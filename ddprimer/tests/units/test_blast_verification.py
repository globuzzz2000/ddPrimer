#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for BlastVerification in the ddPrimer pipeline.

These tests verify the functionality of the BlastVerification class,
which is responsible for verifying BLAST databases before pipeline execution.
"""

import os
import pytest
import tempfile
from unittest.mock import patch, MagicMock, call

# Import the module to test
from ...utils import BlastVerification
from ...config import Config


class TestBlastVerification:
    """Test class for BlastVerification functions."""
    
    @pytest.fixture
    def mock_logger(self):
        """Create a mock logger."""
        return MagicMock()
    
    @pytest.fixture
    def setup_config(self):
        """Set up Config with a test database path."""
        # Save the original DB_PATH
        original_db_path = Config.DB_PATH
        
        # Set a test DB_PATH
        Config.DB_PATH = "/path/to/test_db"
        
        yield
        
        # Restore the original DB_PATH
        Config.DB_PATH = original_db_path
    
    @patch('ddprimer.utils.blast_verification.subprocess.run')
    @patch('os.path.exists')
    def test_verify_blast_database_success(self, mock_exists, mock_run, mock_logger, setup_config):
        """Test successful BLAST database verification."""
        # Setup mocks
        mock_exists.return_value = True  # All files exist
        mock_run.return_value = MagicMock(returncode=0, stdout="BLAST search completed", stderr="")
        
        # Verify the database
        result = BlastVerification.verify_blast_database(mock_logger)
        
        # Verify the result
        assert result is True, "Should return True for valid database"
        assert mock_exists.called, "Should check if files exist"
        assert mock_run.called, "Should call subprocess.run for BLAST test"
    
    @patch('ddprimer.utils.blast_verification.subprocess.run')
    @patch('os.path.exists')
    def test_verify_blast_database_missing_files(self, mock_exists, mock_run, mock_logger, setup_config):
        """Test BLAST database verification with missing files."""
        # Setup mocks
        mock_exists.return_value = False  # Files don't exist
        
        # Mock _handle_failed_verification to return False
        with patch('ddprimer.utils.blast_verification.BlastVerification._handle_failed_verification', return_value=False):
            # Verify the database
            result = BlastVerification.verify_blast_database(mock_logger)
            
            # Verify the result
            assert result is False, "Should return False for missing files"
            assert mock_exists.called, "Should check if files exist"
            assert not mock_run.called, "Should not call subprocess.run for BLAST test"
    
    @patch('ddprimer.utils.blast_verification.subprocess.run')
    @patch('os.path.exists')
    def test_verify_blast_database_blast_failure(self, mock_exists, mock_run, mock_logger, setup_config):
        """Test BLAST database verification with BLAST command failure."""
        # Setup mocks
        mock_exists.return_value = True  # All files exist
        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="Error: blast failed")
        
        # Mock _handle_failed_verification to return False
        with patch('ddprimer.utils.blast_verification.BlastVerification._handle_failed_verification', return_value=False):
            # Mock _retry_verify_with_blastdbcmd to also return False
            with patch('ddprimer.utils.blast_verification.BlastVerification._retry_verify_with_blastdbcmd', return_value=False):
                # Verify the database
                result = BlastVerification.verify_blast_database(mock_logger)
                
                # Verify the result
                assert result is False, "Should return False for BLAST failure"
                assert mock_exists.called, "Should check if files exist"
                assert mock_run.called, "Should call subprocess.run for BLAST test"
    
    @patch('ddprimer.utils.blast_verification.subprocess.run')
    def test_retry_verify_with_blastdbcmd_success(self, mock_run, mock_logger):
        """Test successful retry verification with blastdbcmd."""
        # Setup mock
        mock_run.return_value = MagicMock(returncode=0, stdout="Database info", stderr="")
        
        # Retry verification
        result = BlastVerification._retry_verify_with_blastdbcmd("/path/to/db", mock_logger)
        
        # Verify the result
        assert result is True, "Should return True for successful retry"
        assert mock_run.called, "Should call subprocess.run for blastdbcmd"
    
    @patch('ddprimer.utils.blast_verification.subprocess.run')
    def test_retry_verify_with_blastdbcmd_failure(self, mock_run, mock_logger):
        """Test failed retry verification with blastdbcmd."""
        # Setup mock
        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="Error: invalid database")
        
        # Mock _handle_failed_verification to return False
        with patch('ddprimer.utils.blast_verification.BlastVerification._handle_failed_verification', return_value=False):
            # Retry verification
            result = BlastVerification._retry_verify_with_blastdbcmd("/path/to/db", mock_logger)
            
            # Verify the result
            assert result is False, "Should return False for failed retry"
            assert mock_run.called, "Should call subprocess.run for blastdbcmd"
    
    @patch('builtins.input')
    def test_handle_failed_verification_exit(self, mock_input, mock_logger):
        """Test handling failed verification with exit choice."""
        # Setup mock
        mock_input.return_value = "4"  # Exit option
        
        # Handle failed verification
        result = BlastVerification._handle_failed_verification(mock_logger)
        
        # Verify the result
        assert result is False, "Should return False for exit choice"
        assert mock_input.called, "Should prompt for input"
        assert mock_logger.error.called, "Should log error message"
    
    @patch('builtins.input')
    def test_handle_failed_verification_invalid_choice(self, mock_input, mock_logger):
        """Test handling failed verification with invalid choice."""
        # Setup mocks
        mock_input.return_value = "invalid"
        
        # Mock _handle_failed_verification to return False on recursive call
        with patch('ddprimer.utils.blast_verification.BlastVerification._handle_failed_verification', return_value=False):
            # Handle failed verification
            result = BlastVerification._handle_failed_verification(mock_logger)
            
            # Verify the result
            assert result is False, "Should return False after recursive call"
            assert mock_input.called, "Should prompt for input"
            assert mock_logger.error.called, "Should log error message"
    
    @patch('builtins.input')
    @patch('ddprimer.utils.blast_verification.ModelOrganismManager.select_model_organism')
    @patch('ddprimer.utils.blast_verification.BlastDBCreator')
    def test_handle_failed_verification_create_from_model_organism(self, mock_blast_db_creator, mock_select, mock_input, mock_logger):
        """Test handling failed verification with creating from model organism."""
        # Setup mocks
        mock_input.return_value = "1"  # Create from model organism
        mock_select.return_value = ("key", "name", "/path/to/fasta")
        mock_blast_db_creator.return_value.create_database.return_value = "/path/to/new_db"
        
        # Mock Config.save_database_config to do nothing
        with patch('ddprimer.utils.blast_verification.Config.save_database_config'):
            # Handle failed verification
            result = BlastVerification._handle_failed_verification(mock_logger)
            
            # Verify the result
            assert result is True, "Should return True for successful DB creation"
            assert mock_input.called, "Should prompt for input"
            assert mock_select.called, "Should call select_model_organism"
            assert mock_blast_db_creator.return_value.create_database.called, "Should call create_database"
    
    @patch('builtins.input')
    @patch('ddprimer.utils.blast_verification.BlastDBCreator')
    def test_handle_failed_verification_create_from_custom_file(self, mock_blast_db_creator, mock_input, mock_logger):
        """Test handling failed verification with creating from custom file."""
        # Setup mocks
        mock_input.return_value = "2"  # Create from custom file
        mock_blast_db_creator.return_value.create_database.return_value = "/path/to/new_db"
        
        # Mock Config.save_database_config to do nothing
        with patch('ddprimer.utils.blast_verification.Config.save_database_config'):
            # Handle failed verification
            result = BlastVerification._handle_failed_verification(mock_logger)
            
            # Verify the result
            assert result is True, "Should return True for successful DB creation"
            assert mock_input.called, "Should prompt for input"
            assert mock_blast_db_creator.return_value.create_database.called, "Should call create_database"
    
    @patch('builtins.input')
    def test_handle_failed_verification_try_again(self, mock_input, mock_logger):
        """Test handling failed verification with try again choice."""
        # Setup mock
        mock_input.return_value = "3"  # Try again option
        
        # Handle failed verification
        result = BlastVerification._handle_failed_verification(mock_logger)
        
        # Verify the result
        assert result is False, "Should return False for try again choice"
        assert mock_input.called, "Should prompt for input"
        assert mock_logger.info.called, "Should log info message"