#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for BLAST database utilities in the ddPrimer pipeline.

These tests verify the functionality of the BlastDBCreator and BlastVerification
classes, which are responsible for creating and verifying BLAST databases.
"""

import os
import sys
import pytest
import tempfile
import shutil
import subprocess
from unittest.mock import patch, MagicMock, call

# Import the modules to test
from ...utils import BlastDBCreator
from ...utils import BlastVerification
from ...config import Config


class TestBlastDBCreator:
    """Test class for BlastDBCreator functions."""
    
    @pytest.fixture
    def mock_logger(self):
        """Create a mock logger."""
        return MagicMock()
    
    @pytest.fixture
    def temp_fasta_file(self):
        """Create a temporary FASTA file with test sequences."""
        temp_dir = tempfile.mkdtemp()
        fasta_path = os.path.join(temp_dir, "test_genome.fasta")
        
        # Create a small test FASTA file with realistic sequences
        with open(fasta_path, 'w') as f:
            f.write(">chr1\n")
            f.write("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")
            f.write(">chr2\n")
            f.write("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
        
        yield fasta_path
        
        # Clean up
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @patch('ddprimer.utils.blast_db_creator.subprocess.run')
    def test_create_db_success(self, mock_run, temp_fasta_file, temp_output_dir, mock_logger):
        """Test successful creation of a BLAST database."""
        # Setup mock subprocess.run to return success
        mock_run.return_value = MagicMock(returncode=0, stdout="Database created successfully", stderr="")
        
        # Mock _verify_db to return True
        with patch('ddprimer.utils.blast_db_creator.BlastDBCreator._verify_db', return_value=True):
            # Create the database
            db_creator = BlastDBCreator()
            db_path = db_creator.create_db(temp_fasta_file, temp_output_dir, "test_db", mock_logger)
            
            # Verify the result
            assert db_path == os.path.join(temp_output_dir, "test_db"), "Should return correct DB path"
            
            # Verify subprocess.run was called correctly
            mock_run.assert_called_once()
            args, kwargs = mock_run.call_args
            assert args[0][0] == "makeblastdb", "Should call makeblastdb"
            
            # Modified: Check that the command uses a file with the same basename,
            # since the actual file is copied to a temporary directory
            assert os.path.basename(args[0][2]) == os.path.basename(temp_fasta_file), "Should use file with same name"
            assert args[0][4] == "nucl", "Should create nucleotide database"
    
    @patch('ddprimer.utils.blast_db_creator.subprocess.run')
    def test_create_db_failure(self, mock_run, temp_fasta_file, temp_output_dir, mock_logger):
        """Test failure in creating a BLAST database."""
        # Setup mock subprocess.run to return failure
        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="Error: invalid input file")
        
        # Create the database - should raise an exception
        db_creator = BlastDBCreator()
        with pytest.raises(Exception, match="BLAST database creation failed"):
            db_creator.create_db(temp_fasta_file, temp_output_dir, "test_db", mock_logger)
    
    def test_create_db_nonexistent_file(self, temp_output_dir, mock_logger):
        """Test creating a database from a non-existent file."""
        # Try to create a database from a non-existent file
        db_creator = BlastDBCreator()
        with pytest.raises(FileNotFoundError):
            db_creator.create_db("nonexistent_file.fasta", temp_output_dir, "test_db", mock_logger)
    
    @patch('ddprimer.utils.blast_db_creator.subprocess.run')
    @patch('ddprimer.utils.blast_db_creator.ModelOrganismManager.select_model_organism')
    def test_create_database_with_model_organism(self, mock_select, mock_run, temp_fasta_file, temp_output_dir, mock_logger):
        """Test creating a database from a model organism."""
        # Setup mocks
        mock_select.return_value = ("Thale cress", "Arabidopsis thaliana", temp_fasta_file)
        mock_run.return_value = MagicMock(returncode=0, stdout="Database created successfully", stderr="")
        
        # Mock _verify_db to return True
        with patch('ddprimer.utils.blast_db_creator.BlastDBCreator._verify_db', return_value=True):
            # Create the database
            db_creator = BlastDBCreator()
            db_path = db_creator.create_database(None, "arabidopsis_db", temp_output_dir, mock_logger)
            
            # Verify the result
            assert db_path == os.path.join(temp_output_dir, "arabidopsis_db"), "Should return correct DB path"
            assert mock_select.called, "Should call select_model_organism"
    
    @patch('ddprimer.utils.blast_db_creator.ModelOrganismManager.select_model_organism')
    def test_create_database_selection_canceled(self, mock_select, mock_logger):
        """Test canceling model organism selection."""
        # Setup mock to simulate canceled selection
        mock_select.return_value = (None, None, None)
        
        # Create the database
        db_creator = BlastDBCreator()
        db_path = db_creator.create_database(None, None, None, mock_logger)
        
        # Verify the result
        assert db_path is None, "Should return None when selection is canceled"
        assert mock_select.called, "Should call select_model_organism"
    
    def test_verify_db_success(self, mock_logger):
        """Test successful database verification."""
        # Mock os.path.exists to return True for all required files
        with patch('os.path.exists', return_value=True):
            # Mock subprocess.run to return success
            with patch('subprocess.run') as mock_run:
                mock_run.return_value = MagicMock(returncode=0, stdout="Database info", stderr="")
                
                # Verify the database
                result = BlastDBCreator._verify_db("valid_db_path", mock_logger)
                
                # Verify the result
                assert result is True, "Should return True for valid database"
                assert mock_run.called, "Should call subprocess.run"
    
    def test_verify_db_missing_files(self, mock_logger):
        """Test database verification with missing files."""
        # Setup mock for checking if files exist
        with patch('os.path.exists', return_value=False):
            # Verify the database
            result = BlastDBCreator._verify_db("invalid_db_path", mock_logger)
            
            # Verify the result
            assert result is False, "Should return False for missing files"
    
    def test_verify_db_command_failure(self, mock_logger):
        """Test database verification with command failure."""
        # First mock os.path.exists to return True for all files
        with patch('os.path.exists', return_value=True):
            # Then mock subprocess.run to return failure
            with patch('subprocess.run') as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="Error: invalid database")
                
                # Verify the database
                result = BlastDBCreator._verify_db("invalid_db_path", mock_logger)
                
                # Verify the result
                assert result is False, "Should return False for command failure"
                assert mock_run.called, "Should call subprocess.run"


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
        # First call returns invalid, second call returns 4 (exit)
        mock_input.side_effect = ["invalid", "4"]
        
        # Handle failed verification
        result = BlastVerification._handle_failed_verification(mock_logger)
        
        # Verify the result
        assert result is False, "Should return False after recursive call"
        assert mock_input.call_count == 2, "Should prompt for input twice"
        assert mock_logger.error.called, "Should log error message"