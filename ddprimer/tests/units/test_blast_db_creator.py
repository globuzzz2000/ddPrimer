#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for BlastDBCreator in the ddPrimer pipeline.

These tests verify the functionality of the BlastDBCreator class,
which is responsible for creating BLAST databases from FASTA files.
"""

import os
import pytest
import tempfile
import shutil
from unittest.mock import patch, MagicMock, call

# Import the module to test
from ...utils import BlastDBCreator


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
            assert args[0][2] == temp_fasta_file, "Should use correct input file"
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