"""
Tests for blast_db_creator module
"""
import os
import pytest
from unittest.mock import patch, MagicMock
import subprocess
import tempfile

from ...utils.blast_db_creator import BlastDBCreator


class TestBlastDBCreator:
    """Tests for BlastDBCreator class."""

    @patch('subprocess.run')
    @patch('os.makedirs')
    @patch('shutil.copy2')
    @patch('os.path.exists')
    @patch('shutil.rmtree')
    @patch('tempfile.mkdtemp')
    def test_create_db(self, mock_mkdtemp, mock_rmtree, mock_exists, mock_copy, mock_makedirs, mock_run, tmp_path):
        """Test creating a BLAST database."""
        # Setup
        mock_exists.return_value = True
        mock_mkdtemp.return_value = os.path.join(tmp_path, "blastdb_temp")
        
        # Configure the subprocess.run mock to simulate successful execution
        # First mock call for makeblastdb
        mock_makeblastdb_process = MagicMock()
        mock_makeblastdb_process.returncode = 0
        mock_makeblastdb_process.stdout = "Database created successfully"
        mock_makeblastdb_process.stderr = ""
        
        # Second mock call for blastdbcmd (verify_db)
        mock_verify_process = MagicMock()
        mock_verify_process.returncode = 0
        mock_verify_process.stdout = "Database verified"
        mock_verify_process.stderr = ""
        
        # Set up the side_effect to return different values for each call
        mock_run.side_effect = [mock_makeblastdb_process, mock_verify_process]
        
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Create a temp fasta file path
        fasta_file = os.path.join(tmp_path, "test.fasta")
        
        # Test creating a database
        output_dir = os.path.join(tmp_path, "blast_db")
        db_path = BlastDBCreator.create_db(fasta_file, output_dir=output_dir, logger=mock_logger)
        
        # Verify
        assert db_path is not None
        assert db_path.endswith("test")
        assert mock_makedirs.called
        assert mock_copy.called
        assert mock_run.called
        assert mock_rmtree.called
        
        # Check that the makeblastdb command was called with correct args
        # Get first call args (makeblastdb)
        first_call_args = mock_run.call_args_list[0][0][0]
        assert "makeblastdb" in first_call_args
        assert "-in" in first_call_args
        assert "-dbtype" in first_call_args
        assert "nucl" in first_call_args
        assert "-out" in first_call_args

    @patch('subprocess.run')
    @patch('os.makedirs')
    @patch('os.path.exists')
    @patch('tempfile.mkdtemp')
    @patch('shutil.copy2')
    def test_create_db_error_handling(self, mock_copy2, mock_mkdtemp, mock_exists, mock_makedirs, mock_run):
        """Test error handling in database creation."""
        # Setup to correctly pass the file existence check
        mock_exists.return_value = True
        mock_mkdtemp.return_value = "/tmp/blastdb_temp"
        
        # Configure the subprocess.run mock to simulate execution failure
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_process.stdout = ""
        mock_process.stderr = "Error creating database"
        mock_run.return_value = mock_process
        
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Test creating a database with failure
        with pytest.raises(Exception) as excinfo:
            BlastDBCreator.create_db("/path/to/fasta.fasta", logger=mock_logger)
        
        # Verify the exception was raised with the correct message
        assert "BLAST database creation failed" in str(excinfo.value)
        assert mock_run.called
        
        # Test file not found case separately
        mock_exists.return_value = False
        
        with pytest.raises(FileNotFoundError):
            BlastDBCreator.create_db("/path/to/nonexistent.fasta", logger=mock_logger)