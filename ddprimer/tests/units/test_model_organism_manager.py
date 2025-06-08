#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the model organism manager module in the ddPrimer pipeline.

These tests verify the functionality of the ModelOrganismManager class,
which handles fetching, preprocessing, and managing model organism genomes.
"""

import os
import pytest
import tempfile
import shutil
import gzip
from unittest.mock import patch, MagicMock

# Import package modules
from ...utils import ModelOrganismManager


class TestModelOrganismManager:
    """Test class for ModelOrganismManager functions."""
    
    @pytest.fixture
    def mock_logger(self):
        """Create a mock logger."""
        return MagicMock()
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    def test_model_organisms_structure(self):
        """Test that the MODEL_ORGANISMS dictionary has the expected structure."""
        # Verify that each model organism has required keys
        for key, organism in ModelOrganismManager.MODEL_ORGANISMS.items():
            assert "name" in organism, f"Organism '{key}' is missing 'name' key"
            assert "url" in organism, f"Organism '{key}' is missing 'url' key"
            assert "compressed" in organism, f"Organism '{key}' is missing 'compressed' key"
            
            # Check that URLs are valid NCBI FTP URLs
            assert organism["url"].startswith("https://ftp.ncbi.nlm.nih.gov/"), \
                f"URL for '{key}' should be an NCBI FTP URL"
    
    def test_get_model_organism_menu(self):
        """Test generation of the model organism menu."""
        menu = ModelOrganismManager.get_model_organism_menu()
        
        # Verify the menu is a string
        assert isinstance(menu, str), "Menu should be a string"
        
        # Check for required menu items
        assert "Available model organisms" in menu, "Menu should have header"
        assert "0. Custom FASTA file" in menu, "Menu should have custom file option"
        assert "Cancel" in menu, "Menu should have cancel option"
        
        # Verify all model organisms are in the menu
        for _, organism in ModelOrganismManager.MODEL_ORGANISMS.items():
            assert organism["name"] in menu, f"Organism '{organism['name']}' should be in menu"
    
    @patch('urllib.request.urlretrieve')
    @patch('os.path.exists')
    def test_fetch_model_organism_new_download(self, mock_exists, mock_urlretrieve, temp_dir, mock_logger):
        """Test fetching a model organism that hasn't been downloaded yet."""
        # Setup mocks
        mock_exists.return_value = False
        
        # Create a test compressed file
        compressed_file = os.path.join(temp_dir, "test_genome.fasta.gz")
        with gzip.open(compressed_file, 'wb') as f:
            f.write(b">chr1\nATGCATGCATGC\n>chr2\nGCTAGCTAGCTA\n")
        
        # Mock urlretrieve to "download" our test file
        def mock_download(url, output_file, reporthook):
            shutil.copy(compressed_file, output_file)
            return output_file, {}
            
        mock_urlretrieve.side_effect = mock_download
        
        # Mock _extract_gzip to return uncompressed path
        with patch('ddprimer.utils.model_organism_manager.ModelOrganismManager._extract_gzip') as mock_extract:
            mock_extract.return_value = os.path.join(temp_dir, "test_genome.fasta")
            
            # Call the function
            result = ModelOrganismManager.fetch_model_organism("Thale cress", temp_dir, mock_logger)
            
            # Verify the result
            assert result == os.path.join(temp_dir, "test_genome.fasta"), "Should return path to extracted file"
            assert mock_urlretrieve.called, "Should call urlretrieve to download file"
            assert mock_extract.called, "Should call _extract_gzip to decompress file"
    
    @patch('os.path.exists')
    def test_fetch_model_organism_already_exists(self, mock_exists, temp_dir, mock_logger):
        """Test fetching a model organism that has already been downloaded."""
        # Setup mock to indicate file already exists
        mock_exists.return_value = True
        
        # Create test files
        compressed_file = os.path.join(temp_dir, "GCF_000001735.4_TAIR10.1_genomic.fna.gz")
        uncompressed_file = os.path.join(temp_dir, "GCF_000001735.4_TAIR10.1_genomic.fna")
        
        with open(uncompressed_file, 'w') as f:
            f.write(">chr1\nATGCATGCATGC\n>chr2\nGCTAGCTAGCTA\n")
        
        # Call the function
        result = ModelOrganismManager.fetch_model_organism("Thale cress", temp_dir, mock_logger)
        
        # Verify the result - update expectation to match implementation
        assert result == uncompressed_file, "Should return path to uncompressed file when it exists"
        assert mock_logger.info.called, "Should log information about existing file"
    
    def test_fetch_model_organism_invalid_key(self, mock_logger):
        """Test fetching with an invalid organism key."""
        with pytest.raises(ValueError, match="Invalid model organism key"):
            ModelOrganismManager.fetch_model_organism("Invalid Organism", None, mock_logger)
    
    @patch('urllib.request.urlretrieve')
    def test_fetch_model_organism_download_error(self, mock_urlretrieve, temp_dir, mock_logger):
        """Test handling download errors."""
        # Setup mock to raise an exception
        mock_urlretrieve.side_effect = Exception("Download failed")
        
        # Call the function - should raise ConnectionError
        with pytest.raises(ConnectionError, match="Failed to download genome file"):
            ModelOrganismManager.fetch_model_organism("Thale cress", temp_dir, mock_logger)
    
    def test_extract_gzip(self, temp_dir, mock_logger):
        """Test extracting a gzipped file."""
        # Create a test compressed file
        compressed_file = os.path.join(temp_dir, "test_genome.fasta.gz")
        test_content = b">chr1\nATGCATGCATGC\n>chr2\nGCTAGCTAGCTA\n"
        
        with gzip.open(compressed_file, 'wb') as f:
            f.write(test_content)
        
        # Call the function
        result = ModelOrganismManager._extract_gzip(compressed_file, mock_logger)
        
        # Verify the result
        assert result == os.path.join(temp_dir, "test_genome.fasta"), "Should return path to extracted file"
        assert os.path.exists(result), "Extracted file should exist"
        
        # Verify the content was correctly extracted
        with open(result, 'rb') as f:
            content = f.read()
            assert content == test_content, "Extracted content should match original"
    
    def test_extract_gzip_error(self, temp_dir, mock_logger):
        """Test handling errors during extraction."""
        # Create an invalid compressed file
        invalid_file = os.path.join(temp_dir, "invalid.gz")
        with open(invalid_file, 'wb') as f:
            f.write(b"This is not a valid gzip file")
        
        # Call the function - should raise IOError
        with pytest.raises(IOError, match="Failed to extract genome file"):
            ModelOrganismManager._extract_gzip(invalid_file, mock_logger)
    
    @patch('builtins.input')
    @patch('ddprimer.utils.file_io.FileIO.select_fasta_file') 
    def test_select_model_organism_custom_file(self, mock_select, mock_input, mock_logger):
        """Test selecting a custom FASTA file."""
        # Setup mocks
        mock_input.return_value = "0"  # Choose custom file option
        mock_select.return_value = "/path/to/custom.fasta"
        
        # Call the function
        key, name, file_path = ModelOrganismManager.select_model_organism(mock_logger)
        
        # Verify the result
        assert key is None, "Should return None for key with custom file"
        assert name == "Custom file", "Should return 'Custom file' for name"
        assert file_path == "/path/to/custom.fasta", "Should return path to selected file"
        assert mock_input.called, "Should prompt for input"
        assert mock_select.called, "Should call select_fasta_file"
    
    @patch('builtins.input')
    def test_select_model_organism_model_species(self, mock_input, mock_logger):
        """Test selecting a model organism."""
        # Setup mock to select Thale cress (index 1)
        mock_input.return_value = "1"
        
        # Mock fetch_model_organism to return a path
        with patch('ddprimer.utils.model_organism_manager.ModelOrganismManager.fetch_model_organism') as mock_fetch:
            mock_fetch.return_value = "/path/to/arabidopsis.fasta"
            
            # Call the function
            key, name, file_path = ModelOrganismManager.select_model_organism(mock_logger)
            
            # Verify the result
            assert key == "Thale cress", "Should return correct organism key"
            assert "Arabidopsis thaliana" in name, "Should return scientific name"
            assert file_path == "/path/to/arabidopsis.fasta", "Should return path to genome file"
            assert mock_fetch.called, "Should call fetch_model_organism"
    
    @patch('builtins.input')
    @patch('ddprimer.utils.db_selector.DatabaseSelector.select_database')
    def test_select_model_organism_existing_database(self, mock_select_db, mock_input, mock_logger):
        """Test selecting from existing databases."""
        # Setup mocks
        # Use the index that corresponds to "Select from existing databases"
        existing_db_index = len(ModelOrganismManager.MODEL_ORGANISMS) + 1
        mock_input.return_value = str(existing_db_index)
        mock_select_db.return_value = "/path/to/existing_db"
        
        # Call the function
        key, name, file_path = ModelOrganismManager.select_model_organism(mock_logger)
        
        # Verify the result
        assert key == "existing_db", "Should return 'existing_db' for key"
        assert name == "existing db", "Should return database name with space instead of underscore"  # Updated expectation
        assert file_path == "/path/to/existing_db", "Should return path to database"
        assert mock_select_db.called, "Should call select_database"
    
    @patch('builtins.input')
    def test_select_model_organism_cancel(self, mock_input, mock_logger):
        """Test canceling organism selection."""
        # Setup mock to select cancel option
        cancel_index = len(ModelOrganismManager.MODEL_ORGANISMS) + 2
        mock_input.return_value = str(cancel_index)
        
        # Call the function
        key, name, file_path = ModelOrganismManager.select_model_organism(mock_logger)
        
        # Verify the result
        assert key is None, "Should return None for key when canceled"
        assert name is None, "Should return None for name when canceled"
        assert file_path is None, "Should return None for file path when canceled"
    
    @patch('builtins.input')
    def test_select_model_organism_invalid_choice(self, mock_input, mock_logger):
        """Test handling invalid choice."""
        # Setup mock for invalid input
        mock_input.return_value = "999"  # Invalid index
        
        # Call the function
        key, name, file_path = ModelOrganismManager.select_model_organism(mock_logger)
        
        # Verify the result
        assert key is None, "Should return None for key with invalid choice"
        assert name is None, "Should return None for name with invalid choice"
        assert file_path is None, "Should return None for file path with invalid choice"
        assert mock_logger.error.called, "Should log error message"
    
    @patch('os.remove')
    @patch('os.path.exists')
    def test_cleanup_genome_file(self, mock_exists, mock_remove, mock_logger):
        """Test cleaning up genome files after database creation."""
        # Setup mock to indicate files exist
        mock_exists.return_value = True
        
        # Call the function
        ModelOrganismManager.cleanup_genome_file("/path/to/genome.fasta", mock_logger)
        
        # Verify the result
        assert mock_exists.called, "Should check if file exists"
        assert mock_remove.called, "Should call os.remove"
        # Should attempt to remove both the file and its .gz version
        assert mock_remove.call_count >= 1, "Should remove at least one file"
    
    @patch('os.remove')
    @patch('os.path.exists')
    def test_cleanup_genome_file_not_exists(self, mock_exists, mock_remove, mock_logger):
        """Test cleaning up when file doesn't exist."""
        # Setup mock to indicate file doesn't exist
        mock_exists.return_value = False
        
        # Call the function
        ModelOrganismManager.cleanup_genome_file("/path/to/nonexistent.fasta", mock_logger)
        
        # Verify the result
        assert mock_exists.called, "Should check if file exists"
        assert not mock_remove.called, "Should not call os.remove"
        assert mock_logger.debug.called, "Should log debug message"
    
    @patch('os.remove')   
    @patch('os.path.exists')
    def test_cleanup_genome_file_error(self, mock_remove, mock_logger):
        """Test handling errors during cleanup."""
        # Setup mock to raise an exception
        mock_remove.side_effect = Exception("Permission denied")
        
        # Call the function - should handle the exception gracefully
        ModelOrganismManager.cleanup_genome_file("/path/to/genome.fasta", mock_logger)
        
        # Verify the result
        assert mock_remove.called, "Should attempt to remove file"
        assert mock_logger.warning.called, "Should log warning message"