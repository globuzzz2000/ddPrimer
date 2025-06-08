#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the database selector module in the ddPrimer pipeline.

These tests verify the functionality of the DatabaseSelector class,
which is responsible for finding and selecting existing BLAST databases.
"""

import pytest
from unittest.mock import patch, MagicMock

# Import package modules
from ...utils import DatabaseSelector


class TestDatabaseSelector:
    """Test class for DatabaseSelector functions."""
    
    @pytest.fixture
    def mock_logger(self):
        """Create a mock logger."""
        return MagicMock()
    
    @patch('os.path.exists')
    @patch('os.listdir')
    def test_find_existing_databases_success(self, mock_listdir, mock_exists, mock_logger):
        """Test finding existing databases when they exist."""
        # Setup mocks
        mock_exists.return_value = True
        
        # Different .nhr files for different locations
        mock_listdir.side_effect = [
            ["db1.nhr", "db2.nhr", "not_a_db.txt"],  # First directory
            ["db3.nhr", "db4.nhr"]                    # Second directory
        ]
        
        # Call the function
        databases = DatabaseSelector.find_existing_databases(mock_logger)
        
        # Verify the result
        assert len(databases) == 4, "Should find 4 databases"
        
        # Check specific paths and display names
        for db_path, display_name in databases.items():
            assert db_path.endswith("db1") or db_path.endswith("db2") or \
                   db_path.endswith("db3") or db_path.endswith("db4"), \
                   "Should have correct database paths"
            
            assert display_name in ["db1", "db2", "db3", "db4"], \
                   "Should have correct display names"
    
    @patch('os.path.exists')
    def test_find_existing_databases_no_directories(self, mock_exists, mock_logger):
        """Test finding databases when no directories exist."""
        # Setup mock to return False for all directories
        mock_exists.return_value = False
        
        # Call the function
        databases = DatabaseSelector.find_existing_databases(mock_logger)
        
        # Verify the result
        assert len(databases) == 0, "Should find no databases"
    
    @patch('os.path.exists')
    @patch('os.listdir')
    def test_find_existing_databases_error(self, mock_listdir, mock_exists, mock_logger):
        """Test handling errors when finding databases."""
        # Setup mocks
        mock_exists.return_value = True
        mock_listdir.side_effect = PermissionError("Permission denied")
        
        # Call the function - should handle the error gracefully
        databases = DatabaseSelector.find_existing_databases(mock_logger)
        
        # Verify the result
        assert len(databases) == 0, "Should find no databases on error"
        assert mock_logger.warning.called, "Should log a warning message"
    
    @patch('ddprimer.utils.db_selector.DatabaseSelector.find_existing_databases')
    @patch('builtins.input')
    def test_select_database_success(self, mock_input, mock_find, mock_logger):
        """Test successful database selection."""
        # Setup mock to return some databases
        mock_find.return_value = {
            "/path/to/db1": "Database 1",
            "/path/to/db2": "Database 2"
        }
        
        # Mock user input to select the first database
        mock_input.return_value = "1"
        
        # Call the function
        result = DatabaseSelector.select_database(mock_logger)
        
        # Verify the result
        assert result == "/path/to/db1", "Should return the path of the selected database"
        assert mock_find.called, "Should call find_existing_databases"
        assert mock_input.called, "Should prompt for input"
        assert mock_logger.info.called, "Should log information messages"
    
    @patch('ddprimer.utils.db_selector.DatabaseSelector.find_existing_databases')
    @patch('builtins.input')
    def test_select_database_cancel(self, mock_input, mock_find, mock_logger):
        """Test canceling database selection."""
        # Setup mock to return some databases
        mock_find.return_value = {
            "/path/to/db1": "Database 1",
            "/path/to/db2": "Database 2"
        }
        
        # Mock user input to select the cancel option (3 in this case)
        mock_input.return_value = "3"
        
        # Call the function
        result = DatabaseSelector.select_database(mock_logger)
        
        # Verify the result
        assert result is None, "Should return None when canceled"
        assert mock_find.called, "Should call find_existing_databases"
        assert mock_input.called, "Should prompt for input"
        assert mock_logger.info.called, "Should log information messages"
    
    @patch('ddprimer.utils.db_selector.DatabaseSelector.find_existing_databases')
    @patch('builtins.input')
    def test_select_database_invalid_choice(self, mock_input, mock_find, mock_logger):
        """Test invalid selection."""
        # Setup mock to return some databases
        mock_find.return_value = {
            "/path/to/db1": "Database 1",
            "/path/to/db2": "Database 2"
        }
        
        # Mock user input to select an invalid option
        mock_input.return_value = "5"  # Out of range
        
        # Call the function
        result = DatabaseSelector.select_database(mock_logger)
        
        # Verify the result
        assert result is None, "Should return None for invalid choice"
        assert mock_find.called, "Should call find_existing_databases"
        assert mock_input.called, "Should prompt for input"
        assert mock_logger.error.called, "Should log error message"
    
    @patch('ddprimer.utils.db_selector.DatabaseSelector.find_existing_databases')
    @patch('builtins.input')
    def test_select_database_non_numeric_choice(self, mock_input, mock_find, mock_logger):
        """Test non-numeric input."""
        # Setup mock to return some databases
        mock_find.return_value = {
            "/path/to/db1": "Database 1",
            "/path/to/db2": "Database 2"
        }
        
        # Mock user input to provide non-numeric input
        mock_input.return_value = "not a number"
        
        # Call the function
        result = DatabaseSelector.select_database(mock_logger)
        
        # Verify the result
        assert result is None, "Should return None for non-numeric choice"
        assert mock_find.called, "Should call find_existing_databases"
        assert mock_input.called, "Should prompt for input"
        assert mock_logger.error.called, "Should log error message"
    
    @patch('ddprimer.utils.db_selector.DatabaseSelector.find_existing_databases')
    def test_select_database_no_databases(self, mock_find, mock_logger):
        """Test when no databases are found."""
        # Setup mock to return empty dictionary
        mock_find.return_value = {}
        
        # Call the function
        result = DatabaseSelector.select_database(mock_logger)
        
        # Verify the result
        assert result is None, "Should return None when no databases are found"
        assert mock_find.called, "Should call find_existing_databases"
        assert mock_logger.info.called, "Should log information message"
    
    @patch('ddprimer.utils.db_selector.DatabaseSelector.find_existing_databases')
    @patch('builtins.input')
    def test_select_database_keyboard_interrupt(self, mock_input, mock_find, mock_logger):
        """Test handling keyboard interrupt during selection."""
        # Setup mock to return some databases
        mock_find.return_value = {
            "/path/to/db1": "Database 1",
            "/path/to/db2": "Database 2"
        }
        
        # Mock input to raise KeyboardInterrupt
        mock_input.side_effect = KeyboardInterrupt()
        
        # Call the function
        result = DatabaseSelector.select_database(mock_logger)
        
        # Verify the result
        assert result is None, "Should return None on keyboard interrupt"
        assert mock_find.called, "Should call find_existing_databases"
        assert mock_input.called, "Should prompt for input"
        assert mock_logger.info.called, "Should log information message"