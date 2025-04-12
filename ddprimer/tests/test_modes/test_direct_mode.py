"""
Tests for direct mode functionality
"""
import os
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

from ...modes.direct_mode import run_direct_mode, load_direct_sequences
from ...config import Config


class TestDirectMode:
    """Tests for direct mode workflow."""

    def test_load_direct_sequences_excel(self, create_mini_direct_input):
        """Test loading sequences from Excel file."""
        # Create a test Excel file
        excel_path = create_mini_direct_input()
        
        # Run the function
        sequences = load_direct_sequences(excel_path)
        
        # Validate the results
        assert isinstance(sequences, dict)
        assert len(sequences) == 3  # Should have 3 sequences
        assert 'seq1' in sequences
        assert 'seq2' in sequences
        assert 'seq3' in sequences
        
        # Check sequence content
        assert sequences['seq1'].startswith('ATGCGGCCACTG')
        assert sequences['seq2'].startswith('CTAGCCGGTAAT')
        assert sequences['seq3'].startswith('ACGTACGTACGT')

    def test_load_direct_sequences_csv(self, tmp_path):
        """Test loading sequences from CSV file."""
        # Create a test CSV file
        csv_path = tmp_path / "test_sequences.csv"
        df = pd.DataFrame({
            'Sequence Name': ['seq1', 'seq2', 'seq3'],
            'Sequence': [
                'ATGCGGCCACTGGCTTAA',
                'CTAGCCGGTAATCGGGCC',
                'ACGTACGTACGTGGCCAA'
            ]
        })
        df.to_csv(csv_path, index=False)
        
        # Run the function
        sequences = load_direct_sequences(str(csv_path))
        
        # Validate the results
        assert isinstance(sequences, dict)
        assert len(sequences) == 3  # Should have 3 sequences
        assert 'seq1' in sequences
        assert 'seq2' in sequences
        assert 'seq3' in sequences
        
        # Check sequence content
        assert sequences['seq1'] == 'ATGCGGCCACTGGCTTAA'
        assert sequences['seq2'] == 'CTAGCCGGTAATCGGGCC'
        assert sequences['seq3'] == 'ACGTACGTACGTGGCCAA'

    @patch('ddPrimer.utils.FileUtils.check_and_create_output_dir')
    @patch('ddPrimer.modes.direct_mode.load_direct_sequences')
    @patch('ddPrimer.modes.common.run_primer_design_workflow')
    def test_direct_mode_workflow(
        self,
        mock_run_workflow,
        mock_load_sequences,
        mock_create_output_dir,
        create_mini_direct_input,
        tmp_path
    ):
        """Test the direct mode workflow."""
        # Create a test file
        direct_input_path = create_mini_direct_input()
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure the mocks
        mock_create_output_dir.return_value = output_dir
        mock_load_sequences.return_value = {
            'seq1': 'ATGCGGCCACTGGCTTAA',
            'seq2': 'CTAGCCGGTAATCGGGCC',
            'seq3': 'ACGTACGTACGTGGCCAA'
        }
        mock_run_workflow.return_value = True
        
        # Create a mock args object
        args = MagicMock()
        args.direct = direct_input_path
        args.output = output_dir
        args.cli = True
        
        # Run the direct mode workflow
        result = run_direct_mode(args)
        
        # Verify the results
        assert result is True
        
        # Check that each step of the workflow was called
        mock_create_output_dir.assert_called_once()
        mock_load_sequences.assert_called_once_with(direct_input_path)
        mock_run_workflow.assert_called_once()
        
        # Check that run_primer_design_workflow was called with correct parameters
        workflow_args = mock_run_workflow.call_args[0]
        assert isinstance(workflow_args[0], dict)  # sequences
        assert workflow_args[1] == output_dir  # output_dir
        assert workflow_args[2] == direct_input_path  # reference_file
        assert workflow_args[3] == 'direct'  # mode
        assert workflow_args[4] is None  # genes should be None for direct mode
        
    def test_direct_mode_error_handling(self):
        """Test that direct mode handles errors gracefully."""
        # Create a mock args object with invalid inputs
        args = MagicMock()
        args.direct = "nonexistent.xlsx"
        args.output = "/path/to/output"
        args.cli = True
        
        # Run the direct mode workflow
        result = run_direct_mode(args)
        
        # The workflow should handle errors and return False
        assert result is False