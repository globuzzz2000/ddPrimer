"""
Tests for file_utils module
"""
import os
import sys
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock, mock_open

from ...utils.file_utils import FileUtils


class TestFileUtils:
    """Tests for FileUtils class."""

    def test_check_and_create_output_dir(self, tmp_path):
        """Test checking and creating output directory."""
        # Create a base path
        base_path = str(tmp_path / "output")
        
        # Test with an existing directory
        os.makedirs(base_path, exist_ok=True)
        result_dir = FileUtils.check_and_create_output_dir(base_path)
        
        # Verify the result
        assert result_dir == base_path
        assert os.path.isdir(result_dir)
        
        # Test with a new directory
        new_path = str(tmp_path / "new_output")
        result_dir = FileUtils.check_and_create_output_dir(new_path)
        
        # Verify the result
        assert result_dir == new_path
        assert os.path.isdir(result_dir)
        
        # Test with a timestamp directory
        timestamp_path = os.path.join(base_path, "results_")
        result_dir = FileUtils.check_and_create_output_dir(timestamp_path, add_timestamp=True)
        
        # Verify the result contains a timestamp
        assert result_dir.startswith(timestamp_path)
        assert os.path.isdir(result_dir)

    @patch('ddprimer.utils.file_utils.FileUtils.get_last_directory')
    @patch('ddprimer.utils.file_utils.FileUtils.save_last_directory')
    @patch('ddprimer.utils.file_utils.FileUtils.use_cli', False)
    @patch('tkinter.filedialog.askopenfilename')
    @patch('tkinter.Tk')
    def test_get_file_gui(self, mock_tk, mock_askopenfilename, mock_save_dir, mock_get_dir):
        """Test getting a file with GUI."""
        # Setup
        mock_get_dir.return_value = "/last/directory"
        mock_askopenfilename.return_value = "/path/to/selected/file.txt"
        
        # Test
        result = FileUtils.get_file("Select file", [("Text Files", "*.txt")])
        
        # Verify
        assert result == "/path/to/selected/file.txt"
        mock_askopenfilename.assert_called_once()
        mock_save_dir.assert_called_once_with(os.path.dirname("/path/to/selected/file.txt"))

    @patch('ddprimer.utils.file_utils.FileUtils.use_cli', True)
    @patch('builtins.input')
    def test_get_file_cli(self, mock_input):
        """Test getting a file with CLI."""
        # Setup
        mock_input.return_value = "/path/to/cli/file.txt"
        
        # Test
        result = FileUtils.get_file("Select file", [("Text Files", "*.txt")])
        
        # Verify
        assert result == "/path/to/cli/file.txt"
        mock_input.assert_called_once()

    @patch('ddprimer.utils.file_utils.FileUtils.use_cli', True)
    @patch('builtins.input')
    def test_get_files_cli(self, mock_input):
        """Test getting multiple files with CLI."""
        # Setup
        mock_input.return_value = "/path/file1.txt /path/file2.txt"
        
        # Test
        result = FileUtils.get_files("Select files", [("Text Files", "*.txt")])
        
        # Verify
        assert result == ["/path/file1.txt", "/path/file2.txt"]
        mock_input.assert_called_once()

    @patch('builtins.open', new_callable=mock_open, read_data=">chr1\nATGCGGCCACTGGCTTAA\n>chr2\nCTAGCCGGTAATCG\n")
    def test_load_fasta(self, mock_file):
        """Test loading a FASTA file."""
        # Test
        sequences = FileUtils.load_fasta("/path/to/fasta.fa")
        
        # Verify
        assert isinstance(sequences, dict)
        assert 'chr1' in sequences
        assert 'chr2' in sequences
        assert sequences['chr1'] == 'ATGCGGCCACTGGCTTAA'
        assert sequences['chr2'] == 'CTAGCCGGTAATCG'
        mock_file.assert_called_once_with("/path/to/fasta.fa", 'r')

    @patch('pandas.read_csv')
    @patch('pandas.read_excel')
    def test_load_sequences_from_table_csv(self, mock_read_excel, mock_read_csv):
        """Test loading sequences from CSV table."""
        # Setup - create a mock DataFrame for CSV
        mock_df = pd.DataFrame({
            'Sequence Name': ['seq1', 'seq2'],
            'Sequence': ['ATGCGGCC', 'CTAGGCCT']
        })
        mock_read_csv.return_value = mock_df
        
        # Test
        sequences = FileUtils.load_sequences_from_table("/path/to/sequences.csv")
        
        # Verify
        assert isinstance(sequences, dict)
        assert 'seq1' in sequences
        assert 'seq2' in sequences
        assert sequences['seq1'] == 'ATGCGGCC'
        assert sequences['seq2'] == 'CTAGGCCT'
        mock_read_csv.assert_called_once_with("/path/to/sequences.csv")
        mock_read_excel.assert_not_called()

    @patch('pandas.read_csv')
    @patch('pandas.read_excel')
    def test_load_sequences_from_table_excel(self, mock_read_excel, mock_read_csv):
        """Test loading sequences from Excel table."""
        # Setup - create a mock DataFrame for Excel
        mock_df = pd.DataFrame({
            'Sequence Name': ['seq1', 'seq2'],
            'Sequence': ['ATGCGGCC', 'CTAGGCCT']
        })
        mock_read_excel.return_value = mock_df
        
        # Test
        sequences = FileUtils.load_sequences_from_table("/path/to/sequences.xlsx")
        
        # Verify
        assert isinstance(sequences, dict)
        assert 'seq1' in sequences
        assert 'seq2' in sequences
        assert sequences['seq1'] == 'ATGCGGCC'
        assert sequences['seq2'] == 'CTAGGCCT'
        mock_read_excel.assert_called_once_with("/path/to/sequences.xlsx")
        mock_read_csv.assert_not_called()

    @patch('pandas.DataFrame.to_excel')
    def test_save_formatted_excel_simple(self, mock_to_excel):
        """Test saving basic Excel file."""
        # Setup
        data = {
            'Gene': ['Gene1', 'Gene2'],
            'Primer F': ['ATGCATGCAT', 'CGCGCGCGCG'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC']
        }
        df = pd.DataFrame(data)
        mock_logger = MagicMock()
        
        # Patch openpyxl availability
        with patch('ddprimer.utils.file_utils.HAS_OPENPYXL', False):
            # Test
            output_path = FileUtils.save_formatted_excel(df, "/path/to/output.xlsx", logger=mock_logger)
            
            # Verify
            assert output_path == "/path/to/output.xlsx"
            mock_to_excel.assert_called_once()
            mock_logger.warning.assert_called_once()

    @patch('pandas.DataFrame.to_excel')
    @patch('ddprimer.utils.file_utils.HAS_OPENPYXL', True)
    @patch('openpyxl.load_workbook')
    def test_save_formatted_excel_with_openpyxl(self, mock_load_workbook, mock_to_excel):
        """Test saving formatted Excel file with openpyxl."""
        # Setup
        data = {
            'Gene': ['Gene1', 'Gene2'],
            'Primer F': ['ATGCATGCAT', 'CGCGCGCGCG'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC']
        }
        df = pd.DataFrame(data)
        mock_logger = MagicMock()
        
        # Create mock workbook and worksheet
        mock_workbook = MagicMock()
        mock_worksheet = MagicMock()
        mock_worksheet.max_row = 3
        mock_worksheet.max_column = 3
        mock_workbook.active = mock_worksheet
        mock_load_workbook.return_value = mock_workbook
        
        # Mock the cell values to simulate headers
        def mock_cell(row, column):
            cell_mock = MagicMock()
            if row == 2:  # Header row
                if column == 1:
                    cell_mock.value = "Gene"
                elif column == 2:
                    cell_mock.value = "Primer F"
                elif column == 3:
                    cell_mock.value = "Primer R"
            return cell_mock
            
        mock_worksheet.cell.side_effect = mock_cell
        
        # Test
        output_path = FileUtils.save_formatted_excel(df, "/path/to/output.xlsx", logger=mock_logger)
        
        # Verify
        assert output_path == "/path/to/output.xlsx"
        mock_to_excel.assert_called_once()
        mock_load_workbook.assert_called_once()
        mock_workbook.save.assert_called_once()

    @patch('ddprimer.utils.file_utils.FileUtils.save_formatted_excel')
    def test_save_results(self, mock_save_formatted, tmp_path):
        """Test saving results to an Excel file."""
        # Create a test DataFrame
        data = {
            'Gene': ['Gene1', 'Gene2'],
            'Primer F': ['ATGCATGCAT', 'CGCGCGCGCG'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC'],
            'Amplicon': ['ATGCATGCATTTAAGCCAGT', 'CGCGCGCGCGTGGCTAGGCC'],
            'Length': [20, 20]
        }
        df = pd.DataFrame(data)
        
        # Set up the output directory
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Mock logger
        mock_logger = MagicMock()
        mock_save_formatted.return_value = os.path.join(output_dir, "Primers_file.xlsx")
        
        # Save the results
        input_file = "/path/to/input/file.fasta"
        output_path = FileUtils.save_results(df, output_dir, input_file, 'standard', logger=mock_logger)
        
        # Verify the results
        assert output_path is not None
        assert 'Primers_' in output_path
        assert output_path.endswith('.xlsx')
        mock_save_formatted.assert_called_once()

    @patch('tempfile.mkdtemp')
    @patch('os.makedirs')
    def test_setup_output_directories(self, mock_makedirs, mock_mkdtemp, tmp_path):
        """Test setting up output and temporary directories."""
        # Setup
        mock_args = MagicMock()
        mock_args.output = str(tmp_path / "output")
        mock_mkdtemp.return_value = str(tmp_path / "temp_dir")
        
        # Test
        output_dir, temp_dir = FileUtils.setup_output_directories(mock_args)
        
        # Verify
        assert output_dir == str(tmp_path / "output")
        assert temp_dir == str(tmp_path / "temp_dir")
        mock_makedirs.assert_called_once_with(output_dir, exist_ok=True)
        mock_mkdtemp.assert_called_once()

    @patch('shutil.rmtree')
    def test_cleanup_temp_directory(self, mock_rmtree, tmp_path):
        """Test cleaning up temporary directory."""
        # Setup
        temp_dir = str(tmp_path / "temp_dir")
        os.makedirs(temp_dir, exist_ok=True)
        mock_logger = MagicMock()
        
        # Test
        with patch('os.path.exists', return_value=True):
            result = FileUtils.cleanup_temp_directory(temp_dir)
        
        # Verify
        assert result is True
        mock_rmtree.assert_called_once_with(temp_dir)

    @patch('os.path.isdir')
    @patch('builtins.open', new_callable=mock_open)
    def test_save_last_directory(self, mock_file, mock_isdir):
        """Test saving last directory to config file."""
        # Setup
        mock_isdir.return_value = True
        
        # Test
        FileUtils.save_last_directory("/path/to/last/dir")
        
        # Verify
        mock_file.assert_called_once()
        mock_file().write.assert_called_once_with("/path/to/last/dir")

    @patch('os.path.exists')
    @patch('os.path.isdir')
    @patch('builtins.open', new_callable=mock_open, read_data="/saved/directory")
    def test_get_last_directory(self, mock_file, mock_isdir, mock_exists):
        """Test getting last directory from config file."""
        # Setup
        mock_exists.return_value = True
        mock_isdir.return_value = True
        
        # Reset class variable to force loading from file
        FileUtils._last_directory = None
        
        # Test
        result = FileUtils.get_last_directory()
        
        # Verify
        assert result == "/saved/directory"
        mock_file.assert_called_once()
        
        # Test caching
        result2 = FileUtils.get_last_directory()
        
        # Should not open the file again
        mock_file.assert_called_once()
        assert result2 == "/saved/directory"