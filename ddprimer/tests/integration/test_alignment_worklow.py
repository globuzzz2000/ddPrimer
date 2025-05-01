#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the alignment_workflow module.

These tests verify that the alignment workflow functions correctly,
including handling MAF files, processing alignments, extracting variants,
loading reference sequences, and creating masked sequences.
"""

import os
import pytest
import tempfile
import shutil
import glob
from unittest.mock import patch, MagicMock, mock_open

from ...helpers.alignment_workflow import (
    run_alignment_workflow,
    _handle_maf_file,
    _run_lastz_alignment,
    _process_maf_file,
    _extract_variant_positions,
    _load_reference_sequences,
    _create_masked_sequences,
    _apply_snp_masking,
    _map_second_variants_to_reference,
    _clean_up_intermediate_files,
)
from ...config import Config, AlignmentError, FileFormatError


# Define base test directories for better control
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEMP_BASE = os.path.join(TEST_DIR, 'temp_test_alignment')


# Session-wide fixtures for cleanup
@pytest.fixture(scope="session", autouse=True)
def cleanup_test_dirs():
    """Clean up any test directories before and after running tests."""
    # Remove any leftover test directories before starting
    if os.path.exists(TEMP_BASE):
        shutil.rmtree(TEMP_BASE)
    os.makedirs(TEMP_BASE, exist_ok=True)
    
    # Run the tests
    yield
    
    # Clean up after tests
    if os.path.exists(TEMP_BASE):
        shutil.rmtree(TEMP_BASE)


# Fixtures for common test data
@pytest.fixture
def mock_args():
    """Create a mock args object with necessary attributes."""
    args = MagicMock()
    args.fasta = "test_reference.fasta"
    args.second_fasta = "test_second.fasta"
    args.vcf = "test_reference.vcf"
    args.second_vcf = "test_second.vcf"
    args.maf = None
    args.snp = True
    args.debug = False  # Add debug flag (controls temp file cleanup)
    return args


@pytest.fixture
def mock_output_dir():
    """Create a temporary directory for output files."""
    temp_dir = tempfile.mkdtemp(dir=TEMP_BASE)
    yield temp_dir
    # No need to clean up here, will be handled by the session fixture


@pytest.fixture
def mock_maf_parser():
    """Create a mock MAF parser."""
    mock_parser = MagicMock()
    
    # Setup mock return values
    mock_parser.analyze_maf_file.return_value = {
        'alignment_count': 10,
        'ref_seq_ids': ['chr1', 'chr2'],
        'query_seq_ids': ['qry1', 'qry2']
    }
    
    mock_parser.parse_maf_file.return_value = {'chr1': [{'ref_start': 100, 'ref_end': 200}]}
    
    mock_parser.identify_conserved_regions.return_value = {
        'chr1': [
            {
                'start': 100, 
                'end': 200, 
                'identity': 90.0,
                'qry_src': 'qry1',
                'qry_start': 300,
                'qry_end': 400,
                'qry_strand': '+'
            }
        ]
    }
    
    mock_parser.generate_coordinate_map.return_value = {
        'chr1': {
            150: {'qry_src': 'qry1', 'qry_pos': 350, 'qry_strand': '+'}
        }
    }
    
    mock_parser.extract_reference_sequences_from_maf.return_value = {
        'chr1': 'ATGCATGCATGCATGCATGC'
    }
    
    return mock_parser


@pytest.fixture
def mock_snp_processor():
    """Create a mock SNP masking processor."""
    mock_processor = MagicMock()
    
    # Setup mock return values
    mock_processor.extract_variants_by_regions.return_value = {
        'chr1': {150, 155, 160}
    }
    
    mock_processor.mask_sequences_for_primer_design.return_value = {
        'chr1': 'ATGCATGCNNNCATGCATGC'
    }
    
    mock_processor.combine_masked_sequences.return_value = {
        'chr1': 'ATGCATGCNNNCATGCATGC'
    }
    
    return mock_processor


@pytest.fixture
def mock_lastz_runner():
    """Create a mock LastZ runner."""
    mock_runner = MagicMock()
    mock_runner.run_parallel_alignment.return_value = "test_output.maf"
    return mock_runner


# Test run_alignment_workflow
@patch('ddprimer.helpers.alignment_workflow.MAFParser')
@patch('ddprimer.helpers.alignment_workflow.SNPMaskingProcessor')
@patch('ddprimer.helpers.alignment_workflow._handle_maf_file')
@patch('ddprimer.helpers.alignment_workflow._process_maf_file')
@patch('ddprimer.helpers.alignment_workflow._extract_variant_positions')
@patch('ddprimer.helpers.alignment_workflow._load_reference_sequences')
@patch('ddprimer.helpers.alignment_workflow._create_masked_sequences')
@patch('ddprimer.helpers.alignment_workflow._clean_up_intermediate_files')
def test_run_alignment_workflow_success(
    mock_clean_up, mock_create_masked, mock_load_ref, mock_extract_variants,
    mock_process_maf, mock_handle_maf, mock_snp_processor_class, mock_maf_parser_class,
    mock_args, mock_output_dir
):
    """Test successful execution of run_alignment_workflow."""
    # Setup mocks
    mock_maf_parser_instance = mock_maf_parser_class.return_value
    mock_snp_processor_instance = mock_snp_processor_class.return_value
    
    # Define return values
    mock_handle_maf.return_value = ("test.maf", True, True)
    mock_process_maf.return_value = ({"chr1": []}, {"chr1": {}})
    mock_extract_variants.return_value = ({"chr1": set()}, {"chr1": set()})
    mock_load_ref.return_value = ({"chr1": "ATGC"}, None)
    mock_create_masked.return_value = {"chr1": "ATGC"}
    
    # Run the function
    result = run_alignment_workflow(mock_args, mock_output_dir)
    
    # Verify function calls
    mock_handle_maf.assert_called_once_with(mock_args, mock_output_dir)
    mock_process_maf.assert_called_once()
    mock_extract_variants.assert_called_once()
    mock_load_ref.assert_called_once()
    mock_create_masked.assert_called_once()
    mock_clean_up.assert_called_once()
    
    # Check result
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert result[0] == {"chr1": "ATGC"}  # masked_sequences
    assert result[1] == {"chr1": {}}      # coordinate_map


@patch('ddprimer.helpers.alignment_workflow.MAFParser')
@patch('ddprimer.helpers.alignment_workflow.SNPMaskingProcessor')
def test_run_alignment_workflow_error(
    mock_snp_processor_class, mock_maf_parser_class,
    mock_args, mock_output_dir
):
    """Test error handling in run_alignment_workflow."""
    # Setup maf_parser to raise an exception
    mock_maf_parser_instance = mock_maf_parser_class.return_value
    mock_maf_parser_instance.parse_maf_file.side_effect = Exception("MAF parsing error")
    
    # Verify the exception is wrapped and re-raised
    with pytest.raises(AlignmentError, match="Alignment workflow failed"):
        run_alignment_workflow(mock_args, mock_output_dir)


# Test _handle_maf_file
def test_handle_maf_file_with_maf(mock_args, mock_output_dir):
    """Test _handle_maf_file when a MAF file is provided."""
    # Set up mock_args with a MAF file
    mock_args.maf = "test.maf"
    
    # Run the function
    maf_file, ref_fasta_required, second_fasta_required = _handle_maf_file(mock_args, mock_output_dir)
    
    # Verify results
    assert maf_file == "test.maf"
    assert ref_fasta_required is False
    assert second_fasta_required is False


@patch('ddprimer.helpers.alignment_workflow._run_lastz_alignment')
def test_handle_maf_file_without_maf(mock_run_lastz, mock_args, mock_output_dir):
    """Test _handle_maf_file when no MAF file is provided."""
    # Set up mock_args without a MAF file
    mock_args.maf = None
    mock_run_lastz.return_value = "generated.maf"
    
    # Run the function
    maf_file, ref_fasta_required, second_fasta_required = _handle_maf_file(mock_args, mock_output_dir)
    
    # Verify results
    assert maf_file == "generated.maf"
    assert ref_fasta_required is True
    assert second_fasta_required is True
    mock_run_lastz.assert_called_once_with(mock_args, mock_output_dir)


def test_handle_maf_file_missing_fasta(mock_args, mock_output_dir):
    """Test _handle_maf_file when FASTA files are missing."""
    # Set up mock_args without a MAF file and missing FASTA files
    mock_args.maf = None
    mock_args.fasta = None
    
    # Verify the exception is raised
    with pytest.raises(AlignmentError, match="Reference genome FASTA .* is required"):
        _handle_maf_file(mock_args, mock_output_dir)
    
    # Set reference FASTA but not second FASTA
    mock_args.fasta = "test.fasta"
    mock_args.second_fasta = None
    
    # Verify the exception is raised
    with pytest.raises(AlignmentError, match="Second FASTA .* is required"):
        _handle_maf_file(mock_args, mock_output_dir)


# Test _run_lastz_alignment
@patch('ddprimer.helpers.alignment_workflow.os.makedirs')
@patch('ddprimer.helpers.alignment_workflow.os.path.dirname')
@patch('ddprimer.helpers.alignment_workflow.os.path.abspath')
@patch('ddprimer.helpers.alignment_workflow.LastZRunner')
def test_run_lastz_alignment_success(
    mock_lastz_runner_class, mock_abspath, mock_dirname, mock_makedirs,
    mock_args, mock_output_dir
):
    """Test successful execution of _run_lastz_alignment."""
    # Set up known problematic directory path for targeted cleanup
    problematic_dir = os.path.join(os.getcwd(), "Alignments")
    ddprimer_dir = os.path.join(os.path.expanduser("~"), "ddPrimer", "Alignments")
    
    try:
        # Setup mocks
        mock_abspath.side_effect = lambda path: os.path.join(mock_output_dir, os.path.basename(path))
        mock_dirname.return_value = mock_output_dir
        
        # Create a mock runner instance
        mock_lastz_runner_instance = mock_lastz_runner_class.return_value
        test_maf_file = os.path.join(mock_output_dir, "output.maf") 
        mock_lastz_runner_instance.run_parallel_alignment.return_value = test_maf_file
        
        # Run the function
        result = _run_lastz_alignment(mock_args, mock_output_dir)
        
        # Verify function calls - all paths should be within mock_output_dir
        alignments_dir = os.path.join(mock_output_dir, "Alignments")
        mock_makedirs.assert_called_once_with(alignments_dir, exist_ok=True)
        
        # Verify the lastz runner was called
        mock_lastz_runner_instance.run_parallel_alignment.assert_called_once()
        
        # Check result
        assert result == test_maf_file
        
    finally:
        # Targeted cleanup of only the specific problematic directories
        import shutil
        
        # Clean up the specific directory we know is problematic
        if os.path.exists(problematic_dir):
            shutil.rmtree(problematic_dir)
            print(f"Cleaned up orphaned directory: {problematic_dir}")
            
        # Also check for the user's home ddPrimer directory
        if os.path.exists(ddprimer_dir):
            shutil.rmtree(ddprimer_dir)
            print(f"Cleaned up orphaned directory: {ddprimer_dir}")

@patch('ddprimer.helpers.alignment_workflow.os.makedirs')
@patch('ddprimer.helpers.alignment_workflow.os.path.dirname')
@patch('ddprimer.helpers.alignment_workflow.os.path.abspath')
@patch('ddprimer.helpers.alignment_workflow.LastZRunner')
def test_run_lastz_alignment_error(
    mock_lastz_runner_class, mock_abspath, mock_dirname, mock_makedirs,
    mock_args, mock_output_dir
):
    """Test error handling in _run_lastz_alignment."""
    # Setup mocks
    mock_abspath.side_effect = lambda path: os.path.join(mock_output_dir, os.path.basename(path))
    mock_dirname.return_value = mock_output_dir
    mock_lastz_runner_instance = mock_lastz_runner_class.return_value
    mock_lastz_runner_instance.run_parallel_alignment.side_effect = Exception("LastZ failed")
    
    # Verify the exception is wrapped and re-raised
    with pytest.raises(AlignmentError, match="LastZ alignment failed"):
        _run_lastz_alignment(mock_args, mock_output_dir)


# Test _process_maf_file
def test_process_maf_file_success(mock_maf_parser, mock_args):
    """Test successful execution of _process_maf_file."""
    # Run the function
    conserved_regions, coordinate_map = _process_maf_file("test.maf", mock_maf_parser, mock_args)
    
    # Verify function calls
    mock_maf_parser.analyze_maf_file.assert_called_once_with("test.maf")
    mock_maf_parser.parse_maf_file.assert_called_once_with("test.maf")
    mock_maf_parser.identify_conserved_regions.assert_called_once()
    mock_maf_parser.generate_coordinate_map.assert_called_once_with(
        mock_maf_parser.identify_conserved_regions.return_value
    )
    
    # Check results
    assert conserved_regions == mock_maf_parser.identify_conserved_regions.return_value
    assert coordinate_map == mock_maf_parser.generate_coordinate_map.return_value


def test_process_maf_file_error(mock_maf_parser, mock_args):
    """Test error handling in _process_maf_file."""
    # Setup mock to raise an exception
    mock_maf_parser.analyze_maf_file.side_effect = Exception("MAF analysis error")
    
    # Verify the exception is wrapped and re-raised
    with pytest.raises(AlignmentError, match="MAF file processing failed"):
        _process_maf_file("test.maf", mock_maf_parser, mock_args)


# Test _extract_variant_positions
def test_extract_variant_positions_success(mock_snp_processor, mock_args):
    """Test successful execution of _extract_variant_positions."""
    # Define test data
    conserved_regions = {
        'chr1': [
            {
                'start': 100, 
                'end': 200,
                'identity': 90.0,
                'qry_src': 'qry1',
                'qry_start': 300,
                'qry_end': 400,
                'qry_strand': '+'
            }
        ]
    }
    coordinate_map = {
        'chr1': {
            150: {'qry_src': 'qry1', 'qry_pos': 350, 'qry_strand': '+'}
        }
    }
    
    # Run the function
    ref_variants, second_variants = _extract_variant_positions(
        mock_args, mock_snp_processor, conserved_regions, coordinate_map
    )
    
    # Verify function calls
    mock_snp_processor.extract_variants_by_regions.assert_called()
    assert mock_snp_processor.extract_variants_by_regions.call_count >= 1
    
    # Check results
    assert ref_variants == mock_snp_processor.extract_variants_by_regions.return_value
    assert isinstance(second_variants, dict)


def test_extract_variant_positions_with_errors(mock_snp_processor, mock_args):
    """Test _extract_variant_positions handling errors gracefully."""
    # Define test data
    conserved_regions = {'chr1': [{'start': 100, 'end': 200, 'qry_src': 'qry1', 'qry_start': 300, 'qry_end': 400, 'qry_strand': '+'}]}
    coordinate_map = {'chr1': {150: {'qry_src': 'qry1', 'qry_pos': 350, 'qry_strand': '+'}}}
    
    # Make the first call succeed but the second fail
    mock_snp_processor.extract_variants_by_regions.side_effect = [
        {'chr1': {150}},  # First call succeeds
        Exception("VCF extraction error")  # Second call fails
    ]
    
    # Run the function
    ref_variants, second_variants = _extract_variant_positions(
        mock_args, mock_snp_processor, conserved_regions, coordinate_map
    )
    
    # The function should continue despite the error
    assert ref_variants == {'chr1': {150}}
    assert second_variants == {}  # Empty because of the error


# Test _load_reference_sequences
@patch('ddprimer.helpers.alignment_workflow.tempfile.mkdtemp')
@patch('ddprimer.helpers.alignment_workflow.FileIO')
def test_load_reference_sequences_from_maf(
    mock_file_io, mock_mkdtemp, mock_maf_parser, mock_args, mock_output_dir
):
    """Test loading reference sequences from MAF file."""
    # Set up mock for temp directory - use controlled path
    temp_dir = os.path.join(TEMP_BASE, "maf_extraction_temp")
    mock_mkdtemp.return_value = temp_dir
    
    # Set up mock for file writing
    mock_open_obj = mock_open()
    with patch('builtins.open', mock_open_obj):
        # Run the function with ref_fasta_required=False
        reference_sequences, temp_dir_result = _load_reference_sequences(
            mock_args, mock_maf_parser, False, mock_output_dir
        )
    
    # Verify that reference sequences were extracted from MAF
    mock_maf_parser.extract_reference_sequences_from_maf.assert_called_once()
    mock_open_obj.assert_called_once()
    
    # Check results
    assert reference_sequences == mock_maf_parser.extract_reference_sequences_from_maf.return_value
    assert temp_dir_result == temp_dir


@patch('ddprimer.helpers.alignment_workflow.FileIO')
def test_load_reference_sequences_from_fasta(
    mock_file_io, mock_maf_parser, mock_args, mock_output_dir
):
    """Test loading reference sequences from FASTA file."""
    # Set up mock for FileIO.load_fasta
    mock_file_io.load_fasta.return_value = {'chr1': 'ATGCATGC'}
    
    # Run the function with ref_fasta_required=True
    reference_sequences, temp_dir = _load_reference_sequences(
        mock_args, mock_maf_parser, True, mock_output_dir
    )
    
    # Verify that sequences were loaded from FASTA
    mock_file_io.load_fasta.assert_called_once_with(mock_args.fasta)
    mock_maf_parser.extract_reference_sequences_from_maf.assert_not_called()
    
    # Check results
    assert reference_sequences == {'chr1': 'ATGCATGC'}
    assert temp_dir is None


@patch('ddprimer.helpers.alignment_workflow.FileIO')
def test_load_reference_sequences_error(
    mock_file_io, mock_maf_parser, mock_args, mock_output_dir
):
    """Test error handling in _load_reference_sequences."""
    # Set up mock to raise an exception
    mock_file_io.load_fasta.side_effect = Exception("FASTA loading error")
    
    # Verify the exception is wrapped and re-raised
    with pytest.raises(FileFormatError, match="Failed to load reference FASTA"):
        _load_reference_sequences(mock_args, mock_maf_parser, True, mock_output_dir)


# Test _create_masked_sequences
@patch('ddprimer.helpers.alignment_workflow.FileIO')
def test_create_masked_sequences_no_snp(
    mock_file_io, mock_snp_processor, mock_maf_parser, mock_args, mock_output_dir
):
    """Test creating masked sequences without SNP masking."""
    # Set up args to disable SNP masking
    mock_args.snp = False
    
    # Set up mock for FileIO.load_fasta
    mock_file_io.load_fasta.return_value = {'chr1': 'ATGCATGC'}
    
    # Set up mock for open
    mock_open_obj = mock_open()
    with patch('builtins.open', mock_open_obj):
        # Run the function
        result = _create_masked_sequences(
            mock_args, mock_snp_processor, mock_maf_parser,
            {'chr1': 'ATGCATGC'}, {'chr1': []}, {'chr1': set()}, {'chr1': set()},
            mock_output_dir, {'chr1': {}}
        )
    
    # Verify MAF parser was called to mask non-conserved regions
    mock_maf_parser.mask_non_conserved_regions.assert_called_once()
    
    # Verify that SNP masking was not applied
    mock_snp_processor.mask_sequences_for_primer_design.assert_not_called()
    
    # Check results
    assert result == mock_file_io.load_fasta.return_value


@patch('ddprimer.helpers.alignment_workflow.FileIO')
@patch('ddprimer.helpers.alignment_workflow._apply_snp_masking')
def test_create_masked_sequences_with_snp(
    mock_apply_snp_masking, mock_file_io, mock_snp_processor, mock_maf_parser, 
    mock_args, mock_output_dir
):
    """Test creating masked sequences with SNP masking."""
    # Set up args to enable SNP masking
    mock_args.snp = True
    
    # Set up mock for FileIO.load_fasta
    mock_file_io.load_fasta.return_value = {'chr1': 'ATGCATGC'}
    
    # Set up mock for _apply_snp_masking
    mock_apply_snp_masking.return_value = {'chr1': 'ATGCNNGC'}
    
    # Set up mock for open
    mock_open_obj = mock_open()
    with patch('builtins.open', mock_open_obj):
        # Run the function
        result = _create_masked_sequences(
            mock_args, mock_snp_processor, mock_maf_parser,
            {'chr1': 'ATGCATGC'}, {'chr1': []}, {'chr1': set()}, {'chr1': set()},
            mock_output_dir, {'chr1': {}}
        )
    
    # Verify MAF parser was called to mask non-conserved regions
    mock_maf_parser.mask_non_conserved_regions.assert_called_once()
    
    # Verify that SNP masking was applied
    mock_apply_snp_masking.assert_called_once()
    
    # Check results
    assert result == {'chr1': 'ATGCNNGC'}


@patch('ddprimer.helpers.alignment_workflow.FileIO')
def test_create_masked_sequences_error(
    mock_file_io, mock_snp_processor, mock_maf_parser, mock_args, mock_output_dir
):
    """Test error handling in _create_masked_sequences."""
    # Set up mock to raise an exception
    mock_maf_parser.mask_non_conserved_regions.side_effect = Exception("Masking error")
    
    # Verify the exception is wrapped and re-raised
    with pytest.raises(AlignmentError, match="Failed to mask non-conserved regions"):
        _create_masked_sequences(
            mock_args, mock_snp_processor, mock_maf_parser,
            {'chr1': 'ATGCATGC'}, {'chr1': []}, {'chr1': set()}, {'chr1': set()},
            mock_output_dir, {'chr1': {}}
        )


# Test _apply_snp_masking
def test_apply_snp_masking_success(mock_snp_processor):
    """Test successful execution of _apply_snp_masking."""
    # Define test data
    alignment_masked_sequences = {'chr1': 'ATGCATGC'}
    ref_variants = {'chr1': {2, 4}}
    second_variants = {'qry1': {352, 354}}
    coordinate_map = {'chr1': {2: {'qry_src': 'qry1', 'qry_pos': 352, 'qry_strand': '+'}}}
    
    # Run the function
    result = _apply_snp_masking(
        mock_snp_processor, alignment_masked_sequences, ref_variants, second_variants, coordinate_map
    )
    
    # Verify function calls
    assert mock_snp_processor.mask_sequences_for_primer_design.call_count >= 1
    mock_snp_processor.combine_masked_sequences.assert_called_once()
    
    # Check results
    assert result == mock_snp_processor.combine_masked_sequences.return_value


def test_apply_snp_masking_error(mock_snp_processor):
    """Test error handling in _apply_snp_masking."""
    # Set up mock to raise an exception
    mock_snp_processor.mask_sequences_for_primer_design.side_effect = Exception("SNP masking error")
    
    # Verify the exception is wrapped and re-raised
    with pytest.raises(AlignmentError, match="Failed to apply SNP masking"):
        _apply_snp_masking(
            mock_snp_processor, {'chr1': 'ATGCATGC'}, {'chr1': {2}}, {}, {}
        )


# Test _map_second_variants_to_reference
def test_map_second_variants_to_reference():
    """Test mapping second genome variants to reference coordinates."""
    # Define test data
    second_variants = {'qry1': {350, 355, 360}}
    coordinate_map = {
        'chr1': {
            100: {'qry_src': 'qry1', 'qry_pos': 350, 'qry_strand': '+'},
            105: {'qry_src': 'qry1', 'qry_pos': 355, 'qry_strand': '+'},
            110: {'qry_src': 'qry1', 'qry_pos': 360, 'qry_strand': '+'},
        }
    }
    
    # Run the function
    result = _map_second_variants_to_reference(second_variants, coordinate_map)
    
    # Check results
    assert 'chr1' in result
    assert result['chr1'] == {100, 105, 110}


# Test _clean_up_intermediate_files
@patch('ddprimer.helpers.alignment_workflow.os.path.exists')
@patch('ddprimer.helpers.alignment_workflow.os.remove')
@patch('ddprimer.helpers.alignment_workflow.shutil.rmtree')
def test_clean_up_intermediate_files(
    mock_rmtree, mock_remove, mock_exists, mock_args, mock_output_dir
):
    """Test cleaning up intermediate files."""
    # Set up mocks
    mock_exists.return_value = True
    
    # First test normal mode (not debug)
    mock_args.debug = False
    
    # Run the function with a temporary directory
    temp_dir = os.path.join(TEMP_BASE, "custom_temp_dir")
    _clean_up_intermediate_files(mock_args, mock_output_dir, temp_dir)
    
    # Verify function calls
    mock_rmtree.assert_called_once_with(temp_dir)
    assert mock_remove.call_count == 2  # Two files to remove
    
    # Reset all mocks to prepare for debug mode test
    mock_rmtree.reset_mock()
    mock_remove.reset_mock()
    
    # Test debug mode (should keep files)
    # In alignment_workflow.py, debug mode check needs to be fixed!
    mock_args.debug = True
    Config.DEBUG_MODE = True  # Also need to set Config.DEBUG_MODE
    
    _clean_up_intermediate_files(mock_args, mock_output_dir, temp_dir)
    
    # Verify nothing was removed in debug mode
    mock_rmtree.assert_not_called()
    mock_remove.assert_not_called()
    
    # Reset Config.DEBUG_MODE for other tests
    Config.DEBUG_MODE = False


@pytest.fixture(scope="session", autouse=True)
def check_cleanup():
    """Verify no orphaned test directories exist after tests."""
    yield
    
    # Check for orphaned ddPrimer directories
    user_home = os.path.expanduser("~")
    ddprimer_dirs = glob.glob(os.path.join(user_home, "**/ddPrimer/Alignments"), recursive=True)
    temp_dirs = glob.glob(os.path.join(user_home, "**/ddprimer_temp_*"), recursive=True)
    
    # Print warning for orphaned directories
    if ddprimer_dirs:
        print(f"WARNING: Found orphaned ddPrimer directories: {ddprimer_dirs}")
    if temp_dirs:
        print(f"WARNING: Found orphaned temp directories: {temp_dirs}")


if __name__ == "__main__":
    pytest.main(["-v"])