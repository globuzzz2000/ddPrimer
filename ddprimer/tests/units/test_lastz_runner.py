#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the lastz_runner module.

These tests verify that the LastZRunner class functions correctly,
including directory creation, command execution, job creation,
and running parallel alignments.
"""

import os
import pytest
import tempfile
from unittest.mock import patch, MagicMock, mock_open, call

# Import package modules
from ...helpers import LastZRunner
from ...config import AlignmentError


# Fixtures
@pytest.fixture
def lastz_runner():
    """Create a LastZRunner instance for testing."""
    return LastZRunner()


@pytest.fixture
def mock_config():
    """Create a mock Config object."""
    config = MagicMock()
    config.NUM_PROCESSES = 4
    config.LASTZ_OPTIONS = "--format=maf --step=10"
    return config


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


# Test create_directory
def test_create_directory_success(lastz_runner, temp_dir):
    """Test successful directory creation."""
    test_dir = os.path.join(temp_dir, "test_dir")
    
    # Call the function
    result = lastz_runner.create_directory(test_dir)
    
    # Verify the directory was created
    assert os.path.exists(test_dir)
    assert os.path.isdir(test_dir)
    assert result == os.path.abspath(test_dir)
    
    # Test with existing directory
    result2 = lastz_runner.create_directory(test_dir)
    assert result2 == os.path.abspath(test_dir)


@patch('ddprimer.helpers.lastz_runner.os.makedirs')
def test_create_directory_error(mock_makedirs, lastz_runner):
    """Test error handling in directory creation."""
    # Setup mock to raise an exception other than EEXIST
    mock_makedirs.side_effect = OSError(13, "Permission denied")
    
    # Call the function and verify it raises the expected exception
    with pytest.raises(AlignmentError, match="Cannot create directory"):
        lastz_runner.create_directory("/nonexistent/directory")


# Test run_command
@patch('ddprimer.helpers.lastz_runner.subprocess.Popen')
def test_run_command_success(mock_popen, lastz_runner):
    """Test successful command execution."""
    # Setup mock process
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_process.communicate.return_value = (b"Command output", None)
    mock_popen.return_value = mock_process
    
    # Call the function
    result = lastz_runner.run_command("echo test")
    
    # Verify subprocess.Popen was called correctly
    mock_popen.assert_called_once()
    assert mock_popen.call_args[0][0] == "echo test"
    
    # Verify the result
    assert result == 0


@patch('ddprimer.helpers.lastz_runner.subprocess.Popen')
def test_run_command_nonzero_exit(mock_popen, lastz_runner):
    """Test command execution with non-zero exit code."""
    # Setup mock process
    mock_process = MagicMock()
    mock_process.returncode = 1
    mock_process.communicate.return_value = (b"Command failed", None)
    mock_popen.return_value = mock_process
    
    # Call the function
    result = lastz_runner.run_command("invalid command")
    
    # Verify the result
    assert result == 1


@patch('ddprimer.helpers.lastz_runner.subprocess.Popen')
def test_run_command_exception(mock_popen, lastz_runner):
    """Test command execution with exception."""
    # Setup mock to raise an exception
    mock_popen.side_effect = Exception("Command execution failed")
    
    # Call the function
    result = lastz_runner.run_command("echo test")
    
    # Verify the result indicates failure
    assert result == -1


# Test run_commands_parallel
@patch('ddprimer.helpers.lastz_runner.Pool')
def test_run_commands_parallel_success(mock_pool, lastz_runner):
    """Test successful parallel command execution."""
    # Setup mock pool
    mock_pool_instance = MagicMock()
    mock_pool_instance.map.return_value = [0, 0, 0]
    mock_pool.return_value.__enter__.return_value = mock_pool_instance
    
    # Call the function
    commands = ["cmd1", "cmd2", "cmd3"]
    lastz_runner.run_commands_parallel(commands, processes=2)
    
    # Verify Pool was created with the correct number of processes
    mock_pool.assert_called_once_with(processes=2)
    
    # Verify the commands were executed
    mock_pool_instance.map.assert_called_once()
    assert mock_pool_instance.map.call_args[0][0] == lastz_runner.run_command
    assert list(mock_pool_instance.map.call_args[0][1]) == commands


@patch('ddprimer.helpers.lastz_runner.Pool')
def test_run_commands_parallel_with_failures(mock_pool, lastz_runner):
    """Test parallel command execution with some failures."""
    # Setup mock pool
    mock_pool_instance = MagicMock()
    mock_pool_instance.map.return_value = [0, 1, 0]  # One command failed
    mock_pool.return_value.__enter__.return_value = mock_pool_instance
    
    # Call the function
    commands = ["cmd1", "cmd2", "cmd3"]
    lastz_runner.run_commands_parallel(commands)
    
    # Verify execution completes despite failures
    mock_pool_instance.map.assert_called_once()


@patch('ddprimer.helpers.lastz_runner.Pool')
def test_run_commands_parallel_empty_commands(mock_pool, lastz_runner):
    """Test parallel command execution with empty command list."""
    # Call the function with empty commands
    lastz_runner.run_commands_parallel([])
    
    # Verify Pool was not created
    mock_pool.assert_not_called()


@patch('ddprimer.helpers.lastz_runner.Pool')
def test_run_commands_parallel_exception(mock_pool, lastz_runner):
    """Test error handling in parallel command execution."""
    # Setup mock to raise an exception
    mock_pool.side_effect = Exception("Parallel execution failed")
    
    # Call the function and verify it raises the expected exception
    with pytest.raises(AlignmentError, match="Parallel command execution failed"):
        lastz_runner.run_commands_parallel(["cmd1", "cmd2"])


# Test create_alignment_jobs
@patch('ddprimer.helpers.lastz_runner.SeqIO.parse')
def test_create_alignment_jobs_success(mock_seqio_parse, lastz_runner, temp_dir):
    """Test successful creation of alignment jobs."""
    # Setup mock sequences
    mock_seqio_parse.side_effect = [
        [MagicMock(id="chr1"), MagicMock(id="chr2")],  # Query sequences
        [MagicMock(id="scaffold1"), MagicMock(id="scaffold2")]  # Target sequences
    ]
    
    # Call the function
    extract_jobs, align_jobs = lastz_runner.create_alignment_jobs(
        "query.fasta", "target.fasta", temp_dir, "--format=maf"
    )
    
    # Verify the extraction jobs
    assert len(extract_jobs) == 4  # 2 query + 2 target sequences
    for job in extract_jobs:
        assert job.startswith("samtools faidx")
    
    # Verify the alignment jobs
    assert len(align_jobs) == 4  # 2x2 pairs
    for job in align_jobs:
        assert job.startswith("lastz")
        assert "--format=maf" in job


@patch('ddprimer.helpers.lastz_runner.SeqIO.parse')
def test_create_alignment_jobs_error(mock_seqio_parse, lastz_runner, temp_dir):
    """Test error handling in alignment job creation."""
    # Setup mock to raise an exception
    mock_seqio_parse.side_effect = Exception("Failed to parse FASTA")
    
    # Call the function and verify it raises the expected exception
    with pytest.raises(AlignmentError, match="Failed to create alignment jobs"):
        lastz_runner.create_alignment_jobs("query.fasta", "target.fasta", temp_dir)


# Test run_parallel_alignment
@patch('ddprimer.helpers.lastz_runner.LastZRunner.create_directory')
@patch('ddprimer.helpers.lastz_runner.LastZRunner.create_alignment_jobs')
@patch('ddprimer.helpers.lastz_runner.LastZRunner.run_commands_parallel')
@patch('ddprimer.helpers.lastz_runner.glob.glob')
@patch('ddprimer.helpers.lastz_runner.os.path.abspath')
@patch('ddprimer.helpers.lastz_runner.os.path.basename')
@patch('ddprimer.helpers.lastz_runner.os.path.splitext')
def test_run_parallel_alignment_success(
    mock_splitext, mock_basename, mock_abspath, mock_glob,
    mock_run_commands, mock_create_jobs, mock_create_dir,
    lastz_runner, temp_dir
):
    """Test successful parallel alignment execution."""
    # The core issue is in mocking the os.path functions
    # For one function call with "ref.fa", the basename should return "ref.fa" 
    # For another call with "qry.fa", the basename should return "qry.fa"
    
    # Set up basename mock to return appropriate values
    def basename_side_effect(path):
        if "ref.fa" in str(path):
            return "ref.fa"
        elif "qry.fa" in str(path):
            return "qry.fa"
        else:
            return "unknown.fa"

    mock_basename.side_effect = basename_side_effect
    
    # Set up abspath mock to return fully qualified paths
    def abspath_side_effect(path):
        if isinstance(path, str):
            if "Alignments" in path:
                return "/abs/path/Alignments"
            elif "temp" in path:
                return temp_dir
            elif "ref.fa" in path:
                return "/abs/path/ref.fa"
            elif "qry.fa" in path:
                return "/abs/path/qry.fa"
        return "/abs/path/unknown"
    
    mock_abspath.side_effect = abspath_side_effect
    
    # Set up splitext mock to return the correct base name and extension
    def splitext_side_effect(path):
        basename = basename_side_effect(path)
        if basename == "ref.fa":
            return ("ref", ".fa")
        elif basename == "qry.fa":
            return ("qry", ".fa")
        else:
            return ("unknown", ".fa")
    
    mock_splitext.side_effect = splitext_side_effect
    
    # Set up output directories
    mock_create_dir.side_effect = [
        "/abs/path/Alignments",
        temp_dir
    ]
    
    # Set up job creation
    mock_create_jobs.return_value = (
        ["samtools faidx ref.fa chr1", "samtools faidx qry.fa scaffold1"],
        ["lastz chr1.fa scaffold1.fa --format=maf > output.TMP"]
    )
    
    # Set up glob return
    mock_glob.return_value = [f"{temp_dir}/output.TMP"]
    
    # Mock open for file writing
    mock_open_obj = mock_open()
    with patch('builtins.open', mock_open_obj):
        # Call the function
        result = lastz_runner.run_parallel_alignment(
            "ref.fa", "qry.fa", temp_dir, "--format=maf", processes=2, keep_temp=False
        )
    
    # Verify directory creation
    assert mock_create_dir.call_count == 2
    
    # Verify job creation
    mock_create_jobs.assert_called_once()
    
    # Verify commands were executed
    assert mock_run_commands.call_count == 2
    
    # Check if the output file path is correct
    assert result == "/abs/path/Alignments/ref_vs_qry.maf"
    
    # Verify that open was called with the correct output path
    mock_open_obj.assert_any_call("/abs/path/Alignments/ref_vs_qry.maf", 'w')
@patch('ddprimer.helpers.lastz_runner.LastZRunner.create_directory')
@patch('ddprimer.helpers.lastz_runner.LastZRunner.create_alignment_jobs')
def test_run_parallel_alignment_error(
    mock_create_jobs, mock_create_dir, lastz_runner, temp_dir
):
    """Test error handling in parallel alignment execution."""
    # Setup mock to raise an exception
    mock_create_jobs.side_effect = Exception("Job creation failed")
    
    # Call the function and verify it raises the expected exception
    with pytest.raises(AlignmentError, match="LastZ alignment failed"):
        lastz_runner.run_parallel_alignment("ref.fa", "qry.fa", temp_dir)


# Test _clean_up_temp_files
@patch('ddprimer.helpers.lastz_runner.os.path.exists')
@patch('ddprimer.helpers.lastz_runner.shutil.rmtree')
@patch('ddprimer.helpers.lastz_runner.os.remove')
def test_clean_up_temp_files(mock_remove, mock_rmtree, mock_exists, lastz_runner):
    """Test temporary file cleanup."""
    # Setup mock
    mock_exists.return_value = True
    
    # Call the function
    temp_dir = "/tmp/test_dir"
    fai_files = ["/path/to/ref.fa.fai", "/path/to/qry.fa.fai"]
    lastz_runner._clean_up_temp_files(temp_dir, fai_files)
    
    # Verify temporary directory was removed
    mock_rmtree.assert_called_once_with(temp_dir)
    
    # Verify .fai files were removed
    assert mock_remove.call_count == 2
    mock_remove.assert_has_calls([
        call("/path/to/ref.fa.fai"),
        call("/path/to/qry.fa.fai")
    ])


@patch('ddprimer.helpers.lastz_runner.os.path.exists')
@patch('ddprimer.helpers.lastz_runner.shutil.rmtree')
@patch('ddprimer.helpers.lastz_runner.os.remove')
def test_clean_up_temp_files_with_errors(mock_remove, mock_rmtree, mock_exists, lastz_runner):
    """Test cleanup handling errors gracefully."""
    # Setup mocks to raise exceptions
    mock_exists.return_value = True
    mock_rmtree.side_effect = Exception("rmtree failed")
    mock_remove.side_effect = Exception("remove failed")
    
    # Need to patch the logger to ensure no errors are actually logged
    with patch.object(lastz_runner, 'logger') as mock_logger:
        # Call the function
        temp_dir = "/tmp/test_dir"
        fai_files = ["/path/to/ref.fa.fai"]
        
        # Should not raise exceptions despite cleanup errors
        lastz_runner._clean_up_temp_files(temp_dir, fai_files)
    
    # Verify attempts were made
    mock_rmtree.assert_called_once()
    mock_remove.assert_called_once()
    # Verify warnings were logged
    mock_logger.warning.assert_any_call("Could not remove temporary directory /tmp/test_dir: rmtree failed")
    mock_logger.warning.assert_any_call("Could not remove index file /path/to/ref.fa.fai: remove failed")


if __name__ == "__main__":
    pytest.main(["-v"])