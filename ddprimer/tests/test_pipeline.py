"""
End-to-end tests for the ddPrimer pipeline
"""
import os
import sys
import pytest
from unittest.mock import patch, MagicMock

# Import the pipeline module
from .. import pipeline
from ..config import Config


class TestPipeline:
    """Tests for the main pipeline module."""

    @patch('ddPrimer.modes.run_standard_mode')
    def test_standard_mode_workflow(self, mock_standard_mode, create_mini_fasta, create_mini_vcf, create_mini_gff, tmp_path):
        """Test the pipeline with standard mode."""
        # Create test files
        fasta_path = create_mini_fasta()
        vcf_path = create_mini_vcf()
        gff_path = create_mini_gff()
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure the mock
        mock_standard_mode.return_value = True
        
        # Prepare command line arguments
        test_args = [
            "ddprimer",
            "--fasta", fasta_path,
            "--vcf", vcf_path,
            "--gff", gff_path,
            "--output", output_dir,
            "--cli"  # Force CLI mode for testing
        ]
        
        # Set system args
        with patch.object(sys, 'argv', test_args):
            # Run the pipeline
            result = pipeline.run_pipeline()
            
            # Verify the results
            assert result is True
            
            # Check that the standard mode function was called
            mock_standard_mode.assert_called_once()
            
            # Check that the arguments were passed correctly
            args = mock_standard_mode.call_args[0][0]
            assert args.fasta == fasta_path
            assert args.vcf == vcf_path
            assert args.gff == gff_path
            assert args.output == output_dir
            assert args.cli is True

    @patch('ddPrimer.modes.run_direct_mode')
    def test_direct_mode_workflow(self, mock_direct_mode, create_mini_direct_input, tmp_path):
        """Test the pipeline with direct mode."""
        # Create test files
        direct_input_path = create_mini_direct_input()
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure the mock
        mock_direct_mode.return_value = True
        
        # Prepare command line arguments
        test_args = [
            "ddprimer",
            "--direct", direct_input_path,
            "--output", output_dir,
            "--cli"  # Force CLI mode for testing
        ]
        
        # Set system args
        with patch.object(sys, 'argv', test_args):
            # Run the pipeline
            result = pipeline.run_pipeline()
            
            # Verify the results
            assert result is True
            
            # Check that the direct mode function was called
            mock_direct_mode.assert_called_once()
            
            # Check that the arguments were passed correctly
            args = mock_direct_mode.call_args[0][0]
            assert args.direct == direct_input_path
            assert args.output == output_dir
            assert args.cli is True

    @patch('ddPrimer.modes.run_alignment_mode')
    def test_alignment_mode_workflow(
        self, mock_alignment_mode, create_mini_fasta, create_mini_vcf, create_mini_gff, tmp_path
    ):
        """Test the pipeline with alignment (cross-species) mode."""
        # Create test files
        fasta1_path = create_mini_fasta()
        vcf1_path = create_mini_vcf()
        fasta2_path = create_mini_fasta()  # Second species
        vcf2_path = create_mini_vcf()  # Second species VCF
        gff_path = create_mini_gff()
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure the mock
        mock_alignment_mode.return_value = True
        
        # Prepare command line arguments
        test_args = [
            "ddprimer",
            "--alignment",
            "--fasta", fasta1_path,
            "--vcf", vcf1_path,
            "--second-fasta", fasta2_path,
            "--second-vcf", vcf2_path,
            "--gff", gff_path,
            "--output", output_dir,
            "--min-identity", "85",
            "--cli"  # Force CLI mode for testing
        ]
        
        # Set system args
        with patch.object(sys, 'argv', test_args):
            # Run the pipeline
            result = pipeline.run_pipeline()
            
            # Verify the results
            assert result is True
            
            # Check that the alignment mode function was called
            mock_alignment_mode.assert_called_once()
            
            # Check that the arguments were passed correctly
            args = mock_alignment_mode.call_args[0][0]
            assert args.alignment is True
            assert args.fasta == fasta1_path
            assert args.vcf == vcf1_path
            assert args.second_fasta == fasta2_path
            assert args.second_vcf == vcf2_path
            assert args.gff == gff_path
            assert args.output == output_dir
            assert args.min_identity == 85.0
            assert args.cli is True

    def test_pipeline_error_handling(self):
        """Test that the pipeline handles errors gracefully."""
        # Prepare command line arguments with invalid inputs
        test_args = [
            "ddprimer",
            "--fasta", "nonexistent.fasta",
            "--vcf", "nonexistent.vcf",
            "--gff", "nonexistent.gff",
            "--cli"  # Force CLI mode for testing
        ]
        
        # Set system args
        with patch.object(sys, 'argv', test_args):
            # Run the pipeline
            result = pipeline.run_pipeline()
            
            # The pipeline should handle errors and return False
            assert result is False