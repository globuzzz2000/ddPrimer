"""
Tests for standard mode functionality
"""
import os
import pytest
from unittest.mock import patch, MagicMock

from ...modes.standard_mode import run_standard_mode
from ...config import Config


class TestStandardMode:
    """Tests for standard mode workflow."""

    @patch('ddPrimer.utils.FileUtils.check_and_create_output_dir')
    @patch('ddPrimer.core.SNP_masking_processor.SNPMaskingProcessor.process_vcf')
    @patch('ddPrimer.utils.FileUtils.load_fasta')
    @patch('ddPrimer.core.SNP_masking_processor.SNPMaskingProcessor.mask_sequences')
    @patch('ddPrimer.core.annotation_processor.AnnotationProcessor.load_annotations')
    @patch('ddPrimer.modes.common.run_primer_design_workflow')
    def test_standard_mode_workflow(
        self,
        mock_run_workflow,
        mock_load_annotations,
        mock_mask_sequences,
        mock_load_fasta,
        mock_process_vcf,
        mock_create_output_dir,
        create_mini_fasta,
        create_mini_vcf,
        create_mini_gff,
        tmp_path
    ):
        """Test the standard mode workflow."""
        # Create test files
        fasta_path = create_mini_fasta()
        vcf_path = create_mini_vcf()
        gff_path = create_mini_gff()
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure the mocks
        mock_create_output_dir.return_value = output_dir
        mock_process_vcf.return_value = [{'chr': 'chr1', 'pos': 10, 'ref': 'G', 'alt': 'A'}]
        mock_load_fasta.return_value = {'chr1': 'ATGCGGCCACTGGCTTAA', 'chr2': 'CTAGCCGGTAATCG'}
        mock_mask_sequences.return_value = {'chr1': 'ATGCNGCCACTGGCTTAA', 'chr2': 'CTAGCCGGTAATCG'}
        mock_load_annotations.return_value = {
            'chr1': [{'start': 1, 'end': 18, 'id': 'Gene1', 'strand': '+'}],
            'chr2': [{'start': 1, 'end': 14, 'id': 'Gene2', 'strand': '+'}]
        }
        mock_run_workflow.return_value = True
        
        # Create a mock args object
        args = MagicMock()
        args.fasta = fasta_path
        args.vcf = vcf_path
        args.gff = gff_path
        args.output = output_dir
        args.cli = True
        
        # Run the standard mode workflow
        result = run_standard_mode(args)
        
        # Verify the results
        assert result is True
        
        # Check that each step of the workflow was called
        mock_create_output_dir.assert_called_once()
        mock_process_vcf.assert_called_once()
        mock_load_fasta.assert_called_once()
        mock_mask_sequences.assert_called_once()
        mock_load_annotations.assert_called_once()
        mock_run_workflow.assert_called_once()
        
        # Check that run_primer_design_workflow was called with correct parameters
        workflow_args = mock_run_workflow.call_args[0]
        assert isinstance(workflow_args[0], dict)  # masked_sequences
        assert workflow_args[1] == output_dir  # output_dir
        assert workflow_args[2] == fasta_path  # reference_file
        assert workflow_args[3] == 'standard'  # mode
        assert isinstance(workflow_args[4], dict)  # genes
        
    def test_standard_mode_error_handling(self):
        """Test that standard mode handles errors gracefully."""
        # Create a mock args object with invalid inputs
        args = MagicMock()
        args.fasta = "nonexistent.fasta"
        args.vcf = "nonexistent.vcf"
        args.gff = "nonexistent.gff"
        args.output = "/path/to/output"
        args.cli = True
        
        # Run the standard mode workflow
        result = run_standard_mode(args)
        
        # The workflow should handle errors and return False
        assert result is False