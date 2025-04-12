"""
Tests for MAF (cross-species) mode functionality
"""
import os
import pytest
from unittest.mock import patch, MagicMock

from ...modes.alignment_mode import run_alignment_mode
from ...config import Config


class TestMafMode:
    """Tests for MAF (cross-species) mode workflow."""

    @patch('ddPrimer.utils.FileUtils.check_and_create_output_dir')
    @patch('ddPrimer.core.SNP_masking_processor.SNPMaskingProcessor.process_vcf')
    @patch('ddPrimer.utils.FileUtils.load_fasta')
    @patch('ddPrimer.core.annotation_processor.AnnotationProcessor.load_annotations')
    @patch('ddPrimer.cross_species.cross_species_workflow.run_cross_species_workflow')
    @patch('ddPrimer.modes.common.run_primer_design_workflow')
    def test_maf_mode_workflow_with_existing_maf(
        self,
        mock_run_workflow,
        mock_cross_species,
        mock_load_annotations,
        mock_load_fasta,
        mock_process_vcf,
        mock_create_output_dir,
        create_mini_fasta,
        create_mini_vcf,
        create_mini_gff,
        tmp_path
    ):
        """Test the MAF mode workflow with an existing MAF file."""
        # Create test files
        fasta_path = create_mini_fasta()
        vcf_path = create_mini_vcf()
        vcf2_path = create_mini_vcf()
        gff_path = create_mini_gff()
        maf_path = tmp_path / "test_alignment.maf"
        with open(maf_path, 'w') as f:
            f.write("##maf version=1\n")
            f.write("a score=100\n")
            f.write("s ref.chr1  0 18 + 50 ATGCGGCCACTGGCTTAA\n")
            f.write("s qry.chr1 99 18 + 50 ATGCAGCCACTGGCTTAA\n")
            
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure the mocks
        mock_create_output_dir.return_value = output_dir
        mock_process_vcf.side_effect = [
            [{'chr': 'chr1', 'pos': 10, 'ref': 'G', 'alt': 'A'}],  # First VCF
            [{'chr': 'chr1', 'pos': 12, 'ref': 'C', 'alt': 'T'}]   # Second VCF
        ]
        mock_load_fasta.return_value = {'chr1': 'ATGCGGCCACTGGCTTAA'}
        mock_load_annotations.return_value = {
            'chr1': [{'start': 1, 'end': 18, 'id': 'Gene1', 'strand': '+'}]
        }
        
        # Mock cross-species workflow result
        mock_cross_species.return_value = (
            {'chr1': 'ATGCNGCCACTGGCTTAA'},  # masked sequences
            {  # coordinate map
                'chr1': {
                    1: {'qry_src': 'qry.chr1', 'qry_pos': 100},
                    2: {'qry_src': 'qry.chr1', 'qry_pos': 101},
                    3: {'qry_src': 'qry.chr1', 'qry_pos': 102}
                }
            }
        )
        
        mock_run_workflow.return_value = True
        
        # Create a mock args object
        args = MagicMock()
        args.alignment = True
        args.maf_file = str(maf_path)
        args.vcf = vcf_path
        args.second_vcf = vcf2_path
        args.gff = gff_path
        args.output = output_dir
        args.min_identity = 80.0
        args.min_length = 20
        args.cli = True
        
        # Run the MAF mode workflow
        result = run_alignment_mode(args)
        
        # Verify the results
        assert result is True
        
        # Check that each step of the workflow was called
        mock_create_output_dir.assert_called_once()
        assert mock_load_fasta.call_count == 2  # Should be called twice (one for each FASTA)
        mock_load_annotations.assert_called_once()
        mock_run_lastz.assert_called_once()
        mock_cross_species.assert_called_once()
        mock_run_workflow.assert_called_once()
        
        # Check that run_primer_design_workflow was called with correct parameters
        workflow_args = mock_run_workflow.call_args[0]
        assert isinstance(workflow_args[0], dict)  # masked_sequences
        assert workflow_args[1] == output_dir  # output_dir
        assert workflow_args[2] == str(maf_path)  # reference_file
        assert workflow_args[3] == 'maf'  # mode
        assert isinstance(workflow_args[4], dict)  # genes
        assert isinstance(workflow_args[5], dict)  # coordinate_map
    
    def test_maf_mode_error_handling(self):
        """Test that MAF mode handles errors gracefully."""
        # Create a mock args object with invalid inputs
        args = MagicMock()
        args.alignment = True
        args.maf_file = "nonexistent.maf"
        args.fasta = None
        args.second_fasta = None
        args.vcf = "nonexistent.vcf"
        args.second_vcf = "nonexistent.vcf"
        args.gff = "nonexistent.gff"
        args.output = "/path/to/output"
        args.min_identity = 80.0
        args.min_length = 20
        args.cli = True
        
        # Run the MAF mode workflow
        result = run_alignment_mode(args)
        
        # The workflow should handle errors and return False
        assert result is False_process_vcf.call_count == 2  # Should be called twice (one for each VCF)
        mock_load_annotations.assert_called_once()
        mock_alignment.assert_called_once()
        mock_run_workflow.assert_called_once()
        
        # Check that run_primer_design_workflow was called with correct parameters
        workflow_args = mock_run_workflow.call_args[0]
        assert isinstance(workflow_args[0], dict)  # masked_sequences
        assert workflow_args[1] == output_dir  # output_dir
        assert workflow_args[2] == str(maf_path)  # reference_file
        assert workflow_args[3] == 'maf'  # mode
        assert isinstance(workflow_args[4], dict)  # genes
        assert isinstance(workflow_args[5], dict)  # coordinate_map

    @patch('ddPrimer.utils.FileUtils.check_and_create_output_dir')
    @patch('ddPrimer.core.SNP_masking_processor.SNPMaskingProcessor.process_vcf')
    @patch('ddPrimer.utils.FileUtils.load_fasta')
    @patch('ddPrimer.core.annotation_processor.AnnotationProcessor.load_annotations')
    @patch('ddPrimer.cross_species.lastz_runner.run_lastz')
    @patch('ddPrimer.cross_species.cross_species_workflow.run_cross_species_workflow')
    @patch('ddPrimer.modes.common.run_primer_design_workflow')
    def test_maf_mode_workflow_with_lastz(
        self,
        mock_run_workflow,
        mock_cross_species,
        mock_run_lastz,
        mock_load_annotations,
        mock_load_fasta,
        mock_process_vcf,
        mock_create_output_dir,
        create_mini_fasta,
        create_mini_vcf,
        create_mini_gff,
        tmp_path
    ):
        """Test the MAF mode workflow with automatically generating alignment."""
        # Create test files
        fasta_path = create_mini_fasta()
        fasta2_path = create_mini_fasta()
        vcf_path = create_mini_vcf()
        vcf2_path = create_mini_vcf()
        gff_path = create_mini_gff()
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Configure the mocks
        mock_create_output_dir.return_value = output_dir
        mock_process_vcf.side_effect = [
            [{'chr': 'chr1', 'pos': 10, 'ref': 'G', 'alt': 'A'}],  # First VCF
            [{'chr': 'chr1', 'pos': 12, 'ref': 'C', 'alt': 'T'}]   # Second VCF
        ]
        mock_load_fasta.side_effect = [
            {'chr1': 'ATGCGGCCACTGGCTTAA'},  # First FASTA
            {'chr1': 'ATGCAGCCACTGGCTTAA'}   # Second FASTA
        ]
        mock_load_annotations.return_value = {
            'chr1': [{'start': 1, 'end': 18, 'id': 'Gene1', 'strand': '+'}]
        }
        
        # Mock LastZ run
        maf_path = tmp_path / "alignment.maf"
        with open(maf_path, 'w') as f:
            f.write("##maf version=1\n")
            f.write("a score=100\n")
            f.write("s ref.chr1  0 18 + 50 ATGCGGCCACTGGCTTAA\n")
            f.write("s qry.chr1 99 18 + 50 ATGCAGCCACTGGCTTAA\n")
        mock_run_lastz.return_value = str(maf_path)
        
        # Mock cross-species workflow result
        mock_cross_species.return_value = (
            {'chr1': 'ATGCNGCCACTGGCTTAA'},  # masked sequences
            {  # coordinate map
                'chr1': {
                    1: {'qry_src': 'qry.chr1', 'qry_pos': 100},
                    2: {'qry_src': 'qry.chr1', 'qry_pos': 101},
                    3: {'qry_src': 'qry.chr1', 'qry_pos': 102}
                }
            }
        )
        
        mock_run_workflow.return_value = True
        
        # Create a mock args object
        args = MagicMock()
        args.alignment = True
        args.maf_file = None
        args.fasta = fasta_path
        args.second_fasta = fasta2_path
        args.vcf = vcf_path
        args.second_vcf = vcf2_path
        args.gff = gff_path
        args.output = output_dir
        args.min_identity = 80.0
        args.min_length = 20
        args.lastz_options = "--format=maf"
        args.cli = True
        
        # Run the MAF mode workflow
        result = run_alignment_mode(args)
        
        # Verify the results
        assert result is True
        
        # Check that each step of the workflow was called
        mock_create_output_dir.assert_called_once()
        assert mock_process_vcf.call_count == 2  # Should be called twice (one for each VCF)
        assert mock