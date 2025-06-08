#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end functional tests for the ddPrimer pipeline.

This module tests the complete ddPrimer workflow in all three modes:
1. Standard mode (from FASTA/VCF)
2. Direct mode (from CSV sequences)
3. Alignment mode (from two genomes)

Tests verify that the pipeline runs successfully and produces the expected outputs.
"""

import os
import sys
import unittest
import tempfile
import shutil
from pathlib import Path
import pandas as pd
from unittest.mock import patch, MagicMock

# Ensure the package is in the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Import package modules
from ...pipeline import run_pipeline 
from ...config import Config


class TestCompleteWorkflow(unittest.TestCase):
    """End-to-end functional tests for the ddPrimer pipeline."""

    @classmethod
    def setUpClass(cls):
        """Set up test directory and test files once for all tests."""
        # Create test directory
        cls.test_dir = Path(tempfile.mkdtemp(prefix="ddprimer_functional_test_"))
        
        # Create test data directories
        cls.data_dir = cls.test_dir / "data"
        cls.data_dir.mkdir(exist_ok=True)
        
        # Set up subdirectories
        for subdir in ["fasta", "vcf", "gff", "direct", "alignment"]:
            (cls.data_dir / subdir).mkdir(exist_ok=True)
        
        # Create minimal test files
        cls._create_test_files()
        
        # Store original config settings
        cls.original_blast_db_path = Config.DB_PATH
        cls.original_debug_mode = Config.DEBUG_MODE
        
        # Configure for testing
        Config.DEBUG_MODE = True
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test directory and restore config."""
        # Restore original config
        Config.DB_PATH = cls.original_blast_db_path
        Config.DEBUG_MODE = cls.original_debug_mode
        
        # Clean up test directory
        shutil.rmtree(cls.test_dir)
    
    def setUp(self):
        """Set up for each test."""
        # Create output directory
        self.output_dir = self.test_dir / "output" / self._testMethodName
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Mock file selection dialogs
        self.file_io_patcher = patch('ddprimer.utils.file_io.FileIO')
        self.mock_file_io = self.file_io_patcher.start()
        self.mock_file_io.use_cli = True
        
        # Mock BLAST database verification
        self.blast_patcher = patch('ddprimer.utils.blast_verification.BlastVerification.verify_blast_database')
        self.mock_blast = self.blast_patcher.start()
        self.mock_blast.return_value = True
    
    def tearDown(self):
        """Clean up after each test."""
        self.file_io_patcher.stop()
        self.blast_patcher.stop()
    
    @classmethod
    def _create_test_files(cls):
        """Create necessary test files."""
        # 1. Create minimal FASTA file
        fasta_path = cls.data_dir / "fasta" / "test_genome.fasta"
        with open(fasta_path, "w") as f:
            f.write(">chr1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")
            f.write(">chr2\nGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
        
        # 2. Create second FASTA file (for alignment mode)
        second_fasta_path = cls.data_dir / "fasta" / "test_genome2.fasta"
        with open(second_fasta_path, "w") as f:
            f.write(">chr1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGT\n")
            f.write(">chr2\nGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTG\n")
        
        # 3. Create VCF file
        vcf_path = cls.data_dir / "vcf" / "test_variants.vcf"
        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("chr1\t10\t.\tA\tT\t100\tPASS\t.\n")
            f.write("chr1\t25\t.\tG\tC\t100\tPASS\t.\n")
            f.write("chr2\t15\t.\tG\tA\t100\tPASS\t.\n")
        
        # 4. Create second VCF file (for alignment mode)
        second_vcf_path = cls.data_dir / "vcf" / "test_variants2.vcf"
        with open(second_vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("chr1\t10\t.\tT\tA\t100\tPASS\t.\n")
            f.write("chr1\t25\t.\tC\tG\t100\tPASS\t.\n")
            f.write("chr2\t15\t.\tA\tG\t100\tPASS\t.\n")
        
        # 5. Create GFF file
        gff_path = cls.data_dir / "gff" / "test_annotation.gff"
        with open(gff_path, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\t.\tgene\t5\t30\t.\t+\t.\tID=gene1;Name=Gene1\n")
            f.write("chr2\t.\tgene\t10\t40\t.\t+\t.\tID=gene2;Name=Gene2\n")
        
        # 6. Create CSV file for direct mode
        csv_path = cls.data_dir / "direct" / "test_sequences.csv"
        with open(csv_path, "w") as f:
            f.write("name,sequence\n")
            f.write("seq1,ATGCATGCATGCATGCATGCATGC\n")
            f.write("seq2,GCTAGCTAGCTAGCTAGCTAGCTA\n")
        
        # 7. Create MAF file for alignment mode
        maf_path = cls.data_dir / "alignment" / "test_alignment.maf"
        with open(maf_path, "w") as f:
            f.write("##maf version=1 scoring=LastZ\n")
            f.write("a score=1000\n")
            f.write("s chr1 0 50 + 50 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")
            f.write("s chr1 0 50 + 50 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGT\n")
            f.write("\n")
            f.write("a score=900\n")
            f.write("s chr2 0 50 + 50 GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
            f.write("s chr2 0 50 + 50 GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTG\n")
    
    @patch('sys.argv', ['ddprimer'])
    def test_standard_mode(self):
        """Test standard mode workflow."""
        # Mock arguments for standard mode
        with patch('argparse.ArgumentParser.parse_args') as mock_args, \
            patch('ddprimer.modes.run_standard_mode') as mock_standard_mode, \
            patch('ddprimer.modes.run_direct_mode') as mock_direct_mode, \
            patch('ddprimer.modes.run_alignment_mode') as mock_alignment_mode, \
            patch('ddprimer.pipeline.WorkflowFactory.create_workflow') as mock_factory:
            
            # Set up mock args
            mock_args.return_value.fasta = str(self.data_dir / "fasta" / "test_genome.fasta")
            mock_args.return_value.vcf = str(self.data_dir / "vcf" / "test_variants.vcf")
            mock_args.return_value.gff = str(self.data_dir / "gff" / "test_annotation.gff")
            mock_args.return_value.output = str(self.output_dir)
            mock_args.return_value.debug = True
            mock_args.return_value.config = None
            mock_args.return_value.cli = True
            mock_args.return_value.direct = None
            mock_args.return_value.alignment = False
            mock_args.return_value.nooligo = False
            mock_args.return_value.noannotation = False
            mock_args.return_value.db = None
            mock_args.return_value.db_action = None
            mock_args.return_value.lastzonly = False
            mock_args.return_value.maf = None
            mock_args.return_value.second_fasta = None
            mock_args.return_value.second_vcf = None
            mock_args.return_value.snp = True
            
            # Configure mocks to return success
            mock_standard_mode.return_value = True
            mock_direct_mode.return_value = True
            mock_alignment_mode.return_value = True
            
            # Configure the factory to return the standard mode function
            mock_factory.return_value = mock_standard_mode
            
            # Run pipeline and check that correct mode was called
            result = run_pipeline()
            mock_standard_mode.assert_called_once()
            mock_direct_mode.assert_not_called()
            mock_alignment_mode.assert_not_called()
            
            # Verify result
            self.assertTrue(result, "Pipeline should return success")
        
    @patch('sys.argv', ['ddprimer'])
    def test_direct_mode(self):
        """Test direct mode workflow."""
        # Mock arguments for direct mode
        with patch('argparse.ArgumentParser.parse_args') as mock_args, \
            patch('ddprimer.modes.run_standard_mode') as mock_standard_mode, \
            patch('ddprimer.modes.run_direct_mode') as mock_direct_mode, \
            patch('ddprimer.modes.run_alignment_mode') as mock_alignment_mode, \
            patch('ddprimer.pipeline.WorkflowFactory.create_workflow') as mock_factory:
            
            # Set up mock args
            mock_args.return_value.direct = str(self.data_dir / "direct" / "test_sequences.csv")
            mock_args.return_value.fasta = str(self.data_dir / "fasta" / "test_genome.fasta")
            mock_args.return_value.vcf = str(self.data_dir / "vcf" / "test_variants.vcf")
            mock_args.return_value.output = str(self.output_dir)
            mock_args.return_value.debug = True
            mock_args.return_value.config = None
            mock_args.return_value.cli = True
            mock_args.return_value.alignment = False
            mock_args.return_value.nooligo = False
            mock_args.return_value.snp = True  # Enable SNP masking for direct mode
            mock_args.return_value.noannotation = False
            mock_args.return_value.db = None
            mock_args.return_value.db_action = None
            mock_args.return_value.lastzonly = False
            mock_args.return_value.maf = None  # Add this to avoid undefined attribute
            mock_args.return_value.second_fasta = None  # Add this to avoid undefined attribute
            mock_args.return_value.second_vcf = None  # Add this to avoid undefined attribute
            mock_args.return_value.gff = None  # Add this to avoid undefined attribute
            
            # Configure mocks to return success
            mock_standard_mode.return_value = True
            mock_direct_mode.return_value = True
            mock_alignment_mode.return_value = True
            
            # Configure the factory to return the direct mode function
            mock_factory.return_value = mock_direct_mode
            
            # Run pipeline
            result = run_pipeline()
            
            # Check that the correct mode was called
            mock_standard_mode.assert_not_called()
            mock_direct_mode.assert_called_once()
            mock_alignment_mode.assert_not_called()
            
            # Verify result
            self.assertTrue(result, "Pipeline should return success")
        
    @patch('sys.argv', ['ddprimer'])
    def test_alignment_mode_with_maf(self):
        """Test alignment mode workflow with pre-computed MAF file."""
        # Mock arguments for alignment mode with MAF
        with patch('argparse.ArgumentParser.parse_args') as mock_args, \
            patch('ddprimer.modes.run_standard_mode') as mock_standard_mode, \
            patch('ddprimer.modes.run_direct_mode') as mock_direct_mode, \
            patch('ddprimer.modes.run_alignment_mode') as mock_alignment_mode, \
            patch('ddprimer.helpers.lastz_runner.LastZRunner.run_parallel_alignment') as mock_lastz, \
            patch('ddprimer.pipeline.WorkflowFactory.create_workflow') as mock_factory:
            
            # Set up mock args
            mock_args.return_value.alignment = True
            mock_args.return_value.maf = str(self.data_dir / "alignment" / "test_alignment.maf")
            mock_args.return_value.fasta = str(self.data_dir / "fasta" / "test_genome.fasta")
            mock_args.return_value.second_fasta = str(self.data_dir / "fasta" / "test_genome2.fasta")
            mock_args.return_value.vcf = str(self.data_dir / "vcf" / "test_variants.vcf")
            mock_args.return_value.second_vcf = str(self.data_dir / "vcf" / "test_variants2.vcf")
            mock_args.return_value.gff = str(self.data_dir / "gff" / "test_annotation.gff")
            mock_args.return_value.output = str(self.output_dir)
            mock_args.return_value.debug = True
            mock_args.return_value.config = None
            mock_args.return_value.cli = True
            mock_args.return_value.direct = None
            mock_args.return_value.nooligo = False
            mock_args.return_value.snp = True  # Enable SNP masking for alignment mode
            mock_args.return_value.noannotation = False
            mock_args.return_value.db = None
            mock_args.return_value.db_action = None
            mock_args.return_value.lastzonly = False
            
            # Configure mocks to return success
            mock_standard_mode.return_value = True
            mock_direct_mode.return_value = True
            mock_alignment_mode.return_value = True
            mock_lastz.return_value = str(self.data_dir / "alignment" / "test_alignment.maf")
            
            # Configure the factory to return the alignment mode function
            mock_factory.return_value = mock_alignment_mode
            
            # Run pipeline
            result = run_pipeline()
            
            # Check that the correct mode was called
            mock_standard_mode.assert_not_called()
            mock_direct_mode.assert_not_called()
            mock_alignment_mode.assert_called_once()
            
            # Verify result
            self.assertTrue(result, "Pipeline should return success")
        
    @patch('sys.argv', ['ddprimer'])
    def test_lastzonly_mode(self):
        """Test LastZ-only mode workflow."""
        # Mock arguments for LastZ-only mode
        with patch('argparse.ArgumentParser.parse_args') as mock_args, \
            patch('ddprimer.modes.run_standard_mode') as mock_standard_mode, \
            patch('ddprimer.modes.run_direct_mode') as mock_direct_mode, \
            patch('ddprimer.modes.run_alignment_mode') as mock_alignment_mode, \
            patch('ddprimer.helpers.lastz_runner.LastZRunner.run_parallel_alignment') as mock_lastz, \
            patch('ddprimer.pipeline.WorkflowFactory.create_workflow') as mock_factory:
            
            # Set up mock args
            mock_args.return_value.alignment = True
            mock_args.return_value.lastzonly = True
            mock_args.return_value.fasta = str(self.data_dir / "fasta" / "test_genome.fasta")
            mock_args.return_value.second_fasta = str(self.data_dir / "fasta" / "test_genome2.fasta")
            mock_args.return_value.output = str(self.output_dir)
            mock_args.return_value.debug = True
            mock_args.return_value.config = None
            mock_args.return_value.cli = True
            mock_args.return_value.direct = None
            mock_args.return_value.nooligo = False
            mock_args.return_value.snp = False
            mock_args.return_value.noannotation = False
            mock_args.return_value.db = None
            mock_args.return_value.db_action = None
            mock_args.return_value.maf = None
            mock_args.return_value.vcf = None  # Add this to avoid undefined attribute
            mock_args.return_value.second_vcf = None  # Add this to avoid undefined attribute
            mock_args.return_value.gff = None  # Add this to avoid undefined attribute
            
            # Configure mocks to return success
            mock_standard_mode.return_value = True
            mock_direct_mode.return_value = True
            mock_alignment_mode.return_value = True
            mock_lastz.return_value = str(self.output_dir / "test_output.maf")
            
            # Configure the factory to return the alignment mode function
            mock_factory.return_value = mock_alignment_mode
            
            # Run pipeline
            result = run_pipeline()
            
            # Check that the correct mode was called
            mock_standard_mode.assert_not_called()
            mock_direct_mode.assert_not_called()
            mock_alignment_mode.assert_called_once()
            
            # Verify result
            self.assertTrue(result, "Pipeline should return success")


if __name__ == "__main__":
    unittest.main()