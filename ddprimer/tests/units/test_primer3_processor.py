#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the Primer3Processor module.
"""

import pytest
import re
import primer3
from unittest.mock import patch, MagicMock, call

from ...core import Primer3Processor
from ...config import SequenceProcessingError


class TestPrimer3Processor:
    """Test cases for the Primer3Processor class."""

    @pytest.fixture
    def processor(self, mock_config):
        """Create a Primer3Processor instance with mocked config."""
        return Primer3Processor(config=mock_config)

    @pytest.fixture
    def sample_input_block(self):
        """Create a sample Primer3 input block."""
        return {
            "SEQUENCE_ID": "test_seq_1",
            "SEQUENCE_TEMPLATE": "ATGCATGCATGCGGCCATGCATGCATGCATGC",
            "SEQUENCE_TARGET": [10, 10]
        }

    @pytest.fixture
    def sample_primer3_result(self):
        """Create a sample Primer3 result dictionary."""
        return {
            "PRIMER_PAIR_0_PENALTY": 0.5,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 120,
            "PRIMER_LEFT_0_SEQUENCE": "ATGCATGCATGC",
            "PRIMER_LEFT_0": (0, 12),
            "PRIMER_LEFT_0_TM": 60.0,
            "PRIMER_LEFT_0_PENALTY": 0.2,
            "PRIMER_RIGHT_0_SEQUENCE": "GCATGCATGCAT",
            "PRIMER_RIGHT_0": (130, 12),
            "PRIMER_RIGHT_0_TM": 59.5,
            "PRIMER_RIGHT_0_PENALTY": 0.3,
            "PRIMER_INTERNAL_0_SEQUENCE": "GCGGCCATGCAT",
            "PRIMER_INTERNAL_0": (50, 12),
            "PRIMER_INTERNAL_0_TM": 65.2,
            "PRIMER_INTERNAL_0_PENALTY": 0.1
        }

    def test_init(self, mock_config):
        """Test initialization of Primer3Processor."""
        processor = Primer3Processor(config=mock_config)
        assert processor.config == mock_config
        
        # Test with default config
        with patch('ddprimer.core.primer3_processor.Config') as default_config:
            processor = Primer3Processor()
            assert processor.config == default_config

    def test_format_primer3_result(self, processor, sample_input_block, sample_primer3_result):
        """Test formatting of Primer3 results to match primer3_core output format."""
        formatted = processor._format_primer3_result(
            sample_input_block["SEQUENCE_ID"],
            sample_input_block,
            sample_primer3_result
        )
        
        # Check basic structure
        assert "SEQUENCE_ID=test_seq_1" in formatted
        assert "SEQUENCE_TEMPLATE=" + sample_input_block["SEQUENCE_TEMPLATE"] in formatted
        assert "PRIMER_PAIR_NUM_RETURNED=1" in formatted
        assert "PRIMER_PAIR_0_PENALTY=0.5" in formatted
        assert "PRIMER_LEFT_0_SEQUENCE=ATGCATGCATGC" in formatted
        assert "PRIMER_RIGHT_0_SEQUENCE=GCATGCATGCAT" in formatted
        assert "PRIMER_INTERNAL_0_SEQUENCE=GCGGCCATGCAT" in formatted
        assert formatted.endswith("=")
        
        # Check coordinate conversion is correct
        # In primer3-py, coordinates are 0-based tuples (start, length)
        # Should be converted to "start,length" strings
        assert "PRIMER_LEFT_0=0,12" in formatted
        assert "PRIMER_RIGHT_0=130,12" in formatted
        assert "PRIMER_INTERNAL_0=50,12" in formatted

    @patch('primer3.bindings.design_primers')
    def test_run_primer3_batch(self, mock_design_primers, processor, sample_input_block, sample_primer3_result):
        """Test running Primer3 on a batch of input blocks."""
        # Setup mock
        mock_design_primers.return_value = sample_primer3_result
        
        # Run with a single input block
        result = processor.run_primer3_batch([sample_input_block])
        
        # Check primer3.bindings.design_primers was called correctly
        mock_design_primers.assert_called_once()
        args, kwargs = mock_design_primers.call_args
        assert kwargs["seq_args"] == sample_input_block
        
        # Check output format
        assert "SEQUENCE_ID=test_seq_1" in result
        assert "PRIMER_PAIR_NUM_RETURNED=1" in result
        assert "PRIMER_PAIR_0_PENALTY=0.5" in result

    @patch('primer3.bindings.design_primers')
    def test_run_primer3_batch_parallel(self, mock_design_primers, processor, sample_input_block, sample_primer3_result):
        """Test running Primer3 on multiple input blocks in parallel."""
        # Setup mock
        mock_design_primers.return_value = sample_primer3_result
        
        # Create multiple input blocks
        blocks = [
            {**sample_input_block, "SEQUENCE_ID": f"test_seq_{i}"} 
            for i in range(5)
        ]
        
        # Run with parallel processing
        with patch('ddprimer.core.primer3_processor.ProcessPoolExecutor') as mock_executor:
            # Setup mock executor
            mock_executor_instance = MagicMock()
            mock_executor.return_value.__enter__.return_value = mock_executor_instance
            
            # Setup future results
            mock_future = MagicMock()
            mock_future.result.return_value = "SEQUENCE_ID=test_seq_0\n=\n"
            mock_executor_instance.submit.return_value = mock_future
            
            # Run the parallel batch
            result = processor.run_primer3_batch_parallel(blocks, max_workers=2)
            
            # Check executor was used correctly
            assert mock_executor.call_args[1]["max_workers"] == 2
            assert mock_executor_instance.submit.call_count > 0

    def test_get_amplicon(self, processor):
        """Test extraction of amplicon sequence from template."""
        # Define a test sequence
        seq = "ATGCATGCATGCGGCCATGCATGCATGCATGC"
        
        # Test normal case
        amplicon = processor.get_amplicon(seq, 5, 4, 20, 4)
        # Amplicon should start at position 4 (1-based, adjusted) and end at position 20
        expected = seq[4-1:20]
        assert amplicon == expected
        
        # Test edge cases
        # 1. Left start position is 1 (minimum)
        amplicon = processor.get_amplicon(seq, 1, 4, 20, 4)
        expected = seq[0:20]  # 0-based indexing for Python
        assert amplicon == expected
        
        # 2. Right position is beyond sequence length
        amplicon = processor.get_amplicon(seq, 5, 4, 100, 4)
        # The function adjusts left_start by -1, so we need to start at index 3 (5-1-1)
        expected = seq[4-1:len(seq)]  # Adjusted to match function behavior
        assert amplicon == expected
        
        # 3. Invalid positions
        assert processor.get_amplicon(seq, -5, 4, 20, 4) == seq[0:20]  # Adjusted to valid range
        assert processor.get_amplicon(seq, 25, 4, 20, 4) == ""  # Invalid (start > end)
        
        # 4. None values
        assert processor.get_amplicon(seq, None, 4, 20, 4) == ""
        assert processor.get_amplicon(seq, 5, None, 20, 4) == ""
        assert processor.get_amplicon(seq, 5, 4, None, 4) == ""
        assert processor.get_amplicon(seq, 5, 4, 20, None) == ""

    def test_parse_primer3_batch(self, processor):
        """Test parsing of Primer3 output."""
        # Create sample Primer3 output
        primer3_output = """SEQUENCE_ID=test_seq_1
SEQUENCE_TEMPLATE=ATGCATGCATGCGGCCATGCATGCATGCATGC
PRIMER_PAIR_NUM_RETURNED=2
PRIMER_PAIR_0_PENALTY=0.5
PRIMER_PAIR_0_PRODUCT_SIZE=120
PRIMER_LEFT_0_SEQUENCE=ATGCATGCATGC
PRIMER_LEFT_0=0,12
PRIMER_LEFT_0_TM=60.0
PRIMER_LEFT_0_PENALTY=0.2
PRIMER_RIGHT_0_SEQUENCE=GCATGCATGCAT
PRIMER_RIGHT_0=119,12
PRIMER_RIGHT_0_TM=59.5
PRIMER_RIGHT_0_PENALTY=0.3
PRIMER_INTERNAL_0_SEQUENCE=GCGGCCATGCAT
PRIMER_INTERNAL_0=50,12
PRIMER_INTERNAL_0_TM=65.2
PRIMER_INTERNAL_0_PENALTY=0.1
PRIMER_PAIR_1_PENALTY=0.8
PRIMER_PAIR_1_PRODUCT_SIZE=100
PRIMER_LEFT_1_SEQUENCE=ATGCGGCCATGC
PRIMER_LEFT_1=8,12
PRIMER_LEFT_1_TM=61.0
PRIMER_LEFT_1_PENALTY=0.3
PRIMER_RIGHT_1_SEQUENCE=GCATGCATGCAT
PRIMER_RIGHT_1=107,12
PRIMER_RIGHT_1_TM=58.5
PRIMER_RIGHT_1_PENALTY=0.5
PRIMER_INTERNAL_1_SEQUENCE=ATGCATGCATGC
PRIMER_INTERNAL_1=40,12
PRIMER_INTERNAL_1_TM=63.2
PRIMER_INTERNAL_1_PENALTY=0.2
=
SEQUENCE_ID=test_seq_2
SEQUENCE_TEMPLATE=GCTAGCTAGCTAGCTAGCTA
PRIMER_PAIR_NUM_RETURNED=1
PRIMER_PAIR_0_PENALTY=0.6
PRIMER_PAIR_0_PRODUCT_SIZE=80
PRIMER_LEFT_0_SEQUENCE=GCTAGCTAG
PRIMER_LEFT_0=0,9
PRIMER_LEFT_0_TM=58.0
PRIMER_LEFT_0_PENALTY=0.3
PRIMER_RIGHT_0_SEQUENCE=TAGCTAGCT
PRIMER_RIGHT_0=79,9
PRIMER_RIGHT_0_TM=57.5
PRIMER_RIGHT_0_PENALTY=0.3
PRIMER_INTERNAL_0_SEQUENCE=CTAGCTAGC
PRIMER_INTERNAL_0=30,9
PRIMER_INTERNAL_0_TM=62.2
PRIMER_INTERNAL_0_PENALTY=0.1
="""

        # Create fragment info for mapping sequence IDs to chromosomes
        fragment_info = {
            "test_seq_1": {"chr": "chr1", "start": 1000, "gene": "Gene1"},
            "test_seq_2": {"chr": "chr2", "start": 2000, "gene": "Gene2"}
        }
        
        # Patch config settings
        with patch.object(processor.config, 'DISABLE_INTERNAL_OLIGO', False), \
             patch.object(processor.config, 'PREFER_PROBE_MORE_C_THAN_G', True), \
             patch.object(processor.config, 'MAX_PRIMER_PAIRS_PER_SEGMENT', 2):
            
            # Parse the output
            records = processor.parse_primer3_batch(primer3_output, fragment_info)
            
            # Check basic structure and count
            assert len(records) == 3  # 2 from first sequence, 1 from second
            
            # Check first record details
            assert records[0]["Gene"] == "Gene1"
            assert records[0]["Index"] == 0
            assert records[0]["Primer F"] == "ATGCATGCATGC"
            assert records[0]["Tm F"] == 60.0
            assert records[0]["Primer R"] == "GCATGCATGCAT"
            assert records[0]["Tm R"] == 59.5
            assert records[0]["Probe"] == "GCGGCCATGCAT"
            assert records[0]["Probe Tm"] == 65.2
            assert records[0]["Pair Penalty"] == 0.5
            assert records[0]["Chromosome"] == "chr1"
            
            # Check absolute position calculation (fragment_start + primer_start - 1)
            assert "1000" in records[0]["Location"]
            
            # Check second sequence's record
            assert records[2]["Gene"] == "Gene2"
            assert records[2]["Chromosome"] == "chr2"
            
        # Test with DISABLE_INTERNAL_OLIGO = True
        with patch.object(processor.config, 'DISABLE_INTERNAL_OLIGO', True), \
             patch.object(processor.config, 'MAX_PRIMER_PAIRS_PER_SEGMENT', 2):
            
            records = processor.parse_primer3_batch(primer3_output, fragment_info)
            
            # Check probe fields are not included
            assert "Probe" not in records[0]
            assert "Probe Tm" not in records[0]
            assert "Probe Penalty" not in records[0]

    def test_parse_primer3_batch_empty(self, processor):
        """Test parsing empty Primer3 output."""
        # Test with empty output
        records = processor.parse_primer3_batch("", {})
        assert records == []
        
        # Test with invalid output
        records = processor.parse_primer3_batch("INVALID OUTPUT", {})
        assert records == []