"""
Tests for sequence_utils module
"""
import pytest
from ...utils.sequence_utils import SequenceUtils


class TestSequenceUtils:
    """Tests for SequenceUtils class."""

    def test_has_disallowed_repeats(self):
        """Test checking for disallowed repeats in sequences."""
        # Test sequences with disallowed repeats
        assert SequenceUtils.has_disallowed_repeats("ATGCCCCTGA") is True  # Has CCCC
        assert SequenceUtils.has_disallowed_repeats("ATGGGGGTGA") is True  # Has GGGG
        assert SequenceUtils.has_disallowed_repeats("CCCCATGGGGG") is True  # Has both CCCC and GGGG
        
        # Test sequences without disallowed repeats
        assert SequenceUtils.has_disallowed_repeats("ATGCCTGA") is False
        assert SequenceUtils.has_disallowed_repeats("ATGCGCTGA") is False
        assert SequenceUtils.has_disallowed_repeats("CCCATGGG") is False
        
        # Test edge cases
        assert SequenceUtils.has_disallowed_repeats("") is False
        assert SequenceUtils.has_disallowed_repeats(None) is True
        assert SequenceUtils.has_disallowed_repeats(123) is True  # Non-string input

    def test_calculate_gc(self):
        """Test calculating GC content of sequences."""
        # Test various GC contents
        assert SequenceUtils.calculate_gc("ATGC") == 50.0  # 50% GC
        assert SequenceUtils.calculate_gc("AAAA") == 0.0   # 0% GC
        assert SequenceUtils.calculate_gc("GGCC") == 100.0 # 100% GC
        assert SequenceUtils.calculate_gc("ATCG") == 50.0  # 50% GC
        assert SequenceUtils.calculate_gc("ATGCAT") == 33.33333333333333  # 33.33% GC
        
        # Test case insensitivity
        assert SequenceUtils.calculate_gc("atgc") == 50.0
        assert SequenceUtils.calculate_gc("ATgc") == 50.0
        
        # Test edge cases
        assert SequenceUtils.calculate_gc("") == 0.0
        assert SequenceUtils.calculate_gc(None) == 0.0
        assert SequenceUtils.calculate_gc(123) == 0.0  # Non-string input

    def test_reverse_complement(self):
        """Test generating reverse complement of sequences."""
        # Test standard cases
        assert SequenceUtils.reverse_complement("ATGC") == "GCAT"
        assert SequenceUtils.reverse_complement("AAAA") == "TTTT"
        assert SequenceUtils.reverse_complement("GGCC") == "GGCC"  # Palindromic
        assert SequenceUtils.reverse_complement("GATC") == "GATC"  # Palindromic
        
        # Test with ambiguous bases - fix the expected output
        assert SequenceUtils.reverse_complement("ATGCN") == "NGCAT"
        assert SequenceUtils.reverse_complement("RYSWKM") == "KMWSRY"
        
        # Test case insensitivity
        assert SequenceUtils.reverse_complement("atgc") == "GCAT"
        assert SequenceUtils.reverse_complement("AtGc") == "GCAT"
        
        # Test edge cases
        assert SequenceUtils.reverse_complement("") == ""
        assert SequenceUtils.reverse_complement(None) == ""
        assert SequenceUtils.reverse_complement(123) == ""  # Non-string input

    def test_ensure_more_c_than_g(self):
        """Test ensuring sequence has more C than G."""
        # Test sequences with more C than G
        seq, reversed = SequenceUtils.ensure_more_c_than_g("ACCCGT")  # 3C, 1G
        assert seq == "ACCCGT"
        assert reversed is False
        
        # Test sequences with equal C and G
        seq, reversed = SequenceUtils.ensure_more_c_than_g("ACCGGT")  # 2C, 2G
        assert seq == "ACCGGT"
        assert reversed is False
        
        # Test sequences with more G than C
        # The implementation should reverse complement this sequence
        # and set reversed to True
        seq, reversed = SequenceUtils.ensure_more_c_than_g("AGGGCT")  # 1C, 3G
        expected_rev_comp = "AGGGCT"[::-1].translate(str.maketrans("ATGC", "TACG"))
        assert seq == expected_rev_comp
        assert reversed is True
        
        # Test case insensitivity
        seq, reversed = SequenceUtils.ensure_more_c_than_g("agggct")  # 1c, 3g
        expected_rev_comp = "AGGGCT"[::-1].translate(str.maketrans("ATGC", "TACG"))
        assert seq == expected_rev_comp
        assert reversed is True
        
        # Test edge cases
        seq, reversed = SequenceUtils.ensure_more_c_than_g("")
        assert seq == ""
        assert reversed is False
        
        seq, reversed = SequenceUtils.ensure_more_c_than_g(None)
        assert seq == None
        assert reversed is False
        
        seq, reversed = SequenceUtils.ensure_more_c_than_g(123)
        assert seq == 123
        assert reversed is False