#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for sequence utility functions in the ddPrimer pipeline.

These tests verify the functionality of the SequenceUtils class,
which provides sequence-specific operations for DNA sequences.
"""

import os
import sys
import pytest

# Import the module to test
from ...utils import SequenceUtils
from ...config import Config


class TestSequenceUtils:
    """Test class for SequenceUtils functions."""
    
    def test_has_disallowed_repeats(self):
        """Test detection of disallowed repeats in sequences."""
        # Test sequences with disallowed repeats
        assert SequenceUtils.has_disallowed_repeats("ATCGCCCCATGC") == True, "Should detect CCCC repeat"
        assert SequenceUtils.has_disallowed_repeats("GGGGTATACGTA") == True, "Should detect GGGG repeat"
        assert SequenceUtils.has_disallowed_repeats("ATCGGGGGCTA") == True, "Should detect GGGG repeat"
        
        # Test sequences without disallowed repeats
        assert SequenceUtils.has_disallowed_repeats("ATCGCCCATGC") == False, "Should not detect CCC as disallowed"
        assert SequenceUtils.has_disallowed_repeats("GGGTATACGTA") == False, "Should not detect GGG as disallowed"
        assert SequenceUtils.has_disallowed_repeats("ATCGAATTCGTA") == False, "Should not detect any repeats"
        
        # Test edge cases
        assert SequenceUtils.has_disallowed_repeats("") == False, "Empty sequence should have no repeats"
        assert SequenceUtils.has_disallowed_repeats(None) == True, "None should be considered invalid"
        assert SequenceUtils.has_disallowed_repeats(123) == True, "Non-string should be considered invalid"
    
    def test_calculate_gc(self):
        """Test GC content calculation in DNA sequences."""
        # Test with various realistic sequences
        assert SequenceUtils.calculate_gc("GCGCGCGCGC") == 100.0, "Should be 100% GC"
        assert SequenceUtils.calculate_gc("ATATATATAT") == 0.0, "Should be 0% GC"
        assert SequenceUtils.calculate_gc("ATGC") == 50.0, "Should be 50% GC"
        assert SequenceUtils.calculate_gc("AGTC") == 50.0, "Should be 50% GC"
        
        # Real biological sequences
        promoter_seq = "TACAAAGGACGTGGGTTCGAATCCCACCTCTGCCG"
        assert round(SequenceUtils.calculate_gc(promoter_seq), 1) == 57.1, "Promoter sequence should have ~57.1% GC"
        
        enhancer_seq = "AGATTGCAGATCTGATTGCCAGCTAGCATGCAT"
        assert round(SequenceUtils.calculate_gc(enhancer_seq), 1) == 45.5, "Enhancer sequence should have ~45.5% GC"
        
        # Test case insensitivity
        assert SequenceUtils.calculate_gc("atgc") == 50.0, "Should handle lowercase"
        assert SequenceUtils.calculate_gc("AtGc") == 50.0, "Should handle mixed case"
        
        # Test edge cases
        assert SequenceUtils.calculate_gc("") == 0, "Empty string should return 0"
        assert SequenceUtils.calculate_gc(None) == 0, "None should return 0"
        assert SequenceUtils.calculate_gc(123) == 0, "Non-string should return 0"
    
    def test_reverse_complement(self):
        """Test reverse complement generation for DNA sequences."""
        # Test with simple sequences
        assert SequenceUtils.reverse_complement("ATGC") == "GCAT", "Should correctly reverse complement ATGC"
        assert SequenceUtils.reverse_complement("GCTA") == "TAGC", "Should correctly reverse complement GCTA"
        
        # Test with longer biological sequences
        primer_seq = "ACGTACCTAGCTAGCTAGCTACGT"
        expected_rev_comp = "ACGTAGCTAGCTAGCTAGGTACGT"
        assert SequenceUtils.reverse_complement(primer_seq) == expected_rev_comp, "Should correctly reverse complement primer"
        
        exon_seq = "ATGGCAAGCTATGCAAGCTGACTAGATCGATCGTAGCTAGC"
        expected_rev_comp = "GCTAGCTACGATCGATCTAGTCAGCTTGCATAGCTTGCCAT"
        assert SequenceUtils.reverse_complement(exon_seq) == expected_rev_comp, "Should correctly reverse complement exon"
        
        # Test with ambiguous bases
        assert SequenceUtils.reverse_complement("ATGCRYSWKMBDHVN") == "NBDHVKMWSRYGCAT", "Should handle ambiguous bases"
        
        # Test case preservation (should convert to uppercase)
        assert SequenceUtils.reverse_complement("atgc") == "GCAT", "Should convert to uppercase"
        
        # Test edge cases
        assert SequenceUtils.reverse_complement("") == "", "Empty string should return empty string"
        assert SequenceUtils.reverse_complement(None) == "", "None should return empty string"
        assert SequenceUtils.reverse_complement(123) == "", "Non-string should return empty string"
    
    def test_ensure_more_c_than_g(self):
        """Test ensuring sequences have more Cs than Gs."""
        # Test sequence with more Cs than Gs already
        seq = "ACCTGCCATGC"  # 5 Cs, 3 Gs
        result, was_reversed = SequenceUtils.ensure_more_c_than_g(seq)
        assert result == seq, "Sequence with more Cs should not be reversed"
        assert was_reversed == False, "Should not flag as reversed"
        
        # Test sequence with equal Cs and Gs
        seq = "ACGTACGT"  # 2 Cs, 2 Gs
        result, was_reversed = SequenceUtils.ensure_more_c_than_g(seq)
        assert result == seq, "Sequence with equal Cs and Gs should not be reversed"
        assert was_reversed == False, "Should not flag as reversed"
        
        # Test sequence with more Gs than Cs
        seq = "AGGGTGCATG"  # 2 Cs, 5 Gs
        result, was_reversed = SequenceUtils.ensure_more_c_than_g(seq)
        expected = "CATGCACCCT"  # Reverse complement has 5 Cs, 2 Gs
        assert result == expected, "Sequence with more Gs should be reverse complemented"
        assert was_reversed == True, "Should flag as reversed"
        
        # Test with real biological sequences
        # Sequence with G-rich strand
        g_rich_seq = "GGGCTAGGGCTGGGACTGAGCT"  # 5 Cs, 10 Gs
        result, was_reversed = SequenceUtils.ensure_more_c_than_g(g_rich_seq)
        expected = "AGCTCAGTCCCAGCCCTAGCCC"  # Reverse complement has 10 Cs, 5 Gs
        assert result == expected, "G-rich sequence should be reverse complemented"
        assert was_reversed == True, "Should flag as reversed"
        
        # Sequence with C-rich strand
        c_rich_seq = "CCCGATCCCGACCCTGACGCT"  # 10 Cs, 5 Gs
        result, was_reversed = SequenceUtils.ensure_more_c_than_g(c_rich_seq)
        assert result == c_rich_seq, "C-rich sequence should not be reversed"
        assert was_reversed == False, "Should not flag as reversed"
        
        # Test edge cases
        assert SequenceUtils.ensure_more_c_than_g("") == ("", False), "Empty string should return empty string"
        assert SequenceUtils.ensure_more_c_than_g(None) == (None, False), "None should return None"
        assert SequenceUtils.ensure_more_c_than_g(123) == (123, False), "Non-string should return original"