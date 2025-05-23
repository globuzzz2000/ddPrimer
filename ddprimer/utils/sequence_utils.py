#!/usr/bin/env python3
"""
Sequence utility functions for ddPrimer pipeline.
"""
import re
import logging
from ..config import Config

class SequenceUtils:
    """Sequence-specific utility functions."""
    
    @staticmethod
    def has_disallowed_repeats(seq):
        """
        Check for disallowed repeats in a sequence.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            bool: True if disallowed repeats found, False otherwise
        """
        logger = logging.getLogger("ddPrimer")
        
        if not isinstance(seq, str):
            logger.debug("Invalid sequence type provided to has_disallowed_repeats")
            return True
            
        return "CCCC" in seq or "GGGG" in seq
    
    @staticmethod
    def calculate_gc(seq):
        """
        Calculate GC content of a sequence.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: GC content percentage
        """
        logger = logging.getLogger("ddPrimer")
        
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to calculate_gc")
            return 0
            
        seq = seq.upper()
        gc_count = sum(1 for base in seq if base in "GC")
        return (gc_count / len(seq)) * 100 if seq else 0
    
    @staticmethod
    def reverse_complement(seq):
        """
        Generate the reverse complement of a DNA sequence.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            str: Reverse complement sequence
        """
        logger = logging.getLogger("ddPrimer")
        
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to reverse_complement")
            return ""
        
        seq = seq.upper()
        
        # Define the complement mapping
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                     'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
                     'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                     'D': 'H', 'H': 'D', 'V': 'B'}
        
        # Generate reverse complement
        rev_comp = ''.join(complement.get(base, base) for base in reversed(seq))
        return rev_comp
    
    @staticmethod
    def ensure_more_c_than_g(seq):
        """
        Check if a sequence has more Cs than Gs. If not, return the reverse complement.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            tuple: (possibly_reversed_sequence, was_reversed)
        """
        logger = logging.getLogger("ddPrimer")
        
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to ensure_more_c_than_g")
            return seq, False
        
        seq = seq.upper()
        c_count = seq.count('C')
        g_count = seq.count('G')
        
        if c_count >= g_count:
            return seq, False  # No need to reverse
        
        # Need to reverse complement - set was_reversed to True
        rev_comp = SequenceUtils.reverse_complement(seq)
        return rev_comp, True