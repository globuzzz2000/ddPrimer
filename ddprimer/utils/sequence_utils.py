#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence utility functions for ddPrimer pipeline.

Contains functionality for:
1. DNA sequence analysis and validation
2. GC content calculation
3. Reverse complement generation
4. Repeat sequence detection

This module provides sequence-specific utility functions for analyzing
and manipulating DNA sequences throughout the ddPrimer pipeline.
"""

import logging

# Set up module logger
logger = logging.getLogger(__name__)


class SequenceUtils:
    """
    Sequence-specific utility functions for DNA analysis.
    
    This class provides static methods for common DNA sequence operations
    including GC content calculation, reverse complement generation,
    and repeat sequence detection.
    
    Example:
        >>> gc_content = SequenceUtils.calculate_gc("ATCGCGTA")
        >>> print(f"GC content: {gc_content:.1f}%")
        >>> rev_comp = SequenceUtils.reverse_complement("ATCG")
        >>> print(f"Reverse complement: {rev_comp}")
    """
    
    @staticmethod
    def has_disallowed_repeats(seq):
        """
        Check for disallowed repeats in a DNA sequence.
        
        Identifies sequences containing runs of four or more consecutive
        G or C nucleotides, which can cause issues in PCR amplification
        and primer binding.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            True if disallowed repeats found, False otherwise
            
        Example:
            >>> SequenceUtils.has_disallowed_repeats("ATCGGGGATC")
            True
            >>> SequenceUtils.has_disallowed_repeats("ATCGATCG")
            False
        """
        if not isinstance(seq, str):
            logger.debug("Invalid sequence type provided to has_disallowed_repeats")
            return True
            
        if not seq:
            logger.debug("Empty sequence provided to has_disallowed_repeats")
            return False
            
        seq_upper = seq.upper()
        has_repeats = "CCCC" in seq_upper or "GGGG" in seq_upper
        
        if logger.isEnabledFor(logging.DEBUG) and has_repeats:
            logger.debug(f"Disallowed repeats found in sequence: {seq[:20]}...")
            
        return has_repeats
    
    @staticmethod
    def calculate_gc(seq):
        """
        Calculate GC content of a DNA sequence.
        
        Computes the percentage of guanine and cytosine nucleotides
        in the provided DNA sequence, which is important for primer
        design and PCR optimization.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            GC content as a percentage (0-100)
            
        Example:
            >>> SequenceUtils.calculate_gc("ATCGGCTA")
            37.5
        """
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to calculate_gc")
            return 0.0
            
        if not seq.strip():
            logger.debug("Empty sequence provided to calculate_gc")
            return 0.0
            
        seq_upper = seq.upper()
        gc_count = sum(1 for base in seq_upper if base in "GC")
        total_bases = len(seq_upper)
        
        if total_bases == 0:
            return 0.0
            
        gc_percentage = (gc_count / total_bases) * 100
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"GC content calculation: {gc_count}/{total_bases} = {gc_percentage:.1f}%")
            
        return gc_percentage
    
    @staticmethod
    def reverse_complement(seq):
        """
        Generate the reverse complement of a DNA sequence.
        
        Creates the reverse complement by reversing the sequence and
        replacing each nucleotide with its complement (A↔T, G↔C).
        Supports IUPAC ambiguous nucleotide codes.
        
        Args:
            seq: DNA sequence to reverse complement
            
        Returns:
            Reverse complement sequence
            
        Raises:
            ValueError: If sequence contains invalid characters
            
        Example:
            >>> SequenceUtils.reverse_complement("ATCG")
            'CGAT'
        """
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to reverse_complement")
            return ""
        
        if not seq.strip():
            logger.debug("Empty sequence provided to reverse_complement")
            return ""
        
        seq_upper = seq.upper()
        
        # Define the complement mapping including IUPAC codes
        complement = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
            'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
            'D': 'H', 'H': 'D', 'V': 'B'
        }
        
        # Check for invalid characters
        invalid_chars = set(seq_upper) - set(complement.keys())
        if invalid_chars:
            error_msg = f"Invalid nucleotide characters found: {', '.join(invalid_chars)}"
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        # Generate reverse complement
        try:
            rev_comp = ''.join(complement.get(base, base) for base in reversed(seq_upper))
            
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(f"Reverse complement: {seq[:20]}... -> {rev_comp[:20]}...")
                
            return rev_comp
        except Exception as e:
            error_msg = f"Error generating reverse complement for sequence: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ValueError(error_msg) from e
    
    @staticmethod
    def ensure_more_c_than_g(seq):
        """
        Ensure a sequence has more Cs than Gs, reversing if necessary.
        
        Checks the C/G ratio in a sequence and returns the reverse
        complement if there are more Gs than Cs. This is useful for
        probe design where C-rich sequences are preferred.
        
        Args:
            seq: DNA sequence to analyze
            
        Returns:
            Tuple of (possibly_reversed_sequence, was_reversed)
            
        Example:
            >>> seq, reversed = SequenceUtils.ensure_more_c_than_g("GGGATC")
            >>> print(f"Sequence: {seq}, Was reversed: {reversed}")
        """
        if not seq or not isinstance(seq, str):
            logger.debug("Invalid sequence provided to ensure_more_c_than_g")
            return seq, False
        
        if not seq.strip():
            logger.debug("Empty sequence provided to ensure_more_c_than_g")
            return seq, False
        
        seq_upper = seq.upper()
        c_count = seq_upper.count('C')
        g_count = seq_upper.count('G')
        
        logger.debug(f"C/G analysis: C={c_count}, G={g_count}")
        
        if c_count >= g_count:
            logger.debug("C count >= G count, no reversal needed")
            return seq, False
        
        # Need to reverse complement
        try:
            rev_comp = SequenceUtils.reverse_complement(seq)
            logger.debug("Sequence reversed to ensure more Cs than Gs")
            return rev_comp, True
        except ValueError as e:
            error_msg = f"Failed to reverse complement sequence: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            # Return original sequence if reverse complement fails
            return seq, False
        except Exception as e:
            error_msg = f"Unexpected error in ensure_more_c_than_g: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return seq, False