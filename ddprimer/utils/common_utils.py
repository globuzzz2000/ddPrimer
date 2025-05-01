#!/usr/bin/env python3
"""
Common utility functions for ddPrimer pipeline.
"""

import logging

class CommonUtils:
    """General utility functions used across the pipeline."""
    
    @staticmethod
    def chunks(lst, n):
        """
        Split a list into chunks of size n.
        
        Args:
            lst (list): List to split
            n (int): Chunk size
            
        Returns:
            list: List of chunks from the original list
        """
        if not isinstance(lst, list):
            logger = logging.getLogger("ddPrimer")
            logger.warning("Input to chunks must be a list")
            return []
            
        if not isinstance(n, int) or n <= 0:
            logger = logging.getLogger("ddPrimer")
            logger.warning("Chunk size must be a positive integer")
            return [lst]
            
        chunks_list = []
        for i in range(0, len(lst), n):
            chunks_list.append(lst[i:i + n])
        return chunks_list
    
    @staticmethod
    def flatten_list(list_of_lists):
        """
        Flatten a list of lists into a single list.
        
        Args:
            list_of_lists (list): List of lists
                
        Returns:
            list: Flattened list
            
        Raises:
            TypeError: If any element in list_of_lists is not a list
        """
        if not isinstance(list_of_lists, list):
            logger = logging.getLogger("ddPrimer")
            logger.warning("Input to flatten_list must be a list")
            return []
            
        # Check that each element is also a list
        for element in list_of_lists:
            if not isinstance(element, list):
                raise TypeError(f"Expected list of lists, but found element of type {type(element).__name__}")
                
        return [item for sublist in list_of_lists for item in sublist]