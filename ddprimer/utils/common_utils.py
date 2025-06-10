#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common utility functions for ddPrimer pipeline.

Contains functionality for:
1. List processing and manipulation utilities
2. Data structure transformation functions
3. General-purpose helper functions
4. Input validation and error handling

This module provides shared utility functions that can be used
across different components of the ddPrimer pipeline.
"""

import logging

# Set up module logger
logger = logging.getLogger(__name__)


class CommonUtils:
    """
    General utility functions used across the pipeline.
    
    This class provides static methods for common operations like
    list chunking, flattening, and other data manipulation tasks
    that are frequently needed throughout the pipeline.
    
    Example:
        >>> chunks = CommonUtils.chunks([1, 2, 3, 4, 5], 2)
        >>> print(chunks)  # [[1, 2], [3, 4], [5]]
    """
    
    @staticmethod
    def chunks(lst, n):
        """
        Split a list into chunks of specified size.
        
        Divides a list into smaller sublists of the specified size,
        with the last chunk potentially being smaller if the list
        length is not evenly divisible by the chunk size.
        
        Args:
            lst: List to split into chunks
            n: Chunk size (must be positive integer)
            
        Returns:
            List of chunks from the original list
            
        Raises:
            TypeError: If lst is not a list
            ValueError: If n is not a positive integer
            
        Example:
            >>> CommonUtils.chunks([1, 2, 3, 4, 5], 2)
            [[1, 2], [3, 4], [5]]
        """
        if not isinstance(lst, list):
            error_msg = f"Input to chunks must be a list, got {type(lst).__name__}"
            logger.error(error_msg)
            raise TypeError(error_msg)
            
        if not isinstance(n, int) or n <= 0:
            error_msg = f"Chunk size must be a positive integer, got {n}"
            logger.error(error_msg)
            raise ValueError(error_msg)
            
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Splitting list of {len(lst)} items into chunks of size {n}")
            
        chunks_list = []
        for i in range(0, len(lst), n):
            chunk = lst[i:i + n]
            chunks_list.append(chunk)
            
        logger.debug(f"Created {len(chunks_list)} chunks")
        return chunks_list
    
    @staticmethod
    def flatten_list(list_of_lists):
        """
        Flatten a list of lists into a single list.
        
        Takes a nested list structure and converts it into a flat
        single-level list by concatenating all sublists.
        
        Args:
            list_of_lists: List containing other lists
                
        Returns:
            Flattened list containing all elements from sublists
            
        Raises:
            TypeError: If input is not a list or contains non-list elements
            
        Example:
            >>> CommonUtils.flatten_list([[1, 2], [3, 4], [5]])
            [1, 2, 3, 4, 5]
        """
        if not isinstance(list_of_lists, list):
            error_msg = f"Input to flatten_list must be a list, got {type(list_of_lists).__name__}"
            logger.error(error_msg)
            raise TypeError(error_msg)
            
        # Check that each element is also a list
        for i, element in enumerate(list_of_lists):
            if not isinstance(element, list):
                error_msg = f"Expected list of lists, but found element of type {type(element).__name__} at index {i}"
                logger.error(error_msg)
                raise TypeError(error_msg)
                
        if logger.isEnabledFor(logging.DEBUG):
            total_sublists = len(list_of_lists)
            total_elements = sum(len(sublist) for sublist in list_of_lists)
            logger.debug(f"Flattening {total_sublists} sublists containing {total_elements} total elements")
            
        flattened = [item for sublist in list_of_lists for item in sublist]
        logger.debug(f"Flattened to list of {len(flattened)} elements")
        return flattened