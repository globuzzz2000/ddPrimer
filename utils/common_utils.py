#!/usr/bin/env python3
"""
Common utility functions for ddPrimer pipeline.
"""


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
        """
        return [item for sublist in list_of_lists for item in sublist]
