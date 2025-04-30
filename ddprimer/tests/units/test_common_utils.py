#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for common utility functions in the ddPrimer pipeline.

These tests verify the basic functionality of utility classes like CommonUtils.
"""

import os
import sys
import pytest
from pathlib import Path

# Import the module to test
from ...utils import CommonUtils


class TestCommonUtils:
    """Test class for CommonUtils functions."""
    
    def test_chunks_normal_case(self):
        """Test chunking a list with normal inputs."""
        # Test with a list of 10 items and chunk size 3
        input_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = CommonUtils.chunks(input_list, 3)
        
        assert len(result) == 4, "Should create 4 chunks"
        assert result[0] == [1, 2, 3], "First chunk should contain first 3 elements"
        assert result[1] == [4, 5, 6], "Second chunk should contain next 3 elements"
        assert result[2] == [7, 8, 9], "Third chunk should contain next 3 elements"
        assert result[3] == [10], "Last chunk should contain remaining elements"
    
    def test_chunks_empty_list(self):
        """Test chunking an empty list."""
        result = CommonUtils.chunks([], 5)
        assert result == [], "Should return empty list for empty input"
    
    def test_chunks_chunk_size_equal_to_list_length(self):
        """Test chunking when chunk size equals list length."""
        input_list = [1, 2, 3, 4, 5]
        result = CommonUtils.chunks(input_list, 5)
        
        assert len(result) == 1, "Should create 1 chunk"
        assert result[0] == input_list, "The chunk should be the entire list"
    
    def test_chunks_chunk_size_greater_than_list_length(self):
        """Test chunking when chunk size is greater than list length."""
        input_list = [1, 2, 3]
        result = CommonUtils.chunks(input_list, 5)
        
        assert len(result) == 1, "Should create 1 chunk"
        assert result[0] == input_list, "The chunk should be the entire list"
    
    def test_chunks_invalid_inputs(self):
        """Test chunking with invalid inputs."""
        # Test with non-list input
        result = CommonUtils.chunks("not a list", 3)
        assert result == [], "Should return empty list for non-list input"
        
        # Test with invalid chunk size
        input_list = [1, 2, 3, 4, 5]
        result = CommonUtils.chunks(input_list, 0)
        assert result == [input_list], "Should return original list for chunk size <= 0"
        
        result = CommonUtils.chunks(input_list, -1)
        assert result == [input_list], "Should return original list for chunk size <= 0"
        
        result = CommonUtils.chunks(input_list, "not a number")
        assert result == [input_list], "Should return original list for non-integer chunk size"
    
    def test_flatten_list_normal_case(self):
        """Test flattening a list of lists with normal inputs."""
        input_list = [[1, 2, 3], [4, 5], [6, 7, 8, 9], [10]]
        result = CommonUtils.flatten_list(input_list)
        
        assert result == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "Should flatten all sublists"
    
    def test_flatten_list_empty_list(self):
        """Test flattening an empty list."""
        result = CommonUtils.flatten_list([])
        assert result == [], "Should return empty list for empty input"
    
    def test_flatten_list_single_item_sublists(self):
        """Test flattening list of single-item sublists."""
        input_list = [[1], [2], [3], [4], [5]]
        result = CommonUtils.flatten_list(input_list)
        
        assert result == [1, 2, 3, 4, 5], "Should flatten single-item sublists"
    
    def test_flatten_list_with_empty_sublists(self):
        """Test flattening list with some empty sublists."""
        input_list = [[1, 2], [], [3, 4], [], [5]]
        result = CommonUtils.flatten_list(input_list)
        
        assert result == [1, 2, 3, 4, 5], "Should handle empty sublists"
    
    def test_flatten_list_invalid_input(self):
        """Test flattening with invalid input."""
        # Test with non-list input
        result = CommonUtils.flatten_list("not a list")
        assert result == [], "Should return empty list for non-list input"
        
        # Test with non-list items in the list
        input_list = [[1, 2], "not a list", [3, 4]]
        with pytest.raises(TypeError):
            CommonUtils.flatten_list(input_list)