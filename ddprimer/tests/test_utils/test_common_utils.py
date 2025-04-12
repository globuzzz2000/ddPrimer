"""
Tests for common_utils module
"""
import pytest
from ...utils.common_utils import CommonUtils


class TestCommonUtils:
    """Tests for CommonUtils class."""

    def test_chunks(self):
        """Test splitting a list into chunks."""
        # Test with even division
        test_list = [1, 2, 3, 4, 5, 6]
        chunk_size = 2
        result = CommonUtils.chunks(test_list, chunk_size)
        
        assert len(result) == 3
        assert result == [[1, 2], [3, 4], [5, 6]]
        
        # Test with uneven division
        test_list = [1, 2, 3, 4, 5, 6, 7]
        chunk_size = 2
        result = CommonUtils.chunks(test_list, chunk_size)
        
        assert len(result) == 4
        assert result == [[1, 2], [3, 4], [5, 6], [7]]
        
        # Test with chunk size larger than list
        test_list = [1, 2, 3]
        chunk_size = 5
        result = CommonUtils.chunks(test_list, chunk_size)
        
        assert len(result) == 1
        assert result == [[1, 2, 3]]
        
        # Test with empty list
        test_list = []
        chunk_size = 2
        result = CommonUtils.chunks(test_list, chunk_size)
        
        assert len(result) == 0
        assert result == []

    def test_flatten_list(self):
        """Test flattening a list of lists."""
        # Test with regular nested list
        nested_list = [[1, 2], [3, 4], [5, 6]]
        result = CommonUtils.flatten_list(nested_list)
        
        assert len(result) == 6
        assert result == [1, 2, 3, 4, 5, 6]
        
        # Test with uneven nested list
        nested_list = [[1], [2, 3, 4], [5, 6]]
        result = CommonUtils.flatten_list(nested_list)
        
        assert len(result) == 6
        assert result == [1, 2, 3, 4, 5, 6]
        
        # Test with empty lists
        nested_list = [[], [1, 2], []]
        result = CommonUtils.flatten_list(nested_list)
        
        assert len(result) == 2
        assert result == [1, 2]
        
        # Test with completely empty list
        nested_list = []
        result = CommonUtils.flatten_list(nested_list)
        
        assert len(result) == 0
        assert result == []
        
        # Test with nested empty lists
        nested_list = [[], []]
        result = CommonUtils.flatten_list(nested_list)
        
        assert len(result) == 0
        assert result == []
        
        # Test with strings
        nested_list = [["a", "b"], ["c", "d"]]
        result = CommonUtils.flatten_list(nested_list)
        
        assert len(result) == 4
        assert result == ["a", "b", "c", "d"]