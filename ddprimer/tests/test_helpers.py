"""
Helper functions for ddPrimer tests.
"""
import os
import pandas as pd
from io import StringIO


def get_test_data_path(filename=None):
    """
    Get the absolute path to the test_data directory or a file within it.
    
    Args:
        filename (str, optional): Name of file within test_data directory.
        
    Returns:
        str: Absolute path to the test_data directory or file.
    """
    # Get the absolute path to the test_data directory
    test_data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_data')
    
    if filename:
        # Return the path to the specified file
        return os.path.join(test_data_dir, filename)
    else:
        # Return the path to the directory
        return test_data_dir


def create_test_csv(data_str):
    """
    Create a pandas DataFrame from a CSV string.
    
    Args:
        data_str (str): CSV data as a string.
        
    Returns:
        pandas.DataFrame: DataFrame created from the CSV data.
    """
    return pd.read_csv(StringIO(data_str))


def compare_dataframes(df1, df2, columns=None):
    """
    Compare two DataFrames for equality, optionally limiting to specific columns.
    
    Args:
        df1 (pandas.DataFrame): First DataFrame.
        df2 (pandas.DataFrame): Second DataFrame.
        columns (list, optional): List of columns to compare. If None, compare all columns.
        
    Returns:
        bool: True if DataFrames are equal for the specified columns, False otherwise.
    """
    if columns:
        return df1[columns].equals(df2[columns])
    else:
        return df1.equals(df2)