"""
Unit tests for the primer_processor module
"""
import pytest
import pandas as pd
import numpy as np
from ...core.primer_processor import PrimerProcessor
from ...config import Config


class TestPrimerProcessor:
    """Tests for PrimerProcessor class."""
    
    def test_filter_by_penalty(self):
        """Test filtering primers by penalty threshold."""
        # Create a test DataFrame
        data = {
            'Gene': ['Gene1', 'Gene2', 'Gene3'],
            'Primer F': ['ATGCGGCCAT', 'CTAGCCGGTA', 'ACGTACGTAC'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC', 'GTACGTACGT'],
            'Penalty': [0.2, 1.5, 0.8]
        }
        df = pd.DataFrame(data)
        
        # Set the penalty threshold
        Config.PENALTY_MAX = 1.0
        
        # Run the filtering function
        filtered_df = PrimerProcessor.filter_by_penalty(df)
        
        # Check the results - only primers with penalty < 1.0 should remain
        assert len(filtered_df) == 2
        assert 'Gene1' in filtered_df['Gene'].values
        assert 'Gene3' in filtered_df['Gene'].values
        assert 'Gene2' not in filtered_df['Gene'].values
    
    def test_filter_by_repeats(self):
        """Test filtering primers with repeat sequences."""
        # Create a test DataFrame with repeats
        data = {
            'Gene': ['Gene1', 'Gene2', 'Gene3'],
            'Primer F': ['ATGCATGCAT', 'AAAAAAATCG', 'ACGTACGTAC'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC', 'GGGGGGACTG']
        }
        df = pd.DataFrame(data)
        
        # Run the filtering function
        filtered_df = PrimerProcessor.filter_by_repeats(df)
        
        # Check the results - primers with repeats should be removed
        assert len(filtered_df) == 1
        assert 'Gene1' in filtered_df['Gene'].values
        assert 'Gene2' not in filtered_df['Gene'].values
        assert 'Gene3' not in filtered_df['Gene'].values
    
    def test_filter_by_gc_content(self):
        """Test filtering primers by GC content."""
        # Create a test DataFrame with different GC contents
        data = {
            'Gene': ['Gene1', 'Gene2', 'Gene3'],
            'Primer F': ['ATGCATGCAT', 'CGCGCGCGCG', 'ATATAGAGTA'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC', 'TAAATAAATA'],
            'Primer F GC': [40.0, 100.0, 20.0],
            'Primer R GC': [50.0, 60.0, 0.0]
        }
        df = pd.DataFrame(data)
        
        # Set GC content limits
        Config.SEQUENCE_MIN_GC = 30.0
        Config.SEQUENCE_MAX_GC = 70.0
        
        # Run the filtering function
        filtered_df = PrimerProcessor.filter_by_gc_content(df)
        
        # Check the results - only primers with GC content within limits should remain
        assert len(filtered_df) == 1
        assert 'Gene1' in filtered_df['Gene'].values
        assert 'Gene2' not in filtered_df['Gene'].values
        assert 'Gene3' not in filtered_df['Gene'].values
    
    def test_process_internal_oligos(self):
        """Test processing of internal oligos."""
        # Create a test DataFrame with probes that need processing
        data = {
            'Gene': ['Gene1', 'Gene2'],
            'Primer F': ['ATGCATGCAT', 'CGCGCGCGCG'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC'],
            'Probe': ['ATGCTAGCAT', 'CGCTAAAGCG'],
            'Probe Reverse Complement': [True, False]
        }
        df = pd.DataFrame(data)
        
        # Run the processing function
        processed_df = PrimerProcessor.process_internal_oligos(df)
        
        # Check the results - probes with Reverse Complement = True should be reverse complemented
        assert processed_df.loc[0, 'Probe'] == 'ATGCTAGCAT'[::-1].translate(str.maketrans('ATGC', 'TACG'))
        assert processed_df.loc[1, 'Probe'] == 'CGCTAAAGCG'
    
    def test_filter_by_blast(self):
        """Test filtering primers by BLAST results."""
        # Create a test DataFrame with BLAST results
        data = {
            'Gene': ['Gene1', 'Gene2', 'Gene3'],
            'Primer F': ['ATGCATGCAT', 'CGCGCGCGCG', 'ACGTACGTAC'],
            'Primer R': ['TTAAGCCAGT', 'TGGCTAGGCC', 'GTACGTACGT'],
            'Primer F BLAST1': [1, 2, 5],
            'Primer F BLAST2': [0, 1, 3],
            'Primer R BLAST1': [1, 3, 1],
            'Primer R BLAST2': [0, 2, 0]
        }
        df = pd.DataFrame(data)
        
        # Set the BLAST filter factor
        Config.BLAST_FILTER_FACTOR = 2
        
        # Run the filtering function
        filtered_df = PrimerProcessor.filter_by_blast(df)
        
        # Check the results - primers with BLAST hits > threshold should be removed
        assert len(filtered_df) == 1
        assert 'Gene1' in filtered_df['Gene'].values
        assert 'Gene2' not in filtered_df['Gene'].values
        assert 'Gene3' not in filtered_df['Gene'].values