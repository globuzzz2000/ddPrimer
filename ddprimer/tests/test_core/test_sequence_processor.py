"""
Unit tests for the sequence_processor module
"""
import pytest
from ...core.sequence_processor import SequenceProcessor
from ...config import Config


class TestSequenceProcessor:
    """Tests for SequenceProcessor class."""
    
    def test_cut_at_restriction_sites(self, mock_sequence_dict):
        """Test cutting sequences at restriction sites."""
        # Set the restriction site to "GGCC" for testing
        Config.RESTRICTION_SITE = "GGCC"
        
        # Run the function
        fragments = SequenceProcessor.cut_at_restriction_sites(mock_sequence_dict)
        
        # Check the results
        assert isinstance(fragments, list)
        assert len(fragments) > 0
        
        # Check specific fragments for seq1 which has 2 "GGCC" sites
        seq1_fragments = [f for f in fragments if f['id'].startswith('seq1')]
        assert len(seq1_fragments) == 3  # Should be split into 3 fragments
        
        # Verify fragment contents
        frag_seqs = [f['sequence'] for f in seq1_fragments]
        assert 'ATGC' in frag_seqs
        assert 'ACTGGCTTAA' in frag_seqs
        assert 'TTATAGGATTCGCCATAGC' in frag_seqs
        assert 'TTAA' in frag_seqs
        
        # Check that each fragment has the necessary keys
        for fragment in fragments:
            assert 'id' in fragment
            assert 'sequence' in fragment
            assert 'chr' in fragment
            assert 'start' in fragment
            assert 'end' in fragment
    
    def test_filter_by_gene_overlap(self, mock_restriction_fragments, mock_gene_annotations):
        """Test filtering fragments by gene overlap."""
        # Configure the gene overlap margin
        Config.GENE_OVERLAP_MARGIN = 5
        
        # Run the function
        filtered_fragments = SequenceProcessor.filter_by_gene_overlap(mock_restriction_fragments, mock_gene_annotations)
        
        # Check the results
        assert isinstance(filtered_fragments, list)
        assert len(filtered_fragments) > 0
        
        # All fragments should be within a single gene boundary
        for fragment in filtered_fragments:
            gene_id = fragment.get('Gene')
            chrom = fragment.get('chr')
            start = fragment.get('start')
            end = fragment.get('end')
            
            # Check if this fragment is within a gene's boundaries
            if chrom in mock_gene_annotations:
                matching_genes = [g for g in mock_gene_annotations[chrom] 
                                if g['id'] == gene_id and 
                                   start >= g['start'] - Config.GENE_OVERLAP_MARGIN and
                                   end <= g['end'] + Config.GENE_OVERLAP_MARGIN]
                                   
                assert len(matching_genes) > 0, f"Fragment {fragment['id']} does not fit within gene {gene_id}"
    
    def test_mask_variants(self):
        """Test masking variants in a sequence."""
        # Create a test sequence and variants
        sequence = "ATGCGGCCACTGGCTTAA"
        variants = [
            {'pos': 5, 'ref': 'G', 'alt': 'A'},
            {'pos': 10, 'ref': 'C', 'alt': 'T'}
        ]
        
        # Run the masking function
        masked_sequence = SequenceProcessor.mask_variants(sequence, variants)
        
        # Check that the variants were masked correctly (replaced with N)
        expected_masked = "ATGCNGGCCANTGGCTTAA"
        assert masked_sequence == expected_masked
    
    def test_no_variants_to_mask(self):
        """Test masking with no variants."""
        sequence = "ATGCGGCCACTGGCTTAA"
        variants = []
        
        masked_sequence = SequenceProcessor.mask_variants(sequence, variants)
        
        # No masking should occur
        assert masked_sequence == sequence