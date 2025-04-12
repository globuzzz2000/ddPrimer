"""
Test fixtures for ddPrimer testing
"""
import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

@pytest.fixture
def test_data_dir():
    """Return the test data directory path."""
    # Assuming tests are in ddPrimer/tests and test_data is in ddPrimer/test_data
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_data')


@pytest.fixture
def mock_sequence_dict():
    """Create a mock sequence dictionary."""
    return {
        'seq1': 'ATGCGGCCACTGGCTTAAGGCCTTATAGGATTCGCCATAGCGGCCTTAA',
        'seq2': 'CTAGCCGGTAATCGGGCCGAATTCGATTCGGCCTAGCCA',
        'seq3': 'ACGTACGTACGTGGCCAATCGATCGTACGTAC'
    }


@pytest.fixture
def mock_restriction_fragments():
    """Create a mock list of restriction fragments."""
    return [
        {
            'id': 'seq1_frag1',
            'sequence': 'ATGCGGCCACTGGCTTAA',
            'chr': 'chr1',
            'start': 1,
            'end': 18,
            'Gene': 'Gene1'
        },
        {
            'id': 'seq1_frag2',
            'sequence': 'GGCCTTATAGGATTCGCCATAGCGGCCTTAA',
            'chr': 'chr1',
            'start': 19,
            'end': 49,
            'Gene': 'Gene1'
        },
        {
            'id': 'seq2_frag1',
            'sequence': 'CTAGCCGGTAATCGGGCCGAATTCGATTCGGCCTAGCCA',
            'chr': 'chr2',
            'start': 1,
            'end': 39,
            'Gene': 'Gene2'
        }
    ]


@pytest.fixture
def mock_gene_annotations():
    """Create mock gene annotation data."""
    return {
        'chr1': [
            {'start': 1, 'end': 30, 'id': 'Gene1', 'strand': '+'},
            {'start': 35, 'end': 60, 'id': 'Gene3', 'strand': '-'}
        ],
        'chr2': [
            {'start': 1, 'end': 40, 'id': 'Gene2', 'strand': '+'}
        ]
    }


@pytest.fixture
def mock_primer3_inputs():
    """Create mock Primer3 input data."""
    return [
        {
            'SEQUENCE_ID': 'seq1_frag1',
            'SEQUENCE_TEMPLATE': 'ATGCGGCCACTGGCTTAA'
        },
        {
            'SEQUENCE_ID': 'seq2_frag1',
            'SEQUENCE_TEMPLATE': 'CTAGCCGGTAATCGGGCCGAATTCGATTCGGCCTAGCCA',
            'SEQUENCE_TARGET': [10, 20]
        }
    ]


@pytest.fixture
def mock_primer3_results():
    """Create mock primer results from Primer3."""
    return [
        {
            'Gene': 'Gene1',
            'Fragment ID': 'seq1_frag1',
            'Chromosome': 'chr1',
            'Location': '1-18',
            'Primer F': 'ATGCGGCCACTGGCT',
            'Primer R': 'TTAAGCCAGTGGCCGCAT',
            'Primer F TM': 60.2,
            'Primer R TM': 61.5,
            'Amplicon': 'ATGCGGCCACTGGCTTAA',
            'Amplicon Length': 18,
            'Penalty': 0.2,
            'Probe': 'GGCCACTGGCT'
        },
        {
            'Gene': 'Gene2',
            'Fragment ID': 'seq2_frag1',
            'Chromosome': 'chr2',
            'Location': '1-39',
            'Primer F': 'CTAGCCGGTAATCGGGC',
            'Primer R': 'TGGCTAGGCCGAATCGAAT',
            'Primer F TM': 59.8,
            'Primer R TM': 60.3,
            'Amplicon': 'CTAGCCGGTAATCGGGCCGAATTCGATTCGGCCTAGCCA',
            'Amplicon Length': 39,
            'Penalty': 0.5,
            'Probe': 'ATCGGGCCGAATTCG'
        }
    ]


@pytest.fixture
def mock_coordinate_map():
    """Create a mock coordinate mapping for cross-species testing."""
    coordinate_map = {
        'chr1': {
            1: {'qry_src': 'qry_chr1', 'qry_pos': 100},
            2: {'qry_src': 'qry_chr1', 'qry_pos': 101},
            3: {'qry_src': 'qry_chr1', 'qry_pos': 102}
        },
        'chr2': {
            1: {'qry_src': 'qry_chr2', 'qry_pos': 200},
            2: {'qry_src': 'qry_chr2', 'qry_pos': 201}
        }
    }
    return coordinate_map


@pytest.fixture
def create_mini_fasta(tmp_path):
    """Create a mini FASTA file for testing."""
    def _create_mini_fasta(sequences=None):
        if sequences is None:
            sequences = {
                'chr1': 'ATGCGGCCACTGGCTTAAGGCCTTATAGGATTCGCCATAGCGGCCTTAA',
                'chr2': 'CTAGCCGGTAATCGGGCCGAATTCGATTCGGCCTAGCCA'
            }
        
        fasta_path = tmp_path / "mini_genome.fasta"
        
        with open(fasta_path, 'w') as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n{seq}\n")
        
        return str(fasta_path)
    
    return _create_mini_fasta


@pytest.fixture
def create_mini_vcf(tmp_path):
    """Create a mini VCF file for testing."""
    def _create_mini_vcf():
        vcf_path = tmp_path / "mini_variants.vcf"
        
        with open(vcf_path, 'w') as f:
            f.write("##fileformat=VCFv4.1\n")
            f.write("##contig=<ID=chr1,length=50>\n")
            f.write("##contig=<ID=chr2,length=40>\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("chr1\t10\t.\tG\tA\t100\tPASS\t.\n")
            f.write("chr1\t25\t.\tA\tT\t100\tPASS\t.\n")
            f.write("chr2\t15\t.\tC\tG\t100\tPASS\t.\n")
        
        return str(vcf_path)
    
    return _create_mini_vcf


@pytest.fixture
def create_mini_gff(tmp_path):
    """Create a mini GFF file for testing."""
    def _create_mini_gff():
        gff_path = tmp_path / "mini_annotations.gff"
        
        with open(gff_path, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("chr1\t.\tgene\t1\t30\t.\t+\t.\tID=Gene1;Name=Gene1\n")
            f.write("chr1\t.\tgene\t35\t60\t.\t-\t.\tID=Gene3;Name=Gene3\n")
            f.write("chr2\t.\tgene\t1\t40\t.\t+\t.\tID=Gene2;Name=Gene2\n")
        
        return str(gff_path)
    
    return _create_mini_gff


@pytest.fixture
def create_mini_direct_input(tmp_path):
    """Create a mini direct input Excel file for testing."""
    def _create_mini_direct_input():
        import pandas as pd
        
        df = pd.DataFrame({
            'Sequence Name': ['seq1', 'seq2', 'seq3'],
            'Sequence': [
                'ATGCGGCCACTGGCTTAAGGCCTTATAGGATTCGCCATAGCGGCCTTAA',
                'CTAGCCGGTAATCGGGCCGAATTCGATTCGGCCTAGCCA',
                'ACGTACGTACGTGGCCAATCGATCGTACGTAC'
            ]
        })
        
        excel_path = tmp_path / "mini_direct.xlsx"
        df.to_excel(excel_path, index=False)
        
        return str(excel_path)
    
    return _create_mini_direct_input