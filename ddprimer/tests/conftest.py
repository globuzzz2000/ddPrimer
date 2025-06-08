#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared pytest fixtures for ddPrimer tests.

This file contains fixtures that can be reused across multiple test modules.
"""

import os
import sys
import pytest
import tempfile
import shutil
import logging
import pandas as pd
from pathlib import Path
import warnings
from unittest import mock



# Import package modules
from ..config import Config

# Add the parent directory to the path so we can import the package
# This is needed when running tests from the tests directory
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)


# ============== Suppress all logging completely ===============
# Aggressively suppress ALL logging from third-party packages
# This creates a null handler that will be applied to all loggers
class NullHandler(logging.Handler):
    def emit(self, record):
        pass


# Completely silence loggers by disabling propagation and adding null handler
def silence_logger(logger_name):
    """Completely silence a logger by name."""
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.CRITICAL)
    logger.propagate = False
    logger.handlers = []
    logger.addHandler(NullHandler())


# Silence all nupack loggers and subloggers
silence_logger("nupack")
silence_logger("nupack.rebind.render")
silence_logger("nupack.design")
silence_logger("nupack.constants")
silence_logger("nupack.model")
silence_logger("nupack.derivations")
silence_logger("nupack.core")
silence_logger("nupack.analysis")
silence_logger("nupack.rebind")

# Silence other potentially noisy loggers
silence_logger("matplotlib")
silence_logger("PIL")
silence_logger("tqdm")
silence_logger("fontTools")

# Set root logger to only show errors or higher
root_logger = logging.getLogger()
root_logger.setLevel(logging.ERROR)

# Set ddPrimer logger to show only warnings and above
ddprimer_logger = logging.getLogger("ddPrimer")
ddprimer_logger.setLevel(logging.WARNING)
ddprimer_logger.propagate = False


# ============== Filter out all warnings ===============
# Filter out all warnings to reduce output noise
warnings.filterwarnings("ignore")


def pytest_configure(config):
    """
    Configure pytest settings at startup.
    """
    # Disable all logging messages during testing
    logging.disable(logging.CRITICAL)
    
    # Disable all warnings during testing
    warnings.simplefilter("ignore")


def pytest_runtest_setup(item):
    """
    Set up logging before each test runs.
    """
    # Reset log levels before each test to ensure consistency
    logging.disable(logging.CRITICAL)
    
    # Re-silence all loggers before each test
    for logger_name in ["nupack", "matplotlib", "PIL", "tqdm", "fontTools"]:
        silence_logger(logger_name)


# Redirect stdout/stderr for extremely noisy libraries
@pytest.fixture(autouse=True)
def redirect_output():
    """Redirect stdout and stderr during test execution."""
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    
    try:
        with open(os.devnull, 'w') as devnull:
            sys.stdout = devnull
            sys.stderr = devnull
            yield
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr


# ============== Test data directories ===============

@pytest.fixture(scope="session")
def test_data_dir():
    """Return the path to the test data directory."""
    return os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(scope="session")
def fasta_data_dir(test_data_dir):
    """Return the path to the test FASTA data directory."""
    return os.path.join(test_data_dir, "fasta")


@pytest.fixture(scope="session")
def vcf_data_dir(test_data_dir):
    """Return the path to the test VCF data directory."""
    return os.path.join(test_data_dir, "vcf")


@pytest.fixture(scope="session")
def gff_data_dir(test_data_dir):
    """Return the path to the test GFF data directory."""
    return os.path.join(test_data_dir, "gff")


@pytest.fixture(scope="session")
def direct_data_dir(test_data_dir):
    """Return the path to the test direct mode data directory."""
    return os.path.join(test_data_dir, "direct")


@pytest.fixture(scope="session")
def alignment_data_dir(test_data_dir):
    """Return the path to the test alignment mode data directory."""
    return os.path.join(test_data_dir, "alignment")


# ============== Temporary files and directories ===============

@pytest.fixture(scope="function")
def temp_dir():
    """Create a temporary directory and clean it up after the test."""
    temp_dir = tempfile.mkdtemp(prefix="ddprimer_test_")
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def test_fasta_file():
    """
    Fixture that creates a temporary FASTA file.
    Returns the path to the file, which is automatically cleaned up.
    """
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as tmp_file:
        # Write some test sequences
        tmp_file.write(">seq1\n")
        tmp_file.write("ATCGATCGATCGTAGCTAGCTAGC\n")
        tmp_file.write(">seq2\n")
        tmp_file.write("GCTAGCTAGCTAGCTAGCTA\n")
        tmp_file.write(">seq3\n")
        tmp_file.write("ATGCATGCATGCATGCATGC\n")
        tmp_file.flush()
        tmp_path = tmp_file.name
    
    yield tmp_path
    
    # Clean up
    if os.path.exists(tmp_path):
        os.unlink(tmp_path)


@pytest.fixture
def test_gff_file():
    """
    Fixture that creates a temporary GFF file.
    Returns the path to the file, which is automatically cleaned up.
    """
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.gff', delete=False) as tmp_file:
        # Write some test GFF content
        tmp_file.write("##gff-version 3\n")
        tmp_file.write("seq1\tGenbankParser\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=GENE1\n")
        tmp_file.write("seq1\tGenbankParser\texon\t10\t90\t.\t+\t.\tID=exon1;Parent=gene1\n")
        tmp_file.write("seq2\tGenbankParser\tgene\t1\t200\t.\t-\t.\tID=gene2;Name=GENE2\n")
        tmp_file.write("seq3\tGenbankParser\tgene\t50\t150\t.\t+\t.\tID=gene3;Name=AT1G12345\n")
        tmp_file.flush()
        tmp_path = tmp_file.name
    
    yield tmp_path
    
    # Clean up
    if os.path.exists(tmp_path):
        os.unlink(tmp_path)


@pytest.fixture
def test_vcf_file():
    """
    Fixture that creates a temporary VCF file.
    Returns the path to the file, which is automatically cleaned up.
    """
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.vcf', delete=False) as tmp_file:
        # Write some test VCF content
        tmp_file.write("##fileformat=VCFv4.2\n")
        tmp_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        tmp_file.write("seq1\t10\t.\tA\tG\t100\tPASS\t.\n")
        tmp_file.write("seq1\t50\t.\tC\tT\t100\tPASS\t.\n")
        tmp_file.write("seq2\t30\t.\tG\tA\t100\tPASS\t.\n")
        tmp_file.flush()
        tmp_path = tmp_file.name
    
    yield tmp_path
    
    # Clean up
    if os.path.exists(tmp_path):
        os.unlink(tmp_path)


# ============== Test data and mock objects ===============

@pytest.fixture(scope="function")
def config():
    """Return a fresh instance of the Config class for testing."""
    # Reset config to defaults
    Config._instance = None
    config_instance = Config.get_instance()
    
    # Override settings for testing
    config_instance.DEBUG_MODE = True
    config_instance.NUM_PROCESSES = 2  # Use fewer processes for testing
    config_instance.SHOW_PROGRESS = False  # Disable progress bars in tests
    
    return config_instance


@pytest.fixture(scope="function")
def mock_config():
    """
    Fixture to mock the Config object with a more direct approach.
    Resets any changes made during tests.
    """
    # Store original values
    original_values = {}
    for attr in dir(Config):
        if not attr.startswith('__'):
            original_values[attr] = getattr(Config, attr)
    
    # Set common test values
    Config.SHOW_PROGRESS = False
    Config.NUM_PROCESSES = 2
    Config.RETAIN_TYPES = ['gene', 'exon', 'cds']
    Config.FILTER_MEANINGFUL_NAMES = True
    Config.RESTRICTION_SITE = "GGCC"
    Config.MIN_SEGMENT_LENGTH = 20
    Config.GENE_OVERLAP_MARGIN = 50
    Config.THERMO_TEMPERATURE = 37
    Config.THERMO_SODIUM = 0.05
    Config.THERMO_MAGNESIUM = 0.002
    Config.THERMO_BACKEND = 'vienna'
    Config.THERMO_AUTO_FALLBACK = True
    Config.BLAST_WORD_SIZE = 7
    Config.BLAST_EVALUE = 1000
    Config.BLAST_FILTER_FACTOR = 2
    Config.DB_PATH = "test_blast_db"
    
    yield Config
    
    # Restore original values
    for attr, value in original_values.items():
        setattr(Config, attr, value)


@pytest.fixture(scope="function")
def simple_sequences():
    """Return a dictionary of simple test sequences."""
    return {
        "seq1": "ATGCATGCATGCGGCCATGCATGCATGC",
        "seq2": "GTACGTACGTACGGCCGTACGTACGTAC",
        "seq3": "GGCCTTTGGCCAAGGCCTTTGGCCAAA"
    }


@pytest.fixture
def test_sequences():
    """
    Fixture that returns a dictionary of test sequences.
    """
    return {
        "seq1": "ATCGATCGATCGTAGCTAGCTAGCGGCCATGCATGC",
        "seq2": "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG",
        "seq3": "ATGCATGCATGCGGCCATGCATGCATGCATGC"
    }


@pytest.fixture(scope="function")
def mock_genes():
    """Return a dictionary of mock gene annotations."""
    return {
        "chr1": [
            {"gene": "geneA", "start": 100, "end": 500},
            {"gene": "geneB", "start": 700, "end": 900}
        ],
        "chr2": [
            {"gene": "geneC", "start": 200, "end": 600}
        ]
    }


@pytest.fixture
def test_genes():
    """
    Fixture that returns a list of test genes.
    """
    return [
        {"id": "gene1", "chr": "seq1", "start": 10, "end": 25, "strand": "+"},
        {"id": "gene2", "chr": "seq2", "start": 5, "end": 15, "strand": "-"},
        {"id": "gene3", "chr": "seq4", "start": 1, "end": 10, "strand": "+"}
    ]


@pytest.fixture
def test_fragments():
    """
    Fixture that returns a list of test sequence fragments.
    """
    return [
        {"id": "seq1_frag0", "chr": "seq1", "start": 1, "end": 8, "sequence": "ATCGATCG"},
        {"id": "seq1_frag1", "chr": "seq1", "start": 13, "end": 20, "sequence": "TAGCTAG"},
        {"id": "seq2", "chr": "seq2", "start": 1, "end": 32, "sequence": "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"}
    ]


@pytest.fixture(scope="function")
def mock_primer3_output():
    """Return mock Primer3 output."""
    return {
        "PRIMER_PAIR_0_PENALTY": "0.5",
        "PRIMER_LEFT_0_SEQUENCE": "ATGCATGCATGC",
        "PRIMER_RIGHT_0_SEQUENCE": "GTACGTACGTAC",
        "PRIMER_INTERNAL_0_SEQUENCE": "TTTGGGCCCAAA",
        "PRIMER_LEFT_0": "0,12",
        "PRIMER_RIGHT_0": "27,12",
        "PRIMER_INTERNAL_0": "13,12",
        "PRIMER_PAIR_0_PRODUCT_SIZE": "28",
        "PRIMER_LEFT_0_TM": "60.0",
        "PRIMER_RIGHT_0_TM": "60.0",
        "PRIMER_INTERNAL_0_TM": "65.0",
        "PRIMER_LEFT_0_GC_PERCENT": "50.0",
        "PRIMER_RIGHT_0_GC_PERCENT": "50.0",
        "PRIMER_INTERNAL_0_GC_PERCENT": "55.0",
    }


@pytest.fixture(scope="function")
def mock_blast_output():
    """Return mock BLAST output."""
    return [
        {
            "target": "chr1",
            "start": 100,
            "end": 112,
            "evalue": 1e-10,
            "identity": 100.0
        },
        {
            "target": "chr2",
            "start": 300,
            "end": 312,
            "evalue": 1e-5,
            "identity": 90.0
        }
    ]


@pytest.fixture
def test_primers_df():
    """
    Fixture that returns a test DataFrame of primer candidates.
    """
    return pd.DataFrame({
        'ID': ['primer1', 'primer2', 'primer3', 'primer4'],
        'Forward Primer': ['ATCGATCGATCG', 'GCTAGCTAGCTA', 'ATGCATGCATGC', 'CGTAGCTAGCTA'],
        'Forward Tm': [60.5, 58.2, 62.1, 59.8],
        'Reverse Primer': ['GCTAGCTAGCTA', 'ATCGATCGATCG', 'GCATGCATGCAT', 'TAGCTAGCTACG'],
        'Reverse Tm': [58.2, 60.5, 61.9, 59.8],
        'Probe': ['CGATCGATCGTA', 'TAGCTAGCTAGC', 'CATGCATGCATG', 'TAGCTAGCTACG'],
        'Probe Tm': [65.3, 63.8, 66.2, 64.5],
        'Amplicon': ['ATCGATCGATCGTAGCTAGCTAGCTA', 'GCTAGCTAGCTATAGCTAGCTAGC', 
                    'ATGCATGCATGCATGCATGCATGCAT', 'CGTAGCTAGCTATAGCTAGCTACG'],
        'Amplicon Length': [26, 24, 27, 24],
        'Penalty': [0.3, 0.5, 0.2, 0.6]
    })