# ddPrimer: Advanced Droplet Digital PCR Primer Design

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)

A comprehensive pipeline for designing primers and probes specifically optimized for droplet digital PCR (ddPCR).

## Key Features

- **Complete End-to-End Pipeline**: Design primers from genome sequences through a streamlined workflow using Primer3
- **Smart SNP Masking**: Avoid designing primers across variant positions using VCF files with intelligent AF-based processing
- **Thermodynamic Optimization**: Calculate ΔG values using ViennaRNA to prevent unwanted secondary structures
- **Specificity Verification**: Integrated BLAST validation for both primers and probes
- **File Preparation**: Automatic VCF normalization, chromosome mapping, and file indexing
- **Standard Workflow**: Design from genome, variant and annotation files with comprehensive processing
- **Comprehensive Results**: Detailed Excel output with key metrics and sequence information

## Project Structure

```
ddPrimer/
├── ddprimer/                      # Main package directory
│   ├── __init__.py
│   ├── main.py                    # Main entry point for the pipeline
│   ├── core/                      # Core processing modules
│   │   ├── __init__.py
│   │   ├── snp_processor.py       # VCF-based variant processing
│   │   ├── annotation_processor.py
│   │   ├── blast_processor.py
│   │   ├── vienna_processor.py    # ViennaRNA thermodynamic processor
│   │   ├── primer3_processor.py
│   │   ├── primer_processor.py
│   │   └── sequence_processor.py
│   ├── utils/                     # Utility functions
│   │   ├── __init__.py
│   │   ├── blast_db_manager.py    # Unified BLAST database management
│   │   ├── file_io.py             # File I/O and Excel formatting
│   │   └── file_preparator.py     # File validation and preparation
│   └── config/                    # Configuration and settings
│       ├── __init__.py
│       ├── config.py              # Core configuration settings
│       ├── config_display.py      # Configuration display utilities
│       ├── exceptions.py          # Custom exceptions
│       ├── logging_config.py      # Logging setup
│       └── template_generator.py  # Configuration template generation
├── pyproject.toml                 # Package configuration and dependencies
└── README.md                      # This file
```

## Installation

### Quick Install with Conda (Recommended)

```bash
# Clone and install
git clone https://github.com/jakobmueller/ddPrimer
cd ddPrimer

# Create environment with all dependencies
conda create -n ddprimer python=3.8
conda activate ddprimer
conda install -c bioconda -c conda-forge blast bcftools samtools viennarna
pip install -e .
```

### Alternative: pip install

```bash
git clone https://github.com/jakobmueller/ddPrimer
cd ddPrimer
pip install -e .

# Then install external tools via system package manager:
# macOS: brew install blast bcftools samtools viennarna
# Linux: sudo apt-get install ncbi-blast+ bcftools samtools
```

### Required External Tools

- **NCBI BLAST+**: For specificity checking
- **bcftools/samtools**: For file processing  
- **ViennaRNA**: For thermodynamic calculations

Python dependencies are automatically installed via pip.

## Quick Start

### Command Line Usage

```bash
# Basic primer design with file preparation
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff

# Basic primer design without gene filtering
ddprimer --noannotation --fasta genome.fasta --vcf variants.vcf
```

### Interactive Mode

Simply run `ddprimer` without arguments to launch the interactive mode, which will guide you through file selection with a graphical interface.

## Workflow Overview

1. **Input Selection**: Choose genome FASTA, variant VCF, and annotation GFF files
2. **File Preparation**: Validate and prepare files (bgzip compression, indexing, normalization, chromosome mapping)
3. **Variant Processing**: Apply VCF variants to sequences with intelligent AF-based masking/substitution
4. **Sequence Preparation**: Filter sequences based on restriction sites and gene boundaries
5. **Primer Design**: Design primer and probe candidates using Primer3
6. **Quality Filtering**: Apply filters for penalties, repeats, GC content, and more
7. **Thermodynamic Analysis**: Calculate secondary structure stability using ViennaRNA
8. **Specificity Checking**: Validate specificity using BLAST
9. **Result Export**: Generate comprehensive Excel output with all design information

## Configuration

Customize the pipeline behavior with a JSON configuration file:

```bash
ddprimer --config config.json
```

Example configuration:

```json
{
  "DEBUG_MODE": true,
  "NUM_PROCESSES": 8,
  "PRIMER_MIN_SIZE": 19,
  "PRIMER_OPT_SIZE": 22,
  "PRIMER_MAX_SIZE": 25,
  "PRIMER_PRODUCT_SIZE_RANGE": [[90, 200]],
  "BLAST_WORD_SIZE": 7,
  "RESTRICTION_SITE": "GGCC",
  "THERMO_TEMPERATURE": 37,
  "THERMO_SODIUM": 0.05,
  "THERMO_MAGNESIUM": 0.002,
  "VCF_ALLELE_FREQUENCY_THRESHOLD": 0.05,
  "VCF_USE_SOFT_MASKING": false
}
```

## Additional Utilities

### Managing a BLAST Database

```bash
# Create a BLAST database from a genome file
ddprimer --db genome.fasta

# Create a BLAST database with a custom name
ddprimer --db genome.fasta my_custom_name

# Launch the model organism selection menu
ddprimer --db
```

### Configuration Management

```bash
# Display current configuration
ddprimer --config

# Display Primer3 settings
ddprimer --config primer3

# Generate configuration template
ddprimer --config template
```

## Output Format

The pipeline generates an Excel file with comprehensive information including:

- **Primer Sequences**: Forward, reverse, and probe sequences
- **Thermodynamic Properties**: Melting temperatures and ΔG values calculated with ViennaRNA
- **Amplicon Details**: Sequence, length, and GC content
- **Location Data**: Genomic coordinates of primers
- **Specificity Results**: Two best BLAST hits for all oligonucleotides

## File Preparation

ddPrimer automatically analyzes and prepares input files:

- **VCF Processing**: bgzip compression, tabix indexing, AF field addition, normalization
- **FASTA Indexing**: samtools faidx indexing for efficient access
- **Chromosome Mapping**: Intelligent chromosome name harmonization between files
- **GFF Processing**: Sorting and indexing for optimal performance

## SNP Processing

The pipeline intelligently handles variants:

- **Fixed Variants (AF=1.0)**: Substituted into reference sequence
- **Variable Variants (AF<1.0)**: Masked to avoid primer placement
- **Quality Filtering**: QUAL score thresholds for variant inclusion
- **Flanking Masking**: Optional masking around variant positions

## Troubleshooting

Common issues and solutions:

- **Missing BLAST database**: Run with `--db` to create or select a database
- **Memory errors**: Try reducing `NUM_PROCESSES` in your configuration file
- **GUI errors**: Use `--cli` to force command-line mode
- **macOS GUI issues**: Ensure pyobjc-core and pyobjc-framework-Cocoa are installed
- **VCF processing errors**: Verify bcftools is correctly installed
- **ViennaRNA installation issues**: 
  - Try installing via conda: `conda install -c bioconda viennarna`
  - If pip fails, install from source following ViennaRNA documentation
  - Ensure ViennaRNA is properly linked to your Python environment
- **File compatibility errors**: The pipeline will attempt automatic file preparation

For more help, run `ddprimer --help` or check the logs in `~/.ddPrimer/logs/`.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Development Installation

```bash
```bash
# Clone the repository
git clone https://github.com/jakobmueller/ddPrimer.git
cd ddPrimer

# Create and activate a conda environment
conda create -n ddprimer-dev python=3.8
conda activate ddprimer-dev

# Install external tools via conda
conda install -c bioconda -c conda-forge blast lastz bcftools samtools viennarna

# Install the package in development mode with ALL dependencies
pip install -e ".[dev]"
```

### Running Tests

```bash
# Run tests with coverage report
pytest --cov=ddprimer tests/
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.