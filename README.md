# ddPrimer: Advanced Droplet Digital PCR Primer Design

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)

A comprehensive pipeline for designing primers and probes specifically optimized for droplet digital PCR (ddPCR).

## Key Features

- **Complete End-to-End Pipeline**: Design primers from genome sequences through a streamlined workflow using Primer3
- **Smart SNP Masking**: Avoid designing primers across variant positions using VCF files
- **Thermodynamic Optimization**: Calculate ΔG values to prevent unwanted secondary structures
- **Specificity Verification**: Integrated BLAST validation for both primers and probes
- **Multiple Workflow Modes**:
  - **Standard Mode**: Design from genome, variant and annotation files
  - **Direct Mode**: Design from desired sequences in CSV/Excel files
  - **Alignment Mode**: Design from genome alignments
- **Comprehensive Results**: Detailed Excel output with key metrics and sequence information

## Project Structure

```
ddPrimer/
├── ddprimer/                      # Main package directory
│   ├── __init__.py
│   ├── pipeline.py                # Main entry point for the pipeline
│   ├── modes/                     # Pipeline operation modes
│   │   ├── __init__.py
│   │   ├── alignment_mode.py      # Alignment-based design
│   │   ├── common.py              # Shared workflow components
│   │   ├── direct_mode.py         # Direct sequence design
│   │   └── standard_mode.py       # Standard pipeline
│   ├── core/                      # Core processing modules
│   │   ├── __init__.py
│   │   ├── SNP_masking_processor.py
│   │   ├── annotation_processor.py
│   │   ├── blast_processor.py
│   │   ├── nupack_processor.py
│   │   ├── vienna_processor.py    # ViennaRNA thermodynamic processor
│   │   ├── thermo_processor.py    # Flexible thermodynamics processor
│   │   ├── primer3_processor.py
│   │   ├── primer_processor.py
│   │   └── sequence_processor.py
│   ├── helpers/                   # Helper functions and classes
│   │   ├── __init__.py
│   │   ├── alignment_workflow.py
│   │   ├── direct_masking.py
│   │   ├── lastz_runner.py
│   │   └── maf_parser.py
│   ├── config/                    # Configuration and settings
│   │   ├── __init__.py
│   │   ├── config.py              # Core configuration settings
│   │   ├── config_display.py      # Configuration display utilities
│   │   ├── exceptions.py          # Custom exceptions
│   │   ├── logging_config.py      # Logging setup
│   │   └── template_generator.py  # Configuration template generation
│   ├── tests/                     # Test data and fixtures
│   │   ├── Alignment.maf
│   │   ├── Annotations.gff
│   │   ├── Sequences.xlsx
│   │   ├── VCF1.vcf.gz
│   │   ├── VCF2.vcf.gz
│   │   ├── fasta1.fna
│   │   └── fasta2.fna
│   └── utils/                     # Utility functions
│       ├── __init__.py
│       ├── blast_db_creator.py
│       ├── blast_verification.py
│       ├── common_utils.py
│       ├── file_io.py
│       └── sequence_utils.py
├── pyproject.toml                 # Package configuration and dependencies
└── README.md                      # This file
```

## Installation

### Using Conda (Recommended)

The easiest way to install ddPrimer with all dependencies is using conda:

```bash
# Clone the repository
git clone https://github.com/globuzzz2000/ddPrimer
cd ddPrimer

# Create and activate a conda environment
conda create -n ddprimer python=3.8
conda activate ddprimer

# Install external tools via conda
conda install -c bioconda -c conda-forge blast lastz

# Install GUI dependency (if needed separately)
conda install -c conda-forge wxpython

# Install the package with all Python dependencies
pip install -e ".[thermo]"
```

### Manual Installation

If you prefer not to use conda, you can install the package and dependencies manually:

```bash
# Clone the repository
git clone https://github.com/globuzzz2000/ddPrimer
cd ddPrimer

# Install the package
pip install -e .

# Install BLAST and LastZ externally using a package manager
# On macOS:
brew install blast lastz primer3
# On Linux:
sudo apt-get install ncbi-blast+ lastz
```

### Dependencies

The following tools are required:

- **Python 3.7+**: For core functionality
- **Primer3**: For primer design (core engine)
- **NCBI BLAST+**: For specificity checking
- **LastZ**: For alignment-based mode

For thermodynamic calculations (one of the following):
- **ViennaRNA**: Default option for thermodynamic calculations
- **NUPACK 4.0+**: Alternative option for thermodynamic calculations with advanced features

### Installing Thermodynamics Engines

ddPrimer supports two thermodynamics engines:

#### ViennaRNA (Default)

ViennaRNA is the default thermodynamics engine and can be installed via pip:

```bash
pip install viennarna
```

#### NUPACK (Alternative)

NUPACK is an alternative thermodynamics engine with more advanced features:

```bash
pip install nupack
```

If the pip installation of NUPACK fails, you may need to install it manually. Follow the installation instructions at [nupack.org](http://www.nupack.org/).


To use NUPACK instead of ViennaRNA, set in your configuration:
```json
{
  "THERMO_BACKEND": "nupack"
}
```

## Quick Start

### Command Line Usage

```bash
# Direct sequence mode
ddprimer --direct sequences.csv

# Direct sequence mode with SNP filtering
ddprimer --direct sequences.csv --snp --fasta genome.fasta --vcf variants.vcf
```

The input csv file should contain:
- Column 1: Sequence IDs or names
- Column 2: DNA sequences

```bash
# Basic primer design
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff

# Basic primer design without gene filtering
ddprimer --noannotation --fasta genome.fasta --vcf variants.vcf

# Alignment based design
ddprimer --alignment --fasta genome1.fasta --second-fasta genome2.fasta --gff annotations.gff

# Alignment based design with SNP filtering 
ddprimer --alignment --snp --fasta genome1.fasta --second-fasta genome2.fasta \
         --vcf genome1.vcf --second-vcf genome2.vcf
```

### Interactive Mode

Simply run `ddprimer` without arguments to launch the interactive mode, which will guide you through file selection with a graphical interface.

## Workflow Overview

1. **Input Selection**: Choose genome FASTA, variant VCF, and annotation GFF files
2. **Variant Masking**: Identify and mask all variant positions in the genome
3. **Alignment** (Alignment mode only): Align genomes and identify conserved regions
4. **Sequence Preparation**: Filter sequences based on restriction sites and gene boundaries
5. **Primer Design**: Design primer and probe candidates using Primer3
6. **Quality Filtering**: Apply filters for penalties, repeats, GC content, and more
7. **Thermodynamic Analysis**: Calculate secondary structure stability via ViennaRNA or NUPACK
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
  "THERMO_BACKEND": "vienna"
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

### Running Only Alignment (without primer design)

```bash
ddprimer --alignment --lastzonly --fasta genome1.fasta --second-fasta genome2.fasta
```


## Output Format

The pipeline generates an Excel file with comprehensive information including:

- **Primer Sequences**: Forward, reverse, and probe sequences
- **Thermodynamic Properties**: Melting temperatures and ΔG values
- **Amplicon Details**: Sequence, length, and GC content
- **Location Data**: Genomic coordinates of primers
- **Specificity Results**: Two best BLAST hits for all oligonucleotides
- **Alignment Mapping**: Coordinates in both reference genomes (for alignment mode)

## Troubleshooting

Common issues and solutions:

- **Missing BLAST database**: Run with `--createdb` to create a new database
- **Memory errors**: Try reducing `NUM_PROCESSES` in your configuration file
- **GUI errors**: Use `--cli` to force command-line mode
- **Thermodynamic calculation issues with NUPACK**: Use ViennaRNA instead by setting `THERMO_BACKEND: "vienna"` in your config file (this is now the default)

For more help, run `ddprimer --help` or check the logs in `~/.ddPrimer/logs/`.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Development Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/ddPrimer.git
cd ddPrimer

# Create and activate a conda environment
conda create -n ddprimer-dev python=3.8
conda activate ddprimer-dev

# Install external tools via conda
conda install -c bioconda -c conda-forge blast lastz

# Install the package in development mode with ALL dependencies
# This will automatically install all Python dependencies defined in pyproject.toml
pip install -e ".[dev,thermo]"
```

### Running Tests

```bash
# Run tests with coverage report
pytest --cov=ddprimer tests/
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.