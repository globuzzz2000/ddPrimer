# ddPrimer: Advanced Droplet Digital PCR Primer Design

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)

A comprehensive pipeline for designing primers and probes specifically optimized for droplet digital PCR (ddPCR).


## Key Features

- **Complete End-to-End Pipeline**: Design primers from genome sequences through a streamlined workflow using Primer3
- **ddPCR-specific utilities**: Restriction cutting, GC% filtering, and optimized design parameters for ddPCR applications
- **Gene annotation filtering**: Filter fragments based on gene overlap using GFF files for targeted design
- **SNP Masking**: Avoid designing primers across variant positions using VCF files with AF-based processing
- **Thermodynamic Optimization**: Calculate ΔG values using ViennaRNA to prevent unwanted secondary structures
- **Specificity Verification**: Integrated BLAST validation for both primers and probes
- **File Preparation**: Automatic VCF normalization, chromosome mapping, and file indexing
- **Output format**: Results saved as Excel file with primer sequences, coordinates, thermodynamics, and quality metrics


## Installation

### Quick Install with Conda

```bash
# Clone and install
git clone https://github.com/globuzzz2000/ddPrimer
cd ddPrimer

# Create environment with all dependencies
conda create -n ddprimer python=3.8
conda activate ddprimer
conda install -c bioconda -c conda-forge blast bcftools samtools
pip install -e .

# Alternatively install external tools via system package manager:
# macOS: brew install blast bcftools samtools
# Linux: sudo apt-get install ncbi-blast+ bcftools samtools
```

### Required External Tools

- **NCBI BLAST+**: For specificity checking
- **bcftools and samtools**: For file processing  

Python dependencies are automatically installed via pip.


## Quick Start

### Command Line Usage

```bash
# Basic primer design with file preparation
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff

# Basic primer design without annotation filtering
ddprimer --noannotation --fasta genome.fasta --vcf variants.vcf
```

### Interactive Mode

Simply run `ddprimer` without arguments to launch the interactive mode, which will guide you through file selection with a graphical interface.


## Project Structure

```
ddPrimer/
├── ddprimer/                         # Main package directory
│   ├── __init__.py
│   ├── main.py                       # Main entry point for the pipeline
│   ├── core/                         # Core processing modules
│   │   ├── __init__.py
│   │   ├── annotation_processor.py   # GFF-based annotation filtering
│   │   ├── sequence_processor.py     # Sequence processing
│   │   ├── snp_processor.py          # VCF-based variant masking
│   │   ├── primer3_processor.py      # Primer3-based primer design
│   │   ├── primer_processor.py       # Primer parsing
│   │   ├── filter_processor.py       # Primer quality filtering
│   │   ├── thermo_processor.py       # ViennaRNA thermodynamic processor
│   │   └── blast_processor.py        # Primer specificity filtering
│   ├── utils/                        # Utility functions
│   │   ├── __init__.py
│   │   ├── file_preparator.py        # File validation and preparation
│   │   ├── file_io.py                # File I/O and Excel formatting
│   │   ├── blast_db_manager.py       # Unified BLAST database management
│   │   ├── direct_mode.py            # Target List-based Primer design
│   │   └── primer_remapper           # Primer coordinate remapping to different genome
│   └── config/                       # Configuration and settings
│       ├── __init__.py
│       ├── config.py                 # Core configuration settings
│       ├── config_display.py         # Configuration display
│       ├── exceptions.py             # Error handling
│       ├── logging_config.py         # Logging setup
│       └── template_generator.py     # Configuration template generation
├── pyproject.toml                    # Package configuration and dependencies
└── README.md                  
```


## Workflow Overview

1. **Input Selection**: Choose genome FASTA, variant VCF, and annotation GFF files
2. **File Preparation**: Validate and prepare files (bgzip compression, indexing, normalization, chromosome mapping)
3. **Sequence Preparation**: Filter sequences based on restriction sites and gene boundaries
4. **Variant Processing**: Apply VCF variants to sequences with intelligent AF-based masking/substitution
5. **Primer Design**: Design primer and probe candidates using Primer3
6. **Quality Filtering**: Apply filters for penalties, repeats, GC content, and more
7. **Thermodynamic Analysis**: Calculate secondary structure stability using ViennaRNA
8. **Specificity Checking**: Validate specificity using BLAST
9. **Result Export**: Generate comprehensive Excel output


## Configuration

Customize the pipeline behavior with a JSON configuration file:

```bash
ddprimer --config config.json
```

Example configuration:

```json
{
  "NUM_PROCESSES": 6,
  "SHOW_PROGRESS": true,
  "PRIMER_MIN_SIZE": 18,
  "PRIMER_OPT_SIZE": 20,
  "PRIMER_MAX_SIZE": 23,
  "PRIMER_MIN_GC": 50.0,
  "PRIMER_MAX_GC": 60.0,
  "MIN_SEGMENT_LENGTH": 90,
  "RETAIN_TYPES": "['gene']",
  "RESTRICTION_SITE": "GGCC",
  "PENALTY_MAX": 5.0,
}
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


## Additional Utilities

### Direct Mode
Target-sequence based primer design using CSV/Excel input, bypassing genome-based processing:

**CSV/Excel format**: Two-column table with sequence ID and DNA sequence

```bash
# Basic direct mode with sequence table
ddprimer --direct sequences.csv

# Direct mode with SNP masking (Target sequences should exactly match reference sequence)
ddprimer --direct sequences.xlsx --snp --vcf variants.vcf --fasta reference.fa

```

### BLAST Database Management
Create and manage BLAST databases for primer specificity checking:

```bash
# Create BLAST database from FASTA file
ddprimer --db genome.fasta

# Select from existing databases model organisms (E. coli, S. cerevisiae, etc.)
ddprimer --db

# Use custom database name
ddprimer --db genome.fasta my_custom_db
```

### Primer Remapping
Update existing ddprimer output coordinates and annotations against a new reference genome:

```bash
# With gene annotation updates
ddprimer --remap primers.csv --fasta ref.fa --gff annotations.gff

# Skip annotation updates
ddprimer --remap primers.xlsx --fasta ref.fa --noannotation
```


## Troubleshooting

Common issues and solutions:

- **Missing BLAST database**: Run with `--db` to create or select a database
- **GUI errors**: Use `--cli` to force command-line mode
- **File compatibility errors**: The pipeline will attempt automatic file preparation

For more detailed output, run `ddprimer --debug` or check the logs in `~/.ddPrimer/logs/`.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.