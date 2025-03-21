# ddPrimer Pipeline

A comprehensive pipeline for primer design and filtering that integrates multiple tools and techniques for optimal primer selection.

## Overview

The ddPrimer pipeline automates the process of designing and filtering PCR primers based on various criteria. It takes into account genetic variants, restriction sites, gene annotations, and performs thorough specificity checks to ensure high-quality, specific primers.

## Features

- Variant-aware primer design (avoids SNPs and other variants)
- Restriction site detection and filtering
- Gene annotation integration
- Primer filtering based on:
  - Penalty scores
  - GC content
  - Repeated sequences (e.g., GGGG, CCCC)
  - BLAST specificity
- Thermodynamic property calculation
- Probe/internal oligo processing and optimization
- Multiple output formats

## Package Structure

```
ddprimer/
├── __init__.py                # Package initialization, version info
├── pipeline.py                # Main pipeline script
├── config.py                  # Configuration settings
├── core/                      # Core processing modules
│   ├── __init__.py            # Core subpackage initialization
│   ├── sequence_processor.py  # Sequence filtering and processing
│   ├── primer_processor.py    # Primer filtering and processing
│   ├── annotation_processor.py # GFF file processing
│   ├── primer3_processor.py   # Primer3 interface
│   ├── blast_processor.py     # BLAST interface
│   ├── nupack_processor.py    # NUPACK interface
│   └── vcf_processor.py       # VCF processing and SNP masking
└── utils/                     # Utility modules
    ├── __init__.py            # Utils subpackage initialization
    ├── sequence_utils.py      # Sequence-specific utilities
    ├── file_utils.py          # File I/O utilities
    ├── ui_utils.py            # UI-related utilities
    └── common_utils.py        # General utilities
```

## Module Descriptions

### Core Modules

- **sequence_processor.py**: Handles sequence filtering by restriction sites and gene overlap. Also includes the `MatchProcessor` class for checking sequence matches against genomes.

- **primer_processor.py**: Contains filtering methods for primers based on various criteria like penalty scores, repeats, GC content, and BLAST specificity.

- **annotation_processor.py**: Processes GFF annotation files to extract gene information, with optimized parallel processing.

- **primer3_processor.py**: Interface to Primer3 for designing primers on target sequences.

- **blast_processor.py**: Runs BLAST for checking primer specificity against reference genomes.

- **nupack_processor.py**: Calculates thermodynamic properties like minimum free energy using NUPACK.

- **vcf_processor.py**: Extracts variants from VCF files and masks sequences at variant positions.

### Utility Modules

- **sequence_utils.py**: Contains functions for working with DNA/RNA sequences, such as GC content calculation and repeat detection.

- **file_utils.py**: Handles file operations including UI dialogs and file parsing, with fallback to CLI mode.

- **ui_utils.py**: UI-related utility functions like the progress spinner.

- **common_utils.py**: General utility functions used across the pipeline, like list chunking.

## Pipeline Workflow

1. **Input Selection**: The user provides FASTA, VCF, and GFF files (either via command line or graphical file picker).

2. **Variant Masking**: Variants from the VCF file are extracted and the corresponding positions in sequences are masked.

3. **Sequence Filtering**: Sequences are cut at restriction sites and filtered based on gene overlap.

4. **Primer Design**: Primer3 is used to design primers on the filtered sequences.

5. **Primer Filtering**: Primers are filtered based on penalty, repeats, and GC content.

6. **Thermodynamic Analysis**: NUPACK is used to calculate thermodynamic properties of primers and amplicons.

7. **Specificity Checking**: BLAST is used to check primer specificity against reference genomes.

8. **Results Output**: Final primer pairs are saved to an Excel file.

## Usage

### Command Line Interface

```bash
python -m ddprimer.pipeline --fasta genome.fa --vcf variants.vcf --gff annotations.gff
```

### Optional Arguments

- `--fasta`: Input FASTA file (if not provided, will prompt)
- `--vcf`: VCF file with variants (if not provided, will prompt)
- `--gff`: GFF annotation file (if not provided, will prompt)
- `--output`: Output directory (defaults to "primers" subdirectory)
- `--config`: Configuration file
- `--cli`: Force CLI mode (no GUI file dialogs)

### As a Python Package

```python
from ddprimer.pipeline import run_pipeline

# Run with default settings (will prompt for files)
output_file = run_pipeline()

# Or specify files programmatically
import sys
sys.argv = ['pipeline.py', '--fasta', 'genome.fa', '--vcf', 'variants.vcf', '--gff', 'annotations.gff']
output_file = run_pipeline()
```

## Dependencies

- Python 3.6+
- pandas
- BioPython
- tkinter (for GUI file dialogs)
- primer3-py
- tqdm (for progress bars)

## Installation

```bash
# From source
git clone https://github.com/yourusername/ddprimer.git
cd ddprimer
pip install -e .

# Or directly from PyPI (if published)
pip install ddprimer
```
