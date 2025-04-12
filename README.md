# ddPrimer: Droplet Digital PCR Primer Design Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive and streamlined pipeline for designing primers and probes specifically optimized for droplet digital PCR (ddPCR), with support for both single-species and cross-species primer design.

## Features

- **Complete Primer Design Pipeline**: From genome to validated primers in one workflow
- **Variant-Aware Design**: Masks variants from VCF files to avoid problematic regions
- **Cross-Species Design**: Designs primers that work across multiple species
- **Restriction Site Handling**: Cuts sequences at restriction sites for fragment analysis
- **Gene Overlap Filtering**: Prevents primers from spanning gene boundaries
- **Thermodynamic Analysis**: Calculates ΔG values using NUPACK
- **Specificity Checking**: BLAST validation for primer and probe specificity
- **Batch Processing**: Parallel execution for faster performance
- **Probe Design**: Optimized internal oligo design for TaqMan assays
- **Comprehensive Filtering**: Quality filters for penalties, repeats, GC content and more
- **Visualization**: Clear progress indicators during processing
- **Excel Export**: Comprehensive results in a formatted Excel file

## Project Structure

The ddPrimer package is organized as follows:

```
ddPrimer/
├── __init__.py
├── pipeline.py                # Main pipeline script
├── config/                    # Configuration settings
│   ├── __init__.py
│   ├── config.py              # Main Config class 
│   ├── exceptions.py          # Custom exceptions
│   └── logging_config.py      # Logging and progress reporting
├── core/                      # Core processing modules
│   ├── __init__.py
│   ├── annotation_processor.py
│   ├── blast_processor.py
│   ├── nupack_processor.py
│   ├── primer3_processor.py
│   ├── primer_processor.py
│   ├── sequence_processor.py
│   └── SNP_masking_processor.py
├── cross_species/             # Cross-species alignment modules
│   ├── __init__.py
│   ├── cross_species_workflow.py
│   ├── lastz_runner.py
│   └── maf_parser.py
├── modes/                     # Pipeline workflow modes
│   ├── __init__.py
│   ├── common.py              # Shared functionality
│   ├── direct_mode.py         # Direct sequence mode
│   ├── maf_mode.py            # Cross-species mode
│   └── standard_mode.py       # Standard mode
└── utils/                     # Utility functions
    ├── __init__.py
    ├── blast_db_creator.py
    ├── common_utils.py
    ├── file_utils.py
    └── sequence_utils.py

```

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/ddPrimer.git
cd ddPrimer

# Create and activate a conda environment
conda create -n ddprimer python=3.8
conda activate ddprimer

# Install the package
pip install -e .
```

### Dependencies

The following tools are required:

- **Python 3.7+**: [python.org](https://www.python.org/downloads/)
- **Biopython**: `pip install biopython`
- **Pandas**: `pip install pandas`
- **NUPACK 4.0+**: [nupack.org](http://www.nupack.org/)
  - Used for thermodynamic calculations
- **NCBI BLAST+**: [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  - Required for specificity checking
- **Primer3**: [Primer3](https://github.com/primer3-org/primer3)
  - Core primer design engine
- **LastZ**: [LastZ](https://github.com/lastz/lastz)
  - Used for cross-species genome alignment
- **PLastZ**: [PLastZ](https://github.com/AntoineHo/PLastZ)
  - Parallel implementation of LastZ for improved performance
- **tqdm**: `pip install tqdm`
  - For progress bars

#### Installing using Conda

Most dependencies can be installed using conda:

```bash
conda install -c bioconda -c conda-forge python=3.8 biopython pandas primer3 blast tqdm
conda install -c bioconda lastz
```

For NUPACK, follow installation instructions at [nupack.org](http://www.nupack.org/).

## Quick Start

### Command Line Usage

```bash
# Basic usage
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff --output results/

# Cross-species primer design
ddprimer --alignment --fasta species1.fasta --vcf species1.vcf \
         --second-fasta species2.fasta --second-vcf species2.vcf \
         --gff annotations.gff
```

### Interactive Mode

```bash
# Launch in interactive mode (will prompt for file selection)
ddprimer
```

## Workflow Overview

1. **Inputs**: FASTA genome(s), VCF file(s) with variants, and GFF annotation file
2. **Variant Masking**: Identifies and masks all variant positions
3. **Cross-Species Alignment** (optional): Aligns genomes and identifies conserved regions
4. **Restriction Site Processing**: Cuts sequences at restriction enzyme recognition sites
5. **Gene Boundary Filtering**: Removes sequences that span gene boundaries
6. **Primer3 Design**: Designs primer and probe candidates
7. **Filtering**: Applies quality filters (penalties, repeats, GC content)
8. **Thermodynamic Analysis**: Calculates secondary structure stability
9. **BLAST Specificity**: Checks for off-target binding
10. **Results**: Exports all data in a formatted Excel file

## Pipeline Modes

The pipeline supports three main operational modes:

### Standard Mode

Standard single-species primer design using FASTA, VCF, and GFF files.

```bash
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff
```

### Direct Mode

Design primers from sequences directly (without coordinate information), using Excel or CSV input.

```bash
ddprimer --direct sequences.csv
```

### Alignmnet Mode

Design primers using genome alignment.

```bash
# Using pre-computed MAF alignment file
ddprimer --alignment --maf-file alignment.maf --vcf species1.vcf \
         --second-vcf species2.vcf --gff annotations.gff

# Generating alignment automatically
ddprimer --alignment --fasta species1.fasta --vcf species1.vcf \
         --second-fasta species2.fasta --second-vcf species2.vcf \
         --gff annotations.gff --min-identity 85
```

## Configuration

Configuration settings can be provided in a JSON or Primer3 format:

```bash
# Using JSON config
ddprimer --config config.json

# Using Primer3 style config
ddprimer --config primer3_settings.txt
```

Example JSON configuration:

```json
{
  "DEBUG_MODE": true,
  "NUM_PROCESSES": 8,
  "PRIMER_MIN_SIZE": 19,
  "PRIMER_OPT_SIZE": 22,
  "PRIMER_MAX_SIZE": 25,
  "PRIMER_PRODUCT_SIZE_RANGE": [[90, 200]],
  "BLAST_WORD_SIZE": 7,
  "RESTRICTION_SITE": "GGCC"
}
```

## Advanced Usage

### Creating a BLAST Database

```bash
ddprimer --dbfasta genome.fasta --dbname my_genome
```

### Customizing LastZ Alignment

```bash
ddprimer --alignment --fasta species1.fasta --second-fasta species2.fasta \
         --lastz-options "--format=maf --ambiguous=iupac --chain"
```

## Logging and Debugging

Logs are stored in `~/.ddPrimer/logs/` with timestamps. Use `--debug` for verbose output:

```bash
ddprimer --debug
```

## Output Format

The pipeline generates an Excel file with comprehensive primer information including:

- Forward and reverse primer sequences with thermodynamic properties
- Amplicon details (sequence, length, GC content)
- Probe information when internal oligos are designed
- BLAST specificity results
- Genomic coordinates of primers
- Cross-species mapping information (for MAF mode)

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

# Install in development mode with testing dependencies
pip install -e ".[dev]"
```

### Running Tests

```bash
# Run tests with coverage report
pytest --cov=ddprimer tests/
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.