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

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/ddPrimer.git
cd ddPrimer

# Create and activate a conda environment
conda create -n ddprimer python=3.8
conda activate ddprimer

# Install the package and dependencies
pip install -e .
```

### Dependencies

The following tools are required:

- **Python 3.7+**: For core functionality
- **Primer3**: For primer design (core engine)
- **NUPACK 4.0+**: For thermodynamic calculations
- **NCBI BLAST+**: For specificity checking
- **LastZ**: For alignment-based mode

Most dependencies can be installed through conda:

```bash
conda install -c bioconda -c conda-forge python=3.8 biopython pandas primer3 blast tqdm
conda install -c bioconda lastz
```

For NUPACK, follow installation instructions at [nupack.org](http://www.nupack.org/).

## Quick Start

### Command Line Usage

```bash
# Direct sequence mode
ddprimer --direct sequences.csv

# Direct sequence mode with SNP filtering
ddprimer --direct sequences.csv --snp --fasta genome.fasta --vcf variants.vcf
```

The input file should contain:
- Column 1: Sequence IDs or names
- Column 2: DNA sequences

```bash
# Basic primer design
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff

# Basic primer design without gene filtering
ddprimer --noannotation --fasta genome.fasta --vcf variants.vcf

# Alignment based design
ddprimer --alignment --fasta species1.fasta --second-fasta species2.fasta --gff annotations.gff

# Alignment based design with SNP filtering 
ddprimer --alignment --snp --fasta species1.fasta --second-fasta species2.fasta \
         --vcf species1.vcf --second-vcf species2.vcf
```

### Interactive Mode

Simply run `ddprimer` without arguments to launch the interactive mode, which will guide you through file selection with a graphical interface.

## Workflow Overview

1. **Input Selection**: Choose genome FASTA, variant VCF, and annotation GFF files
2. **Variant Masking**: Identify and mask all variant positions in the genome
3. **Alignment** (Alignmnet mode only): Align genomes and identify conserved regions
4. **Sequence Preparation**: Filter sequences based on restriction sites and gene boundaries
5. **Primer Design**: Design primer and probe candidates using Primer3
6. **Quality Filtering**: Apply filters for penalties, repeats, GC content, and more
7. **Thermodynamic Analysis**: Calculate secondary structure stability via NUPACK
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
  "RESTRICTION_SITE": "GGCC"
}
```

## Additional Utilities

### Creating a BLAST Database

```bash
ddprimer --createdb genome.fasta --dbname my_genome
```

### Running Only Alignment (without primer design)

```bash
ddprimer --alignment --lastzonly --fasta species1.fasta --second-fasta species2.fasta
```


## Output Format

The pipeline generates an Excel file with comprehensive information including:

- **Primer Sequences**: Forward, reverse, and probe sequences
- **Thermodynamic Properties**: Melting temperatures and ΔG values
- **Amplicon Details**: Sequence, length, and GC content
- **Location Data**: Genomic coordinates of primers
- **Specificity Results**: Two best BLAST hits for all oligonucleotides
- **Cross-Species Mapping**: Coordinates in both reference genomes (for alignment mode)

## Troubleshooting

Common issues and solutions:

- **Missing BLAST database**: Run with `--createdb` to create a new database
- **Memory errors**: Try reducing `NUM_PROCESSES` in your configuration file
- **GUI errors**: Use `--cli` to force command-line mode

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
