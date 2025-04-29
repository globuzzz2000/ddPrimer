# ddPrimer: Advanced Droplet Digital PCR Primer Design

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)

A comprehensive pipeline for designing primers and probes specifically optimized for droplet digital PCR (ddPCR), with support for both single-species and cross-species design workflows.

## Key Features

- **Complete End-to-End Pipeline**: Design primers from genome sequences through a streamlined workflow
- **Smart SNP Masking**: Avoid designing primers across variant positions using VCF files
- **Cross-Species Design**: Create primers that work consistently across multiple genomes 
- **Thermodynamic Optimization**: Calculate ΔG values to prevent unwanted secondary structures
- **Specificity Verification**: Integrated BLAST validation for both primers and probes
- **Multiple Workflow Modes**:
  - **Standard Mode**: Design from genome and annotation files
  - **Direct Mode**: Design from raw sequences in CSV/Excel files
  - **Alignment Mode**: Design from genome alignments
- **User-Friendly Interface**: GUI file selection with clear progress indicators
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
- **NCBI BLAST+**: For specificity checking
- **NUPACK 4.0+**: For thermodynamic calculations
- **LastZ/PLastZ**: For cross-species alignments

Most dependencies can be installed through conda:

```bash
conda install -c bioconda -c conda-forge python=3.8 biopython pandas primer3 blast tqdm
conda install -c bioconda lastz
```

For NUPACK, follow installation instructions at [nupack.org](http://www.nupack.org/).

## Quick Start

### Command Line Usage

```bash
# Basic primer design
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff

# Cross-species design
ddprimer --alignment --fasta species1.fasta --second-fasta species2.fasta \
         --vcf species1.vcf --second-vcf species2.vcf

# Direct sequence mode
ddprimer --direct sequences.csv
```

### Interactive Mode

Simply run `ddprimer` without arguments to launch the interactive mode, which will guide you through file selection with a graphical interface.

## Workflow Overview

1. **Input Selection**: Choose genome FASTA, variant VCF, and annotation GFF files
2. **Variant Masking**: Identify and mask all variant positions in the genome
3. **Cross-Species Alignment** (optional): Align genomes and identify conserved regions
4. **Sequence Preparation**: Filter sequences based on restriction sites and gene boundaries
5. **Primer Design**: Design primer and probe candidates using Primer3
6. **Quality Filtering**: Apply filters for penalties, repeats, GC content, and more
7. **Thermodynamic Analysis**: Calculate secondary structure stability via NUPACK
8. **Specificity Checking**: Validate specificity using BLAST
9. **Result Export**: Generate comprehensive Excel output with all design information

## Pipeline Modes

### Standard Mode

Design primers from a single genome:

```bash
ddprimer --fasta genome.fasta --vcf variants.vcf --gff annotations.gff --output results/
```

### Direct Mode

Design primers directly from sequences in CSV or Excel format:

```bash
ddprimer --direct sequences.xlsx
```

The input file should contain:
- Column 1: Sequence IDs or names
- Column 2: DNA sequences

### Alignment Mode

Design primers that work across multiple species:

```bash
# Using two genomes (automatically computes alignment)
ddprimer --alignment --fasta species1.fasta --second-fasta species2.fasta \
         --vcf species1.vcf --second-vcf species2.vcf

# Using a pre-computed MAF alignment file
ddprimer --alignment --maf-file alignment.maf --vcf species1.vcf --second-vcf species2.vcf
```

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

## Advanced Usage

### Creating a BLAST Database

```bash
ddprimer --createdb genome.fasta --dbname my_genome
```

### Running Only Alignment (without primer design)

```bash
ddprimer --alignment --lastzonly --fasta species1.fasta --second-fasta species2.fasta
```

### Customizing LastZ Alignment

```bash
ddprimer --alignment --fasta species1.fasta --second-fasta species2.fasta \
         --lastz-options "--format=maf --ambiguous=iupac --chain"
```

## Output Format

The pipeline generates an Excel file with comprehensive information including:

- **Primer Sequences**: Forward, reverse, and probe sequences
- **Thermodynamic Properties**: Melting temperatures and ΔG values
- **Amplicon Details**: Sequence, length, and GC content
- **Location Data**: Genomic coordinates of primers
- **Specificity Results**: BLAST hits for all oligonucleotides
- **Cross-Species Mapping**: Coordinates in both reference genomes (for alignment mode)

## Troubleshooting

Common issues and solutions:

- **Missing BLAST database**: Run with `--createdb` to create a new database
- **Memory errors**: Try reducing `NUM_PROCESSES` in your configuration file
- **GUI errors**: Use `--cli` to force command-line mode
- **Alignment failures**: Check input sequence quality and try adjusting `--min-identity`

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
