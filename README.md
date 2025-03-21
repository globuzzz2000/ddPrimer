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
├── config.py           # Configuration settings
├── pipeline.py         # Main pipeline script
├── core/               # Core processing modules
│   ├── __init__.py
│   ├── annotation_processor.py
│   ├── blast_processor.py
│   ├── nupack_processor.py
│   ├── primer3_processor.py
│   ├── primer_processor.py
│   ├── sequence_processor.py
│   └── SNP_masking_processor.py
├── cross_species/      # Cross-species alignment modules
│   ├── __init__.py
│   ├── cross_species_workflow.py
│   ├── lastz_runner.py
│   └── maf_parser.py
└── utils/              # Utility functions
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

## Cross-Species Design Mode

The cross-species workflow uses LastZ for genome alignment and identifies conserved regions suitable for primer design. It maps variants from both species and ensures primers will work across multiple genomes.

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

Configuration settings can be provided in a YAML file:

```bash
ddprimer --config config.yaml
```

Example configuration:

```yaml
#############################################################################
#                           Pipeline Mode Options
#############################################################################
ONLY_LASTZ_MULTIZ: False           # Run only LastZ/MultiZ alignments
RUN_MATCH_CHECKER: False           # Run only the Match Checker
RUN_MATCH_CHECKER_ON_OUTPUT: False # Run Match Checker on pipeline output
DEBUG_MODE: False                  # Debug logging mode

#############################################################################
#                           Performance Settings
#############################################################################
NUM_PROCESSES: 8                   # Number of cores to use for parallel processing
BATCH_SIZE: 100                    # Batch size for parallel processing
MAF_CHUNK_SIZE: 10000              # Chunk size for MAF parsing
SHOW_PROGRESS: True                # Show progress bars

#############################################################################
#                           Design Parameters
#############################################################################
# Basic primer design constraints
PRIMER_MIN_SIZE: 18
PRIMER_OPT_SIZE: 20
PRIMER_MAX_SIZE: 23
PRIMER_MIN_TM: 50.0
PRIMER_OPT_TM: 57.5
PRIMER_MAX_TM: 65.0
PRIMER_MIN_GC: 50.0
PRIMER_MAX_GC: 60.0

# Product size constraints
PRIMER_PRODUCT_SIZE_RANGE: [[90, 200]]

# Pipeline parameters
MIN_SEGMENT_LENGTH: 100
RETAIN_TYPES: "gene"               # gff filtering: "gene", "mRNA", "CDS", "exon", etc.
FILTER_MEANINGFUL_NAMES: True      # only use named genes from gff
COUNT_AMBIGUOUS_AS_MISMATCH: False
GENE_OVERLAP_MARGIN: 25
RESTRICTION_SITE: "GGCC"
PENALTY_MAX: 5.0
MAX_PRIMER_PAIRS_PER_SEGMENT: 3
PREFER_PROBE_MORE_C_THAN_G: True   # For probe design

#############################################################################
#                           BLAST Database Options
#############################################################################
DB_FASTA: null                     # Path to a FASTA file to create a BLAST database from
DB_OUTPUT_DIR: null                # Custom output directory for the BLAST database
DB_NAME: null                      # Custom name for the BLAST database
USE_CUSTOM_DB: False               # Whether to use a custom database or the default
DB_PATH: "/path/to/blast/database" # Path to the existing BLAST database

# BLAST+ parameters
BLAST_WORD_SIZE: 7
BLAST_EVALUE: 10
BLAST_MAX_TARGET_SEQS: 100
BLAST_FILTER_FACTOR: 100           # E-value filtering factor

# NUPACK parameters
NUPACK_TEMPERATURE: 37             # Celsius
NUPACK_SODIUM: 0.05                # Molar
NUPACK_MAGNESIUM: 0.0              # Molar
```

The configuration also supports extensive Primer3 settings, which can be found in the full `config.py` file.

## Output

The pipeline generates an Excel file with comprehensive primer information:

- **Primers**: Forward and reverse primer sequences
- **Thermodynamics**: Tm, ΔG, and GC content
- **Amplicon**: Size, sequence, and properties
- **Genomic Location**: Chromosome and position
- **Quality Metrics**: Penalty scores and BLAST results
- **Probe Information**: Probe sequence and properties (if designed)
- **Cross-Species Data**: Mapped coordinates in second genome (if applicable)

### Output Columns Explained

- **Gene**: Target gene name (extracted from sequence ID or GFF)
- **Primer F/R**: Forward and reverse primer sequences
- **Tm F/R**: Melting temperature for each primer (°C)
- **Penalty F/R**: Primer3 penalty scores (lower is better)
- **Primer F/R dG**: Free energy (ΔG) for secondary structures
- **Primer F/R BLAST**: BLAST specificity results
- **Pair Penalty**: Combined penalty for the primer pair
- **Amplicon**: Full amplicon sequence
- **Length**: Amplicon length in base pairs
- **Amplicon GC%**: GC content percentage of the amplicon
- **Amplicon dG**: Free energy of the amplicon
- **Chromosome**: Source chromosome or contig
- **Location**: Genomic coordinates of the amplicon
- **Probe**: Probe sequence (if internal oligo was designed)
- **Probe Tm/Penalty/dG/BLAST**: Probe-specific metrics

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

### Custom Pipeline Modes

The pipeline supports several specialized execution modes:

```bash
# Run only the LastZ/MultiZ alignment step
ddprimer --fasta species1.fasta --second-fasta species2.fasta --config custom_config.yaml

# Run only the Match Checker to validate primers against specific genomes
ddprimer --config match_checker_config.yaml
```

In your configuration file, enable these special modes:

```yaml
# Pipeline mode options
ONLY_LASTZ_MULTIZ: True            # Run only alignment step
RUN_MATCH_CHECKER: False           # Run only validation step
RUN_MATCH_CHECKER_ON_OUTPUT: True  # Run validation on pipeline output
```

## Setting Up BLAST

For BLAST to work correctly, you need to create a BLAST database:

### Option 1: Use ddPrimer's Built-in Database Creation

```bash
# Create a BLAST database from a genome FASTA file
ddprimer --dbfasta /path/to/genome.fasta --dbname MyGenome
```

### Option 2: Manual BLAST Database Creation

```bash
# Install BLAST+
conda install -c bioconda blast

# Make a BLAST database
makeblastdb -in genome.fasta -dbtype nucl -out MyGenome

# Set the path in config.py or using --config
DB_PATH = "/path/to/MyGenome"
```

## Cross-Species Primer Design Requirements

For cross-species primer design, LastZ must be properly installed:

```bash
# Install LastZ with conda
# Install LastZ with conda
conda install -c bioconda lastz

# Verify installations
lastz --version
```

LastZ is used to generate Multiple Alignment Format (MAF) files that identify conserved regions between genomes, enabling the design of primers that work across multiple species.

## Logging and Debugging

Logs are stored in `~/.ddPrimer/logs/` with timestamps. Use `--debug` for verbose output:

```bash
ddprimer --debug
```

Log files contain detailed information about each step of the pipeline execution, including:
- Input file validation
- Variant extraction from VCF files
- Sequence processing steps
- Primer3 parameters and output
- BLAST results
- Error messages and warnings

This information is valuable for troubleshooting and ensuring that the pipeline is functioning as expected.

## Troubleshooting

If you encounter issues:

1. **Check Log Files**: Enable debug mode with `--debug` and examine the log file
2. **Input File Format**: Ensure FASTA, VCF, and GFF files are properly formatted
3. **Dependency Issues**: Verify all required tools are correctly installed
4. **Memory Issues**: For large genomes, try increasing available memory or processing in batches
5. **BLAST Database**: Confirm your BLAST database path is correctly set in the configuration
6. **Cross-Species Alignment**: For LastZ errors, check that input sequences have consistent chromosome/contig naming

Common errors and solutions:

| Error | Possible Solution |
|-------|-------------------|
| FASTA parsing error | Check for invalid characters or line endings in FASTA file |
| VCF parsing error | Ensure VCF file is properly formatted and bgzipped if compressed |
| No primers found | Try relaxing primer constraints in the configuration |
| LastZ alignment fails | Verify LastZ installation and input sequence validity |
| BLAST not working | Check BLAST+ installation and database path |

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citations

If you use ddPrimer in your research, please cite the following tools and resources that make this pipeline possible:

### Main Tools

- **Primer3**:  
  Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. (2012) Primer3 - new capabilities and interfaces. *Nucleic Acids Research* 40(15):e115  
  DOI: [10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596)

- **BLAST+**:  
  Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. (2009) BLAST+: architecture and applications. *BMC Bioinformatics* 10:421  
  DOI: [10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)

- **LastZ**:  
  Harris RS. (2007) Improved pairwise alignment of genomic DNA. Ph.D. Thesis, The Pennsylvania State University  
  [GitHub Repository](https://github.com/lastz/lastz)

- **PLastZ**:  
  Ho, A. (2020) PLastZ: Parallel version of the LastZ aligner.  
  [GitHub Repository](https://github.com/AntoineHo/PLastZ)

- **NUPACK**:  
  Dirks RM, Bois JS, Schaeffer JM, Winfree E, Pierce NA. (2007) Thermodynamic Analysis of Interacting Nucleic Acid Strands. *SIAM Review* 49(1):65-88  
  DOI: [10.1137/060651100](https://doi.org/10.1137/060651100)

### Python Libraries

- **Biopython**:  
  Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422-1423  
  DOI: [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

- **Pandas**:  
  McKinney W. (2010) Data Structures for Statistical Computing in Python. *Proceedings of the 9th Python in Science Conference* 51-56  
  DOI: [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)

### File Formats

- **VCF Format**:  
  Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group. (2011) The variant call format and VCFtools. *Bioinformatics* 27(15):2156-2158  
  DOI: [10.1093/bioinformatics/btr330](https://doi.org/10.1093/bioinformatics/btr330)

- **GFF Format**:  
  The Sequence Ontology Consortium. Generic Feature Format Version 3 (GFF3).  
  [Specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

- **MAF Format**:  
  UCSC Genome Browser. Multiple Alignment Format.  
  [Specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format5)

## Package Development

### Installation for Development

For developers who want to contribute to the project:

```bash
# Clone the repository
git clone https://github.com/yourusername/ddPrimer.git
cd ddPrimer

# Create and activate a conda environment
conda create -n ddprimer-dev python=3.8
conda activate ddprimer-dev

# Install in development mode
pip install -e ".[dev]"
```

### Running Tests

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run tests with coverage report
pytest --cov=ddPrimer tests/
```

### Building Documentation

```bash
# Install documentation dependencies
pip install sphinx sphinx-rtd-theme

# Build documentation
cd docs
make html
```