# ddPrimer

A comprehensive, object-oriented pipeline for designing specific primers across genome comparisons.

## Overview

This tool automates a complete primer design pipeline, integrating PLastZ genome alignment, GFF gene annotation parsing, sequence filtering, Primer3 primer design, BLAST specificity filtering, and NUPACK thermodynamic evaluation.

## Features

- **Flexible Genome Comparison**: Compare two genomes with PLastZ or multiple genomes with MultiZ
- **Gene-focused Design**: Targets primers within or near annotated genes
- **Restriction Site Filtering**: Excludes sequences containing specified restriction sites
- **Specificity Validation**: Uses BLAST to ensure primer uniqueness
- **Thermodynamic Analysis**: Calculates structure and stability using NUPACK
- **Comprehensive Processing**: Automates the entire primer design workflow
- **Comprehensive Output**: Detailed Excel report with all primer properties
- **Optimized Parallel Processing**: Uses PLastZ for efficient parallel genome comparisons

## System Requirements

- Python 3.7+
- 8GB+ RAM recommended (16GB+ for large genomes)
- 10GB+ free disk space
- External dependencies (see Setup section)

## Installation

### 1. Clone/Download the Package

Download the ddPrimer package which includes:
- `ddPrimer.py` - Main script
- `Setup.py` - Setup helper
- `PLastZ/` - Directory containing PLastZ scripts:
  - `PLastZ.py` - The parallel processing script for LastZ

### 2. Run the Setup Script

The Setup.py script will guide you through the installation of dependencies:

```bash
python Setup.py
```

For systems without a graphical interface (e.g., remote servers, SSH sessions):

```bash
python Setup.py --headless
```

The setup script will:
- Install required Python packages
- Help you install and configure external tools
- Guide you through BLAST database setup
- Compile MultiZ if needed
- Configure paths for your system
- Set up the PLastZ implementation

### 3. Manual Installation of Dependencies

If the setup script doesn't work for your system, follow these manual steps:

#### Python Dependencies

```bash
pip install pandas tqdm nupack
```

#### LastZ

LastZ is used for genome-to-genome alignments.

**macOS (with Homebrew):**
```bash
brew install lastz
```

**Linux:**
```bash
# Debian/Ubuntu
apt-get install lastz

# From source
wget https://github.com/lastz/lastz/archive/refs/tags/1.04.15.tar.gz
tar -xzf 1.04.15.tar.gz
cd lastz-1.04.15
make
sudo make install
```

**Windows:**
Download from [LastZ GitHub](https://github.com/lastz/lastz/releases) and add to PATH.

#### Samtools

Samtools is required for PLastZ functionality.

**macOS (with Homebrew):**
```bash
brew install samtools
```

**Linux:**
```bash
# Debian/Ubuntu
apt-get install samtools

# From source
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
tar -xjf samtools-1.17.tar.bz2
cd samtools-1.17
./configure
make
sudo make install
```

**Windows:**
Download from [Samtools](http://www.htslib.org/download/) and add to PATH.

#### MultiZ

MultiZ is used for multiple sequence alignment when comparing more than two genomes.

**macOS/Linux:**
```bash
# From source
git clone https://github.com/multiz/multiz.git
cd multiz
make
sudo cp multiz /usr/local/bin/
```

**Windows:**
Compile from source or use WSL to run the Linux version.

#### Primer3

Used for designing PCR primers based on specified parameters.

**macOS (with Homebrew):**
```bash
brew install primer3
```

**Linux:**
```bash
# Debian/Ubuntu
apt-get install primer3

# From source
wget https://github.com/primer3-org/primer3/archive/refs/tags/v2.6.1.tar.gz
tar -xzf v2.6.1.tar.gz
cd primer3-2.6.1/src
make
sudo make install
```

**Windows:**
Download from [Primer3 Releases](https://github.com/primer3-org/primer3/releases) and add to PATH.

#### BLAST+

Used for primer specificity checking.

**macOS (with Homebrew):**
```bash
brew install blast
```

**Linux:**
```bash
# Debian/Ubuntu
apt-get install ncbi-blast+

# From NCBI
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.13.0+-x64-linux.tar.gz
# Add to PATH or move binaries to /usr/local/bin
```

**Windows:**
Download installer from [NCBI BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and install.

#### NUPACK

The script uses the Python package `nupack` which can be installed via pip:

```bash
pip install nupack
```

Note: On some systems, you may need additional compilers for NUPACK installation.

### 4. Configure BLAST Database

You need to create a BLAST database from your reference genome:

```bash
makeblastdb -in reference_genome.fasta -dbtype nucl -out my_reference_db
```

Then update the `DB_PATH` variable in the script to point to your database.

## Configuration

Before running, adjust parameters in the script to match your system and requirements:

```python
# Tool paths
PLASTZ_PATH = "/Library/Application Support/ddPrimer/PLastZ/PLastZ.py"
PRIMER3_CORE_PATH = "/opt/homebrew/bin/primer3_core"
MULTIZ_PATH = "/Library/Application Support/ddPrimer/Multiz/multiz"
BLASTN_PATH = "/Library/Application Support/ddPrimer/Blast+/bin/blastn"
DB_PATH = "/Library/Application Support/ddPrimer/Tair DB/TAIR10_db"
```

## Usage

1. Run the script:
   ```bash
   python ddPrimer.py
   ```

2. Follow the interactive prompts:
   - Choose between comparing two genomes or multiple genomes
   - Choose to run LastZ or use existing MAF files
   - Select reference FASTA file (reference genome)
   - Select query FASTA file (query genome)
   - Select additional genomes if using multiple genome comparison
   - Select GFF3 annotation file

3. The script will run through all analysis steps with progress indicators:
   - For multiple genome comparison, it will run PLastZ pairwise between all genomes
   - If using multiple genomes, it will run MultiZ to combine alignments
   - Process gene annotations from GFF
   - Extract conserved segments from alignments
   - Design primers with Primer3
   - Filter for specificity using BLAST
   - Calculate thermodynamic properties with NUPACK

4. Results will be saved to an Excel file named `Primers_RefGenome_+_QueryGenome.xlsx`.

## About PLastZ

This version of ddPrimer includes PLastZ that offers improved performance over standard LastZ:

- Parallel processing of genome comparisons
- Better handling of files with spaces in paths
- MAF output formatting
- Robust sequence processing
- Optimized for multi-core processing
- Speeds up genome-to-genome alignments significantly
- Bundled with the ddPrimer package

## Cross-Platform Compatibility

While optimized for M1/M2 Macs, this script works on:

### macOS
- Works on both Apple Silicon and Intel Macs with appropriate path settings
- Includes special handling for macOS UI via Cocoa integration when available

### Linux
- Functions with path adjustments
- For Ubuntu/Debian, ensure build essentials are installed: `apt-get install build-essential`

### Windows
- Requires path adjustments for all tools
- Use forward slashes (/) or escaped backslashes (\\\\) in paths
- Install tools via Windows installers or WSL
- May require modifications to subprocess calls

To adapt for Windows:
```python
# Example path modifications for Windows in the script
BLASTN_PATH = "C:/Program Files/NCBI/blast/bin/blastn.exe"
MULTIZ_PATH = "C:/path/to/multiz.exe"
```

## Output Format

The Excel output includes columns for:
- Sequence identifiers
- Primer and probe sequences
- Melting temperatures
- Penalty scores
- BLAST specificity results
- ΔG values
- Amplicon details
- Genomic coordinates in both reference and query genomes

## Troubleshooting

### Common Issues

1. **LastZ errors**: Ensure reference and query FASTA files are properly formatted.

2. **PLastZ errors**: Make sure Samtools is installed and in your PATH.

3. **MultiZ errors**: If using multiple genomes, ensure MultiZ is properly installed and accessible.

4. **BLAST errors**: Verify BLAST installation and database paths.

5. **Memory errors**: Make sure you have sufficient RAM, especially for large genomes.

6. **No primers found**: Check Primer3 settings or adjust `PENALTY_MAX`.

## License

This script is provided for academic and research purposes. Please cite the corresponding tools (LastZ, MultiZ, Primer3, BLAST, NUPACK) in your publications.

## Citation

If you use this tool in your research, please cite:
- LastZ: Harris, R.S. (2007) Improved pairwise alignment of genomic DNA. Ph.D. Thesis
- PLastZ: Ho, A. (2023). PLastZ: Run lastz jobs in parallel. GitHub. https://github.com/AntoineHo/PLastZ
- MultiZ: Blanchette M, et al. (2004). Aligning multiple genomic sequences with the threaded blockset aligner. Genome Res. 14:708-715
- Primer3: Untergasser et al. (2012). Primer3—new capabilities and interfaces. Nucleic Acids Research, 40(15), e115
- BLAST: Camacho C et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10:421
- NUPACK: Dirks RM, Pierce NA (2004). An algorithm for computing nucleic acid base-pairing probabilities including pseudoknots. J Comput Chem, 25(10):1295-1304
- Samtools: Li H, et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16):2078-9

## Contact

For questions or support, please contact globus1333@gmail.com.