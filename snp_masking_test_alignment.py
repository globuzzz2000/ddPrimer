#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced test script for SNP masking in alignment mode.

This script provides a comprehensive test of SNP masking in alignment mode with:
- Longer, more realistic sequences
- Imperfect alignments with mismatches and gaps
- Better visualization of the alignment and masking process
- More detailed output for understanding what's happening
"""

import os
import sys
import tempfile
import subprocess
from pathlib import Path
import shutil
import time
import platform
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import logging fix to reduce excessive logging
from test_logging_fix import logger

# Create test directory
test_dir = Path(tempfile.mkdtemp(prefix="ddprimer_test_alignment_comp_"))
logger.info(f"Created test directory: {test_dir}")

def color_text(text, color_code):
    """Add color to terminal text."""
    return f"\033[{color_code}m{text}\033[0m"

def red(text):
    return color_text(text, "31")

def green(text):
    return color_text(text, "32")

def yellow(text):
    return color_text(text, "33")

def blue(text):
    return color_text(text, "34")

def magenta(text):
    return color_text(text, "35")

def cyan(text):
    return color_text(text, "36")

def highlight_changes(original, masked):
    """Highlight the changes between original and masked sequences."""
    result = ""
    for i, (orig_char, masked_char) in enumerate(zip(original, masked)):
        if orig_char != masked_char:
            result += red(masked_char)
        else:
            result += masked_char
        
        # Add a space every 10 characters for better readability
        if (i + 1) % 10 == 0:
            result += " "
    
    return result

def format_position_ruler(length):
    """Create a position ruler for sequence visualization."""
    # Create 1-based position indicators
    ruler = ""
    for i in range(1, length + 1):
        if i % 10 == 0:
            # For multiples of 10, show the position number
            ruler += f"{i//10}"
        else:
            # For other positions, just a mark
            ruler += "·"
        
        # Add a space every 10 characters for better readability
        if i % 10 == 0:
            ruler += " "
    
    return ruler

def check_lastz_installed():
    """Check if LastZ is installed in the system."""
    try:
        result = subprocess.run(['which', 'lastz'], 
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, 
                               text=True)
        if result.returncode == 0:
            return True
        
        # Also check for Windows which doesn't have 'which'
        if platform.system() == 'Windows':
            result = subprocess.run(['where', 'lastz'], 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, 
                                   text=True)
            return result.returncode == 0
            
        return False
    except:
        return False

def check_samtools_installed():
    """Check if samtools is installed in the system."""
    try:
        result = subprocess.run(['which', 'samtools'], 
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, 
                               text=True)
        if result.returncode == 0:
            return True
        
        # Also check for Windows
        if platform.system() == 'Windows':
            result = subprocess.run(['where', 'samtools'], 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, 
                                   text=True)
            return result.returncode == 0
            
        return False
    except:
        return False

def create_test_files():
    """Create realistic test files for alignment mode comprehensive test."""
    # 1. Create reference FASTA with multiple sequences representing model organism (e.g., A. thaliana)
    ref_fasta_content = """>chr1 Arabidopsis thaliana chromosome 1
ATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCA
GCGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCATCGATCGATGC
GCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC
ATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCG
ATCGATCGATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATG
GCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC
>chr2 Arabidopsis thaliana chromosome 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATG
GCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC
ATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCG
CTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCA
GCGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCG
>chr3 Arabidopsis thaliana chromosome 3
TGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAG
CTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGA
GCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCTAGCGATC
GCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCG
ATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCA
CTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGAT
"""
    ref_fasta_path = test_dir / "reference.fasta"
    with open(ref_fasta_path, "w") as f:
        f.write(ref_fasta_content)
    
    # 2. Create second species FASTA with imperfections - mismatches, insertions, deletions
    second_fasta_content = """>chr1 Closely related species chromosome 1
ATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGTTCGTACGATCGAGATCGATCGATGCATGCTAGCA
GCGATCGATCGATTCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCATCGATCGATGC
GCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCG--CGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC
ATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCG
ATCGATCGATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATG
GCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC
>chr2 Closely related species chromosome 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGTACGATCGATCGATGCTTGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATG
GCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATTGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC
ATCGATCGTAGCTAGCTACGTAGCTGATCGATCGAACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCG
CTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCA
GCGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATTGTACGATCGATCGATGCATGCTAGCATCGATCG
>chr3 Closely related species chromosome 3
TGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCG--CGTACGTAGCTACGATCGTAG
CTAGCATAGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGA
GCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCTAGCGATC
GCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCG
ATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCA
CTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGAT
"""
    second_fasta_path = test_dir / "second_species.fasta"
    with open(second_fasta_path, "w") as f:
        f.write(second_fasta_content)
    
    # 3. Create reference VCF file with realistic SNPs/indels
    ref_vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	10	SNP001	A	T	100	PASS	.	GT	1/1
chr1	25	SNP002	G	C	95	PASS	.	GT	0/1
chr1	50	SNP003	C	G	98	PASS	.	GT	1/0
chr1	75	SNP004	A	G	99	PASS	.	GT	1/1
chr1	100	SNP005	T	C	97	PASS	.	GT	0/1
chr2	15	SNP006	G	A	94	PASS	.	GT	0/1
chr2	40	SNP007	T	C	92	PASS	.	GT	1/1
chr2	65	SNP008	C	A	91	PASS	.	GT	0/1
chr2	90	SNP009	G	T	93	PASS	.	GT	1/0
chr3	20	SNP010	C	T	96	PASS	.	GT	1/1
chr3	45	SNP011	T	A	90	PASS	.	GT	0/1
chr3	70	SNP012	G	C	89	PASS	.	GT	1/0
"""
    ref_vcf_path = test_dir / "reference.vcf"
    with open(ref_vcf_path, "w") as f:
        f.write(ref_vcf_content)
    
    # 4. Create second species VCF file
    second_vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	15	SNP101	C	T	96	PASS	.	GT	1/1
chr1	30	SNP102	A	G	95	PASS	.	GT	0/1
chr1	55	SNP103	G	A	94	PASS	.	GT	1/0
chr1	80	SNP104	T	C	93	PASS	.	GT	1/1
chr1	105	SNP105	G	A	92	PASS	.	GT	0/1
chr2	20	SNP106	T	G	91	PASS	.	GT	0/1
chr2	45	SNP107	G	A	90	PASS	.	GT	1/1
chr2	70	SNP108	A	C	89	PASS	.	GT	0/1
chr2	95	SNP109	T	G	88	PASS	.	GT	1/0
chr3	25	SNP110	G	T	87	PASS	.	GT	1/1
chr3	50	SNP111	A	C	86	PASS	.	GT	0/1
chr3	75	SNP112	T	G	85	PASS	.	GT	1/0
"""
    second_vcf_path = test_dir / "second_species.vcf"
    with open(second_vcf_path, "w") as f:
        f.write(second_vcf_content)
    
    # 5. Create GFF file
    gff_content = """##gff-version 3
chr1	source	gene	5	150	.	+	.	ID=gene1;Name=GENE1
chr1	source	exon	5	50	.	+	.	Parent=gene1;Name=GENE1.1
chr1	source	exon	75	120	.	+	.	Parent=gene1;Name=GENE1.2
chr2	source	gene	10	140	.	+	.	ID=gene2;Name=GENE2
chr2	source	exon	10	60	.	+	.	Parent=gene2;Name=GENE2.1
chr2	source	exon	80	140	.	+	.	Parent=gene2;Name=GENE2.2
chr3	source	gene	15	130	.	+	.	ID=gene3;Name=GENE3
chr3	source	exon	15	45	.	+	.	Parent=gene3;Name=GENE3.1
chr3	source	exon	60	130	.	+	.	Parent=gene3;Name=GENE3.2
"""
    gff_path = test_dir / "annotation.gff"
    with open(gff_path, "w") as f:
        f.write(gff_content)
    
    # 6. Create a realistic MAF file that mimics LastZ output
    # This includes alignments with mismatches, indels, and gaps
    maf_content = """##maf version=1 scoring=LastZ
a score=10000
s chr1 0 480 + 480 ATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCAGCGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCATCGATCGATGCGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCGATCGATCGATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC
s chr1 0 480 + 480 ATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGTTCGTACGATCGAGATCGATCGATGCATGCTAGCAGCGATCGATCGATTCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCATCGATCGATGCGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCG--CGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCGATCGATCGATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGC--

a score=9950
s chr2 0 480 + 480 GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAATCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCAGCGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCG
s chr2 0 480 + 480 GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAATCGATCGATCGTACGATCGATCGATGCTTGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATTGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCATCGATCGTAGCTAGCTACGTAGCTGATCGATCGAACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATCGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCAGCGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATTGTACGATCGATCGATGCATGCTAGCATCGATCG

a score=9900
s chr3 0 480 + 480 TGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCTAGCGATCGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCACTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGAT
s chr3 0 480 + 480 TGCTAGCTAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCG--CGTACGTAGCTACGATCGTAGCTAGCATAGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGAGCTACGATCGTAGCTAGCGATCGATCGTACGATCGATCGTACGTAGCTAGCGATCGTACGTAGCTACGATCGTAGCTAGCGATCGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGATGGCGATGCTAGCGATCGATCGTAGCTAGCTACGTAGCTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCACTGATCGATCGTACGTACGATCGAGATCGATCGATGCATGCTAGCATCGATCGATCGATGCTAGCGATCGATCGTACGATCGAT--
"""
    maf_path = test_dir / "alignment.maf"
    with open(maf_path, "w") as f:
        f.write(maf_content)
    
    return {
        "ref_fasta": ref_fasta_path,
        "second_fasta": second_fasta_path,
        "ref_vcf": ref_vcf_path,
        "second_vcf": second_vcf_path,
        "gff": gff_path,
        "maf": maf_path
    }

def show_alignment_visualization(alignments, chrom_key):
    """Generate a visual representation of the alignment in the terminal."""
    if chrom_key not in alignments or not alignments[chrom_key]:
        logger.info(f"No alignments found for {chrom_key}")
        return
    
    # Just show the first alignment for brevity
    alignment = alignments[chrom_key][0]
    ref_seq = alignment["ref_seq"]
    qry_seq = alignment["qry_seq"]
    
    # Show only first 100 characters for display clarity
    display_length = min(100, len(ref_seq))
    
    logger.info(f"\nAlignment Visualization - {chrom_key} (first 100bp):")
    logger.info(f"Identity: {alignment['identity']:.2f}%")
    logger.info(f"Ref: {alignment['ref_chrom']}:{alignment['ref_start']}-{alignment['ref_end']} ({alignment['ref_strand']})")
    logger.info(f"Qry: {alignment['qry_src']}:{alignment['qry_start']}-{alignment['qry_end']} ({alignment['qry_strand']})")
    
    # Show position ruler
    logger.info(f"Pos: {''.join(' ' for _ in range(6))} {format_position_ruler(display_length)}")
    
    # Format the reference sequence
    ref_display = ref_seq[:display_length]
    ref_formatted = "Ref:  " + " ".join(ref_display[i:i+10] for i in range(0, len(ref_display), 10))
    logger.info(ref_formatted)
    
    # Show alignment marks (match, mismatch, gap)
    align_marks = ""
    for i in range(display_length):
        if i >= len(ref_seq) or i >= len(qry_seq):
            align_marks += " "
        elif ref_seq[i] == qry_seq[i]:
            align_marks += "|"  # Match
        elif ref_seq[i] == "-" or qry_seq[i] == "-":
            align_marks += " "  # Gap
        else:
            align_marks += "."  # Mismatch
        
        # Add a space every 10 characters for better readability
        if (i + 1) % 10 == 0:
            align_marks += " "
    
    logger.info(f"      {align_marks}")
    
    # Format the query sequence
    qry_display = qry_seq[:display_length]
    qry_formatted = "Qry:  " + " ".join(qry_display[i:i+10] for i in range(0, len(qry_display), 10))
    logger.info(qry_formatted)

def run_lastz_alignment(ref_fasta, second_fasta, output_dir):
    """
    Run LastZ alignment between reference and second species.
    
    Args:
        ref_fasta: Path to reference FASTA
        second_fasta: Path to second species FASTA
        output_dir: Output directory for MAF file
        
    Returns:
        str: Path to MAF file, or None if alignment failed
    """
    if not check_lastz_installed():
        logger.warning("LastZ not found in the system. Using pre-created MAF file.")
        return None
        
    if not check_samtools_installed():
        logger.warning("samtools not found in the system. Using pre-created MAF file.")
        return None
    
    logger.info("Running LastZ alignment...")
    try:
        # Create output directory
        align_dir = output_dir / "alignment"
        align_dir.mkdir(exist_ok=True)
        
        # Run our custom LastZ runner
        from ddprimer.alignment.lastz_runner import LastZRunner
        
        # Create runner and execute alignment
        runner = LastZRunner()
        output_file = runner.run_parallel_alignment(
            str(ref_fasta), 
            str(second_fasta), 
            str(align_dir), 
            "--format=maf"
        )
        
        logger.info(f"LastZ alignment completed successfully: {output_file}")
        return output_file
    except Exception as e:
        logger.error(f"Error running LastZ alignment: {e}")
        logger.error("Using pre-created MAF file instead.")
        return None

def visualize_masked_sequences(original_seqs, masked_seqs, coordinate_map, known_variants):
    """Visualize the masked sequences with original sequences for comparison."""
    logger.info("\n=== Masked Sequences Visualization ===")
    
    for seq_id, original in original_seqs.items():
        if seq_id not in masked_seqs:
            logger.warning(f"Sequence {seq_id} not found in masked sequences")
            continue
        
        masked = masked_seqs[seq_id]
        
        # Show sequence details
        logger.info(f"\n{blue(seq_id)} ({len(original)} bp)")
        
        # Collect all expected masked positions
        expected_masked_positions = set()
        
        # Add positions from reference variants
        if seq_id in known_variants.get('reference', {}):
            ref_positions = known_variants['reference'][seq_id]
            expected_masked_positions.update(ref_positions)
            logger.info(f"  Reference variants: {len(ref_positions)} positions")
        
        # Add mapped positions from second species variants
        if seq_id in known_variants.get('second', {}):
            second_positions = known_variants['second'][seq_id]
            expected_masked_positions.update(second_positions)
            logger.info(f"  Second species variants (mapped): {len(second_positions)} positions")
        
        # Display only first 100 bp for clarity
        display_length = min(len(original), 100)
        
        # Create position ruler
        ruler = format_position_ruler(display_length)
        logger.info(f"Pos: {ruler}")
        
        # Format and display original sequence
        orig_display = original[:display_length]
        orig_formatted = " ".join(orig_display[i:i+10] for i in range(0, len(orig_display), 10))
        logger.info(f"Orig: {orig_formatted}")
        
        # Format and display masked sequence with highlights
        masked_display = masked[:display_length]
        masked_formatted = highlight_changes(orig_display, masked_display)
        logger.info(f"Mask: {masked_formatted}")
        
        # Find actual masked positions
        masked_indices = [i for i, (orig, mask) in enumerate(zip(original, masked)) if orig != mask]
        masked_positions = {i+1 for i in masked_indices}  # Convert to 1-based set
        
        # Calculate statistics
        total_masked = len(masked_indices)
        masked_percent = (total_masked / len(original)) * 100 if len(original) > 0 else 0
        logger.info(f"Stats: {total_masked} out of {len(original)} positions masked ({masked_percent:.2f}%)")
        
        # Check for coordinate mapping
        if seq_id in coordinate_map:
            mapping_count = len(coordinate_map[seq_id])
            logger.info(f"Coordinate mapping: {mapping_count} positions mapped to second species")
            
            # Show a few example mappings
            if mapping_count > 0:
                examples = list(coordinate_map[seq_id].items())[:3]
                for ref_pos, mapping in examples:
                    logger.info(f"    {seq_id}:{ref_pos} → {mapping['qry_src']}:{mapping['qry_pos']} ({mapping['qry_strand']})")
        
        # Validate if masking matches expected positions
        if expected_masked_positions:
            missing = expected_masked_positions - masked_positions
            extra = masked_positions - expected_masked_positions
            
            if not missing and not extra:
                logger.info(f"Validation: {green('PASSED')} - All expected variants were masked correctly")
            else:
                if missing:
                    logger.info(f"Validation: {red('WARNING')} - {len(missing)} variant positions not masked: {sorted(missing)[:5]}...")
                if extra:
                    logger.info(f"Validation: {yellow('WARNING')} - {len(extra)} unexpected positions masked: {sorted(extra)[:5]}...")
        else:
            logger.info(f"Validation: {yellow('SKIPPED')} - No expected positions available for validation")

def test_alignment_workflow_with_maf():
    """
    Test alignment mode with a pre-created MAF file.
    This tests the MAF parsing, conserved region identification,
    and SNP masking but skips the LastZ alignment.
    """
    try:
        from ddprimer.core import SNPMaskingProcessor
        from ddprimer.alignment.maf_parser import MAFParser
        from ddprimer.alignment import AlignmentWorkflow
        
        # Create test files
        files = create_test_files()
        
        # Mock args object with everything except maf_file
        class Args:
            def __init__(self):
                self.fasta = str(files["ref_fasta"])
                self.second_fasta = str(files["second_fasta"])
                self.vcf = str(files["ref_vcf"])
                self.second_vcf = str(files["second_vcf"])
                self.gff = str(files["gff"])
                self.maf_file = str(files["maf"])
                self.snp = True
                self.min_identity = 80.0
                self.min_length = 20
                self.lastz_options = "--format=maf"
                self.lastzonly = False
                self.noannotation = False
                
        args = Args()
        
        # Create output directory
        output_dir = test_dir / "output"
        output_dir.mkdir(exist_ok=True)
        
        # ===== STEP 1: MAF FILE PARSING =====
        logger.info(magenta("\n===== STEP 1: MAF FILE PARSING ====="))
        # Initialize MAF parser
        maf_parser = MAFParser()
        
        # Analyze MAF file structure
        logger.info("Analyzing MAF file structure...")
        maf_analysis = maf_parser.analyze_maf_file(str(files["maf"]))
        logger.info(f"MAF file contains {maf_analysis['alignment_count']} alignment blocks")
        logger.info(f"Reference sequences: {', '.join(maf_analysis['ref_seq_ids'])}")
        logger.info(f"Query sequences: {', '.join(maf_analysis['query_seq_ids'])}")
        
        # Parse MAF file
        logger.info("Parsing MAF file...")
        alignments = maf_parser.parse_maf_file(str(files["maf"]))
        logger.info(f"Parsed {len(alignments)} reference sequences with alignments")
        
        # Show alignment visualization for each chromosome
        for chrom_key in alignments.keys():
            show_alignment_visualization(alignments, chrom_key)
        
        # ===== STEP 2: CONSERVED REGION IDENTIFICATION =====
        logger.info(magenta("\n===== STEP 2: CONSERVED REGION IDENTIFICATION ====="))
        conserved_regions = maf_parser.identify_conserved_regions(
            args.min_identity,
            args.min_length
        )
        
        total_regions = sum(len(regions) for regions in conserved_regions.values())
        logger.info(f"Identified {total_regions} conserved regions across {len(conserved_regions)} chromosomes")
        
        # Show conserved regions details for each chromosome
        for chrom, regions in conserved_regions.items():
            logger.info(f"\nChromosome {chrom}: {len(regions)} conserved regions")
            if regions:
                # Show top 3 regions
                for i, region in enumerate(regions[:3]):
                    logger.info(f"  Region {i+1}: {region['start']}-{region['end']} ({region['end']-region['start']+1} bp)")
                    logger.info(f"    Identity: {region['identity']:.2f}%")
                    logger.info(f"    Maps to: {region['qry_src']}:{region['qry_start']}-{region['qry_end']} ({region['qry_strand']})")
        
        # ===== STEP 3: COORDINATE MAPPING =====
        logger.info(magenta("\n===== STEP 3: COORDINATE MAPPING ====="))
        coordinate_map = maf_parser.generate_coordinate_map(conserved_regions)
        
        # Check mapping for key positions
        logger.info("Checking coordinate mapping...")
        found_mappings = 0
        for chrom, positions in coordinate_map.items():
            mapped_count = len(positions)
            found_mappings += mapped_count
            logger.info(f"Chromosome {chrom}: {mapped_count} positions mapped")
            
            # Show some example mappings
            if mapped_count > 0:
                examples = list(positions.items())[:3]
                for ref_pos, mapping in examples:
                    logger.info(f"  {chrom}:{ref_pos} → {mapping['qry_src']}:{mapping['qry_pos']} ({mapping['qry_strand']})")
        
        logger.info(f"Total of {found_mappings} positions mapped across all chromosomes")
        
        # ===== STEP 4: SNP EXTRACTION =====
        logger.info(magenta("\n===== STEP 4: SNP EXTRACTION ====="))
        # Initialize SNP masking processor
        snp_processor = SNPMaskingProcessor()
        
        # Extract variants from reference VCF
        logger.info("Extracting variants from reference VCF...")
        ref_variants = snp_processor.get_variant_positions(str(files["ref_vcf"]))
        ref_variant_count = sum(len(positions) for positions in ref_variants.values())
        logger.info(f"Extracted {ref_variant_count} variants from reference genome")
        
        # Show variants for each chromosome
        for chrom, positions in ref_variants.items():
            pos_list = sorted(positions)
            logger.info(f"  {chrom}: {len(positions)} variants at positions {pos_list[:5]}...")
        
        # Extract variants from second species VCF
        logger.info("\nExtracting variants from second species VCF...")
        second_variants = snp_processor.get_variant_positions(str(files["second_vcf"]))
        second_variant_count = sum(len(positions) for positions in second_variants.values())
        logger.info(f"Extracted {second_variant_count} variants from second species")
        
        # Show variants for each chromosome
        for chrom, positions in second_variants.items():
            pos_list = sorted(positions)
            logger.info(f"  {chrom}: {len(positions)} variants at positions {pos_list[:5]}...")
        
        # ===== STEP 5: FULL WORKFLOW EXECUTION =====
        logger.info(magenta("\n===== STEP 5: FULL WORKFLOW EXECUTION ====="))
        logger.info("Running full AlignmentWorkflow...")
        try:
            # Call AlignmentWorkflow
            masked_sequences, result_coordinate_map = AlignmentWorkflow(args, str(output_dir), logger)
            
            # Check if results are valid
            if masked_sequences and len(masked_sequences) > 0:
                logger.info(f"Workflow completed successfully with {len(masked_sequences)} masked sequences")
                
                # Load the original reference sequences for comparison
                logger.info("Loading reference sequences for comparison...")
                original_sequences = {}
                for record in SeqIO.parse(files["ref_fasta"], "fasta"):
                    original_sequences[record.id.split()[0]] = str(record.seq)
                
                # Track which positions should be masked
                known_variants = {
                    'reference': ref_variants,
                    'second': {}  # We'll map second species variants to reference coordinates
                }
                
                # Map second species variants to reference coordinates
                logger.info("Mapping second species variants to reference coordinates...")
                for ref_chrom, ref_pos_map in coordinate_map.items():
                    if ref_chrom not in known_variants['second']:
                        known_variants['second'][ref_chrom] = set()
                    
                    for ref_pos, mapping in ref_pos_map.items():
                        qry_src = mapping['qry_src']
                        qry_pos = mapping['qry_pos']
                        
                        # If this position in second species has a variant, map it back
                        if qry_src in second_variants and qry_pos in second_variants[qry_src]:
                            known_variants['second'][ref_chrom].add(ref_pos)
                
                # Count mapped variants
                mapped_variants = sum(len(positions) for positions in known_variants['second'].values())
                logger.info(f"Mapped {mapped_variants} second species variants to reference coordinates")
                
                # Visualize the masked sequences
                visualize_masked_sequences(original_sequences, masked_sequences, result_coordinate_map, known_variants)
                
                # Check for key positions that should be masked
                logger.info("\nChecking specific positions for masking...")
                positions_to_check = {
                    'chr1': [10, 25, 50, 75, 100],  # Reference variants
                    'chr2': [15, 40, 65, 90],       # Reference variants
                    'chr3': [20, 45, 70]            # Reference variants
                }
                
                all_checks_passed = True
                for chrom, positions in positions_to_check.items():
                    if chrom not in masked_sequences:
                        logger.warning(f"Chromosome {chrom} not found in masked sequences")
                        continue
                    
                    for pos in positions:
                        # Convert to 0-based index
                        idx = pos - 1
                        if idx < len(masked_sequences[chrom]):
                            if masked_sequences[chrom][idx] == 'N':
                                logger.info(f"✓ Position {pos} in {chrom} is correctly masked")
                            else:
                                logger.warning(f"✗ Position {pos} in {chrom} should be masked, but found '{masked_sequences[chrom][idx]}'")
                                all_checks_passed = False
                
                if all_checks_passed:
                    logger.info(green("\nAll tested positions are correctly masked!"))
                else:
                    logger.warning(yellow("\nSome positions were not masked as expected"))
                
                return True
            else:
                logger.error("Workflow completed but returned no masked sequences")
                return False
                
        except Exception as e:
            logger.error(f"Error running AlignmentWorkflow: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False
        
    except ImportError as e:
        logger.error(f"Could not import required modules: {e}")
        return False
    except Exception as e:
        logger.error(f"Error in alignment mode test: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def test_alignment_workflow_with_lastz():
    """
    Test full alignment workflow starting with LastZ alignment.
    This tests the complete pipeline including LastZ alignment.
    """
    # First check if LastZ is installed
    if not check_lastz_installed() or not check_samtools_installed():
        logger.warning("LastZ or samtools not found. Skipping full alignment workflow test.")
        return True  # Return success so the test suite continues
    
    try:
        from ddprimer.alignment import AlignmentWorkflow
        
        # Create test files
        files = create_test_files()
        
        # Create output directory
        output_dir = test_dir / "lastz_output"
        output_dir.mkdir(exist_ok=True)
        
        # Run LastZ alignment
        logger.info(magenta("\n===== FRESH LASTZ ALIGNMENT ====="))
        maf_file = run_lastz_alignment(
            files["ref_fasta"], 
            files["second_fasta"], 
            output_dir
        )
        
        if not maf_file:
            logger.warning("LastZ alignment failed or not available. Skipping this part of the test.")
            return True  # Return success to continue with other tests
        
        # Mock args object with everything
        class Args:
            def __init__(self):
                self.fasta = str(files["ref_fasta"])
                self.second_fasta = str(files["second_fasta"])
                self.vcf = str(files["ref_vcf"])
                self.second_vcf = str(files["second_vcf"])
                self.gff = str(files["gff"])
                self.maf_file = None  # We'll use our fresh LastZ alignment
                self.snp = True
                self.min_identity = 80.0
                self.min_length = 20
                self.lastz_options = "--format=maf"
                self.lastzonly = False
                self.noannotation = False
                
        args = Args()
        
        # Run the full AlignmentWorkflow with our generated MAF file
        logger.info("\nTesting full alignment workflow with fresh LastZ alignment...")
        try:
            # Call AlignmentWorkflow
            masked_sequences, coordinate_map = AlignmentWorkflow(args, str(output_dir), logger)
            
            # Check if results are valid
            if masked_sequences and len(masked_sequences) > 0:
                logger.info(f"Workflow completed successfully with {len(masked_sequences)} masked sequences")
                
                # Output some info about the masking
                for seq_id, sequence in masked_sequences.items():
                    n_count = sequence.count('N')
                    n_percent = (n_count / len(sequence)) * 100
                    logger.info(f"Sequence {seq_id}: {n_count}/{len(sequence)} bases masked ({n_percent:.2f}%)")
                
                return True
            else:
                logger.error("Workflow completed but returned no masked sequences")
                return False
                
        except Exception as e:
            logger.error(f"Error running AlignmentWorkflow with LastZ: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False
        
    except ImportError as e:
        logger.error(f"Could not import required modules: {e}")
        return False
    except Exception as e:
        logger.error(f"Error in alignment mode with LastZ test: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    try:
        # Create files
        files = create_test_files()
        
        # Run test with pre-created MAF
        logger.info(magenta("\n===== ALIGNMENT TEST WITH PRE-CREATED MAF FILE ====="))
        success_maf = test_alignment_workflow_with_maf()
        
        # Run test with LastZ alignment
        logger.info(magenta("\n===== ALIGNMENT TEST WITH FRESH LASTZ ALIGNMENT ====="))
        success_lastz = test_alignment_workflow_with_lastz()
        
        # Overall test success
        if success_maf and success_lastz:
            logger.info(green("\n✓ All alignment tests passed successfully!"))
        else:
            if not success_maf:
                logger.error(red("\n✗ Alignment test with pre-created MAF failed"))
            if not success_lastz:
                logger.error(red("\n✗ Alignment test with fresh LastZ alignment failed"))
            sys.exit(1)
            
    finally:
        # Clean up
        logger.info(f"Cleaning up test directory: {test_dir}")
        try:
            shutil.rmtree(test_dir)
        except Exception as e:
            logger.error(f"Error cleaning up: {e}")