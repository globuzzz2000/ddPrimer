#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced test script for SNP masking in direct mode.

This script provides a comprehensive test of SNP masking in direct mode with:
- Longer, more varied sequences
- Multiple sequence variants
- More realistic sequence names and contexts
- Improved human-readable output
"""

import os
import sys
import tempfile
from pathlib import Path
import shutil
import csv
import re

# Import logging fix
from test_logging_fix import logger

# Create test directory
test_dir = Path(tempfile.mkdtemp(prefix="ddprimer_test_direct_"))
logger.info(f"Created test directory: {test_dir}")

def create_test_files():
    """Create enhanced test files for direct mode SNP masking test."""
    # 1. Create reference FASTA with longer, more varied sequences
    ref_fasta_content = """>chr1 Arabidopsis thaliana chromosome 1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
GCTAGCTAGCTAGCTCGATCGATCGTAGCTAGCTAGCTGACTCGATCGATCGTAGCTGGCTAGCTCGATCGATCGTAGCTCGAT
TATCCGATCGTAGCTGGCTATCCGATCGTAGCACTGCATACGTACGTCAGTAGCATCGATGCATACGTACGTCAGTAGCATCGA
GACTGCATACGTACGTCAGATCGTAGCTGGCTCCGAACGTACGTCAGATCGTAGCTCGATCGATGCATACGTACGTCAGATCGT
TATTGCTACGTCAGATCGTAGCTCGATCGATGCATACGTACTCAGATCGTAGCTCGACGTCAGTAGCATCGATGCATACGTACG
GACTACGTCAGATCGTAGCTCGATCACGTACGTCAGTAGCATCGATGACGTACGTCAGTAGCATCGATGCATACGCCAGTAGCA
>chr2 Arabidopsis thaliana chromosome 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
TGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCA
ACGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGT
CTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCA
TACGCGATCGACTCGATCATACGCGATCGACTCGATCACACGTCAGTAGCATCGACGTACGCGATCGACACGTCAGTAGCATCG
TGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCAT
>chr3 Arabidopsis thaliana chromosome 3
CGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTC
TCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCAT
TGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
TACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGAT
"""
    ref_fasta_path = test_dir / "reference.fasta"
    with open(ref_fasta_path, "w") as f:
        f.write(ref_fasta_content)
    
    # 2. Create more comprehensive VCF file with variants on different chromosomes
    vcf_content = """##fileformat=VCFv4.2
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
    vcf_path = test_dir / "variants.vcf"
    with open(vcf_path, "w") as f:
        f.write(vcf_content)
    
    # 3. Create more realistic CSV file with varied sequences
    # Include naming patterns that match chromosomes for matching logic testing
    csv_content = [
        ["name", "sequence"],
        # Gene from chr1
        ["AT1G00630_CDS_chr1", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCGCTAGCTAGCTAGCTCGATCGATCGTAGCTAG"],
        # Gene from chr2
        ["AT2G15690_CDS_chromosome2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTTGCATGCATGCATGCATGCTGCATGCATGCAT"],
        # Gene from chr3
        ["AT3G44870_CDS_Chr3", "CGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTTCGATCATACGCGATCGACTCGATCATACGC"],
        # Gene with no clear chromosome match
        ["AT4G12345_fragment", "TGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGC"],
        # Longer sequence spanning multiple variants
        ["AT1G76540_promoter_region_chr1", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCGCTAGCTAGCTAGCTCGATCGATCGTAGCTAGCTAGCTGACTCGATCGATCGTAGCTGGCTAGCT"],
    ]
    
    csv_path = test_dir / "sequences.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(csv_content)
    
    return {
        "ref_fasta": ref_fasta_path,
        "vcf": vcf_path,
        "csv": csv_path
    }

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

def test_direct_mode_masking():
    """Test SNP masking in direct mode."""
    try:
        from ddprimer.core import SNPMaskingProcessor
        
        # Create test files
        files = create_test_files()
        
        # Initialize SNP masking processor
        processor = SNPMaskingProcessor()
        
        # Get variant positions from VCF
        variants = processor.get_variant_positions(str(files["vcf"]))
        logger.info(f"Extracted variants from VCF: {sum(len(pos) for pos in variants.values())} positions across {len(variants)} chromosomes")
        for chrom, positions in variants.items():
            logger.info(f"  {chrom}: {len(positions)} variants at positions {sorted(positions)[:5]}...")
        
        # Load sequences from CSV
        sequences = {}
        with open(files["csv"], "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sequences[row["name"]] = row["sequence"]
        
        logger.info(f"Loaded {len(sequences)} sequences from CSV:")
        for seq_id, seq in sequences.items():
            logger.info(f"  {seq_id}: {len(seq)} bp")
        
        # Mock the direct mode workflow - first try to match sequences to chromosomes
        logger.info("\n=== Testing Chromosome Matching Logic ===")
        seq_to_chrom = {}
        for seq_id in sequences.keys():
            # Try to match the sequence to a chromosome
            matched_chrom = None
            
            # Method 1: Check for explicit chromosome indicators in sequence ID
            for chrom in variants.keys():
                # Look for patterns like "chr1", "Chr1", "chromosome1"
                patterns = [
                    f"{chrom}$",                 # Ends with exact chrom name
                    f"{chrom}_",                 # Contains chrom name followed by underscore
                    f"chromosome{chrom[-1]}",    # "chromosome" followed by chrom number
                    f"Chr{chrom[-1]}",           # "Chr" followed by chrom number
                    f"chr{chrom[-1]}"            # "chr" followed by chrom number
                ]
                
                for pattern in patterns:
                    if re.search(pattern, seq_id, re.IGNORECASE):
                        matched_chrom = chrom
                        logger.info(f"  {seq_id} → {green(chrom)} (pattern match: {pattern})")
                        break
                
                if matched_chrom:
                    break
            
            if not matched_chrom:
                logger.info(f"  {seq_id} → {yellow('No chromosome match found')}")
            
            seq_to_chrom[seq_id] = matched_chrom
        
        # Now mask each sequence with variants from its matched chromosome
        logger.info("\n=== Testing SNP Masking with Chromosome Matching ===")
        masked_sequences = {}
        masked_positions = {}
        
        for seq_id, sequence in sequences.items():
            logger.info(f"\nProcessing sequence: {blue(seq_id)}")
            chrom = seq_to_chrom[seq_id]
            
            if not chrom:
                logger.info(f"  No chromosome match for {seq_id}, skipping masking")
                masked_sequences[seq_id] = sequence
                masked_positions[seq_id] = []
                continue
            
            # Get variants for this chromosome
            chrom_variants = variants.get(chrom, set())
            logger.info(f"  Using variants from {chrom}: {len(chrom_variants)} positions")
            
            # Find which variant positions are within sequence length
            applicable_variants = {pos for pos in chrom_variants if pos <= len(sequence)}
            logger.info(f"  {len(applicable_variants)} variant positions fall within sequence length: {sorted(applicable_variants)}")
            
            # Mask the sequence
            masked_sequence = processor.mask_variants(sequence, applicable_variants)
            masked_sequences[seq_id] = masked_sequence
            
            # Record which positions were actually masked
            mask_indices = [i for i, (orig, mask) in enumerate(zip(sequence, masked_sequence)) if orig != mask]
            masked_pos = [i+1 for i in mask_indices]  # Convert to 1-based
            masked_positions[seq_id] = masked_pos
            
            logger.info(f"  Masked {len(masked_pos)} positions: {masked_pos}")
        
        # Verify and display results
        logger.info("\n=== Visualization of Masked Sequences ===")
        for seq_id, original in sequences.items():
            masked = masked_sequences[seq_id]
            logger.info(f"\n{blue(seq_id)} (Chromosome: {seq_to_chrom.get(seq_id, 'Unknown')})")
            
            # Show position ruler for reference
            display_length = min(len(original), 100)  # Show at most 100 bp for clarity
            ruler = format_position_ruler(display_length)
            logger.info(f"Pos: {ruler}")
            
            # Show original sequence
            orig_display = original[:display_length]
            orig_formatted = " ".join(orig_display[i:i+10] for i in range(0, len(orig_display), 10))
            logger.info(f"Orig: {orig_formatted}")
            
            # Show masked sequence with highlights
            masked_display = masked[:display_length]
            masked_formatted = highlight_changes(orig_display, masked_display)
            logger.info(f"Mask: {masked_formatted}")
            
            # Show statistics
            total_masked = sum(1 for a, b in zip(original, masked) if a != b)
            masked_percent = (total_masked / len(original)) * 100 if len(original) > 0 else 0
            logger.info(f"Stats: {total_masked} out of {len(original)} positions masked ({masked_percent:.2f}%)")
            
            # Verify that masked positions match expected variants
            chrom = seq_to_chrom.get(seq_id)
            if chrom and chrom in variants:
                expected_positions = {pos for pos in variants[chrom] if pos <= len(original)}
                actual_positions = set(masked_positions[seq_id])
                
                if expected_positions == actual_positions:
                    logger.info(f"Validation: {green('PASSED')} - All expected variants were masked")
                else:
                    # Check if any expected positions were missed
                    missed = expected_positions - actual_positions
                    extra = actual_positions - expected_positions
                    
                    if missed:
                        logger.info(f"Validation: {red('FAILED')} - Missed positions: {missed}")
                    if extra:
                        logger.info(f"Validation: {red('FAILED')} - Extra positions masked: {extra}")
            else:
                logger.info(f"Validation: {yellow('SKIPPED')} - No chromosome match for validation")
        
        # Overall validation
        all_passed = True
        for seq_id in sequences.keys():
            chrom = seq_to_chrom.get(seq_id)
            if chrom and chrom in variants:
                expected_positions = {pos for pos in variants[chrom] if pos <= len(sequences[seq_id])}
                actual_positions = set(masked_positions[seq_id])
                
                if expected_positions != actual_positions:
                    all_passed = False
                    break
        
        if all_passed:
            logger.info(f"\n{green('✓ All sequences were correctly masked according to their chromosome matching')}")
            return True
        else:
            logger.error(f"\n{red('✗ Some sequences were not masked correctly')}")
            return False
        
    except ImportError:
        logger.error("Could not import required modules. Make sure ddprimer is installed.")
        return False
    except AssertionError as e:
        logger.error(f"Verification failed: {e}")
        return False
    except Exception as e:
        logger.error(f"Error in direct mode test: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    try:
        success = test_direct_mode_masking()
        if not success:
            sys.exit(1)
    finally:
        # Clean up
        logger.info(f"Cleaning up test directory: {test_dir}")
        shutil.rmtree(test_dir)