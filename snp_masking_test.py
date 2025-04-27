#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced test script for basic SNP masking functionality in ddPrimer pipeline.

This script provides comprehensive testing of the SNP masking core functionality with:
- Longer, more varied sequences
- Multiple variant types (SNPs at different positions)
- Better visualization of masking results
- More detailed output for understanding the process
"""

import os
import sys
import tempfile
from pathlib import Path
import random

# Import logging fix
from test_logging_fix import logger

# Create test directory
test_dir = Path(tempfile.mkdtemp(prefix="ddprimer_test_basic_"))
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

def create_test_files():
    """Create more comprehensive test files for SNP masking test."""
    # Create a reference FASTA file with longer sequences
    ref_fasta_content = """>chr1 Arabidopsis thaliana chromosome 1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
GCTAGCTAGCTAGCTCGATCGATCGTAGCTAGCTAGCTGACTCGATCGATCGTAGCTGGCTAGCTCGATCGATCGTAGCTCGAT
>chr2 Arabidopsis thaliana chromosome 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
TGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCATGCATGCATGCTGCATGCA
>chr3 Arabidopsis thaliana chromosome 3
CGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTCAGTAGCATCGACGTACGTC
TCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCATACGCGATCGACTCGATCAT
"""
    ref_fasta_path = test_dir / "reference.fasta"
    with open(ref_fasta_path, "w") as f:
        f.write(ref_fasta_content)
    
    # Create a more comprehensive VCF file with variants
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
    
    return {
        "ref_fasta": ref_fasta_path,
        "vcf": vcf_path
    }

def create_test_sequences():
    """Create a variety of test sequences for masking tests."""
    # Create some varied test sequences
    sequences = {
        "simple_sequence": "ATGCATGCATGCATGCATGC",
        "long_sequence": "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
        "repetitive_sequence": "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG",
        "gc_rich_sequence": "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
        "mixed_sequence": "ATGCTAGCTAGCTATGCTAGCTATGCGATCGATCGTAGCTAGCTACGTAGCTAGCTAGCTACGTAGCTAGCTAGCTAGCTAG"
    }
    
    return sequences

def test_snp_masking():
    """Test basic SNP masking functionality."""
    # Create test files
    files = create_test_files()
    
    # Import the SNP masking processor
    try:
        from ddprimer.core import SNPMaskingProcessor
    except ImportError:
        logger.error("Could not import SNPMaskingProcessor. Make sure ddprimer is installed.")
        return False
    
    # Initialize SNP masking processor
    processor = SNPMaskingProcessor()
    
    # Test 1: Extract variants from VCF
    logger.info(magenta("\n===== TEST 1: VARIANT EXTRACTION ====="))
    variants = processor.get_variant_positions(str(files["vcf"]))
    
    # Verify expected positions were extracted
    total_variants = sum(len(positions) for positions in variants.values())
    logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes:")
    
    # Show variants for each chromosome
    for chrom, positions in variants.items():
        logger.info(f"  {chrom}: {len(positions)} variants at positions {sorted(positions)}")
    
    # Test if all expected chromosomes and positions are present
    expected_chroms = {'chr1', 'chr2', 'chr3'}
    expected_positions = {
        'chr1': {10, 25, 50, 75, 100},
        'chr2': {15, 40, 65, 90},
        'chr3': {20, 45, 70}
    }
    
    all_chroms_present = set(variants.keys()) == expected_chroms
    all_positions_present = True
    
    for chrom, expected_pos in expected_positions.items():
        if chrom not in variants or not expected_pos.issubset(variants[chrom]):
            all_positions_present = False
            break
    
    if all_chroms_present and all_positions_present:
        logger.info(green("✓ All expected variants were correctly extracted from VCF"))
    else:
        logger.error(red("✗ Some expected variants were not extracted correctly"))
        if not all_chroms_present:
            logger.error(f"  Missing chromosomes: {expected_chroms - set(variants.keys())}")
        if not all_positions_present:
            for chrom, expected_pos in expected_positions.items():
                if chrom in variants:
                    missing = expected_pos - variants[chrom]
                    if missing:
                        logger.error(f"  Missing positions in {chrom}: {missing}")
        return False
    
    # Test 2: Simple masking test with different sequences
    logger.info(magenta("\n===== TEST 2: BASIC MASKING TESTS ====="))
    
    test_sequences = create_test_sequences()
    for name, sequence in test_sequences.items():
        logger.info(f"\nTesting with {name} ({len(sequence)} bp):")
        logger.info(f"Original: {sequence[:50]}..." if len(sequence) > 50 else f"Original: {sequence}")
        
        # Create a set of positions to test
        # For each sequence, we'll create positions that are multiples of 10, up to its length
        positions = {pos for pos in range(10, len(sequence) + 1, 10)}
        logger.info(f"Variant positions to mask: {positions}")
        
        # Mask the sequence
        masked_sequence = processor.mask_variants(sequence, positions)
        
        # Verify masking was applied correctly
        expected_indices = [pos - 1 for pos in positions]  # Convert 1-based to 0-based
        
        # Show position ruler for reference
        display_length = min(len(sequence), 100)  # Show at most 100 bp for clarity
        ruler = format_position_ruler(display_length)
        logger.info(f"Pos: {ruler}")
        
        # Show original and masked sequences for display length
        orig_display = sequence[:display_length]
        masked_display = masked_sequence[:display_length]
        
        # Format for display with spaces every 10 chars
        orig_formatted = " ".join(orig_display[i:i+10] for i in range(0, len(orig_display), 10))
        logger.info(f"Orig: {orig_formatted}")
        
        # Show masked sequence with highlighted changes
        masked_formatted = highlight_changes(orig_display, masked_display)
        logger.info(f"Mask: {masked_formatted}")
        
        # Verify each N is in the expected position and count them
        masked_indices = [i for i, char in enumerate(masked_sequence) if char == 'N']
        expected_masked_count = len([pos for pos in positions if pos <= len(sequence)])
        actual_masked_count = masked_sequence.count('N')
        
        # Check if masking was applied correctly
        masked_positions = {i+1 for i in masked_indices}  # Convert to 1-based for logging
        if set(masked_indices) == set(expected_indices) and actual_masked_count == expected_masked_count:
            logger.info(green(f"✓ All positions correctly masked: {masked_positions}"))
        else:
            missing = {pos for pos in positions if pos <= len(sequence)} - masked_positions
            extra = masked_positions - {pos for pos in positions if pos <= len(sequence)}
            
            if missing:
                logger.error(red(f"✗ Missing masked positions: {missing}"))
            if extra:
                logger.error(red(f"✗ Extra masked positions: {extra}"))
            if actual_masked_count != expected_masked_count:
                logger.error(red(f"✗ Expected {expected_masked_count} masked positions, found {actual_masked_count}"))
            return False
    
    # Test 3: Edge cases
    logger.info(magenta("\n===== TEST 3: EDGE CASES ====="))
    
    # Test with empty variant set
    test_sequence = "ATGCATGCATGCATGCATGC"
    logger.info("\nTesting with empty variant set:")
    masked = processor.mask_variants(test_sequence, set())
    if masked == test_sequence:
        logger.info(green("✓ Empty variant set - sequence remained unchanged"))
    else:
        logger.error(red(f"✗ Empty variant set should not change sequence: {masked}"))
        return False
    
    # Test with positions outside sequence length
    test_sequence = "ATGCATGCATGCATGCATGC"
    out_of_range = {5, 25, 50, 100}  # Position 25, 50, 100 are beyond sequence length
    logger.info(f"\nTesting with positions outside sequence length: {out_of_range}")
    masked = processor.mask_variants(test_sequence, out_of_range)
    
    # Only position 5 should be masked (1-based, so index 4)
    expected = list(test_sequence)
    expected[4] = 'N'
    expected = ''.join(expected)
    
    logger.info(f"Original: {test_sequence}")
    logger.info(f"Masked:   {highlight_changes(test_sequence, masked)}")
    logger.info(f"Expected: {highlight_changes(test_sequence, expected)}")
    
    if masked == expected:
        logger.info(green("✓ Out-of-range positions were correctly ignored"))
    else:
        logger.error(red(f"✗ Out-of-range positions not handled correctly"))
        return False
    
    # Test with empty sequence
    logger.info("\nTesting with empty sequence:")
    masked = processor.mask_variants("", {10, 20})
    if masked == "":
        logger.info(green("✓ Empty sequence - result remained empty"))
    else:
        logger.error(red(f"✗ Empty sequence should return empty result: {masked}"))
        return False
    
    # Test with many positions
    test_sequence = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    # Create a large set of positions (1 to length of sequence)
    many_positions = set(range(1, len(test_sequence) + 1))
    logger.info(f"\nTesting with all positions (1-{len(test_sequence)}):")
    masked = processor.mask_variants(test_sequence, many_positions)
    
    # All positions should be masked
    expected = 'N' * len(test_sequence)
    
    logger.info(f"Original: {test_sequence}")
    logger.info(f"Masked:   {highlight_changes(test_sequence, masked)}")
    
    if masked == expected:
        logger.info(green("✓ All positions were correctly masked"))
    else:
        logger.error(red(f"✗ All positions masking failed"))
        unmasked_count = sum(1 for c in masked if c != 'N')
        logger.error(f"  {unmasked_count} positions were not masked")
        return False
    
    # Test 4: Real chromosome sequences
    logger.info(magenta("\n===== TEST 4: REAL CHROMOSOME SEQUENCES ====="))
    
    # Load sequences from FASTA
    chromosomes = {}
    current_chrom = None
    with open(files["ref_fasta"], 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_chrom = line[1:].split()[0]  # Get chromosome name
                chromosomes[current_chrom] = ""
            elif current_chrom:
                chromosomes[current_chrom] += line
    
    # Test masking on each chromosome with its specific variants
    all_chrom_tests_passed = True
    for chrom, sequence in chromosomes.items():
        if chrom not in variants:
            logger.warning(f"No variants found for chromosome {chrom}")
            continue
        
        logger.info(f"\nTesting chromosome {blue(chrom)} ({len(sequence)} bp):")
        chrom_variants = variants[chrom]
        logger.info(f"Variant positions: {sorted(chrom_variants)}")
        
        # Mask the sequence
        masked_sequence = processor.mask_variants(sequence, chrom_variants)
        
        # Verify masking
        display_length = min(100, len(sequence))  # Show at most 100 bp
        
        # Position ruler
        ruler = format_position_ruler(display_length)
        logger.info(f"Pos: {ruler}")
        
        # Show original and masked sequences for display length
        orig_display = sequence[:display_length]
        masked_display = masked_sequence[:display_length]
        
        # Format for display with spaces every 10 chars
        orig_formatted = " ".join(orig_display[i:i+10] for i in range(0, len(orig_display), 10))
        logger.info(f"Orig: {orig_formatted}")
        
        # Show masked sequence with highlighted changes
        masked_formatted = highlight_changes(orig_display, masked_display)
        logger.info(f"Mask: {masked_formatted}")
        
        # Check if all expected positions are masked
        masked_indices = [i for i, char in enumerate(masked_sequence) if char == 'N']
        masked_positions = {i+1 for i in masked_indices}  # Convert to 1-based
        
        # Only consider positions within sequence length
        expected_positions = {pos for pos in chrom_variants if pos <= len(sequence)}
        
        if masked_positions == expected_positions:
            logger.info(green(f"✓ All {len(expected_positions)} expected positions correctly masked"))
        else:
            all_chrom_tests_passed = False
            missing = expected_positions - masked_positions
            extra = masked_positions - expected_positions
            
            if missing:
                logger.error(red(f"✗ Missing masked positions: {missing}"))
            if extra:
                logger.error(red(f"✗ Extra masked positions: {extra}"))
    
    if not all_chrom_tests_passed:
        logger.error(red("✗ Some chromosome masking tests failed"))
        return False
    
    # Test 5: Performance with large sequences
    logger.info(magenta("\n===== TEST 5: PERFORMANCE TEST WITH LARGE SEQUENCE ====="))
    
    # Create a large sequence (1 million bp)
    logger.info("Creating 1 million bp sequence...")
    bases = ['A', 'T', 'G', 'C']
    large_sequence = ''.join(random.choice(bases) for _ in range(1_000_000))
    
    # Create sparse variants (every 10,000 positions)
    sparse_variants = set(range(10000, 1_000_001, 10000))
    logger.info(f"Testing with {len(sparse_variants)} sparse variants (every 10,000 bp)")
    
    # Time the masking operation
    import time
    start_time = time.time()
    masked_large = processor.mask_variants(large_sequence, sparse_variants)
    elapsed_time = time.time() - start_time
    
    # Verify masking
    masked_count = masked_large.count('N')
    
    logger.info(f"Masked {masked_count} positions in {elapsed_time:.3f} seconds")
    
    if masked_count == len(sparse_variants):
        logger.info(green("✓ Large sequence masking succeeded"))
    else:
        logger.error(red(f"✗ Large sequence masking failed: Expected {len(sparse_variants)} masked positions, found {masked_count}"))
        return False
    
    # All tests passed!
    logger.info(green("\n✓ All SNP masking tests passed successfully!"))
    return True

if __name__ == "__main__":
    try:
        success = test_snp_masking()
        if success:
            logger.info(green("\n=== All tests passed successfully! ==="))
        else:
            logger.error(red("\n=== Some tests failed! ==="))
            sys.exit(1)
        
    except Exception as e:
        logger.error(red(f"\nError during testing: {e}"))
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    finally:
        # Clean up
        logger.info(f"Cleaning up test directory: {test_dir}")
        import shutil
        shutil.rmtree(test_dir)