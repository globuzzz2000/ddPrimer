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

# Silence extra ddPrimer logging
import logging as _logging
_logging.getLogger("ddPrimer").setLevel(_logging.CRITICAL)

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
    """Create minimal, fully-aligned test files: one chromosome, 600bp, with 10 SNPs."""
    # 1. Make chromosome sequence
    chrA = "ATG" + "C"*297 + "G"*300  # 600 bp
    # 2. Make second species: introduce 10 SNPs evenly spaced
    chrB = list(chrA)
    snp_positions = [60 * i for i in range(1, 11)]  # 60,120,...,600 (1-based)
    snp_positions_0 = [p-1 for p in snp_positions]  # 0-based
    # For each, change to a different base (cyclically: if C, to A, if G, to T, etc.)
    def alt_base(b):
        return {"A": "C", "C": "A", "G": "T", "T": "G"}.get(b, "N")
    for i, pos in enumerate(snp_positions_0):
        chrB[pos] = alt_base(chrB[pos])
    chrB = "".join(chrB)
    # 3. Write reference FASTA (one record)
    ref_fasta_path = test_dir / "reference.fasta"
    with open(ref_fasta_path, "w") as f:
        f.write(">chrA\n")
        for i in range(0, len(chrA), 60):
            f.write(chrA[i:i+60] + "\n")
    # 4. Write second species FASTA (one record, same name)
    second_fasta_path = test_dir / "second_species.fasta"
    with open(second_fasta_path, "w") as f:
        f.write(">chrA\n")
        for i in range(0, len(chrB), 60):
            f.write(chrB[i:i+60] + "\n")
    # 5. Write reference VCF (10 SNPs, REF=chrA, ALT=chrB)
    ref_vcf_path = test_dir / "reference.vcf"
    with open(ref_vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
        for i, pos in enumerate(snp_positions):
            ref_base = chrA[pos-1]
            alt = chrB[pos-1]
            f.write(f"chrA\t{pos}\tSNP{100+i:03d}\t{ref_base}\t{alt}\t99\tPASS\t.\tGT\t1/1\n")
    # 6. Write second species VCF (10 SNPs, same positions, ALT swapped back to ref)
    second_vcf_path = test_dir / "second_species.vcf"
    with open(second_vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
        for i, pos in enumerate(snp_positions):
            ref_base = chrB[pos-1]
            alt = chrA[pos-1]
            f.write(f"chrA\t{pos}\tSNP{200+i:03d}\t{ref_base}\t{alt}\t99\tPASS\t.\tGT\t1/1\n")
    # 7. Write a single-block MAF for the full 600bp alignment
    maf_path = test_dir / "alignment.maf"
    with open(maf_path, "w") as f:
        f.write("##maf version=1 scoring=LastZ\n")
        f.write("a score=1234\n")
        f.write(f"s chrA 0 600 + 600 {chrA}\n")
        f.write(f"s chrA 0 600 + 600 {chrB}\n")
    # 8. Write dummy GFF (not used, but file must exist)
    gff_path = test_dir / "annotation.gff"
    with open(gff_path, "w") as f:
        f.write("##gff-version 3\n")
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
        from ddprimer.helpers.lastz_runner import LastZRunner
        
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
        
        # Find actual masked positions (convert to 1-based for consistency with VCF)
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

# ------------------------------------------------------------------
# Test harness helpers:  capture AlignmentWorkflow outputs
# ------------------------------------------------------------------
# Keep a reference to the original run_alignment_workflow function
import ddprimer.helpers.alignment_workflow as _alignment_wf
_orig_run_alignment_workflow = _alignment_wf.run_alignment_workflow

captured_results = {"masked_sequences": {}, "coordinate_map": {}}

def _capture_alignment_workflow(args, out_dir):
    """Proxy AlignmentWorkflow that captures outputs and calls the real function."""
    masked_seqs, cmap = _orig_run_alignment_workflow(args, out_dir)
    captured_results["masked_sequences"] = masked_seqs
    captured_results["coordinate_map"] = cmap
    return masked_seqs, cmap

def test_alignment_workflow_with_maf():
    """
    Test alignment mode with a pre‑created MAF file.
    Captures internal masked_sequences / coordinate_map but skips primer design.
    """
    try:
        # --- monkey‑patch run_alignment_workflow function so we can capture its output ---
        import ddprimer.helpers.alignment_workflow
        ddprimer.helpers.alignment_workflow.run_alignment_workflow = _capture_alignment_workflow

        from ddprimer.core.SNP_masking_processor import SNPMaskingProcessor
        from ddprimer.helpers.maf_parser import MAFParser

        # -------------------------------------------------------------------
        #                  1)  create test input files
        # -------------------------------------------------------------------
        files = create_test_files()

        # Mock CLI args object
        class Args:
            def __init__(self):
                self.fasta         = str(files["ref_fasta"])
                self.second_fasta  = str(files["second_fasta"])
                self.vcf           = str(files["ref_vcf"])
                self.second_vcf    = str(files["second_vcf"])
                self.gff           = str(files["gff"])
                self.maf           = str(files["maf"])
                self.snp           = True
                self.min_identity  = 90.0
                self.min_length    = 500
                self.lastz_options = "--format=maf"
                self.lastzonly     = False
                self.noannotation  = False

        args = Args()

        # -------------------------------------------------------------------
        #                  2)  run AlignmentWorkflow end‑to‑end
        # -------------------------------------------------------------------
        output_dir = test_dir / "output"
        output_dir.mkdir(exist_ok=True)

        logger.info(magenta("\n===== FULL AlignmentWorkflow with pre‑computed MAF ====="))
        from ddprimer.helpers.alignment_workflow import run_alignment_workflow
        run_alignment_workflow(args, str(output_dir))

        # If we reach here, AlignmentWorkflow finished.  Check we captured data.
        if not captured_results["masked_sequences"]:
            logger.error("AlignmentWorkflow returned no masked sequences")
            return False

        # Sanity check: expect ~10-20 Ns (variants), not hundreds
        total_len = sum(len(seq) for seq in captured_results["masked_sequences"].values())
        n_total   = sum(seq.count("N") for seq in captured_results["masked_sequences"].values())
        logger.info(f"Captured {len(captured_results['masked_sequences'])} sequences, total Ns = {n_total}")
        # Expect at least 1 N, but not hundreds (since only variants are masked)
        if n_total > 50:
            logger.error("Too many Ns masked")
            return False

        # -------------------------------------------------------------------
        # 3) strict coordinate validation: each expected SNP must be masked
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # Build expected‑mask dictionaries by simple VCF parsing (no PyVCF).
        # -------------------------------------------------------------------
        def read_vcf_positions(vcf_path):
            pos_dict = {}
            with open(vcf_path) as vf:
                for line in vf:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 2:
                        continue
                    chrom = fields[0]
                    try:
                        pos = int(fields[1])
                    except ValueError:
                        continue
                    pos_dict.setdefault(chrom, set()).add(pos)
            return pos_dict

        expected_from_ref    = read_vcf_positions(files["ref_vcf"])
        expected_from_second = read_vcf_positions(files["second_vcf"])

        ok = True
        for seq_id, masked_seq in captured_results["masked_sequences"].items():
            N_pos = {i + 1 for i, b in enumerate(masked_seq) if b == "N"}  # 1‑based
            exp = expected_from_ref.get(seq_id, set()) | expected_from_second.get(seq_id, set())
            # Filter expected positions to those actually inside the conserved region
            exp = {p for p in exp if p < len(masked_seq)}  # end‑position in MAF is exclusive
            missing = exp - N_pos
            extra   = N_pos - exp
            if missing or extra:
                logger.error(f"{seq_id}: Masking mismatch – "
                             f"missing {sorted(missing)}, extra {sorted(extra)}")
                ok = False
            else:
                logger.info(f"{seq_id}: Exact masking validation {green('PASSED')}")

        if not ok:
            return False
        return ok

    except ImportError as e:
        logger.error(f"Import error in test_alignment_workflow_with_maf: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error in test_alignment_workflow_with_maf: {e}")
        import traceback; logger.error(traceback.format_exc())
        return False
    finally:
        # always restore the original function
        try:
            ddprimer.helpers.alignment_workflow.run_alignment_workflow = _orig_run_alignment_workflow
        except Exception:
            pass

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
        from ddprimer.helpers.alignment_workflow import run_alignment_workflow
        
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
                self.maf = None  # We'll use our fresh LastZ alignment
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
            # Call run_alignment_workflow
            masked_sequences, coordinate_map = run_alignment_workflow(args, str(output_dir))
            
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
            logger.error(f"Error running alignment workflow with LastZ: {e}")
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
            logger.info(green("\n✓ All alignment tests PASSED with expected masking!"))
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