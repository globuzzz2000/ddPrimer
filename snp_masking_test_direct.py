#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated test script for SNP masking in direct mode with exact-position assertions,
while preserving the original visual log output.
"""
import os
import sys
import tempfile
import logging
from pathlib import Path
import shutil
import csv

# Configure test logger
logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')
logger = logging.getLogger("test")
logger.setLevel(logging.INFO)
# Clean default handlers
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
logger.addHandler(handler)

# Silence other loggers
override = logging.CRITICAL
logging.getLogger("nupack").setLevel(override)
logging.getLogger("ddPrimer").setLevel(override)
# Disable BLAST usage report
os.environ["BLAST_USAGE_REPORT"] = "false"
# Disable NUPACK verbose logging
os.environ["NUPACK_DISABLE_LOGGING"] = "1"

# Create test directory
test_dir = Path(tempfile.mkdtemp(prefix="ddprimer_test_direct_"))
logger.info(f"Created test directory: {test_dir}")


def create_test_files():
    """Create reference FASTA, VCF, and CSV for testing."""
    # Chromosome sequences
    chrom1_seq = ("TGTCCCACTTCTGCTCTAGCAGAGCACTTGTTTCGGTAACAGCTCCTTCGACAACTCGCACCGCAGAAGGATCAAGATACTGTTCCAGTAACTTAGTGAGCAGAGCTGACGGATACAAGTGGTGCTGAGTCTGACCCGATGTCTAACGTCGGGTTTAGTGGTTTGTCGTCTTTTAATGCGTCGAGTATGGTGTCTCCGCGCTCATCAGGTCAAGTTCAGGGCAAGAGGACTAGCTTCACGGGATACAATGCCTTGAAGACTAATGACCGCGGAAGAAGCTGCAAGAAGGAGATTAGCGGCGCCACGGCTATCTTCCGAAGTGATAACCGCCCGGAAAATTCCACATCCCGCGAGCAAAGGACCACTTCTTGATCCAATCACACCCGTTGGATGGATATGGGAC")
    chrom2_seq = ("ACTCTTCCCTATCCCTTCTGATGACAATGGCGTAACAGGCGACACAGCACCAACGCCTAAAGCTCCAGGAATTGCAGCAGAAGCACCTTGGACAAAGCCCACATTGTTCGTTAAAGACTGATCTCCAAGTCCCACAAGACCACCACCTCTTATCCCAGGACTCCATGAAGGTGTCTCCAACGTTGACGACAAGAGCACCAGGACGAGGGCGAACGGTCTGCCAACTACCTGCGGCAAAGACCTCTAAGCCCACAACATCGAGTCGGAATCAGTGAAGGAAGTAGAAGAAGACACTGCGTCAGAGTTGATTTGCGTCGCCGTVTCATGTTGCCGCACCGGTGACTTATCCGACGACTTC")
    chrom3_seq = ("GCATTGGAAGATCGAGGAGGGCGGCTAACACCGGTTTAAAGGATTTGACCATTTGGTGCATCCGTTTTATAGCACCGTGTCCAGCGGATTGAGCCCACGCGAGGTCAAAGCCATTGGAGAAGAACTTTCCGTGACCGGTAGTGATGAGGACGGATCCGGAGTATCAAGCATGTCGGATATGAACATGGAGAACCTTATGGAGGACTCTGTTGCTTTTAGGGTTCGGGCTAAACGTGGTTGCGCAACTCATCCCCGCAGCATTGCCGAGAGGGTAAACCTTCTTAAGGGACCTGCTGGTCTGGAATATCTGAAGAAGTCATTTTCCAGTCGGCATGGTTCCCCAGACCAAGCGTCTTCCTCGCTACCATTGCAGTGGTCTACTTCAGTGCCCCACCATTTGAGGTCGTGTCCAGCATTCTGCAGGGCGAAGAACAAGAGCACCCCCATGAACGCGGTCCCTGCATCGAGCGCTGCAGAGAGTACGTAATTGTACTTCTGCCACCATCTCTTGTGGTAATTGAACACAAAGTAGTTGAAGATGGTTCCTGTGACC")
    # Write FASTA
    ref = test_dir / "reference.fasta"
    with open(ref, 'w') as f:
        f.write(f">chr1\n{chrom1_seq}\n>chr2\n{chrom2_seq}\n>chr3\n{chrom3_seq}\n")
    # Write VCF
    vcf = test_dir / "variants.vcf"
    vcf_data = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	T	G	.	PASS	.
chr1	20	.	C	A	.	PASS	.
chr1	30	.	T	C	.	PASS	.
chr2	40	.	A	T	.	PASS	.
chr2	85	.	C	G	.	PASS	.
chr3	90	.	C	A	.	PASS	.
chr3	150	.	A	G	.	PASS	.
"""
    with open(vcf, 'w') as f:
        f.write(vcf_data)
    # Write CSV
    csvp = test_dir / "sequences.csv"
    rows = [["name", "sequence"],
            ["seq_from_chr1", chrom1_seq[5:55]],
            ["seq_from_chr2", chrom2_seq[55:100]],
            ["seq_from_chr3", chrom3_seq[100:200]],
            ["seq_no_match", "G"*50]]
    with open(csvp, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rows)
    return {'fasta': str(ref), 'vcf': str(vcf), 'csv': str(csvp)}

# Capture state for assertions
test_results = {'matching_status': {}, 'masked_sequences': {}, 'all_sequences': {}, 'adjusted_variants': {}}


def test_direct_mode_masking():
    try:
        from ddprimer.modes import direct_mode
        from ddprimer.core import SNPMaskingProcessor
    except ImportError as e:
        logger.error(f"Import error: {e}")
        return False

    files = create_test_files()
    logger.info(f"Created test files: FASTA={files['fasta']}, VCF={files['vcf']}, CSV={files['csv']}")
    logger.info("\n=== Testing Direct Mode SNP Masking ===")

    # Monkey-patch workflow to skip downstream
    orig_workflow = direct_mode.common.run_primer_design_workflow
    direct_mode.common.run_primer_design_workflow = lambda *a, **k: True
    orig_blast = direct_mode.find_sequence_location_with_blast
    orig_run = direct_mode.run

    def patched_run(args):
        sequences = direct_mode.FileUtils.load_sequences_from_table(args.direct)
        test_results['all_sequences'] = sequences.copy()
        sp = SNPMaskingProcessor()
        variants = sp.get_variant_positions(args.vcf)
        matching, masked = {}, {}

        for seq_id, seq in sequences.items():
            logger.info(f"\nProcessing sequence {seq_id} ({len(seq)} bp)")
            chrom, sstart, send, pid = orig_blast(seq, args.fasta)
            if chrom:
                matching[seq_id] = 'Success'
                logger.info(f"Matched {seq_id} to {chrom}:{sstart}-{send} with {pid}% identity")
                vars_chr = variants.get(chrom, set())
                logger.info(f"Variants for {chrom}: {vars_chr}")
                adjusted = set()
                for vp in vars_chr:
                    if sstart <= vp <= send:
                        adjusted.add(vp - sstart + 1)
                test_results['adjusted_variants'][seq_id] = adjusted
                logger.info(f"Adjusted variant positions: {adjusted}")
                if adjusted:
                    logger.info(f"Masking sequence {seq_id} with {len(adjusted)} variants")
                    masked_seq = sp.mask_variants(seq, adjusted)
                else:
                    masked_seq = seq
                # Visual output
                logger.info(f"Original: {seq}")
                logger.info(f"Masked:   {masked_seq}")
                arrow = ''.join('↑' if i+1 in adjusted else ' ' for i in range(len(seq)))
                logger.info(f"Mask pos: {arrow}")
                masked[seq_id] = masked_seq
            else:
                matching[seq_id] = 'Failure'
                logger.info(f"Could not match sequence {seq_id} to reference genome")
        test_results['matching_status'] = matching
        test_results['masked_sequences'] = masked
        return True

    direct_mode.run = patched_run

    class Args: pass
    args = Args()
    args.direct = files['csv']; args.snp = True; args.fasta = files['fasta']; args.vcf = files['vcf']
    args.output = str(test_dir); args.noannotation = False; args.cli = True; args.debug = True

    try:
        if not direct_mode.run(args):
            logger.error("Direct mode failed")
            return False
        logger.info("\n✓ Direct mode completed successfully")

        # Verifications
        all_ok = True
        logger.info("\n=== Verifying Sequence Matching ===")
        for seq_id, status in test_results['matching_status'].items():
            exp = 'Failure' if seq_id=='seq_no_match' else 'Success'
            if status != exp:
                logger.error(f"✗ {seq_id}: expected {exp}, got {status}")
                all_ok = False
            else:
                logger.info(f"✓ {seq_id}: {status}")

        logger.info("\n=== Verifying SNP Masking Counts ===")
        for seq_id, seq in test_results['masked_sequences'].items():
            cnt = seq.count('N')
            exp_cnt = len(test_results['adjusted_variants'].get(seq_id, []))
            if cnt != exp_cnt:
                logger.error(f"✗ {seq_id}: expected {exp_cnt} Ns, found {cnt}")
                all_ok = False
            else:
                logger.info(f"✓ {seq_id}: {cnt} Ns as expected")

        logger.info("\n=== Verifying Exact Mask Positions ===")
        for seq_id, seq in test_results['masked_sequences'].items():
            exp_pos = test_results['adjusted_variants'].get(seq_id, set())
            for pos in exp_pos:
                idx = pos-1
                if seq[idx] != 'N':
                    logger.error(f"✗ {seq_id}: expected 'N' at {pos}, found '{seq[idx]}'")
                    all_ok = False
                else:
                    logger.info(f"✓ {seq_id}: correct 'N' at {pos}")
            extras = [i+1 for i,b in enumerate(seq) if b=='N' and (i+1) not in exp_pos]
            if extras:
                logger.error(f"✗ {seq_id}: unexpected Ns at {extras}")
                all_ok = False

        if all_ok:
            logger.info("\n✓ TEST PASSED: SNP masking in direct mode is correct")
            return True
        else:
            logger.error("\n✗ TEST FAILED: Issues detected in SNP masking")
            return False

    finally:
        direct_mode.run = orig_run
        direct_mode.common.run_primer_design_workflow = orig_workflow
        direct_mode.find_sequence_location_with_blast = orig_blast
        logger.info(f"\nCleaning up test directory: {test_dir}")
        shutil.rmtree(test_dir)

if __name__ == "__main__":
    success = test_direct_mode_masking()
    sys.exit(0 if success else 1)
