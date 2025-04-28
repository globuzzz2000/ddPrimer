#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for SNP masking in direct mode, mimicking real-file I/O.
This version:
 1. Creates reference FASTA, VCF, CSV, and mapping.csv
 2. Uses real SNPMaskingProcessor.get_region_variants
 3. Stubs DirectMasking.find_location based on mapping.csv
 4. Captures masked_sequences and matching_status via a stubbed run_primer_design_workflow
 5. Runs direct_mode.run(args) unmodified
 6. Verifies both count and exact positions of Ns
"""
import os
import sys
import tempfile
import logging
from pathlib import Path
import shutil
import csv

# Configure test logger
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger('test')

# Silence heavy modules
import logging as _logging
_logging.getLogger('Bio').setLevel(_logging.CRITICAL)
_logging.getLogger('nupack').setLevel(_logging.CRITICAL)

# Create test directory
TEST_DIR = Path(tempfile.mkdtemp(prefix='ddprimer_test_realio_'))
logger.info(f'Created test directory: {TEST_DIR}')


def create_test_files(test_dir):
    """Create FASTA, VCF, CSV, and mapping.csv files."""
    chrom_seqs = {
        'chr1': 'TGTCCCACTTCTGCTCTAGCAGAGCACTTGTTTCGGTAACAGCTCCTTCGACAACTCGCACCGCAGAAGGATCAAGATACTGTTCCAGTAACTTAGTGAGCAGAGCTGACGGATACAAGTGGTGCTGAGTCTGACCCGATGTCTAACGTCGGGTTTAGTGGTTTGTCGTCTTTTAATGCGTCGAGTATGGTGTCTCCGCGCTCATCAGGTCAAGTTCAGGGCAAGAGGACTAGCTTCACGGGATACAATGCCTTGAAGACTAATGACCGCGGAAGAAGCTGCAAGAAGGAGATTAGCGGCGCCACGGCTATCTTCCGAAGTGATAACCGCCCGGAAAATTCCACATCCCGCGAGCAAAGGACCACTTCTTGATCCAATCACACCCGTTGGATGGATATGGGAC',
        'chr2': 'ACTCTTCCCTATCCCTTCTGATGACAATGGCGTAACAGGCGACACAGCACCAACGCCTAAAGCTCCAGGAATTGCAGCAGAAGCACCTTGGACAAAGCCCACATTGTTCGTTAAAGACTGATCTCCAAGTCCCACAAGACCACCACCTCTTATCCCAGGACTCCATGAAGGTGTCTCCAACGTTGACGACAAGAGCACCAGGACGAGGGCGAACGGTCTGCCAACTACCTGCGGCAAAGACCTCTAAGCCCACAACATCGAGTCGGAATCAGTGAAGGAAGTAGAAGAAGACACTGCGTCAGAGTTGATTTGCGTCGCCGTGTCATGTTGCCGCACCGGTGACTTATCCGACGACTTC',
        'chr3': 'GCATTGGAAGATCGAGGAGGGCGGCTAACACCGGTTTAAAGGATTTGACCATTTGGTGCATCCGTTTTATAGCACCGTGTCCAGCGGATTGAGCCCACGCGAGGTCAAAGCCATTGGAGAAGAACTTTCCGTGACCGGTAGTGATGAGGACGGATCCGGAGTATCAAGCATGTCGGATATGAACATGGAGAACCTTATGGAGGACTCTGTTGCTTTTAGGGTTCGGGCTAAACGTGGTTGCGCAACTCATCCCCGCAGCATTGCCGAGAGGGTAAACCTTCTTAAGGGACCTGCTGGTCTGGAATATCTGAAGAAGTCATTTTCCAGTCGGCATGGTTCCCCAGACCAAGCGTCTTCCTCGCTACCATTGCAGTGGTCTACTTCAGTGCCCCACCATTTGAGGTCGTGTCCAGCATTCTGCAGGGCGAAGAACAAGAGCACCCCCATGAACGCGGTCCCTGCATCGAGCGCTGCAGAGAGTACGTAATTGTACTTCTGCCACCATCTCTTGTGGTAATTGAACACAAAGTAGTTGAAGATGGTTCCTGTGACC'
    }
    # Write FASTA
    fasta = test_dir / 'reference.fasta'
    with open(fasta, 'w') as f:
        for chrom, seq in chrom_seqs.items():
            f.write(f'>{chrom}\n{seq}\n')

    # Write VCF
    vcf = test_dir / 'variants.vcf'
    vcf_records = [
        ('chr1', 10), ('chr1', 20), ('chr1', 30),
        ('chr2', 40), ('chr2', 85),
        ('chr3', 90), ('chr3', 150)
    ]
    with open(vcf, 'w') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos in vcf_records:
            f.write(f'{chrom}\t{pos}\t.\tN\tN\t.\tPASS\t.\n')

    # Prepare CSV and mapping
    csvp = test_dir / 'sequences.csv'
    mapping = test_dir / 'mapping.csv'
    rows = [['name', 'sequence']]
    map_rows = [['name', 'chrom', 'start', 'end']]
    # seq_from_chr1: genome 6-55
    rows.append(['seq_from_chr1', chrom_seqs['chr1'][5:55]])
    map_rows.append(['seq_from_chr1', 'chr1', '6', '55'])
    # seq_from_chr2: genome 56-100
    rows.append(['seq_from_chr2', chrom_seqs['chr2'][55:100]])
    map_rows.append(['seq_from_chr2', 'chr2', '56', '100'])
    # seq_from_chr3: genome 101-200
    rows.append(['seq_from_chr3', chrom_seqs['chr3'][100:200]])
    map_rows.append(['seq_from_chr3', 'chr3', '101', '200'])
    # non-match
    rows.append(['seq_no_match', 'G'*50])
    map_rows.append(['seq_no_match', 'none', '0', '0'])
    # Write files
    with open(csvp, 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    with open(mapping, 'w', newline='') as f:
        csv.writer(f).writerows(map_rows)
    return str(fasta), str(vcf), str(csvp), str(mapping)

# Generate test files
fasta_path, vcf_path, csv_path, map_path = create_test_files(TEST_DIR)

# Import production code
try:
    from ddprimer.modes import direct_mode
    from ddprimer.core import SNPMaskingProcessor
    from ddprimer.helpers import DirectMasking
    from ddprimer.utils import FileUtils
except ImportError as e:
    logger.error(f"Import failed: {e}")
    sys.exit(1)

# Load mapping.csv
mapping_dict = {}
with open(map_path) as mf:
    for r in csv.DictReader(mf):
        mapping_dict[r['name']] = (r['chrom'], int(r['start']), int(r['end']))

# Load sequences as pipeline will
test_sequences = FileUtils.load_sequences_from_table(csv_path)

# Prepare capture container
captured = {'matching_status': {}, 'masked_sequences': {}}

# Define stub for DirectMasking.find_location
orig_find = DirectMasking.find_location

def patched_find_location(sequence, ref_fasta, min_identity=90, min_coverage=90):
    for sid, seq in test_sequences.items():
        if seq == sequence:
            chrom, start, end = mapping_dict[sid]
            if chrom == 'none':
                return None, None, None, None
            return chrom, start, end, 100.0
    return None, None, None, None
DirectMasking.find_location = patched_find_location

# Define stub to capture workflow inputs
def capture_workflow(**kwargs):
    captured['matching_status'] = kwargs.get('matching_status', {})
    captured['masked_sequences'] = kwargs.get('masked_sequences', {})
    return True
# Monkey-patch workflow
orig_workflow = direct_mode.common.run_primer_design_workflow
direct_mode.common.run_primer_design_workflow = capture_workflow

# Build args
tmp = {}
class Args: pass
args = Args()
args.direct = csv_path
args.snp = True
args.fasta = fasta_path
args.vcf = vcf_path
args.output = str(TEST_DIR)
args.noannotation = False
args.cli = True
args.debug = True

# Run pipeline and tests
def main():
    if not direct_mode.run(args):
        logger.error('Pipeline run failed')
        return 1
    logger.info('Pipeline completed, verifying results...')

    sp = SNPMaskingProcessor()
    all_ok = True

    # Check matching statuses
    logger.info('=== Sequence matching ===')
    for name in test_sequences:
        expected = 'Failure' if name=='seq_no_match' else 'Success'
        got = captured['matching_status'].get(name)
        if got != expected:
            logger.error(f"{name}: expected {expected}, got {got}")
            all_ok = False
        else:
            logger.info(f"{name}: {got}")

    # Check SNP masking
    logger.info('=== SNP masking ===')
    for name, seq in test_sequences.items():
        chrom, start, end = mapping_dict[name]
        expected_positions = set()
        if chrom!='none':
            vars_in_region = sp.get_region_variants(vcf_path, chrom, start, end)
            for vp in vars_in_region:
                expected_positions.add(vp - start + 1)
        masked_seq = captured['masked_sequences'].get(name, '')
        # Visual comparison
        logger.info(f"\n{name}:")
        logger.info(f"  Original: {seq}")
        logger.info(f"  Masked:   {masked_seq}")
        # Build arrow indicator for expected mask positions
        arrow_line = ''.join('↑' if (i + 1) in expected_positions else ' ' for i in range(len(seq)))
        logger.info(f"  Mask pos: {arrow_line}")
        # Count
        cnt = masked_seq.count('N')
        if cnt != len(expected_positions):
            logger.error(f"{name}: expected {len(expected_positions)} Ns, got {cnt}")
            all_ok = False
        # Exact positions
        for pos in expected_positions:
            idx = pos-1
            val = masked_seq[idx] if 0<=idx<len(masked_seq) else None
            if val != 'N':
                logger.error(f"{name}: expected N at pos {pos}, got '{val}'")
                all_ok=False
        # Unexpected Ns
        extras = [i+1 for i,b in enumerate(masked_seq) if b=='N' and (i+1) not in expected_positions]
        if extras:
            logger.error(f"{name}: unexpected Ns at positions {extras}")
            all_ok=False

    return 0 if all_ok else 1

if __name__=='__main__':
    exit_code = main()
    # Restore
    DirectMasking.find_location = orig_find
    direct_mode.common.run_primer_design_workflow = orig_workflow
    # Cleanup
    shutil.rmtree(TEST_DIR)
    sys.exit(exit_code)