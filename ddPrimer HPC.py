#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 21:30:15 2025

@author: jakob
"""

import os
import sys
import subprocess
import threading
import time
import re
import tempfile
import tkinter as tk
from tkinter import filedialog
from io import StringIO
import pandas as pd
from tqdm import tqdm
from multiprocessing import cpu_count
import concurrent.futures
import functools
import nupack

##############################################################################
#                           Variables
##############################################################################
# Design Parameters
MIN_SEGMENT_LENGTH = 100
RETAIN_TYPES = "gene"  ## gff filtering, "gene", "mRNA", "CDS", "exon",...
FILTER_MEANINGFUL_NAMES = True # will only use named genes from gff
COUNT_AMBIGUOUS_AS_MISMATCH = False
GENE_OVERLAP_MARGIN = 25
restriction_site = "GGCC"
PENALTY_MAX = 5.0
MAX_PRIMER_PAIRS_PER_SEGMENT = 3

# Tool paths
PARALLEL_LASTZ_PATH = "/cluster/home/jwimmer/ddPrimer/Tools/parallelLastZ/parallelLastz.pl"
PARALLEL_LASTZ_CONFIG = "/cluster/home/jwimmer/ddPrimer/Tools/parallelLastZ/parallelLastz.config"
PRIMER3_CORE_PATH = "/opt/homebrew/bin/primer3_core"
MULTIZ_PATH = "/cluster/home/jwimmer/ddPrimer/ToolsMultiz/multiz"
BLASTN_PATH = "/cluster/home/jwimmer/ddPrimer/Tools/Blast+/bin/blastn"
DB_PATH = "/cluster/home/jwimmer/ddPrimer/ToolsTair DB/TAIR10_db"

# Primer3 settings
PRIMER3_SETTINGS_DICTIONARY = {
    "P3_FILE_TYPE": "settings", "P3_FILE_ID": "User settings", "P3P_DEBUG_MODE": 0, "P3P_GB_ORIENTATION": "+",
    "P3P_PRIMER_NAME_ACRONYM_INTERNAL": "IN", "P3P_PRIMER_NAME_ACRONYM_LEFT": "F", "P3P_PRIMER_NAME_ACRONYM_RIGHT": "R",
    "P3P_PRIMER_NAME_ACRONYM_SPACER": "_", "PRIMER_ANNEALING_TEMP": 52.0, "PRIMER_DMSO_CONC": 0.0,
    "PRIMER_DMSO_FACTOR": 0.6, "PRIMER_DNA_CONC": 50.0, "PRIMER_DNTP_CONC": 0.8, "PRIMER_FIRST_BASE_INDEX": 1,
    "PRIMER_FORMAMIDE_CONC": 0.0, "PRIMER_GC_CLAMP": 1, "PRIMER_INSIDE_PENALTY": -1.0, "PRIMER_INTERNAL_DMSO_CONC": 0.0,
    "PRIMER_INTERNAL_DMSO_FACTOR": 0.6, "PRIMER_INTERNAL_DNA_CONC": 50.0, "PRIMER_INTERNAL_DNTP_CONC": 0.0,
    "PRIMER_INTERNAL_FORMAMIDE_CONC": 0.0, "PRIMER_INTERNAL_MAX_BOUND": 110.0, "PRIMER_INTERNAL_MAX_GC": 80.0,
    "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 47.0, "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB": 12.0, "PRIMER_INTERNAL_MAX_NS_ACCEPTED": 0,
    "PRIMER_INTERNAL_MAX_POLY_X": 4, "PRIMER_INTERNAL_MAX_SELF_ANY": 12.0, "PRIMER_INTERNAL_MAX_SELF_ANY_TH": 47.0,
    "PRIMER_INTERNAL_MAX_SELF_END": 12.0, "PRIMER_INTERNAL_MAX_SELF_END_TH": 47.0, "PRIMER_INTERNAL_MAX_SIZE": 27,
    "PRIMER_INTERNAL_MAX_TM": 70, "PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4,
    "PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7, "PRIMER_INTERNAL_MIN_BOUND": -10.0, "PRIMER_INTERNAL_MIN_GC": 30.0,
    "PRIMER_INTERNAL_MIN_QUALITY": 0, "PRIMER_INTERNAL_MIN_SIZE": 15, "PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE": -1,
    "PRIMER_INTERNAL_MIN_TM": 64, "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME": "hnnnn", "PRIMER_INTERNAL_OPT_BOUND": 97.0,
    "PRIMER_INTERNAL_OPT_GC_PERCENT": 50.0, "PRIMER_INTERNAL_OPT_SIZE": 20, "PRIMER_INTERNAL_OPT_TM": 65,
    "PRIMER_INTERNAL_SALT_DIVALENT": 0.0, "PRIMER_INTERNAL_SALT_MONOVALENT": 50.0, "PRIMER_LIBERAL_BASE": 1,
    "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS": 0, "PRIMER_LOWERCASE_MASKING": 0, "PRIMER_MAX_BOUND": 110.0,
    "PRIMER_MAX_END_GC": 3, "PRIMER_MAX_END_STABILITY": 9.0, "PRIMER_MAX_GC": 60.0, "PRIMER_MAX_HAIRPIN_TH": 47.0,
    "PRIMER_MAX_LIBRARY_MISPRIMING": 12.0, "PRIMER_MAX_NS_ACCEPTED": 0, "PRIMER_MAX_POLY_X": 4, "PRIMER_MAX_SELF_ANY": 8.0,
    "PRIMER_MAX_SELF_ANY_TH": 47.0, "PRIMER_MAX_SELF_END": 3.0, "PRIMER_MAX_SELF_END_TH": 47.0, "PRIMER_MAX_SIZE": 23,
    "PRIMER_MAX_TEMPLATE_MISPRIMING": 12.0, "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": 47.0, "PRIMER_MAX_TM": 65.0,
    "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4, "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7, "PRIMER_MIN_BOUND": -10.0,
    "PRIMER_MIN_END_QUALITY": 0, "PRIMER_MIN_GC": 50.0, "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": 3,
    "PRIMER_MIN_QUALITY": 0, "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": 3, "PRIMER_MIN_SIZE": 18, "PRIMER_MIN_TM": 50.0,
    "PRIMER_NUM_RETURN": 10, "PRIMER_OPT_BOUND": 97.0, "PRIMER_OPT_GC_PERCENT": 52.5, "PRIMER_OPT_SIZE": 20,
    "PRIMER_OPT_TM": 60, "PRIMER_OUTSIDE_PENALTY": 0.0, "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0, "PRIMER_PAIR_MAX_COMPL_END": 3.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 1, "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": 24.0, "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": 24.0,
    "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": 47.0, "PRIMER_PAIR_WT_COMPL_ANY": 0.0, "PRIMER_PAIR_WT_COMPL_ANY_TH": 0.0,
    "PRIMER_PAIR_WT_COMPL_END": 0.0, "PRIMER_PAIR_WT_COMPL_END_TH": 0.0, "PRIMER_PAIR_WT_DIFF_TM": 0.0,
    "PRIMER_PAIR_WT_IO_PENALTY": 0.0, "PRIMER_PAIR_WT_LIBRARY_MISPRIMING": 0.0, "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.0,
    "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.0, "PRIMER_PAIR_WT_PRODUCT_TM_GT": 0.0, "PRIMER_PAIR_WT_PRODUCT_TM_LT": 0.0,
    "PRIMER_PAIR_WT_PR_PENALTY": 1.0, "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING": 0.0, "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH": 0.0,
    "PRIMER_PICK_ANYWAY": 0, "PRIMER_PICK_INTERNAL_OLIGO": 1, "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PRODUCT_MAX_TM": 1000000.0, "PRIMER_PRODUCT_MIN_TM": -1000000.0, "PRIMER_PRODUCT_OPT_SIZE": 0,
    "PRIMER_PRODUCT_OPT_TM": 0.0, "PRIMER_PRODUCT_SIZE_RANGE": "90-200", "PRIMER_QUALITY_RANGE_MAX": 100,
    "PRIMER_QUALITY_RANGE_MIN": 0, "PRIMER_SALT_CORRECTIONS": 1, "PRIMER_SALT_DIVALENT": 3.8, "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT": 1, "PRIMER_SEQUENCING_ACCURACY": 20, "PRIMER_SEQUENCING_INTERVAL": 250,
    "PRIMER_SEQUENCING_LEAD": 50, "PRIMER_SEQUENCING_SPACING": 500, "PRIMER_TASK": "generic",
    "PRIMER_THERMODYNAMIC_PARAMETERS_PATH": "/opt/homebrew/Cellar/primer3/2.4.0/share/primer3/primer3_config/",
    "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1, "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": 0, "PRIMER_TM_FORMULA": 1,
    "PRIMER_WT_BOUND_GT": 0.0, "PRIMER_WT_BOUND_LT": 0.0, "PRIMER_WT_END_QUAL": 0.0, "PRIMER_WT_END_STABILITY": 0.0,
    "PRIMER_WT_GC_PERCENT_GT": 0.5, "PRIMER_WT_GC_PERCENT_LT": 0.5, "PRIMER_WT_HAIRPIN_TH": 0.0,
    "PRIMER_WT_LIBRARY_MISPRIMING": 0.0, "PRIMER_WT_NUM_NS": 0.0, "PRIMER_WT_POS_PENALTY": 0.0, "PRIMER_WT_SELF_ANY": 0.0,
    "PRIMER_WT_SELF_ANY_TH": 0.0, "PRIMER_WT_SELF_END": 0.0, "PRIMER_WT_SELF_END_TH": 0.0, "PRIMER_WT_SEQ_QUAL": 0.0,
    "PRIMER_WT_SIZE_GT": 1.0, "PRIMER_WT_SIZE_LT": 1.0, "PRIMER_WT_TEMPLATE_MISPRIMING": 0.0,
    "PRIMER_WT_TEMPLATE_MISPRIMING_TH": 0.0, "PRIMER_WT_TM_GT": 1.0, "PRIMER_WT_TM_LT": 1.0
}

def format_primer3_settings(settings_dict):
    return "\n".join(f"{key}={value}" for key, value in settings_dict.items()) + "\n"

PRIMER3_SETTINGS = format_primer3_settings(PRIMER3_SETTINGS_DICTIONARY)

# Blast+ parameters
BLAST_WORD_SIZE = 7
BLAST_EVALUE = 10
BLAST_MAX_TARGET_SEQS = 100
BLAST_REWARD = 2
BLAST_PENALTY = -3
BLAST_GAPOPEN = 5
BLAST_GAPEXTEND = 2
# best E-value must be at least BLAST_FILTER_FACTOR times smaller than the second-best E-value.
BLAST_FILTER_FACTOR = 100

# NUPACK parameters
NUPACK_TEMPERATURE = 37  # Celsius
NUPACK_SODIUM = 0.05     # Molar
NUPACK_MAGNESIUM = 0.0   # Molar

# Optimization Parameters
# Use 75% of available cores by default, but allow at least 1
NUM_PROCESSES = max(1, int(cpu_count() * 0.75))
# Batch size for multiprocessing tasks
BATCH_SIZE = 100
# MAF processing chunk size (number of lines to read at once)
MAF_CHUNK_SIZE = 10000  # Adjust based on your available memory

# Progress display
SHOW_PROGRESS = True

##############################################################################
#                          File Loader
##############################################################################

try:
    from Cocoa import NSApp
    macos_fix = True
except ImportError:
    macos_fix = False
# Detect if running in a headless (GUI-less) environment
euler_fix = "euler" in os.uname().nodename.lower()  # Euler cluster detection
headless_fix = (
    not sys.stdout.isatty() or  # Running in a non-interactive terminal
    "DISPLAY" not in os.environ or  # No X server available (Linux)
    "WSL_DISTRO_NAME" in os.environ or  # Windows Subsystem for Linux
    "SSH_CONNECTION" in os.environ  # Remote SSH session
)
def _tk_get_file(prompt, filetypes):
    """ macOS file dialog wrapper. """
    os.environ['TK_SILENCE_DEPRECATION'] = '1'
    root = tk.Tk()
    root.withdraw()
    if macos_fix:
        NSApp().setActivationPolicy_(1)
    file_path = filedialog.askopenfilename(title=prompt, filetypes=filetypes)
    root.destroy()
    return file_path

def _tk_get_files(prompt, filetypes):
    """ macOS multi-file dialog wrapper. """
    os.environ['TK_SILENCE_DEPRECATION'] = '1'
    root = tk.Tk()
    root.withdraw()
    if macos_fix:
        NSApp().setActivationPolicy_(1)
    file_paths = filedialog.askopenfilenames(title=prompt, filetypes=filetypes)
    root.destroy()
    return file_paths

def get_file(prompt, filetypes):
    """
    Show a file dialog on Mac or prompt for input on Euler.
    """
    if macos_fix and not headless_fix:
        return _tk_get_file(prompt, filetypes)
    else:
        file_path = input(f"{prompt} (Enter path manually): ").strip()
        if not file_path:
            print(f"[ERROR] {prompt} - No file provided. Exiting.")
            sys.exit(1)
        return file_path

def get_files(prompt, filetypes):
    """
    Show a file dialog on Mac or prompt for input on Euler.
    """
    if macos_fix and not headless_fix:
        return _tk_get_files(prompt, filetypes)
    else:
        paths = input(f"{prompt} (Enter comma-separated paths): ").strip()
        file_paths = paths.split(",") if paths else []
        if not file_paths:
            print(f"[ERROR] {prompt} - No files provided. Exiting.")
            sys.exit(1)
        return file_paths

def load_fasta(filepath):
    """
    Load sequences from a FASTA file into a dict: {header_without_gt: sequence}.
    Optimized for memory efficiency.
    """
    sequences = {}
    name = None
    seq_chunks = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous
                if name:
                    sequences[name] = "".join(seq_chunks).upper()
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # Last record
        if name:
            sequences[name] = "".join(seq_chunks).upper()
    
    return sequences

##############################################################################
#                          GFF Filtering
##############################################################################

def parse_gff_attributes(attribute_str):
    """
    Convert GFF attribute string (key1=val1;key2=val2) -> dict.
    Keys forced to lower case.
    """
    attr_dict = {}
    for attr in attribute_str.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attr_dict[key.strip().lower()] = value.strip()
    return attr_dict

# Patterns for general placeholder detection across species
PLACEHOLDER_PATTERNS = [
    re.compile(r'^AT\dG\d{5}$', re.IGNORECASE),   # Arabidopsis
    re.compile(r'^LOC_.*$', re.IGNORECASE),       # NCBI-style, rice, etc.
    re.compile(r'^GRMZM.*$', re.IGNORECASE),      # Maize
    re.compile(r'^ENS.*$', re.IGNORECASE)         # Ensembl IDs (humans, etc.)
]

def is_meaningful_name(name, gene_id=None, locus_tag=None):
    """
    Determine whether the 'name' is a meaningful gene symbol,
    excluding placeholders, technical locus tags, or accessions.
    """
    if not name:
        return False

    name_lower = name.lower()

    if gene_id and name_lower == gene_id.lower():
        return False

    if locus_tag and name_lower == locus_tag.lower():
        return False

    for pattern in PLACEHOLDER_PATTERNS:
        if pattern.match(name):
            return False

    return True

def process_gff_chunk(chunk):
    """
    Process a chunk of GFF file lines.
    Used for parallel processing.
    """
    chunk_genes = []
    for line in chunk:
        if line.startswith('#'):
            continue

        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue

        ftype = parts[2].lower()
        if ftype not in RETAIN_TYPES:
            continue

        seqname, _, _, start, end, _, strand, _, attributes = parts
        attr_dict = parse_gff_attributes(attributes)

        name = attr_dict.get('name')
        gene_id = attr_dict.get('id')
        locus_tag = attr_dict.get('locus_tag')

        if FILTER_MEANINGFUL_NAMES and not is_meaningful_name(name, gene_id, locus_tag):
            continue

        try:
            chunk_genes.append({
                "chr": seqname,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "id": name
            })
        except ValueError:
            pass
    
    return chunk_genes

def load_genes_from_gff(gff_path):
    """
    Return a list of gene dicts: {chr, start, end, strand, id}.
    Uses global RETAIN_TYPES and FILTER_MEANINGFUL_NAMES for filtering.
    Optimized version that processes the file in parallel chunks.
    """
    # Read all lines from the file
    with open(gff_path, 'r') as f:
        all_lines = f.readlines()
    
    # Calculate chunk size for parallel processing
    chunk_size = max(1, len(all_lines) // NUM_PROCESSES)
    chunks = [all_lines[i:i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
    
    # Process chunks in parallel
    genes = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_PROCESSES) as executor:
        futures = [executor.submit(process_gff_chunk, chunk) for chunk in chunks]
        
        if SHOW_PROGRESS:
            for future in tqdm(concurrent.futures.as_completed(futures), 
                              total=len(futures), 
                              desc="Processing GFF file"):
                genes.extend(future.result())
        else:
            for future in concurrent.futures.as_completed(futures):
                genes.extend(future.result())

    print(f"[INFO] Loaded {len(genes)} entries from GFF.")
    return genes

##############################################################################
#                           Progression Spinner
##############################################################################

def spinner_task(stop_event):
    """
    Print a rotating spinner while stop_event is not set.
    """
    spinner = ['|', '/', '-', '\\']
    idx = 0
    sys.stdout.write(spinner[idx])
    sys.stdout.flush()
    while not stop_event.is_set():
        time.sleep(0.1)
        idx = (idx + 1) % 4
        sys.stdout.write("\b" + spinner[idx])
        sys.stdout.flush()
    sys.stdout.write("\b")
    sys.stdout.flush()

##############################################################################
#                           MAF Alignment Processing
##############################################################################

def parse_maf_chunk(chunk_lines):
    """
    Parse a chunk of MAF format lines into alignment records.
    """
    alignments = []
    current_alignment = {}

    for line in chunk_lines:
        line = line.strip()
        if line.startswith('a '):
            # Start of a new alignment block
            if current_alignment and 'ref_name' in current_alignment and 'query_name' in current_alignment:
                alignments.append(current_alignment)
            current_alignment = {}
        elif line.startswith('s '):
            fields = line.split()
            if len(fields) >= 7:  # Make sure we have enough fields
                seq_name = fields[1]
                start = int(fields[2])
                sequence = fields[6]

                if 'ref_name' not in current_alignment:
                    # First sequence in block = reference
                    current_alignment['ref_name'] = seq_name
                    current_alignment['ref_start'] = start + 1  # Convert to 1-based
                    current_alignment['ref_sequence'] = sequence
                else:
                    # Second sequence in block = query
                    current_alignment['query_name'] = seq_name
                    current_alignment['query_start'] = start + 1
                    current_alignment['query_sequence'] = sequence

    # Don't forget the last alignment
    if current_alignment and 'ref_name' in current_alignment and 'query_name' in current_alignment:
        alignments.append(current_alignment)

    return alignments

def process_maf_file_in_chunks(maf_path, chunk_processor, chunk_size=10000):
    """
    Process a MAF file in chunks to avoid loading it all into memory.
    The chunk_processor function is called for each chunk of alignments.
    """
    chunk_lines = []
    segments = []
    
    with open(maf_path, 'r') as maf_file:
        for line in maf_file:
            chunk_lines.append(line)
            
            # Process chunk when it reaches the desired size or contains complete alignment blocks
            if len(chunk_lines) >= chunk_size:
                # Ensure we don't split in the middle of an alignment block
                # Find the last occurrence of a line starting with 'a '
                last_a_idx = None
                for i in range(len(chunk_lines) - 1, -1, -1):
                    if chunk_lines[i].startswith('a '):
                        last_a_idx = i
                        break
                
                # If we found an 'a ' line and it's not the last line, process up to that point
                if last_a_idx is not None and last_a_idx < len(chunk_lines) - 1:
                    process_chunk = chunk_lines[:last_a_idx]
                    chunk_lines = chunk_lines[last_a_idx:]
                else:
                    process_chunk = chunk_lines
                    chunk_lines = []
                
                # Parse and process this chunk
                alignments = parse_maf_chunk(process_chunk)
                new_segments = chunk_processor(alignments)
                segments.extend(new_segments)
    
    # Process any remaining lines
    if chunk_lines:
        alignments = parse_maf_chunk(chunk_lines)
        new_segments = chunk_processor(alignments)
        segments.extend(new_segments)
    
    return segments

##############################################################################
#                           Mismatch identifier
##############################################################################

def bases_match(ref_base, qry_base):
    """
    Checks if two bases match, with optional ambiguity handling.
    If COUNT_AMBIGUOUS_AS_MISMATCH is True, any ambiguous base is treated as mismatch.
    """
    if ref_base == qry_base:
        return True
    if COUNT_AMBIGUOUS_AS_MISMATCH:
        return False

    # Basic IUPAC handling
    ref_base = ref_base.upper()
    qry_base = qry_base.upper()

    ambiguity_map = {
        'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
        'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
    }

    ref_options = ambiguity_map.get(ref_base, {ref_base})
    qry_options = ambiguity_map.get(qry_base, {qry_base})

    # They match if their sets of possible bases overlap
    return not ref_options.isdisjoint(qry_options)

def split_on_mismatches(ref_seq, qry_seq):
    """
    Split ref_seq into perfect-match segments based on mismatches and indels.
    Return only segments of length >= MIN_SEGMENT_LENGTH.
    """
    segments = []
    curr_segment = []

    for r, q in zip(ref_seq, qry_seq):
        if r == '-' or q == '-' or not bases_match(r, q):
            if len(curr_segment) >= MIN_SEGMENT_LENGTH:
                segments.append("".join(curr_segment))
            curr_segment = []
        else:
            curr_segment.append(r)

    if len(curr_segment) >= MIN_SEGMENT_LENGTH:
        segments.append("".join(curr_segment))

    return segments

##############################################################################
#                     Segment Processing
##############################################################################

def collect_segments_from_alignments(alignments, genes):
    """
    Process a batch of alignments and collect segments meeting criteria.
    This function is used as the chunk processor for MAF file processing.
    """
    chunk_segments = []
    
    for aln in alignments:
        chr_name = aln['ref_name']
        ref_start = aln['ref_start']
        ref_seq = aln['ref_sequence']
        qry_seq = aln['query_sequence']

        # ref_start is the 1-based coordinate of the first base in ref_seq
        # so the ith base of ref_seq is at genome position ref_start + i - 1

        for gene in genes:
            if gene['chr'] != chr_name:
                continue

            # Overlap range
            gene_start = gene['start'] - GENE_OVERLAP_MARGIN
            gene_end = gene['end'] + GENE_OVERLAP_MARGIN

            # If this alignment doesn't overlap the gene region, skip
            aln_ref_end = ref_start + len(ref_seq.replace('-', '')) - 1
            if aln_ref_end < gene_start or ref_start > gene_end:
                continue

            # Overlapping portion
            overlap_start = max(ref_start, gene_start)
            overlap_end = min(aln_ref_end, gene_end)

            # Indices relative to the alignment's first base (0-based)
            aln_offset_start = overlap_start - ref_start
            aln_offset_end = aln_offset_start + (overlap_end - overlap_start) + 1

            # Make sure we don't go out of bounds
            if aln_offset_start < 0:
                aln_offset_start = 0
            if aln_offset_end > len(ref_seq):
                aln_offset_end = len(ref_seq)

            # Extract the overlapping fragments
            ref_fragment = ref_seq[aln_offset_start:aln_offset_end]
            qry_fragment = qry_seq[aln_offset_start:aln_offset_end]

            # Split into perfect match segments
            perfect_segments = split_on_mismatches(ref_fragment, qry_fragment)

            # Further split perfect segments on restriction sites
            for segment in perfect_segments:
                # Use re.split to handle overlapping restriction sites
                subfragments = re.split(f'(?={restriction_site})', segment)
                for i, frag in enumerate(subfragments):
                    # Remove the restriction site prefix if present (except for first fragment)
                    if i > 0 and frag.startswith(restriction_site):
                        frag = frag[len(restriction_site):]
                    
                    if len(frag) >= MIN_SEGMENT_LENGTH:
                        n_existing = sum(1 for s in chunk_segments if s['gene'] == gene['id'])
                        chunk_segments.append({
                            "gene": gene['id'],
                            "index": n_existing + 1,
                            "sequence": frag
                        })
    
    return chunk_segments

##############################################################################
#                        Sequence Validation
##############################################################################

def exact_match_in_genome(sequence, genome_dict):
    """
    Check if a given sequence is found exactly in any of the sequences in genome_dict.
    """
    return any(sequence in seq for seq in genome_dict.values())

def validate_sequence(seq_info, ref_genome, qry_genome, extra_genomes=[]):
    """
    Validate a single segment sequence against multiple genomes.
    """
    seq = seq_info["sequence"]
    return (
        exact_match_in_genome(seq, ref_genome) and 
        exact_match_in_genome(seq, qry_genome) and
        all(exact_match_in_genome(seq, g) for g in extra_genomes)
    )

def validate_sequences_batch(batch, ref_genome, qry_genome, extra_genomes=[]):
    """
    Process a batch of sequences in parallel.
    """
    results = []
    for seg in batch:
        if validate_sequence(seg, ref_genome, qry_genome, extra_genomes):
            results.append(seg)
    return results

##############################################################################
#                          BLAST Functions
##############################################################################

def blast_short_seq(seq, db=f'"{DB_PATH}"'):
    """Runs BLASTn for short sequences and returns the two lowest e-values separately."""
    if not seq or not isinstance(seq, str) or not seq.strip():
        return None, None  # Ensuring consistency in output

    # Use a temporary file for the query sequence
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_query:
        tmp_query.write(f">seq\n{seq}\n")
        tmp_query.flush()
        tmp_filename = tmp_query.name

    # Run BLASTn
    try:
        result = subprocess.run(
            [
                BLASTN_PATH,
                "-task", "blastn-short",
                "-db", db,
                "-query", tmp_filename,
                "-word_size", str(BLAST_WORD_SIZE),
                "-evalue", str(BLAST_EVALUE),
                "-reward", str(BLAST_REWARD),
                "-penalty", str(BLAST_PENALTY),
                "-gapopen", str(BLAST_GAPOPEN),
                "-gapextend", str(BLAST_GAPEXTEND),
                "-max_target_seqs", str(BLAST_MAX_TARGET_SEQS),
                "-outfmt", "6 evalue"
            ],
            text=True,
            capture_output=True
        )
    finally:
        # Ensure we clean up even if an exception occurs
        try:
            os.remove(tmp_filename)  # Clean up temp file
        except OSError:
            pass

    # If BLAST fails, print the error and return None values
    if result.returncode != 0:
        print(f"[ERROR] BLAST failed: {result.stderr}")
        return None, None

    # Parse BLAST output correctly
    try:
        evalues = sorted([float(line.strip()) for line in result.stdout.strip().split("\n") if line.strip()])
    except ValueError:
        evalues = []

    if not evalues:
        return None, None

    # Return the two lowest e-values as separate values
    best = evalues[0] if len(evalues) > 0 else None
    second = evalues[1] if len(evalues) > 1 else None

    return best, second  # Separate columns

def process_blast_batch(batch_data):
    """
    Process a batch of sequences for BLAST in parallel.
    """
    batch, col_name = batch_data
    results = []
    
    for seq in batch:
        if pd.notnull(seq):
            blast1, blast2 = blast_short_seq(seq)
            results.append((blast1, blast2))
        else:
            results.append((None, None))
    
    return results

##############################################################################
#                        Primer3 Processing
##############################################################################

def run_primer3_batch(input_blocks):
    """
    Run primer3_core on a batch of input blocks.
    Returns the stdout output.
    """
    primer3_input_str = "".join(input_blocks)

    try:
        primer3_proc = subprocess.Popen(
            [PRIMER3_CORE_PATH],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True
        )
        stdout_data, stderr_data = primer3_proc.communicate(input=primer3_input_str)
        
        if primer3_proc.returncode != 0:
            print("[ERROR] Primer3 failed:")
            print(stderr_data)
            return ""
        
        if stderr_data and SHOW_PROGRESS:
            print("\n[WARNING] Primer3 warnings or info for batch:")
            print(stderr_data)
            
        return stdout_data
    except Exception as e:
        print(f"[ERROR] Failed to run Primer3: {e}")
        return ""

def parse_primer3_batch(stdout_data):
    """
    Parse primer3 output for a batch.
    Returns a list of records.
    """
    records = []
    current_id = None
    sequence_template = ""
    pairs = []
    
    def get_amplicon(seq, left_start, left_len, right_start, right_len):
        """
        Return the substring from the 5' end of the forward primer
        to the 3' end of the reverse primer, 1-based coords in the template.
        """
        if None in (left_start, left_len, right_start, right_len):
            return ""
        amp_start = left_start
        amp_end = right_start
        if amp_start < 1 or amp_end > len(seq) or amp_start > amp_end:
            return ""
        return seq[amp_start - 1 : amp_end]

    # Helper function: finalize the current record
    def finalize_record():
        """
        Save up to 3 pairs that meet penalty requirements.
        """
        nonlocal current_id, pairs, sequence_template

        if not current_id or not pairs:
            return

        acceptable = []
        for p in pairs:
            pair_penalty = p.get("pair_penalty", 999)
            probe_penalty = p.get("internal_penalty", 999)  # if no probe
            if pair_penalty <= PENALTY_MAX and probe_penalty <= PENALTY_MAX:
                acceptable.append(p)

        # take up to 3
        acceptable = acceptable[:MAX_PRIMER_PAIRS_PER_SEGMENT]
        for p in acceptable:
            left_seq = p.get("left_sequence", "")
            right_seq = p.get("right_sequence", "")
            probe_seq = p.get("internal_sequence", "")
            ls, ll = p.get("left_start"), p.get("left_len")
            rs, rl = p.get("right_start"), p.get("right_len")

            ampseq = get_amplicon(sequence_template, ls, ll, rs, rl)
            rec = {
                "Sequence": current_id,
                "Index": p["idx"],
                "Template": sequence_template,
                "Primer F": left_seq,
                "Tm F": p.get("left_tm", None),
                "Penalty F": p.get("left_penalty", None),
                "Primer R": right_seq,
                "Tm R": p.get("right_tm", None),
                "Penalty R": p.get("right_penalty", None),
                "Probe": probe_seq,
                "Probe Tm": p.get("internal_tm", None),
                "Probe Penalty": p.get("internal_penalty", None),
                "Pair Penalty": p.get("pair_penalty", None),
                "Amplicon": ampseq,
                "Length": p.get("product_size", None)
            }
            records.append(rec)

        current_id = None
        sequence_template = ""
        pairs = []

    lines = stdout_data.splitlines()
    for line in lines:
        line = line.strip()
        if line.startswith("SEQUENCE_ID="):
            # finalize any existing record first
            finalize_record()
            current_id = line.split("=", 1)[1]
            sequence_template = ""
            pairs = []

        elif line.startswith("SEQUENCE_TEMPLATE="):
            sequence_template = line.split("=", 1)[1].upper()

        elif re.match(r'^PRIMER_PAIR_(\d+)_PENALTY=', line):
            match = re.match(r'^PRIMER_PAIR_(\d+)_PENALTY=(.*)', line)
            idx, val = int(match.group(1)), float(match.group(2))
            pair = next((p for p in pairs if p['idx'] == idx), None)
            if not pair:
                pair = {"idx": idx}
                pairs.append(pair)
            pair["pair_penalty"] = val

        elif re.match(r'^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=', line):
            match = re.match(r'^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(.*)', line)
            idx, val = int(match.group(1)), int(match.group(2))
            pair = next((p for p in pairs if p['idx'] == idx), None)
            if pair:
                pair["product_size"] = val

        elif re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)_(SEQUENCE|TM|PENALTY)=', line):
            match = re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)_(SEQUENCE|TM|PENALTY)=(.*)', line)
            side, idx, attr, val = match.groups()
            idx = int(idx)
            pair = next((p for p in pairs if p['idx'] == idx), None)
            if not pair:
                pair = {"idx": idx}
                pairs.append(pair)
            attr_key = f"{side.lower()}_{attr.lower()}"
            if attr in ["TM", "PENALTY"]:
                pair[attr_key] = float(val)
            else:
                pair[attr_key] = val.upper()

        elif re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)=(\d+),(\d+)', line):
            match = re.match(r'^PRIMER_(LEFT|RIGHT|INTERNAL)_(\d+)=(\d+),(\d+)', line)
            side, idx, start, length = match.groups()
            idx = int(idx)
            start = int(start) + 1  # Convert to 1-based for consistency
            length = int(length)
            pair = next((p for p in pairs if p['idx'] == idx), None)
            if not pair:
                pair = {"idx": idx}
                pairs.append(pair)
            pair[f"{side.lower()}_start"] = start
            pair[f"{side.lower()}_len"] = length

        elif line == "=":
            # end of current record
            finalize_record()

    # If the last record never ended with '='
    finalize_record()
    
    return records

##############################################################################
#                     Amplicon Location
##############################################################################

def find_amplicon_locations_batch(batch_data):
    """
    Find the locations of a batch of amplicons in parallel.
    """
    amplicon_seqs, genome_dict = batch_data
    results = []
    
    for seq in amplicon_seqs:
        if pd.notnull(seq):
            for chrom, ref_seq in genome_dict.items():
                start_pos = ref_seq.find(seq)
                if start_pos != -1:
                    results.append((chrom, start_pos + 1))  # Convert to 1-based
                    break
            else:
                results.append((None, None))
        else:
            results.append((None, None))
    
    return results

##############################################################################
#                     NUPACK ΔG Calculation
##############################################################################

def calc_deltaG(seq):
    """
    Calculate the minimum free energy of a sequence using NUPACK.
    """
    if not isinstance(seq, str) or seq == "":
        return None
        
    dna_pattern = re.compile(r'^[ACGTNacgtn]+$')
    if not dna_pattern.match(seq):
        return None
        
    try:
        model = nupack.Model(
            material='dna',
            celsius=NUPACK_TEMPERATURE,
            sodium=NUPACK_SODIUM,
            magnesium=NUPACK_MAGNESIUM
        )
        result = nupack.mfe(seq, model=model)
        if result:
            # result[0].energy is the MFE structure's free energy
            return result[0].energy
    except Exception as e:
        if SHOW_PROGRESS:
            preview = seq[:50] + ("..." if len(seq) > 50 else "")
            print(f"[WARNING] NUPACK error for sequence {preview}: {e}")
        return None
    return None

def calc_deltaG_batch(seqs):
    """
    Calculate ΔG for a batch of sequences in parallel.
    """
    results = []
    for seq in seqs:
        if pd.notnull(seq):
            results.append(calc_deltaG(seq))
        else:
            results.append(None)
    return results

##############################################################################
#                          Utility Functions
##############################################################################

def chunks(lst, n):
    """
    Split a list into chunks of size n.
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def flatten_list(list_of_lists):
    """
    Flatten a list of lists into a single list.
    """
    return [item for sublist in list_of_lists for item in sublist]

def create_alignment_directory(reference_path):
    """
    Create an 'Alignments' directory in the same folder as the reference genome.
    Return the path to the alignments directory.
    """
    # Get the directory containing the reference genome
    ref_dir = os.path.dirname(os.path.abspath(reference_path))
    
    # Define alignments directory path
    alignments_dir = os.path.join(ref_dir, "Alignments")
    
    # Create the directory if it doesn't exist
    if not os.path.exists(alignments_dir):
        os.makedirs(alignments_dir)
        print(f"[INFO] Created Alignments directory: {alignments_dir}")
    
    return alignments_dir

def has_disallowed_repeats(seq):
    """
    Check for disallowed repeats in a sequence.
    """
    if not isinstance(seq, str):
        return True
    return "CCCC" in seq or "GGGG" in seq

def calculate_gc(seq):
    """
    Calculate GC content of a sequence.
    """
    if not seq or not isinstance(seq, str):
        return 0
    seq = seq.upper()
    gc_count = sum(1 for base in seq if base in "GC")
    return (gc_count / len(seq)) * 100 if seq else 0

def passes_blast_filter(row, col):
    """
    Check if a sequence passes the BLAST specificity filter.
    """
    best = row[f"{col} BLAST1"]
    second = row[f"{col} BLAST2"]
    # If best is None, no hits => discard or keep? (Often you'd discard, but user preference may vary.)
    if pd.isna(best):
        return False
    # If no second best, it's effectively unique
    if pd.isna(second):
        return True
    # best must be at least FILTER_FACTOR times smaller => best * FACTOR <= second
    return best * BLAST_FILTER_FACTOR <= second

##############################################################################
#                           Main Program
##############################################################################

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    print("\n========== ddPrimer: Combined Primer Design Pipeline ==========")
    print("\n[INPUT] Do you want to compare two or more genomes?")
    print("1. Two genomes (skip MultiZ)")
    print("2. More than two genomes (use MultiZ)")
    comparison_choice = input("Enter 1 or 2: ").strip()
    
    use_multiz = comparison_choice == "2"
    num_genomes = 2 if not use_multiz else int(input("[INPUT] How many genomes? (≥3): ").strip())

    print("\n[INPUT] Do you want to run LastZ, or do you already have .maf files?")
    print("1. Run LastZ")
    print("2. Use existing .maf files")
    maf_choice = input("Enter 1 or 2: ").strip()

    maf_files = []

    if maf_choice == "1":
        print("[INFO] Proceeding with LastZ alignments...")
    else:
        print("[INFO] Using existing MAF files...")
        print("[INPUT] Select MAF files for analysis")
        maf_files = get_files("Select MAF files", [("MAF Files", "*.maf"), ("All Files", "*.*")])

    ##############################################################################
    #                          1. File Import
    ##############################################################################

    if euler_fix:
        print("[INFO] Running in Euler cluster mode - using manual path input.")
    print("[INPUT] Select the Reference FASTA File (reference genome)")
    ref_path = get_file("Select Reference FASTA File", [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])

    print("[INPUT] Select the Query FASTA File (query genome)")
    qry_path = get_file("Select Query FASTA File", [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])

    # Determine additional genomes beyond the first two (already provided)
    extra_genome_paths = []
    extra_genomes = []
    num_extra = max(0, num_genomes - 2)  # Ensure we only ask for additional genomes beyond the first two
    if num_extra > 0:
        print(f"[INPUT] Select {num_extra} additional genomes for filtering:")
        for i in range(num_extra):
            extra_path = get_file(f"Select {i+1}. Additional FASTA File", [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])
            extra_genome_paths.append(extra_path)
            extra_genomes.append(load_fasta(extra_path))

    print("[INPUT] Select the GFF3 Annotation File")
    gff_path = get_file("Select GFF3 Annotation File", [("GFF3 Files", "*.gff3 *.gff"), ("All Files", "*.*")])

    ##############################################################################
    #                          2. Run LastZ
    ##############################################################################

    genes = load_genes_from_gff(gff_path)

    # Store output MAF files
    if maf_choice == "1":
        print("[INFO] Running LastZ alignments...")
        lastz_output_files = []
        
        # Create the Alignments directory
        alignments_dir = create_alignment_directory(ref_path)
        
        # Use already selected genomes: reference, query, and any additional genomes
        genome_paths = [ref_path, qry_path] + extra_genome_paths
        
        # Run ParallelLastZ pairwise for all genome combinations using all possible pairs
        from itertools import combinations
        for ref_genome_path, qry_genome_path in combinations(genome_paths, 2):
            # Extract base names for reference and query genomes
            ref_basename = os.path.splitext(os.path.basename(ref_genome_path))[0]
            qry_basename = os.path.splitext(os.path.basename(qry_genome_path))[0]
            
            # Create a specific filename based on your desired naming convention
            alignment_name = f"{ref_basename}_vs_{qry_basename}"
            
            # Set the output file path in the Alignments directory
            lastz_output_file = os.path.join(alignments_dir, f"{alignment_name}.maf")
            lastz_output_files.append(lastz_output_file)
            
            # Create a temporary directory for ParallelLastZ
            temp_output_dir = os.path.join(alignments_dir, f"temp_{alignment_name}")
            os.makedirs(temp_output_dir, exist_ok=True)
            
            # Construct the parallelLastz command
            lastz_cmd = [
                "perl", PARALLEL_LASTZ_PATH,
                "--qfile", qry_genome_path,
                "--tfile", ref_genome_path,
                "--cfile", PARALLEL_LASTZ_CONFIG,
                "--speedup", str(NUM_PROCESSES),  # Use the same number of processes as elsewhere in the pipeline
                "--length", str(MIN_SEGMENT_LENGTH),  # Use the minimum segment length
                "--output", temp_output_dir
            ]
            
            print(f"\n[INFO] Running LastZ for {ref_basename} vs {qry_basename}... ")
            stop_event = threading.Event()
            thread = threading.Thread(target=spinner_task, args=(stop_event,))
            thread.start()
            
            proc = subprocess.Popen(lastz_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            _, lastz_error = proc.communicate()
            
            stop_event.set()
            thread.join()
            
            if proc.returncode != 0:
                print(f"\n[ERROR] ParallelLastZ failed with error:\n{lastz_error}")
                sys.exit(1)
                
            # Find the MAF file in the temporary output directory and move it to the desired location
            maf_files_in_dir = [os.path.join(temp_output_dir, f) for f in os.listdir(temp_output_dir) if f.endswith('.maf')]
            if maf_files_in_dir:
                # Copy/move the first MAF file to our target location
                # Using subprocess.run instead of shutil to better handle larger files
                subprocess.run(['cp', maf_files_in_dir[0], lastz_output_file], check=True)
                print(f"[INFO] LastZ alignment saved to: {lastz_output_file}")
            else:
                print(f"[WARNING] No MAF files found in output directory: {temp_output_dir}")
                # Keep looking for any maf file with a more flexible approach
                # This is a fallback in case the file naming is different than expected
                lastz_output_files[-1] = None  # Mark as not found in case we don't find anything
            
            # Clean up temp directory after processing (optional)
            import shutil
            shutil.rmtree(temp_output_dir, ignore_errors=True)
            
        # Filter out any None values from the output files list
        lastz_output_files = [f for f in lastz_output_files if f is not None]
        maf_files = lastz_output_files  # Store the generated MAF files

    ##############################################################################
    #                          3. MultiZ (If More Than Two Genomes)
    ##############################################################################

    if use_multiz and len(maf_files) > 1:  # Only run MultiZ if we have multiple MAF files
        import tempfile
        
        # Create a temporary file that will be automatically deleted when closed
        with tempfile.NamedTemporaryFile(suffix='.maf', delete=False) as temp_file:
            multiz_output_file = temp_file.name
            
        command = f'"{MULTIZ_PATH}" ' + " ".join(f'"{file}"' for file in maf_files) + f' > "{multiz_output_file}"'
        print("\n[INFO] Running MultiZ for multiple sequence alignment... ")
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"\n[ERROR] MultiZ failed with error:\n{e}")
            sys.exit(1)
        
        print("[INFO] MultiZ alignment completed successfully")
        maf_files = [multiz_output_file]  # Replace MAF files with MultiZ output
        
        # Register the temporary file for deletion when the program exits
        import atexit
        atexit.register(lambda file=multiz_output_file: os.remove(file) if os.path.exists(file) else None)

    ##############################################################################
    #                          4. Filtering
    ##############################################################################

    # Load the reference and query genomes
    print("[INFO] Loading genomes...")
    ref_genome = load_fasta(ref_path)
    qry_genome = load_fasta(qry_path)

    # Validate extracted sequences in parallel batches
    print("[INFO] Validating sequences...")
    validated_segments = []

    if not maf_files:
        print("[ERROR] No MAF files provided for sequence extraction.")
        sys.exit(1)

    segments = []

    for maf_file in maf_files:
        # Process each MAF file to extract segments
        maf_segments = process_maf_file_in_chunks(
            maf_file, 
            lambda alns: collect_segments_from_alignments(alns, genes), 
            MAF_CHUNK_SIZE
        )
        segments.extend(maf_segments)

    if not segments:
        print("[ERROR] No valid segments extracted from MAF.")
        sys.exit(1)

    print(f"[INFO] Extracted {len(segments)} initial segments from MAF files.")

    # Process in batches for better memory management
    batch_size = BATCH_SIZE
    seq_batches = list(chunks(segments, batch_size))

    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_PROCESSES) as executor:
        # Create a partial function with fixed genome parameters
        validate_batch = functools.partial(
            validate_sequences_batch, 
            ref_genome=ref_genome, 
            qry_genome=qry_genome, 
            extra_genomes=extra_genomes
        )
        
        # Submit all batches for processing
        futures = [executor.submit(validate_batch, batch) for batch in seq_batches]
        
        # Collect results as they complete
        if SHOW_PROGRESS:
            for future in tqdm(concurrent.futures.as_completed(futures), 
                            total=len(futures), 
                            desc="Verifying sequences"):
                validated_segments.extend(future.result())
        else:
            for future in concurrent.futures.as_completed(futures):
                validated_segments.extend(future.result())

    # Replace original segments list with the validated one
    segments = validated_segments

    print(f"[INFO] Retained {len(segments)} sequences after exact-match verification.")
    if not segments:
        print("[ERROR] No sequences remained after verification. Exiting.")
        sys.exit(1)

    ##############################################################################
    #                          5. Run Primer3
    ##############################################################################

    primer3_settings_content = PRIMER3_SETTINGS.strip()

    # Create primer3 input blocks
    primer3_input_blocks = []
    for seg in segments:
        seq_id = f"{seg['gene']}_{seg['index']}"
        block = (
            f"SEQUENCE_ID={seq_id}\n"
            f"SEQUENCE_TEMPLATE={seg['sequence']}\n"
            f"{primer3_settings_content}\n"
            "=\n"
        )
        primer3_input_blocks.append(block)

    # Run primer3 in batches for better memory management and parallel processing
    print("\n[INFO] Running Primer3...")
    batch_size = max(1, len(primer3_input_blocks) // (NUM_PROCESSES * 2))  # Smaller batches for memory management
    batches = list(chunks(primer3_input_blocks, batch_size))
    
    all_records = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_PROCESSES) as executor:
        # Submit all batches for processing
        batch_futures = []
        for i, batch in enumerate(batches):
            future = executor.submit(run_primer3_batch, batch)
            batch_futures.append(future)
        
        # Process outputs as they complete
        if SHOW_PROGRESS:
            for future in tqdm(concurrent.futures.as_completed(batch_futures), 
                             total=len(batch_futures), 
                             desc="Processing Primer3 batches"):
                stdout_data = future.result()
                if stdout_data:
                    batch_records = parse_primer3_batch(stdout_data)
                    all_records.extend(batch_records)
        else:
            for future in concurrent.futures.as_completed(batch_futures):
                stdout_data = future.result()
                if stdout_data:
                    batch_records = parse_primer3_batch(stdout_data)
                    all_records.extend(batch_records)

    # Create dataframe from records
    df = pd.DataFrame(all_records)
    if df.empty:
        print("[ERROR] No suitable primer pairs found (no pairs with penalty ≤ 5). Exiting.")
        sys.exit(0)

    ##############################################################################
    #                           6. Filter Primers - Repeats and GC%
    ##############################################################################

    print("[INFO] Filtering primers based on repeats and GC content...")
    # Apply filtering on the main thread (DataFrame operations aren't always faster in parallel)
    df["Amplicon GC%"] = df["Amplicon"].apply(calculate_gc)

    df = df[
        ~(
            df["Primer F"].apply(has_disallowed_repeats) |
            df["Primer R"].apply(has_disallowed_repeats) |
            (df["Amplicon GC%"] < 50) |
            (df["Amplicon GC%"] > 60)
        )
    ].reset_index(drop=True)

    if df.empty:
        print("[ERROR] After repeat/GC% filtering, no primers remain.")
        sys.exit(0)
    print(f"[INFO] Retained {len(df)} primer pairs after filtering.")

    ##############################################################################
    #                           7. Run BLASTn
    ##############################################################################

    db_name = f'"{DB_PATH}"'.split("/")[-1].split("_db")[0]
    print(f"\n[INFO] Running BLAST against {db_name}...")
    
    # Run BLAST in parallel for each column
    blast_cols = {"Primer F": "Primer F", "Primer R": "Primer R", "Probe": "Probe"}
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_PROCESSES) as executor:
        for col in blast_cols.keys():
            # Convert Series to list and split into batches
            seq_list = df[col].tolist()
            batch_size = BATCH_SIZE  # Adjust this value based on your sequences
            seq_batches = list(chunks(seq_list, batch_size))
            
            # Prepare batch data with column name
            batch_data = [(batch, col) for batch in seq_batches]
            
            # Submit all batches for processing
            futures = [executor.submit(process_blast_batch, data) for data in batch_data]
            
            # Collect results as they complete
            all_results = []
            if SHOW_PROGRESS:
                for future in tqdm(concurrent.futures.as_completed(futures), 
                                 total=len(futures), 
                                 desc=f"Blasting {col}"):
                    all_results.extend(future.result())
            else:
                for future in concurrent.futures.as_completed(futures):
                    all_results.extend(future.result())
                    
            # Update DataFrame with results
            base_col = blast_cols[col]
            df[f"{base_col} BLAST1"], df[f"{base_col} BLAST2"] = zip(*all_results)
    
    # Filter by BLAST specificity
    print("[INFO] Filtering primers based on BLAST specificity...")
    keep_indices = []
    for i in range(len(df)):
        row = df.iloc[i]
        keepF = passes_blast_filter(row, "Primer F")
        keepR = passes_blast_filter(row, "Primer R")
        keepP = pd.isna(row["Probe"]) or row["Probe"] == "" or passes_blast_filter(row, "Probe")
        if keepF and keepR and keepP:
            keep_indices.append(i)

    df = df.loc[keep_indices].reset_index(drop=True)
    print(f"[INFO] Retained {len(df)} primer pairs after BLAST specificity filtering.")
    if df.empty:
        print("[ERROR] No primers passed the BLAST specificity filter.")
        sys.exit(0)
        
    ##############################################################################
    #                           8. Find location
    ##############################################################################        

    # Process reference genome locations
    print("[INFO] Finding amplicon locations in genomes...")
    amplicon_seqs = df["Amplicon"].tolist()
    batch_size = BATCH_SIZE
    seq_batches = list(chunks(amplicon_seqs, batch_size))
    
    ref_results = []
    qry_results = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_PROCESSES) as executor:
        # Process reference genome
        print("[INFO] Processing reference genome locations...")
        ref_futures = [executor.submit(find_amplicon_locations_batch, (batch, ref_genome)) 
                     for batch in seq_batches]
        
        for future in concurrent.futures.as_completed(ref_futures):
            ref_results.extend(future.result())
        
        # Process query genome
        print("[INFO] Processing query genome locations...")
        qry_futures = [executor.submit(find_amplicon_locations_batch, (batch, qry_genome)) 
                      for batch in seq_batches]
        
        for future in concurrent.futures.as_completed(qry_futures):
            qry_results.extend(future.result())
    
    # Update DataFrame
    df["Ref Chromosome"], df["Ref Start"] = zip(*ref_results)
    df["Qry Chromosome"], df["Qry Start"] = zip(*qry_results)

    ##############################################################################
    #                           9. Run NUPACK dG
    ##############################################################################

    print("\n[INFO] Calculating NUPACK ΔG values...")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_PROCESSES) as executor:
        for col in ["Primer F", "Primer R", "Probe", "Amplicon"]:
            # Skip if no sequences
            if df[col].isna().all():
                continue
                
            # Convert Series to list and split into batches
            seq_list = df[col].tolist()
            batch_size = BATCH_SIZE 
            seq_batches = list(chunks(seq_list, batch_size))
            
            # Submit all batches for processing
            futures = [executor.submit(calc_deltaG_batch, batch) for batch in seq_batches]
            
            # Collect results as they complete
            all_results = []
            if SHOW_PROGRESS:
                for future in tqdm(concurrent.futures.as_completed(futures), 
                                 total=len(futures), 
                                 desc=f"Calculating ΔG for {col}"):
                    all_results.extend(future.result())
            else:
                for future in concurrent.futures.as_completed(futures):
                    all_results.extend(future.result())
                    
            # Update DataFrame with results
            dG_col = f"{col} dG"
            df[dG_col] = all_results

    ##############################################################################
    #                           10. Save Results
    ##############################################################################

    final_cols = [
        "Sequence", "Index", "Template",
        "Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST1", "Primer F BLAST2",
        "Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST1", "Primer R BLAST2",
        "Pair Penalty", "Probe", "Probe Tm", "Probe Penalty", "Probe dG", "Probe BLAST1", 
        "Probe BLAST2", "Amplicon", "Length", "Amplicon GC%", "Amplicon dG",
        "Ref Chromosome", "Ref Start", "Qry Chromosome", "Qry Start"
    ]
    
    # Reorder columns
    df = df[final_cols]

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    output_file = os.path.join(script_dir, f"Primer_Results_{timestamp}.xlsx")
    
    df.to_excel(output_file, index=False)
    
    print(f"\n[SUCCESS] Analysis complete! Found {len(df)} primer pairs.")
    print(f"[SUCCESS] Results saved to: {output_file}")
    print("\n=============================================================")

if __name__ == "__main__":
    main()