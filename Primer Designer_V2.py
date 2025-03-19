#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 21:30:15 2025

@author: jakob
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
tqdm.pandas()
import nupack

##############################################################################
#                           Variables
##############################################################################
MIN_SEGMENT_LENGTH = 100
RETAIN_TYPES = "gene"  ## gff filtering, "gene", "mRNA", "CDS", "exon",...
FILTER_MEANINGFUL_NAMES = True # will only use named genes from gff
COUNT_AMBIGUOUS_AS_MISMATCH = False
GENE_OVERLAP_MARGIN = 25
restriction_site = "GGCC"
PENALTY_MAX = 5.0
MAX_PRIMER_PAIRS_PER_SEGMENT = 3

# BLAST parameters
BLASTN_PATH = "/Library/Blast/ncbi-blast-2.16.0+/bin/blastn"
DB_PATH = "/Library/Blast/TAIR10_db"

BLAST_WORD_SIZE = 7
BLAST_EVALUE = 10
BLAST_MAX_TARGET_SEQS = 100
BLAST_REWARD = 2
BLAST_PENALTY = -3
BLAST_GAPOPEN = 5
BLAST_GAPEXTEND = 2

# BLAST specificity filter:
# best E-value must be at least BLAST_FILTER_FACTOR times smaller than the second-best E-value.
# If you literally want "50× smaller," set this to 50. If you only want best <= second, set it to 1.
BLAST_FILTER_FACTOR = 100

# NUPACK parameters
NUPACK_TEMPERATURE = 37  # Celsius
NUPACK_SODIUM = 0.05     # Molar
NUPACK_MAGNESIUM = 0.0   # Molar

##############################################################################
#                          File Loader
##############################################################################

def get_file(prompt, filetypes):
    """
    Show a file dialog to pick one file, or exit if none selected.
    """
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title=prompt, filetypes=filetypes)
    if not file_path:
        print(f"{prompt} - No file selected. Exiting.")
        sys.exit(1)
    return file_path

def load_fasta(filepath):
    """
    Load sequences from a FASTA file into a dict: {header_without_gt: sequence}.
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

def load_genes_from_gff(gff_path):
    """
    Return a list of gene dicts: {chr, start, end, strand, id}.
    Uses global RETAIN_TYPES and FILTER_MEANINGFUL_NAMES for filtering.
    """
    genes = []
    with open(gff_path, 'r') as f:
        for line in f:
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
                genes.append({
                    "chr": seqname,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "id": name
                })
            except ValueError:
                pass

    print(f"\nLoaded {len(genes)} entries from GFF.")
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
#                           MAF Alingment
##############################################################################

def parse_maf(maf_stream):
    alignments = []
    current_alignment = {}

    for line in maf_stream:
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
#                           Main Program
##############################################################################

def main():
    print("=== Combined Primer Design Pipeline ===")

##############################################################################
#                          1. File Import
##############################################################################

    print("Select the Reference FASTA File (reference genome)")
    ref_path = get_file("Select Reference FASTA File", [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])

    print("Select the Query FASTA File (query genome)")
    qry_path = get_file("Select Query FASTA File", [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])

    print("Select the GFF3 Annotation File")
    gff_path = get_file("Select GFF3 Annotation File", [("GFF3 Files", "*.gff3 *.gff"), ("All Files", "*.*")])

    print("Select the Primer3 Settings File")
    primer3_settings_path = get_file("Select Primer3 Settings File", [("Text Files", "*.txt"), ("All Files", "*.*")])

##############################################################################
#                          2. Run LastZ
##############################################################################
  
    genes = load_genes_from_gff(gff_path)

    #    Ensure `lastz` is in PATH, or replace `lastz_cmd = ["/path/to/lastz", ...]`
    lastz_cmd = [
        "lastz", f"{ref_path}[multiple]", f"{qry_path}[multiple]",
        "--format=maf", "--ambiguous=iupac"
    ]
    print("Running LastZ... ", end="", flush=True)
    stop_event = threading.Event()
    thread = threading.Thread(target=spinner_task, args=(stop_event,))
    thread.start()

    proc = subprocess.Popen(lastz_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    lastz_output, lastz_error = proc.communicate()

    stop_event.set()
    thread.join()

    if proc.returncode != 0:
        print(f"\nLastZ failed with error:\n{lastz_error}")
        sys.exit(1)

    maf_alignments = parse_maf(StringIO(lastz_output))
    print(f"Parsed {len(maf_alignments)} alignments from MAF.")
    
##############################################################################
#                          3. Filter - Gene Overlap,
#                          Mismatches and Restriction Site                      
##############################################################################  
    
    print("\nFiltering sequences...")

    segments = []
    stop_event = threading.Event()
    spinner_thread = threading.Thread(target=spinner_task, args=(stop_event,))
    spinner_thread.start()

    for aln in maf_alignments:
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
                        n_existing = sum(1 for s in segments if s['gene'] == gene['id'])
                        segments.append({
                            "gene": gene['id'],
                            "index": n_existing + 1,
                            "sequence": frag
                        })
    stop_event.set()
    spinner_thread.join()

##############################################################################
#                          4. Filtering - Doublecheck
##############################################################################

    def exact_match_in_genome(sequence, genome_dict):
        """
        Check if a given sequence is found exactly in any of the sequences in genome_dict.
        """
        return any(sequence in seq for seq in genome_dict.values())

    # Load the reference and query genomes
    ref_genome = load_fasta(ref_path)
    qry_genome = load_fasta(qry_path)

    # Validate extracted sequences before proceeding with a progress bar
    validated_segments = []
    for seg in tqdm(segments, desc="Verifying sequences"):
        seq = seg["sequence"]
        if exact_match_in_genome(seq, ref_genome) and exact_match_in_genome(seq, qry_genome):
            validated_segments.append(seg)

    # Replace original segments list with the validated one
    segments = validated_segments

    print(f"\nRetained {len(segments)} sequences after exact-match verification.")
    if not segments:
        print("No sequences remained after verification. Exiting.")
        sys.exit(1)

##############################################################################
#                          5. Run Primer3
##############################################################################

    with open(primer3_settings_path, 'r') as f:
        primer3_settings_content = f.read().strip()

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

    primer3_input_str = "".join(primer3_input_blocks)

    #    Make sure `primer3_core` is in PATH, or replace with the correct full path
    print("\nRunning Primer3...")
    
    def run_primer3(primer3_input_str):
        """
        Runs primer3_core using the Homebrew-installed location.
        """
        primer3_core_path = "/opt/homebrew/bin/primer3_core"  # Adjust if `which primer3_core` gives something else
        
        stop_event = threading.Event()
        spinner_thread = threading.Thread(target=spinner_task, args=(stop_event,))
        spinner_thread.start()

        primer3_proc = subprocess.Popen(
            [primer3_core_path],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True
        )
        stdout_data, stderr_data = primer3_proc.communicate(input=primer3_input_str)
        
        stop_event.set()
        spinner_thread.join()

        if primer3_proc.returncode != 0:
            print("\nPrimer3 Error:")
            print(stderr_data)
            sys.exit(1)

        if stderr_data:
            print("\n⚠️ Primer3 warnings or info:")
            print(stderr_data)

        print("Primer3 completed.", end=" ")
        return stdout_data

    # Run primer3 and get the output
    stdout_data = run_primer3(primer3_input_str)

##############################################################################
#                          6. Parse Primer 3 Output
##############################################################################

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

    df = pd.DataFrame(records)
    if df.empty:
        print("No suitable primer pairs found (no pairs with penalty ≤ 5). Exiting.")
        sys.exit(0)

##############################################################################
#                           7. Filter Primers - Repeats and GC%
##############################################################################

    def has_disallowed_repeats(seq):
        if not isinstance(seq, str):
            return True
        return "CCCC" in seq or "GGGG" in seq

    def calculate_gc(seq):
        if not seq or not isinstance(seq, str):
            return 0
        seq = seq.upper()
        gc_count = sum(1 for base in seq if base in "GC")
        return (gc_count / len(seq)) * 100 if seq else 0

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
        print("After repeat/GC% filtering, no primers remain.")
        sys.exit(0)
    print(f"Retained {len(df)} primer pairs after penatly filtering.")

##############################################################################
#                           8. Run BLASTn
##############################################################################

    db_name = DB_PATH.split("/")[-1].split("_db")[0]
    print(f"\nBlasting versus {db_name}")
    
    def blast_short_seq(seq, db=DB_PATH):
        """Runs BLASTn for short sequences and returns the two lowest e-values separately."""
        if not seq or not isinstance(seq, str) or not seq.strip():
            return None, None  # Ensuring consistency in output

        # Use a temporary file for the query sequence
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_query:
            tmp_query.write(f">seq\n{seq}\n")
            tmp_query.flush()
            tmp_filename = tmp_query.name

        # Run BLASTn
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

        os.remove(tmp_filename)  # Clean up temp file

        # If BLAST fails, print the error and return None values
        if result.returncode != 0:
            print(f"BLAST Error: {result.stderr}")
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

    # Update column names for BLAST results
    blast_cols = {"Primer F": "Primer F", "Primer R": "Primer R", "Probe": "Probe"}

    # Create new columns in the DataFrame
    for col in blast_cols.values():
        df[f"{col} BLAST1"] = None
        df[f"{col} BLAST2"] = None

    # Run BLAST and store results in separate columns
    for col, base_col in blast_cols.items():
        results = [blast_short_seq(seq) if pd.notnull(seq) else (None, None) for seq in tqdm(df[col], desc=f"Blasting {col}")]

        # Unpack tuples into two columns
        df[f"{base_col} BLAST1"], df[f"{base_col} BLAST2"] = zip(*results)

    # Filter by BLAST specificity
    def passes_blast_filter(row, col):
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

    keep_indices = []
    for i in range(len(df)):
        row = df.iloc[i]
        keepF = passes_blast_filter(row, "Primer F")
        keepR = passes_blast_filter(row, "Primer R")
        keepP = row["Probe"] == "" or passes_blast_filter(row, "Probe")  # Only check probes if they exist
        if keepF and keepR and keepP:
            keep_indices.append(i)

    df = df.loc[keep_indices].reset_index(drop=True)
    print(f"Retained {len(df)} primer pairs after BLAST specificity filtering.")
    if df.empty:
        print("No primers passed the BLAST specificity filter.")
        sys.exit(0)
        
##############################################################################
#                           9. Find location
##############################################################################        
        
    def find_amplicon_location(amplicon_seq, genome_dict):
        """
        Search for an exact match of amplicon_seq in genome_dict.
        Returns: (chromosome_name, start_position) or (None, None) if not found.
        """
        for chrom, seq in genome_dict.items():
            start_pos = seq.find(amplicon_seq)
            if start_pos != -1:
                return chrom, start_pos + 1  # Convert to 1-based index
        return None, None

    # Find locations of amplicons in both genomes
    df["Ref Chromosome"], df["Ref Start"] = zip(
        *df["Amplicon"].progress_apply(lambda seq: find_amplicon_location(seq, ref_genome) if pd.notnull(seq) else (None, None))
    )

    df["Qry Chromosome"], df["Qry Start"] = zip(
        *df["Amplicon"].progress_apply(lambda seq: find_amplicon_location(seq, qry_genome) if pd.notnull(seq) else (None, None))
    )

##############################################################################
#                           10. Run NUPACK dG
##############################################################################

    print("\nCalculating NUPACK ΔG for primers, probes, and amplicons...")
    model = nupack.Model(
        material='dna',
        celsius=NUPACK_TEMPERATURE,
        sodium=NUPACK_SODIUM,
        magnesium=NUPACK_MAGNESIUM
    )
    dna_pattern = re.compile(r'^[ACGTNacgtn]+$')

    def calc_deltaG(seq):
        if not isinstance(seq, str) or seq == "":
            return None
        if not dna_pattern.match(seq):
            return None
        try:
            result = nupack.mfe(seq, model=model)
            if result:
                # result[0].energy is the MFE structure's free energy
                return result[0].energy
        except Exception as e:
            print(f"NUPACK error for sequence {seq[:20]}...: {e}")
            return None
        return None

    for col in ["Primer F", "Primer R", "Probe", "Amplicon"]:
        dG_col = f"{col} dG"
        df[dG_col] = df[col].progress_apply(calc_deltaG)

##############################################################################
#                           11. Save Results
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
    output_file = f"Primer_Results_{timestamp}.xlsx"
    
    print()
    print(f"\nSaving results to {output_file}...")
    df.to_excel(output_file, index=False)
    
    print(f"Done! Analysis complete. Found {len(df)} primer pairs.")

if __name__ == "__main__":
    main()