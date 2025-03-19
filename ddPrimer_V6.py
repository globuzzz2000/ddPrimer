#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ddPrimer - Primer Design Pipeline for Multiple Genomes

This script designs primers that work across multiple genomes by:
1. Processing GFF files to identify regions of interest
2. Analyzing MAF alignment files to find conserved regions
3. Designing primers using Primer3
4. Validating primers with BLAST and NUPACK

Created on Fri Mar 7 21:30:15 2025
@author: jakob
"""

import os
import sys
import subprocess
import threading
import time
import re
import tempfile
import shutil
import tkinter as tk
from tkinter import filedialog
from io import StringIO
import pandas as pd
from tqdm import tqdm
from multiprocessing import cpu_count
import concurrent.futures
import functools
import nupack
from itertools import combinations

##############################################################################
#                           Configuration
##############################################################################
class Config:
    """Central configuration settings for the ddPrimer pipeline."""
    
    # Design Parameters
    MIN_SEGMENT_LENGTH = 100
    RETAIN_TYPES = "gene"  # gff filtering: "gene", "mRNA", "CDS", "exon", etc.
    FILTER_MEANINGFUL_NAMES = True  # only use named genes from gff
    COUNT_AMBIGUOUS_AS_MISMATCH = False
    GENE_OVERLAP_MARGIN = 25
    RESTRICTION_SITE = "GGCC"
    PENALTY_MAX = 5.0
    MAX_PRIMER_PAIRS_PER_SEGMENT = 3

    # Tool paths
    PLASTZ_PATH = "/Library/Application Support/ddPrimer/PLastZ/PLastZ.py"
    PRIMER3_CORE_PATH = "/opt/homebrew/bin/primer3_core"
    MULTIZ_PATH = "/Library/Application Support/ddPrimer/Multiz/multiz"
    BLASTN_PATH = "/Library/Application Support/ddPrimer/Blast+/bin/blastn"
    DB_PATH = "/Library/Application Support/ddPrimer/Tair DB/TAIR10_db"

    # Blast+ parameters
    BLAST_WORD_SIZE = 7
    BLAST_EVALUE = 10
    BLAST_MAX_TARGET_SEQS = 100
    BLAST_REWARD = 2
    BLAST_PENALTY = -3
    BLAST_GAPOPEN = 5
    BLAST_GAPEXTEND = 2
    BLAST_FILTER_FACTOR = 100  # E-value filtering factor

    # NUPACK parameters
    NUPACK_TEMPERATURE = 37  # Celsius
    NUPACK_SODIUM = 0.05     # Molar
    NUPACK_MAGNESIUM = 0.0   # Molar

    # Performance settings
    NUM_PROCESSES = max(1, int(cpu_count() * 0.75))  # Use 75% of cores
    BATCH_SIZE = 100
    MAF_CHUNK_SIZE = 10000
    SHOW_PROGRESS = True

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
        "PRIMER_OPT_TM": 57.5, "PRIMER_OUTSIDE_PENALTY": 0.0, "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
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

    @classmethod
    def format_primer3_settings(cls):
        """Format the Primer3 settings dictionary as a string."""
        return "\n".join(f"{key}={value}" for key, value in cls.PRIMER3_SETTINGS_DICTIONARY.items()) + "\n"

    @classmethod
    def get_primer3_settings(cls):
        """Get formatted Primer3 settings string."""
        return cls.format_primer3_settings()


##############################################################################
#                          File I/O Utilities
##############################################################################
class FileUtilities:
    """Handles file operations including UI dialogs and file parsing."""

    @staticmethod
    def get_file(prompt, filetypes):
        """Show a file dialog to pick one file, or exit if none selected."""
        os.environ['TK_SILENCE_DEPRECATION'] = '1'  # Suppress macOS warnings

        # Check for macOS-specific fix
        try:
            from Cocoa import NSApp
            macos_fix = True
        except ImportError:
            macos_fix = False

        root = tk.Tk()
        root.withdraw()  # Hide the Tk window

        # Apply macOS-specific fix
        if macos_fix:
            NSApp().setActivationPolicy_(1)  # Prevent Tk from appearing in the menu bar

        file_path = filedialog.askopenfilename(title=prompt, filetypes=filetypes)
        root.destroy()  # Ensure Tk properly closes after selection

        if not file_path:
            print(f"{prompt} - No file selected. Exiting.")
            sys.exit(1)

        return file_path

    @staticmethod
    def get_files(prompt, filetypes):
        """Show a file dialog to pick multiple files, or exit if none selected."""
        os.environ['TK_SILENCE_DEPRECATION'] = '1'  # Suppress macOS warnings

        # Check for macOS-specific fix
        try:
            from Cocoa import NSApp
            macos_fix = True
        except ImportError:
            macos_fix = False

        root = tk.Tk()
        root.withdraw()  # Hide the Tk window

        # Apply macOS-specific fix
        if macos_fix:
            NSApp().setActivationPolicy_(1)  # Prevent Tk from appearing in the menu bar

        file_paths = filedialog.askopenfilenames(title=prompt, filetypes=filetypes)
        root.destroy()  # Ensure Tk properly closes after selection

        if not file_paths:
            print(f"{prompt} - No files selected. Exiting.")
            sys.exit(1)

        return file_paths

    @staticmethod
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
#                          GFF Processing
##############################################################################
class GFFProcessor:
    """Processes GFF annotation files to extract gene information."""
    
    # Patterns for general placeholder detection across species
    PLACEHOLDER_PATTERNS = [
        re.compile(r'^AT\dG\d{5}$', re.IGNORECASE),   # Arabidopsis
        re.compile(r'^LOC_.*$', re.IGNORECASE),       # NCBI-style, rice, etc.
        re.compile(r'^GRMZM.*$', re.IGNORECASE),      # Maize
        re.compile(r'^ENS.*$', re.IGNORECASE)         # Ensembl IDs (humans, etc.)
    ]
    
    @staticmethod
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
    
    @classmethod
    def is_meaningful_name(cls, name, gene_id=None, locus_tag=None):
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

        for pattern in cls.PLACEHOLDER_PATTERNS:
            if pattern.match(name):
                return False

        return True
    
    @classmethod
    def process_gff_chunk(cls, chunk):
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
            if ftype not in Config.RETAIN_TYPES:
                continue

            seqname, _, _, start, end, _, strand, _, attributes = parts
            attr_dict = cls.parse_gff_attributes(attributes)

            name = attr_dict.get('name')
            gene_id = attr_dict.get('id')
            locus_tag = attr_dict.get('locus_tag')

            if Config.FILTER_MEANINGFUL_NAMES and not cls.is_meaningful_name(name, gene_id, locus_tag):
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
    
    @classmethod
    def load_genes_from_gff(cls, gff_path):
        """
        Return a list of gene dicts: {chr, start, end, strand, id}.
        Uses global RETAIN_TYPES and FILTER_MEANINGFUL_NAMES for filtering.
        Optimized version that processes the file in parallel chunks.
        """
        # Read all lines from the file
        with open(gff_path, 'r') as f:
            all_lines = f.readlines()
        
        # Calculate chunk size for parallel processing
        chunk_size = max(1, len(all_lines) // Config.NUM_PROCESSES)
        chunks = [all_lines[i:i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
        
        # Process chunks in parallel
        genes = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            futures = [executor.submit(cls.process_gff_chunk, chunk) for chunk in chunks]
            
            if Config.SHOW_PROGRESS:
                for future in tqdm(concurrent.futures.as_completed(futures), 
                                  total=len(futures), 
                                  desc="Processing GFF file"):
                    genes.extend(future.result())
            else:
                for future in concurrent.futures.as_completed(futures):
                    genes.extend(future.result())

        print(f"\nLoaded {len(genes)} entries from GFF.")
        return genes


##############################################################################
#                           UI Utilities
##############################################################################
class UIUtilities:
    """UI-related utility functions."""
    
    @staticmethod
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
class MAFProcessor:
    """Processes MAF alignment files to identify conserved segments."""
    
    @staticmethod
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
    
    @classmethod
    def process_maf_file_in_chunks(cls, maf_path, chunk_processor, chunk_size=Config.MAF_CHUNK_SIZE):
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
                    alignments = cls.parse_maf_chunk(process_chunk)
                    new_segments = chunk_processor(alignments)
                    segments.extend(new_segments)
        
        # Process any remaining lines
        if chunk_lines:
            alignments = cls.parse_maf_chunk(chunk_lines)
            new_segments = chunk_processor(alignments)
            segments.extend(new_segments)
        
        return segments
    
    @staticmethod
    def bases_match(ref_base, qry_base):
        """
        Checks if two bases match, with optional ambiguity handling.
        If COUNT_AMBIGUOUS_AS_MISMATCH is True, any ambiguous base is treated as mismatch.
        """
        if ref_base == qry_base:
            return True
        if Config.COUNT_AMBIGUOUS_AS_MISMATCH:
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
    
    @classmethod
    def split_on_mismatches(cls, ref_seq, qry_seq):
        """
        Split ref_seq into perfect-match segments based on mismatches and indels.
        Return only segments of length >= MIN_SEGMENT_LENGTH.
        """
        segments = []
        curr_segment = []

        for r, q in zip(ref_seq, qry_seq):
            if r == '-' or q == '-' or not cls.bases_match(r, q):
                if len(curr_segment) >= Config.MIN_SEGMENT_LENGTH:
                    segments.append("".join(curr_segment))
                curr_segment = []
            else:
                curr_segment.append(r)

        if len(curr_segment) >= Config.MIN_SEGMENT_LENGTH:
            segments.append("".join(curr_segment))

        return segments
    
    @classmethod
    def collect_segments_from_alignments(cls, alignments, genes):
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
                gene_start = gene['start'] - Config.GENE_OVERLAP_MARGIN
                gene_end = gene['end'] + Config.GENE_OVERLAP_MARGIN

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
                perfect_segments = cls.split_on_mismatches(ref_fragment, qry_fragment)

                # Further split perfect segments on restriction sites
                for segment in perfect_segments:
                    # Use re.split to handle overlapping restriction sites
                    subfragments = re.split(f'(?={Config.RESTRICTION_SITE})', segment)
                    for i, frag in enumerate(subfragments):
                        # Remove the restriction site prefix if present (except for first fragment)
                        if i > 0 and frag.startswith(Config.RESTRICTION_SITE):
                            frag = frag[len(Config.RESTRICTION_SITE):]
                        
                        if len(frag) >= Config.MIN_SEGMENT_LENGTH:
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
class SequenceValidator:
    """Validates sequences against multiple genomes."""
    
    @staticmethod
    def exact_match_in_genome(sequence, genome_dict):
        """
        Check if a given sequence is found exactly in any of the sequences in genome_dict.
        """
        return any(sequence in seq for seq in genome_dict.values())

    @classmethod
    def validate_sequence(cls, seq_info, ref_genome, qry_genome, extra_genomes=None):
        """
        Validate a single segment sequence against multiple genomes.
        """
        if extra_genomes is None:
            extra_genomes = []
            
        seq = seq_info["sequence"]
        return (
            cls.exact_match_in_genome(seq, ref_genome) and 
            cls.exact_match_in_genome(seq, qry_genome) and
            all(cls.exact_match_in_genome(seq, g) for g in extra_genomes)
        )

    @classmethod
    def validate_sequences_batch(cls, batch, ref_genome, qry_genome, extra_genomes=None):
        """
        Process a batch of sequences in parallel.
        """
        if extra_genomes is None:
            extra_genomes = []
            
        results = []
        for seg in batch:
            if cls.validate_sequence(seg, ref_genome, qry_genome, extra_genomes):
                results.append(seg)
        return results


##############################################################################
#                          BLAST Functions
##############################################################################
class BlastProcessor:
    """Handles BLAST operations for primer specificity checking."""
    
    @staticmethod
    def blast_short_seq(seq, db=f'"{Config.DB_PATH}"'):
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
                    Config.BLASTN_PATH,
                    "-task", "blastn-short",
                    "-db", db,
                    "-query", tmp_filename,
                    "-word_size", str(Config.BLAST_WORD_SIZE),
                    "-evalue", str(Config.BLAST_EVALUE),
                    "-reward", str(Config.BLAST_REWARD),
                    "-penalty", str(Config.BLAST_PENALTY),
                    "-gapopen", str(Config.BLAST_GAPOPEN),
                    "-gapextend", str(Config.BLAST_GAPEXTEND),
                    "-max_target_seqs", str(Config.BLAST_MAX_TARGET_SEQS),
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

    @classmethod
    def process_blast_batch(cls, batch_data):
        """
        Process a batch of sequences for BLAST in parallel.
        """
        batch, col_name = batch_data
        results = []
        
        for seq in batch:
            if pd.notnull(seq):
                blast1, blast2 = cls.blast_short_seq(seq)
                results.append((blast1, blast2))
            else:
                results.append((None, None))
        
        return results
    
    @staticmethod
    def passes_blast_filter(row, col):
        """
        Check if a sequence passes the BLAST specificity filter.
        """
        best = row[f"{col} BLAST1"]
        second = row[f"{col} BLAST2"]
        # If best is None, no hits => discard
        if pd.isna(best):
            return False
        # If no second best, it's effectively unique
        if pd.isna(second):
            return True
        # best must be at least FILTER_FACTOR times smaller => best * FACTOR <= second
        return best * Config.BLAST_FILTER_FACTOR <= second


##############################################################################
#                        Primer3 Processing
##############################################################################
class Primer3Processor:
    """Handles Primer3 operations for primer design."""
    
    @staticmethod
    def run_primer3_batch(input_blocks):
        """
        Run primer3_core on a batch of input blocks.
        Returns the stdout output.
        """
        primer3_input_str = "".join(input_blocks)

        try:
            primer3_proc = subprocess.Popen(
                [Config.PRIMER3_CORE_PATH],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True
            )
            stdout_data, stderr_data = primer3_proc.communicate(input=primer3_input_str)
            
            if primer3_proc.returncode != 0:
                print("\nPrimer3 Error:")
                print(stderr_data)
                return ""
            
            if stderr_data and Config.SHOW_PROGRESS:
                print("\n⚠️ Primer3 warnings or info for batch:")
                print(stderr_data)
                
            return stdout_data
        except Exception as e:
            print(f"Error running Primer3: {e}")
            return ""

    @staticmethod
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
                if pair_penalty <= Config.PENALTY_MAX and probe_penalty <= Config.PENALTY_MAX:
                    acceptable.append(p)

            # take up to 3
            acceptable = acceptable[:Config.MAX_PRIMER_PAIRS_PER_SEGMENT]
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
class AmpliconLocator:
    """Locates amplicons in genomes."""
    
    @staticmethod
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
class NupackProcessor:
    """Handles thermodynamic calculations using NUPACK."""
    
    @staticmethod
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
                celsius=Config.NUPACK_TEMPERATURE,
                sodium=Config.NUPACK_SODIUM,
                magnesium=Config.NUPACK_MAGNESIUM
            )
            result = nupack.mfe(seq, model=model)
            if result:
                # result[0].energy is the MFE structure's free energy
                return result[0].energy
        except Exception as e:
            if Config.SHOW_PROGRESS:
                preview = seq[:50] + ("..." if len(seq) > 50 else "")
                print(f"NUPACK error for sequence {preview}: {e}")
            return None
        return None

    @classmethod
    def calc_deltaG_batch(cls, seqs):
        """
        Calculate ΔG for a batch of sequences in parallel.
        """
        results = []
        for seq in seqs:
            if pd.notnull(seq):
                results.append(cls.calc_deltaG(seq))
            else:
                results.append(None)
        return results


##############################################################################
#                          Utility Functions
##############################################################################
class Utils:
    """General utility functions."""
    
    @staticmethod
    def chunks(lst, n):
        """Split a list into chunks of size n."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    @staticmethod
    def flatten_list(list_of_lists):
        """Flatten a list of lists into a single list."""
        return [item for sublist in list_of_lists for item in sublist]

    @staticmethod
    def has_disallowed_repeats(seq):
        """Check for disallowed repeats in a sequence."""
        if not isinstance(seq, str):
            return True
        return "CCCC" in seq or "GGGG" in seq

    @staticmethod
    def calculate_gc(seq):
        """Calculate GC content of a sequence."""
        if not seq or not isinstance(seq, str):
            return 0
        seq = seq.upper()
        gc_count = sum(1 for base in seq if base in "GC")
        return (gc_count / len(seq)) * 100 if seq else 0


##############################################################################
#                           Pipeline
##############################################################################
class PrimerDesignPipeline:
    """Main class that orchestrates the primer design pipeline."""
    
    def __init__(self):
        """Initialize the pipeline."""
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.ref_path = None
        self.qry_path = None
        self.extra_genome_paths = []
        self.extra_genomes = []
        self.gff_path = None
        self.maf_files = []
        self.use_multiz = False
        self.num_genomes = 2
        self.genes = []
        self.segments = []
        self.ref_genome = {}
        self.qry_genome = {}
        
    def get_user_input(self):
        """Get initial user inputs for the pipeline."""
        print("=== Combined Primer Design Pipeline ===")
        print("Do you want to compare two or more genomes?")
        print("1. Two genomes (skip MultiZ)")
        print("2. More than two genomes (use MultiZ)")
        comparison_choice = input("Enter 1 or 2: ").strip()
        
        self.use_multiz = comparison_choice == "2"
        self.num_genomes = 2 if not self.use_multiz else int(input("How many genomes? (≥3): ").strip())

        print("\nDo you want to run LastZ, or do you already have .maf files?")
        print("1. Run LastZ")
        print("2. Use existing .maf files")
        maf_choice = input("Enter 1 or 2: ").strip()

        if maf_choice == "1":
            print("Proceeding with LastZ alignments...")
        else:
            print("Select MAF files for analysis")
            self.maf_files = FileUtilities.get_files("Select MAF files", [("MAF Files", "*.maf"), ("All Files", "*.*")])
            
        return maf_choice
    
    def load_input_files(self):
        """Load all required input files."""
        print("Select the Reference FASTA File (reference genome)")
        self.ref_path = FileUtilities.get_file("Select Reference FASTA File", 
                                             [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])

        print("Select the Query FASTA File (query genome)")
        self.qry_path = FileUtilities.get_file("Select Query FASTA File", 
                                             [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])

        # Load additional genomes if needed
        num_extra = max(0, self.num_genomes - 2)
        if num_extra > 0:
            print(f"Select {num_extra} additional genomes for filtering:")
            for i in range(num_extra):
                extra_path = FileUtilities.get_file(f"Select {i+1}. Additional FASTA File", 
                                                  [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])
                self.extra_genome_paths.append(extra_path)
                self.extra_genomes.append(FileUtilities.load_fasta(extra_path))

        print("Select the GFF3 Annotation File")
        self.gff_path = FileUtilities.get_file("Select GFF3 Annotation File", 
                                             [("GFF3 Files", "*.gff3 *.gff"), ("All Files", "*.*")])
        
        # Process GFF file
        self.genes = GFFProcessor.load_genes_from_gff(self.gff_path)
        
        # Load reference and query genomes
        self.ref_genome = FileUtilities.load_fasta(self.ref_path)
        self.qry_genome = FileUtilities.load_fasta(self.qry_path)
    
    def run_lastz(self, maf_choice):
        """Run PLastZ alignments if needed."""
        if maf_choice == "1":
            lastz_output_files = []

            # Use already selected genomes: reference, query, and any additional genomes
            genome_paths = [self.ref_path, self.qry_path] + self.extra_genome_paths
            primary_reference = genome_paths[0]

            # Run PLastZ pairwise for all genome combinations
            for ref_genome_path, qry_genome_path in combinations(genome_paths, 2):
                # Get the directory of the input files
                input_dir = os.path.dirname(ref_genome_path)
                
                # Create an "Alignments" folder in the input directory if it doesn't exist
                alignments_dir = os.path.join(input_dir, "Alignments")
                os.makedirs(alignments_dir, exist_ok=True)
                
                # Simplified naming convention without timestamp or prefix
                alignment_name = (f"{os.path.splitext(os.path.basename(ref_genome_path))[0]}_"
                                f"vs_{os.path.splitext(os.path.basename(qry_genome_path))[0]}.maf")
                lastz_output_file = os.path.join(alignments_dir, alignment_name)
                lastz_output_files.append(lastz_output_file)

                # Modified to use PLastZ instead of LastZ
                # Using -p for processes and -lo for lastz-options (short forms)
                plastz_cmd = [
                    "python", Config.PLASTZ_PATH, 
                    ref_genome_path, qry_genome_path, 
                    self.script_dir,
                    "-lo=--format=maf --ambiguous=iupac",
                    "-p", str(Config.NUM_PROCESSES)
                ]
                
                print(f"\nRunning PLastZ for {os.path.basename(ref_genome_path)} vs {os.path.basename(qry_genome_path)}... ")
                stop_event = threading.Event()
                thread = threading.Thread(target=UIUtilities.spinner_task, args=(stop_event,))
                thread.start()
                
                # PLastZ outputs to a fixed file name "plastz.maf" in the output directory
                temp_output_name = os.path.join(self.script_dir, "plastz.maf")
                
                # Run PLastZ with stdout suppressed
                proc = subprocess.Popen(
                    plastz_cmd, 
                    stderr=subprocess.PIPE, 
                    stdout=subprocess.DEVNULL,
                    text=True
                )
                _, plastz_error = proc.communicate()

                stop_event.set()
                thread.join()

                if proc.returncode != 0:
                    print(f"\nPLastZ failed with error:\n{plastz_error}")
                    sys.exit(1)
                
                # Check if the output file exists and rename it
                if os.path.exists(temp_output_name):
                    try:
                        # Make sure the destination doesn't exist before copying
                        if os.path.exists(lastz_output_file):
                            os.remove(lastz_output_file)
                        
                        # Copy the file to our desired output location and name
                        shutil.copy(temp_output_name, lastz_output_file)
                        
                        # Optionally, remove the original file if you don't need it
                        os.remove(temp_output_name)
                        
                        print(f"\nSaved PLastZ output to {lastz_output_file}")
                    except Exception as e:
                        print(f"\nError renaming PLastZ output file: {e}")
                        sys.exit(1)
                else:
                    print(f"\nPLastZ did not generate the expected output file: {temp_output_name}")
                    print(f"Contents of directory: {os.listdir(self.script_dir)}")
                    sys.exit(1)

            self.maf_files = lastz_output_files
    
    def run_multiz(self):
        """Run MultiZ if needed for multiple genomes."""
        if self.use_multiz and len(self.maf_files) > 1:
            # Get the directory of the input files
            input_dir = os.path.dirname(self.ref_path)
            
            # Create an "Alignments" folder in the input directory if it doesn't exist
            alignments_dir = os.path.join(input_dir, "Alignments")
            os.makedirs(alignments_dir, exist_ok=True)
            
            # Generate filename based on genome names
            ref_name = os.path.splitext(os.path.basename(self.ref_path))[0]
            qry_name = os.path.splitext(os.path.basename(self.qry_path))[0]
            
            filename_parts = [ref_name, qry_name]
            
            # Add additional genomes to the filename if present
            for extra_path in self.extra_genome_paths:
                extra_name = os.path.splitext(os.path.basename(extra_path))[0]
                filename_parts.append(extra_name)
            
            # Join all genome names with + signs for the multiz output
            multiz_filename = "_+_".join(filename_parts) + "_multiz.maf"
            
            # Full path to the output file
            multiz_output_file = os.path.join(alignments_dir, multiz_filename)
            
            command = f'"{Config.MULTIZ_PATH}" ' + " ".join(f'"{file}"' for file in self.maf_files) + f' > "{multiz_output_file}"'
            print("\nRunning MultiZ for multiple sequence alignment... ")
            try:
                subprocess.run(command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"\nMultiZ failed with error:\n{e}")
                sys.exit(1)
            
            print(f"MultiZ alignment saved to: {multiz_output_file}")
            self.maf_files = [multiz_output_file]  # Replace MAF files with MultiZ output
    
    def extract_segments(self):
        """Extract conserved segments from MAF files."""
        if not self.maf_files:
            print("Error: No MAF files provided for sequence extraction.")
            sys.exit(1)

        for maf_file in self.maf_files:
            # Process each MAF file to extract segments
            maf_segments = MAFProcessor.process_maf_file_in_chunks(
                maf_file, 
                lambda alns: MAFProcessor.collect_segments_from_alignments(alns, self.genes), 
                Config.MAF_CHUNK_SIZE
            )
            self.segments.extend(maf_segments)

        if not self.segments:
            print("Error: No valid segments extracted from MAF.")
            sys.exit(1)

        print(f"Extracted {len(self.segments)} initial segments from MAF files.")
    
    def validate_segments(self):
        """Validate segments against all genomes."""
        validated_segments = []
        
        # Process in batches for better memory management
        batch_size = Config.BATCH_SIZE
        seq_batches = list(Utils.chunks(self.segments, batch_size))

        with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            # Create a partial function with fixed genome parameters
            validate_batch = functools.partial(
                SequenceValidator.validate_sequences_batch, 
                ref_genome=self.ref_genome, 
                qry_genome=self.qry_genome, 
                extra_genomes=self.extra_genomes
            )
            
            # Submit all batches for processing
            futures = [executor.submit(validate_batch, batch) for batch in seq_batches]
            
            # Collect results as they complete
            if Config.SHOW_PROGRESS:
                for future in tqdm(concurrent.futures.as_completed(futures), 
                                total=len(futures), 
                                desc="Verifying sequences"):
                    validated_segments.extend(future.result())
            else:
                for future in concurrent.futures.as_completed(futures):
                    validated_segments.extend(future.result())

        # Replace original segments list with the validated one
        self.segments = validated_segments

        print(f"\nRetained {len(self.segments)} sequences after exact-match verification.")
        if not self.segments:
            print("No sequences remained after verification. Exiting.")
            sys.exit(1)
    
    def run_primer3(self):
        """Run Primer3 for primer design."""
        primer3_settings_content = Config.get_primer3_settings().strip()

        # Create primer3 input blocks
        primer3_input_blocks = []
        for seg in self.segments:
            seq_id = f"{seg['gene']}_{seg['index']}"
            block = (
                f"SEQUENCE_ID={seq_id}\n"
                f"SEQUENCE_TEMPLATE={seg['sequence']}\n"
                f"{primer3_settings_content}\n"
                "=\n"
            )
            primer3_input_blocks.append(block)

        # Run primer3 in batches for better memory management and parallel processing
        print("\nRunning Primer3...")
        batch_size = max(1, len(primer3_input_blocks) // (Config.NUM_PROCESSES * 2))
        batches = list(Utils.chunks(primer3_input_blocks, batch_size))
        
        all_records = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            # Submit all batches for processing
            batch_futures = []
            for i, batch in enumerate(batches):
                future = executor.submit(Primer3Processor.run_primer3_batch, batch)
                batch_futures.append(future)
            
            # Process outputs as they complete
            if Config.SHOW_PROGRESS:
                for future in tqdm(concurrent.futures.as_completed(batch_futures), 
                                 total=len(batch_futures), 
                                 desc="Processing Primer3 batches"):
                    stdout_data = future.result()
                    if stdout_data:
                        batch_records = Primer3Processor.parse_primer3_batch(stdout_data)
                        all_records.extend(batch_records)
            else:
                for future in concurrent.futures.as_completed(batch_futures):
                    stdout_data = future.result()
                    if stdout_data:
                        batch_records = Primer3Processor.parse_primer3_batch(stdout_data)
                        all_records.extend(batch_records)

        # Create dataframe from records
        df = pd.DataFrame(all_records)
        if df.empty:
            print("No suitable primer pairs found (no pairs with penalty ≤ 5). Exiting.")
            sys.exit(0)
            
        return df
    
    def filter_primers(self, df):
        """Filter primers based on repeats and GC content."""
        df["Amplicon GC%"] = df["Amplicon"].apply(Utils.calculate_gc)

        df = df[
            ~(
                df["Primer F"].apply(Utils.has_disallowed_repeats) |
                df["Primer R"].apply(Utils.has_disallowed_repeats) |
                (df["Amplicon GC%"] < 50) |
                (df["Amplicon GC%"] > 60)
            )
        ].reset_index(drop=True)

        if df.empty:
            print("After repeat/GC% filtering, no primers remain.")
            sys.exit(0)
            
        print(f"Retained {len(df)} primer pairs after penalty filtering.")
        return df
    
    def run_blast(self, df):
        """Run BLAST for primer specificity checking."""
        db_name = f'"{Config.DB_PATH}"'.split("/")[-1].split("_db")[0]
        print(f"\nBlasting versus {db_name}...")
        
        # Run BLAST in parallel for each column
        blast_cols = {"Primer F": "Primer F", "Primer R": "Primer R", "Probe": "Probe"}
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            for col in blast_cols.keys():
                # Convert Series to list and split into batches
                seq_list = df[col].tolist()
                batch_size = Config.BATCH_SIZE
                seq_batches = list(Utils.chunks(seq_list, batch_size))
                
                # Prepare batch data with column name
                batch_data = [(batch, col) for batch in seq_batches]
                
                # Submit all batches for processing
                futures = [executor.submit(BlastProcessor.process_blast_batch, data) for data in batch_data]
                
                # Collect results as they complete
                all_results = []
                if Config.SHOW_PROGRESS:
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
        keep_indices = []
        for i in range(len(df)):
            row = df.iloc[i]
            keepF = BlastProcessor.passes_blast_filter(row, "Primer F")
            keepR = BlastProcessor.passes_blast_filter(row, "Primer R")
            keepP = pd.isna(row["Probe"]) or row["Probe"] == "" or BlastProcessor.passes_blast_filter(row, "Probe")
            if keepF and keepR and keepP:
                keep_indices.append(i)

        df = df.loc[keep_indices].reset_index(drop=True)
        print(f"Retained {len(df)} primer pairs after BLAST specificity filtering.")
        if df.empty:
            print("No primers passed the BLAST specificity filter.")
            sys.exit(0)
        
        return df
    
    def find_locations(self, df):
        """Find the genomic locations of amplicons."""
        # Process reference genome locations
        amplicon_seqs = df["Amplicon"].tolist()
        batch_size = Config.BATCH_SIZE
        seq_batches = list(Utils.chunks(amplicon_seqs, batch_size))
        
        ref_results = []
        qry_results = []
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            # Process reference genome
            print("Processing reference genome locations...")
            ref_futures = [executor.submit(AmpliconLocator.find_amplicon_locations_batch, (batch, self.ref_genome)) 
                         for batch in seq_batches]
            
            for future in concurrent.futures.as_completed(ref_futures):
                ref_results.extend(future.result())
            
            # Process query genome
            print("Processing query genome locations...")
            qry_futures = [executor.submit(AmpliconLocator.find_amplicon_locations_batch, (batch, self.qry_genome)) 
                          for batch in seq_batches]
            
            for future in concurrent.futures.as_completed(qry_futures):
                qry_results.extend(future.result())
        
        # Update DataFrame
        df["Ref Chromosome"], df["Ref Start"] = zip(*ref_results)
        df["Qry Chromosome"], df["Qry Start"] = zip(*qry_results)
        
        return df
    
    def run_nupack(self, df):
        """Calculate thermodynamic properties using NUPACK."""
        print("\nCalculating NUPACK ΔG...")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            for col in ["Primer F", "Primer R", "Probe", "Amplicon"]:
                # Skip if no sequences
                if df[col].isna().all():
                    continue
                    
                # Convert Series to list and split into batches
                seq_list = df[col].tolist()
                batch_size = Config.BATCH_SIZE 
                seq_batches = list(Utils.chunks(seq_list, batch_size))
                
                # Submit all batches for processing
                futures = [executor.submit(NupackProcessor.calc_deltaG_batch, batch) for batch in seq_batches]
                
                # Collect results as they complete
                all_results = []
                if Config.SHOW_PROGRESS:
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
        
        return df
    
    def save_results(self, df):
        """Save final results to Excel."""
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

        # Create a Primers directory in the script folder if it doesn't exist
        primers_dir = os.path.join(self.script_dir, "Primers")
        os.makedirs(primers_dir, exist_ok=True)
        
        # Generate filename based on genome names
        ref_name = os.path.splitext(os.path.basename(self.ref_path))[0]
        qry_name = os.path.splitext(os.path.basename(self.qry_path))[0]
        
        filename_parts = [ref_name, qry_name]
        
        # Add additional genomes to the filename if present
        for extra_path in self.extra_genome_paths:
            extra_name = os.path.splitext(os.path.basename(extra_path))[0]
            filename_parts.append(extra_name)
        
        # Join all genome names with + signs
        filename_base = "_+_".join(filename_parts)
        
        output_file = os.path.join(primers_dir, f"Primers_{filename_base}.xlsx")
        
        df.to_excel(output_file, index=False)
        
        print(f"Done! Analysis complete. Found {len(df)} primer pairs.")
        print(f"Results saved to: {output_file}")
    
    def run(self):
        """Run the full pipeline."""
        # 1. Get user input
        maf_choice = self.get_user_input()
        
        # 2. Load input files
        self.load_input_files()
        
        # 3. Run LastZ alignments if needed
        self.run_lastz(maf_choice)
        
        # 4. Run MultiZ if needed
        self.run_multiz()
        
        # 5. Extract segments from MAF files
        self.extract_segments()
        
        # 6. Validate segments
        self.validate_segments()
        
        # 7. Run Primer3
        df = self.run_primer3()
        
        # 8. Filter primers
        df = self.filter_primers(df)
        
        # 9. Run BLAST
        df = self.run_blast(df)
        
        # 10. Find locations
        df = self.find_locations(df)
        
        # 11. Run NUPACK
        df = self.run_nupack(df)
        
        # 12. Save results
        self.save_results(df)


##############################################################################
#                           Main Program
##############################################################################
def main():
    """Main entry point for the program."""
    pipeline = PrimerDesignPipeline()
    pipeline.run()


if __name__ == "__main__":
    main()