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
import regex  # Fast approximate matching
from Levenshtein import distance as levenshtein_distance
import tempfile
import shutil
import platform
import tkinter as tk
from tkinter import filedialog
from tqdm import tqdm
from multiprocessing import cpu_count
from itertools import combinations
import concurrent.futures
import functools
import pandas as pd
import nupack
import argparse

##############################################################################
#                           Configuration
##############################################################################
class Config:
    """Central configuration settings for the ddPrimer pipeline."""
    # Pipeline Mode Options
    ONLY_LASTZ_MULTIZ = False          # Run only LastZ/MultiZ alignments
    RUN_MATCH_CHECKER = False          # Run only the Match Checker
    RUN_MATCH_CHECKER_ON_OUTPUT = False # Run Match Checker on pipeline output
    DEBUG_MODE = False                 # Debug logging mode (enable with --debug flag))
    
    # Performance settings
    NUM_PROCESSES = max(1, int(cpu_count() * 0.75))  # Use 75% of cores
    BATCH_SIZE = 100
    MAF_CHUNK_SIZE = 10000
    SHOW_PROGRESS = True

    # Design Parameters
    MIN_SEGMENT_LENGTH = 100
    RETAIN_TYPES = "gene"  # gff filtering: "gene", "mRNA", "CDS", "exon", etc.
    FILTER_MEANINGFUL_NAMES = True  # only use named genes from gff
    COUNT_AMBIGUOUS_AS_MISMATCH = False
    GENE_OVERLAP_MARGIN = 25
    RESTRICTION_SITE = "GGCC"
    PENALTY_MAX = 5.0
    MAX_PRIMER_PAIRS_PER_SEGMENT = 3
    PREFER_PROBE_MORE_C_THAN_G = True  # Set to False to disable
    
    # Validation options
    VALIDATION_MODE = "TOLERANT"  # "STRICT" or "TOLERANT"
    ALLOW_AMP_MISMATCHES = 2  # Number of mismatches allowed in amplicon for TOLERANT mode
    ALLOW_AMP_MISMATCH_PERCENT = 5  # Percentage of mismatches allowed (0 = use absolute count)
    MAX_SEARCH_LENGTH = 1000000  # Limit search in large chromosomes to improve performance

    # Tool paths
    MULTIZ_PATH = "/Library/Application Support/ddPrimer/Multiz/multiz"
    DB_PATH = "/Library/Application Support/ddPrimer/Tair DB/TAIR10_db"

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

    @staticmethod
    def debug(message):
        """
        Print debug messages if debug mode is enabled.
        
        Args:
            message (str): The debug message to print
        """
        if Config.DEBUG_MODE:
            print(f"[DEBUG] {message}")


##############################################################################
#                          File I/O Utilities
##############################################################################
class FileUtilities:
    """Handles file operations including UI dialogs and file parsing."""
    
    # Default to GUI mode unless explicitly detected as headless
    use_cli = False

    # Detect headless environments explicitly with refined logic
    if platform.system() == "Linux":
        if not os.environ.get("DISPLAY", ""):
            use_cli = True
            reason = "No DISPLAY environment variable on Linux."
        else:
            reason = "DISPLAY variable found on Linux."
    elif platform.system() == "Windows":
        if not sys.stdout.isatty():
            use_cli = True
            reason = "Non-interactive terminal on Windows."
        else:
            reason = "Interactive terminal on Windows."
    elif platform.system() == "Darwin":
        # Force GUI mode by default on macOS unless explicitly headless
        use_cli = False
        reason = "Defaulting to GUI mode on macOS."
    else:
        reason = "Unknown platform, defaulting to GUI mode."

    @staticmethod
    def get_file(prompt, filetypes):
        """
        Show a file dialog or prompt for a file path if CLI mode is enabled.
        Falls back to CLI input if GUI initialization fails.
        """
        if FileUtilities.use_cli:
            file_path = input(f"{prompt} (enter file path manually): ").strip()
            if not file_path:
                print(f"{prompt} - No file selected. Exiting.")
                sys.exit(1)
            return file_path

        try:
            root = tk.Tk()
        except Exception as e:
            print(f"GUI initialization failed ({e}), falling back to CLI mode.")
            FileUtilities.use_cli = True
            return FileUtilities.get_file(prompt, filetypes)

        root.withdraw()  # Hide the Tk window

        # Apply macOS-specific fix if available
        try:
            from Cocoa import NSApp
            NSApp().setActivationPolicy_(1)
        except ImportError:
            pass

        # Redirect stderr to suppress NSOpenPanel warnings
        original_stderr = sys.stderr
        sys.stderr = open(os.devnull, "w")

        file_path = filedialog.askopenfilename(title=prompt, filetypes=filetypes)

        # Restore stderr
        sys.stderr = original_stderr

        root.destroy()

        if not file_path:
            print(f"{prompt} - No file selected. Exiting.")
            sys.exit(1)

        return file_path

    @staticmethod
    def get_files(prompt, filetypes):
        """
        Show a file dialog for multiple files or prompt if CLI mode is enabled.
        Falls back to CLI input if GUI initialization fails.
        """
        if FileUtilities.use_cli:
            paths = input(f"{prompt} (enter file paths separated by spaces): ").strip()
            file_paths = paths.split() if paths else []
            if not file_paths:
                print(f"{prompt} - No files selected. Exiting.")
                sys.exit(1)
            return file_paths

        try:
            root = tk.Tk()
        except Exception as e:
            print(f"GUI initialization failed ({e}), falling back to CLI mode.")
            FileUtilities.use_cli = True
            return FileUtilities.get_files(prompt, filetypes)

        root.withdraw()

        try:
            from Cocoa import NSApp
            NSApp().setActivationPolicy_(1)
        except ImportError:
            pass

        original_stderr = sys.stderr
        sys.stderr = open(os.devnull, "w")

        file_paths = filedialog.askopenfilenames(title=prompt, filetypes=filetypes)

        sys.stderr = original_stderr
        root.destroy()

        if not file_paths:
            print(f"{prompt} - No files selected. Exiting.")
            sys.exit(1)

        return file_paths

    @staticmethod
    def load_fasta(filepath):
        """
        Load sequences from a FASTA file into a dict: {header_without_gt: sequence}.
        Optimized for memory efficiency.
        
        Args:
            filepath (str): Path to the FASTA file
            
        Returns:
            dict: Dictionary mapping sequence headers to sequences
        """
        sequences = {}
        name = None
        seq_chunks = []
        
        Config.debug(f"Loading FASTA file: {filepath}")
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if name:
                            sequences[name] = "".join(seq_chunks).upper()
                        name = line[1:].split()[0]
                        seq_chunks = []
                    else:
                        seq_chunks.append(line)
                if name:
                    sequences[name] = "".join(seq_chunks).upper()
        except Exception as e:
            print(f"Error reading FASTA file {os.path.abspath(filepath)}: {e}")
            sys.exit(1)
        
        Config.debug(f"Loaded {len(sequences)} sequences from FASTA file")
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
        
        Args:
            attribute_str (str): GFF attribute string
            
        Returns:
            dict: Dictionary of attribute key-value pairs
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
        
        Args:
            name (str): Gene name to check
            gene_id (str): Gene ID for comparison
            locus_tag (str): Locus tag for comparison
            
        Returns:
            bool: True if name is meaningful, False otherwise
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
        
        Args:
            chunk (list): List of GFF file lines
            
        Returns:
            list: List of gene dictionaries
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
        
        Args:
            gff_path (str): Path to the GFF file
            
        Returns:
            list: List of gene dictionaries
        """
        Config.debug(f"Loading genes from GFF file: {gff_path}")
        # Read all lines from the file
        with open(gff_path, 'r') as f:
            all_lines = f.readlines()
        
        # Calculate chunk size for parallel processing
        chunk_size = max(1, len(all_lines) // Config.NUM_PROCESSES)
        chunks = [all_lines[i:i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
        
        Config.debug(f"Processing GFF in {len(chunks)} chunks with {Config.NUM_PROCESSES} processes")
        
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

        print(f"Loaded {len(genes)} entries from GFF.")
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
        Ensures proper cleanup when stopping.
        
        Args:
            stop_event (threading.Event): Event to signal when to stop the spinner
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
        
        # Make sure we clear the spinner character
        sys.stdout.write("\b \b")  # backspace, space, backspace to fully erase
        sys.stdout.flush()


##############################################################################
#                           PLastZ Runner
##############################################################################
# This code is based on PLastZ.py, originally created by Antoine Ho (AntoineHo)
# and obtained from GitHub (https://github.com/AntoineHo/PLastZ).
# It has been modified and integrated into ddPrimer for direct integration
# and improved compatibility.
#
# Original code licensed under Apache License, Version 2.0
# Original copyright information preserved as per license requirements.
# Modified for integration with ddPrimer by Jakob (2025)

import os
import errno
import shutil
import subprocess
import shlex
from multiprocessing import Pool
from Bio import SeqIO

class LastZRunner:
    """
    Integrated parallel LastZ runner for sequence alignment.
    Handles splitting, parallel alignment, and result aggregation.
    """
    
    def __init__(self, config=None):
        """
        Initialize with configuration settings.
        
        Args:
            config: Configuration object (defaults to global Config)
        """
        self.config = config if config else Config
    
    def create_directory(self, path):
        """
        Create a directory if it doesn't exist.
        
        Args:
            path (str): Directory path to create
            
        Returns:
            str: Absolute path to the directory
            
        Raises:
            Exception: If directory cannot be created
        """
        try:
            os.makedirs(path)
            return os.path.abspath(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise Exception(f"ERROR: Cannot create directory: {path}")
            else:
                return os.path.abspath(path)
    
    def run_command(self, cmd):
        """
        Execute a shell command with error handling.
        """
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            proc.communicate()
            return proc.returncode
        except Exception as e:
            print(f"Error executing command: {cmd}\n{e}")
            return -1
    
    def run_commands_parallel(self, commands, processes=None):
        """
        Run multiple commands in parallel.
        
        Args:
            commands (list): List of commands to execute
            processes (int): Number of parallel processes (default: from config)
        """
        if processes is None:
            processes = self.config.NUM_PROCESSES
        
        with Pool(processes=processes) as pool:
            pool.map(self.run_command, commands)
    
    def create_alignment_jobs(self, query_path, target_path, temp_dir, lastz_options=None):
        """
        Create extraction and alignment jobs.
        
        Args:
            query_path (str): Path to query FASTA file
            target_path (str): Path to target FASTA file
            temp_dir (str): Directory to store temporary files
            lastz_options (str): Additional options for LastZ
            
        Returns:
            tuple: Lists of extraction and alignment job commands
        """
        extract_jobs = []
        align_jobs = []
        extracted_sequences = []
        pairs_done = []
        
        Config.debug(f"Creating alignment jobs for {query_path} vs {target_path}")
        
        # Extract unique sequence IDs from query and target
        query_contigs = []
        for qry in SeqIO.parse(query_path, "fasta"):
            query_contigs.append(qry.id)
            
        target_contigs = []
        for tgt in SeqIO.parse(target_path, "fasta"):
            target_contigs.append(tgt.id)
        
        # Create jobs for each sequence pair
        for query_id in query_contigs:
            if query_id not in extracted_sequences:
                # Extract query sequence
                out_fa_path = os.path.join(temp_dir, f"{query_id}.fa")
                query_quoted = shlex.quote(query_path)
                out_quoted = shlex.quote(out_fa_path)
                extract_jobs.append(f"samtools faidx {query_quoted} {query_id} > {out_quoted}")
                extracted_sequences.append(query_id)
            
            for target_id in target_contigs:
                # Skip if we've already processed this pair (in either direction)
                if [query_id, target_id] in pairs_done:
                    continue
                
                # Mark both directions as processed
                pairs_done.append([query_id, target_id])
                pairs_done.append([target_id, query_id])
                
                # Extract target sequence if not already done
                if target_id not in extracted_sequences:
                    out_fa_path = os.path.join(temp_dir, f"{target_id}.fa")
                    target_quoted = shlex.quote(target_path)
                    out_quoted = shlex.quote(out_fa_path)
                    extract_jobs.append(f"samtools faidx {target_quoted} {target_id} > {out_quoted}")
                    extracted_sequences.append(target_id)
                
                # Create alignment job
                tmp_output = os.path.join(temp_dir, f"{query_id}_V_{target_id}.TMP")
                query_fa = os.path.join(temp_dir, f"{query_id}.fa")
                target_fa = os.path.join(temp_dir, f"{target_id}.fa")
                
                # Quote paths for shell safety
                query_fa_quoted = shlex.quote(query_fa)
                target_fa_quoted = shlex.quote(target_fa)
                tmp_output_quoted = shlex.quote(tmp_output)
                
                # Build LastZ command
                if lastz_options:
                    cmd = f"lastz {query_fa_quoted} {target_fa_quoted} {lastz_options} > {tmp_output_quoted}"
                else:
                    cmd = f"lastz {query_fa_quoted} {target_fa_quoted} > {tmp_output_quoted}"
                
                align_jobs.append(cmd)
        
        Config.debug(f"Created {len(extract_jobs)} extraction jobs and {len(align_jobs)} alignment jobs")
        return extract_jobs, align_jobs
    
    def run_parallel_alignment(self, ref_path, qry_path, output_dir, lastz_options=None, processes=None, keep_temp=False):
        """
        Run LastZ alignment between reference and query genomes in parallel.
        
        Args:
            ref_path (str): Path to reference FASTA file
            qry_path (str): Path to query FASTA file
            output_dir (str): Directory to store results
            lastz_options (str): Additional options for LastZ
            processes (int): Number of parallel processes
            keep_temp (bool): Whether to keep temporary files
            
        Returns:
            str: Path to the output MAF file
        """
        # Set default processes if not provided
        if processes is None:
            processes = self.config.NUM_PROCESSES
        
        # Ensure paths are absolute
        ref_path = os.path.abspath(ref_path)
        qry_path = os.path.abspath(qry_path)
        
        # Create output directory
        output_dir = self.create_directory(output_dir)
        temp_dir = self.create_directory(os.path.join(output_dir, "temp"))
        
        # Track .fai files created by samtools
        fai_files = []
        if ref_path != qry_path:
            fai_files = [f"{ref_path}.fai", f"{qry_path}.fai"]
        else:
            fai_files = [f"{ref_path}.fai"]
        
        Config.debug("\nPreparing LastZ alignment...")
        extract_jobs, align_jobs = self.create_alignment_jobs(ref_path, qry_path, temp_dir, lastz_options)
        
        Config.debug(f"Extracting sequences from reference and query genomes...")
        self.run_commands_parallel(extract_jobs, processes)
        
        print("Running LastZ alignments...")
        self.run_commands_parallel(align_jobs, processes)
        
        # Combine all alignment results
        output_file = os.path.join(output_dir, "plastz.maf")
        Config.debug("Combining alignment results...")
        
        # Use Python's file I/O capabilities for a more elegant solution
        with open(output_file, 'w') as outfile:
            # Get all TMP files in the temp directory
            import glob
            temp_files = glob.glob(os.path.join(temp_dir, "*.TMP"))
            for temp_file in temp_files:
                with open(temp_file, 'r') as infile:
                    outfile.write(infile.read())
        
        # Clean up temporary files
        if not keep_temp:
            Config.debug("Cleaning up temporary files...")
            shutil.rmtree(temp_dir)
            
            # Clean up .fai files
            for fai_file in fai_files:
                if os.path.exists(fai_file):
                    try:
                        os.remove(fai_file)
                    except OSError as e:
                        print(f"Warning: Could not remove index file {fai_file}: {e}")
        
        return output_file


##############################################################################
#                           MAF Alignment Processing
##############################################################################
class MAFProcessor:
    """Processes MAF alignment files to identify conserved segments."""
    
    # Precomputed lookup table for IUPAC compatibility
    IUPAC_COMPATIBILITY = {}
    # Populate the lookup table for IUPAC codes
    IUPAC_BASE_MAP = {
        'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
        'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
        'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
    }
    # For each pair of bases, precompute compatibility
    for base1 in IUPAC_BASE_MAP:
        for base2 in IUPAC_BASE_MAP:
            IUPAC_COMPATIBILITY[(base1, base2)] = not IUPAC_BASE_MAP[base1].isdisjoint(IUPAC_BASE_MAP[base2])

    @staticmethod
    def parse_maf_chunk(chunk_lines):
        """
        Parse a chunk of MAF format lines into alignment records.
        Improved to handle multiple sequences per alignment block.
        
        Args:
            chunk_lines (list): List of MAF file lines
            
        Returns:
            list: List of alignment dictionaries
        """
        alignments = []
        current_alignment = {"sequences": []}

        for line in chunk_lines:
            line = line.strip()
            if line.startswith('a '):
                # Start of a new alignment block
                if current_alignment["sequences"]:
                    alignments.append(current_alignment)
                current_alignment = {"sequences": []}
            elif line.startswith('s '):
                fields = line.split()
                if len(fields) >= 7:  # Make sure we have enough fields
                    seq_entry = {
                        "name": fields[1],
                        "start": int(fields[2]) + 1,  # Convert to 1-based
                        "sequence": fields[6]
                    }
                    current_alignment["sequences"].append(seq_entry)

        # Don't forget the last alignment
        if current_alignment["sequences"]:
            alignments.append(current_alignment)

        return alignments
    
    @classmethod
    def process_maf_file_in_chunks(cls, maf_path, chunk_processor, chunk_size=None):
        """
        Process a MAF file in chunks to avoid loading it all into memory.
        The chunk_processor function is called for each chunk of alignments.
        Pre-filters candidate regions and builds a lookup table for fast sequence matching.
        
        Args:
            maf_path (str): Path to the MAF file
            chunk_processor (function): Function to process each chunk of alignments
            chunk_size (int): Size of chunks to process (default: from Config)
            
        Returns:
            list: Combined list of processed segments
        """
        if chunk_size is None:
            chunk_size = Config.MAF_CHUNK_SIZE
            
        start_time = time.time()
        chunk_lines = []
        segments = []
        
        Config.debug(f"Processing MAF file in chunks: {maf_path}")
        
        with open(maf_path, 'r') as maf_file:
            for line in maf_file:
                chunk_lines.append(line)
                if len(chunk_lines) >= chunk_size:
                    last_a_idx = None
                    for i in range(len(chunk_lines) - 1, -1, -1):
                        if chunk_lines[i].startswith('a '):
                            last_a_idx = i
                            break
                    if last_a_idx is not None and last_a_idx < len(chunk_lines) - 1:
                        process_chunk = chunk_lines[:last_a_idx]
                        chunk_lines = chunk_lines[last_a_idx:]
                    else:
                        process_chunk = chunk_lines
                        chunk_lines = []
                    alignments = cls.parse_maf_chunk(process_chunk)
                    new_segments = chunk_processor(alignments)
                    segments.extend(new_segments)
        
        # Process any remaining lines
        if chunk_lines:
            alignments = cls.parse_maf_chunk(chunk_lines)
            new_segments = chunk_processor(alignments)
            segments.extend(new_segments)
        
        # Build lookup table for sequence validation
        for seg in segments:
            SequenceValidator.add_to_maf_lookup(seg["sequence"], ("MAF", seg["index"]))
        
        elapsed = time.time() - start_time
        Config.debug(f"MAF processing completed: {len(segments)} segments in {elapsed:.2f} seconds")
        return segments
    
    @classmethod
    def bases_match(cls, ref_base, qry_base):
        """
        Checks if two bases match, with optional ambiguity handling.
        Uses precomputed lookup for efficiency.
        
        Args:
            ref_base (str): Reference base
            qry_base (str): Query base
            
        Returns:
            bool: True if bases match, False otherwise
        """
        if ref_base == qry_base:
            return True
        if Config.COUNT_AMBIGUOUS_AS_MISMATCH:
            return False

        # Use precomputed lookup table for IUPAC compatibility
        ref_base, qry_base = ref_base.upper(), qry_base.upper()
        key = (ref_base, qry_base)
        if key in cls.IUPAC_COMPATIBILITY:
            return cls.IUPAC_COMPATIBILITY[key]
        
        # Fallback for unexpected bases
        return False
    
    @classmethod
    def split_on_mismatches(cls, ref_seq, qry_seq):
        """
        Split ref_seq into perfect-match segments based on mismatches and indels.
        Return only segments of length >= MIN_SEGMENT_LENGTH.
        
        Args:
            ref_seq (str): Reference sequence
            qry_seq (str): Query sequence
            
        Returns:
            list: List of matching segments
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
        Handles multiple genomes in MAF files by comparing all sequence pairs.
        """
        chunk_segments = []
        
        for aln in alignments:
            # Skip alignment blocks with fewer than 2 sequences
            if len(aln["sequences"]) < 2:
                continue
                
            # Process each gene
            for gene in genes:
                # Find all sequences in this alignment that overlap with the gene
                matching_seqs = []
                for seq in aln["sequences"]:
                    chr_name = seq["name"]
                    
                    # Skip sequences that don't match the gene's chromosome
                    if gene['chr'] != chr_name:
                        continue
                        
                    seq_start = seq["start"]
                    seq_sequence = seq["sequence"]
                    
                    # Calculate sequence end position
                    seq_end = seq_start + len(seq_sequence.replace('-', '')) - 1
                    
                    # Overlap range for the gene
                    gene_start = gene['start'] - Config.GENE_OVERLAP_MARGIN
                    gene_end = gene['end'] + Config.GENE_OVERLAP_MARGIN
                    
                    # If this sequence doesn't overlap the gene region, skip
                    if seq_end < gene_start or seq_start > gene_end:
                        continue
                        
                    # Overlapping portion
                    overlap_start = max(seq_start, gene_start)
                    overlap_end = min(seq_end, gene_end)
                    
                    # Indices relative to the alignment's first base (0-based)
                    aln_offset_start = overlap_start - seq_start
                    aln_offset_end = aln_offset_start + (overlap_end - overlap_start) + 1
                    
                    # Make sure we don't go out of bounds
                    if aln_offset_start < 0:
                        aln_offset_start = 0
                    if aln_offset_end > len(seq_sequence):
                        aln_offset_end = len(seq_sequence)
                    
                    # This sequence overlaps with the gene, add it to matching sequences
                    matching_seqs.append({
                        "name": chr_name,
                        "start": seq_start,
                        "sequence": seq_sequence,
                        "aln_offset_start": aln_offset_start,
                        "aln_offset_end": aln_offset_end
                    })
                
                # If we found at least one matching sequence, process all pairs
                if len(matching_seqs) >= 1:
                    # For each matching sequence, compare it with all other sequences in the alignment
                    for i, ref_seq in enumerate(matching_seqs):
                        ref_fragment = ref_seq["sequence"][ref_seq["aln_offset_start"]:ref_seq["aln_offset_end"]]
                        
                        # Compare with all other sequences in the alignment
                        for qry_seq in aln["sequences"]:
                            # Skip self-comparison
                            if qry_seq["name"] == ref_seq["name"]:
                                continue
                                
                            # Extract corresponding fragment from the query sequence
                            # Make sure query sequence is long enough
                            if len(qry_seq["sequence"]) <= ref_seq["aln_offset_end"]:
                                continue
                                
                            qry_fragment = qry_seq["sequence"][ref_seq["aln_offset_start"]:ref_seq["aln_offset_end"]]
                            
                            # Get segments based on validation mode
                            if Config.VALIDATION_MODE == "STRICT":
                                segments_to_process = cls.split_on_mismatches(ref_fragment, qry_fragment)
                            else:
                                if len(ref_fragment.replace('-', '')) >= Config.MIN_SEGMENT_LENGTH:
                                    cleaned_fragment = ''.join([r for r, q in zip(ref_fragment, qry_fragment) 
                                                            if r != '-' and q != '-'])
                                    segments_to_process = [cleaned_fragment] if len(cleaned_fragment) >= Config.MIN_SEGMENT_LENGTH else []
                                else:
                                    segments_to_process = []
                            
                            # Process segments
                            for segment in segments_to_process:
                                # Split on restriction sites
                                subfragments = re.split(f'(?={Config.RESTRICTION_SITE})', segment)
                                for j, frag in enumerate(subfragments):
                                    # Remove the restriction site prefix if present
                                    if j > 0 and frag.startswith(Config.RESTRICTION_SITE):
                                        frag = frag[len(Config.RESTRICTION_SITE):]
                                    
                                    if len(frag) >= Config.MIN_SEGMENT_LENGTH:
                                        # Find where this segment appears in the original fragment
                                        segment_start_in_fragment = ref_fragment.find(frag)
                                        
                                        if segment_start_in_fragment != -1:
                                            # Calculate actual genomic positions
                                            real_offset = ref_seq["aln_offset_start"] + segment_start_in_fragment
                                            
                                            # Calculate positions in both genomes
                                            ref_bases_count = len(ref_seq["sequence"][:real_offset].replace('-', ''))
                                            ref_pos = ref_seq["start"] + ref_bases_count - 1
                                            
                                            qry_bases_count = len(qry_seq["sequence"][:real_offset].replace('-', ''))
                                            qry_pos = qry_seq["start"] + qry_bases_count - 1
                                            
                                            n_existing = sum(1 for s in chunk_segments 
                                                        if s['gene'] == gene['id'] and 
                                                            s['ref_chr'] == ref_seq["name"] and
                                                            s['qry_chr'] == qry_seq["name"])
                                            
                                            # Create an identifier that records both genomes
                                            pair_id = f"{ref_seq['name']}_{qry_seq['name']}"
                                            
                                            chunk_segments.append({
                                                "gene": gene['id'],
                                                "index": n_existing + 1,
                                                "sequence": frag,
                                                "ref_chr": ref_seq["name"],
                                                "ref_pos": ref_pos,
                                                "qry_chr": qry_seq["name"],
                                                "qry_pos": qry_pos,
                                                "genome_pair": pair_id  # Add this to track which genome pair this segment belongs to
                                            })
        
        return chunk_segments

##############################################################################
#                        Sequence Validation 
##############################################################################
class SequenceValidator:
    """Handles sequence validation against multiple genomes with optimized performance."""
    # Lookup table built from MAF alignment data for fast sequence matching
    maf_lookup = {}
    
    # Initialize sequence index cache for faster lookups
    sequence_index = {}
    
    # Cache of previously validated sequences for reuse
    validation_cache = {}
    
    @classmethod
    def add_to_maf_lookup(cls, sequence, location):
        """
        Add a sequence to the MAF lookup table.
        
        Args:
            sequence (str): Sequence to add
            location (tuple): Location tuple (type, index)
        """
        if sequence:
            cls.maf_lookup[sequence] = location
    
    @classmethod
    def build_sequence_index(cls, genome_dict, k=10):
        """
        Build a k-mer index for faster sequence lookups.
        Index structure: {k-mer: [(chrom, position), ...]}
        
        Args:
            genome_dict (dict): Dictionary mapping chromosome names to sequences
            k (int): k-mer size
        """
        start_time = time.time()
        genome_id = hash(tuple(sorted(genome_dict.keys())))
        
        total_size = sum(len(seq) for seq in genome_dict.values())
        Config.debug(f"Building k-mer index for genome (total size: {total_size} bases)")
        
        if genome_id in cls.sequence_index:
            Config.debug(f"Using existing k-mer index for genome with id {genome_id}")
            return  # Already indexed
            
        index = {}
        # Use larger stride for very large genomes to reduce index size
        stride = 1
        if total_size > 100_000_000:  # 100MB
            stride = 2
        elif total_size > 500_000_000:  # 500MB
            stride = 3
            
        for chrom, seq in genome_dict.items():
            for i in range(0, len(seq) - k + 1, stride):
                kmer = seq[i:i+k]
                index.setdefault(kmer, []).append((chrom, i))
        
        cls.sequence_index[genome_id] = {'k': k, 'index': index, 'genome': genome_dict}
        
        total_kmers = sum(len(v) for v in index.values())
        elapsed = time.time() - start_time
        Config.debug(f"k-mer index built: {total_kmers} k-mers in {elapsed:.2f} seconds")
    
    @classmethod
    def matches_with_tolerance(cls, sequence, genome_dict):
        """
        Fast search for matches in genome with allowed mismatches.
        Uses multiple optimized strategies based on sequence and genome properties.
        
        Args:
            sequence (str): Sequence to search for
            genome_dict (dict): Dictionary mapping chromosome names to sequences
            
        Returns:
            tuple: (chromosome, position) if match found, (None, None) otherwise
        """
        # Skip invalid sequences
        if not sequence or len(sequence) < Config.MIN_SEGMENT_LENGTH:
            return None, None
        if Config.VALIDATION_MODE == "STRICT":
            pass
        
        # Generate cache key and check cache first (fast path)
        genome_id = hash(tuple(sorted(genome_dict.keys())))
        cache_key = (sequence, genome_id)
        if cache_key in cls.validation_cache:
            return cls.validation_cache[cache_key]
        
        # Check MAF lookup (second fastest path)
        if sequence in cls.maf_lookup:
            Config.debug(f"Using MAF lookup for sequence: {sequence[:30]}...")
            return cls.maf_lookup[sequence]
        
        # Calculate allowed mismatches
        if Config.ALLOW_AMP_MISMATCH_PERCENT > 0:
            allowed_mismatches = max(1, int(len(sequence) * (Config.ALLOW_AMP_MISMATCH_PERCENT / 100.0)))
        else:
            allowed_mismatches = Config.ALLOW_AMP_MISMATCHES
        
        # Get genome size for strategy selection
        total_genome_size = sum(len(seq) for seq in genome_dict.values())
        
        # STRATEGY 1: Direct search for small genomes/short sequences
        if total_genome_size < 1_000_000 or len(sequence) < 30:
            Config.debug(f"Using direct search strategy for small genome/short sequence")
            result = cls._direct_search(sequence, genome_dict, allowed_mismatches)
            if result != (None, None):
                cls.validation_cache[cache_key] = result
                return result
        
        # STRATEGY 2: K-mer index search for larger genomes/sequences
        else:
            # Build index if needed
            if genome_id not in cls.sequence_index:
                # Choose k-mer size based on sequence and genome properties
                k_size = 8 if len(sequence) < 50 else 10
                if total_genome_size > 500_000_000:  # Very large genomes
                    k_size = 12
                cls.build_sequence_index(genome_dict, k=k_size)
            
            Config.debug(f"Using k-mer index search strategy")
            result = cls._kmer_search(sequence, genome_dict, genome_id, allowed_mismatches)
            if result != (None, None):
                cls.validation_cache[cache_key] = result
                return result
            
            # STRATEGY 3: Targeted substring search if k-mer fails
            Config.debug(f"K-mer search failed, trying targeted substring search")
            result = cls._targeted_search(sequence, genome_dict, allowed_mismatches)
            if result != (None, None):
                cls.validation_cache[cache_key] = result
                return result
        
        # Cache and return no match
        cls.validation_cache[cache_key] = (None, None)
        return None, None
    
    @classmethod
    def _direct_search(cls, sequence, genome_dict, allowed_mismatches):
        """Direct regex search for small genomes or short sequences."""
        # First try direct exact matching (much faster)
        for chrom, chrom_seq in genome_dict.items():
            if sequence in chrom_seq:
                pos = chrom_seq.find(sequence)
                return (chrom, pos + 1)  # 1-based position
        
        # Calculate max length for safe regex search
        seq_len = len(sequence)
        
        # For very long sequences, try substring search instead of regex
        if seq_len > 200:
            # Take a few strategic substrings from the sequence
            fragments = [
                sequence[:min(50, seq_len)],          # Start
                sequence[seq_len//2-25:seq_len//2+25] if seq_len > 50 else "", # Middle
                sequence[-min(50, seq_len):]          # End
            ]
            
            for frag in fragments:
                if not frag or len(frag) < 20:
                    continue
                    
                for chrom, chrom_seq in genome_dict.items():
                    frag_pos = chrom_seq.find(frag)
                    if frag_pos >= 0:
                        # Calculate where the full sequence would start
                        seq_start = frag_pos - sequence.find(frag)
                        if seq_start < 0 or seq_start + seq_len > len(chrom_seq):
                            continue
                            
                        # Extract and compare the full sequence
                        candidate = chrom_seq[seq_start:seq_start+seq_len]
                        dist = levenshtein_distance(sequence, candidate)
                        if dist <= allowed_mismatches:
                            return (chrom, seq_start + 1)  # 1-based
        
        # If fragments didn't work, try regex with longer timeout for shorter sequences
        if seq_len < 100:  # Only use regex for shorter sequences
            for chrom, chrom_seq in genome_dict.items():
                try:
                    Config.debug(f"Attempting regex search for sequence length {seq_len}")
                    matches = list(regex.finditer(
                        fr'({sequence}){{e<={allowed_mismatches}}}', 
                        chrom_seq, overlapped=True, timeout=5  # Increased timeout
                    ))
                    if matches:
                        return (chrom, matches[0].start() + 1)  # 1-based
                except (regex.error, TimeoutError) as e:
                    Config.debug(f"Regex search error: {str(e)[:100]}...")
                    continue  # Try next chromosome
        
        return None, None
    
    @classmethod
    def _kmer_search(cls, sequence, genome_dict, genome_id, allowed_mismatches):
        """Use k-mer index to find candidate positions, then verify."""
        index_data = cls.sequence_index[genome_id]
        k = index_data['k']
        index = index_data['index']
        
        # Extract strategic k-mers (beginning, middle, end)
        kmers_to_check = []
        
        # Beginning k-mers (more likely to match)
        for i in range(0, min(30, len(sequence) - k + 1), k):
            kmers_to_check.append((i, sequence[i:i+k]))
            
        # Middle k-mer
        if len(sequence) > 60:
            mid = len(sequence) // 2
            mid_pos = max(0, mid - k//2)
            kmers_to_check.append((mid_pos, sequence[mid_pos:mid_pos+k]))
            
        # Last k-mer (if not already covered)
        if len(sequence) > k*2:
            last_pos = len(sequence) - k
            kmers_to_check.append((last_pos, sequence[last_pos:last_pos+k]))
        
        # Find candidates from index
        candidates = {}
        for start_pos, kmer in kmers_to_check:
            if kmer in index:
                for chrom, pos in index[kmer]:
                    adj_pos = pos - start_pos
                    if adj_pos < 0:
                        continue
                    key = (chrom, adj_pos)
                    candidates[key] = candidates.get(key, 0) + 1
        
        # Check most promising positions first
        top_candidates = sorted(candidates.items(), key=lambda x: x[1], reverse=True)
        for (chrom, pos), _ in top_candidates[:5]:  # Check top 5 candidates
            if pos + len(sequence) <= len(genome_dict[chrom]):
                target = genome_dict[chrom][pos:pos+len(sequence)]
                
                # Try exact match first
                if sequence == target:
                    return (chrom, pos + 1)  # Exact match
                
                # Then check edit distance
                dist = levenshtein_distance(sequence, target)
                if dist <= allowed_mismatches:
                    return (chrom, pos + 1)  # 1-based position
        
        return None, None
    
    @classmethod
    def _targeted_search(cls, sequence, genome_dict, allowed_mismatches):
        """Search using strategic substrings of the sequence for when k-mer and direct search fail."""
        import random
        
        # Create several non-overlapping substrings to search for
        # Use longer substrings for better specificity
        substrings = []
        seq_len = len(sequence)
        
        # Extract exact substrings from beginning, middle, end with more coverage
        if seq_len >= 100:
            # For longer sequences, use more substrings
            fifth = seq_len // 5
            substrings.append(sequence[:min(30, fifth)])         # First segment
            substrings.append(sequence[fifth:fifth*2])           # Second segment
            substrings.append(sequence[fifth*2:fifth*3])         # Middle segment
            substrings.append(sequence[fifth*3:fifth*4])         # Fourth segment
            substrings.append(sequence[max(0, seq_len-30):])     # Last segment
        elif seq_len >= 60:
            third = seq_len // 3
            substrings.append(sequence[:min(25, third)])
            substrings.append(sequence[third:third*2])
            substrings.append(sequence[max(0, seq_len-25):])
        elif seq_len >= 30:
            half = seq_len // 2
            substrings.append(sequence[:min(20, half)])
            substrings.append(sequence[max(0, seq_len-20):])
        else:
            # For short sequences, use the whole sequence
            substrings.append(sequence)
        
        # Filter out any substrings that are too short to be useful
        substrings = [s for s in substrings if len(s) >= 15]
        if not substrings:
            substrings = [sequence]  # Fallback to using the full sequence
        
        # Select chromosomes to search (stratified approach)
        search_chroms = list(genome_dict.keys())
        if len(search_chroms) > 5:
            # Take 5 chromosomes: first, last, middle, and two random ones
            mid_idx = len(search_chroms) // 2
            sampled = [search_chroms[0], search_chroms[-1], search_chroms[mid_idx]]
            
            # Add two random chromosomes that aren't already selected
            remaining = [c for c in search_chroms if c not in sampled]
            if remaining:
                sampled.extend(random.sample(remaining, min(2, len(remaining))))
            search_chroms = sampled
        
        # Search each chromosome
        for chrom in search_chroms:
            chrom_seq = genome_dict[chrom]
            
            # For very large chromosomes, search strategic windows
            if len(chrom_seq) > Config.MAX_SEARCH_LENGTH:
                windows = []
                total_length = len(chrom_seq)
                
                # Create windows at the beginning, middle, and end of chromosome
                windows.append((0, min(Config.MAX_SEARCH_LENGTH, total_length)))
                
                if total_length > Config.MAX_SEARCH_LENGTH:
                    mid_start = max(0, (total_length // 2) - (Config.MAX_SEARCH_LENGTH // 2))
                    mid_end = min(total_length, mid_start + Config.MAX_SEARCH_LENGTH)
                    windows.append((mid_start, mid_end))
                
                if total_length > 2 * Config.MAX_SEARCH_LENGTH:
                    end_start = max(0, total_length - Config.MAX_SEARCH_LENGTH)
                    windows.append((end_start, total_length))
                
                # Search each window
                for start, end in windows:
                    window = chrom_seq[start:end]
                    
                    # Search each substring in window
                    for substr in substrings:
                        pos = window.find(substr)
                        if pos >= 0:
                            substr_offset = sequence.find(substr)
                            seq_start = start + pos - substr_offset
                            
                            # Check boundaries
                            if seq_start < 0 or seq_start + len(sequence) > len(chrom_seq):
                                continue
                            
                            potential_match = chrom_seq[seq_start:seq_start + len(sequence)]
                            dist = levenshtein_distance(sequence, potential_match)
                            if dist <= allowed_mismatches:
                                return (chrom, seq_start + 1)  # 1-based
            else:
                # For smaller chromosomes, check each substring directly
                for substr in substrings:
                    pos = chrom_seq.find(substr)
                    if pos >= 0:
                        substr_offset = sequence.find(substr)
                        seq_start = pos - substr_offset
                        
                        if seq_start < 0 or seq_start + len(sequence) > len(chrom_seq):
                            continue
                        
                        potential_match = chrom_seq[seq_start:seq_start + len(sequence)]
                        dist = levenshtein_distance(sequence, potential_match)
                        if dist <= allowed_mismatches:
                            return (chrom, seq_start + 1)  # 1-based
        
        # If all else fails, try limited regex with much higher timeout on a small sample
        try:
            # Choose one random chromosome
            if search_chroms:
                random_chrom = random.choice(search_chroms)
                chrom_seq = genome_dict[random_chrom]
                
                # For large chromosomes, use a random window
                if len(chrom_seq) > Config.MAX_SEARCH_LENGTH:
                    start = random.randint(0, max(0, len(chrom_seq) - Config.MAX_SEARCH_LENGTH))
                    search_seq = chrom_seq[start:start + Config.MAX_SEARCH_LENGTH]
                    
                    Config.debug(f"Last resort regex search with longer timeout")
                    matches = list(regex.finditer(
                        fr'({sequence}){{e<={allowed_mismatches}}}', 
                        search_seq, overlapped=True, timeout=10  # Much longer timeout
                    ))
                    if matches:
                        return (random_chrom, start + matches[0].start() + 1)
                elif len(chrom_seq) < 10_000_000:  # Only try regex on manageable sequences
                    matches = list(regex.finditer(
                        fr'({sequence}){{e<={allowed_mismatches}}}', 
                        chrom_seq, overlapped=True, timeout=10  # Longer timeout
                    ))
                    if matches:
                        return (random_chrom, matches[0].start() + 1)
        except (regex.error, IndexError, TimeoutError) as e:
            Config.debug(f"Final fallback regex search failed: {str(e)[:100]}...")
        
        return None, None

    @classmethod
    def find_exact_match(cls, sequence, genome_dict):
        """
        Optimized method to find exact matches of a sequence in a genome.
        
        Args:
            sequence (str): Sequence to search for
            genome_dict (dict): Dictionary mapping chromosome names to sequences
            
        Returns:
            bool: True if found, False otherwise
        """
        if not sequence or len(sequence) < 10:
            return False
        
        # Generate cache key
        genome_id = hash(tuple(sorted(genome_dict.keys())))
        cache_key = ("exact", sequence, genome_id)
        
        # Check cache
        if cache_key in cls.validation_cache:
            return cls.validation_cache[cache_key]
        
        # Check if we have an index for this genome
        if genome_id in cls.sequence_index:
            index_data = cls.sequence_index[genome_id]
            k = index_data['k']
            index = index_data['index']
            
            # For sequences shorter than k, use direct search
            if len(sequence) < k:
                for chrom_seq in genome_dict.values():
                    if sequence in chrom_seq:
                        cls.validation_cache[cache_key] = True
                        return True
                        
                cls.validation_cache[cache_key] = False
                return False
            
            # For longer sequences, use the index
            kmer = sequence[:k]
            if kmer not in index:
                cls.validation_cache[cache_key] = False
                return False
                
            # Check each potential location
            for chrom, pos in index[kmer]:
                # Ensure we don't go out of bounds
                if pos + len(sequence) <= len(genome_dict[chrom]):
                    if genome_dict[chrom][pos:pos+len(sequence)] == sequence:
                        cls.validation_cache[cache_key] = True
                        return True
            
            cls.validation_cache[cache_key] = False
            return False
        
        # Fall back to direct search if no index
        for chrom_seq in genome_dict.values():
            if sequence in chrom_seq:
                cls.validation_cache[cache_key] = True
                return True
                
        cls.validation_cache[cache_key] = False
        return False

    @classmethod
    def validate_sequence(cls, seq_info, ref_genome, qry_genome, extra_genomes=None):
        """
        Validate a sequence against multiple genomes with optimized performance.
        Adapts validation approach based on configuration.
        
        Args:
            seq_info (dict): Sequence information dictionary
            ref_genome (dict): Reference genome dictionary
            qry_genome (dict): Query genome dictionary
            extra_genomes (list): List of additional genome dictionaries
            
        Returns:
            bool: True if sequence is valid across all genomes, False otherwise
        """
        if extra_genomes is None:
            extra_genomes = []
        
        # Extract sequence info
        seq = seq_info["sequence"]
        primer_f = seq_info.get("Primer F", "")
        primer_r = seq_info.get("Primer R", "")
        probe = seq_info.get("Probe", "")
        
        # Generate a stable cache key that works across processes
        validation_key = (seq, hash(tuple([hash(tuple(g.keys())) for g in [ref_genome, qry_genome] + extra_genomes])))
        
        # Quick check in cache
        if validation_key in cls.validation_cache:
            return cls.validation_cache[validation_key]
        
        # Flag to track validation progress
        is_valid = True
        
        # Check mode-specific requirements
        if Config.VALIDATION_MODE == "STRICT":
            # In strict mode, everything must match exactly
            # Check primers first (cheapest validation)
            if primer_f and not cls.find_exact_match(primer_f, ref_genome):
                is_valid = False
            
            if is_valid and primer_r and not cls.find_exact_match(primer_r, qry_genome):
                is_valid = False
                
            # Check amplicon next
            if is_valid:
                for genome in [ref_genome, qry_genome] + extra_genomes:
                    if not cls.find_exact_match(seq, genome):
                        is_valid = False
                        break
            
            # Check probe last (if we're still valid)
            if is_valid and probe:
                for genome in [ref_genome, qry_genome] + extra_genomes:
                    if not cls.find_exact_match(probe, genome):
                        is_valid = False
                        break
        else:  # TOLERANT mode
            # Primers and probes must match exactly, amplicons can have mismatches
            # Check primers first
            if primer_f:
                primer_f_valid = cls.find_exact_match(primer_f, ref_genome) and \
                                cls.find_exact_match(primer_f, qry_genome) and \
                                all(cls.find_exact_match(primer_f, g) for g in extra_genomes)
                if not primer_f_valid:
                    is_valid = False
            
            if is_valid and primer_r:
                primer_r_valid = cls.find_exact_match(primer_r, ref_genome) and \
                                cls.find_exact_match(primer_r, qry_genome) and \
                                all(cls.find_exact_match(primer_r, g) for g in extra_genomes)
                if not primer_r_valid:
                    is_valid = False
            
            # Check probe
            if is_valid and probe:
                probe_valid = cls.find_exact_match(probe, ref_genome) and \
                            cls.find_exact_match(probe, qry_genome) and \
                            all(cls.find_exact_match(probe, g) for g in extra_genomes)
                if not probe_valid:
                    is_valid = False
            
            # Check amplicon with mismatches
            if is_valid:
                ref_match = cls.matches_with_tolerance(seq, ref_genome) != (None, None)
                if not ref_match:
                    is_valid = False
                
                if is_valid:
                    qry_match = cls.matches_with_tolerance(seq, qry_genome) != (None, None)
                    if not qry_match:
                        is_valid = False
                
                if is_valid and extra_genomes:
                    for genome in extra_genomes:
                        if cls.matches_with_tolerance(seq, genome) == (None, None):
                            is_valid = False
                            break
        
        # Cache and return result
        cls.validation_cache[validation_key] = is_valid
        return is_valid

    @classmethod
    def validate_sequences_batch(cls, batch, ref_genome, qry_genome, extra_genomes=None):
        """
        Process a batch of sequences in parallel with optimized validation.
        
        Args:
            batch (list): List of sequence dictionaries to validate
            ref_genome (dict): Reference genome dictionary
            qry_genome (dict): Query genome dictionary
            extra_genomes (list): List of additional genome dictionaries
            
        Returns:
            list: List of valid sequence dictionaries
        """
        if extra_genomes is None:
            extra_genomes = []
        
        # Build indices if needed (once per genome)
        if Config.VALIDATION_MODE == "TOLERANT":
            cls.build_sequence_index(ref_genome)
            cls.build_sequence_index(qry_genome)
            for genome in extra_genomes:
                cls.build_sequence_index(genome)
        
        # Process in smaller sub-batches to increase reliability
        sub_batch_size = min(25, len(batch))
        sub_batches = list(Utils.chunks(batch, sub_batch_size))
        
        start_time = time.time()
        Config.debug(f"Validating {len(batch)} sequences in {len(sub_batches)} sub-batches")
        
        results = []
        for i, sub_batch in enumerate(sub_batches):
            Config.debug(f"Processing sub-batch {i+1}/{len(sub_batches)} with {len(sub_batch)} sequences")
            
            # Use ThreadPoolExecutor for shared memory advantage with the indices
            with concurrent.futures.ThreadPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
                # Process each sequence with robust error handling
                futures = {}
                for seg in sub_batch:
                    future = executor.submit(cls._safe_validate_sequence, 
                                           seg, ref_genome, qry_genome, extra_genomes)
                    futures[future] = seg
                
                # Collect valid sequences with proper error handling
                for future in concurrent.futures.as_completed(futures):
                    segment = futures[future]
                    try:
                        is_valid = future.result()
                        seq_preview = segment['sequence'][:20] + "..." if len(segment['sequence']) > 20 else segment['sequence']
                        
                        if is_valid:
                            results.append(segment)
                            Config.debug(f"Sequence valid: {seq_preview}")
                        else:
                            Config.debug(f"Sequence invalid: {seq_preview}")
                    except Exception as e:
                        Config.debug(f"Error validating sequence: {str(e)[:100]}... - skipping")
                
        elapsed = time.time() - start_time
        Config.debug(f"Batch validation completed: {len(results)}/{len(batch)} valid in {elapsed:.2f}s")
        return results
    
    @classmethod
    def _safe_validate_sequence(cls, seq_info, ref_genome, qry_genome, extra_genomes=None):
        """Wrapper for validate_sequence with exception handling."""
        try:
            return cls.validate_sequence(seq_info, ref_genome, qry_genome, extra_genomes)
        except Exception as e:
            Config.debug(f"Validation error: {str(e)[:100]}...")
            return False

    @classmethod
    def clear_caches(cls):
        """Clear all caches to free memory."""
        cls.validation_cache.clear()
        Config.debug("Validation caches cleared")


##############################################################################
#                          Primer3 Processing
##############################################################################
class Primer3Processor:
    """Handles Primer3 operations for primer design."""
    
    @staticmethod
    def run_primer3_batch(input_blocks):
        """
        Run primer3_core on a batch of input blocks.
        Checks for primer3_core availability before execution.
        """
        primer3_input_str = "".join(input_blocks)

        if shutil.which("primer3_core") is None:
            print("Error: primer3_core is not available in PATH.")
            return ""

        try:
            primer3_proc = subprocess.Popen(
                ["primer3_core"],
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
        
        Args:
            stdout_data (str): Primer3 stdout output
            
        Returns:
            list: List of primer record dictionaries
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
                
                # Check if we need to reverse complement the probe based on C/G content
                probe_reversed = False
                if Config.PREFER_PROBE_MORE_C_THAN_G and probe_seq:
                    probe_seq, probe_reversed = Utils.ensure_more_c_than_g(probe_seq)
                    
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
                    "Probe Reversed": probe_reversed,  # Add this field to track reversals
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
        Find the locations of a batch of amplicons in multiple genomes.
        Uses improved search strategies for better matching.
        
        Args:
            batch_data (tuple): Tuple of (amplicon_seqs, genome_dicts)
            
        Returns:
            list: List of location tuples
        """
        amplicon_seqs, genome_dicts = batch_data
        results = []

        for seq in amplicon_seqs:
            if pd.notnull(seq):
                # Since we're processing one genome at a time
                if isinstance(genome_dicts, dict):
                    if len(genome_dicts) == 1:
                        genome_dict = next(iter(genome_dicts.values()))  # Use the single genome dictionary
                    else:
                        genome_dict = genome_dicts  # Pass all genome dictionaries
                else:
                    # Ensure genome_dicts is treated as a dictionary
                    genome_dict = {"default": genome_dicts} if isinstance(genome_dicts, str) else genome_dicts
                
                # Try primary exact match first (faster)
                primary_match = False
                for chrom, chrom_seq in genome_dict.items():
                    if seq in chrom_seq:
                        pos = chrom_seq.find(seq) + 1  # Convert to 1-based position
                        results.append((chrom, pos))
                        primary_match = True
                        break
                
                if not primary_match:
                    # Try alternate strategies if exact match fails
                    # Strategy 1: Check both ends of the amplicon (for primer regions)
                    if len(seq) >= 40:
                        left_end = seq[:20]
                        right_end = seq[-20:]
                        
                        best_match = None
                        
                        for chrom, chrom_seq in genome_dict.items():
                            left_pos = chrom_seq.find(left_end)
                            if left_pos >= 0:
                                # Check if a matching region of appropriate length follows
                                potential_region = chrom_seq[left_pos:left_pos + len(seq)]
                                if len(potential_region) >= len(seq) * 0.9:  # Allow for some mismatch
                                    # Count mismatches
                                    mismatches = sum(1 for a, b in zip(seq, potential_region) if a != b)
                                    if mismatches <= Config.ALLOW_AMP_MISMATCHES:
                                        best_match = (chrom, left_pos + 1)
                                        break  # Good enough match
                            
                            # Try right end if left end didn't work
                            right_pos = chrom_seq.find(right_end)
                            if right_pos >= 0:
                                # Check if a matching region of appropriate length precedes
                                start_pos = max(0, right_pos - (len(seq) - len(right_end)))
                                potential_region = chrom_seq[start_pos:right_pos + len(right_end)]
                                if len(potential_region) >= len(seq) * 0.9:  # Allow for some mismatch
                                    # Count mismatches
                                    aligned_region = potential_region[-len(seq):] if len(potential_region) > len(seq) else potential_region
                                    padded_region = aligned_region.ljust(len(seq), 'X')  # Pad if needed
                                    mismatches = sum(1 for a, b in zip(seq, padded_region) if a != b)
                                    if mismatches <= Config.ALLOW_AMP_MISMATCHES:
                                        best_match = (chrom, start_pos + 1)
                                        break  # Good enough match
                        
                        if best_match:
                            results.append(best_match)
                        else:
                            # Strategy 2: Try the Levenshtein fuzzy search as fallback
                            location = SequenceValidator.matches_with_tolerance(seq, genome_dict)
                            results.append(location)
                    else:
                        # For short sequences, just use the original fuzzy search
                        location = SequenceValidator.matches_with_tolerance(seq, genome_dict)
                        results.append(location)
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
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: Minimum free energy, or None if calculation fails
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
        
        Args:
            seqs (list): List of DNA sequences
            
        Returns:
            list: List of ΔG values
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
        """
        Split a list into chunks of size n.
        
        Args:
            lst (list): List to split
            n (int): Chunk size
            
        Yields:
            list: Chunks of the original list
        """
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    @staticmethod
    def flatten_list(list_of_lists):
        """
        Flatten a list of lists into a single list.
        
        Args:
            list_of_lists (list): List of lists
            
        Returns:
            list: Flattened list
        """
        return [item for sublist in list_of_lists for item in sublist]

    @staticmethod
    def has_disallowed_repeats(seq):
        """
        Check for disallowed repeats in a sequence.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            bool: True if disallowed repeats found, False otherwise
        """
        if not isinstance(seq, str):
            return True
        return "CCCC" in seq or "GGGG" in seq

    @staticmethod
    def calculate_gc(seq):
        """
        Calculate GC content of a sequence.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: GC content percentage
        """
        if not seq or not isinstance(seq, str):
            return 0
        seq = seq.upper()
        gc_count = sum(1 for base in seq if base in "GC")
        return (gc_count / len(seq)) * 100 if seq else 0

    @staticmethod
    def ensure_more_c_than_g(seq):
        """
        Check if a sequence has more Cs than Gs. If not, return the reverse complement.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            tuple: (possibly_reversed_sequence, was_reversed)
        """
        if not seq or not isinstance(seq, str):
            return seq, False
            
        seq = seq.upper()
        c_count = seq.count('C')
        g_count = seq.count('G')
        
        if c_count >= g_count:
            return seq, False  # No need to reverse
        
        # Need to reverse complement
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                    'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 
                    'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 
                    'D': 'H', 'H': 'D', 'V': 'B'}
        rev_comp = ''.join(complement.get(base, base) for base in reversed(seq))
        return rev_comp, True


##############################################################################
#                          BLAST Functions
##############################################################################
class BlastProcessor:
    """Handles BLAST operations for primer specificity checking."""
    
    @staticmethod
    def blast_short_seq(seq, db=None):
        """
        Runs BLASTn for short sequences and returns the two lowest e-values separately.
        
        Args:
            seq (str): Sequence to BLAST
            db (str): BLAST database path (default: from Config)
            
        Returns:
            tuple: (best_evalue, second_best_evalue)
        """
        if db is None:
            db = f'"{Config.DB_PATH}"'
            
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
                    "blastn",
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
        
        Args:
            batch_data (tuple): Tuple of (batch, col_name)
            
        Returns:
            list: List of (best_evalue, second_best_evalue) tuples
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
        
        Args:
            row (pandas.Series): DataFrame row
            col (str): Column prefix
            
        Returns:
            bool: True if passes filter, False otherwise
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


    @classmethod
    def validate_sequences_batch(cls, batch, ref_genome, qry_genome, extra_genomes=None):
        """
        Process a batch of sequences in parallel with optimized validation.
        
        Args:
            batch (list): List of sequence dictionaries to validate
            ref_genome (dict): Reference genome dictionary
            qry_genome (dict): Query genome dictionary
            extra_genomes (list): List of additional genome dictionaries
            
        Returns:
            list: List of valid sequence dictionaries
        """
        if extra_genomes is None:
            extra_genomes = []
        
        # Build indices if needed (once per genome)
        if Config.VALIDATION_MODE == "TOLERANT":
            cls.build_sequence_index(ref_genome)
            cls.build_sequence_index(qry_genome)
            for genome in extra_genomes:
                cls.build_sequence_index(genome)
            
        start_time = time.time()
        Config.debug(f"Validating batch of {len(batch)} sequences...")
        
        # Use ThreadPoolExecutor for shared memory advantage with the indices
        with concurrent.futures.ThreadPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
            # Process each sequence
            future_to_seg = {
                executor.submit(cls.validate_sequence, seg, ref_genome, qry_genome, extra_genomes): seg 
                for seg in batch
            }
            
            # Collect valid sequences
            valid_sequences = []
            for future in concurrent.futures.as_completed(future_to_seg):
                segment = future_to_seg[future]
                try:
                    is_valid = future.result()
                    if is_valid:
                        valid_sequences.append(segment)
                        
                    seq_preview = segment['sequence'][:30] + "..." if len(segment['sequence']) > 30 else segment['sequence']
                    Config.debug(f"Sequence validation result: {seq_preview} - {'Valid' if is_valid else 'Invalid'}")
                except Exception as e:
                    Config.debug(f"Error validating sequence: {e}")
        
        elapsed = time.time() - start_time
        Config.debug(f"Batch validation completed: {len(valid_sequences)}/{len(batch)} valid in {elapsed:.2f}s")
        return valid_sequences


##############################################################################
#                          Match Checker Implementation  
##############################################################################
class MatchChecker:
    """Handles sequence matching between primers and genomes based on the Match checker.py logic."""
    
    @staticmethod
    def run_match_checker(primer_file, genome_file, min_match_length=None):
        """
        Run the Match Checker functionality.
        
        Args:
            primer_file (str): Path to the Excel or CSV file with primers
            genome_file (str): Path to the FASTA file with genome sequences
            min_match_length (int): Minimum length of match to report
            
        Returns:
            str: Path to the output Excel file
        """
        if min_match_length is None:
            min_match_length = Config.MIN_SEGMENT_LENGTH
            
        print(f"Running Match Checker with minimum match length: {min_match_length}")
        
        # Load genome sequences
        genome_seqs = {}
        for record in SeqIO.parse(genome_file, "fasta"):
            genome_seqs[record.id] = str(record.seq).upper()
        
        # Load query sequences
        if (primer_file.endswith(".csv")):
            df = pd.read_csv(primer_file)
        else:
            df = pd.read_excel(primer_file)
        
        # Assume we need the first two columns (gene name and sequence)
        if len(df.columns) < 2:
            print("Error: Primer file must have at least two columns (gene name and sequence).")
            sys.exit(1)
            
        # Create query_data as a list of (gene_name, sequence) tuples
        query_data = df.iloc[:, :2].dropna().astype(str).values.tolist()
        
        results = []
        for gene_name, seq in tqdm(query_data, desc="Processing Sequences"):
            chrom, pos, match_len = MatchChecker.find_exact_subsequence(seq, genome_seqs, min_match_length)
            if match_len > 0:  # Only store results where a match is found
                results.append((gene_name, seq, chrom, pos, match_len, match_len == len(seq)))
        
        # Create a dataframe to store results
        df_results = pd.DataFrame(results, columns=["Gene Name", "Sequence", "Chromosome", "Start", "Match Length", "Full Match"])
        
        # Output to Excel in the Primers directory (same as main pipeline)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        primers_dir = os.path.join(script_dir, "Primers")
        os.makedirs(primers_dir, exist_ok=True)
        output_file = os.path.join(primers_dir, "Subsequence_Match_Results.xlsx")
        df_results.to_excel(output_file, index=False)
        
        print(f"Subsequence matching completed. Results saved to {output_file}")
        return output_file
    
    @staticmethod
    def find_exact_subsequence(seq, genome_seqs, min_match_length=100):
        """
        Finds the longest exact subsequence match for a given query sequence in the genome sequences.
        
        Parameters:
        - seq (str): The query sequence.
        - genome_seqs (dict): Dictionary of genome sequences with chromosome names as keys.
        - min_match_length (int): Minimum length of the matching subsequence.
        
        Returns:
        - chrom (str): The chromosome where the match is found.
        - pos (int): The start position of the match (1-based).
        - match_len (int): Length of the matched sequence.
        """
        best_match = (None, None, 0)  # (chromosome, position, match length)
        seq_len = len(seq)

        for length in range(seq_len, min_match_length - 1, -1):  # Start from full sequence length
            for start in range(seq_len - length + 1):  # Iterate over possible start positions
                sub_seq = seq[start:start + length]
                for chrom, genome_seq in genome_seqs.items():
                    pos = genome_seq.find(sub_seq)
                    if pos != -1:
                        return chrom, pos + 1, length  # Return immediately when the longest match is found

        return best_match  # Return None if no match is found

    @staticmethod
    def run_from_ddprimer_output(output_file, genome_file):
        """
        Run Match Checker on the output of the ddPrimer pipeline.
        
        Args:
            output_file (str): Path to the Excel file from ddPrimer
            genome_file (str): Path to the FASTA file to check against
            
        Returns:
            str: Path to the output Excel file
        """
        print(f"Running Match Checker on ddPrimer output: {output_file}")
        
        # Load ddPrimer output
        df = pd.read_excel(output_file)
        
        # Check if this is a direct pipeline output or already preprocessed
        if "Amplicon" in df.columns:
            # Direct pipeline output - create a new dataframe with gene names and amplicon sequences
            # Use both Sequence and Index to create a unique identifier
            query_df = pd.DataFrame({
                "Gene Name": df["Sequence"] + "_" + df["Index"].astype(str),
                "Sequence": df["Amplicon"]
            })
        else:
            # Assume it's already in the right format (e.g., from a previous Match Checker run)
            query_df = df.copy()
            if "Gene Name" not in query_df.columns or "Sequence" not in query_df.columns:
                print("Error: Input file doesn't have required columns (Gene Name, Sequence)")
                return None
        
        # Filter out rows with empty sequences
        query_df = query_df[~query_df["Sequence"].isna()]
        
        # Make sure we have data to process
        if len(query_df) == 0:
            print("Error: No valid sequences found in the input file.")
            return None
        
        # Save this as a temporary file
        temp_file = os.path.join(os.path.dirname(output_file), "temp_match_checker_input.xlsx")
        query_df.to_excel(temp_file, index=False)
        
        # Run the match checker
        result = MatchChecker.run_match_checker(temp_file, genome_file)
        
        # Clean up the temporary file
        try:
            os.remove(temp_file)
        except:
            print(f"Warning: Could not remove temporary file {temp_file}")
        
        return result


##############################################################################
#                          PLastZ Only Implementation
##############################################################################
class LastZOnlyRunner:
    """Implements PLastZ-only functionality based on the provided script."""
    
    @staticmethod
    def run_lastz_only(new_sequence, existing_sequences, output_dir="Alignments", num_processes=None):
        """
        Run PLastZ alignments between a new sequence and a list of existing sequences.
        
        Args:
            new_sequence (str): Path to the new FASTA file
            existing_sequences (list): List of paths to existing FASTA files
            output_dir (str): Directory to store alignments
            num_processes (int): Number of processes to use (defaults to Config.NUM_PROCESSES)
            
        Returns:
            list: Paths to the generated MAF files
        """
        if num_processes is None:
            num_processes = Config.NUM_PROCESSES
            
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Create an instance of LastZRunner from the main pipeline
        lastz_runner = LastZRunner(config=Config)
        
        output_files = []
        
        for existing_seq in existing_sequences:
            # Create alignment output name
            alignment_name = (f"{os.path.splitext(os.path.basename(new_sequence))[0]}_"
                             f"vs_{os.path.splitext(os.path.basename(existing_seq))[0]}.maf")
            output_path = os.path.join(output_dir, alignment_name)
            
            # Prepare LastZ options
            lastz_options = "--format=maf --ambiguous=iupac"
            
            print(f"\nAligning {os.path.basename(new_sequence)} vs {os.path.basename(existing_seq)}... ")
            stop_event = threading.Event()
            thread = threading.Thread(target=UIUtilities.spinner_task, args=(stop_event,))
            thread.start()
            
            try:
                # Create a temporary directory for this pair
                pair_temp_dir = os.path.join(output_dir, "temp_" + os.path.basename(output_path))
                
                # Run the alignment
                temp_output = lastz_runner.run_parallel_alignment(
                    existing_seq, 
                    new_sequence,
                    pair_temp_dir,
                    lastz_options=lastz_options,
                    processes=num_processes
                )
                
                # Rename and move the output file
                if os.path.exists(temp_output):
                    if os.path.exists(output_path):
                        os.remove(output_path)
                    shutil.copy(temp_output, output_path)
                    
                    # Clean up the temporary directory
                    if os.path.exists(pair_temp_dir):
                        shutil.rmtree(pair_temp_dir)
                    
                    output_files.append(output_path)
                    print(f"Saved output to {output_path}")
                else:
                    print(f"\nPLastZ did not generate the expected output file")
                    
            except Exception as e:
                print(f"\nPLastZ alignment failed with error: {e}")
            finally:
                stop_event.set()
                thread.join()
        
        print("All alignments completed.")
        return output_files


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
        """
        Get initial user inputs for the pipeline.
        
        Returns:
            str: MAF choice (1 for run LastZ, 2 for use existing MAF files)
        """
        print("=== Combined Primer Design Pipeline ===")
        if FileUtilities.use_cli:
            print("Detected CLI environment. Running in text-based mode.")
        
        # Handle Match Checker mode separately
        if Config.RUN_MATCH_CHECKER:
            print("Running in Match Checker mode...")
            print("Select the genome sequence FASTA File")
            genome_file = FileUtilities.get_file("Select Genome FASTA File", 
                                            [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])
            
            print("Select the primer/sequence file (Excel or CSV)")
            primer_file = FileUtilities.get_file("Select Primer/Sequence File", 
                                            [("Excel Files", "*.xlsx"), ("CSV Files", "*.csv"), ("All Files", "*.*")])
            
            # Run Match Checker
            output_file = MatchChecker.run_match_checker(primer_file, genome_file, Config.MIN_SEGMENT_LENGTH)
            
            print(f"Match Checker completed. Results saved to {output_file}")
            sys.exit(0)
        
        # For both full pipeline and LastZ/MultiZ-only mode, get the same input initially
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
            print(f"Select {num_extra} additional genomes for filtering")
            for i in range(num_extra):
                extra_path = FileUtilities.get_file(f"Select {i+1}. Additional FASTA File", 
                                                    [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])
                self.extra_genome_paths.append(extra_path)
                self.extra_genomes.append(FileUtilities.load_fasta(extra_path))

        # Load GFF3 annotation file
        print("Select the GFF3 Annotation File")
        self.gff_path = FileUtilities.get_file("Select GFF3 Annotation File", 
                                             [("GFF3 Files", "*.gff3 *.gff"), ("All Files", "*.*")])
        
        # Process GFF file
        self.genes = GFFProcessor.load_genes_from_gff(self.gff_path)
        
        # Load reference and query genomes
        self.ref_genome = FileUtilities.load_fasta(self.ref_path)
        self.qry_genome = FileUtilities.load_fasta(self.qry_path)
    
    def run_lastz(self, maf_choice):
        """
        Run parallel LastZ alignments using the integrated LastZRunner.
        
        Args:
            maf_choice (str): MAF choice (1 for run LastZ, 2 for use existing MAF files)
        """
        if (maf_choice == "1"):
            # Use the integrated LastZRunner for alignments
            lastz_runner = LastZRunner(config=Config)
            lastz_output_files = []

            # Use already selected genomes: reference, query, and any additional genomes
            genome_paths = [self.ref_path, self.qry_path] + self.extra_genome_paths
            
            # Run pairwise alignments for all genome combinations
            for ref_genome_path, qry_genome_path in combinations(genome_paths, 2):
                # Get the directory of the input files
                input_dir = os.path.dirname(ref_genome_path)
                
                # Create an "Alignments" folder in the input directory if it doesn't exist
                alignments_dir = os.path.join(input_dir, "Alignments")
                os.makedirs(alignments_dir, exist_ok=True)
                
                # Create alignment output name
                alignment_name = (f"{os.path.splitext(os.path.basename(ref_genome_path))[0]}_"
                                f"vs_{os.path.splitext(os.path.basename(qry_genome_path))[0]}.maf")
                lastz_output_file = os.path.join(alignments_dir, alignment_name)
                
                # Prepare LastZ options
                lastz_options = "--format=maf --ambiguous=iupac"
                
                print(f"\nAligning {os.path.basename(ref_genome_path)} vs {os.path.basename(qry_genome_path)}... ")
                stop_event = threading.Event()
                thread = threading.Thread(target=UIUtilities.spinner_task, args=(stop_event,))
                thread.start()
                
                try:
                    # Create a temporary directory for this pair
                    pair_temp_dir = os.path.join(alignments_dir, "temp_" + os.path.basename(lastz_output_file))
                    
                    # Run the alignment
                    temp_output = lastz_runner.run_parallel_alignment(
                        ref_genome_path, 
                        qry_genome_path,
                        pair_temp_dir,
                        lastz_options=lastz_options,
                        processes=Config.NUM_PROCESSES
                    )
                    
                    # Rename and move the output file
                    if os.path.exists(temp_output):
                        if os.path.exists(lastz_output_file):
                            os.remove(lastz_output_file)
                        shutil.copy(temp_output, lastz_output_file)
                        
                        # Clean up the temporary directory
                        if os.path.exists(pair_temp_dir):
                            shutil.rmtree(pair_temp_dir)
                        
                        lastz_output_files.append(lastz_output_file)
                        print(f"Saved output to {lastz_output_file}")
                    else:
                        print(f"\nPLastZ did not generate the expected output file")
                        sys.exit(1)
                        
                except Exception as e:
                    print(f"\nPLastZ alignment failed with error: {e}")
                    sys.exit(1)
                finally:
                    stop_event.set()
                    thread.join()

            self.maf_files = lastz_output_files
        else:
            # Use existing MAF files
            print("Using existing MAF files for analysis")
    
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
            if os.path.exists(multiz_output_file):
                print(f"⚠️ Warning: Output file '{multiz_output_file}' already exists and may cause issues with MultiZ.")
                print("Consider deleting it manually or renaming it before running MultiZ.")
            
            # Order MAF files to match the same order they would have been generated
            # Generate list of expected genome pairs based on combinations
            genome_paths = [self.ref_path, self.qry_path] + self.extra_genome_paths
            expected_pairs = list(combinations(genome_paths, 2))
            
            # Function to get base filename without extension
            def get_base_filename(path):
                return os.path.splitext(os.path.basename(path))[0]
            
            # Create dictionary of expected filename patterns
            expected_patterns = {}
            for ref_genome_path, qry_genome_path in expected_pairs:
                ref_name = get_base_filename(ref_genome_path)
                qry_name = get_base_filename(qry_genome_path)
                # Both possible orderings in the filename
                expected_patterns[(ref_genome_path, qry_genome_path)] = [
                    f"{ref_name}_vs_{qry_name}",
                    f"{qry_name}_vs_{ref_name}"
                ]
            
            # Try to match each MAF file to a genome pair based on filename
            ordered_maf_files = []
            unmatched_maf_files = list(self.maf_files)  # Copy to track unmatched files
            
            # First pass: match by filename patterns
            for pair in expected_pairs:
                patterns = expected_patterns[pair]
                match_found = False
                
                for maf_file in list(unmatched_maf_files):  # Use a copy for safe iteration
                    maf_basename = get_base_filename(maf_file)
                    
                    # Check if this MAF file matches any expected pattern
                    for pattern in patterns:
                        if pattern in maf_basename:
                            ordered_maf_files.append(maf_file)
                            unmatched_maf_files.remove(maf_file)
                            match_found = True
                            break
                    
                    if match_found:
                        break
            
            # If any files remain unmatched, add them at the end
            # In a real implementation, we would analyze MAF contents here
            if unmatched_maf_files:
                print("  Warning: Some MAF files could not be matched by name pattern.")
                print("  Appending unmatched files at the end. Results may be inconsistent.")
                ordered_maf_files.extend(unmatched_maf_files)
            
            # Use the ordered files for MultiZ
            command = f'"{Config.MULTIZ_PATH}" ' + " ".join(f'"{file}"' for file in ordered_maf_files) + f' > "{multiz_output_file}"'
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
            for seg in maf_segments:
                # Standardize sequence ID and index format
                seg['id'] = f"{seg['gene']}_{seg['index']}"
            self.segments.extend(maf_segments)

        if not self.segments:
            print("Error: No valid segments extracted from MAF.")
            sys.exit(1)

        print(f"\nExtracted {len(self.segments)} initial segments from MAF files.")
        
        # Debug: Print examples of segments with their coordinates
        if Config.DEBUG_MODE:
            print("\nDEBUG: Examples of extracted segments with coordinates:")
            for i, seg in enumerate(self.segments[:5]):  # Print first 5 segments
                print(f"Segment {i+1}:")
                print(f"  Gene: {seg['gene']}, Index: {seg['index']}")
                print(f"  Sequence start: {seg['sequence'][:30]}...")
                if 'ref_chr' in seg and 'ref_pos' in seg:
                    print(f"  Ref location: {seg['ref_chr']}:{seg['ref_pos']}")
                    print(f"  Qry location: {seg['qry_chr']}:{seg['qry_pos']}")
                else:
                    print("  No location information stored for this segment")
            
            # Check how many segments have location information
            segments_with_locations = sum(1 for seg in self.segments if 'ref_chr' in seg and 'ref_pos' in seg)
            print(f"\nDEBUG: {segments_with_locations}/{len(self.segments)} segments have location information")
    
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

        print(f"Retained {len(self.segments)} sequences after sequence verification.")
        if not self.segments:
            print("No sequences remained after verification. Exiting.")
            sys.exit(1)
        
        # Clear validation caches to free memory
        SequenceValidator.clear_caches()
    
    def run_primer3(self):
        """
        Run Primer3 for primer design.
        
        Returns:
            pandas.DataFrame: DataFrame with primer design results
        """
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
                        all_records.extend(Primer3Processor.parse_primer3_batch(stdout_data))
            else:
                for future in concurrent.futures.as_completed(batch_futures):
                    stdout_data = future.result()
                    if stdout_data:
                        all_records.extend(Primer3Processor.parse_primer3_batch(stdout_data))

        # Create dataframe from records
        df = pd.DataFrame(all_records)
        if df.empty:
            print("No suitable primer pairs found (no pairs with penalty ≤ 5). Exiting.")
            sys.exit(0)
            
        return df
    
    def filter_primers(self, df):
        """
        Filter primers based on repeats and GC content.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer design results
            
        Returns:
            pandas.DataFrame: Filtered DataFrame
        """
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
        """
        Run BLAST for primer specificity checking.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer design results
            
        Returns:
            pandas.DataFrame: DataFrame with BLAST results
        """
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
        """Find the genomic locations of amplicons across all genomes."""
        # First, try to get locations from stored segment positions
        df["Ref Chromosome"] = None
        df["Ref Start"] = None
        df["Qry Chromosome"] = None
        df["Qry Start"] = None
        
        # Create a segment location map
        segment_locations = {
            f"{seg['gene']}_{seg['index']}": {
                "ref_chr": seg["ref_chr"],
                "ref_pos": seg["ref_pos"],
                "qry_chr": seg["qry_chr"],
                "qry_pos": seg["qry_pos"]
            }
            for seg in self.segments if "ref_chr" in seg and "ref_pos" in seg
        }
        
        if Config.DEBUG_MODE:
            print(f"\nDEBUG: Built location map with {len(segment_locations)} entries")
            print("\nDEBUG: Examples of segment location keys:")
            for key in list(segment_locations.keys())[:5]:  # Print first 5 keys
                print(f"  {key}")
            
            print("\nDEBUG: Examples of primer dataframe sequence IDs:")
            for i, row in df.head(5).iterrows():
                print(f"  {row['Sequence']}_{row['Index']}")
        
        # Try to use stored locations first
        match_count = 0
        for i, row in df.iterrows():
            seq_id = f"{row['Sequence']}_{row['Index']}"
            if seq_id in segment_locations:
                loc = segment_locations[seq_id]
                df.at[i, "Ref Chromosome"] = loc["ref_chr"]
                df.at[i, "Ref Start"] = loc["ref_pos"]
                df.at[i, "Qry Chromosome"] = loc["qry_chr"]
                df.at[i, "Qry Start"] = loc["qry_pos"]
                match_count += 1
            elif Config.DEBUG_MODE:
                print(f"DEBUG: Could not match ID: {seq_id}")
        
        if Config.DEBUG_MODE:
            print(f"\nDEBUG: Found locations for {match_count} out of {len(df)} primers")
        
        # For rows without location information, fall back to sequence search
        missing_ref = df["Ref Chromosome"].isnull()
        missing_qry = df["Qry Chromosome"].isnull()
        
        if Config.DEBUG_MODE:
            print(f"\nDEBUG: Total sequences: {len(df)}")
            print(f"DEBUG: Sequences missing Ref location: {missing_ref.sum()}")
            print(f"DEBUG: Sequences missing Qry location: {missing_qry.sum()}")
        
        # Only do sequence search if there are missing locations
        if missing_ref.any() or missing_qry.any():
            print(f"\nFalling back to sequence search for {missing_ref.sum()} reference and {missing_qry.sum()} query locations")
            amplicon_seqs = df.loc[missing_ref | missing_qry, "Amplicon"].tolist()
            batch_size = Config.BATCH_SIZE
            seq_batches = list(Utils.chunks(amplicon_seqs, batch_size))
            
            # Initialize storage for results
            ref_results, qry_results = [], []

            with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
                # Process reference genome locations
                if missing_ref.any():
                    ref_futures = [executor.submit(AmpliconLocator.find_amplicon_locations_batch, (batch, self.ref_genome))
                                for batch in seq_batches]
                    for future in concurrent.futures.as_completed(ref_futures):
                        ref_results.extend(future.result())

                # Process query genome locations
                if missing_qry.any():
                    qry_futures = [executor.submit(AmpliconLocator.find_amplicon_locations_batch, (batch, self.qry_genome))
                                for batch in seq_batches]
                    for future in concurrent.futures.as_completed(qry_futures):
                        qry_results.extend(future.result())
            
            # Update DataFrame with sequence search results
            missing_indices = df.index[missing_ref | missing_qry].tolist()
            
            if Config.DEBUG_MODE:
                found_ref = sum(1 for r in ref_results if r[0] is not None)
                found_qry = sum(1 for r in qry_results if r[0] is not None)
                print(f"DEBUG: Sequence search found {found_ref}/{len(ref_results)} ref locations")
                print(f"DEBUG: Sequence search found {found_qry}/{len(qry_results)} qry locations")
            
            for i, idx in enumerate(missing_indices):
                if i < len(ref_results) and missing_ref.iloc[missing_indices.index(idx)]:
                    df.at[idx, "Ref Chromosome"], df.at[idx, "Ref Start"] = ref_results[i]
                if i < len(qry_results) and missing_qry.iloc[missing_indices.index(idx)]:
                    df.at[idx, "Qry Chromosome"], df.at[idx, "Qry Start"] = qry_results[i]

        # Process extra genomes if provided
        if self.extra_genomes:
            extra_results = {f"Extra_{i+1}": [] for i in range(len(self.extra_genomes))}
            
            amplicon_seqs = df["Amplicon"].tolist()
            batch_size = Config.BATCH_SIZE
            seq_batches = list(Utils.chunks(amplicon_seqs, batch_size))
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
                for i, genome in enumerate(self.extra_genomes):
                    genome_name = f"Extra_{i+1}"
                    extra_futures = [executor.submit(AmpliconLocator.find_amplicon_locations_batch, (batch, genome))
                                    for batch in seq_batches]
                    for future in concurrent.futures.as_completed(extra_futures):
                        extra_results[genome_name].extend(future.result())

            # Add extra genome results dynamically
            for genome_name, results in extra_results.items():
                df[f"{genome_name} Chromosome"], df[f"{genome_name} Start"] = zip(*results)

        # Final summary with improved error reporting
        if Config.DEBUG_MODE:
            print(f"\nDEBUG: Final location statistics:")
            print(f"DEBUG: Sequences with Ref location: {df['Ref Chromosome'].notnull().sum()}/{len(df)}")
            print(f"DEBUG: Sequences with Qry location: {df['Qry Chromosome'].notnull().sum()}/{len(df)}")
            for i, genome in enumerate(self.extra_genomes):
                genome_name = f"Extra_{i+1}"
                col = f"{genome_name} Chromosome"
                print(f"DEBUG: Sequences with {genome_name} location: {df[col].notnull().sum()}/{len(df)}")
            
            # Identify sequences that couldn't be located in any genome
            missing_both = df[df['Ref Chromosome'].isnull() & df['Qry Chromosome'].isnull()]
            
            if not missing_both.empty:
                print(f"\nDEBUG: {len(missing_both)} sequences could not be located in any genome:")
                for _, row in missing_both.iterrows():
                    amplicon = row['Amplicon']
                    if pd.notnull(amplicon):
                        print(f"  Sequence: {row['Sequence']}, Index: {row['Index']}")
                        print(f"  Amplicon length: {len(amplicon)}")
                        
                        # Check for potential issues that might prevent matching
                        gc_content = Utils.calculate_gc(amplicon)
                        repeats = Utils.has_disallowed_repeats(amplicon)
                        ambiguous = 'N' in amplicon
                        
                        print(f"  GC content: {gc_content:.1f}%, Has disallowed repeats: {repeats}, Contains N: {ambiguous}")
                        
                        # Check 20bp on each end of the amplicon for features that might hinder mapping
                        start_seq = amplicon[:20] if len(amplicon) >= 20 else amplicon
                        end_seq = amplicon[-20:] if len(amplicon) >= 20 else amplicon
                        print(f"  Starts with: {start_seq}, Ends with: {end_seq}")
                        
                        # Check exact match for 30bp subsections to help diagnose partial matches
                        if len(amplicon) >= 60:
                            subsections = [amplicon[:30], amplicon[len(amplicon)//2-15:len(amplicon)//2+15], amplicon[-30:]]
                            for i, subsec in enumerate(subsections):
                                ref_match = any(subsec in seq for seq in self.ref_genome.values())
                                qry_match = any(subsec in seq for seq in self.qry_genome.values())
                                sec_name = "start" if i == 0 else "middle" if i == 1 else "end"
                                print(f"  {sec_name} 30bp: Exists in Ref: {ref_match}, Exists in Qry: {qry_match}")

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
        """
        Save final results to Excel.
        
        Args:
            df (pandas.DataFrame): DataFrame with final results
            
        Returns:
            str: Path to the output Excel file
        """
        final_cols = [
            "Sequence", "Index", "Template",
            "Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST1", "Primer F BLAST2",
            "Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST1", "Primer R BLAST2",
            "Pair Penalty", "Probe", "Probe Tm", "Probe Penalty", "Probe Reversed", "Probe dG", "Probe BLAST1", 
            "Probe BLAST2", "Amplicon", "Length", "Amplicon GC%", "Amplicon dG"
        ]
        # Add location columns in alternating order
        location_cols = []
        for label in ["Ref", "Qry"] + [f"Extra_{i+1}" for i in range(len(self.extra_genomes))]:
            location_cols.append(f"{label} Chromosome")
            location_cols.append(f"{label} Start")

        final_cols += location_cols
        
        # Reorder columns
        df = df[final_cols]

        # Create a Primers directory in the parent folder of the input files if it doesn't exist
        input_dir = os.path.dirname(self.ref_path)
        primers_dir = os.path.join(input_dir, "primers")
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
        
        print(f"\nDone! Analysis complete. Found {len(df)} primer pairs.")
        print(f"Results saved to: {output_file}")
        
        return output_file  # Return the path to the output file
    
    def run(self):
        """Run the full pipeline."""
        # 1. Get user input - this will handle Match Checker mode and exit if needed
        maf_choice = self.get_user_input()
        
        # 2. Load input files
        self.load_input_files()
        
        # 3. Run LastZ alignments if needed
        self.run_lastz(maf_choice)
        
        # 4. Run MultiZ if needed
        if self.use_multiz and len(self.maf_files) > 1:
            self.run_multiz()
        
        # If running in LastZ/MultiZ-only mode, exit here
        if Config.ONLY_LASTZ_MULTIZ:
            print("\nLastZ/MultiZ processing completed successfully.")
            if self.maf_files:
                print("Output MAF files:")
                for maf_file in self.maf_files:
                    print(f" - {maf_file}")
            sys.exit(0)
        
        # Only continue with the rest of the pipeline if not in ONLY_LASTZ_MULTIZ mode
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
        output_file = self.save_results(df)
        
        # 13. Optionally run Match Checker on output
        if Config.RUN_MATCH_CHECKER_ON_OUTPUT:
            print("\nWould you like to run Match Checker on the output against a specific genome? (y/n)")
            choice = input().strip().lower()
            if choice == 'y':
                print("Select the target genome FASTA file to check matches against:")
                target_genome = FileUtilities.get_file("Select Target Genome", 
                                                    [("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")])
                                                    
                # Run Match Checker on the output
                MatchChecker.run_from_ddprimer_output(output_file, target_genome)


##############################################################################
#                           Main Program
##############################################################################
def main():
    """Main entry point for the program."""
    parser = argparse.ArgumentParser(description='ddPrimer - Primer Design Pipeline for Multiple Genomes')
    parser.add_argument('--lastz-only', action='store_true', 
                       help='Run only LastZ/MultiZ alignments (skip primer design)')
    parser.add_argument('--match-checker', action='store_true', 
                       help='Run only the Match Checker')
    parser.add_argument('--no-match-checker-on-output', action='store_true', 
                       help='Disable Match Checker on pipeline output')
    parser.add_argument('--debug', action='store_true', 
                       help='Enable debug mode for verbose logging')
    args = parser.parse_args()
    
    # Set configuration based on command line arguments
    if args.lastz_only:
        Config.ONLY_LASTZ_MULTIZ = True
    elif args.match_checker:
        Config.RUN_MATCH_CHECKER = True
    
    if args.no_match_checker_on_output:
        Config.RUN_MATCH_CHECKER_ON_OUTPUT = False
        
    if args.debug:
        Config.DEBUG_MODE = True
    
    pipeline = PrimerDesignPipeline()
    pipeline.run()

if __name__ == "__main__":
    main()