#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotation processing module for ddPrimer pipeline.

Handles GFF annotation file parsing and processing.
"""

import re
import concurrent.futures
from tqdm import tqdm
from ..config import Config
from ..utils.common_utils import CommonUtils


class AnnotationProcessor:
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
        chunks = CommonUtils.chunks(all_lines, chunk_size)
        
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
