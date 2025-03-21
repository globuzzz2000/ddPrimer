#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence processing module for ddPrimer pipeline.

Contains functionality for:
1. Restriction site cutting
2. Gene overlap filtering
3. Sequence matching against genomes
"""

import re
import os
import sys
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from ..config import Config
from ..utils.sequence_utils import SequenceUtils


class SequenceProcessor:
    """Handles sequence processing and filtering operations."""
    
    @staticmethod
    def cut_at_restriction_sites(sequences, restriction_site=None):
        """
        Cut sequences at restriction sites.
        
        Args:
            sequences (dict): Dictionary mapping sequence IDs to sequences
            restriction_site (str): Restriction site pattern (default: from Config)
            
        Returns:
            list: List of fragment dictionaries
        """
        if restriction_site is None:
            restriction_site = Config.RESTRICTION_SITE
        
        # No cutting if no restriction site is defined
        if not restriction_site:
            fragments = []
            for seq_id, sequence in sequences.items():
                fragment = {
                    "id": seq_id,
                    "chr": seq_id,
                    "start": 1,  # 1-based coordinates
                    "end": len(sequence),
                    "sequence": sequence
                }
                fragments.append(fragment)
            return fragments
            
        # Compile regex for restriction site
        restriction_pattern = re.compile(restriction_site, re.IGNORECASE)
        
        # Process each sequence
        fragments = []
        
        for seq_id, sequence in sequences.items():
            # Find all restriction sites
            matches = list(restriction_pattern.finditer(sequence))
            
            if not matches:
                # No restriction sites - keep the entire sequence
                fragment = {
                    "id": seq_id,
                    "chr": seq_id,
                    "start": 1,  # 1-based coordinates
                    "end": len(sequence),
                    "sequence": sequence
                }
                fragments.append(fragment)
                continue
            
            # Create fragments between restriction sites
            last_end = 0
            for i, match in enumerate(matches):
                start = last_end
                end = match.start()
                
                if end - start >= Config.MIN_SEGMENT_LENGTH:
                    fragment = {
                        "id": f"{seq_id}_frag{i}",
                        "chr": seq_id,
                        "start": start + 1,  # Convert to 1-based
                        "end": end,
                        "sequence": sequence[start:end]
                    }
                    fragments.append(fragment)
                
                last_end = match.end()
            
            # Add final fragment after last restriction site
            if len(sequence) - last_end >= Config.MIN_SEGMENT_LENGTH:
                fragment = {
                    "id": f"{seq_id}_frag{len(matches)}",
                    "chr": seq_id,
                    "start": last_end + 1,  # Convert to 1-based
                    "end": len(sequence),
                    "sequence": sequence[last_end:]
                }
                fragments.append(fragment)
        
        return fragments
    
    @staticmethod
    def filter_by_gene_overlap(fragments, genes, margin=None):
        """
        Filter fragments by gene overlap.
        
        Args:
            fragments (list): List of fragment dictionaries
            genes (list): List of gene dictionaries from GFF
            margin (int): Gene overlap margin in base pairs (default: from Config)
            
        Returns:
            list: List of filtered fragment dictionaries
        """
        if margin is None:
            margin = Config.GENE_OVERLAP_MARGIN
        
        filtered_fragments = []
        
        for fragment in fragments:
            # Check if fragment overlaps with any gene
            overlapping_gene = None
            
            for gene in genes:
                if gene["chr"] == fragment["chr"]:
                    # Check for overlap
                    if (fragment["start"] <= gene["end"] and 
                        fragment["end"] >= gene["start"]):
                        
                        overlapping_gene = gene
                        break
            
            if overlapping_gene:
                # Truncate fragment if it extends beyond gene boundaries
                new_start = max(fragment["start"], 
                               overlapping_gene["start"] - margin)
                new_end = min(fragment["end"], 
                             overlapping_gene["end"] + margin)
                
                # Check if truncated fragment is still valid
                if new_end - new_start + 1 >= Config.MIN_SEGMENT_LENGTH:
                    # Adjust sequence to match truncated coordinates
                    seq_start = new_start - fragment["start"]
                    seq_end = new_end - fragment["start"] + 1
                    
                    # Ensure valid sequence indices
                    seq_start = max(0, seq_start)
                    seq_end = min(len(fragment["sequence"]), seq_end)
                    
                    # Create updated fragment
                    filtered_fragment = {
                        "id": f"{fragment['id']}_{overlapping_gene['id']}",
                        "chr": fragment["chr"],
                        "start": new_start,
                        "end": new_end,
                        "gene": overlapping_gene["id"],
                        "sequence": fragment["sequence"][seq_start:seq_end]
                    }
                    
                    filtered_fragments.append(filtered_fragment)
        
        return filtered_fragments
