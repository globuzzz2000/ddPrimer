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
import logging
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Import package modules
from ..config import Config, SequenceProcessingError


class SequenceProcessor:
    """Handles sequence processing and filtering operations."""
    
    # Get module logger
    logger = logging.getLogger("ddPrimer.sequence_processor")

    @staticmethod
    def cut_at_restriction_sites(sequences, restriction_site=None):
        """
        Cut sequences at restriction sites.
        
        Args:
            sequences (dict): Dictionary mapping sequence IDs to sequences
            restriction_site (str, optional): Restriction site pattern. Defaults to Config.RESTRICTION_SITE.
            
        Returns:
            list: List of fragment dictionaries
            
        Raises:
            SequenceProcessingError: If there's an issue with sequence processing
        """
        logger = logging.getLogger("ddPrimer.sequence_processor")
        logger.debug(f"Processing {len(sequences)} sequences for restriction site cutting")
        
        if restriction_site is None:
            restriction_site = Config.RESTRICTION_SITE
        
        # No cutting if no restriction site is defined
        if not restriction_site:
            logger.debug("No restriction site defined - keeping sequences intact")
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
        try:
            restriction_pattern = re.compile(restriction_site, re.IGNORECASE)
            logger.debug(f"Using restriction site pattern: {restriction_site}")
        except re.error as e:
            logger.error(f"Invalid restriction site regex pattern: {e}")
            raise SequenceProcessingError(f"Invalid restriction site pattern: {e}")
        
        # Process each sequence
        fragments = []
        
        for seq_id, sequence in sequences.items():
            try:
                # Find all restriction sites
                matches = list(restriction_pattern.finditer(sequence))
                
                if not matches:
                    # No restriction sites - keep the entire sequence
                    logger.debug(f"No restriction sites found in {seq_id} - keeping entire sequence ({len(sequence)} bp)")
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
                logger.debug(f"Found {len(matches)} restriction sites in {seq_id}")
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
                        logger.debug(f"Created fragment {fragment['id']}: {end-start} bp")
                    else:
                        logger.debug(f"Skipped fragment {i} in {seq_id} - too short ({end-start} bp)")
                    
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
                    logger.debug(f"Created final fragment {fragment['id']}: {len(sequence)-last_end} bp")
                else:
                    logger.debug(f"Skipped final fragment in {seq_id} - too short ({len(sequence)-last_end} bp)")
                    
            except Exception as e:
                logger.error(f"Error processing sequence {seq_id}: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                continue
        
        logger.debug(f"Created {len(fragments)} fragments from {len(sequences)} sequences")
        return fragments
    
    @staticmethod
    def filter_by_gene_overlap(fragments, genes, margin=None):
        """
        Filter fragments by gene overlap.
        
        Args:
            fragments (list): List of fragment dictionaries
            genes (list): List of gene dictionaries from GFF
            margin (int, optional): Gene overlap margin in base pairs. Defaults to Config.GENE_OVERLAP_MARGIN.
            
        Returns:
            list: List of filtered fragment dictionaries
            
        Raises:
            SequenceProcessingError: If there's an issue with sequence processing
        """
        logger = logging.getLogger("ddPrimer.sequence_processor")
        logger.debug(f"Filtering {len(fragments)} fragments by gene overlap")
        
        if margin is None:
            margin = Config.GENE_OVERLAP_MARGIN
        
        logger.debug(f"Using gene overlap margin: {margin} bp")
        
        if not genes:
            logger.warning("No genes provided for filtering - returning original fragments")
            return fragments
            
        filtered_fragments = []
        
        try:
            # Create a mapping of chromosomes to genes for faster lookup
            gene_map = {}
            for gene in genes:
                chr_name = gene["chr"]
                if chr_name not in gene_map:
                    gene_map[chr_name] = []
                gene_map[chr_name].append(gene)
            
            logger.debug(f"Created gene map with {len(gene_map)} chromosomes")
            
            # Process each fragment
            for fragment in fragments:
                # Check if fragment overlaps with any gene
                fragment_chr = fragment["chr"]
                overlapping_gene = None
                
                # Only check genes on the same chromosome
                chr_genes = gene_map.get(fragment_chr, [])
                
                for gene in chr_genes:
                    # Check for overlap
                    if (fragment["start"] <= gene["end"] and 
                        fragment["end"] >= gene["start"]):
                        
                        overlapping_gene = gene
                        logger.debug(f"Fragment {fragment['id']} overlaps with gene {gene['id']}")
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
                        logger.debug(f"Added filtered fragment {filtered_fragment['id']}: {seq_end-seq_start} bp")
                    else:
                        logger.debug(f"Filtered fragment too short after truncation: {new_end-new_start+1} bp")
                else:
                    logger.debug(f"Fragment {fragment['id']} does not overlap with any gene")
        
        except Exception as e:
            logger.error(f"Error during gene overlap filtering: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(f"Gene overlap filtering failed: {e}")
        
        logger.debug(f"Filtered to {len(filtered_fragments)} fragments that overlap with genes")
        return filtered_fragments