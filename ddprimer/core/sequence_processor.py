#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence processing module for ddPrimer pipeline.

Contains functionality for:
1. Restriction site cutting with configurable patterns
2. Gene overlap filtering and fragment generation
3. Sequence validation and coordinate handling
4. Fragment length filtering and optimization
"""

import re
import logging

from ..config import Config, SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)


class SequenceProcessor:
    """
    Handles sequence processing and filtering operations.
    
    This class provides comprehensive sequence manipulation capabilities
    including restriction site analysis, fragment generation, and
    coordinate-based filtering for primer design workflows.
    
    Features:
        - Configurable restriction site cutting
        - Fragment length validation
        - Coordinate system handling (1-based genomic coordinates)
        - Gene overlap detection and filtering
        
    Example:
        >>> fragments = SequenceProcessor.cut_at_restriction_sites(sequences)
        >>> filtered = SequenceProcessor.filter_by_gene_overlap(fragments, genes)
    """

    @staticmethod
    def cut_at_restriction_sites(sequences, restriction_site=None):
        """
        Cut sequences at restriction sites to generate fragments.
        
        Processes DNA sequences by identifying restriction enzyme cut sites
        and generating fragments suitable for primer design, with proper
        coordinate tracking and length validation.
        
        Args:
            sequences: Dictionary mapping sequence IDs to DNA sequences
            restriction_site: Restriction site pattern, defaults to Config.RESTRICTION_SITE
            
        Returns:
            List of fragment dictionaries with coordinates and sequences
            
        Raises:
            SequenceProcessingError: If restriction site pattern is invalid or processing fails
            
        Example:
            >>> sequences = {"chr1": "GAATTCATCGAATTCGCTA"}
            >>> fragments = SequenceProcessor.cut_at_restriction_sites(sequences, "GAATTC")
            >>> print(f"Generated {len(fragments)} fragments")
        """
        logger.debug(f"=== RESTRICTION SITE CUTTING ===")
        logger.debug(f"Processing {len(sequences)} sequences for restriction site cutting")
        
        if restriction_site is None:
            restriction_site = Config.RESTRICTION_SITE
        
        # Handle case where no restriction site is defined
        if not restriction_site:
            logger.debug("No restriction site pattern defined - keeping sequences intact")
            fragments = []
            
            for seq_id, sequence in sequences.items():
                fragment = {
                    "id": seq_id,
                    "chr": seq_id,
                    "start": 1,  # 1-based genomic coordinates
                    "end": len(sequence),
                    "sequence": sequence
                }
                fragments.append(fragment)
                logger.debug(f"Created intact fragment: {seq_id} ({len(sequence)} bp)")
            
            logger.debug(f"Created {len(fragments)} intact fragments")
            return fragments
            
        # Compile and validate restriction site pattern
        try:
            restriction_pattern = re.compile(restriction_site, re.IGNORECASE)
            logger.debug(f"Using restriction site pattern: {restriction_site}")
        except re.error as e:
            error_msg = f"Invalid restriction site regex pattern: {restriction_site}"
            logger.error(error_msg)
            logger.debug(f"Regex compilation error: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
        
        # Process each sequence for restriction site cutting
        fragments = []
        total_sites_found = 0
        
        for seq_id, sequence in sequences.items():
            try:
                logger.debug(f"Processing sequence {seq_id} ({len(sequence)} bp)")
                
                # Find all restriction sites in sequence
                matches = list(restriction_pattern.finditer(sequence))
                sites_in_sequence = len(matches)
                total_sites_found += sites_in_sequence
                
                if not matches:
                    # No restriction sites found - keep entire sequence
                    logger.debug(f"No restriction sites found in {seq_id} - keeping entire sequence")
                    fragment = {
                        "id": seq_id,
                        "chr": seq_id,
                        "start": 1,
                        "end": len(sequence),
                        "sequence": sequence
                    }
                    fragments.append(fragment)
                    continue
                
                logger.debug(f"Found {sites_in_sequence} restriction sites in {seq_id}")
                
                # Generate fragments between restriction sites
                fragment_count = 0
                last_end = 0
                
                for i, match in enumerate(matches):
                    start_pos = last_end
                    end_pos = match.start()
                    fragment_length = end_pos - start_pos
                    
                    # Check minimum fragment length requirement
                    if fragment_length >= Config.MIN_SEGMENT_LENGTH:
                        fragment = {
                            "id": f"{seq_id}_frag{fragment_count}",
                            "chr": seq_id,
                            "start": start_pos + 1,  # Convert to 1-based coordinates
                            "end": end_pos,
                            "sequence": sequence[start_pos:end_pos]
                        }
                        fragments.append(fragment)
                        fragment_count += 1
                        
                        logger.debug(f"Created fragment {fragment['id']}: {fragment_length} bp "
                                   f"({fragment['start']}-{fragment['end']})")
                    else:
                        logger.debug(f"Skipped fragment {i} in {seq_id}: {fragment_length} bp < {Config.MIN_SEGMENT_LENGTH} bp minimum")
                    
                    last_end = match.end()
                
                # Handle final fragment after last restriction site
                final_start = last_end
                final_length = len(sequence) - final_start
                
                if final_length >= Config.MIN_SEGMENT_LENGTH:
                    fragment = {
                        "id": f"{seq_id}_frag{fragment_count}",
                        "chr": seq_id,
                        "start": final_start + 1,  # Convert to 1-based coordinates
                        "end": len(sequence),
                        "sequence": sequence[final_start:]
                    }
                    fragments.append(fragment)
                    fragment_count += 1
                    
                    logger.debug(f"Created final fragment {fragment['id']}: {final_length} bp "
                               f"({fragment['start']}-{fragment['end']})")
                else:
                    logger.debug(f"Skipped final fragment in {seq_id}: {final_length} bp < {Config.MIN_SEGMENT_LENGTH} bp minimum")
                
                logger.debug(f"Generated {fragment_count} valid fragments from {seq_id}")
                    
            except Exception as e:
                error_msg = f"Error processing sequence {seq_id} for restriction sites"
                logger.error(error_msg)
                logger.debug(f"Processing error: {str(e)}", exc_info=True)
                # Continue with other sequences rather than failing completely
                continue
        
        # Log comprehensive cutting results
        logger.debug(f"=== RESTRICTION SITE CUTTING RESULTS ===")
        logger.debug(f"Input sequences: {len(sequences)}")
        logger.debug(f"Total restriction sites found: {total_sites_found}")
        logger.debug(f"Output fragments: {len(fragments)}")
        logger.debug(f"Minimum fragment length: {Config.MIN_SEGMENT_LENGTH} bp")
        
        if logger.isEnabledFor(logging.DEBUG):
            # Fragment length statistics
            fragment_lengths = [len(frag['sequence']) for frag in fragments]
            if fragment_lengths:
                avg_length = sum(fragment_lengths) / len(fragment_lengths)
                min_length = min(fragment_lengths)
                max_length = max(fragment_lengths)
                logger.debug(f"Fragment length statistics: avg={avg_length:.0f}, min={min_length}, max={max_length}")
        
        logger.info(f"Restriction site cutting completed: {len(fragments)} fragments generated")
        logger.debug(f"=== END RESTRICTION SITE CUTTING ===")
        
        return fragments
    
    @staticmethod
    def filter_by_gene_overlap(fragments, genes, margin=None):
        """
        Filter fragments by gene overlap with coordinate-based truncation.
        
        Identifies fragments that overlap with annotated genes and truncates
        them to gene boundaries with optional margin extension for improved
        primer design targeting.
        
        Args:
            fragments: List of fragment dictionaries with coordinates
            genes: List of gene dictionaries from GFF annotations
            margin: Gene overlap margin in base pairs, defaults to Config.GENE_OVERLAP_MARGIN
            
        Returns:
            List of filtered fragment dictionaries truncated to gene regions
            
        Raises:
            SequenceProcessingError: If gene overlap processing fails
            
        Example:
            >>> filtered = SequenceProcessor.filter_by_gene_overlap(fragments, genes, margin=100)
            >>> print(f"Filtered to {len(filtered)} gene-overlapping fragments")
        """
        logger.debug(f"=== GENE OVERLAP FILTERING ===")
        logger.debug(f"Filtering {len(fragments)} fragments by gene overlap")
        
        if margin is None:
            margin = Config.GENE_OVERLAP_MARGIN
        
        logger.debug(f"Using gene overlap margin: {margin} bp")
        
        if not genes:
            logger.warning("No gene annotations provided for filtering - returning original fragments")
            return fragments
            
        filtered_fragments = []
        
        try:
            # Create chromosome-based gene mapping for efficient lookup
            gene_map = {}
            for gene in genes:
                chr_name = gene["chr"]
                if chr_name not in gene_map:
                    gene_map[chr_name] = []
                gene_map[chr_name].append(gene)
            
            logger.debug(f"Created gene map with {len(gene_map)} chromosomes")
            
            # Log gene statistics
            total_genes = sum(len(chr_genes) for chr_genes in gene_map.values())
            logger.debug(f"Total genes available for overlap detection: {total_genes}")
            
            # Process each fragment for gene overlap
            fragments_with_overlaps = 0
            
            for fragment in fragments:
                fragment_chr = fragment["chr"]
                fragment_start = fragment["start"]
                fragment_end = fragment["end"]
                overlapping_gene = None
                
                logger.debug(f"Checking fragment {fragment['id']} ({fragment_chr}:{fragment_start}-{fragment_end})")
                
                # Get genes on the same chromosome
                chr_genes = gene_map.get(fragment_chr, [])
                
                if not chr_genes:
                    logger.debug(f"No genes found on chromosome {fragment_chr}")
                    continue
                
                # Find first overlapping gene
                for gene in chr_genes:
                    gene_start = gene["start"]
                    gene_end = gene["end"]
                    
                    # Check for overlap using standard interval logic
                    if fragment_start <= gene_end and fragment_end >= gene_start:
                        overlapping_gene = gene
                        logger.debug(f"Fragment {fragment['id']} overlaps with gene {gene['id']} "
                                   f"({gene_start}-{gene_end})")
                        break
                
                if overlapping_gene:
                    fragments_with_overlaps += 1
                    
                    # Calculate truncated coordinates with margin
                    new_start = max(fragment_start, overlapping_gene["start"] - margin)
                    new_end = min(fragment_end, overlapping_gene["end"] + margin)
                    
                    # Ensure coordinates remain within original fragment bounds
                    new_start = max(new_start, fragment_start)
                    new_end = min(new_end, fragment_end)
                    
                    # Check if truncated fragment meets minimum length requirement
                    truncated_length = new_end - new_start + 1
                    
                    if truncated_length >= Config.MIN_SEGMENT_LENGTH:
                        # Calculate sequence indices for truncation
                        seq_start = new_start - fragment_start
                        seq_end = new_end - fragment_start + 1
                        
                        # Validate sequence boundaries
                        seq_start = max(0, seq_start)
                        seq_end = min(len(fragment["sequence"]), seq_end)
                        
                        # Create filtered fragment
                        filtered_fragment = {
                            "id": f"{fragment['id']}_{overlapping_gene['id']}",
                            "chr": fragment["chr"],
                            "start": new_start,
                            "end": new_end,
                            "gene": overlapping_gene["id"],
                            "sequence": fragment["sequence"][seq_start:seq_end]
                        }
                        
                        filtered_fragments.append(filtered_fragment)
                        
                        logger.debug(f"Added filtered fragment {filtered_fragment['id']}: "
                                   f"{len(filtered_fragment['sequence'])} bp "
                                   f"({new_start}-{new_end})")
                    else:
                        logger.debug(f"Truncated fragment too short: {truncated_length} bp < {Config.MIN_SEGMENT_LENGTH} bp")
                else:
                    logger.debug(f"Fragment {fragment['id']} does not overlap with any genes")
        
        except Exception as e:
            error_msg = f"Gene overlap filtering failed"
            logger.error(error_msg)
            logger.debug(f"Filtering error: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
        
        # Log filtering results
        logger.debug(f"=== GENE OVERLAP FILTERING RESULTS ===")
        logger.debug(f"Input fragments: {len(fragments)}")
        logger.debug(f"Fragments with gene overlaps: {fragments_with_overlaps}")
        logger.debug(f"Output filtered fragments: {len(filtered_fragments)}")
        
        if logger.isEnabledFor(logging.DEBUG):
            # Analyze filtered fragment statistics
            if filtered_fragments:
                filtered_lengths = [len(frag['sequence']) for frag in filtered_fragments]
                avg_filtered_length = sum(filtered_lengths) / len(filtered_lengths)
                logger.debug(f"Filtered fragment length statistics: avg={avg_filtered_length:.0f}")
                
                # Count unique genes
                unique_genes = set(frag.get('gene', 'unknown') for frag in filtered_fragments)
                logger.debug(f"Unique genes represented: {len(unique_genes)}")
        
        logger.info(f"Gene overlap filtering completed: {len(filtered_fragments)} gene-overlapping fragments")
        logger.debug(f"=== END GENE OVERLAP FILTERING ===")
        
        return filtered_fragments