#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotation processing module for ddPrimer pipeline.

Handles GFF annotation file parsing and processing to extract gene information
from prepared GFF files for use in primer design workflows.
Contains functionality for:
1. GFF file parsing with parallel processing support
2. Gene overlap detection and fragment creation
3. Gene overlap extraction for comprehensive coverage
4. Integration with prepared/indexed GFF files

This module integrates with the broader ddPrimer pipeline to provide
robust gene annotation capabilities for primer design workflows. Works
with GFF files prepared by FilePreparator (sorted, compressed, indexed).
"""

import concurrent.futures
import gzip
import logging
import os
from tqdm import tqdm
from typing import List, Dict, Optional

# Import package modules
from ..config import Config, FileError, SequenceProcessingError

# Set up module logger
logger = logging.getLogger(__name__)


class AnnotationProcessor:
    """
    Processes GFF annotation files to extract gene information.
    
    This class provides methods for parsing prepared GFF files, extracting gene overlaps,
    and creating gene-specific fragments for primer design workflows. Works with both
    prepared (sorted/indexed) and standard GFF files.
    
    Example:
        >>> processor = AnnotationProcessor()
        >>> genes = processor.load_genes_from_gff("annotations.gff")
        >>> gene_fragments = processor.extract_all_gene_overlaps(fragments, genes)
    """
    
    @staticmethod
    def parse_gff_attributes(attribute_str: str) -> Dict[str, str]:
        """
        Convert GFF attribute string (key1=val1;key2=val2) -> dict.
        
        Args:
            attribute_str: GFF attribute string in format key1=val1;key2=val2
            
        Returns:
            Dictionary of attribute key-value pairs with keys forced to lowercase
        """
        attr_dict = {}
        for attr in attribute_str.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key.strip().lower()] = value.strip()
        return attr_dict
    
    @classmethod
    def process_gff_chunk(cls, chunk: List[str]) -> List[Dict]:
        """
        Process a chunk of GFF file lines.
        Used for parallel processing.
        
        Args:
            chunk: List of GFF file lines
            
        Returns:
            List of gene dictionaries with extracted information
        """
        chunk_genes = []
        
        # Parse RETAIN_TYPES - handle both string and list formats
        if isinstance(Config.RETAIN_TYPES, str):
            # Split by comma and clean whitespace
            retain_types = [t.strip().lower() for t in Config.RETAIN_TYPES.split(',')]
        elif isinstance(Config.RETAIN_TYPES, list):
            retain_types = [t.lower() for t in Config.RETAIN_TYPES]
        else:
            retain_types = [str(Config.RETAIN_TYPES).lower()]
        
        for line in chunk:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqname, _, ftype, start, end, _, strand, _, attributes = parts
            
            # Check if feature type is in our retain list
            if ftype.lower() not in retain_types:
                continue

            attr_dict = cls.parse_gff_attributes(attributes)

            # Get the best available identifier
            name = attr_dict.get('name')
            gene_id = attr_dict.get('id')
            locus_tag = attr_dict.get('locus_tag')
            
            # Use the first available identifier (prefer name, then id, then locus_tag)
            identifier = name or gene_id or locus_tag or f"feature_{start}_{end}"

            try:
                chunk_genes.append({
                    "chr": seqname,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "id": identifier
                })
            except ValueError as e:
                logger.debug(f"Error converting position data: {e}")
        
        return chunk_genes
    
    @classmethod
    def load_genes_from_gff(cls, gff_path: str) -> List[Dict]:
        """
        Extract gene information from a GFF file.
        
        Handles both compressed (.gz) and uncompressed GFF files.
        Uses parallel processing for large files.
        
        Args:
            gff_path: Path to the GFF file (prepared or standard)
            
        Returns:
            List of gene dictionaries with extracted information
            
        Raises:
            FileError: If GFF file cannot be read
            SequenceProcessingError: If GFF processing fails
        """
        if not os.path.exists(gff_path):
            error_msg = f"GFF file not found: {gff_path}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        try:
            # Determine if file is compressed
            is_compressed = gff_path.endswith('.gz')
            opener = gzip.open if is_compressed else open
            mode = 'rt' if is_compressed else 'r'
            
            # Read all lines from the file
            with opener(gff_path, mode) as f:
                all_lines = f.readlines()
            
            logger.debug(f"Read {len(all_lines)} lines from GFF file")
            
            # Calculate chunk size for parallel processing
            chunk_size = max(1, len(all_lines) // Config.NUM_PROCESSES)
            chunks = [all_lines[i:i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
            
            logger.debug(f"Split into {len(chunks)} chunks for parallel processing")
            
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
            
            logger.info(f"Extracted {len(genes)} genes from GFF file")
            
            # Log sample of genes if debug enabled (only first few)
            if logger.isEnabledFor(logging.DEBUG) and genes:
                logger.debug("Sample genes extracted:")
                for i, gene in enumerate(genes[:5]):  # Only show first 5
                    logger.debug(f"  {i+1}. {gene['id']} at {gene['chr']}:{gene['start']}-{gene['end']}")
                if len(genes) > 5:
                    logger.debug(f"  ... and {len(genes) - 5} more genes")
            
            return genes
        
        except Exception as e:
            error_msg = f"Failed to load genes from GFF: {str(e)}"
            logger.error(error_msg)
            raise SequenceProcessingError(error_msg) from e

    @staticmethod
    def extract_gene_name(sequence_id: str) -> str:
        """
        Extract just the gene name from a sequence identifier.
        
        Handles format Chr_Fragment_Gene, returning the gene name component.
        
        Args:
            sequence_id: The full sequence identifier
            
        Returns:
            The extracted gene name
        """
        # Split by underscore
        parts = sequence_id.split("_")
        
        # For format Chr_Fragment_Gene, return the third part if it exists
        if len(parts) >= 3:
            return parts[2]  # Return the Gene part
        
        # If the ID doesn't have enough parts, return the original
        return sequence_id

    @classmethod
    def filter_by_gene_overlap(cls, restriction_fragments: List[Dict], 
                              genes: List[Dict]) -> List[Dict]:
        """
        Filter restriction fragments by gene overlap and extract gene regions.
        
        This is the main function for extracting gene-overlapping regions from
        restriction fragments. Creates separate fragments for each gene overlap.
        
        Args:
            restriction_fragments: List of restriction fragments
            genes: List of gene annotations
            
        Returns:
            List of all gene-overlapping fragments
        """
        overlap_margin = getattr(Config, 'GENE_OVERLAP_MARGIN', 0)
        
        # Extract all gene overlaps
        gene_fragments = cls.extract_all_gene_overlaps(
            restriction_fragments, genes, overlap_margin
        )
        
        logger.info(f"Extracted {len(gene_fragments)} gene fragments from {len(restriction_fragments)} restriction fragments")
        return gene_fragments

    @classmethod
    def extract_all_gene_overlaps(cls, restriction_fragments: List[Dict], genes: List[Dict], 
                                 overlap_margin: int = 0) -> List[Dict]:
        """
        Extract ALL gene-overlapping regions from each restriction fragment.
        
        This function processes each restriction fragment and identifies all genes
        that overlap with it, creating separate fragments for each gene overlap.
        
        Args:
            restriction_fragments: List of restriction fragments with chr, start, end, sequence
            genes: List of gene annotations with chr, start, end, strand, id
            overlap_margin: Additional margin around genes (bp)
            
        Returns:
            List of gene-overlapping fragments
        """
        gene_fragments = []
        processed_fragments = 0
        fragments_with_overlaps = 0
        
        for fragment in restriction_fragments:
            fragment_chr = fragment.get("chr")
            fragment_start = fragment.get("start")
            fragment_end = fragment.get("end")
            fragment_seq = fragment.get("sequence", "")
            
            if not all([fragment_chr, fragment_start is not None, fragment_end is not None]):
                # Only log missing data occasionally
                if processed_fragments % 100 == 0:
                    logger.warning(f"Fragment {fragment.get('id', 'unknown')} missing position data, skipping")
                processed_fragments += 1
                continue
            
            # Find all genes that overlap with this fragment
            overlapping_genes = cls.find_overlapping_genes(
                fragment_chr, fragment_start, fragment_end, genes, overlap_margin
            )
            
            if overlapping_genes:
                fragments_with_overlaps += 1
                # Create a separate fragment for each overlapping gene
                for gene_idx, gene in enumerate(overlapping_genes):
                    gene_fragment = cls.create_gene_fragment(
                        fragment, gene, gene_idx, overlap_margin
                    )
                    
                    if gene_fragment:
                        gene_fragments.append(gene_fragment)
                
                # Sample logging - only log every 100th fragment OR those with many overlaps
                if (logger.isEnabledFor(logging.DEBUG) and 
                    (processed_fragments % 100 == 0 or len(overlapping_genes) > 3)):
                    logger.debug(f"Fragment {fragment['id']} ({fragment_chr}:{fragment_start}-{fragment_end}) "
                               f"overlaps with {len(overlapping_genes)} genes")
            else:
                # Only log fragments without overlaps occasionally
                if logger.isEnabledFor(logging.DEBUG) and processed_fragments % 100 == 0:
                    logger.debug(f"Fragment {fragment['id']} ({fragment_chr}:{fragment_start}-{fragment_end}) "
                               f"has no gene overlaps")
            
            processed_fragments += 1
        
        logger.debug(f"Gene overlap extraction complete: {processed_fragments} fragments processed, "
                   f"{fragments_with_overlaps} had overlaps, generated {len(gene_fragments)} gene fragments")
        
        return gene_fragments

    @classmethod
    def find_overlapping_genes(cls, fragment_chr: str, fragment_start: int, fragment_end: int,
                              genes: List[Dict], overlap_margin: int = 0) -> List[Dict]:
        """
        Find all genes that overlap with a given genomic region.
        
        Args:
            fragment_chr: Chromosome/sequence name of the fragment
            fragment_start: Start position of the fragment
            fragment_end: End position of the fragment
            genes: List of gene annotations
            overlap_margin: Additional margin around genes
            
        Returns:
            List of overlapping genes
        """
        overlapping_genes = []
        genes_checked = 0
        
        for gene in genes:
            # Check if gene is on the same chromosome/sequence
            if gene.get("chr") != fragment_chr:
                continue
            
            gene_start = gene.get("start")
            gene_end = gene.get("end")
            
            if gene_start is None or gene_end is None:
                continue
            
            # Apply overlap margin to gene coordinates
            gene_start_with_margin = max(1, gene_start - overlap_margin)
            gene_end_with_margin = gene_end + overlap_margin
            
            # Check for overlap: genes overlap if they don't NOT overlap
            # No overlap if: fragment_end < gene_start OR fragment_start > gene_end
            if not (fragment_end < gene_start_with_margin or fragment_start > gene_end_with_margin):
                overlapping_genes.append(gene)
                
                # Only log details for first few overlaps or every 100th gene checked
                if (logger.isEnabledFor(logging.DEBUG) and 
                    (len(overlapping_genes) <= 3 or genes_checked % 100 == 0)):
                    logger.debug(f"Gene {gene['id']} ({gene_start}-{gene_end}) overlaps with fragment "
                               f"({fragment_start}-{fragment_end})")
            
            genes_checked += 1
        
        return overlapping_genes

    @classmethod
    def create_gene_fragment(cls, fragment: Dict, gene: Dict, gene_idx: int, 
                            overlap_margin: int = 0) -> Optional[Dict]:
        """
        Create a new fragment for a specific gene overlap region.
        
        Args:
            fragment: Original restriction fragment
            gene: Gene annotation that overlaps with the fragment
            gene_idx: Index of this gene (for unique naming)
            overlap_margin: Overlap margin used
            
        Returns:
            New fragment dictionary for the gene region, or None if invalid
        """
        fragment_chr = fragment.get("chr")
        fragment_start = fragment.get("start")
        fragment_end = fragment.get("end")
        fragment_seq = fragment.get("sequence", "")
        
        gene_start = gene.get("start")
        gene_end = gene.get("end")
        gene_id = gene.get("id", "unknown_gene")
        
        if not all([fragment_start is not None, fragment_end is not None, 
                    gene_start is not None, gene_end is not None]):
            return None
        
        # Calculate the overlap region with margin
        gene_start_with_margin = max(1, gene_start - overlap_margin)
        gene_end_with_margin = gene_end + overlap_margin
        
        # Find the actual overlap region within the fragment
        overlap_start = max(fragment_start, gene_start_with_margin)
        overlap_end = min(fragment_end, gene_end_with_margin)
        
        if overlap_start >= overlap_end:
            return None
        
        # Calculate positions within the fragment sequence
        seq_start = overlap_start - fragment_start
        seq_end = overlap_end - fragment_start
        
        # Validate sequence boundaries
        if seq_start < 0 or seq_end > len(fragment_seq) or seq_start >= seq_end:
            # Only log validation warnings occasionally
            if gene_idx % 100 == 0:
                logger.warning(f"Invalid sequence boundaries for gene {gene_id}: "
                              f"seq_start={seq_start}, seq_end={seq_end}, fragment_len={len(fragment_seq)}")
            return None
        
        # Extract the overlapping sequence
        gene_sequence = fragment_seq[seq_start:seq_end]
        
        # Use Config minimum segment length
        if len(gene_sequence) < Config.MIN_SEGMENT_LENGTH:
            # Only log length warnings occasionally
            if gene_idx % 100 == 0:
                logger.debug(f"Gene fragment too short for {gene_id}: {len(gene_sequence)} < {Config.MIN_SEGMENT_LENGTH}")
            return None
        
        # Create new fragment ID
        original_id = fragment.get("id", "unknown")
        new_id = f"{original_id}_gene{gene_idx}_{gene_id}"
        
        # Create the new fragment
        gene_fragment = {
            "id": new_id,
            "sequence": gene_sequence,
            "chr": fragment_chr,
            "start": overlap_start,
            "end": overlap_end,
            "Gene": gene_id,
            "original_fragment": original_id,
            "gene_start": gene_start,
            "gene_end": gene_end,
            "overlap_margin": overlap_margin
        }
        
        # Only log successful fragment creation occasionally
        if logger.isEnabledFor(logging.DEBUG) and gene_idx % 100 == 0:
            logger.debug(f"Created gene fragment {new_id}: {len(gene_sequence)} bp "
                       f"({overlap_start}-{overlap_end})")
        
        return gene_fragment