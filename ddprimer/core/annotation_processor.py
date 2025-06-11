#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotation processing module for ddPrimer pipeline.

Handles GFF annotation file parsing and processing to extract gene information
from GFF files for use in primer design workflows.
Contains functionality for:
1. GFF file parsing with parallel processing support
2. Gene overlap detection and fragment creation
3. Meaningful gene name filtering with species-specific patterns
4. Enhanced gene overlap extraction for comprehensive coverage

This module integrates with the broader ddPrimer pipeline to provide
robust gene annotation capabilities for primer design workflows.
"""

import re
import concurrent.futures
import logging
from tqdm import tqdm
from typing import List, Dict

# Import package modules
from ..config import Config
from ..utils import CommonUtils

# Set up module logger
logger = logging.getLogger(__name__)


class AnnotationProcessor:
    """
    Processes GFF annotation files to extract gene information.
    
    This class provides methods for parsing GFF files, extracting gene overlaps,
    and creating gene-specific fragments for primer design workflows.
    
    Attributes:
        PLACEHOLDER_PATTERNS: Regex patterns for detecting placeholder gene names
        
    Example:
        >>> processor = AnnotationProcessor()
        >>> genes = processor.load_genes_from_gff("annotations.gff")
        >>> gene_fragments = processor.extract_all_gene_overlaps(fragments, genes)
    """
    
    # Patterns for general placeholder detection across species
    PLACEHOLDER_PATTERNS = [
        re.compile(r'^AT\dG\d{5}$', re.IGNORECASE),   # Arabidopsis
        re.compile(r'^LOC_.*$', re.IGNORECASE),       # NCBI-style, rice, etc.
        re.compile(r'^GRMZM.*$', re.IGNORECASE),      # Maize
        re.compile(r'^ENS.*$', re.IGNORECASE)         # Ensembl IDs (humans, etc.)
    ]
    
    @staticmethod
    def parse_gff_attributes(attribute_str: str) -> dict:
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
    def is_meaningful_name(cls, name: str, gene_id: str = None, locus_tag: str = None) -> bool:
        """
        Determine whether the 'name' is a meaningful gene symbol.
        
        Excludes placeholders, technical locus tags, or accessions.
        
        Args:
            name: Gene name to check
            gene_id: Gene ID for comparison
            locus_tag: Locus tag for comparison
            
        Returns:
            True if name is meaningful, False otherwise
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
            except ValueError as e:
                logger.debug(f"Error converting position data: {e}")
        
        return chunk_genes
    
    @classmethod
    def load_genes_from_gff(cls, gff_path: str) -> List[Dict]:
        """
        Extract gene information from a GFF file.
        
        Uses global RETAIN_TYPES and FILTER_MEANINGFUL_NAMES for filtering.
        Optimized version that processes the file in parallel chunks.
        
        Args:
            gff_path: Path to the GFF file
            
        Returns:
            List of gene dictionaries with extracted information
            
        Raises:
            FileError: If GFF file cannot be read
            SequenceProcessingError: If GFF processing fails
        """
        logger.debug(f"Loading genes from GFF file: {gff_path}")
        
        try:
            # Read all lines from the file
            with open(gff_path, 'r') as f:
                all_lines = f.readlines()
            
            # Calculate chunk size for parallel processing
            chunk_size = max(1, len(all_lines) // Config.NUM_PROCESSES)
            chunks = CommonUtils.chunks(all_lines, chunk_size)
            
            logger.debug(f"Processing GFF in {len(chunks)} chunks with {Config.NUM_PROCESSES} processes")
            
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
            return genes
        
        except Exception as e:
            error_msg = f"Failed to load genes from GFF: {e}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return []

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
        
        logger.debug(f"Processing {len(restriction_fragments)} restriction fragments")
        logger.debug(f"Available genes: {len(genes)}")
        
        # Log restriction fragment statistics
        fragment_lengths = [len(frag.get('sequence', '')) for frag in restriction_fragments]
        if fragment_lengths:
            avg_len = sum(fragment_lengths) / len(fragment_lengths)
            logger.debug(f"Fragment length stats - avg: {avg_len:.0f}, min: {min(fragment_lengths)}, max: {max(fragment_lengths)}")
        
        fragments_with_genes = 0
        total_gene_overlaps = 0
        
        # Track what gets filtered out for debugging
        debug_stats = {
            'fragments_processed': 0,
            'fragments_with_overlaps': 0,
            'gene_fragments_created': 0,
            'gene_fragments_too_short': 0,
            'gene_fragments_invalid': 0,
            'total_overlaps_found': 0
        }
        
        for fragment in restriction_fragments:
            debug_stats['fragments_processed'] += 1
            
            fragment_chr = fragment.get("chr")
            fragment_start = fragment.get("start")
            fragment_end = fragment.get("end")
            fragment_seq = fragment.get("sequence", "")
            
            if not all([fragment_chr, fragment_start is not None, fragment_end is not None]):
                logger.warning(f"Fragment {fragment.get('id', 'unknown')} missing position data, skipping")
                continue
            
            # Find all genes that overlap with this fragment
            overlapping_genes = cls.find_overlapping_genes(
                fragment_chr, fragment_start, fragment_end, genes, overlap_margin
            )
            
            if overlapping_genes:
                debug_stats['fragments_with_overlaps'] += 1
                debug_stats['total_overlaps_found'] += len(overlapping_genes)
                fragments_with_genes += 1
                total_gene_overlaps += len(overlapping_genes)
                
                logger.debug(f"Fragment {fragment.get('id')}: found {len(overlapping_genes)} overlapping genes")
                
                # Create a separate fragment for each overlapping gene
                for gene_idx, gene in enumerate(overlapping_genes):
                    gene_fragment = cls.create_gene_fragment(
                        fragment, gene, gene_idx, overlap_margin
                    )
                    
                    if gene_fragment:
                        debug_stats['gene_fragments_created'] += 1
                        gene_fragments.append(gene_fragment)
                        logger.debug(f"Created gene fragment: {gene_fragment['id']} "
                                f"(Gene: {gene['id']}, Length: {len(gene_fragment['sequence'])} bp)")
                    else:
                        debug_stats['gene_fragments_invalid'] += 1
                        # Try to determine why it was invalid
                        gene_start = gene.get("start")
                        gene_end = gene.get("end")
                        if gene_start and gene_end:
                            potential_length = min(fragment_end, gene_end) - max(fragment_start, gene_start)
                            if potential_length < Config.MIN_SEGMENT_LENGTH:
                                debug_stats['gene_fragments_too_short'] += 1
                                logger.debug(f"Gene fragment for {gene.get('id', 'unknown')} rejected: "
                                              f"potential length {potential_length} < {Config.MIN_SEGMENT_LENGTH} bp minimum")
                        logger.debug(f"Failed to create gene fragment for {gene.get('id', 'unknown')}")
            else:
                logger.debug(f"Fragment {fragment.get('id')}: no overlapping genes found")
        
        # Log detailed statistics
        logger.debug(f"GENE OVERLAP EXTRACTION DETAILED STATS:")
        logger.debug(f"  Restriction fragments processed: {debug_stats['fragments_processed']}")
        logger.debug(f"  Fragments with gene overlaps: {debug_stats['fragments_with_overlaps']}")
        logger.debug(f"  Total gene overlaps found: {debug_stats['total_overlaps_found']}")
        logger.debug(f"  Gene fragments created successfully: {debug_stats['gene_fragments_created']}")
        logger.debug(f"  Gene fragments rejected (too short): {debug_stats['gene_fragments_too_short']}")
        logger.debug(f"  Gene fragments rejected (other reasons): {debug_stats['gene_fragments_invalid'] - debug_stats['gene_fragments_too_short']}")
        logger.debug(f"  Success rate: {debug_stats['gene_fragments_created']}/{debug_stats['total_overlaps_found']} ({100*debug_stats['gene_fragments_created']/max(1, debug_stats['total_overlaps_found']):.1f}%)")
        
        logger.debug(f"GENE OVERLAP RESULTS:")
        logger.debug(f"  Input restriction fragments: {len(restriction_fragments)}")
        logger.debug(f"  Fragments with gene overlaps: {fragments_with_genes}")
        logger.debug(f"  Total gene overlaps found: {total_gene_overlaps}")
        logger.debug(f"  Output gene fragments: {len(gene_fragments)}")
        logger.debug(f"  Average overlaps per fragment: {total_gene_overlaps/len(restriction_fragments):.2f}")
        
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
        
        return overlapping_genes

    @classmethod
    def create_gene_fragment(cls, fragment: Dict, gene: Dict, gene_idx: int, 
                            overlap_margin: int = 0) -> Dict:
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
            logger.warning(f"Invalid overlap region for gene {gene_id}")
            return None
        
        # Calculate positions within the fragment sequence
        seq_start = overlap_start - fragment_start
        seq_end = overlap_end - fragment_start
        
        # Validate sequence boundaries
        if seq_start < 0 or seq_end > len(fragment_seq) or seq_start >= seq_end:
            logger.warning(f"Invalid sequence boundaries for gene {gene_id}: "
                          f"seq_start={seq_start}, seq_end={seq_end}, fragment_len={len(fragment_seq)}")
            return None
        
        # Extract the overlapping sequence
        gene_sequence = fragment_seq[seq_start:seq_end]
        
        # Use Config minimum segment length instead of hardcoded value
        if len(gene_sequence) < Config.MIN_SEGMENT_LENGTH:
            logger.debug(f"Skipping short gene fragment for {gene_id}: {len(gene_sequence)} bp < {Config.MIN_SEGMENT_LENGTH} bp minimum")
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
        
        return gene_fragment

    @classmethod
    def filter_by_gene_overlap_enhanced(cls, restriction_fragments: List[Dict], 
                                       genes: List[Dict]) -> List[Dict]:
        """
        Enhanced version of filter_by_gene_overlap that extracts ALL gene overlaps.
        
        This is the replacement function for the existing filter_by_gene_overlap
        in SequenceProcessor.
        
        Args:
            restriction_fragments: List of restriction fragments
            genes: List of gene annotations
            
        Returns:
            List of all gene-overlapping fragments
        """
        overlap_margin = getattr(Config, 'GENE_OVERLAP_MARGIN', 0)
        
        logger.debug(f"Extracting gene-overlapping regions with margin: {overlap_margin} bp")
        
        # Extract all gene overlaps
        gene_fragments = cls.extract_all_gene_overlaps(
            restriction_fragments, genes, overlap_margin
        )
        
        # Log statistics
        logger.debug(f"Gene overlap extraction results:")
        logger.debug(f"  Input restriction fragments: {len(restriction_fragments)}")
        logger.debug(f"  Output gene fragments: {len(gene_fragments)}")
        
        if gene_fragments:
            # Count unique genes
            unique_genes = set(frag.get("Gene", "unknown") for frag in gene_fragments)
            logger.debug(f"  Unique genes covered: {len(unique_genes)}")
            
            # Log fragment length statistics
            lengths = [len(frag.get("sequence", "")) for frag in gene_fragments]
            if lengths:
                avg_length = sum(lengths) / len(lengths)
                min_length = min(lengths)
                max_length = max(lengths)
                logger.debug(f"  Fragment lengths: avg={avg_length:.0f}, min={min_length}, max={max_length}")
        
        return gene_fragments