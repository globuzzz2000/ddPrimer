#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotation processing module for ddPrimer pipeline.

Handles gene annotation processing and filtering operations including GFF file parsing,
gene overlap detection, and fragment filtering for primer design workflows.

Contains functionality for:
1. GFF file parsing and gene annotation loading
2. Gene overlap detection and coordinate handling
3. Fragment filtering based on gene annotations
4. Coordinate system validation and conversion

COORDINATE SYSTEM:
- Internal coordinates: 0-based half-open intervals [start, end)
- GFF file coordinates: 1-based inclusive intervals [start, end]
- Conversion: GFF [1001, 2000] becomes internal [1000, 2000)
- Half-open overlap: intervals [a,b) and [c,d) overlap if a < d and b > c
"""

import concurrent.futures
import gzip
import logging
import os
from tqdm import tqdm
from typing import List, Dict, Optional

# Import package modules
from ..config import Config, FileError, SequenceProcessingError, CoordinateValidationError, DebugLogLimiter

# Set up module logger
logger = logging.getLogger(__name__)


class AnnotationProcessor:
    """
    Handles gene annotation processing and filtering operations.
    
    This class provides comprehensive gene annotation capabilities
    including GFF parsing, gene overlap detection, and coordinate-based
    filtering for primer design workflows.
    
    Uses 0-based half-open intervals [start, end) internally for consistency
    with Python string indexing and bioinformatics standards.
    """
    
    #############################################################################
    #                           Workflow Wrappers
    #############################################################################
    
    @classmethod
    def filter_fragments_by_gene_overlap_workflow(cls, restriction_fragments: List[Dict], 
                                                genes: List[Dict], 
                                                skip_annotation_filtering: bool = False) -> List[Dict]:
        """
        Filter fragments based on gene overlap for workflow integration.
        
        Removed SNP processing statistics logging - that should be handled by SNP processor.
        """
        logger.info("Extracting annotation-filtered regions...")    
        logger.debug("=== WORKFLOW: GENE OVERLAP FILTERING ===")
        logger.debug(f"Processing {len(restriction_fragments)} fragments for gene overlap")
        logger.debug(f"Skip annotation filtering: {skip_annotation_filtering}")
        
        try:
            if skip_annotation_filtering:
                filtered_fragments = []
                for fragment in restriction_fragments:
                    # Validate fragment coordinates before processing
                    validation_result = cls._validate_fragment_coordinates(fragment)
                    if not validation_result["valid"]:
                        logger.warning(f"Fragment {fragment.get('id', 'unknown')} has invalid coordinates: {validation_result['errors']}")
                        continue
                    
                    # Create simplified fragment, preserving SNP processing stats if present
                    simplified_fragment = {
                        "id": fragment["id"],
                        "sequence": fragment["sequence"],
                        "chromosome": fragment.get("chromosome", ""),
                        "start": fragment.get("start", 0),  # 0-based
                        "end": fragment.get("end", len(fragment["sequence"])),  # 0-based
                        "length": len(fragment["sequence"]),
                        "Gene": fragment["id"].split("_")[-1]
                    }
                    
                    # Preserve SNP processing statistics if they exist (but don't log them here)
                    if '_processing_stats' in fragment:
                        simplified_fragment['_processing_stats'] = fragment['_processing_stats']
                    
                    filtered_fragments.append(simplified_fragment)
                    
                logger.debug(f"Created {len(filtered_fragments)} simplified fragments (no annotation filtering)")
                
            else:
                if not genes:
                    logger.error("Gene annotations not provided for standard mode.")
                    return []

                # Validate all gene coordinates before processing
                valid_genes = []
                for gene in genes:
                    validation_result = cls._validate_gene_coordinates(gene)
                    if validation_result["valid"]:
                        valid_genes.append(gene)
                    else:
                        logger.warning(f"Gene {gene.get('id', 'unknown')} has invalid coordinates: {validation_result['errors']}")

                filtered_fragments = cls.filter_by_gene_overlap(restriction_fragments, valid_genes)
            
            logger.info(f"Generated {len(filtered_fragments)} annotation-filtered fragments\n")
            logger.debug("=== END WORKFLOW: GENE OVERLAP FILTERING ===")
            
            return filtered_fragments
            
        except Exception as e:
            error_msg = f"Error in gene overlap filtering workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: GENE OVERLAP FILTERING ===")
            raise SequenceProcessingError(error_msg) from e
    
    #############################################################################
    #                           Coordinate Validation
    #############################################################################
    
    @classmethod
    def _validate_fragment_coordinates(cls, fragment: Dict) -> Dict[str, any]:
        """
        Validate fragment coordinates for 0-based half-open interval consistency.
        
        Args:
            fragment: Fragment dictionary with 0-based coordinates
            
        Returns:
            Dictionary with validation results
        """
        errors = []
        
        # Check required fields
        start = fragment.get("start")
        end = fragment.get("end") 
        sequence = fragment.get("sequence", "")
        chromosome = fragment.get("chromosome")
        
        if start is None:
            errors.append("Missing 'start' coordinate")
        if end is None:
            errors.append("Missing 'end' coordinate")
        if not chromosome:
            errors.append("Missing 'chromosome' field")
        if not sequence:
            errors.append("Missing or empty 'sequence' field")
            
        # Validate coordinate logic (0-based half-open system)
        if start is not None and end is not None:
            if start < 0:
                errors.append(f"Start coordinate cannot be negative: {start}")
            if end <= start:
                errors.append(f"End coordinate ({end}) must be greater than start ({start})")
            if sequence and (end - start) != len(sequence):
                errors.append(f"Coordinate span ({end - start}) doesn't match sequence length ({len(sequence)})")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "coordinates": {
                "start": start,
                "end": end,
                "length": end - start if (start is not None and end is not None) else None,
                "sequence_length": len(sequence)
            }
        }
    
    @classmethod
    def _validate_gene_coordinates(cls, gene: Dict) -> Dict[str, any]:
        """
        Validate gene coordinates for 0-based half-open interval consistency.
        
        Args:
            gene: Gene dictionary with 0-based coordinates
            
        Returns:
            Dictionary with validation results
        """
        errors = []
        
        # Check required fields
        start = gene.get("start")
        end = gene.get("end")
        chromosome = gene.get("chr")
        gene_id = gene.get("id")
        
        if start is None:
            errors.append("Missing 'start' coordinate")
        if end is None:
            errors.append("Missing 'end' coordinate")
        if not chromosome:
            errors.append("Missing 'chr' field")
        if not gene_id:
            errors.append("Missing 'id' field")
            
        # Validate coordinate logic (0-based half-open system)
        if start is not None and end is not None:
            if start < 0:
                errors.append(f"Start coordinate cannot be negative: {start}")
            if end <= start:
                errors.append(f"End coordinate ({end}) must be greater than start ({start})")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "coordinates": {
                "start": start,
                "end": end,
                "length": end - start if (start is not None and end is not None) else None
            }
        }
    
    #############################################################################
    #                           GFF Processing
    #############################################################################
    
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
    def _convert_gff_to_internal_coordinates(cls, gff_start: int, gff_end: int) -> tuple[int, int]:
        """
        Convert GFF coordinates to internal 0-based half-open coordinates.
        
        GFF format: 1-based inclusive [start, end]
        Internal format: 0-based half-open [start, end)
        
        Args:
            gff_start: GFF start position (1-based, inclusive)
            gff_end: GFF end position (1-based, inclusive)
            
        Returns:
            Tuple of (internal_start, internal_end) in 0-based half-open format
        """
        # COORDINATE CONVERSION: GFF (1-based inclusive) -> Internal (0-based half-open)
        internal_start = gff_start - 1  # Convert to 0-based
        internal_end = gff_end          # Keep same value for half-open interval
        
        if logger.isEnabledFor(logging.DEBUG) and DebugLogLimiter.should_log('gff_coordinate_conversion', interval=1000, max_initial=3):
            logger.debug(f"COORDINATE CONVERSION: GFF [{gff_start}, {gff_end}] -> Internal [{internal_start}, {internal_end})")
        
        return internal_start, internal_end
    
    @classmethod
    def process_gff_chunk(cls, chunk: List[str]) -> List[Dict]:
        """
        Process a chunk of GFF file lines with coordinate conversion.
        Used for parallel processing.
        
        Args:
            chunk: List of GFF file lines
            
        Returns:
            List of gene dictionaries with 0-based internal coordinates
        """
        chunk_genes = []
        
        if isinstance(Config.RETAIN_TYPES, str):
            retain_types = [t.strip().lower() for t in Config.RETAIN_TYPES.split(',')]
        elif isinstance(Config.RETAIN_TYPES, list):
            retain_types = [t.lower() for t in Config.RETAIN_TYPES]
        else:
            retain_types = [str(Config.RETAIN_TYPES).lower()]
        
        processed_count = 0
        valid_features = 0
        conversion_errors = 0
        
        for line in chunk:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqname, _, ftype, start, end, _, strand, _, attributes = parts
            processed_count += 1
            
            if ftype.lower() not in retain_types:
                continue

            attr_dict = cls.parse_gff_attributes(attributes)
            name = attr_dict.get('name')
            gene_id = attr_dict.get('id')
            locus_tag = attr_dict.get('locus_tag')
            identifier = name or gene_id or locus_tag or f"feature_{start}_{end}"

            try:
                # Parse GFF coordinates (1-based inclusive)
                gff_start = int(start)
                gff_end = int(end)
                
                # COORDINATE CONVERSION: Convert to internal 0-based half-open
                internal_start, internal_end = cls._convert_gff_to_internal_coordinates(gff_start, gff_end)
                
                # Create gene record with 0-based coordinates
                gene_record = {
                    "chr": seqname,
                    "start": internal_start,  # 0-based start
                    "end": internal_end,      # 0-based end (exclusive)
                    "strand": strand,
                    "id": identifier,
                    # Keep original GFF coordinates for debugging/validation
                    "_original_gff_start": gff_start,
                    "_original_gff_end": gff_end
                }
                
                # Validate the converted coordinates
                validation_result = cls._validate_gene_coordinates(gene_record)
                if validation_result["valid"]:
                    chunk_genes.append(gene_record)
                    valid_features += 1
                    
                    if (logger.isEnabledFor(logging.DEBUG) and 
                        DebugLogLimiter.should_log('gff_feature_processing', interval=1000, max_initial=3)):
                        logger.debug(f"Processed feature {valid_features}: {identifier} "
                                   f"GFF[{gff_start},{gff_end}] -> Internal[{internal_start},{internal_end})")
                else:
                    conversion_errors += 1
                    if (logger.isEnabledFor(logging.DEBUG) and 
                        DebugLogLimiter.should_log('gff_validation_errors', interval=100, max_initial=3)):
                        logger.debug(f"Invalid coordinates for feature {identifier}: {validation_result['errors']}")
                    
            except ValueError as e:
                conversion_errors += 1
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('gff_processing_errors', interval=100, max_initial=2)):
                    logger.debug(f"Error converting position data for feature {identifier}: {e}")
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Chunk processing complete: {processed_count} lines processed, "
                       f"{valid_features} valid genes extracted, {conversion_errors} conversion errors")
        
        return chunk_genes
    
    @classmethod
    def load_genes_from_gff(cls, gff_path: str) -> List[Dict]:
        """
        Extract gene information from a GFF file with coordinate conversion.
        
        Handles both compressed (.gz) and uncompressed GFF files.
        Uses parallel processing for large files without progress bar.
        """
        if not os.path.exists(gff_path):
            error_msg = f"GFF file not found: {gff_path}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        try:
            is_compressed = gff_path.endswith('.gz')
            opener = gzip.open if is_compressed else open
            mode = 'rt' if is_compressed else 'r'
            
            with opener(gff_path, mode) as f:
                all_lines = f.readlines()
            
            logger.debug(f"Read {len(all_lines)} lines from GFF file")
            logger.debug("COORDINATE SYSTEM: Converting GFF (1-based inclusive) to Internal (0-based half-open)")
            
            chunk_size = max(1, len(all_lines) // Config.NUM_PROCESSES)
            chunks = [all_lines[i:i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
            
            logger.debug(f"Split into {len(chunks)} chunks for parallel processing")
            
            genes = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=Config.NUM_PROCESSES) as executor:
                futures = [executor.submit(cls.process_gff_chunk, chunk) for chunk in chunks]
                
                # Process without progress bar
                for future in concurrent.futures.as_completed(futures):
                    genes.extend(future.result())
            
            logger.debug(f"Extracted {len(genes)} total genes from GFF file (0-based coordinates)")
            
            if logger.isEnabledFor(logging.DEBUG) and genes:
                logger.debug("Sample genes extracted (showing internal 0-based coordinates):")
                for i, gene in enumerate(genes[:5]):
                    original_start = gene.get('_original_gff_start', 'N/A')
                    original_end = gene.get('_original_gff_end', 'N/A')
                    logger.debug(f"  {i+1}. {gene['id']} at {gene['chr']}:{gene['start']}-{gene['end']} "
                            f"(was GFF:{original_start}-{original_end})")
                if len(genes) > 5:
                    logger.debug(f"  ... and {len(genes) - 5} more genes")
            
            return genes
        
        except Exception as e:
            error_msg = f"Failed to load genes from GFF: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e

    #############################################################################
    #                           Gene Overlap Processing
    #############################################################################
    
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
        parts = sequence_id.split("_")
        
        if len(parts) >= 3:
            return parts[2]
        
        return sequence_id

    @classmethod
    def filter_by_gene_overlap(cls, restriction_fragments: List[Dict], 
                              genes: List[Dict]) -> List[Dict]:
        """
        Filter restriction fragments by gene overlap and extract gene regions.
        
        This is the main function for extracting gene-overlapping regions from
        restriction fragments. Creates separate fragments for each gene overlap.
        All coordinates are in 0-based format.
        
        Args:
            restriction_fragments: List of restriction fragments (0-based coordinates)
            genes: List of gene annotations (0-based coordinates)
            
        Returns:
            List of all gene-overlapping fragments (0-based coordinates)
        """
        overlap_margin = getattr(Config, 'GENE_OVERLAP_MARGIN', 0)
        
        gene_fragments = cls.extract_all_gene_overlaps(
            restriction_fragments, genes, overlap_margin
        )
        
        logger.debug(f"Extracted {len(gene_fragments)} gene fragments from {len(restriction_fragments)} restriction fragments")
        return gene_fragments

    @classmethod
    def extract_all_gene_overlaps(cls, restriction_fragments: List[Dict], genes: List[Dict], 
                                 overlap_margin: int = 0) -> List[Dict]:
        """
        Extract ALL gene-overlapping regions from each restriction fragment.
        
        This function processes each restriction fragment and identifies all genes
        that overlap with it, creating separate fragments for each gene overlap.
        All coordinates are in 0-based half-open format.
        
        Args:
            restriction_fragments: List of restriction fragments (0-based coordinates)
            genes: List of gene annotations (0-based coordinates)
            overlap_margin: Additional margin around genes (bp)
            
        Returns:
            List of gene-overlapping fragments (0-based coordinates)
        """
        logger.debug(f"=== GENE OVERLAP EXTRACTION DEBUG ===")
        logger.debug(f"Processing {len(restriction_fragments)} fragments for gene overlaps with {overlap_margin}bp margin")
        logger.debug("COORDINATE SYSTEM: All coordinates are 0-based half-open intervals")
        
        gene_fragments = []
        processed_fragments = 0
        fragments_with_overlaps = 0
        validation_errors = 0
        
        for fragment in restriction_fragments:
            # Validate fragment coordinates first
            fragment_validation = cls._validate_fragment_coordinates(fragment)
            if not fragment_validation["valid"]:
                validation_errors += 1
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('fragment_validation_errors', interval=100, max_initial=3)):
                    logger.debug(f"Fragment {fragment.get('id', 'unknown')} validation failed: {fragment_validation['errors']}")
                processed_fragments += 1
                continue
            
            fragment_chr = fragment.get("chromosome")
            fragment_start = fragment.get("start")  # 0-based
            fragment_end = fragment.get("end")      # 0-based
            
            overlapping_genes = cls.find_overlapping_genes(
                fragment_chr, fragment_start, fragment_end, genes, overlap_margin
            )
            
            if overlapping_genes:
                fragments_with_overlaps += 1
                for gene_idx, gene in enumerate(overlapping_genes):
                    gene_fragment = cls.create_gene_fragment(
                        fragment, gene, gene_idx, overlap_margin
                    )
                    
                    if gene_fragment:
                        # Validate the created gene fragment
                        gene_frag_validation = cls._validate_fragment_coordinates(gene_fragment)
                        if gene_frag_validation["valid"]:
                            gene_fragments.append(gene_fragment)
                        else:
                            validation_errors += 1
                            if (logger.isEnabledFor(logging.DEBUG) and 
                                DebugLogLimiter.should_log('gene_fragment_validation_errors', interval=100, max_initial=3)):
                                logger.debug(f"Created gene fragment validation failed: {gene_frag_validation['errors']}")
                
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('fragment_overlap_success', interval=200, max_initial=3)):
                    logger.debug(f"Fragment {fragment['id']} ({fragment_chr}:{fragment_start}-{fragment_end}) "
                               f"overlaps with {len(overlapping_genes)} genes")
            else:
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('fragment_no_overlap', interval=500, max_initial=2)):
                    logger.debug(f"Fragment {fragment['id']} ({fragment_chr}:{fragment_start}-{fragment_end}) "
                               f"has no gene overlaps")
            
            processed_fragments += 1
        
        logger.debug(f"Gene overlap extraction complete: {processed_fragments} fragments processed, "
                   f"{fragments_with_overlaps} had overlaps, generated {len(gene_fragments)} gene fragments")
        if validation_errors > 0:
            logger.debug(f"Validation errors encountered: {validation_errors}")
        logger.debug(f"=== END GENE OVERLAP EXTRACTION DEBUG ===")
        
        return gene_fragments

    @classmethod
    def find_overlapping_genes(cls, fragment_chr: str, fragment_start: int, fragment_end: int,
                              genes: List[Dict], overlap_margin: int = 0) -> List[Dict]:
        """
        Find all genes that overlap with a given genomic region.
        
        Uses 0-based half-open coordinate system for all calculations.
        Two intervals [a,b) and [c,d) overlap if: a < d and b > c
        
        Args:
            fragment_chr: Chromosome/sequence name of the fragment
            fragment_start: Start position of the fragment (0-based)
            fragment_end: End position of the fragment (0-based, exclusive)
            genes: List of gene annotations (0-based coordinates)
            overlap_margin: Additional margin around genes
            
        Returns:
            List of overlapping genes (0-based coordinates)
        """
        overlapping_genes = []
        genes_checked = 0
        overlap_calculations = 0
        
        for gene in genes:
            if gene.get("chr") != fragment_chr:
                continue
            
            gene_start = gene.get("start")  # 0-based
            gene_end = gene.get("end")      # 0-based
            
            if gene_start is None or gene_end is None:
                continue
            
            # Apply margin to gene coordinates (still in 0-based space)
            gene_start_with_margin = max(0, gene_start - overlap_margin)
            gene_end_with_margin = gene_end + overlap_margin
            
            # HALF-OPEN INTERVAL OVERLAP: [a,b) and [c,d) overlap if a < d and b > c
            overlap_detected = fragment_start < gene_end_with_margin and fragment_end > gene_start_with_margin
            
            if overlap_detected:
                # Calculate overlap region for debugging
                overlap_start = max(fragment_start, gene_start_with_margin)
                overlap_end = min(fragment_end, gene_end_with_margin)
                overlap_length = overlap_end - overlap_start
                
                overlapping_genes.append(gene)
                    
                if (logger.isEnabledFor(logging.DEBUG) and 
                    DebugLogLimiter.should_log('gene_overlap_details', interval=1000, max_initial=2)):
                    logger.debug(f"Gene {gene['id']} ({gene_start}-{gene_end}) overlaps with fragment "
                               f"({fragment_start}-{fragment_end}), overlap region: "
                               f"{{'start': {overlap_start}, 'end': {overlap_end}, 'length': {overlap_length}}}")
                
                overlap_calculations += 1
            
            genes_checked += 1
        
        if (logger.isEnabledFor(logging.DEBUG) and 
            DebugLogLimiter.should_log('gene_overlap_summary', interval=500, max_initial=2)):
            logger.debug(f"Checked {genes_checked} genes, found {len(overlapping_genes)} overlaps, "
                       f"performed {overlap_calculations} overlap calculations")
        
        return overlapping_genes

    @classmethod
    def create_gene_fragment(cls, fragment: Dict, gene: Dict, gene_idx: int, 
                            overlap_margin: int = 0) -> Optional[Dict]:
        """
        Create a new fragment for a specific gene overlap region.

        Uses 0-based half-open coordinates throughout. Creates a fragment that represents
        the intersection of the original fragment and the gene region (with margin).
        
        Args:
            fragment: Original restriction fragment (0-based coordinates, may have SNP stats)
            gene: Gene annotation that overlaps with the fragment (0-based coordinates)
            gene_idx: Index of this gene (for unique naming)
            overlap_margin: Overlap margin used
            
        Returns:
            New fragment dictionary for the gene region (0-based coordinates), or None if invalid
        """
        fragment_chr = fragment.get("chromosome")
        fragment_start = fragment.get("start")  # 0-based
        fragment_end = fragment.get("end")      # 0-based
        fragment_seq = fragment.get("sequence", "")
        
        gene_start = gene.get("start")  # 0-based
        gene_end = gene.get("end")      # 0-based
        gene_id = gene.get("id", "unknown_gene")
        
        if not all([fragment_start is not None, fragment_end is not None, 
                    gene_start is not None, gene_end is not None]):
            return None
        
        # Apply margin to gene coordinates
        gene_start_with_margin = max(0, gene_start - overlap_margin)
        gene_end_with_margin = gene_end + overlap_margin
        
        # Calculate overlap region in 0-based half-open coordinates
        overlap_start = max(fragment_start, gene_start_with_margin)
        overlap_end = min(fragment_end, gene_end_with_margin)
        
        if overlap_start >= overlap_end:
            return None
        
        # Calculate sequence boundaries (0-based indexing into sequence)
        seq_start = overlap_start - fragment_start
        seq_end = overlap_end - fragment_start
        
        if seq_start < 0 or seq_end > len(fragment_seq) or seq_start >= seq_end:
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('gene_fragment_validation_errors', interval=200, max_initial=2)):
                logger.debug(f"Invalid sequence boundaries for gene {gene_id}: "
                            f"seq_start={seq_start}, seq_end={seq_end}, fragment_len={len(fragment_seq)}")
            return None
        
        gene_sequence = fragment_seq[seq_start:seq_end]
        
        if len(gene_sequence) < Config.MIN_SEGMENT_LENGTH:
            if (logger.isEnabledFor(logging.DEBUG) and 
                DebugLogLimiter.should_log('gene_fragment_too_short', interval=200, max_initial=2)):
                logger.debug(f"Gene fragment too short for {gene_id}: {len(gene_sequence)} < {Config.MIN_SEGMENT_LENGTH}")
            return None
        
        original_id = fragment.get("id", "unknown")
        new_id = f"{original_id}_gene{gene_idx}_{gene_id}"
        
        # Create gene fragment with 0-based coordinates
        gene_fragment = {
            "id": new_id,
            "sequence": gene_sequence,
            "chromosome": fragment_chr,
            "start": overlap_start,    # 0-based start
            "end": overlap_end,        # 0-based end (exclusive)
            "length": len(gene_sequence),
            "Gene": gene_id,
            "original_fragment": original_id,
            "gene_start": gene_start,  # Original gene coordinates (0-based)
            "gene_end": gene_end,      # Original gene coordinates (0-based)
            "overlap_margin": overlap_margin,
            # Validation metadata
            "_fragment_coordinates": {
                "original_start": fragment_start,
                "original_end": fragment_end,
                "seq_extraction": f"[{seq_start}:{seq_end}]"
            }
        }
        
        # UPDATED: Preserve SNP processing statistics if they exist in the original fragment
        if '_processing_stats' in fragment:
            # Calculate proportional statistics for the gene fragment
            original_length = fragment_end - fragment_start
            gene_fragment_length = len(gene_sequence)
            
            if original_length > 0:
                # Scale the statistics proportionally to the fragment size
                proportion = gene_fragment_length / original_length
                original_stats = fragment['_processing_stats']
                
                scaled_stats = {
                    'applied': int(original_stats.get('applied', 0) * proportion),
                    'substituted': int(original_stats.get('substituted', 0) * proportion),
                    'masked': int(original_stats.get('masked', 0) * proportion),
                    'validation_failures': int(original_stats.get('validation_failures', 0) * proportion),
                    'inherited_from': original_id,
                    'scaling_factor': proportion
                }
                
                gene_fragment['_processing_stats'] = scaled_stats
        
        if (logger.isEnabledFor(logging.DEBUG) and 
            DebugLogLimiter.should_log('gene_fragment_creation_success', interval=500, max_initial=3)):
            stats_info = ""
            if '_processing_stats' in gene_fragment:
                stats = gene_fragment['_processing_stats']
                stats_info = f" (inherited {stats['applied']} SNP variants)"
            logger.debug(f"Created gene fragment {new_id}: {len(gene_sequence)} bp "
                    f"(coordinates: {overlap_start}-{overlap_end}){stats_info}")
        
        return gene_fragment