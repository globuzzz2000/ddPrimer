# MAFParser module for ddPrimer pipeline
# For parsing MAF alignment files and masking sequences
# Integrated with ddPrimer by Jakob (2025)

"""
MAF Parser module for ddPrimer pipeline.

This module provides functionality for parsing MAF files from LastZ alignments
and masking sequences based on alignment data for primer design.
"""

import os
import re
import logging
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from ..config import Config
from ..config.exceptions import AlignmentError


class MAFParser:
    """
    Parser for Multiple Alignment Format (MAF) files.
    Handles parsing alignment data and sequence masking for primer design.
    """
    
    def __init__(self, config=None):
        """
        Initialize with configuration settings.
        
        Args:
            config: Configuration object (defaults to global Config)
        """
        self.config = config if config else Config
        self.logger = logging.getLogger("ddPrimer.helpers")
        self.alignments = defaultdict(list)  # Store alignments by reference sequence
        self.masked_regions = defaultdict(list)  # Store regions to mask
    
    def parse_maf_file(self, maf_file):
        """
        Parse a MAF file and extract alignment information.
        
        Args:
            maf_file (str): Path to MAF file
            
        Returns:
            dict: Parsed alignment data organized by reference sequence
            
        Raises:
            AlignmentError: If there's an error parsing the MAF file
        """
        self.logger.debug(f"Parsing MAF file: {maf_file}")
        
        if not os.path.exists(maf_file):
            self.logger.error(f"MAF file not found: {maf_file}")
            raise AlignmentError(f"MAF file not found: {maf_file}")
            
        try:
            current_alignment = []
            with open(maf_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue
                    
                    # New alignment block
                    if line.startswith('a'):
                        # Process previous alignment if it exists
                        if current_alignment:
                            self._process_alignment_block(current_alignment)
                            current_alignment = []
                        # Start new alignment with header line
                        current_alignment.append(line)
                    
                    # Sequence lines
                    elif line.startswith('s'):
                        current_alignment.append(line)
                    
                    # End of alignment block (blank line or new 'a' line)
                    elif current_alignment and not line:
                        self._process_alignment_block(current_alignment)
                        current_alignment = []
                
                # Process the last alignment block if it exists
                if current_alignment:
                    self._process_alignment_block(current_alignment)
            
            self.logger.debug(f"Parsed {len(self.alignments)} reference sequences with alignments")
            return self.alignments
            
        except Exception as e:
            self.logger.error(f"Error parsing MAF file: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Failed to parse MAF file: {str(e)}")
    
    def _process_alignment_block(self, alignment_lines):
        """
        Process a single alignment block from the MAF file.
        
        Args:
            alignment_lines (list): Lines from a single alignment block
        """
        if len(alignment_lines) < 3:  # Need at least a header and two sequences
            return
        
        # Extract sequence information from the alignment block
        sequences = []
        for line in alignment_lines[1:]:  # Skip the 'a' line
            if line.startswith('s'):
                parts = re.split(r'\s+', line.strip())
                if len(parts) >= 7:
                    seq_info = {
                        'src': parts[1],        # Use full sequence name (e.g., "CP002688.1")
                        'start': int(parts[2]),
                        'size': int(parts[3]),
                        'strand': parts[4],
                        'src_size': int(parts[5]),
                        'text': parts[6]
                    }
                    sequences.append(seq_info)
        
        # Need at least two sequences to form an alignment
        if len(sequences) < 2:
            return
        
        # First sequence is always the reference sequence in our case
        ref_seq = sequences[0]
        
        # Use the full sequence name as the reference chromosome
        ref_chrom = ref_seq['src']
        
        # Record the alignment for this reference sequence
        for qry_seq in sequences[1:]:
            alignment = {
                'ref_chrom': ref_chrom,      # Full sequence name
                'ref_start': ref_seq['start'],
                'ref_end': ref_seq['start'] + ref_seq['size'],
                'ref_strand': ref_seq['strand'],
                'qry_src': qry_seq['src'],    # Full query sequence name
                'qry_start': qry_seq['start'],
                'qry_end': qry_seq['start'] + qry_seq['size'],
                'qry_strand': qry_seq['strand'],
                'ref_seq': ref_seq['text'],
                'qry_seq': qry_seq['text']
            }
            
            # Calculate alignment identity
            alignment['identity'] = self._calculate_identity(ref_seq['text'], qry_seq['text'])
            
            # Add to alignments dictionary using full sequence name
            self.alignments[ref_chrom].append(alignment)
    
    def _calculate_identity(self, seq1, seq2):
        """
        Calculate sequence identity between two aligned sequences.
        
        Args:
            seq1 (str): First aligned sequence with gaps
            seq2 (str): Second aligned sequence with gaps
            
        Returns:
            float: Percentage identity (0-100)
        """
        min_len = min(len(seq1), len(seq2))
        
        matches = 0
        valid_positions = 0
        
        for i in range(min_len):
            # Skip positions with gaps in either sequence
            if seq1[i] == '-' or seq2[i] == '-':
                continue
            
            valid_positions += 1
            if seq1[i].upper() == seq2[i].upper():
                matches += 1
        
        if valid_positions == 0:
            return 0.0
            
        # Use exact multiplication before division to maintain precision
        return (matches * 100.0) / valid_positions
    
    def identify_conserved_regions(self, min_identity=80, min_length=20):
        """
        Identify regions that are conserved between genomes.
        
        Args:
            min_identity (float): Minimum sequence identity percentage (default: 80)
            min_length (int): Minimum length of conserved region (default: 20)
            
        Returns:
            dict: Conserved regions by chromosome
            
        Raises:
            AlignmentError: If no alignments are loaded or conservation analysis fails
        """
        if not self.alignments:
            self.logger.error("No alignments loaded. Call parse_maf_file() first.")
            raise AlignmentError("No alignments loaded. Call parse_maf_file() first.")
            
        self.logger.debug(f"Identifying conserved regions (min identity: {min_identity}%, min length: {min_length}bp)")
        
        try:
            conserved_regions = defaultdict(list)
            
            for chrom, alignments in self.alignments.items():
                for alignment in alignments:
                    # Skip alignments below the identity threshold
                    if alignment['identity'] < min_identity:
                        continue
                    
                    # Find conserved blocks within this alignment
                    blocks = self._find_conserved_blocks(
                        alignment['ref_seq'], 
                        alignment['qry_seq'],
                        alignment['ref_start'],
                        min_length
                    )
                    
                    for block in blocks:
                        # Map query coordinates for this block
                        qry_start = self._map_to_query_coords(
                            block['start'] - alignment['ref_start'],
                            alignment['ref_seq'],
                            alignment['qry_seq'],
                            alignment['qry_start']
                        )
                        
                        qry_end = self._map_to_query_coords(
                            block['end'] - alignment['ref_start'],
                            alignment['ref_seq'],
                            alignment['qry_seq'],
                            alignment['qry_start']
                        )
                        
                        conserved_regions[chrom].append({
                            'start': block['start'],
                            'end': block['end'],
                            'identity': block['identity'],
                            'qry_src': alignment['qry_src'],
                            'qry_start': qry_start,
                            'qry_end': qry_end,
                            'qry_strand': alignment['qry_strand']
                        })
            
            total_regions = sum(len(regions) for regions in conserved_regions.values())
            self.logger.debug(f"Identified {total_regions} conserved regions across {len(conserved_regions)} chromosomes")
            
            return conserved_regions
            
        except Exception as e:
            self.logger.error(f"Error identifying conserved regions: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Failed to identify conserved regions: {str(e)}")
    
    def _find_conserved_blocks(self, ref_seq, qry_seq, ref_start, min_length):
        """
        Find conserved blocks within an alignment.
        
        Args:
            ref_seq (str): Reference sequence with gaps
            qry_seq (str): Query sequence with gaps
            ref_start (int): Start position of the reference sequence
            min_length (int): Minimum length of conserved block
            
        Returns:
            list: Conserved blocks with start, end, and identity
        """
        blocks = []
        current_block = {'start': None, 'matches': 0, 'length': 0, 'ref_positions': 0}
        ref_pos = ref_start
        
        for i in range(len(ref_seq)):
            # Skip positions with gaps in either sequence
            if ref_seq[i] == '-' or qry_seq[i] == '-':
                # If we have an active block, check if it meets criteria
                if current_block['start'] is not None and current_block['ref_positions'] >= min_length:
                    identity = (current_block['matches'] / current_block['length']) * 100
                    blocks.append({
                        'start': current_block['start'],
                        'end': ref_pos - 1,  # Last position before the gap
                        'identity': identity
                    })
                
                # Reset current block
                current_block = {'start': None, 'matches': 0, 'length': 0, 'ref_positions': 0}
                
                # Only advance reference position if there's no gap in the reference
                if ref_seq[i] != '-':
                    ref_pos += 1
                
                continue
            
            # Start a new block if needed
            if current_block['start'] is None:
                current_block['start'] = ref_pos
            
            # Update block statistics
            current_block['length'] += 1
            current_block['ref_positions'] += 1
            if ref_seq[i].upper() == qry_seq[i].upper():
                current_block['matches'] += 1
            
            # Advance reference position
            ref_pos += 1
        
        # Check if the last block meets criteria
        if current_block['start'] is not None and current_block['ref_positions'] >= min_length:
            identity = (current_block['matches'] / current_block['length']) * 100
            blocks.append({
                'start': current_block['start'],
                'end': ref_pos - 1,
                'identity': identity
            })
        
        return blocks
    
    def _map_to_query_coords(self, ref_offset, ref_seq, qry_seq, qry_start):
        """
        Map a position in the reference sequence to the query sequence.
        
        Args:
            ref_offset (int): Offset in the reference sequence (ungapped)
            ref_seq (str): Reference sequence with gaps
            qry_seq (str): Query sequence with gaps
            qry_start (int): Start position of the query sequence
            
        Returns:
            int: Corresponding position in query sequence
        """
        ref_ungapped_pos = 0
        qry_ungapped_pos = 0
        
        for i in range(len(ref_seq)):
            if ref_ungapped_pos == ref_offset:
                return qry_start + qry_ungapped_pos
            
            # Only advance reference position if there's no gap
            if ref_seq[i] != '-':
                ref_ungapped_pos += 1
            
            # Only advance query position if there's no gap
            if qry_seq[i] != '-':
                qry_ungapped_pos += 1
        
        # If we reach this point, return the last position
        return qry_start + qry_ungapped_pos
    
    def extract_reference_sequences_from_maf(self):
        """
        Extract reference sequences directly from a MAF file.
        This is useful when no reference FASTA file is provided.
        
        Returns:
            dict: Dictionary mapping reference sequence IDs to sequences
            
        Raises:
            AlignmentError: If no alignments are loaded or extraction fails
        """
        self.logger.debug("Extracting reference sequences from MAF file")
        
        if not self.alignments:
            self.logger.error("No alignments loaded. Call parse_maf_file() first.")
            raise AlignmentError("No alignments loaded. Call parse_maf_file() first.")
        
        try:
            # Dictionary to store reference sequences
            reference_sequences = {}
            
            # Keep track of sequence lengths
            sequence_lengths = {}
            
            # For each chromosome/reference sequence in the alignments
            for chrom, alignments in self.alignments.items():
                self.logger.debug(f"Processing reference sequence: {chrom}")
                
                # Find the maximum coordinate to determine sequence length
                max_end = 0
                for alignment in alignments:
                    max_end = max(max_end, alignment['ref_end'])
                
                # Create a blank sequence of appropriate length
                sequence_lengths[chrom] = max_end
                sequence = ['N'] * max_end
                
                # Fill in the sequence from alignments
                for alignment in alignments:
                    ref_start = alignment['ref_start']
                    ref_seq = alignment['ref_seq']
                    
                    # Remove gaps from reference sequence
                    ref_seq_ungapped = ref_seq.replace('-', '')
                    
                    # Fill in the sequence at the appropriate position
                    for i, base in enumerate(ref_seq_ungapped):
                        pos = ref_start + i
                        if pos < len(sequence):
                            sequence[pos] = base
                
                # Convert to string and store
                reference_sequences[chrom] = ''.join(sequence)
                
                self.logger.debug(f"Extracted reference sequence {chrom}: {len(reference_sequences[chrom])} bp")
            
            # Check if we have empty sequences
            for chrom, seq in reference_sequences.items():
                if set(seq) == {'N'}:
                    self.logger.warning(f"Reference sequence {chrom} contains only N's")
            
            self.logger.debug(f"Extracted {len(reference_sequences)} reference sequences from MAF file")
            return reference_sequences
            
        except Exception as e:
            self.logger.error(f"Error extracting reference sequences from MAF: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Failed to extract reference sequences: {str(e)}")
    
    def mask_non_conserved_regions(self, input_file_or_dict, output_file, conserved_regions, min_identity=80):
        """
        Create a masked version of the reference genome keeping only conserved regions.
        
        Args:
            input_file_or_dict (str/dict): Input FASTA file or dictionary of sequences
            output_file (str): Output masked FASTA file
            conserved_regions (dict): Conserved regions by chromosome
            min_identity (float): Minimum identity to consider a region conserved
            
        Returns:
            str: Path to masked FASTA file
            
        Raises:
            AlignmentError: If masking fails
        """
        self.logger.debug(f"Masking non-conserved regions in reference genome")
        
        try:
            # Get sequences either from file or from provided dictionary
            sequences = self._load_input_sequences(input_file_or_dict)
            
            # Get the list of FASTA sequence IDs
            fasta_ids = {seq_id: len(sequence) for seq_id, sequence in sequences.items()}
            
            # Get the list of alignment chromosome names
            align_chroms = list(conserved_regions.keys())
            
            self.logger.debug(f"FASTA IDs ({len(fasta_ids)}): {', '.join(list(fasta_ids.keys())[:5])}...")
            self.logger.debug(f"Alignment sequence IDs ({len(align_chroms)}): {', '.join(align_chroms[:5])}...")
            
            # Create chromosome mapping between alignment and FASTA sequences
            chrom_map = self._map_chromosomes(align_chroms, fasta_ids.keys())
            
            # Log the mapping results
            self.logger.debug("Chromosome mapping results:")
            for ac, fid in chrom_map.items():
                self.logger.debug(f"  Mapped: '{ac}' => '{fid}'")
            
            # Initialize masks with all positions masked (value of 0)
            masks = {}
            for record_id, length in fasta_ids.items():
                masks[record_id] = [0] * length
            
            # Create reverse mapping from FASTA IDs to alignment chromosomes
            reverse_map = {}
            for ac, fid in chrom_map.items():
                if fid not in reverse_map:
                    reverse_map[fid] = []
                reverse_map[fid].append(ac)
            
            # Apply masking to sequences
            masked_records, stats = self._apply_region_masks(
                sequences, masks, reverse_map, conserved_regions, min_identity
            )
            
            # Write masked sequences to output file
            SeqIO.write(masked_records, output_file, "fasta")
            
            # Calculate and log overall statistics
            total_n_count = stats["total_n_count"]
            total_bases = stats["total_bases"]
            applied_regions = stats["applied_regions"]
            total_regions = stats["total_regions"]
            
            mask_percent = (total_n_count / total_bases) * 100 if total_bases > 0 else 0
            unmask_percent = 100 - mask_percent
            
            self.logger.debug(f"Created masked FASTA file: {output_file}")
            self.logger.debug(f"Total statistics: {applied_regions}/{total_regions} regions applied")
            self.logger.debug(f"Masking summary: {total_n_count}/{total_bases} bases masked to N ({mask_percent:.2f}%)")
            self.logger.debug(f"Conserved bases: {total_bases - total_n_count}/{total_bases} bases kept ({unmask_percent:.2f}%)")
            
            # Check if everything was masked
            if total_n_count == total_bases:
                self.logger.warning("ALL bases in output FASTA are masked (100% Ns)")
                self.logger.warning("This indicates no regions were successfully mapped between alignment and FASTA")
            
            return output_file
            
        except Exception as e:
            self.logger.error(f"Error masking non-conserved regions: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Failed to mask non-conserved regions: {str(e)}")
    
    def _load_input_sequences(self, input_file_or_dict):
        """
        Load sequences from a file or use provided dictionary.
        
        Args:
            input_file_or_dict (str/dict): Input FASTA file or dictionary of sequences
            
        Returns:
            dict: Dictionary of sequences
        """
        if isinstance(input_file_or_dict, str):
            # Load sequences from FASTA file
            sequences = {}
            for record in SeqIO.parse(input_file_or_dict, "fasta"):
                sequences[record.id] = str(record.seq)
            self.logger.debug(f"Loaded {len(sequences)} sequences from {input_file_or_dict}")
            return sequences
        else:
            # Use provided dictionary
            self.logger.debug(f"Using provided dictionary with {len(input_file_or_dict)} sequences")
            return input_file_or_dict
            
    def _apply_region_masks(self, sequences, masks, reverse_map, conserved_regions, min_identity):
        """
        Apply masking to sequences based on conserved regions.
        
        Args:
            sequences (dict): Dictionary of sequences
            masks (dict): Dictionary of mask arrays
            reverse_map (dict): Mapping from FASTA IDs to alignment chromosomes
            conserved_regions (dict): Conserved regions by chromosome
            min_identity (float): Minimum identity to consider a region conserved
            
        Returns:
            tuple: (masked_records, stats) - List of masked sequence records and statistics
        """
        # Count total regions for statistics
        total_regions = sum(len(regions) for regions in conserved_regions.values())
        applied_regions = 0
        unmasked_bases = 0
        
        # Analyze each FASTA sequence
        for fasta_id, sequence_length in {id: len(seq) for id, seq in sequences.items()}.items():
            # Get all alignment chromosomes that map to this FASTA ID
            corresponding_align_chroms = reverse_map.get(fasta_id, [])
            
            if not corresponding_align_chroms:
                self.logger.debug(f"No alignment chromosomes map to FASTA ID '{fasta_id}'")
                continue
                
            self.logger.debug(f"Processing FASTA ID '{fasta_id}', mapped to alignment chromosomes: {corresponding_align_chroms}")
            
            # Apply all relevant conserved regions
            regions_applied = 0
            bases_unmasked = 0
            
            for align_chrom in corresponding_align_chroms:
                regions = conserved_regions.get(align_chrom, [])
                self.logger.debug(f"  - Applying {len(regions)} regions from alignment chromosome '{align_chrom}'")
                
                # Use progress bar for large numbers of regions
                region_iter = regions
                if Config.SHOW_PROGRESS and len(regions) > 100:
                    region_iter = tqdm(regions, desc=f"Processing {align_chrom}", leave=False)
                
                for region in region_iter:
                    if region['identity'] >= min_identity:
                        regions_applied += 1
                        start = region['start']
                        end = min(region['end'], len(masks[fasta_id]) - 1)  # Ensure within bounds
                        
                        if start >= len(masks[fasta_id]):
                            continue  # Skip if start is beyond sequence length
                        
                        for i in range(start, end + 1):
                            if 0 <= i < len(masks[fasta_id]):
                                if masks[fasta_id][i] == 0:  # Only count if not already unmasked
                                    bases_unmasked += 1
                                masks[fasta_id][i] = 1
            
            applied_regions += regions_applied
            unmasked_bases += bases_unmasked
            self.logger.debug(f"Applied {regions_applied} conserved regions to {fasta_id}, unmasked {bases_unmasked} bases")
        
        # Create masked sequences
        masked_records = []
        total_n_count = 0
        total_bases = 0
        
        for seq_id, sequence in sequences.items():
            total_bases += len(sequence)
            
            if seq_id in masks:
                masked_seq = ""
                n_count = 0
                
                for i, base in enumerate(sequence):
                    if i < len(masks[seq_id]) and masks[seq_id][i] == 1:
                        masked_seq += base  # Keep original base
                    else:
                        masked_seq += "N"  # Mask with N
                        n_count += 1
                
                total_n_count += n_count
                n_percent = (n_count / len(sequence)) * 100 if len(sequence) > 0 else 0
                
                self.logger.debug(f"Sequence {seq_id}: {n_count}/{len(sequence)} bases masked to N ({n_percent:.2f}%)")
                
                # Create new record with masked sequence
                masked_record = SeqRecord(
                    Seq(masked_seq),
                    id=seq_id,
                    name=seq_id,
                    description=f"Masked sequence with {n_percent:.2f}% masked"
                )
                masked_records.append(masked_record)
            else:
                # No mask for this chromosome - include unmasked
                self.logger.debug(f"No mask created for {seq_id}, adding unmasked sequence")
                masked_record = SeqRecord(
                    Seq(sequence),
                    id=seq_id,
                    name=seq_id,
                    description="Unmasked sequence (no alignment data)"
                )
                masked_records.append(masked_record)
        
        # Return records and statistics
        stats = {
            "total_n_count": total_n_count,
            "total_bases": total_bases,
            "applied_regions": applied_regions,
            "total_regions": total_regions
        }
        
        return masked_records, stats
    
    def generate_coordinate_map(self, conserved_regions):
        """
        Generate a mapping between coordinates in reference and query genomes.
        
        Args:
            conserved_regions (dict): Conserved regions with mapping information
            
        Returns:
            dict: Coordinate mapping from reference to query
        """
        self.logger.debug("Generating coordinate map between reference and query genomes")
        
        coordinate_map = {}
        
        for chrom, regions in conserved_regions.items():
            coordinate_map[chrom] = {}
            
            for region in regions:
                qry_src = region['qry_src']
                qry_strand = region['qry_strand']
                
                # Create mapping for each position in this region
                ref_range = range(region['start'], region['end'] + 1)
                qry_range = range(region['qry_start'], region['qry_end'] + 1)
                
                # Handle reverse strand
                if qry_strand == '-':
                    qry_range = reversed(list(qry_range))
                
                # Map each reference position to query position
                for ref_pos, qry_pos in zip(ref_range, qry_range):
                    coordinate_map[chrom][ref_pos] = {
                        'qry_src': qry_src,
                        'qry_pos': qry_pos,
                        'qry_strand': qry_strand
                    }
        
        # Log statistics
        total_mappings = sum(len(pos_map) for pos_map in coordinate_map.values())
        self.logger.debug(f"Generated coordinate map with {total_mappings} position mappings across {len(coordinate_map)} chromosomes")
        
        return coordinate_map
    
    def map_second_variants_to_reference(self, second_variants, coordinate_map):
        """
        Map variant positions from second genome to reference genome coordinates.
        
        Args:
            second_variants (dict): Dictionary mapping second genomme chromosomes to sets of variant positions
            coordinate_map (dict): Coordinate mapping from reference to second genome
            
        Returns:
            dict: Dictionary mapping reference chromosomes to sets of mapped variant positions
            
        Raises:
            AlignmentError: If mapping fails
        """
        self.logger.debug("Mapping second genome variants to reference genome coordinates")
        
        try:
            # Dictionary to store mapped variants for each reference chromosome
            mapped_variants = {}
            
            # Create reverse mapping (from second genome to reference)
            reverse_map = self._create_reverse_coordinate_map(coordinate_map)
            
            # Count statistics
            total_second_variants = sum(len(pos_set) for pos_set in second_variants.values())
            mapped_count = 0
            
            # For each second genome chromosome and its variants
            for qry_chrom, positions in second_variants.items():
                if qry_chrom not in reverse_map:
                    self.logger.debug(f"No mapping found for second genome chromosome: {qry_chrom}")
                    continue
                    
                # For each variant position
                for pos in positions:
                    # Check if this position is in our mapping
                    if pos in reverse_map[qry_chrom]:
                        # Get the corresponding reference position
                        mapping = reverse_map[qry_chrom][pos]
                        ref_chrom = mapping['ref_chrom']
                        ref_pos = mapping['ref_pos']
                        
                        # Add to mapped variants
                        if ref_chrom not in mapped_variants:
                            mapped_variants[ref_chrom] = set()
                            
                        mapped_variants[ref_chrom].add(ref_pos)
                        mapped_count += 1
            
            # Log mapping statistics
            self.logger.debug(f"Successfully mapped {mapped_count}/{total_second_variants} second genome variants to reference genome")
            for ref_chrom, variants in mapped_variants.items():
                self.logger.debug(f"Reference chromosome {ref_chrom}: {len(variants)} mapped variants")
            
            return mapped_variants
            
        except Exception as e:
            self.logger.error(f"Error mapping second genome variants: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Failed to map second genome variants: {str(e)}")
    
    
    def _create_reverse_coordinate_map(self, coordinate_map):
        """
        Create a reverse mapping from query coordinates to reference coordinates.
        
        Args:
            coordinate_map (dict): Forward coordinate mapping (ref -> query)
            
        Returns:
            dict: Reverse coordinate mapping (query -> ref)
        """
        reverse_map = {}
        
        # Build the reverse map
        for ref_chrom, positions in coordinate_map.items():
            for ref_pos, mapping in positions.items():
                qry_src = mapping['qry_src']
                qry_pos = mapping['qry_pos']
                qry_strand = mapping['qry_strand']
                
                # Initialize the nested dictionary if needed
                if qry_src not in reverse_map:
                    reverse_map[qry_src] = {}
                    
                # Store the mapping: qry_pos -> ref_chrom, ref_pos
                reverse_map[qry_src][qry_pos] = {
                    'ref_chrom': ref_chrom,
                    'ref_pos': ref_pos,
                    'strand': qry_strand
                }
        
        return reverse_map
    
    def analyze_maf_file(self, maf_file):
        """
        Analyze a MAF file to identify sequence naming patterns and structure.
        
        Args:
            maf_file (str): Path to MAF file
            
        Returns:
            dict: Information about the MAF file structure
            
        Raises:
            AlignmentError: If analysis fails
        """
        self.logger.debug(f"Analyzing MAF file structure: {maf_file}")
        
        if not os.path.exists(maf_file):
            self.logger.error(f"MAF file not found: {maf_file}")
            raise AlignmentError(f"MAF file not found: {maf_file}")
            
        try:
            seq_ids = set()
            ref_seq_ids = set()
            query_seq_ids = set()
            alignment_count = 0
            
            with open(maf_file, 'r') as f:
                current_alignment = []
                
                for line in f:
                    line = line.strip()
                    
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue
                    
                    # New alignment block
                    if line.startswith('a'):
                        # Process previous alignment if it exists
                        if current_alignment:
                            alignment_count += 1
                            if len(current_alignment) >= 2:
                                # First sequence after the 'a' line is reference
                                parts = re.split(r'\s+', current_alignment[0].strip())
                                if len(parts) >= 7:
                                    ref_id = parts[1]
                                    ref_seq_ids.add(ref_id)
                                    seq_ids.add(ref_id)
                                
                                # Other sequences are queries
                                for seq_line in current_alignment[1:]:
                                    parts = re.split(r'\s+', seq_line.strip())
                                    if len(parts) >= 7:
                                        query_id = parts[1]
                                        query_seq_ids.add(query_id)
                                        seq_ids.add(query_id)
                            
                            current_alignment = []
                        
                    # Sequence lines
                    elif line.startswith('s'):
                        current_alignment.append(line)
                
                # Process the last alignment
                if current_alignment:
                    alignment_count += 1
                    if len(current_alignment) >= 2:
                        # First sequence after the 'a' line is reference
                        parts = re.split(r'\s+', current_alignment[0].strip())
                        if len(parts) >= 7:
                            ref_id = parts[1]
                            ref_seq_ids.add(ref_id)
                            seq_ids.add(ref_id)
                        
                        # Other sequences are queries
                        for seq_line in current_alignment[1:]:
                            parts = re.split(r'\s+', seq_line.strip())
                            if len(parts) >= 7:
                                query_id = parts[1]
                                query_seq_ids.add(query_id)
                                seq_ids.add(query_id)
            
            # Log results
            self.logger.debug(f"MAF file contains {alignment_count} alignment blocks")
            self.logger.debug(f"Found {len(seq_ids)} unique sequence IDs: {', '.join(list(seq_ids)[:5])}...")
            self.logger.debug(f"Reference sequences ({len(ref_seq_ids)}): {', '.join(list(ref_seq_ids)[:5])}...")
            self.logger.debug(f"Query sequences ({len(query_seq_ids)}): {', '.join(list(query_seq_ids)[:5])}...")
            
            return {
                'alignment_count': alignment_count,
                'seq_ids': sorted(seq_ids),
                'ref_seq_ids': sorted(ref_seq_ids),
                'query_seq_ids': sorted(query_seq_ids)
            }
            
        except Exception as e:
            self.logger.error(f"Error analyzing MAF file: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"MAF file analysis failed: {str(e)}")
    
    def _map_chromosomes(self, align_chroms, fasta_ids):
        """
        Create mapping between alignment chromosome names and FASTA IDs.
        
        Args:
            align_chroms (list): Chromosome names from alignments
            fasta_ids (list): Sequence IDs from FASTA file
            
        Returns:
            dict: Mapping from alignment chromosome names to FASTA IDs
        """
        chrom_map = {}
        fasta_ids = list(fasta_ids)
        
        # Try multiple strategies to match chromosome names
        for align_chrom in align_chroms:
            # Skip if already mapped
            if align_chrom in chrom_map:
                continue
                
            # Strategy 1: Direct match
            if align_chrom in fasta_ids:
                chrom_map[align_chrom] = align_chrom
                continue
                
            # Strategy 2: FASTA ID contains alignment chromosome
            matching_ids = [fid for fid in fasta_ids if align_chrom in fid]
            if matching_ids:
                chrom_map[align_chrom] = matching_ids[0]
                continue
                
            # Strategy 3: Try numeric chromosome matching
            if align_chrom.isdigit():
                # Look for FASTA IDs that start with this number
                matching_ids = [fid for fid in fasta_ids if fid.startswith(f"Chr{align_chrom}") or 
                                                            fid.startswith(f"chr{align_chrom}") or 
                                                            fid.startswith(f"chromosome{align_chrom}") or
                                                            f"_{align_chrom}" in fid]
                if matching_ids:
                    chrom_map[align_chrom] = matching_ids[0]
                    continue
            
            # Strategy 4: CP00268x.1 format matching
            if align_chrom.isdigit():
                cp_matches = [fid for fid in fasta_ids if f"CP00268{align_chrom}" in fid]
                if cp_matches:
                    chrom_map[align_chrom] = cp_matches[0]
                    continue
        
        # FALLBACK: If no matches found but we have same number of chromosomes,
        # use positional mapping (this is risky but better than nothing)
        if len(chrom_map) < len(align_chroms) and len(align_chroms) <= len(fasta_ids):
            unmapped_align = [ac for ac in align_chroms if ac not in chrom_map]
            unmapped_fasta = [fid for fid in fasta_ids if fid not in chrom_map.values()]
            
            if unmapped_align and unmapped_fasta:
                self.logger.debug("Using fallback positional mapping for unmapped chromosomes")
                
                # Try to sort both lists in a meaningful way
                if all(c.isdigit() for c in unmapped_align):
                    unmapped_align.sort(key=int)
                
                # Sort FASTA IDs by numeric part if they have the CP00268x.1 format
                def extract_num(fid):
                    match = re.search(r'CP00268(\d)\.1', fid)
                    return int(match.group(1)) if match else 999
                    
                unmapped_fasta.sort(key=extract_num)
                
                # Create mappings
                for i, (ac, fid) in enumerate(zip(unmapped_align, unmapped_fasta)):
                    chrom_map[ac] = fid
        
        return chrom_map