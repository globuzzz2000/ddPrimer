#!/usr/bin/env python3
"""
SNP processing using VCF normalization approach with intelligent chromosome mapping.

This module implements the industry-standard approach for handling VCF variants:
1. Normalize VCF using bcftools to eliminate coordinate ambiguity
2. Apply sequence modifications in a single pass
3. Handle both masking and substitution operations cleanly
4. Intelligent chromosome name mapping between VCF and FASTA files
5. Scalable processing using bcftools streaming for large files

Requires bcftools to be installed and available in PATH.
"""

import subprocess
import tempfile
import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import from your exceptions module
from ..config import VCFNormalizationError

logger = logging.getLogger(__name__)


class SNPMaskingProcessor:
    """
    SNP processor using VCF normalization approach with intelligent chromosome mapping.
    
    This class handles SNP processing using the industry-standard VCF normalization
    workflow to eliminate coordinate drift issues, combined with intelligent
    chromosome name mapping between VCF and FASTA files. Uses scalable bcftools
    operations for efficient processing of large VCF files.
    """
    
    def __init__(self, reference_fasta: str):
        """
        Initialize SNP processor.
        
        Args:
            reference_fasta: Path to reference FASTA file
            
        Raises:
            FileNotFoundError: If reference FASTA not found
            VCFNormalizationError: If bcftools not available
        """
        self.reference_fasta = Path(reference_fasta)
        if not self.reference_fasta.exists():
            raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")
        
        # Check if bcftools is available
        try:
            subprocess.run(['bcftools', '--version'], 
                         capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise VCFNormalizationError(
                "bcftools not found. Please install bcftools and ensure it's in PATH"
            )
        
        # Initialize chromosome mapping (will be built when needed)
        self.chromosome_mapping = None
    
    def build_chromosome_mapping(self, vcf_path: str):
        """
        Build chromosome mapping between VCF and FASTA using ChromosomeMapper.
        
        Args:
            vcf_path: Path to VCF file for chromosome analysis
            
        Returns:
            Dictionary mapping VCF chromosomes to FASTA sequence names
        """
        if self.chromosome_mapping is not None:
            return self.chromosome_mapping
        
        try:
            # Import the ChromosomeMapper from your existing module
            from ..utils import ChromosomeMapper
            
            mapper = ChromosomeMapper()
            
            # Check compatibility and get mapping
            analysis = mapper.check_chromosome_compatibility(
                vcf_path, str(self.reference_fasta)
            )
            
            if analysis['compatible'] and not analysis['needs_mapping']:
                # Direct mapping (chromosome names match exactly)
                logger.debug("VCF and FASTA chromosome names are compatible")
                self.chromosome_mapping = {chrom: chrom for chrom in analysis['exact_matches']}
            
            elif analysis['needs_mapping'] and analysis['suggested_mapping']:
                # Use suggested automatic mapping
                self.chromosome_mapping = analysis['suggested_mapping']
                logger.info(f"Using automatic chromosome mapping for {len(self.chromosome_mapping)} chromosomes")
                
                if logger.isEnabledFor(logging.DEBUG):
                    for vcf_chr, fasta_chr in self.chromosome_mapping.items():
                        logger.debug(f"  VCF '{vcf_chr}' -> FASTA '{fasta_chr}'")
            
            else:
                logger.warning("No chromosome mapping could be established")
                self.chromosome_mapping = {}
            
            return self.chromosome_mapping
            
        except Exception as e:
            logger.error(f"Error building chromosome mapping: {e}")
            # Fall back to empty mapping
            self.chromosome_mapping = {}
            return self.chromosome_mapping
    
    def map_chromosome_name(self, vcf_chromosome: str) -> Optional[str]:
        """
        Map VCF chromosome name to FASTA sequence name.
        
        Args:
            vcf_chromosome: Chromosome name from VCF file
            
        Returns:
            Corresponding FASTA sequence name, or None if no mapping exists
        """
        if self.chromosome_mapping is None:
            logger.warning("Chromosome mapping not initialized")
            return vcf_chromosome  # Fall back to original name
        
        return self.chromosome_mapping.get(vcf_chromosome, None)
    
    def _create_chromosome_mapping_file(self) -> str:
        """
        Create a chromosome mapping file for bcftools --rename-chrs.
        
        Format is "old_name new_name" separated by tab, one per line.
        
        Returns:
            Path to the mapping file
            
        Raises:
            VCFNormalizationError: If mapping file creation fails
        """
        try:
            # Create temporary mapping file
            fd, mapping_file = tempfile.mkstemp(suffix='.txt', prefix='chr_mapping_')
            
            with os.fdopen(fd, 'w') as f:
                for vcf_chr, fasta_chr in self.chromosome_mapping.items():
                    f.write(f"{vcf_chr}\t{fasta_chr}\n")
                    logger.debug(f"Mapping: {vcf_chr} -> {fasta_chr}")
            
            logger.debug(f"Created chromosome mapping file with {len(self.chromosome_mapping)} mappings")
            return mapping_file
            
        except Exception as e:
            error_msg = f"Failed to create chromosome mapping file: {str(e)}"
            logger.error(error_msg)
            raise VCFNormalizationError(error_msg) from e
    
    def normalize_vcf(self, vcf_path: str, output_path: Optional[str] = None) -> str:
        """
        Normalize VCF file using bcftools with scalable chromosome mapping.
        
        Uses bcftools annotate --rename-chrs for efficient chromosome mapping
        and bcftools norm for variant normalization. This approach is scalable
        for large VCF files but requires temporary disk space.
        
        Args:
            vcf_path: Path to input VCF file
            output_path: Path for output normalized VCF (if None, uses temp file)
            
        Returns:
            Path to normalized VCF file
            
        Raises:
            VCFNormalizationError: If normalization fails
        """
        if output_path is None:
            # Create temporary file for normalized VCF
            fd, output_path = tempfile.mkstemp(suffix='.vcf.gz', prefix='normalized_')
            os.close(fd)
        
        # Build chromosome mapping if not already done
        if self.chromosome_mapping is None:
            self.build_chromosome_mapping(vcf_path)
        
        try:
            logger.info(f"Normalizing VCF with scalable approach: {vcf_path}")
            
            # Step 1: Create chromosome mapping file for bcftools
            mapping_file = None
            intermediate_vcf = None
            
            if self.chromosome_mapping:
                mapping_file = self._create_chromosome_mapping_file()
                logger.debug(f"Created chromosome mapping file: {mapping_file}")
                
                # Create intermediate file with renamed chromosomes
                fd, intermediate_vcf = tempfile.mkstemp(suffix='.vcf.gz', prefix='renamed_')
                os.close(fd)
                
                # Rename chromosomes using bcftools (memory-efficient)
                rename_cmd = [
                    'bcftools', 'annotate',
                    '--rename-chrs', mapping_file,
                    '--output-type', 'z',  # Compressed output
                    '--output', intermediate_vcf,
                    vcf_path
                ]
                
                logger.debug(f"Running chromosome rename: {' '.join(rename_cmd)}")
                result = subprocess.run(rename_cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    raise VCFNormalizationError(
                        f"Chromosome renaming failed:\n"
                        f"STDOUT: {result.stdout}\n"
                        f"STDERR: {result.stderr}"
                    )
                
                input_for_norm = intermediate_vcf
                logger.debug("Chromosome renaming completed successfully")
            else:
                logger.warning("No chromosome mapping available - using original VCF")
                input_for_norm = vcf_path
            
            # Step 2: Normalize the VCF (now with correct chromosome names)
            norm_cmd = [
                'bcftools', 'norm',
                '-m-both',  # Split multiallelic sites into biallelic records
                '-f', str(self.reference_fasta),  # Reference FASTA
                '--output-type', 'z',  # Compressed output
                '--output', output_path,
                input_for_norm
            ]
            
            logger.debug(f"Running normalization: {' '.join(norm_cmd)}")
            result = subprocess.run(norm_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise VCFNormalizationError(
                    f"bcftools normalization failed:\n"
                    f"STDOUT: {result.stdout}\n"
                    f"STDERR: {result.stderr}"
                )
            
            # Index the output file for efficient access
            try:
                index_cmd = ['bcftools', 'index', output_path]
                subprocess.run(index_cmd, capture_output=True)
                logger.debug(f"Indexed normalized VCF: {output_path}")
            except:
                logger.debug("Could not index VCF (non-critical)")
            
            logger.info(f"VCF normalized successfully: {output_path}")
            return output_path
            
        except subprocess.CalledProcessError as e:
            raise VCFNormalizationError(f"bcftools execution failed: {e}")
        except Exception as e:
            # Clean up output file if we created it and there was an error
            if output_path and os.path.exists(output_path):
                try:
                    os.unlink(output_path)
                except:
                    pass
            raise VCFNormalizationError(f"Normalization failed: {e}")
        finally:
            # Clean up intermediate files
            if intermediate_vcf and os.path.exists(intermediate_vcf):
                try:
                    os.unlink(intermediate_vcf)
                    # Also clean up index if it exists
                    if os.path.exists(intermediate_vcf + '.csi'):
                        os.unlink(intermediate_vcf + '.csi')
                    logger.debug(f"Cleaned up intermediate VCF: {intermediate_vcf}")
                except:
                    pass
            
            if mapping_file and os.path.exists(mapping_file):
                try:
                    os.unlink(mapping_file)
                    logger.debug(f"Cleaned up mapping file: {mapping_file}")
                except:
                    pass
    
    def parse_normalized_vcf(self, vcf_path: str, 
                           chromosome: str) -> List[Dict]:
        """
        Parse normalized VCF file and extract variants for a specific chromosome.
        
        Since we now remap chromosomes before normalization, the normalized VCF
        already contains the correct FASTA chromosome names. Uses bcftools query
        for memory-efficient extraction of variants.
        
        Args:
            vcf_path: Path to normalized VCF file
            chromosome: FASTA chromosome name to extract variants for
            
        Returns:
            List of variant dictionaries with keys: pos, ref, alt, qual, info
        """
        logger.debug(f"Parsing normalized VCF for chromosome: {chromosome}")
        
        try:
            # Use bcftools query for memory-efficient variant extraction
            query_cmd = [
                'bcftools', 'query',
                '-r', chromosome,  # Only this chromosome
                '-f', '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO\n',
                vcf_path
            ]
            
            logger.debug(f"Querying variants: {' '.join(query_cmd)}")
            result = subprocess.run(query_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.warning(f"bcftools query failed for {chromosome}: {result.stderr}")
                return []
            
            variants = []
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                    
                fields = line.split('\t')
                if len(fields) < 6:
                    continue
                
                chrom, pos_str, ref, alt, qual_str, info = fields
                
                # Skip if alt is '.' (no alternative)
                if alt == '.':
                    continue
                
                # Parse position and quality
                try:
                    pos = int(pos_str)
                except ValueError:
                    continue
                
                qual_value = None
                if qual_str != "." and qual_str != "":
                    try:
                        qual_value = float(qual_str)
                    except ValueError:
                        qual_value = None
                
                variants.append({
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'qual': qual_value,
                    'info': info
                })
            
            logger.info(f"Parsed {len(variants)} variants for chromosome {chromosome}")
            return variants
            
        except Exception as e:
            logger.error(f"Error parsing VCF file {vcf_path} for chromosome {chromosome}: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return []
    
    def apply_variants_to_sequence(self, 
                                 sequence: str, 
                                 variants: List[Dict],
                                 snp_allele_frequency_threshold: Optional[float] = None,
                                 snp_quality_threshold: Optional[float] = None,
                                 snp_flanking_mask_size: int = 0,
                                 snp_use_soft_masking: bool = False) -> str:
        """
        Apply variants to sequence using single-pass approach.
        
        Since the VCF is normalized, we can apply all variants in a single pass
        without coordinate drift issues.
        
        Args:
            sequence: Input DNA sequence
            variants: List of normalized variants
            snp_allele_frequency_threshold: Minimum AF to apply masking
            snp_quality_threshold: Minimum QUAL to include variants
            snp_flanking_mask_size: Bases to mask around each variant
            snp_use_soft_masking: Use lowercase instead of 'N'
            
        Returns:
            Modified sequence with variants applied
        """
        if not variants:
            return sequence
        
        # Convert to mutable list for in-place modifications
        seq_list = list(sequence)
        
        # Sort variants by position (descending) to avoid coordinate drift
        sorted_variants = sorted(variants, key=lambda x: x['pos'], reverse=True)
        
        logger.info(f"Applying {len(sorted_variants)} variants to sequence")
        
        masked_count = 0
        substituted_count = 0
        
        for variant in sorted_variants:
            pos = variant['pos'] - 1  # Convert to 0-based indexing
            ref = variant['ref']
            alt = variant['alt']
            qual = variant['qual']
            info = variant['info']
            
            # Apply quality threshold
            if snp_quality_threshold is not None and qual is not None:
                if qual < snp_quality_threshold:
                    logger.debug(f"Skipping variant at {pos+1}: QUAL {qual} < {snp_quality_threshold}")
                    continue
            
            # Extract allele frequency
            af = self._extract_allele_frequency(info)
            
            # Determine if we should mask or substitute
            should_mask = (
                snp_allele_frequency_threshold is not None and
                af is not None and
                af >= snp_allele_frequency_threshold
            )
            
            # Validate reference sequence matches
            if pos + len(ref) <= len(seq_list):
                ref_in_seq = ''.join(seq_list[pos:pos+len(ref)])
                if ref_in_seq.upper() != ref.upper():
                    logger.warning(
                        f"Reference mismatch at position {pos+1}: "
                        f"expected {ref}, found {ref_in_seq}"
                    )
                    continue
            else:
                logger.warning(f"Variant at position {pos+1} extends beyond sequence")
                continue
            
            if should_mask:
                # Apply masking
                self._apply_masking(
                    seq_list, pos, len(ref), 
                    snp_flanking_mask_size, snp_use_soft_masking
                )
                masked_count += 1
                logger.debug(f"Masked variant at {pos+1}: {ref} (AF={af})")
            else:
                # Apply substitution
                self._apply_substitution(seq_list, pos, ref, alt)
                substituted_count += 1
                logger.debug(f"Substituted variant at {pos+1}: {ref}->{alt}")
        
        logger.info(f"Applied {masked_count} maskings and {substituted_count} substitutions")
        return ''.join(seq_list)
    
    def _extract_allele_frequency(self, info: str) -> Optional[float]:
        """Extract allele frequency from INFO field."""
        # Look for common AF tags
        af_tags = ['AF=', 'AC=', 'MAF=']
        for tag in af_tags:
            for field in info.split(';'):
                if field.startswith(tag):
                    try:
                        value = field.split('=')[1]
                        # Handle comma-separated values (take first)
                        if ',' in value:
                            value = value.split(',')[0]
                        return float(value)
                    except (ValueError, IndexError):
                        continue
        return None
    
    def _apply_masking(self, seq_list: List[str], pos: int, length: int,
                      flanking_size: int, use_soft_masking: bool):
        """Apply masking to sequence at specified position."""
        mask_char = 'n' if use_soft_masking else 'N'
        
        # Calculate masking range
        start = max(0, pos - flanking_size)
        end = min(len(seq_list), pos + length + flanking_size)
        
        # Apply masking
        for i in range(start, end):
            seq_list[i] = mask_char
    
    def _apply_substitution(self, seq_list: List[str], pos: int, ref: str, alt: str):
        """Apply substitution to sequence at specified position."""
        # Handle different variant types
        if len(ref) == len(alt):
            # SNP or MNP
            for i, base in enumerate(alt):
                if pos + i < len(seq_list):
                    seq_list[pos + i] = base
        elif len(ref) > len(alt):
            # Deletion
            # Replace ref with alt, then remove extra bases
            for i, base in enumerate(alt):
                if pos + i < len(seq_list):
                    seq_list[pos + i] = base
            # Remove the extra bases
            del seq_list[pos + len(alt):pos + len(ref)]
        else:
            # Insertion
            # Replace ref with alt (which is longer)
            for i, base in enumerate(ref):
                if pos + i < len(seq_list):
                    seq_list[pos + i] = base
            # Insert additional bases
            insertion = alt[len(ref):]
            for i, base in enumerate(insertion):
                seq_list.insert(pos + len(ref) + i, base)
    
    def process_sequence_with_vcf(self, 
                                sequence: Union[str, SeqRecord],
                                vcf_path: str,
                                chromosome: str,
                                **kwargs) -> str:
        """
        Complete workflow: normalize VCF and apply variants to sequence.
        
        This method uses the scalable bcftools-based approach for processing
        large VCF files efficiently.
        
        Args:
            sequence: Input sequence (string or SeqRecord)
            vcf_path: Path to VCF file
            chromosome: FASTA chromosome identifier (will be mapped to VCF chromosome)
            **kwargs: Additional arguments for apply_variants_to_sequence
            
        Returns:
            Modified sequence string
        """
        # Extract sequence string if SeqRecord
        if isinstance(sequence, SeqRecord):
            seq_str = str(sequence.seq)
        else:
            seq_str = sequence
        
        # Normalize VCF using scalable approach
        normalized_vcf = None
        try:
            normalized_vcf = self.normalize_vcf(vcf_path)
            
            # Parse variants for this chromosome using bcftools query
            variants = self.parse_normalized_vcf(normalized_vcf, chromosome)
            
            # Apply variants
            modified_sequence = self.apply_variants_to_sequence(
                seq_str, variants, **kwargs
            )
            
            return modified_sequence
            
        finally:
            # Clean up temporary normalized VCF
            if normalized_vcf and os.path.exists(normalized_vcf):
                try:
                    os.unlink(normalized_vcf)
                    # Also clean up index if it exists
                    if os.path.exists(normalized_vcf + '.csi'):
                        os.unlink(normalized_vcf + '.csi')
                    logger.debug(f"Cleaned up normalized VCF: {normalized_vcf}")
                except:
                    pass


def process_fasta_with_vcf(fasta_path: str, 
                          vcf_path: str, 
                          output_path: str,
                          reference_fasta: str,
                          **kwargs):
    """
    Process FASTA file with VCF variants using scalable normalization approach.
    
    Args:
        fasta_path: Input FASTA file
        vcf_path: Input VCF file
        output_path: Output FASTA file
        reference_fasta: Reference FASTA for VCF normalization
        **kwargs: Additional arguments for variant processing
    """
    from Bio import SeqIO
    
    processor = SNPMaskingProcessor(reference_fasta)
    
    with open(output_path, 'w') as output_handle:
        for record in SeqIO.parse(fasta_path, 'fasta'):
            logger.info(f"Processing sequence: {record.id}")
            
            # Process sequence with VCF
            modified_seq = processor.process_sequence_with_vcf(
                record, vcf_path, record.id, **kwargs
            )
            
            # Create new record with modified sequence
            new_record = SeqRecord(
                Seq(modified_seq),
                id=record.id,
                description=record.description + " [VCF_processed]"
            )
            
            SeqIO.write(new_record, output_handle, 'fasta')
    
    logger.info(f"Processed sequences written to: {output_path}")