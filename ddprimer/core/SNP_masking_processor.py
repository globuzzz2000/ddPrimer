import os
import logging
import subprocess
import gzip
from Bio import SeqIO
from tqdm import tqdm

# Import package modules
from ..config import Config, FileError, SequenceProcessingError

# Set up logging
logger = logging.getLogger("ddPrimer.snp_masking_processor")

class SNPMaskingProcessor:
    """Handles masking of SNPs in sequences to prepare them for primer design."""
    
    def __init__(self):
        """Initialize SNP masking processor."""
        logger.debug("Initialized SNPMaskingProcessor")
    
    def get_variant_positions(self, vcf_file, chromosome=None):
        """
        Extract variant positions from the VCF file using bcftools.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        logger.debug(f"\nFetching variant positions from {vcf_file}")
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
        
        # Base command
        command = f'bcftools query -f "%CHROM\\t%POS\\n" "{vcf_file}"'
        
        # Add chromosome filter if specified
        if chromosome:
            command += f' -r "{chromosome}"'
            logger.debug(f"Filtering variants for chromosome: {chromosome}")
        
        # Run command
        logger.debug(f"Running command: {command}")
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"Error running bcftools: {result.stderr}")
                raise FileError(f"bcftools query failed: {result.stderr}")
            
            # Parse results
            variants = {}
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                    
                parts = line.split("\t")
                if len(parts) == 2:
                    chrom, pos = parts
                    try:
                        pos = int(pos)
                        
                        if chrom not in variants:
                            variants[chrom] = set()
                            
                        variants[chrom].add(pos)
                    except ValueError:
                        logger.warning(f"Invalid position value in VCF: {pos}")
            
            total_variants = sum(len(positions) for positions in variants.values())
            logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes")
            return variants
            
        except subprocess.SubprocessError as e:
            logger.error(f"Failed to run bcftools: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Failed to run bcftools: {e}")
    
    def get_filtered_variants(self, vcf_file, chromosome=None, min_af=None, min_qual=None):
        """
        Extract variant positions from VCF file with allele frequency and quality filtering.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            min_af (float, optional): Minimum allele frequency threshold (0.0-1.0)
            min_qual (float, optional): Minimum QUAL score threshold
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        logger.debug(f"Fetching filtered variants from {vcf_file}")
        logger.debug(f"Filters: min_af={min_af}, min_qual={min_qual}, chromosome={chromosome}")
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
        
        # If no filters specified, use basic extraction
        if min_af is None and min_qual is None:
            logger.debug("No filters specified, using basic variant extraction")
            return self.get_variant_positions(vcf_file, chromosome)
        
        # First, try bcftools approach with filtering
        try:
            return self._extract_with_bcftools(vcf_file, chromosome, min_af, min_qual)
        except Exception as e:
            logger.warning(f"bcftools filtering failed: {str(e)}")
            logger.info("Falling back to manual VCF parsing with filtering")
            # Use manual parsing but with filtering applied
            return self._extract_variants_manually_global(vcf_file, chromosome, min_af, min_qual)
    
    def _extract_with_bcftools(self, vcf_file, chromosome=None, min_af=None, min_qual=None):
        """
        Extract variants using bcftools with proper filtering.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            min_af (float, optional): Minimum allele frequency threshold
            min_qual (float, optional): Minimum QUAL score threshold
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
        """
        # Build filters carefully
        filters = []
        
        # Quality filter
        if min_qual is not None:
            filters.append(f"QUAL>={min_qual}")
        
        # For AF filtering, we'll include it in the filter but handle errors gracefully
        if min_af is not None:
            filters.append(f"INFO/AF>={min_af}")
        
        # Build command
        if filters:
            filter_string = " && ".join(filters)
            command = f'bcftools query -i "{filter_string}" -f "%CHROM\\t%POS\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        else:
            command = f'bcftools query -f "%CHROM\\t%POS\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        
        # Add chromosome filter if specified
        if chromosome:
            command += f' -r "{chromosome}"'
            logger.debug(f"Filtering variants for chromosome: {chromosome}")
        
        logger.debug(f"Running bcftools command: {command}")
        
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                # If AF filtering fails, try without AF filter but with manual AF filtering
                if min_af is not None and "INFO/AF" in result.stderr:
                    logger.warning("AF field might not exist in VCF INFO, trying alternative approach")
                    return self._extract_with_bcftools_fallback(vcf_file, chromosome, min_af, min_qual)
                else:
                    raise subprocess.SubprocessError(f"bcftools failed: {result.stderr}")
            
            # Parse results and apply manual filtering as backup
            variants = {}
            total_processed = 0
            total_kept = 0
            
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                    
                parts = line.split("\t")
                if len(parts) >= 2:
                    chrom = parts[0]
                    try:
                        pos = int(parts[1])
                        qual = float(parts[2]) if len(parts) > 2 and parts[2] != "." else None
                        af = None
                        if len(parts) > 3 and parts[3] not in [".", ""]:
                            try:
                                af_values = [float(x.replace(",", ".")) for x in parts[3].split(",")]
                                af = sum(af_values)
                            except ValueError:
                                logger.warning(f"Invalid AF value: {parts[3]}")
                        
                        total_processed += 1
                        
                        # Apply manual filters as backup (bcftools should have done this, but double-check)
                        passes_filter = True
                        
                        if min_qual is not None and qual is not None and qual < min_qual:
                            passes_filter = False
                        
                        if min_af is not None and af is not None and af < min_af:
                            passes_filter = False
                        
                        if passes_filter:
                            if chrom not in variants:
                                variants[chrom] = set()
                            variants[chrom].add(pos)
                            total_kept += 1
                            
                    except ValueError:
                        logger.warning(f"BInvalid position value in VCF: {parts[1]}")
                        continue
            
            logger.info(f"bcftools processing: {total_processed} variants examined, {total_kept} kept after filtering")
            if min_af is not None:
                logger.debug(f"AF threshold applied: >= {min_af}")
            if min_qual is not None:
                logger.debug(f"QUAL threshold applied: >= {min_qual}")
            
            return variants
            
        except subprocess.SubprocessError:
            # Re-raise subprocess errors
            raise
        except Exception as e:
            logger.warning(f"Error in bcftools processing: {str(e)}")
            raise
    
    def _extract_with_bcftools_fallback(self, vcf_file, chromosome=None, min_af=None, min_qual=None):
        """
        Fallback bcftools approach when AF field filtering fails.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            min_af (float, optional): Minimum allele frequency threshold
            min_qual (float, optional): Minimum QUAL score threshold
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
        """
        # Try with only QUAL filter first, then manually filter AF
        filters = []
        if min_qual is not None:
            filters.append(f"QUAL>={min_qual}")
        
        if filters:
            filter_string = " && ".join(filters)
            command = f'bcftools query -i "{filter_string}" -f "%CHROM\\t%POS\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        else:
            command = f'bcftools query -f "%CHROM\\t%POS\\t%QUAL\\t%INFO/AF\\n" "{vcf_file}"'
        
        if chromosome:
            command += f' -r "{chromosome}"'
        
        logger.debug(f"Fallback bcftools command: {command}")
        
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise subprocess.SubprocessError(f"Fallback bcftools failed: {result.stderr}")
        
        # Parse and manually filter
        variants = {}
        total_processed = 0
        total_kept = 0
        
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
                
            parts = line.split("\t")
            if len(parts) >= 2:
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                    qual = float(parts[2]) if len(parts) > 2 and parts[2] != "." else None
                    af_str = parts[3] if len(parts) > 3 else ""
                    
                    total_processed += 1
                    
                    # Parse AF from INFO field if needed
                    if min_af is not None:
                        af = self._parse_af_from_string(af_str)
                        if af is None or af < min_af:
                            continue
                    
                    # QUAL should already be filtered by bcftools, but double-check
                    if min_qual is not None and qual is not None and qual < min_qual:
                        continue
                    
                    if chrom not in variants:
                        variants[chrom] = set()
                    variants[chrom].add(pos)
                    total_kept += 1
                    
                except ValueError:
                    logger.warning(f"Invalid position value: {parts[1]}")
                    continue
        
        logger.info(f"bcftools fallback: {total_processed} variants examined, {total_kept} kept after filtering")
        return variants
    
    def _extract_variants_manually_global(self, vcf_file, chromosome=None, min_af=None, min_qual=None):
        """
        Extract variants using manual VCF parsing with filtering applied globally.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str, optional): Specific chromosome to filter
            min_af (float, optional): Minimum allele frequency threshold
            min_qual (float, optional): Minimum QUAL score threshold
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
        """
        variants = {}
        total_processed = 0
        total_kept = 0
        
        try:
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with open_func(vcf_file, mode) as vcf:
                for line in vcf:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue
                    
                    vcf_chrom = fields[0]
                    if chromosome and vcf_chrom != chromosome:
                        continue
                    
                    try:
                        pos = int(fields[1])
                        total_processed += 1
                        
                        # Apply quality filter
                        if min_qual is not None:
                            try:
                                qual = float(fields[5])
                                if qual < min_qual:
                                    continue
                            except (ValueError, IndexError):
                                # Skip if QUAL cannot be parsed and filter is required
                                continue
                        
                        # Apply allele frequency filter
                        if min_af is not None:
                            info_field = fields[7] if len(fields) > 7 else ""
                            af_value = self._parse_af_from_info(info_field)
                            if af_value is None or af_value < min_af:
                                continue
                        
                        if vcf_chrom not in variants:
                            variants[vcf_chrom] = set()
                        variants[vcf_chrom].add(pos)
                        total_kept += 1
                        
                    except ValueError:
                        logger.warning(f"Invalid position value: {fields[1]}")
                        continue
            
            logger.info(f"Manual parsing: {total_processed} variants examined, {total_kept} kept after filtering")
            if min_af is not None:
                logger.info(f"Manual AF filtering applied: >= {min_af}")
            if min_qual is not None:
                logger.info(f"Manual QUAL filtering applied: >= {min_qual}")
            
        except Exception as e:
            logger.error(f"Error in manual VCF parsing: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        
        return variants
    
    def prepare_vcf_file(self, vcf_file):
        """
        Prepare a VCF file for region queries by ensuring it's compressed with bgzip and indexed with tabix.
        
        Args:
            vcf_file (str): Path to VCF file (.vcf or .vcf.gz)
            
        Returns:
            str: Path to the prepared VCF file
            
        Raises:
            FileError: If VCF file cannot be prepared
        """
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
        
        # Check if bgzip and tabix are available
        try:
            subprocess.run(["bgzip", "--version"], capture_output=True, check=False)
            subprocess.run(["tabix", "--version"], capture_output=True, check=False)
        except FileNotFoundError:
            logger.warning("bgzip or tabix not found. Please install them for optimal performance.")
            logger.warning("Falling back to slower VCF parsing method.")
            return vcf_file
        
        # Determine if we need to compress the file
        if vcf_file.endswith('.vcf'):
            compressed_path = vcf_file + '.gz'
            
            # Check if the compressed file already exists and is newer than the original
            if os.path.exists(compressed_path) and os.path.getmtime(compressed_path) > os.path.getmtime(vcf_file):
                logger.debug(f"Using existing compressed VCF: {compressed_path}")
            else:
                logger.info(f"Compressing VCF file with bgzip...")
                try:
                    result = subprocess.run(
                        ["bgzip", "-c", vcf_file], 
                        stdout=open(compressed_path, 'wb'),
                        stderr=subprocess.PIPE,
                        check=True
                    )
                    logger.debug("VCF compression completed successfully")
                except subprocess.CalledProcessError as e:
                    logger.warning(f"Failed to compress VCF file: {e.stderr.decode('utf-8')}")
                    return vcf_file
            
            # Use the compressed file for subsequent operations
            vcf_file = compressed_path
        
        # Ensure the compressed file is indexed
        index_path = vcf_file + '.tbi'
        if not os.path.exists(index_path) or os.path.getmtime(index_path) < os.path.getmtime(vcf_file):
            logger.info(f"Indexing VCF file with tabix...")
            try:
                result = subprocess.run(
                    ["tabix", "-p", "vcf", vcf_file],
                    stderr=subprocess.PIPE,
                    check=True
                )
                logger.debug("VCF indexing completed successfully")
            except subprocess.CalledProcessError as e:
                logger.warning(f"Failed to index VCF file: {e.stderr.decode('utf-8')}")
                # Even if indexing fails, return the compressed file as we can still use it
        
        return vcf_file
    
    def get_region_variants(self, vcf_file, chromosome, start_pos, end_pos, min_af=None, min_qual=None):
        """
        Extract variant positions from the VCF file for a specific genomic region with filtering.
        
        Args:
            vcf_file (str): Path to VCF file
            chromosome (str): Chromosome name
            start_pos (int): Start position of the region
            end_pos (int): End position of the region
            min_af (float, optional): Minimum allele frequency threshold
            min_qual (float, optional): Minimum QUAL score threshold
            
        Returns:
            set: Set of variant positions within the specified region
            
        Raises:
            SequenceProcessingError: If region variants cannot be extracted
        """
        logger.debug(f"Fetching variants for region {chromosome}:{start_pos}-{end_pos}")
        logger.debug(f"Filters: min_af={min_af}, min_qual={min_qual}")
        
        # Prepare the VCF file (compress and index if needed)
        try:
            prepared_vcf = self.prepare_vcf_file(vcf_file)
            
            # Build filters for bcftools
            filters = []
            if min_qual is not None:
                filters.append(f"QUAL>={min_qual}")
            if min_af is not None:
                filters.append(f"INFO/AF>={min_af}")
            
            # Build bcftools command
            if filters:
                filter_string = " && ".join(filters)
                command = f'bcftools query -i "{filter_string}" -f "%POS\\n" -r "{chromosome}:{start_pos}-{end_pos}" "{prepared_vcf}"'
            else:
                command = f'bcftools query -f "%POS\\n" -r "{chromosome}:{start_pos}-{end_pos}" "{prepared_vcf}"'
            
            logger.debug(f"Running command: {command}")
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                # Parse results from bcftools
                variants = set()
                for line in result.stdout.strip().split("\n"):
                    if line:  # Skip empty lines
                        try:
                            pos = int(line.strip())
                            variants.add(pos)
                        except ValueError:
                            logger.warning(f"Invalid position value: {line}")
                            continue
                
                logger.debug(f"Extracted {len(variants)} variants in region {chromosome}:{start_pos}-{end_pos}")
                return variants
            else:
                logger.warning(f"bcftools query failed: {result.stderr}")
                # Fall back to manual parsing
        except Exception as e:
            logger.warning(f"Error running bcftools: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            # Fall back to manual parsing
        
        # Manual parsing approach for when bcftools fails
        logger.debug("Falling back to manual VCF parsing")
        return self._extract_variants_manually(prepared_vcf, chromosome, start_pos, end_pos, min_af, min_qual)
    
    def _extract_variants_manually(self, vcf_file, chrom, start, end, min_af=None, min_qual=None):
        """
        Extract variants from a VCF file for a specific region using manual parsing.
        This is a fallback method when bcftools fails.
        
        Args:
            vcf_file (str): Path to VCF file
            chrom (str): Chromosome name
            start (int): Start position
            end (int): End position
            min_af (float, optional): Minimum allele frequency threshold
            min_qual (float, optional): Minimum QUAL score threshold
            
        Returns:
            set: Set of variant positions in the region
        """
        variants = set()
        
        try:
            # Determine if we need to handle gzip compression
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with open_func(vcf_file, mode) as vcf:
                for line in vcf:
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Parse VCF data line
                    fields = line.strip().split('\t')
                    if len(fields) < 8:  # Minimum valid VCF line (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
                        continue
                    
                    # Check chromosome match
                    vcf_chrom = fields[0]
                    if vcf_chrom != chrom:
                        continue
                    
                    # Check position is in our region
                    try:
                        pos = int(fields[1])
                        if not (start <= pos <= end):
                            continue
                        
                        # Apply quality filter
                        if min_qual is not None:
                            try:
                                qual = float(fields[5])
                                if qual < min_qual:
                                    continue
                            except (ValueError, IndexError):
                                # Skip if QUAL cannot be parsed
                                if min_qual is not None:
                                    continue
                        
                        # Apply allele frequency filter
                        if min_af is not None:
                            info_field = fields[7] if len(fields) > 7 else ""
                            af_value = self._parse_af_from_info(info_field)
                            if af_value is None or af_value < min_af:
                                continue
                        
                        variants.add(pos)
                        
                    except ValueError:
                        logger.warning(f"Invalid position value in VCF: {fields[1]}")
                        continue
            
            logger.debug(f"Manually extracted {len(variants)} variants in region {chrom}:{start}-{end}")
            
        except Exception as e:
            logger.error(f"Error in manual VCF parsing for region {chrom}:{start}-{end}: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        
        return variants
    
    def _parse_af_from_info(self, info_field):
        """
        Parse allele frequency from INFO field.
        
        Args:
            info_field (str): INFO field from VCF line
            
        Returns:
            float or None: Allele frequency if found, None otherwise
        """
        try:
            # Look for AF= in the INFO field
            for info_item in info_field.split(';'):
                if info_item.startswith('AF='):
                    af_str = info_item[3:]  # Remove 'AF=' prefix
                    # Handle multiple values (take the first one)
                    if ',' in af_str:
                        af_str = af_str.split(',')[0]
                    return float(af_str)
        except (ValueError, AttributeError):
            pass
        
        return None
    
    def _parse_af_from_string(self, af_str):
        """
        Parse allele frequency from a string that might be from INFO/AF field.
        
        Args:
            af_str (str): String containing AF value
            
        Returns:
            float or None: Allele frequency if found, None otherwise
        """
        try:
            if af_str and af_str != "." and af_str != "":
                # Handle comma-separated values (take first)
                if ',' in af_str:
                    af_str = af_str.split(',')[0]
                return float(af_str)
        except (ValueError, AttributeError):
            pass
        
        return None
    
    def extract_reference_sequences(self, fasta_file):
        """
        Retrieve all reference sequences from the FASTA file.
        
        Args:
            fasta_file (str): Path to FASTA file
            
        Returns:
            dict: Dictionary mapping sequence IDs to sequences
            
        Raises:
            FileError: If FASTA file cannot be read
        """
        logger.info(f"Extracting sequences from {fasta_file}")
        
        if not os.path.exists(fasta_file):
            logger.error(f"FASTA file not found: {fasta_file}")
            raise FileError(f"FASTA file not found: {fasta_file}")
            
        sequences = {}
        
        try:
            with open(fasta_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequences[record.id] = str(record.seq)
                    
            logger.info(f"Extracted {len(sequences)} sequences from FASTA file")
            return sequences
            
        except Exception as e:
            logger.error(f"Failed to read FASTA file: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(f"Failed to read FASTA file: {e}")
    
    def mask_variants(self, sequence, variant_positions, flanking_size=0, use_soft_masking=False):
        """
        Mask variants in sequence with options for flanking regions and soft masking.
        
        Args:
            sequence (str): Input DNA sequence
            variant_positions (set): Set of variant positions (1-based)
            flanking_size (int): Number of bases to mask around each variant (default: 0)
            use_soft_masking (bool): Use lowercase letters instead of 'N' (default: False)
            
        Returns:
            str: Masked sequence
        """
        if not variant_positions:
            logger.debug("No variant positions to mask")
            return sequence
            
        sequence_list = list(sequence)
        masked_positions = set()
        
        # Calculate all positions to mask (including flanking regions)
        positions_to_mask = set()
        
        variant_pos_list = list(variant_positions)
        logger.debug(f"Processing {len(variant_pos_list)} variant positions with flanking_size={flanking_size}")
        
        if Config.SHOW_PROGRESS and len(variant_pos_list) > 1000:  # Only show progress for large sets
            variant_iter = tqdm(variant_pos_list, desc="Calculating mask positions")
        else:
            variant_iter = variant_pos_list
            
        for pos in variant_iter:
            # Convert 1-based VCF position to 0-based sequence index
            center_idx = pos - 1
            
            # Add flanking positions
            for offset in range(-flanking_size, flanking_size + 1):
                mask_idx = center_idx + offset
                if 0 <= mask_idx < len(sequence_list):
                    positions_to_mask.add(mask_idx)
        
        # Apply masking
        if Config.SHOW_PROGRESS and len(positions_to_mask) > 10000:
            mask_iter = tqdm(positions_to_mask, desc="Applying masks")
        else:
            mask_iter = positions_to_mask
            
        for idx in mask_iter:
            original_base = sequence_list[idx]
            
            if use_soft_masking:
                # Convert to lowercase for soft masking
                sequence_list[idx] = original_base.lower()
            else:
                # Use 'N' for hard masking
                sequence_list[idx] = 'N'
            
            if original_base.upper() in 'ATCG':  # Only count actual nucleotides as masked
                masked_positions.add(idx)
        
        masked_sequence = "".join(sequence_list)
        
        # Calculate masking statistics
        masked_count = len(masked_positions)
        masked_ratio = masked_count / len(sequence) if len(sequence) > 0 else 0
        
        mask_type = "soft" if use_soft_masking else "hard"
        logger.debug(f"Applied {mask_type} masking to {masked_count} positions ({masked_ratio:.2%} of sequence)")
        
        if flanking_size > 0:
            logger.debug(f"Flanking region size: {flanking_size} bases around each variant")
        
        return masked_sequence
    
    def mask_sequences_for_primer_design(self, sequences, variants, flanking_size=0, 
                                       use_soft_masking=False, min_af=None, min_qual=None):
        """
        Mask variants in sequences for primer design with advanced filtering options.
        
        Args:
            sequences (dict): Dictionary of sequences
            variants (dict): Dictionary mapping chromosomes to sets of variant positions
            flanking_size (int): Number of bases to mask around each variant (default: 0)
            use_soft_masking (bool): Use lowercase letters instead of 'N' (default: False)
            min_af (float, optional): Minimum allele frequency threshold for filtering
            min_qual (float, optional): Minimum QUAL score threshold for filtering
            
        Returns:
            dict: Dictionary of masked sequences
        """
        mask_type = "soft" if use_soft_masking else "hard"
        logger.info(f"Masking variants in sequences using {mask_type} masking...")
        
        if flanking_size > 0:
            logger.info(f"Using flanking region size: {flanking_size} bases")
        if min_af is not None:
            logger.info(f"AF filter: variants with AF >= {min_af}")
        if min_qual is not None:
            logger.info(f"QUAL filter: variants with QUAL >= {min_qual}")
        
        masked_sequences = {}
        
        if Config.SHOW_PROGRESS:
            sequence_iter = tqdm(sequences.items(), desc=f"Masking sequences ({mask_type})")
        else:
            sequence_iter = sequences.items()
        
        for seq_id, sequence in sequence_iter:
            # Try to find matching variants for this sequence
            # Assume sequence ID corresponds to chromosome name
            if seq_id in variants:
                variant_positions = variants[seq_id]
                masked_sequence = self.mask_variants(
                    sequence, 
                    variant_positions, 
                    flanking_size=flanking_size,
                    use_soft_masking=use_soft_masking
                )
                masked_sequences[seq_id] = masked_sequence
            else:
                # No variants found for this sequence, use original
                masked_sequences[seq_id] = sequence
                logger.debug(f"No variants found for sequence {seq_id}, using original")
        
        logger.info(f"Completed variant masking for {len(masked_sequences)} sequences")
        return masked_sequences
    
    def extract_variants_by_regions(self, vcf_file, conserved_regions, min_af=None, min_qual=None):
        """
        Extract variants from VCF file only for specific conserved regions with filtering.
        This optimizes memory usage by only loading variants in regions of interest.
        
        Args:
            vcf_file (str): Path to VCF file
            conserved_regions (dict): Dictionary mapping chromosome names to lists of conserved regions
            min_af (float, optional): Minimum allele frequency threshold
            min_qual (float, optional): Minimum QUAL score threshold
            
        Returns:
            dict: Dictionary mapping chromosomes to sets of variant positions
            
        Raises:
            FileError: If VCF file cannot be processed
        """
        logger.debug(f"Extracting variants for specific regions from {vcf_file}")
        logger.debug(f"Filters: min_af={min_af}, min_qual={min_qual}")
        
        if not os.path.exists(vcf_file):
            logger.error(f"VCF file not found: {vcf_file}")
            raise FileError(f"VCF file not found: {vcf_file}")
            
        if not conserved_regions:
            logger.warning("No conserved regions provided")
            return {}
        
        # Prepare the VCF file (compress and index if needed)
        try:
            prepared_vcf = self.prepare_vcf_file(vcf_file)
            
            # Dictionary to store variants by chromosome
            variants = {}
            
            # Process each chromosome and its conserved regions
            for chrom, regions in conserved_regions.items():
                if not regions:
                    continue
                    
                variants[chrom] = set()
                logger.debug(f"Processing {len(regions)} conserved regions for {chrom}")
                
                # Process each conserved region
                for region in regions:
                    start = region['start']
                    end = region['end']
                    
                    # Skip very small regions
                    if end - start < 10:
                        logger.debug(f"Skipping small region {chrom}:{start}-{end} (< 10 bp)")
                        continue
                        
                    # Extract variants for this specific region with filtering
                    region_variants = self.get_region_variants(
                        prepared_vcf, chrom, start, end, min_af=min_af, min_qual=min_qual
                    )
                    variants[chrom].update(region_variants)
            
            # Count total variants
            total_variants = sum(len(positions) for positions in variants.values())
            logger.info(f"Extracted {total_variants} variants across {len(variants)} chromosomes")
            
            return variants
            
        except Exception as e:
            logger.error(f"Failed to extract variants by regions: {e}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(f"Failed to extract variants by regions: {e}")
    
    def combine_masked_sequences(self, masked_sequences1, masked_sequences2):
        """
        Combine two sets of masked sequences, masking a position if it's masked in either set.
        
        Args:
            masked_sequences1 (dict): First set of masked sequences
            masked_sequences2 (dict): Second set of masked sequences
            
        Returns:
            dict: Combined masked sequences
        """
        combined_sequences = {}
        
        # Get all sequence IDs from both sets
        all_seq_ids = set(masked_sequences1.keys()) | set(masked_sequences2.keys())
        logger.info(f"Combining masking from two sets ({len(masked_sequences1)} and {len(masked_sequences2)} sequences)")
        
        for seq_id in all_seq_ids:
            # If sequence only exists in one set, use that version
            if seq_id not in masked_sequences1:
                logger.debug(f"Sequence {seq_id} only exists in second set")
                combined_sequences[seq_id] = masked_sequences2[seq_id]
                continue
            elif seq_id not in masked_sequences2:
                logger.debug(f"Sequence {seq_id} only exists in first set")
                combined_sequences[seq_id] = masked_sequences1[seq_id]
                continue
            
            # If sequence exists in both sets, combine the masking
            seq1 = masked_sequences1[seq_id]
            seq2 = masked_sequences2[seq_id]
            
            # Ensure sequences are the same length
            if len(seq1) != len(seq2):
                logger.warning(f"Sequence length mismatch for {seq_id}. Using first sequence.")
                combined_sequences[seq_id] = seq1
                continue
            
            # Combine masking - preserve masking from either sequence
            combined_seq = ""
            for i in range(len(seq1)):
                char1, char2 = seq1[i], seq2[i]
                
                # If either position is hard masked ('N'), use hard masking
                if char1 == 'N' or char2 == 'N':
                    combined_seq += 'N'
                # If either position is soft masked (lowercase), use soft masking
                elif char1.islower() or char2.islower():
                    # Use the lowercase version, preferring the original base if available
                    if char1.islower():
                        combined_seq += char1
                    else:
                        combined_seq += char2.lower()
                # If neither is masked, use the original base (should be the same)
                else:
                    combined_seq += char1.upper()
            
            # Calculate masking statistics
            n_count = combined_seq.count('N')
            lowercase_count = sum(1 for c in combined_seq if c.islower())
            total_masked = n_count + lowercase_count
            mask_ratio = total_masked / len(combined_seq) if len(combined_seq) > 0 else 0
            
            logger.debug(f"Combined sequence {seq_id}: {total_masked}/{len(combined_seq)} bases masked ({mask_ratio:.2%})")
            logger.debug(f"  Hard masked (N): {n_count}, Soft masked (lowercase): {lowercase_count}")
            
            combined_sequences[seq_id] = combined_seq
        
        return combined_sequences