def AlignmentWorkflow(args, output_dir, logger):
    """
    Run the alignment primer design workflow with improved order of operations:
    1. Run alignment first (or use pre-computed MAF)
    2. Identify conserved regions
    3. Apply SNP/indel masking from both genomes (optional)
    4. Design primers on fully masked sequences
    
    Args:
        args: Command line arguments
        output_dir: Output directory path
        logger: Logger object
    
    Returns:
        tuple: (final_masked_sequences, coordinate_map) - Masked sequences and coordinate mapping
    """
    import os
    import traceback
    from .lastz_runner import LastZRunner
    from .maf_parser import MAFParser  
    from ..utils import FileUtils
    from ..core import SNPMaskingProcessor
    from ..config import Config
    
    # Initialize processors
    maf_parser = MAFParser()
    snp_processor = SNPMaskingProcessor()
    
    # Step 1: Handle MAF file (pre-computed or generate new one)
    if args.maf_file:
        logger.info("\n>>> Using pre-computed MAF file <<<")
        maf_file = args.maf_file
        logger.debug(f"Using MAF file: {maf_file}")
        
        # With pre-computed MAF, we don't need the FASTA files
        ref_fasta_required = False
        second_fasta_required = False
    else:
        # Need FASTA files to run alignment
        ref_fasta_required = True
        second_fasta_required = True
        
        # Check if required files are provided
        if not args.fasta:
            raise ValueError("Reference genome FASTA (--fasta) is required for alignment")
        if not args.second_fasta:
            raise ValueError("Second FASTA (--second-fasta) is required for alignment")
        
        logger.info("\n>>> Running LastZ alignment between genomes <<<")

        # Create alignment directory
        input_dir = os.path.dirname(os.path.abspath(args.fasta))
        alignments_dir = os.path.join(input_dir, "Alignments")
        os.makedirs(alignments_dir, exist_ok=True)
                    
        # Ensure lastz_options include MAF format
        lastz_options = args.lastz_options
        if "--format=maf" not in lastz_options and "-f maf" not in lastz_options:
            lastz_options += " --format=maf"
                    
        # Run LastZ alignment
        try:
            lastz_runner = LastZRunner()
            maf_file = lastz_runner.run_parallel_alignment(
                args.fasta,
                args.second_fasta,
                alignments_dir,  # Use the alignments_dir directly, not alignment_dir
                lastz_options
            )
            logger.info(f"LastZ alignment completed: {maf_file}")
        except Exception as e:
            logger.error(f"Error running LastZ alignment: {e}")
            logger.debug(traceback.format_exc())
            raise
    
    # Step 2: Parse MAF file and identify conserved regions
    logger.info("\n>>> Parsing alignment and identifying conserved regions <<<")
    try:
        # First analyze the MAF file to understand its structure
        maf_analysis = maf_parser.analyze_maf_file(maf_file)
        logger.debug(f"MAF file contains {maf_analysis['alignment_count']} alignment blocks")
        logger.debug(f"Reference sequences: {', '.join(maf_analysis['ref_seq_ids'])}")
        logger.debug(f"Query sequences: {', '.join(maf_analysis['query_seq_ids'])}")
        
        # Parse MAF file
        alignments = maf_parser.parse_maf_file(maf_file)
        
        # Identify conserved regions
        conserved_regions = maf_parser.identify_conserved_regions(
            args.min_identity,
            args.min_length
        )
        
        # Generate coordinate mapping between reference and query genomes
        coordinate_map = maf_parser.generate_coordinate_map(conserved_regions)
        logger.debug("Generated coordinate mapping between reference and second genome")
        
        total_regions = sum(len(regions) for regions in conserved_regions.values())
        logger.info(f"Identified {total_regions} conserved regions across {len(conserved_regions)} chromosomes")
    except Exception as e:
        logger.error(f"Error processing MAF file: {e}")
        logger.debug(traceback.format_exc())
        raise
    
    # Step 3: Prepare for masking by getting variants from both genomes (if SNP masking is enabled)
    ref_variants = {}
    second_variants = {}
    
    # Get reference genome variants if VCF provided and SNP masking is enabled
    if not args.no_snp_masking and args.vcf:
        logger.info("\n>>> Extracting variants from reference genome VCF <<<")
        try:
            ref_variants = snp_processor.get_variant_positions(args.vcf)
            total_ref_variants = sum(len(positions) for positions in ref_variants.values())
            logger.info(f"Extracted {total_ref_variants} variants from reference genome")
        except Exception as e:
            logger.error(f"Error extracting reference variants: {e}")
            logger.debug(traceback.format_exc())
            raise

    # Get second genome variants if VCF provided and SNP masking is enabled
    if not args.no_snp_masking and args.second_vcf:
        logger.info("\n>>> Extracting variants from second VCF <<<")
        try:
            second_variants = snp_processor.get_variant_positions(args.second_vcf)
            total_second_variants = sum(len(positions) for positions in second_variants.values())
            logger.info(f"Extracted {total_second_variants} variants from second genome")
        except Exception as e:
            logger.error(f"Error extracting second genome variants: {e}")
            logger.debug(traceback.format_exc())
            raise
    
    # Step 4: Load reference FASTA if needed for masking
    reference_sequences = {}
    
    # If we have a pre-computed MAF and don't need to run alignment
    if not ref_fasta_required:
        # We need to generate a minimal reference FASTA from the MAF file
        # containing only the conserved regions
        logger.info("\n>>> Generating reference sequences from MAF file <<<")
        try:
            reference_sequences = maf_parser.extract_reference_sequences_from_maf()
            logger.debug(f"Generated {len(reference_sequences)} reference sequences from MAF file")
            
            # Write to a temporary FASTA file for later use
            temp_ref_fasta = os.path.join(output_dir, "ref_from_maf.fasta")
            with open(temp_ref_fasta, 'w') as f:
                for seq_id, seq in reference_sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")
            
            logger.debug(f"Wrote reference sequences to temporary file: {temp_ref_fasta}")
            
        except Exception as e:
            logger.error(f"Error generating reference sequences from MAF: {e}")
            logger.debug(traceback.format_exc())
            raise
    else:
        # Load reference sequences from provided FASTA
        logger.info("\n>>> Loading reference sequences from FASTA file <<<")
        try:
            reference_sequences = FileUtils.load_fasta(args.fasta)
            logger.debug(f"Loaded {len(reference_sequences)} sequences from reference FASTA")
        except Exception as e:
            logger.error(f"Error loading reference FASTA: {e}")
            logger.debug(traceback.format_exc())
            raise
    
    # Step 5: Apply masking to the reference sequences
    logger.info("\n>>> Creating masked genome with conserved regions <<<")
    
    # First create a masked reference with only conserved regions
    masked_fasta_path = os.path.join(output_dir, "masked_reference.fasta")
    try:
        maf_parser.mask_non_conserved_regions(
            reference_sequences if not ref_fasta_required else args.fasta,
            masked_fasta_path,
            conserved_regions,
            args.min_identity
        )
        logger.debug(f"Created alignment-masked reference genome: {masked_fasta_path}")
    except Exception as e:
        logger.error(f"Error masking non-conserved regions: {e}")
        logger.debug(traceback.format_exc())
        raise
    
    # Load the alignment-masked sequences
    alignment_masked_sequences = FileUtils.load_fasta(masked_fasta_path)
    logger.debug(f"Loaded {len(alignment_masked_sequences)} alignment-masked sequences")
    
    # Now mask variants from both genomes if SNP masking is enabled
    if args.no_snp_masking:
        logger.info("\n>>> Skipping SNP masking as requested <<<")
        # Use alignment-masked sequences directly without SNP masking
        final_masked_sequences = alignment_masked_sequences
    else:
        logger.info("\n>>> Applying SNP masking from VCF files <<<")
        final_masked_sequences = {}
        
        # Process each sequence
        for seq_id, sequence in alignment_masked_sequences.items():
            # Get variants for reference genome
            ref_seq_variants = ref_variants.get(seq_id, set())
            
            # Get mapped variants from second genome
            second_seq_variants = set()
            
            # Map second genome variants to reference coordinates
            if seq_id in coordinate_map and second_variants:
                for second_chrom, positions in second_variants.items():
                    # Find all mappings from second genome to reference
                    for pos in positions:
                        # Check all reference positions to find mappings (inefficient but works)
                        for ref_pos, mapping in coordinate_map[seq_id].items():
                            if mapping["qry_src"] == second_chrom and mapping["qry_pos"] == pos:
                                # Found a mapping from second genome variant to reference
                                second_seq_variants.add(ref_pos)
            
            # Combine variants from both genomes
            all_variants = ref_seq_variants.union(second_seq_variants)
            logger.debug(f"Sequence {seq_id}: {len(ref_seq_variants)} reference variants, "
                        f"{len(second_seq_variants)} mapped second genome variants, "
                        f"{len(all_variants)} total variants")
            
            # Mask all variants in the sequence
            if all_variants:
                logger.debug(f"Masking {len(all_variants)} variants in {seq_id}...")
                try:
                    variant_masked_seq = snp_processor.mask_variants(sequence, all_variants)
                    final_masked_sequences[seq_id] = variant_masked_seq
                except Exception as e:
                    logger.error(f"Error masking variants in {seq_id}: {e}")
                    logger.debug(traceback.format_exc())
                    raise
            else:
                logger.debug(f"No variants to mask in {seq_id}")
                final_masked_sequences[seq_id] = sequence
    
    # Write final masked sequences to file for debugging/verification
    final_masked_path = os.path.join(output_dir, "final_masked.fasta")
    with open(final_masked_path, 'w') as f:
        for seq_id, seq in final_masked_sequences.items():
            f.write(f">{seq_id}\n{seq}\n")
    
    if args.no_snp_masking:
        logger.info(f"Created alignment-masked reference genome (no SNP masking): {final_masked_path}")
    else:
        logger.info(f"Created fully masked reference genome (with SNP masking): {final_masked_path}")
    logger.info(f"Ready for primer design on {len(final_masked_sequences)} masked sequences")

    # Clean up intermediate files if not in debug mode
    if not Config.DEBUG_MODE:
        try:
            # List of files to clean up
            files_to_clean = [
                os.path.join(output_dir, "masked_reference.fasta"),
                final_masked_path  # This is final_masked.fasta
            ]
            
            for file_path in files_to_clean:
                if os.path.exists(file_path):
                    logger.debug(f"Cleaning up intermediate file: {file_path}")
                    os.remove(file_path)
        except Exception as e:
            logger.warning(f"Error during cleanup of masked files: {e}")
            logger.debug(traceback.format_exc())
    
    return final_masked_sequences, coordinate_map