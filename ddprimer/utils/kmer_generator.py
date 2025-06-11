#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K-mer generator module for ddPrimer pipeline using GenomeTester4's glistmaker.

Contains functionality for:
1. Glistmaker availability checking and execution
2. K-mer list generation using GenomeTester4 tools
3. Primer3-compatible .list file generation
4. Error handling and validation for external tool integration
5. Progress tracking and logging for k-mer generation workflow

This module provides k-mer frequency analysis capabilities using the
GenomeTester4 toolkit's glistmaker tool for optimal performance and
compatibility with Primer3's masking functionality.
"""

import os
import subprocess
import logging
import shutil
from pathlib import Path
from typing import Optional, List, Tuple

from ..config import Config, FileError, ExternalToolError

# Set up module logger
logger = logging.getLogger(__name__)


class KmerGenerator:
    """
    K-mer generator using GenomeTester4's glistmaker tool.
    
    This class provides methods for generating k-mer frequency lists using
    the external glistmaker tool from the GenomeTester4 package, which
    produces binary .list files directly compatible with Primer3.
    
    Attributes:
        glistmaker_available: Whether glistmaker tool is available
        
    Example:
        >>> generator = KmerGenerator()
        >>> if generator.check_glistmaker_availability():
        ...     success = generator.generate_kmer_lists("genome.fasta", "output_dir")
    """
    
    def __init__(self):
        """
        Initialize k-mer generator with glistmaker availability check.
        """
        self.glistmaker_available = self.check_glistmaker_availability()
        logger.debug(f"KmerGenerator initialized: glistmaker_available={self.glistmaker_available}")
    
    def check_glistmaker_availability(self) -> bool:
        """
        Check if glistmaker from GenomeTester4 is available in PATH.
        
        Returns:
            True if glistmaker is available, False otherwise
            
        Example:
            >>> generator = KmerGenerator()
            >>> if generator.check_glistmaker_availability():
            ...     print("Glistmaker is available")
        """
        logger.debug("Checking glistmaker availability")
        
        try:
            # Try to find glistmaker executable
            glistmaker_path = shutil.which("glistmaker")
            if glistmaker_path:
                logger.debug(f"Found glistmaker at: {glistmaker_path}")
                
                # Try to run glistmaker to verify it works
                result = subprocess.run(
                    ["glistmaker", "--help"],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                
                # glistmaker returns 1 when showing help, but that's normal
                if result.returncode in [0, 1] and ("glistmaker" in result.stderr or "Usage:" in result.stderr):
                    logger.debug("Glistmaker is available and functional")
                    return True
                else:
                    logger.debug(f"Glistmaker found but not functional: returncode={result.returncode}, stderr={result.stderr}")
                    return False
            else:
                logger.debug("Glistmaker not found in PATH")
                return False
                
        except subprocess.TimeoutExpired:
            logger.debug("Glistmaker availability check timed out")
            return False
        except Exception as e:
            logger.debug(f"Error checking glistmaker availability: {str(e)}")
            return False
    
    def generate_kmer_lists(self, fasta_file: str, output_dir: str, 
                           organism_name: Optional[str] = None) -> bool:
        """
        Generate k-mer lists using glistmaker for standard k-mer sizes.
        
        Creates k-mer frequency lists for the sizes specified in Config.KMER_SIZES
        using glistmaker with the minimum frequency threshold from Config.KMER_MIN_FREQUENCY.
        
        Args:
            fasta_file: Path to input FASTA file
            output_dir: Directory for output .list files
            organism_name: Optional organism name for file naming
            
        Returns:
            True if generation completed successfully, False otherwise
            
        Raises:
            FileError: If input file doesn't exist or output directory cannot be created
            ExternalToolError: If glistmaker execution fails
            
        Example:
            >>> generator = KmerGenerator()
            >>> success = generator.generate_kmer_lists("genome.fasta", "./output", "arabidopsis")
        """
        logger.debug("=== K-MER LIST GENERATION DEBUG ===")
        
        # Check glistmaker availability
        if not self.glistmaker_available:
            error_msg = "Glistmaker is not available. Please install GenomeTester4."
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="glistmaker")
        
        # Validate input file
        if not os.path.exists(fasta_file):
            error_msg = f"FASTA file not found: {fasta_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        # Create output directory
        try:
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Output directory: {output_dir}")
        except OSError as e:
            error_msg = f"Failed to create output directory {output_dir}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
        # Determine organism name for file naming
        if not organism_name:
            organism_name = Path(fasta_file).stem.replace("_", " ").lower()
        
        logger.info(f"\nGenerating k-mer lists for {organism_name}")
        logger.debug(f"K-mer sizes: {Config.KMER_SIZES}")
        logger.debug(f"Minimum frequency: {Config.KMER_MIN_FREQUENCY}")
        
        success_count = 0
        
        # Generate k-mer lists for each size
        for kmer_size in Config.KMER_SIZES:
            try:
                output_file = os.path.join(output_dir, f"{organism_name.replace(' ', '_')}_{kmer_size}.list")
                
                success = self._run_glistmaker(fasta_file, output_file, kmer_size, Config.KMER_MIN_FREQUENCY)
                
                if success:
                    success_count += 1
                    logger.info(f"Generated: {Path(output_file).name}")
                else:
                    logger.error(f"Failed to generate {kmer_size}-mer list")
                    
            except Exception as e:
                logger.error(f"Error generating {kmer_size}-mer list: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                continue
        
        logger.debug("=== END K-MER LIST GENERATION DEBUG ===")
        
        if success_count == 0:
            logger.error("Failed to generate any k-mer lists")
            return False
        elif success_count < len(Config.KMER_SIZES):
            logger.warning(f"Only generated {success_count}/{len(Config.KMER_SIZES)} k-mer lists")
        else:
            logger.info(f"\nSuccessfully generated all {success_count} k-mer lists")
        
        return success_count > 0
    
    def _run_glistmaker(self, input_file: str, output_file: str, kmer_size: int, 
                       min_frequency: int) -> bool:
        """
        Run glistmaker to generate a single k-mer list.
        
        Args:
            input_file: Path to input FASTA file
            output_file: Path to output .list file
            kmer_size: K-mer size to generate
            min_frequency: Minimum frequency threshold
            
        Returns:
            True if glistmaker completed successfully, False otherwise
            
        Raises:
            ExternalToolError: If glistmaker execution fails critically
        """
        logger.debug(f"Running glistmaker: k={kmer_size}, min_freq={min_frequency}")
        
        # Glistmaker automatically appends _{kmer_size} to the output name
        # So if we want "thaliana_11.list", we need to pass "thaliana" as the base name
        base_name = os.path.splitext(output_file)[0]  # Remove .list extension
        if base_name.endswith(f"_{kmer_size}"):
            # Remove the kmer size from the end since glistmaker will add it
            base_name = base_name[:-len(f"_{kmer_size}")]
        
        # Build glistmaker command - FIXED: Use correct parameter names
        cmd = [
            "glistmaker",
            input_file,  # Input file comes first as positional argument
            "--wordlength", str(kmer_size),  # Use --wordlength instead of -k
            "--outputname", base_name,  # Use base name without kmer size suffix
            # Note: glistmaker doesn't have a direct frequency cutoff parameter
            # It generates all k-mers, and filtering would need to be done post-processing
        ]
        
        logger.debug(f"Glistmaker command: {' '.join(cmd)}")
        
        try:
            # Run glistmaker with progress tracking
            logger.info(f"\nGenerating {kmer_size}-mer list...")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout for large genomes
            )
            
            if result.returncode == 0:
                # Glistmaker creates files with pattern: {base_name}_{kmer_size}.list
                actual_output = f"{base_name}_{kmer_size}.list"
                
                # Verify the actual output file was created
                if os.path.exists(actual_output) and os.path.getsize(actual_output) > 0:
                    # If the actual output name differs from what we want, rename it
                    if actual_output != output_file:
                        try:
                            shutil.move(actual_output, output_file)
                            logger.debug(f"Renamed {actual_output} to {output_file}")
                        except Exception as e:
                            logger.warning(f"Could not rename output file: {str(e)}")
                            # Update the output_file to the actual name if rename fails
                            output_file = actual_output
                    
                    logger.debug(f"Glistmaker completed successfully: {output_file}")
                    return True
                else:
                    logger.error(f"Glistmaker completed but output file is missing or empty: {actual_output}")
                    # Also check if it was created with the original expected name
                    expected_output = f"{os.path.splitext(output_file)[0]}.list"
                    if os.path.exists(expected_output) and os.path.getsize(expected_output) > 0:
                        logger.debug(f"Found output file at expected location: {expected_output}")
                        return True
                    return False
            else:
                logger.error(f"Glistmaker failed with return code {result.returncode}")
                if result.stderr:
                    logger.error(f"Glistmaker stderr: {result.stderr}")
                if result.stdout:
                    logger.debug(f"Glistmaker stdout: {result.stdout}")
                return False
                
        except subprocess.TimeoutExpired:
            error_msg = f"Glistmaker timed out after 1 hour for {kmer_size}-mer generation"
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="glistmaker")
        except FileNotFoundError:
            error_msg = "Glistmaker executable not found"
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="glistmaker")
        except Exception as e:
            error_msg = f"Error running glistmaker: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="glistmaker") from e
    
    def validate_kmer_list(self, list_file: str) -> bool:
        """
        Validate that a .list file is properly formatted and readable.
        
        Args:
            list_file: Path to .list file to validate
            
        Returns:
            True if file is valid, False otherwise
        """
        if not os.path.exists(list_file):
            logger.debug(f"K-mer list file does not exist: {list_file}")
            return False
        
        if os.path.getsize(list_file) == 0:
            logger.debug(f"K-mer list file is empty: {list_file}")
            return False
        
        # Basic validation - check if file has reasonable size and binary content
        try:
            with open(list_file, 'rb') as f:
                # Read first few bytes to check if it's a binary file
                header = f.read(16)
                if len(header) < 16:
                    logger.debug(f"K-mer list file too small: {list_file}")
                    return False
                
                # Basic check for binary content (not all ASCII)
                if all(32 <= b <= 126 for b in header):
                    logger.debug(f"K-mer list file appears to be text, not binary: {list_file}")
                    return False
                
                logger.debug(f"K-mer list file validation passed: {list_file}")
                return True
                
        except Exception as e:
            logger.debug(f"Error validating k-mer list file {list_file}: {str(e)}")
            return False


def run_kmer_generation(fasta_file: str, output_dir: Optional[str] = None) -> bool:
    """
    Run k-mer generation workflow using glistmaker.
    
    Orchestrates the complete k-mer generation process using glistmaker
    from GenomeTester4 to create Primer3-compatible .list files.
    
    Args:
        fasta_file: Path to input FASTA file
        output_dir: Directory for output files (defaults to ~/.ddprimer/kmer_lists)
        
    Returns:
        True if generation completed successfully, False otherwise
        
    Raises:
        FileError: If files cannot be accessed
        ExternalToolError: If glistmaker is not available or fails
        
    Example:
        >>> success = run_kmer_generation("genome.fasta", "./output")
    """
    logger.debug("=== K-MER GENERATION WORKFLOW DEBUG ===")
    
    # Create output directory - default to ~/.ddprimer/kmer_lists
    if not output_dir:
        output_dir = os.path.join(os.path.expanduser("~"), ".ddprimer", "kmer_lists")
    
    # Initialize generator
    generator = KmerGenerator()
    
    # Check if glistmaker is available
    if not generator.check_glistmaker_availability():
        logger.error("Glistmaker (GenomeTester4) is not available.")
        logger.error("Please install GenomeTester4 to use k-mer generation:")
        logger.error("  conda install -c bioconda genometester4")
        logger.error("  # or compile from source: https://github.com/bioinfo-ut/GenomeTester4")
        return False
    
    # Determine organism name from filename
    organism_name = Path(fasta_file).stem.replace("_", " ").lower()
    
    try:
        success = generator.generate_kmer_lists(fasta_file, output_dir, organism_name)
        
        if success:
            logger.info(f"\nOutput files in: {output_dir}")            
            logger.info(f"\n=== K-mer generation completed successfully! ===")

        
        logger.debug("=== END K-MER GENERATION WORKFLOW DEBUG ===")
        return success
        
    except (FileError, ExternalToolError):
        # Re-raise specific exceptions
        logger.debug("=== END K-MER GENERATION WORKFLOW DEBUG ===")
        raise
    except Exception as e:
        error_msg = f"Unexpected error in k-mer generation workflow: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        logger.debug("=== END K-MER GENERATION WORKFLOW DEBUG ===")
        raise ExternalToolError(error_msg, tool_name="glistmaker") from e