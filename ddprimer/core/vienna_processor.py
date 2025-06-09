#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thermodynamic calculation module using ViennaRNA CLI for DNA oligos.

Provides thermodynamic calculations for DNA sequences using ViennaRNA
command-line interface with DNA-specific parameters and salt corrections.
Contains functionality for:
1. DNA parameter file detection and validation
2. Minimum free energy calculations via RNAfold CLI
3. Batch processing with progress tracking
4. ViennaRNA installation validation and configuration reporting

This module integrates with the broader ddPrimer pipeline to provide
robust thermodynamic analysis capabilities for primer design workflows.
"""

import re
import logging
import subprocess
import os
from pathlib import Path
from typing import List, Optional, Union
import pandas as pd
from tqdm import tqdm

from ..config import Config, SequenceProcessingError, ExternalToolError

# Type alias for path inputs
PathLike = Union[str, Path]

# Set up module logger
logger = logging.getLogger(__name__)


class ViennaRNAProcessor:
    """
    Handles thermodynamic calculations using ViennaRNA CLI for DNA oligos.
    
    This class provides methods for calculating minimum free energy (ΔG)
    of DNA sequences using ViennaRNA command-line interface with DNA-specific parameters.
    
    Attributes:
        _param_file_cache: Class-level cache for parameter file path
        
    Example:
        >>> processor = ViennaRNAProcessor()
        >>> deltaG = processor.calc_deltaG("ATCGATCG")
        >>> batch_results = processor.calc_deltaG_batch(sequences)
    """

    _param_file_cache = None  # Cache for parameter file path

    @classmethod
    def _find_dna_parameter_file(cls) -> Optional[Path]:
        """
        Find the DNA parameter file in standard ViennaRNA locations.
        
        Searches for DNA parameter files in common ViennaRNA installation
        locations with preference for newer parameter sets.
        
        Returns:
            Path to DNA parameter file if found, None otherwise
            
        Example:
            >>> param_file = ViennaRNAProcessor._find_dna_parameter_file()
            >>> if param_file:
            ...     print(f"Found DNA parameters: {param_file}")
        """
        # Return cached result if available
        if cls._param_file_cache is not None:
            return cls._param_file_cache
            
        try:
            # First try to find RNAfold binary path
            try:
                rnafold_path = subprocess.run(
                    ["which", "RNAfold"],
                    capture_output=True,
                    text=True,
                    check=True
                ).stdout.strip()

                if rnafold_path:
                    # Assume param file is in ../share/ViennaRNA/dna_mathews2004.par relative to RNAfold
                    base_dir = os.path.abspath(os.path.join(os.path.dirname(rnafold_path), "..", "share", "ViennaRNA"))
                    param_file = Path(base_dir) / "dna_mathews2004.par"

                    if param_file.exists() and param_file.is_file():
                        logger.debug(f"Found DNA parameter file via RNAfold path: {param_file}")
                        cls._param_file_cache = param_file
                        return param_file
            except subprocess.CalledProcessError:
                logger.debug("Could not find RNAfold using 'which' command")
            
            # Try different DNA parameter files in order of preference
            candidate_files = [
                # System-wide installation paths
                Path("/usr/local/share/ViennaRNA/dna_mathews2004.par"),
                Path("/usr/share/ViennaRNA/dna_mathews2004.par"),
                Path("/opt/ViennaRNA/share/ViennaRNA/dna_mathews2004.par"),
                # Alternative parameter files
                Path("/usr/local/share/ViennaRNA/dna_mathews1999.par"),
                Path("/usr/share/ViennaRNA/dna_mathews1999.par"),
                Path("/opt/ViennaRNA/share/ViennaRNA/dna_mathews1999.par"),
                Path("/usr/local/share/ViennaRNA/dna_mathews.par"),
                Path("/usr/share/ViennaRNA/dna_mathews.par"),
                Path("/opt/ViennaRNA/share/ViennaRNA/dna_mathews.par")
            ]
            
            for param_file in candidate_files:
                if param_file.exists() and param_file.is_file():
                    logger.debug(f"Found DNA parameter file: {param_file}")
                    cls._param_file_cache = param_file
                    return param_file
                    
            logger.debug("No DNA parameter file found, will use default RNA parameters")
            cls._param_file_cache = None
            return None
            
        except Exception as e:
            logger.debug(f"Error finding DNA parameter file: {e}")
            cls._param_file_cache = None
            return None

    @classmethod
    def _run_rnafold_cli(cls, seq: str, use_dna_params: bool = True) -> Optional[float]:
        """
        Run RNAfold via command line interface.
        
        Args:
            seq: DNA sequence (A,C,G,T,N)
            use_dna_params: Whether to use DNA parameter file if available
            
        Returns:
            Minimum free energy in kcal/mol, or None if calculation fails
            
        Raises:
            ExternalToolError: If RNAfold execution fails
        """
        if not seq or not isinstance(seq, str):
            raise SequenceProcessingError("Invalid sequence provided for RNAfold")
            
        try:
            # Convert T to U for ViennaRNA (even when using DNA parameters)
            rna_seq = seq.upper().replace("T", "U")
            
            # Validate converted sequence
            if not cls._is_valid_rna_sequence(rna_seq):
                raise SequenceProcessingError(f"Invalid sequence after T->U conversion: {seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...")

            # Build RNAfold command
            cmd = ["RNAfold"]
            
            # Disable PostScript output to avoid creating rna.ps files
            cmd.append("--noPS")
            
            # Add DNA parameter file if available and requested
            if use_dna_params:
                param_file = cls._find_dna_parameter_file()
                if param_file:
                    cmd.append(f"--paramFile={param_file}")
                    logger.debug(f"Using DNA parameter file: {param_file}")
            
            # Add temperature setting
            cmd.append(f"--temp={Config.THERMO_TEMPERATURE}")
            
            # Run RNAfold
            try:
                result = subprocess.run(
                    cmd,
                    input=rna_seq.encode(),
                    capture_output=True,
                    check=True,
                    timeout=30  # Add timeout to prevent hanging
                )
                
                # Parse output to extract energy value
                output = result.stdout.decode()
                energy = cls._parse_rnafold_output(output)
                
                logger.debug(f"RNAfold CLI successful: {energy:.3f} kcal/mol for {len(seq)} bp sequence")
                return energy
                
            except subprocess.CalledProcessError as e:
                error_msg = f"RNAfold CLI execution failed: {e.stderr.decode() if e.stderr else str(e)}"
                logger.warning(error_msg)
                raise ExternalToolError(error_msg, tool_name="RNAfold") from e
            except subprocess.TimeoutExpired:
                error_msg = "RNAfold CLI execution timed out"
                logger.warning(error_msg)
                raise ExternalToolError(error_msg, tool_name="RNAfold")
                
        except SequenceProcessingError:
            # Re-raise without wrapping
            raise
        except ExternalToolError:
            # Re-raise without wrapping
            raise
        except Exception as e:
            error_msg = f"Error running RNAfold CLI for sequence '{seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...': {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="RNAfold") from e

    @staticmethod
    def _parse_rnafold_output(output: str) -> float:
        """
        Parse RNAfold output to extract energy value.
        
        Args:
            output: Raw output from RNAfold command
            
        Returns:
            Energy value in kcal/mol
            
        Raises:
            ExternalToolError: If output cannot be parsed
        """
        try:
            lines = output.strip().split('\n')
            
            # Look for energy value in parentheses, e.g., ( -1.20)
            for line in lines:
                if '(' in line and ')' in line:
                    # Extract energy value from parentheses
                    energy_match = re.search(r'\(\s*(-?\d+\.?\d*)\s*\)', line)
                    if energy_match:
                        return float(energy_match.group(1))
            
            error_msg = "Could not find energy value in RNAfold output"
            raise ExternalToolError(error_msg, tool_name="RNAfold")
            
        except ValueError as e:
            error_msg = f"Could not parse energy value from RNAfold output: {e}"
            raise ExternalToolError(error_msg, tool_name="RNAfold") from e
        except Exception as e:
            error_msg = f"Error parsing RNAfold output: {e}"
            raise ExternalToolError(error_msg, tool_name="RNAfold") from e

    @staticmethod
    def _is_valid_rna_sequence(seq: str) -> bool:
        """
        Validate that a sequence contains only valid RNA characters.
        
        Args:
            seq: RNA sequence to validate
            
        Returns:
            True if sequence is valid RNA, False otherwise
        """
        if not seq:
            return False
            
        rna_pattern = re.compile(r"^[ACGUNacgun]+$")
        return bool(rna_pattern.match(seq))

    @staticmethod
    def _is_valid_dna_sequence(seq: str) -> bool:
        """
        Validate that a sequence contains only valid DNA characters.
        
        Args:
            seq: Sequence to validate
            
        Returns:
            True if sequence is valid DNA
            
        Example:
            >>> ViennaRNAProcessor._is_valid_dna_sequence("ATCG")
            True
            >>> ViennaRNAProcessor._is_valid_dna_sequence("ATCGX")
            False
        """
        if not seq:
            return False
            
        dna_pattern = re.compile(r"^[ACGTNacgtn]+$")
        return bool(dna_pattern.match(seq))

    @classmethod
    def calc_deltaG(cls, seq: str) -> Optional[float]:
        """
        Calculate the minimum free energy (ΔG, kcal/mol) of a DNA oligo using CLI.
        
        Computes the minimum free energy for DNA secondary structure formation
        using ViennaRNA command-line interface with DNA-specific parameters.
        
        Args:
            seq: DNA sequence string (A, T, C, G, N allowed)
            
        Returns:
            Minimum free energy in kcal/mol, or None if calculation fails
            
        Raises:
            SequenceProcessingError: If sequence is invalid
            
        Example:
            >>> processor = ViennaRNAProcessor()
            >>> deltaG = processor.calc_deltaG("ATCGATCGATCG")
            >>> if deltaG is not None:
            ...     print(f"ΔG = {deltaG:.2f} kcal/mol")
        """
        # Validate input
        if not isinstance(seq, str) or seq == "":
            logger.debug("Empty or invalid sequence provided")
            return None

        # Validate DNA sequence format
        if not cls._is_valid_dna_sequence(seq):
            error_msg = f"Invalid DNA sequence: contains non-DNA characters"
            logger.debug(f"Invalid DNA sequence format: {seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...")
            raise SequenceProcessingError(error_msg)

        try:
            # Use CLI interface with DNA parameters
            energy = cls._run_rnafold_cli(seq, use_dna_params=True)
            return energy
            
        except SequenceProcessingError:
            # Already logged in _run_rnafold_cli
            raise
        except ExternalToolError:
            # Already logged in _run_rnafold_cli
            raise
        except Exception as e:
            error_msg = f"ViennaRNA calculation failed for sequence {seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...: {str(e)}"
            logger.warning(error_msg)
            logger.debug(f"ViennaRNA error details: {str(e)}", exc_info=True)
            return None

    @classmethod
    def calc_deltaG_batch(
        cls, 
        seqs: List[str], 
        description: str = "Calculating ΔG with ViennaRNA CLI"
    ) -> List[Optional[float]]:
        """
        Calculate ΔG for a batch of sequences with progress tracking.
        
        Processes multiple sequences efficiently with optional progress display
        and error handling for individual sequence failures.
        
        Args:
            seqs: List of DNA sequences
            description: Description for the progress bar
            
        Returns:
            List of ΔG values (same length as input, with None for failed calculations)
            
        Raises:
            SequenceProcessingError: If batch processing setup fails
            
        Example:
            >>> sequences = ["ATCG", "GCTA", "AAAA"]
            >>> processor = ViennaRNAProcessor()
            >>> results = processor.calc_deltaG_batch(sequences)
            >>> len(results) == len(sequences)
            True
        """
        if not isinstance(seqs, list):
            error_msg = "Sequences must be provided as a list"
            raise SequenceProcessingError(error_msg)
            
        logger.info(f"Processing batch of {len(seqs)} sequences for ΔG calculation using CLI")
        logger.debug(f"Progress tracking: {Config.SHOW_PROGRESS}")
        
        results = []
        failed_count = 0
        
        try:
            # Apply progress tracking if enabled
            sequence_iter = tqdm(seqs, desc=description) if Config.SHOW_PROGRESS else seqs
                
            for seq_num, seq in enumerate(sequence_iter):
                if pd.notnull(seq) and seq:
                    try:
                        result = cls.calc_deltaG(seq)
                        results.append(result)
                        if result is None:
                            failed_count += 1
                    except Exception as e:
                        logger.debug(f"ΔG calculation failed for sequence {seq_num}: {e}")
                        results.append(None)
                        failed_count += 1
                else:
                    results.append(None)
                    
            logger.debug(f"Completed ΔG calculations: {len(results)} total, {failed_count} failed")
            return results
            
        except Exception as e:
            error_msg = f"Error during batch ΔG calculation: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
        
    @classmethod
    def process_deltaG_series(
        cls, 
        series: pd.Series, 
        description: str = "Processing sequences with ViennaRNA CLI"
    ) -> pd.Series:
        """
        Helper method for pandas.apply() with progress tracking.
        
        Provides a pandas-compatible interface for ΔG calculations with
        optional progress display.
        
        Args:
            series: Series of DNA sequences
            description: Description for the progress bar
            
        Returns:
            Series of ΔG values with same index as input
            
        Raises:
            SequenceProcessingError: If series processing fails
            
        Example:
            >>> import pandas as pd
            >>> sequences = pd.Series(["ATCG", "GCTA", "TTTT"])
            >>> processor = ViennaRNAProcessor()
            >>> deltaG_values = processor.process_deltaG_series(sequences)
            >>> isinstance(deltaG_values, pd.Series)
            True
        """
        if not isinstance(series, pd.Series):
            error_msg = "Input must be a pandas Series"
            raise SequenceProcessingError(error_msg)
            
        logger.info(f"Processing {len(series)} sequences for ΔG with pandas using CLI")
        
        try:
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc=description)
                return series.progress_apply(cls.calc_deltaG)
            else:
                return series.apply(cls.calc_deltaG)
                
        except Exception as e:
            error_msg = f"Error processing ΔG series: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e

    @classmethod
    def validate_vienna_setup(cls) -> bool:
        """
        Validate that ViennaRNA CLI is properly installed and accessible.
        
        Checks ViennaRNA CLI installation, parameter files, and basic functionality
        to ensure thermodynamic calculations will work properly.
        
        Returns:
            True if ViennaRNA setup is valid, False otherwise
            
        Raises:
            ExternalToolError: If ViennaRNA validation fails
            
        Example:
            >>> processor = ViennaRNAProcessor()
            >>> if processor.validate_vienna_setup():
            ...     print("ViennaRNA ready for calculations")
        """
        logger.debug("Validating ViennaRNA CLI setup")
        
        try:
            # Test basic ViennaRNA CLI functionality
            test_sequence = "AUCG"
            
            try:
                # Test basic RNAfold command
                result = subprocess.run(
                    ["RNAfold", "--version"],
                    capture_output=True,
                    text=True,
                    check=True,
                    timeout=10
                )
                logger.debug(f"RNAfold version check successful: {result.stdout.strip()}")
                
            except subprocess.CalledProcessError as e:
                error_msg = f"RNAfold CLI not accessible: {str(e)}"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="RNAfold") from e
            except subprocess.TimeoutExpired:
                error_msg = "RNAfold CLI version check timed out"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="RNAfold")
            
            # Test actual folding
            try:
                energy = cls._run_rnafold_cli(test_sequence, use_dna_params=False)
                logger.debug(f"ViennaRNA CLI basic test successful: energy={energy}")
                
            except Exception as e:
                error_msg = f"ViennaRNA CLI basic functionality test failed: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                raise ExternalToolError(error_msg, tool_name="RNAfold") from e
            
            # Check for DNA parameter files
            dna_param_file = cls._find_dna_parameter_file()
            if dna_param_file:
                logger.debug(f"DNA parameters available: {dna_param_file}")
            else:
                logger.warning("No DNA parameter files found - using RNA parameters")
            
            logger.debug("ViennaRNA CLI setup validation successful")
            return True
            
        except ExternalToolError:
            # Re-raise without wrapping
            raise
        except Exception as e:
            error_msg = f"Unexpected error validating ViennaRNA CLI setup: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="RNAfold") from e

    @classmethod
    def get_vienna_info(cls) -> dict:
        """
        Get information about the ViennaRNA CLI installation and configuration.
        
        Returns:
            Dictionary containing ViennaRNA configuration details
            
        Example:
            >>> processor = ViennaRNAProcessor()
            >>> info = processor.get_vienna_info()
            >>> print(f"RNAfold available: {info.get('rnafold_available', False)}")
        """
        info = {
            'version': 'unknown',
            'rnafold_available': False,
            'rnafold_path': None,
            'dna_params_available': False,
            'dna_params_file': None,
            'temperature': Config.THERMO_TEMPERATURE,
            'sodium_concentration': Config.THERMO_SODIUM,
            'magnesium_concentration': Config.THERMO_MAGNESIUM,
            'interface': 'CLI'
        }
        
        try:
            # Check if RNAfold is available
            try:
                result = subprocess.run(
                    ["which", "RNAfold"],
                    capture_output=True,
                    text=True,
                    check=True
                )
                info['rnafold_available'] = True
                info['rnafold_path'] = result.stdout.strip()
            except subprocess.CalledProcessError:
                info['rnafold_available'] = False
            
            # Get ViennaRNA version if available
            if info['rnafold_available']:
                try:
                    result = subprocess.run(
                        ["RNAfold", "--version"],
                        capture_output=True,
                        text=True,
                        check=True,
                        timeout=5
                    )
                    # Parse version from output
                    version_line = result.stdout.split('\n')[0]
                    info['version'] = version_line.strip()
                except Exception:
                    pass
            
            # Check for DNA parameters
            dna_param_file = cls._find_dna_parameter_file()
            if dna_param_file:
                info['dna_params_available'] = True
                info['dna_params_file'] = str(dna_param_file)
            
            logger.debug(f"ViennaRNA CLI info: {info}")
            
        except Exception as e:
            logger.debug(f"Error getting ViennaRNA CLI info: {e}")
        
        return info