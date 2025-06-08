#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thermodynamic calculation module using ViennaRNA for DNA oligos.

Provides thermodynamic calculations for DNA sequences using ViennaRNA
with DNA-specific parameters and salt corrections.
"""

import re
import logging
from pathlib import Path
from typing import List, Optional, Union
import pandas as pd
from tqdm import tqdm
import RNA

from ..config import Config, SequenceProcessingError, ExternalToolError

# Type alias for path inputs
PathLike = Union[str, Path]


class ViennaRNAProcessor:
    """
    Handles thermodynamic calculations using ViennaRNA for DNA oligos.
    
    This class provides methods for calculating minimum free energy (ΔG)
    of DNA sequences using ViennaRNA with DNA-specific parameters.
    
    Attributes:
        logger: Logger instance for this processor
        
    Example:
        >>> processor = ViennaRNAProcessor()
        >>> deltaG = processor.calc_deltaG("ATCGATCG")
        >>> batch_results = processor.calc_deltaG_batch(sequences)
    """

    # Set up module logger
    logger = logging.getLogger("ddPrimer.vienna_processor")

    @classmethod
    def _make_fold_compound(cls, seq: str) -> RNA.fold_compound:
        """
        Create a ViennaRNA fold_compound object loaded with DNA parameters.
        
        Creates a fold compound configured for DNA thermodynamics with
        appropriate temperature and salt conditions.
        
        Args:
            seq: DNA sequence (A,C,G,T,N)
            
        Returns:
            Ready-to-use fold compound for MFE calculations
            
        Raises:
            SequenceProcessingError: If fold compound creation fails
            ExternalToolError: If ViennaRNA parameters cannot be loaded
            
        Example:
            >>> fc = ViennaRNAProcessor._make_fold_compound("ATCG")
            >>> structure, energy = fc.mfe()
        """
        if not seq or not isinstance(seq, str):
            raise SequenceProcessingError("Invalid sequence provided for fold compound creation")
            
        try:
            # ViennaRNA expects U instead of T even when using DNA parameters
            rna_seq = seq.upper().replace("T", "U")
            
            # Validate converted sequence
            if not cls._is_valid_rna_sequence(rna_seq):
                raise SequenceProcessingError(f"Invalid sequence after T->U conversion: {seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...")

            # Model details object lets us set temperature etc.
            md = RNA.md()
            md.temperature = Config.THERMO_TEMPERATURE
            cls.logger.debug(f"Set ViennaRNA temperature: {Config.THERMO_TEMPERATURE}°C")

            # Build fold compound
            fc = RNA.fold_compound(rna_seq, md, RNA.OPTION_DEFAULT)

            # Load DNA parameter file
            dna_param_file = cls._find_dna_parameter_file()
            if dna_param_file:
                try:
                    fc.params_load(str(dna_param_file))
                    cls.logger.debug(f"Loaded DNA parameters from: {dna_param_file}")
                except Exception as e:
                    cls.logger.warning(f"Could not load DNA parameters from {dna_param_file}: {e}")
                    raise ExternalToolError(
                        f"Failed to load ViennaRNA DNA parameters: {e}",
                        tool_name="ViennaRNA"
                    ) from e

            # Set salt concentrations if the API is available
            cls._set_salt_concentrations(fc)

            return fc
            
        except SequenceProcessingError:
            # Re-raise without wrapping
            raise
        except ExternalToolError:
            # Re-raise without wrapping
            raise
        except Exception as e:
            error_msg = f"Error creating fold compound for sequence '{seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...': {str(e)}"
            cls.logger.error(error_msg)
            cls.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e

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
        try:
            # Get base path from ViennaRNA
            try:
                base_path = Path(RNA.params_path())
                cls.logger.debug(f"ViennaRNA params path: {base_path}")
            except Exception as e:
                cls.logger.debug(f"Could not get ViennaRNA params path: {e}")
                base_path = None
            
            # Try different DNA parameter files in order of preference
            candidate_files = []
            
            # Add base path candidates if available
            if base_path:
                candidate_files.extend([
                    base_path / "dna_mathews2004.par",
                    base_path / "dna_mathews1999.par", 
                    base_path / "dna_mathews.par"
                ])
            
            # Add system-wide installation paths
            candidate_files.extend([
                Path("/usr/local/share/ViennaRNA/dna_mathews2004.par"),
                Path("/usr/share/ViennaRNA/dna_mathews2004.par"),
                Path("/opt/ViennaRNA/share/ViennaRNA/dna_mathews2004.par")
            ])
            
            for param_file in candidate_files:
                if param_file.exists() and param_file.is_file():
                    cls.logger.debug(f"Found DNA parameter file: {param_file}")
                    return param_file
                    
            cls.logger.debug("No DNA parameter file found, using default RNA parameters")
            return None
            
        except Exception as e:
            cls.logger.debug(f"Error finding DNA parameter file: {e}")
            return None

    @classmethod
    def _set_salt_concentrations(cls, fc: RNA.fold_compound) -> None:
        """
        Set salt concentrations on the fold compound if the API supports it.
        
        Configures sodium and magnesium concentrations for more accurate
        thermodynamic calculations.
        
        Args:
            fc: ViennaRNA fold compound object
            
        Example:
            >>> fc = RNA.fold_compound("AUCG")
            >>> ViennaRNAProcessor._set_salt_concentrations(fc)
        """
        # Set sodium concentration if the API is available (ViennaRNA ≥ 2.6)
        if hasattr(fc, "params_set_salt"):
            try:
                fc.params_set_salt(Config.THERMO_SODIUM)
                cls.logger.debug(f"Set sodium concentration: {Config.THERMO_SODIUM} M")
            except Exception as e:
                cls.logger.debug(f"Could not set sodium concentration: {e}")
                
        # Set magnesium concentration if available and configured
        if hasattr(fc, "params_set_salt_MgdefaultK") and Config.THERMO_MAGNESIUM > 0:
            try:
                # Convert from M to mM
                mg_concentration_mm = Config.THERMO_MAGNESIUM * Config.MOLARITY_TO_MILLIMOLAR
                fc.params_set_salt_MgdefaultK(mg_concentration_mm)
                cls.logger.debug(f"Set magnesium concentration: {mg_concentration_mm} mM")
            except Exception as e:
                cls.logger.debug(f"Could not set magnesium concentration: {e}")

    @classmethod
    def calc_deltaG(cls, seq: str) -> Optional[float]:
        """
        Calculate the minimum free energy (ΔG, kcal/mol) of a DNA oligo.
        
        Computes the minimum free energy for DNA secondary structure formation
        using ViennaRNA with DNA-specific parameters and salt corrections.
        
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
            cls.logger.debug("Empty or invalid sequence provided")
            return None

        # Validate DNA sequence format
        if not cls._is_valid_dna_sequence(seq):
            cls.logger.debug(f"Invalid DNA sequence format: {seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...")
            raise SequenceProcessingError(f"Invalid DNA sequence: contains non-DNA characters")

        try:
            fc = cls._make_fold_compound(seq)
            structure, energy = fc.mfe()
            
            cls.logger.debug(f"ΔG calculation successful: {energy:.3f} kcal/mol for {len(seq)} bp sequence")
            return energy  # kcal/mol
            
        except SequenceProcessingError:
            # Already logged in _make_fold_compound
            raise
        except ExternalToolError:
            # Already logged in _make_fold_compound
            raise
        except Exception as e:
            error_msg = f"ViennaRNA calculation failed for sequence {seq[:Config.MAX_SEQUENCE_DISPLAY_LENGTH]}...: {str(e)}"
            cls.logger.warning(error_msg)
            cls.logger.debug(f"ViennaRNA error details: {str(e)}", exc_info=True)
            return None

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
    def calc_deltaG_batch(
        cls, 
        seqs: List[str], 
        description: str = "Calculating ΔG with ViennaRNA"
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
            raise SequenceProcessingError("Sequences must be provided as a list")
            
        cls.logger.info(f"Processing batch of {len(seqs)} sequences for ΔG calculation")
        cls.logger.debug(f"Progress tracking: {Config.SHOW_PROGRESS}")
        
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
                        cls.logger.debug(f"ΔG calculation failed for sequence {seq_num}: {e}")
                        results.append(None)
                        failed_count += 1
                else:
                    results.append(None)
                    
            cls.logger.debug(f"Completed ΔG calculations: {len(results)} total, {failed_count} failed")
            return results
            
        except Exception as e:
            error_msg = f"Error during batch ΔG calculation: {str(e)}"
            cls.logger.error(error_msg)
            raise SequenceProcessingError(error_msg) from e
        
    @classmethod
    def process_deltaG_series(
        cls, 
        series: pd.Series, 
        description: str = "Processing sequences with ViennaRNA"
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
            
        Example:
            >>> import pandas as pd
            >>> sequences = pd.Series(["ATCG", "GCTA", "TTTT"])
            >>> processor = ViennaRNAProcessor()
            >>> deltaG_values = processor.process_deltaG_series(sequences)
            >>> isinstance(deltaG_values, pd.Series)
            True
        """
        if not isinstance(series, pd.Series):
            raise SequenceProcessingError("Input must be a pandas Series")
            
        cls.logger.info(f"Processing {len(series)} sequences for ΔG with pandas")
        
        try:
            if Config.SHOW_PROGRESS:
                tqdm.pandas(desc=description)
                return series.progress_apply(cls.calc_deltaG)
            else:
                return series.apply(cls.calc_deltaG)
                
        except Exception as e:
            error_msg = f"Error processing ΔG series: {str(e)}"
            cls.logger.error(error_msg)
            raise SequenceProcessingError(error_msg) from e

    @classmethod
    def validate_vienna_setup(cls) -> bool:
        """
        Validate that ViennaRNA is properly installed and accessible.
        
        Checks ViennaRNA installation, parameter files, and basic functionality
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
        cls.logger.debug("Validating ViennaRNA setup")
        
        try:
            # Test basic ViennaRNA functionality
            test_sequence = "AUCG"
            
            try:
                # Create a simple fold compound
                md = RNA.md()
                fc = RNA.fold_compound(test_sequence, md, RNA.OPTION_DEFAULT)
                structure, energy = fc.mfe()
                
                cls.logger.debug(f"ViennaRNA basic test successful: energy={energy}")
                
            except Exception as e:
                error_msg = f"ViennaRNA basic functionality test failed: {str(e)}"
                cls.logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="ViennaRNA") from e
            
            # Check for DNA parameter files
            dna_param_file = cls._find_dna_parameter_file()
            if dna_param_file:
                cls.logger.debug(f"DNA parameters available: {dna_param_file}")
            else:
                cls.logger.warning("No DNA parameter files found - using RNA parameters")
            
            # Test temperature setting
            try:
                md = RNA.md()
                md.temperature = Config.THERMO_TEMPERATURE
                cls.logger.debug(f"Temperature setting successful: {Config.THERMO_TEMPERATURE}°C")
            except Exception as e:
                cls.logger.warning(f"Could not set temperature: {e}")
            
            cls.logger.debug("ViennaRNA setup validation successful")
            return True
            
        except ExternalToolError:
            # Re-raise without wrapping
            raise
        except Exception as e:
            error_msg = f"Unexpected error validating ViennaRNA setup: {str(e)}"
            cls.logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="ViennaRNA") from e

    @classmethod
    def get_vienna_info(cls) -> dict:
        """
        Get information about the ViennaRNA installation and configuration.
        
        Returns:
            Dictionary containing ViennaRNA configuration details
            
        Example:
            >>> processor = ViennaRNAProcessor()
            >>> info = processor.get_vienna_info()
            >>> print(f"ViennaRNA version: {info.get('version', 'unknown')}")
        """
        info = {
            'version': 'unknown',
            'params_path': None,
            'dna_params_available': False,
            'dna_params_file': None,
            'temperature': Config.THERMO_TEMPERATURE,
            'sodium_concentration': Config.THERMO_SODIUM,
            'magnesium_concentration': Config.THERMO_MAGNESIUM
        }
        
        try:
            # Get ViennaRNA version if available
            if hasattr(RNA, '__version__'):
                info['version'] = RNA.__version__
            elif hasattr(RNA, 'version'):
                info['version'] = RNA.version()
            
            # Get parameters path
            try:
                info['params_path'] = RNA.params_path()
            except Exception:
                pass
            
            # Check for DNA parameters
            dna_param_file = cls._find_dna_parameter_file()
            if dna_param_file:
                info['dna_params_available'] = True
                info['dna_params_file'] = str(dna_param_file)
            
            cls.logger.debug(f"ViennaRNA info: {info}")
            
        except Exception as e:
            cls.logger.debug(f"Error getting ViennaRNA info: {e}")
        
        return info