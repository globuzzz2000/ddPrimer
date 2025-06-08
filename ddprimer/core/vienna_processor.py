#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thermodynamic calculation module using ViennaRNA for DNA oligos.
"""

import re
import logging
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import RNA

from ..config import Config, SequenceProcessingError

class ViennaRNAProcessor:
    """Handles thermodynamic calculations using ViennaRNA for DNA oligos."""

    # Set up module logger
    logger = logging.getLogger("ddPrimer.vienna_processor")

    @staticmethod
    def _make_fold_compound(seq):
        """
        Create a ViennaRNA fold_compound object loaded with DNA parameters.
        
        Args:
            seq (str): DNA sequence (A,C,G,T,N)
            
        Returns:
            RNA.fold_compound: Ready-to-use fold compound for MFE calculations
            
        Raises:
            SequenceProcessingError: If fold compound creation fails
        """
        try:
            # ViennaRNA expects U instead of T even when using DNA parameters
            rna_seq = seq.upper().replace("T", "U")

            # Model details object lets us set temperature etc.
            md = RNA.md()
            md.temperature = Config.THERMO_TEMPERATURE

            # Build fold compound
            fc = RNA.fold_compound(rna_seq, md, RNA.OPTION_DEFAULT)

            # Load DNA parameter file
            dna_param_file = ViennaRNAProcessor._find_dna_parameter_file()
            if dna_param_file:
                try:
                    fc.params_load(str(dna_param_file))
                    ViennaRNAProcessor.logger.debug(f"Loaded DNA parameters from: {dna_param_file}")
                except Exception as e:
                    ViennaRNAProcessor.logger.warning(f"Could not load DNA parameters from {dna_param_file}: {e}")

            # Set salt concentrations if the API is available
            ViennaRNAProcessor._set_salt_concentrations(fc)

            return fc
            
        except Exception as e:
            ViennaRNAProcessor.logger.error(f"Error creating fold compound for sequence '{seq[:50]}...': {e}")
            raise SequenceProcessingError(f"Failed to create ViennaRNA fold compound: {e}")

    @staticmethod
    def _find_dna_parameter_file():
        """
        Find the DNA parameter file in standard ViennaRNA locations.
        
        Returns:
            Path or None: Path to DNA parameter file if found
        """
        try:
            base_path = Path(RNA.params_path())
            
            # Try different DNA parameter files in order of preference
            candidate_files = [
                base_path / "dna_mathews2004.par",
                base_path / "dna_mathews1999.par", 
                base_path / "dna_mathews.par",
                Path("/usr/local/share/ViennaRNA/dna_mathews2004.par"),
                Path("/usr/share/ViennaRNA/dna_mathews2004.par")
            ]
            
            for param_file in candidate_files:
                if param_file.exists():
                    return param_file
                    
            ViennaRNAProcessor.logger.debug("No DNA parameter file found, using default RNA parameters")
            return None
            
        except Exception as e:
            ViennaRNAProcessor.logger.debug(f"Error finding DNA parameter file: {e}")
            return None

    @staticmethod
    def _set_salt_concentrations(fc):
        """
        Set salt concentrations on the fold compound if the API supports it.
        
        Args:
            fc: ViennaRNA fold compound object
        """
        # Set sodium concentration if the API is available (ViennaRNA ≥ 2.6)
        if hasattr(fc, "params_set_salt"):
            try:
                fc.params_set_salt(Config.THERMO_SODIUM)
                ViennaRNAProcessor.logger.debug(f"Set sodium concentration: {Config.THERMO_SODIUM} M")
            except Exception as e:
                ViennaRNAProcessor.logger.debug(f"Could not set sodium concentration: {e}")
                
        # Set magnesium concentration if available and configured
        if hasattr(fc, "params_set_salt_MgdefaultK") and Config.THERMO_MAGNESIUM > 0:
            try:
                # Convert from M to mM
                mg_concentration_mm = Config.THERMO_MAGNESIUM * 1000
                fc.params_set_salt_MgdefaultK(mg_concentration_mm)
                ViennaRNAProcessor.logger.debug(f"Set magnesium concentration: {mg_concentration_mm} mM")
            except Exception as e:
                ViennaRNAProcessor.logger.debug(f"Could not set magnesium concentration: {e}")

    @classmethod
    def calc_deltaG(cls, seq):
        """
        Calculate the minimum free energy (ΔG, kcal/mol) of a DNA oligo.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float or None: Minimum free energy in kcal/mol, or None if calculation fails
        """
        if not isinstance(seq, str) or seq == "":
            return None

        # Validate DNA sequence format
        if not cls._is_valid_dna_sequence(seq):
            cls.logger.debug(f"Invalid DNA sequence format: {seq[:50]}...")
            return None

        try:
            fc = cls._make_fold_compound(seq)
            structure, energy = fc.mfe()
            return energy  # kcal/mol
            
        except SequenceProcessingError:
            # Already logged in _make_fold_compound
            return None
        except Exception as e:
            cls.logger.warning(f"ViennaRNA calculation failed for sequence {seq[:50]}...: {e}")
            cls.logger.debug(f"ViennaRNA error details: {str(e)}", exc_info=True)
            return None

    @staticmethod
    def _is_valid_dna_sequence(seq):
        """
        Validate that a sequence contains only valid DNA characters.
        
        Args:
            seq (str): Sequence to validate
            
        Returns:
            bool: True if sequence is valid DNA
        """
        dna_pattern = re.compile(r"^[ACGTNacgtn]+$")
        return bool(dna_pattern.match(seq))

    @classmethod
    def calc_deltaG_batch(cls, seqs, description="Calculating ΔG with ViennaRNA"):
        """
        Calculate ΔG for a batch of sequences with progress tracking.
        
        Args:
            seqs (list): List of DNA sequences
            description (str, optional): Description for the progress bar
            
        Returns:
            list: List of ΔG values (same length as input, with None for failed calculations)
        """
        cls.logger.info(f"Processing batch of {len(seqs)} sequences for ΔG calculation")
        results = []
        
        # Apply progress tracking if enabled
        sequence_iter = tqdm(seqs, desc=description) if Config.SHOW_PROGRESS else seqs
            
        for seq in sequence_iter:
            if pd.notnull(seq) and seq:
                results.append(cls.calc_deltaG(seq))
            else:
                results.append(None)
                
        cls.logger.debug(f"Completed ΔG calculations for {len(seqs)} sequences")
        return results
        
    @classmethod
    def process_deltaG_series(cls, series, description="Processing sequences with ViennaRNA"):
        """
        Helper method for pandas.apply() with progress tracking.
        
        Args:
            series (pandas.Series): Series of DNA sequences
            description (str, optional): Description for the progress bar
            
        Returns:
            pandas.Series: Series of ΔG values
        """
        cls.logger.info(f"Processing {len(series)} sequences for ΔG with pandas")
        
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc=description)
            return series.progress_apply(cls.calc_deltaG)
        else:
            return series.apply(cls.calc_deltaG)