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

from ..config import Config
from ..config.exceptions import SequenceProcessingError

class ViennaRNAProcessor:
    """Handles thermodynamic calculations using ViennaRNA for DNA oligos."""

    # Set up module logger
    logger = logging.getLogger("ddPrimer.vienna_processor")
    
    @classmethod
    def _check_vienna_available(cls):
        """
        Check if ViennaRNA package is available.
        
        Returns:
            bool: True if ViennaRNA is available, False otherwise
        """
        try:
            import RNA
            return True
        except ImportError:
            cls.logger.warning("ViennaRNA package not available. Install with pip install ViennaRNA")
            return False

    @staticmethod
    def _make_fold_compound(seq):
        """
        Create a ViennaRNA fold_compound object loaded with DNA parameters.
        
        Args:
            seq (str): DNA sequence (A,C,G,T,N)
            
        Returns:
            RNA.fold_compound: Ready-to-use fold compound for MFE calculations
        """
        try:
            import RNA
            
            # ViennaRNA expects U instead of T even when using DNA parameters
            rna_seq = seq.upper().replace("T", "U")

            # Model details object lets us set temperature etc.
            md = RNA.md()
            md.temperature = Config.THERMO_TEMPERATURE  # Use generalized temperature parameter

            # Build fold compound
            fc = RNA.fold_compound(rna_seq, md, RNA.OPTION_DEFAULT)

            # Try to find DNA parameter file
            dna_param_file = None
            try:
                # Standard path in newer ViennaRNA installations
                dna_param_file = str(Path(RNA.params_path()) / "dna_mathews2004.par")
                if not Path(dna_param_file).exists():
                    # Alternative paths to check
                    alt_paths = [
                        Path(RNA.params_path()) / "dna_mathews1999.par",
                        Path(RNA.params_path()) / "dna_mathews.par",
                        Path("/usr/local/share/ViennaRNA/rna_turner2004.par"),  # Common system paths
                        Path("/usr/share/ViennaRNA/rna_turner2004.par")
                    ]
                    
                    for path in alt_paths:
                        if path.exists():
                            dna_param_file = str(path)
                            break
            except Exception as e:
                ViennaRNAProcessor.logger.debug(f"Error finding DNA parameter file: {e}")
                
            # Load DNA parameter set if found
            if dna_param_file and Path(dna_param_file).exists():
                try:
                    fc.params_load(dna_param_file)
                    ViennaRNAProcessor.logger.debug(f"Loaded DNA parameters from: {dna_param_file}")
                except Exception as e:
                    ViennaRNAProcessor.logger.debug(f"Error loading DNA parameters: {e}")

            # Set salt concentration if the API is available (ViennaRNA ≥ 2.6)
            if hasattr(fc, "params_set_salt"):
                try:
                    fc.params_set_salt(Config.THERMO_SODIUM)  # Use generalized sodium parameter
                except Exception as e:
                    ViennaRNAProcessor.logger.debug(f"Could not set salt concentration: {e}")
                    
            # Set magnesium concentration if available
            if hasattr(fc, "params_set_salt_MgdefaultK") and Config.THERMO_MAGNESIUM > 0:
                try:
                    # Convert from M to mM
                    fc.params_set_salt_MgdefaultK(Config.THERMO_MAGNESIUM * 1000)  # Use generalized magnesium parameter
                except Exception as e:
                    ViennaRNAProcessor.logger.debug(f"Could not set magnesium concentration: {e}")

            return fc
            
        except ImportError:
            ViennaRNAProcessor.logger.error("ViennaRNA package not available")
            return None
        except Exception as e:
            ViennaRNAProcessor.logger.error(f"Error creating fold compound: {e}")
            return None

    @classmethod
    def calc_deltaG(cls, seq):
        """
        Calculate the minimum free energy (ΔG, kcal/mol) of a DNA oligo.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: Minimum free energy, or None if calculation fails
        """
        if not cls._check_vienna_available():
            return None
            
        if not isinstance(seq, str) or seq == "":
            return None

        # Validate DNA sequence format
        dna_pattern = re.compile(r"^[ACGTNacgtn]+$")
        if not dna_pattern.match(seq):
            cls.logger.debug(f"Invalid DNA sequence format: {seq[:50]}")
            return None

        try:
            import RNA
            
            # Create fold compound and calculate MFE
            fc = cls._make_fold_compound(seq)
            if fc is None:
                return None
                
            structure, energy = fc.mfe()
            return energy  # kcal/mol
            
        except ImportError:
            cls.logger.error("ViennaRNA package not available")
            return None
        except Exception as e:
            if Config.SHOW_PROGRESS:
                preview = seq[:50] + ("..." if len(seq) > 50 else "")
                cls.logger.warning(f"ViennaRNA error for sequence {preview}: {e}")
            cls.logger.debug(f"ViennaRNA error details: {str(e)}", exc_info=True)
            return None

    @classmethod
    def calc_deltaG_batch(cls, seqs, description="Calculating ΔG with ViennaRNA"):
        """
        Calculate ΔG for a batch of sequences with progress bar.
        
        Args:
            seqs (list): List of DNA sequences
            description (str, optional): Description for the progress bar
            
        Returns:
            list: List of ΔG values
        """
        if not cls._check_vienna_available():
            return [None] * len(seqs)
            
        cls.logger.info(f"Processing batch of {len(seqs)} sequences for ΔG calculation")
        results = []
        
        # Apply tqdm for progress tracking if enabled
        if Config.SHOW_PROGRESS:
            sequence_iter = tqdm(seqs, desc=description)
        else:
            sequence_iter = seqs
            
        for seq in sequence_iter:
            if pd.notnull(seq):
                results.append(cls.calc_deltaG(seq))
            else:
                results.append(None)
                
        cls.logger.debug(f"Completed ΔG calculations for {len(seqs)} sequences")
        return results
        
    @classmethod
    def process_deltaG(cls, series, description="Processing sequences with ViennaRNA"):
        """
        Helper method for pandas.apply() with progress tracking.
        
        Args:
            series (pandas.Series): Series of DNA sequences
            description (str, optional): Description for the progress bar
            
        Returns:
            pandas.Series: Series of ΔG values
        """
        if not cls._check_vienna_available():
            return pd.Series([None] * len(series), index=series.index)
            
        cls.logger.info(f"Processing {len(series)} sequences for ΔG with pandas")
        
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc=description)
            return series.progress_apply(cls.calc_deltaG)
        else:
            return series.apply(cls.calc_deltaG)