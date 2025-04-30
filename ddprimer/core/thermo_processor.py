#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flexible thermodynamic processing module that can switch between NUPACK
and ViennaRNA packages for calculating thermodynamic properties.
"""

import logging
import pandas as pd
from tqdm import tqdm

from ..config import Config

# Initialize logger
logger = logging.getLogger("ddPrimer.thermo_processor")

class ThermoProcessor:
    """
    Flexible thermodynamic processor that can use either NUPACK or ViennaRNA
    based on availability and configuration.
    """
    
    # Initialize backend preference
    _backend = None
    
    @classmethod
    def check_nupack_available(cls):
        """
        Check if NUPACK is available.
        
        Returns:
            bool: True if NUPACK is available, False otherwise
        """
        # Use the is_available method from NupackProcessor
        try:
            from ..core.nupack_processor import NupackProcessor
            return NupackProcessor.is_available()
        except ImportError:
            logger.debug("Could not import NupackProcessor")
            return False
    
    @classmethod
    def check_vienna_available(cls):
        """
        Check if ViennaRNA is available.
        
        Returns:
            bool: True if ViennaRNA is available, False otherwise
        """
        # Use the is_available method from ViennaRNAProcessor
        try:
            from ..core.vienna_processor import ViennaRNAProcessor
            return ViennaRNAProcessor.is_available()
        except ImportError:
            logger.debug("Could not import ViennaRNAProcessor")
            return False
    
    @classmethod
    def get_backend(cls):
        """
        Determine which thermodynamics backend to use based on availability and configuration.
        
        Returns:
            str: 'nupack', 'vienna', or None if neither is available
        """
        # If already determined, return the cached result
        if cls._backend is not None:
            return cls._backend
            
        # Try to use the configured backend preference first
        preferred_backend = getattr(Config, 'THERMO_BACKEND', 'nupack').lower()
        auto_fallback = getattr(Config, 'THERMO_AUTO_FALLBACK', True)
        
        if preferred_backend == 'nupack' and cls.check_nupack_available():
            cls._backend = 'nupack'
            logger.info("Using NUPACK for thermodynamic calculations")
        elif preferred_backend == 'vienna' and cls.check_vienna_available():
            cls._backend = 'vienna'
            logger.info("Using ViennaRNA for thermodynamic calculations")
        # Fall back to whatever is available if auto fallback is enabled
        elif auto_fallback and cls.check_nupack_available():
            cls._backend = 'nupack'
            logger.info("Falling back to NUPACK for thermodynamic calculations")
        elif auto_fallback and cls.check_vienna_available():
            cls._backend = 'vienna'
            logger.info("Falling back to ViennaRNA for thermodynamic calculations")
        else:
            cls._backend = None
            logger.warning("No thermodynamic calculation package available. "
                         "Install either NUPACK or ViennaRNA.")
            
        return cls._backend
    
    @classmethod
    def set_backend(cls, backend):
        """
        Explicitly set the thermodynamics backend.
        
        Args:
            backend (str): 'nupack' or 'vienna'
            
        Returns:
            bool: True if successful, False otherwise
        """
        if backend.lower() == 'nupack':
            if cls.check_nupack_available():
                cls._backend = 'nupack'
                logger.info("Switched to NUPACK for thermodynamic calculations")
                return True
            else:
                logger.warning("Cannot switch to NUPACK as it is not available")
                return False
        elif backend.lower() == 'vienna':
            if cls.check_vienna_available():
                cls._backend = 'vienna'
                logger.info("Switched to ViennaRNA for thermodynamic calculations")
                return True
            else:
                logger.warning("Cannot switch to ViennaRNA as it is not available")
                return False
        else:
            logger.warning(f"Unknown backend: {backend}. Use 'nupack' or 'vienna'")
            return False
    
    @classmethod
    def calc_deltaG(cls, seq):
        """
        Calculate the minimum free energy of a sequence using the current backend.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: Minimum free energy, or None if calculation fails
        """
        # Get the appropriate backend
        backend = cls.get_backend()
        
        if backend is None:
            logger.debug("No thermodynamic backend available")
            return None
            
        # Calculate using the selected backend
        try:
            if backend == 'nupack':
                from ..core.nupack_processor import NupackProcessor
                return NupackProcessor.calc_deltaG(seq)
            else:  # backend == 'vienna'
                from ..core.vienna_processor import ViennaRNAProcessor
                return ViennaRNAProcessor.calc_deltaG(seq)
        except ImportError as e:
            logger.warning(f"Could not import thermodynamic processor: {e}")
            return None
    
    @classmethod
    def calc_deltaG_batch(cls, seqs, description=None):
        """
        Calculate ΔG for a batch of sequences with progress bar.
        
        Args:
            seqs (list): List of DNA sequences
            description (str, optional): Description for the progress bar
            
        Returns:
            list: List of ΔG values
        """
        # Get the appropriate backend
        backend = cls.get_backend()
        
        if backend is None:
            logger.debug("No thermodynamic backend available")
            return [None] * len(seqs)
            
        # Use appropriate description if not provided
        if description is None:
            backend_name = "NUPACK" if backend == 'nupack' else "ViennaRNA"
            description = f"Calculating ΔG with {backend_name}"
            
        # Calculate using the selected backend
        try:
            if backend == 'nupack':
                from ..core.nupack_processor import NupackProcessor
                return NupackProcessor.calc_deltaG_batch(seqs, description)
            else:  # backend == 'vienna'
                from ..core.vienna_processor import ViennaRNAProcessor
                return ViennaRNAProcessor.calc_deltaG_batch(seqs, description)
        except ImportError as e:
            logger.warning(f"Could not import thermodynamic processor: {e}")
            return [None] * len(seqs)
    
    @classmethod
    def process_deltaG(cls, series, description=None):
        """
        Helper method for pandas.apply() with progress tracking.
        
        Args:
            series (pandas.Series): Series of DNA sequences
            description (str, optional): Description for the progress bar
            
        Returns:
            pandas.Series: Series of ΔG values
        """
        # Get the appropriate backend
        backend = cls.get_backend()
        
        if backend is None:
            logger.debug("No thermodynamic backend available")
            return pd.Series([None] * len(series), index=series.index)
            
        # Use appropriate description if not provided
        if description is None:
            backend_name = "NUPACK" if backend == 'nupack' else "ViennaRNA"
            description = f"Processing sequences with {backend_name}"
            
        # Calculate using the selected backend
        try:
            if backend == 'nupack':
                from ..core.nupack_processor import NupackProcessor
                return NupackProcessor.process_deltaG(series, description)
            else:  # backend == 'vienna'
                from ..core.vienna_processor import ViennaRNAProcessor
                return ViennaRNAProcessor.process_deltaG(series, description)
        except ImportError as e:
            logger.warning(f"Could not import thermodynamic processor: {e}")
            return pd.Series([None] * len(series), index=series.index)