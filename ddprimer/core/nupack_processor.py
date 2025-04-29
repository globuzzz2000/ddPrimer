import re
import logging
import pandas as pd
import nupack
from tqdm import tqdm
from ..config import Config
from ..config.exceptions import SequenceProcessingError

class NupackProcessor:
    """Handles thermodynamic calculations using NUPACK."""
    
    # Get module logger
    logger = logging.getLogger("ddPrimer.nupack_processor")
    
    @staticmethod
    def calc_deltaG(seq):
        """
        Calculate the minimum free energy of a sequence using NUPACK.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: Minimum free energy, or None if calculation fails
            
        Raises:
            SequenceProcessingError: When sequence format is invalid
        """
        if not isinstance(seq, str) or seq == "":
            return None
            
        dna_pattern = re.compile(r'^[ACGTNacgtn]+$')
        if not dna_pattern.match(seq):
            NupackProcessor.logger.debug(f"Invalid DNA sequence format: {seq[:50]}")
            return None
            
        try:
            model = nupack.Model(
                material='dna',
                celsius=Config.NUPACK_TEMPERATURE,
                sodium=Config.NUPACK_SODIUM,
                magnesium=Config.NUPACK_MAGNESIUM
            )
            result = nupack.mfe(seq, model=model)
            if result:
                # result[0].energy is the MFE structure's free energy
                return result[0].energy
        except Exception as e:
            if Config.SHOW_PROGRESS:
                preview = seq[:50] + ("..." if len(seq) > 50 else "")
                NupackProcessor.logger.warning(f"NUPACK error for sequence {preview}: {e}")
            NupackProcessor.logger.debug(f"NUPACK error details: {str(e)}", exc_info=True)
            return None
            
        return None
        
    @classmethod
    def calc_deltaG_batch(cls, seqs, description="Calculating ΔG"):
        """
        Calculate ΔG for a batch of sequences with progress bar.
        
        Args:
            seqs (list): List of DNA sequences
            description (str, optional): Description for the progress bar. Defaults to "Calculating ΔG".
            
        Returns:
            list: List of ΔG values
        """
        cls.logger.info(f"Processing batch of {len(seqs)} sequences for ΔG calculation")
        results = []
        
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
    def process_deltaG(cls, series, description="Processing sequences"):
        """
        Helper method to use with pandas.apply() that handles progress tracking.
        
        Args:
            series (pandas.Series): Series of DNA sequences
            description (str, optional): Description for the progress bar. Defaults to "Processing sequences".
            
        Returns:
            pandas.Series: Series of deltaG values
        """
        cls.logger.info(f"Processing {len(series)} sequences for ΔG with pandas")
        
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc=description)
            return series.progress_apply(cls.calc_deltaG)
        else:
            return series.apply(cls.calc_deltaG)