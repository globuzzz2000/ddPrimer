import re
import pandas as pd
import nupack
from tqdm import tqdm
from ..config import Config

class NupackProcessor:
    """Handles thermodynamic calculations using NUPACK."""
    
    @staticmethod
    def calc_deltaG(seq):
        """
        Calculate the minimum free energy of a sequence using NUPACK.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: Minimum free energy, or None if calculation fails
        """
        if not isinstance(seq, str) or seq == "":
            return None
            
        dna_pattern = re.compile(r'^[ACGTNacgtn]+$')
        if not dna_pattern.match(seq):
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
                print(f"NUPACK error for sequence {preview}: {e}")
            return None
            
        return None
        
    @classmethod
    def calc_deltaG_batch(cls, seqs, description="Calculating ΔG"):
        """
        Calculate ΔG for a batch of sequences with progress bar.
        
        Args:
            seqs (list): List of DNA sequences
            description (str): Description for the progress bar
            
        Returns:
            list: List of ΔG values
        """
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
                
        return results
        
    @classmethod
    def pandas_calc_deltaG(cls, series, description="Processing sequences"):
        """
        Helper method to use with pandas.apply() that handles progress tracking.
        
        Args:
            series (pandas.Series): Series of DNA sequences
            description (str): Description for the progress bar
            
        Returns:
            pandas.Series: Series of deltaG values
        """
        if Config.SHOW_PROGRESS:
            tqdm.pandas(desc=description)
            return series.progress_apply(cls.calc_deltaG)
        else:
            return series.apply(cls.calc_deltaG)