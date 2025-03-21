import re
import pandas as pd
import nupack
from tqdm import tqdm
from ..config import Config

class NupackProcessor:
    """Handles thermodynamic calculations using NUPACK."""
    
    @classmethod
    def calc_deltaG_batch(cls, sequences, description="Calculating ΔG"):
        """
        Calculate deltaG for a batch of sequences with progress bar.
        
        Args:
            sequences (list): List of DNA sequences
            description (str): Description for the progress bar
            
        Returns:
            list: List of deltaG values
        """
        results = []
        
        if Config.SHOW_PROGRESS:
            sequence_iter = tqdm(sequences, desc=description)
        else:
            sequence_iter = sequences
            
        for seq in sequence_iter:
            if seq is not None and isinstance(seq, str) and seq.strip():
                results.append(cls.calc_deltaG(seq))
            else:
                results.append(None)
                
        return results
    
    @staticmethod
    def calc_deltaG(seq):
        """
        Calculate the minimum free energy (ΔG) of a DNA sequence using NUPACK.
        
        Args:
            seq (str): DNA sequence
            
        Returns:
            float: Minimum free energy in kcal/mol, or None on error
        """
        if not seq or not isinstance(seq, str) or not seq.strip():
            return None
            
        try:
            # Set up NUPACK model with our parameters
            model = nupack.Model(
                material='dna',
                temperature=Config.NUPACK_TEMPERATURE + 273.15,  # Kelvin
                sodium=Config.NUPACK_SODIUM,
                magnesium=Config.NUPACK_MAGNESIUM
            )
            
            # Calculate and return minimum free energy
            mfe = nupack.mfe([seq], model=model)
            return float(mfe[0].energy)
        except Exception as e:
            print(f"Error calculating deltaG for sequence {seq}: {e}")
            return None