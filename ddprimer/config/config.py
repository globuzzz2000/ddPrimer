#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration module for the ddPrimer pipeline.
"""

import os
import json
from multiprocessing import cpu_count
from typing import Dict, List, Union, Any, Optional

class Config:
    """Central configuration settings for the ddPrimer pipeline with singleton pattern."""
    
    # Singleton instance
    _instance = None
    
    #############################################################################
    #                           Pipeline Mode Options
    #############################################################################
    DEBUG_MODE = False                   # Debug logging mode (enable with --debug flag)
    DISABLE_INTERNAL_OLIGO = False      # Disable internal oligo design
    
    #############################################################################
    #                           Performance Settings
    #############################################################################
    NUM_PROCESSES = max(1, int(cpu_count() * 0.75))  # Use 75% of cores
    BATCH_SIZE = 100
    MAF_CHUNK_SIZE = 10000
    SHOW_PROGRESS = True
    
    #############################################################################
    #                           Design Parameters
    #############################################################################
    # Basic primer design constraints
    PRIMER_MIN_SIZE = 18
    PRIMER_OPT_SIZE = 20
    PRIMER_MAX_SIZE = 23
    PRIMER_MIN_TM = 50.0
    PRIMER_OPT_TM = 57.5
    PRIMER_MAX_TM = 65.0
    PRIMER_MIN_GC = 50.0
    PRIMER_MAX_GC = 60.0
    
    # Product size constraints
    PRIMER_PRODUCT_SIZE_RANGE = [[90, 200]]
    
    # Pipeline parameters
    MIN_SEGMENT_LENGTH = 90  # Match first entry in PRIMER_PRODUCT_SIZE_RANGE
    RETAIN_TYPES = "gene"  # gff filtering: "gene", "mRNA", "CDS", "exon", etc.
    FILTER_MEANINGFUL_NAMES = True  # only use named genes from gff
    COUNT_AMBIGUOUS_AS_MISMATCH = False
    GENE_OVERLAP_MARGIN = 25
    RESTRICTION_SITE = "GGCC"
    PENALTY_MAX = 5.0
    MAX_PRIMER_PAIRS_PER_SEGMENT = 3
    PREFER_PROBE_MORE_C_THAN_G = True  # Set to False to disable
    SEQUENCE_MIN_GC = 50.0
    SEQUENCE_MAX_GC = 60.0
    
    #############################################################################
    #                           Validation Options
    #############################################################################
    VALIDATION_MODE = "TOLERANT"  # "STRICT" or "TOLERANT"
    ALLOW_AMP_MISMATCHES = 2  # Number of mismatches allowed in amplicon for TOLERANT mode
    ALLOW_AMP_MISMATCH_PERCENT = 5  # Percentage of mismatches allowed (0 = use absolute count)
    MAX_SEARCH_LENGTH = 1000000  # Limit search in large chromosomes to improve performance
    
    #############################################################################
    #                           BLAST Database Options
    #############################################################################
    DB_FASTA = None            # Path to a FASTA file to create a BLAST database from
    DB_OUTPUT_DIR = None       # Custom output directory for the BLAST database
    DB_NAME = None             # Custom name for the BLAST database
    USE_CUSTOM_DB = False      # Whether to use a custom database or the default
    DB_PATH = "/Library/Application Support/Blast_DBs/Tair DB/TAIR10"
    
    # BLAST+ parameters
    BLAST_WORD_SIZE = 7
    BLAST_EVALUE = 10
    BLAST_MAX_TARGET_SEQS = 100
    BLAST_REWARD = 2
    BLAST_PENALTY = -3
    BLAST_GAPOPEN = 5
    BLAST_GAPEXTEND = 2
    BLAST_FILTER_FACTOR = 100  # E-value filtering factor
    
    # NUPACK parameters
    NUPACK_TEMPERATURE = 37  # Celsius
    NUPACK_SODIUM = 0.05     # Molar
    NUPACK_MAGNESIUM = 0.0   # Molar

    # SNP masking settings
    SNP_MASK_PADDING = 3  # Number of bases to mask on each side of a SNP
    
    #############################################################################
    #                           Primer3 Settings
    #############################################################################
    # Complete Primer3 settings dictionary
    PRIMER3_SETTINGS = {
        # General settings
        "P3_FILE_TYPE": "settings", 
        "P3_FILE_ID": "User settings", 
        "P3P_DEBUG_MODE": 0, 
        "P3P_GB_ORIENTATION": "+",
        "P3P_PRIMER_NAME_ACRONYM_INTERNAL": "IN", 
        "P3P_PRIMER_NAME_ACRONYM_LEFT": "F", 
        "P3P_PRIMER_NAME_ACRONYM_RIGHT": "R",
        "P3P_PRIMER_NAME_ACRONYM_SPACER": "_", 
        
        # Primer conditions
        "PRIMER_ANNEALING_TEMP": 52.0, 
        "PRIMER_DMSO_CONC": 0.0,
        "PRIMER_DMSO_FACTOR": 0.6, 
        "PRIMER_DNA_CONC": 50.0, 
        "PRIMER_DNTP_CONC": 0.8, 
        "PRIMER_FIRST_BASE_INDEX": 1,
        "PRIMER_FORMAMIDE_CONC": 0.0, 
        "PRIMER_GC_CLAMP": 1, 
        "PRIMER_INSIDE_PENALTY": -1.0, 
        
        # Internal oligo parameters
        "PRIMER_INTERNAL_DMSO_CONC": 0.0,
        "PRIMER_INTERNAL_DMSO_FACTOR": 0.6, 
        "PRIMER_INTERNAL_DNA_CONC": 50.0, 
        "PRIMER_INTERNAL_DNTP_CONC": 0.0,
        "PRIMER_INTERNAL_FORMAMIDE_CONC": 0.0, 
        "PRIMER_INTERNAL_MAX_BOUND": 110.0, 
        "PRIMER_INTERNAL_MAX_GC": 80.0,
        "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 47.0, 
        "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB": 12.0, 
        "PRIMER_INTERNAL_MAX_NS_ACCEPTED": 0,
        "PRIMER_INTERNAL_MAX_POLY_X": 4, 
        "PRIMER_INTERNAL_MAX_SELF_ANY": 12.0, 
        "PRIMER_INTERNAL_MAX_SELF_ANY_TH": 47.0,
        "PRIMER_INTERNAL_MAX_SELF_END": 12.0, 
        "PRIMER_INTERNAL_MAX_SELF_END_TH": 47.0, 
        "PRIMER_INTERNAL_MAX_SIZE": 27,
        "PRIMER_INTERNAL_MAX_TM": 70, 
        "PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4,
        "PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7, 
        "PRIMER_INTERNAL_MIN_BOUND": -10.0, 
        "PRIMER_INTERNAL_MIN_GC": 30.0,
        "PRIMER_INTERNAL_MIN_QUALITY": 0, 
        "PRIMER_INTERNAL_MIN_SIZE": 15, 
        "PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE": -1,
        "PRIMER_INTERNAL_MIN_TM": 64, 
        "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME": "hnnnn",
        "PRIMER_INTERNAL_OPT_BOUND": 97.0,
        "PRIMER_INTERNAL_OPT_GC_PERCENT": 50.0, 
        "PRIMER_INTERNAL_OPT_SIZE": 20, 
        "PRIMER_INTERNAL_OPT_TM": 65,
        "PRIMER_INTERNAL_SALT_DIVALENT": 0.0, 
        "PRIMER_INTERNAL_SALT_MONOVALENT": 50.0, 
        
        # General primer constraints
        "PRIMER_LIBERAL_BASE": 1,
        "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS": 0, 
        "PRIMER_LOWERCASE_MASKING": 0, 
        "PRIMER_MAX_BOUND": 110.0,
        "PRIMER_MAX_END_GC": 3, 
        "PRIMER_MAX_END_STABILITY": 9.0, 
        "PRIMER_MAX_GC": 60.0, 
        "PRIMER_MAX_HAIRPIN_TH": 47.0,
        "PRIMER_MAX_LIBRARY_MISPRIMING": 12.0, 
        "PRIMER_MAX_NS_ACCEPTED": 0, 
        "PRIMER_MAX_POLY_X": 4, 
        "PRIMER_MAX_SELF_ANY": 8.0,
        "PRIMER_MAX_SELF_ANY_TH": 47.0, 
        "PRIMER_MAX_SELF_END": 3.0, 
        "PRIMER_MAX_SELF_END_TH": 47.0, 
        "PRIMER_MAX_SIZE": 23,
        "PRIMER_MAX_TEMPLATE_MISPRIMING": 12.0, 
        "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": 47.0, 
        "PRIMER_MAX_TM": 65.0,
        "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4, 
        "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7, 
        "PRIMER_MIN_BOUND": -10.0,
        "PRIMER_MIN_END_QUALITY": 0, 
        "PRIMER_MIN_GC": 50.0, 
        "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": 3,
        "PRIMER_MIN_QUALITY": 0, 
        "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": 3, 
        "PRIMER_MIN_SIZE": 18, 
        "PRIMER_MIN_TM": 50.0,
        "PRIMER_NUM_RETURN": 10, 
        "PRIMER_OPT_BOUND": 97.0, 
        "PRIMER_OPT_GC_PERCENT": 52.5, 
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_OPT_TM": 57.5, 
        "PRIMER_OUTSIDE_PENALTY": 0.0, 
        
        # Primer pair parameters
        "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0, 
        "PRIMER_PAIR_MAX_COMPL_END": 3.0, 
        "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,
        "PRIMER_PAIR_MAX_DIFF_TM": 1, 
        "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": 24.0, 
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": 24.0,
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": 47.0, 
        "PRIMER_PAIR_WT_COMPL_ANY": 0.0, 
        "PRIMER_PAIR_WT_COMPL_ANY_TH": 0.0,
        "PRIMER_PAIR_WT_COMPL_END": 0.0, 
        "PRIMER_PAIR_WT_COMPL_END_TH": 0.0, 
        "PRIMER_PAIR_WT_DIFF_TM": 0.0,
        "PRIMER_PAIR_WT_IO_PENALTY": 0.0, 
        "PRIMER_PAIR_WT_LIBRARY_MISPRIMING": 0.0, 
        "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.0,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.0, 
        "PRIMER_PAIR_WT_PRODUCT_TM_GT": 0.0, 
        "PRIMER_PAIR_WT_PRODUCT_TM_LT": 0.0,
        "PRIMER_PAIR_WT_PR_PENALTY": 1.0, 
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING": 0.0, 
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH": 0.0,
        
        # Primer selection and product parameters
        "PRIMER_PICK_ANYWAY": 0, 
        "PRIMER_PICK_INTERNAL_OLIGO": 1, 
        "PRIMER_PICK_LEFT_PRIMER": 1, 
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PRODUCT_MAX_TM": 1000000.0, 
        "PRIMER_PRODUCT_MIN_TM": -1000000.0, 
        "PRIMER_PRODUCT_OPT_SIZE": 0,
        "PRIMER_PRODUCT_OPT_TM": 0.0, 
        "PRIMER_PRODUCT_SIZE_RANGE": "90-200", 
        "PRIMER_QUALITY_RANGE_MAX": 100,
        "PRIMER_QUALITY_RANGE_MIN": 0, 
        
        # Chemistry and thermodynamic parameters
        "PRIMER_SALT_CORRECTIONS": 1, 
        "PRIMER_SALT_DIVALENT": 3.8, 
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT": 1, 
        "PRIMER_SEQUENCING_ACCURACY": 20, 
        "PRIMER_SEQUENCING_INTERVAL": 250,
        "PRIMER_SEQUENCING_LEAD": 50, 
        "PRIMER_SEQUENCING_SPACING": 500, 
        "PRIMER_TASK": "generic",
        "PRIMER_THERMODYNAMIC_PARAMETERS_PATH": "/opt/homebrew/Cellar/primer3/2.4.0/share/primer3/primer3_config/",
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1, 
        "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": 0, 
        "PRIMER_TM_FORMULA": 1,
        
        # Penalty weights
        "PRIMER_WT_BOUND_GT": 0.0, 
        "PRIMER_WT_BOUND_LT": 0.0, 
        "PRIMER_WT_END_QUAL": 0.0, 
        "PRIMER_WT_END_STABILITY": 0.0,
        "PRIMER_WT_GC_PERCENT_GT": 0.5, 
        "PRIMER_WT_GC_PERCENT_LT": 0.5, 
        "PRIMER_WT_HAIRPIN_TH": 0.0,
        "PRIMER_WT_LIBRARY_MISPRIMING": 0.0, 
        "PRIMER_WT_NUM_NS": 0.0, 
        "PRIMER_WT_POS_PENALTY": 0.0, 
        "PRIMER_WT_SELF_ANY": 0.0,
        "PRIMER_WT_SELF_ANY_TH": 0.0, 
        "PRIMER_WT_SELF_END": 0.0, 
        "PRIMER_WT_SELF_END_TH": 0.0, 
        "PRIMER_WT_SEQ_QUAL": 0.0,
        "PRIMER_WT_SIZE_GT": 1.0, 
        "PRIMER_WT_SIZE_LT": 1.0, 
        "PRIMER_WT_TEMPLATE_MISPRIMING": 0.0,
        "PRIMER_WT_TEMPLATE_MISPRIMING_TH": 0.0, 
        "PRIMER_WT_TM_GT": 1.0, 
        "PRIMER_WT_TM_LT": 1.0
    }
    
    def __init__(self):
        """Initialize Config instance with default values."""
        # Implementation left empty as we're using class variables
        pass
    
    @classmethod
    def get_instance(cls) -> 'Config':
        """
        Get the singleton instance of Config.
        
        Returns:
            Config: Singleton instance
        """
        if cls._instance is None:
            cls._instance = cls()
            # Update MIN_SEGMENT_LENGTH based on PRIMER_PRODUCT_SIZE_RANGE
            if cls.PRIMER_PRODUCT_SIZE_RANGE and cls.PRIMER_PRODUCT_SIZE_RANGE[0]:
                cls.MIN_SEGMENT_LENGTH = cls.PRIMER_PRODUCT_SIZE_RANGE[0][0]
        return cls._instance
    
    @classmethod
    def get_primer3_global_args(cls) -> Dict[str, Any]:
        """
        Get global primer3 arguments as a dictionary.
        This combines the simplified settings with the complete settings.
        
        Returns:
            dict: Dictionary of primer3 settings
        """
        # Start with complete settings
        settings = cls.PRIMER3_SETTINGS.copy()
        
        # Override with simplified settings for commonly adjusted parameters
        settings.update({
            "PRIMER_MIN_SIZE": cls.PRIMER_MIN_SIZE,
            "PRIMER_OPT_SIZE": cls.PRIMER_OPT_SIZE,
            "PRIMER_MAX_SIZE": cls.PRIMER_MAX_SIZE,
            "PRIMER_MIN_TM": cls.PRIMER_MIN_TM,
            "PRIMER_OPT_TM": cls.PRIMER_OPT_TM,
            "PRIMER_MAX_TM": cls.PRIMER_MAX_TM,
            "PRIMER_MIN_GC": cls.PRIMER_MIN_GC,
            "PRIMER_MAX_GC": cls.PRIMER_MAX_GC,
            "PRIMER_NUM_RETURN": cls.MAX_PRIMER_PAIRS_PER_SEGMENT,
        })
        
        # Convert the product size range from list to string format if needed
        if isinstance(cls.PRIMER_PRODUCT_SIZE_RANGE, list):
            # Convert [[min1, max1], [min2, max2], ...] to "min1-max1 min2-max2 ..."
            size_range_str = " ".join([f"{r[0]}-{r[1]}" for r in cls.PRIMER_PRODUCT_SIZE_RANGE])
            settings["PRIMER_PRODUCT_SIZE_RANGE"] = size_range_str
        
        return settings
    
    @classmethod
    def format_settings_for_file(cls) -> str:
        """
        Format settings for writing to a primer3 settings file.
        
        Returns:
            str: Formatted settings string
        """
        settings = cls.get_primer3_global_args()
        return "\n".join(f"{key}={value}" for key, value in settings.items()) + "\n"
    
    @classmethod
    def load_from_file(cls, filepath: str) -> bool:
        """
        Load settings from a configuration file.
        Supports both Primer3 format and JSON format.
        
        Args:
            filepath: Path to the settings file
            
        Returns:
            bool: True if settings were loaded successfully
        """
        try:
            # Initialize the singleton if not already done
            cls.get_instance()
            
            # JSON configuration
            if filepath.endswith('.json'):
                return cls._load_from_json(filepath)
            # Primer3 format (key=value)
            else:
                return cls._load_from_primer3_format(filepath)
                
        except Exception as e:
            print(f"Error loading settings from {filepath}: {e}")
            return False
    
    @classmethod
    def _load_from_json(cls, filepath: str) -> bool:
        """
        Load settings from a JSON file.
        
        Args:
            filepath: Path to the JSON file
            
        Returns:
            bool: True if settings were loaded successfully
        """
        try:
            with open(filepath, 'r') as f:
                settings = json.load(f)
            
            # Update class attributes based on JSON
            for key, value in settings.items():
                if hasattr(cls, key):
                    setattr(cls, key, value)
                elif key in cls.PRIMER3_SETTINGS:
                    cls.PRIMER3_SETTINGS[key] = value
            
            # Handle special case for PRIMER_PRODUCT_SIZE_RANGE
            if "PRIMER_PRODUCT_SIZE_RANGE" in settings:
                value = settings["PRIMER_PRODUCT_SIZE_RANGE"]
                if isinstance(value, list):
                    cls.PRIMER_PRODUCT_SIZE_RANGE = value
                    # Also update MIN_SEGMENT_LENGTH
                    if value and value[0]:
                        cls.MIN_SEGMENT_LENGTH = value[0][0]
            
            return True
        except Exception as e:
            print(f"Error loading JSON settings from {filepath}: {e}")
            return False
    
    @classmethod
    def _load_from_primer3_format(cls, filepath: str) -> bool:
        """
        Load settings from a Primer3 settings file (key=value format).
        
        Args:
            filepath: Path to the Primer3 settings file
            
        Returns:
            bool: True if settings were loaded successfully
        """
        try:
            with open(filepath, 'r') as f:
                settings_text = f.read()
                
            # Parse settings
            settings = {}
            for line in settings_text.strip().split('\n'):
                if '=' in line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Try to convert to appropriate type
                    if value.replace('.', '', 1).isdigit():
                        if '.' in value:
                            value = float(value)
                        else:
                            value = int(value)
                    elif value.lower() in ('true', 'false'):
                        value = value.lower() == 'true'
                        
                    settings[key] = value
            
            # Update settings
            cls.PRIMER3_SETTINGS.update(settings)
            
            # Update simplified settings if corresponding keys exist
            simplified_settings_map = {
                "PRIMER_MIN_SIZE": "PRIMER_MIN_SIZE",
                "PRIMER_OPT_SIZE": "PRIMER_OPT_SIZE",
                "PRIMER_MAX_SIZE": "PRIMER_MAX_SIZE",
                "PRIMER_MIN_TM": "PRIMER_MIN_TM",
                "PRIMER_OPT_TM": "PRIMER_OPT_TM",
                "PRIMER_MAX_TM": "PRIMER_MAX_TM",
                "PRIMER_MIN_GC": "PRIMER_MIN_GC",
                "PRIMER_MAX_GC": "PRIMER_MAX_GC",
                "PRIMER_NUM_RETURN": "MAX_PRIMER_PAIRS_PER_SEGMENT",
            }
            
            for primer3_key, config_attr in simplified_settings_map.items():
                if primer3_key in settings:
                    setattr(cls, config_attr, settings[primer3_key])
            
            # Handle product size range
            if "PRIMER_PRODUCT_SIZE_RANGE" in settings:
                value = settings["PRIMER_PRODUCT_SIZE_RANGE"]
                if isinstance(value, str):
                    # Parse "min1-max1 min2-max2 ..." format
                    ranges = []
                    for r in value.split():
                        if '-' in r:
                            min_val, max_val = r.split('-')
                            ranges.append([int(min_val), int(max_val)])
                    cls.PRIMER_PRODUCT_SIZE_RANGE = ranges
                    
                    # Update MIN_SEGMENT_LENGTH
                    if ranges and ranges[0]:
                        cls.MIN_SEGMENT_LENGTH = ranges[0][0]
            
            # Handle BLAST database options
            blast_settings_map = {
                "DB_FASTA": "DB_FASTA",
                "DB_OUTPUT_DIR": "DB_OUTPUT_DIR",
                "DB_NAME": "DB_NAME",
                "DB_PATH": "DB_PATH"
            }
            
            for setting_key, config_attr in blast_settings_map.items():
                if setting_key in settings:
                    setattr(cls, config_attr, settings[setting_key])
                    if setting_key == "DB_FASTA":
                        cls.USE_CUSTOM_DB = True
            
            return True
        except Exception as e:
            print(f"Error loading Primer3 settings from {filepath}: {e}")
            return False
    
    @classmethod
    def save_to_file(cls, filepath: str, format_type: str = "primer3") -> bool:
        """
        Save current settings to a file.
        
        Args:
            filepath: Path to save the settings
            format_type: Format to use ("primer3" or "json")
            
        Returns:
            bool: True if settings were saved successfully
        """
        try:
            if format_type.lower() == "json":
                return cls._save_to_json(filepath)
            else:  # Default to primer3 format
                return cls._save_to_primer3_format(filepath)
        except Exception as e:
            print(f"Error saving settings to {filepath}: {e}")
            return False
    
    @classmethod
    def _save_to_primer3_format(cls, filepath: str) -> bool:
        """
        Save settings in Primer3 format (key=value).
        
        Args:
            filepath: Path to save the settings
            
        Returns:
            bool: True if settings were saved successfully
        """
        try:
            with open(filepath, 'w') as f:
                f.write(cls.format_settings_for_file())
            return True
        except Exception as e:
            print(f"Error saving Primer3 settings to {filepath}: {e}")
            return False
    
    @classmethod
    def _save_to_json(cls, filepath: str) -> bool:
        """
        Save settings in JSON format.
        
        Args:
            filepath: Path to save the settings
            
        Returns:
            bool: True if settings were saved successfully
        """
        try:
            # Get all public class attributes (excluding methods, private attrs, etc.)
            settings = {}
            
            # Add all class variables that don't start with underscore 
            # and aren't methods or callables
            for key in dir(cls):
                if not key.startswith('_') and not callable(getattr(cls, key)):
                    value = getattr(cls, key)
                    # Only include serializable types
                    if isinstance(value, (str, int, float, bool, list, dict, tuple)) or value is None:
                        settings[key] = value
            
            with open(filepath, 'w') as f:
                json.dump(settings, f, indent=4)
            
            return True
        except Exception as e:
            print(f"Error saving JSON settings to {filepath}: {e}")
            return False
    
    @classmethod
    def get_all_settings(cls) -> Dict[str, Any]:
        """
        Get all settings as a dictionary.
        
        Returns:
            dict: Dictionary of all settings
        """
        settings = {}
        
        # Add all class variables that don't start with underscore
        for key in dir(cls):
            if not key.startswith('_') and not callable(getattr(cls, key)):
                settings[key] = getattr(cls, key)
        
        return settings
    
    @staticmethod
    def debug(message: str) -> None:
        """
        Print debug messages if debug mode is enabled.
        
        Args:
            message: The debug message to print
        """
        if Config.DEBUG_MODE:
            print(f"[DEBUG] {message}")