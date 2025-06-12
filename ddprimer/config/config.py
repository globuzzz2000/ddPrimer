#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration module for ddPrimer pipeline.

Contains functionality for:
1. Central configuration settings management with singleton pattern
2. JSON and Primer3 format configuration file loading/saving
3. BLAST database path management and auto-detection
4. Comprehensive Primer3 settings configuration

This module provides centralized configuration management for the ddPrimer
pipeline, supporting both simple parameter overrides and complete
Primer3 settings customization.
"""

import os
import json
from pathlib import Path
import logging
from multiprocessing import cpu_count
from typing import Dict, Any

logger = logging.getLogger(__name__)


class Config:
    """
    Central configuration settings for ddPrimer pipeline with singleton pattern.
    
    This class manages all configuration settings for the ddPrimer pipeline,
    providing a singleton pattern for consistent settings access across
    all modules and comprehensive configuration file support.
    
    Attributes:
        DEBUG_MODE: Enable debug logging mode
        NUM_PROCESSES: Number of parallel processes to use
        PRIMER_MIN_SIZE: Minimum primer length
        PRIMER_MAX_SIZE: Maximum primer length
        
    Example:
        >>> config = Config.get_instance()
        >>> config.PRIMER_MIN_SIZE = 20
        >>> Config.load_from_file("my_config.json")
    """
    
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
    MIN_SEGMENT_LENGTH = 90
    RETAIN_TYPES = ["gene", "mRNA"]  # gff filtering: "gene", "mRNA", "CDS", "exon", etc.
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
    #                           SNP Masking Parameters
    #############################################################################
    SNP_ALLELE_FREQUENCY_THRESHOLD = None    # Minimum allele frequency (AF) to mask (e.g., None, 0.05 for 5%)
    SNP_QUALITY_THRESHOLD = None             # Minimum QUAL score to include variants (e.g., None, 30.0)
    SNP_FLANKING_MASK_SIZE = 0               # Number of bases to mask around each SNP (0 = just the SNP)
    SNP_USE_SOFT_MASKING = False             # Use lowercase letters instead of 'N' characters
    
    #############################################################################
    #                           Thermodynamic Calculation Settings
    #############################################################################
    THERMO_TEMPERATURE = 37  # Celsius
    THERMO_SODIUM = 0.05     # Molar
    THERMO_MAGNESIUM = 0.0   # Molar

    #############################################################################
    #                           BLAST Database Options
    #############################################################################
    DB_FASTA = None            # Path to a FASTA file to create a BLAST database from
    DB_OUTPUT_DIR = None       # Custom output directory for the BLAST database
    DB_NAME = None             # Custom name for the BLAST database
    USE_CUSTOM_DB = False      # Whether to use a custom database or an auto-detected one
    DB_PATH = None
    
    # BLAST+ parameters
    BLAST_WORD_SIZE = 7
    BLAST_EVALUE = 10
    BLAST_MAX_TARGET_SEQS = 100
    BLAST_REWARD = 2
    BLAST_PENALTY = -3
    BLAST_GAPOPEN = 5
    BLAST_GAPEXTEND = 2
    BLAST_FILTER_FACTOR = 100  # E-value filtering factor


    #############################################################################
    #                           Alignment Parameters
    #############################################################################
    # Sequence identity settings for alignment mode
    MIN_IDENTITY = 80.0  # Minimum sequence identity for primer regions (percentage)
    MIN_LENGTH = 20      # Minimum length of conserved regions

    # LastZ alignment options
    LASTZ_OPTIONS = "--format=maf --chain --gapped --step=10 --ambiguous=iupac"
    # Options commonly used:
    # --format=maf      : Output in Multiple Alignment Format (required)
    # --step=10         : Step size for search (smaller = more sensitive but slower)
    # --hspthresh=3000  : Threshold for high-scoring segment pairs (lower = more alignments)
    # --chain           : Enable chaining of HSPs
    # --gapped          : Enable gapped extension of HSPs
    # --inner=2000      : Inner radius for seeding alignments
    # --ydrop=3400      : Y-drop parameter (controls extension termination)
    # --gappedthresh=3000 : Threshold for gapped alignments

        
    #############################################################################
    #                           Unified Configuration Storage
    #############################################################################
    
    @classmethod
    def get_user_config_dir(cls) -> str:
        """
        Get the user configuration directory, creating it if necessary.
        
        Returns:
            Path to user configuration directory
            
        Raises:
            ConfigError: If directory cannot be created
        """
        config_dir = os.path.join(os.path.expanduser("~"), ".ddprimer")
        
        try:
            os.makedirs(config_dir, exist_ok=True)
            logger.debug(f"User config directory: {config_dir}")
            return config_dir
        except OSError as e:
            error_msg = f"Failed to create user config directory {config_dir}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def save_user_setting(cls, key: str, value: Any) -> bool:
        """
        Save a user setting to persistent storage.
        
        Args:
            key: Setting key name
            value: Setting value (must be JSON-serializable)
            
        Returns:
            True if setting was saved successfully
            
        Raises:
            ConfigError: If setting cannot be saved
        """
        logger.debug(f"Saving user setting: {key} = {value}")
        
        try:
            config_dir = cls.get_user_config_dir()
            settings_file = os.path.join(config_dir, "user_settings.json")
            
            # Load existing settings
            settings = {}
            if os.path.exists(settings_file):
                try:
                    with open(settings_file, 'r') as f:
                        settings = json.load(f)
                except (json.JSONDecodeError, OSError) as e:
                    logger.warning(f"Error reading existing settings, starting fresh: {str(e)}")
                    settings = {}
            
            # Update with new setting
            settings[key] = value
            
            # Save back to file
            with open(settings_file, 'w') as f:
                json.dump(settings, f, indent=2)
            
            logger.debug(f"Successfully saved user setting: {key}")
            return True
            
        except Exception as e:
            error_msg = f"Failed to save user setting {key}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def load_user_setting(cls, key: str, default: Any = None) -> Any:
        """
        Load a user setting from persistent storage.
        
        Args:
            key: Setting key name
            default: Default value if setting doesn't exist
            
        Returns:
            Setting value or default if not found
            
        Raises:
            ConfigError: If there's a critical error reading settings
        """
        logger.debug(f"Loading user setting: {key}")
        
        try:
            config_dir = cls.get_user_config_dir()
            settings_file = os.path.join(config_dir, "user_settings.json")
            
            if not os.path.exists(settings_file):
                logger.debug(f"Settings file doesn't exist, returning default for {key}")
                return default
            
            with open(settings_file, 'r') as f:
                settings = json.load(f)
            
            value = settings.get(key, default)
            logger.debug(f"Loaded user setting {key} = {value}")
            return value
            
        except (json.JSONDecodeError, OSError) as e:
            logger.warning(f"Error reading user settings for {key}, returning default: {str(e)}")
            return default
        except Exception as e:
            error_msg = f"Critical error loading user setting {key}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e


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
            
        Example:
            >>> config = Config.get_instance()
            >>> config.PRIMER_MIN_SIZE = 20
        """
        if cls._instance is None:
            logger.debug("Creating new Config singleton instance")
            cls._instance = cls()
            # Update MIN_SEGMENT_LENGTH based on PRIMER_PRODUCT_SIZE_RANGE
            if cls.PRIMER_PRODUCT_SIZE_RANGE and cls.PRIMER_PRODUCT_SIZE_RANGE[0]:
                cls.MIN_SEGMENT_LENGTH = cls.PRIMER_PRODUCT_SIZE_RANGE[0][0]
                logger.debug(f"Updated MIN_SEGMENT_LENGTH to {cls.MIN_SEGMENT_LENGTH}")
            
            # Initialize the BLAST database path if not already set
            if cls.DB_PATH is None:
                cls.DB_PATH = cls._initialize_db_path()
                logger.debug(f"Initialized DB_PATH to {cls.DB_PATH}")
        return cls._instance
    
    @classmethod
    def get_primer3_global_args(cls) -> Dict[str, Any]:
        """
        Get global primer3 arguments as a dictionary.
        
        This combines the simplified settings with the complete settings,
        allowing for easy parameter override while maintaining full
        Primer3 compatibility.
        
        Returns:
            dict: Dictionary of primer3 settings ready for use
            
        Example:
            >>> args = Config.get_primer3_global_args()
            >>> args['PRIMER_MIN_SIZE']
            18
        """
        logger.debug("=== PRIMER3 ARGS DEBUG ===")
        logger.debug("Generating Primer3 global arguments")
        
        try:
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
                logger.debug(f"Converted product size range to: {size_range_str}")
            
            logger.debug(f"Generated {len(settings)} Primer3 arguments")
            return settings
            
        except Exception as e:
            error_msg = f"Failed to generate Primer3 global arguments"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
        
    logger.debug("=== END PRIMER3 ARGS DEBUG ===")
    
    @classmethod
    def format_settings_for_file(cls) -> str:
        """
        Format settings for writing to a primer3 settings file.
        
        Returns:
            str: Formatted settings string in key=value format
            
        Raises:
            ConfigError: If formatting fails
            
        Example:
            >>> formatted = Config.format_settings_for_file()
            >>> print(formatted.split('\\n')[0])
            P3_FILE_TYPE=settings
        """
        try:
            settings = cls.get_primer3_global_args()
            formatted = "\n".join(f"{key}={value}" for key, value in settings.items()) + "\n"
            logger.debug(f"Formatted {len(settings)} settings for file output")
            return formatted
        except Exception as e:
            error_msg = f"Failed to format settings for file"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def load_from_file(cls, filepath: str) -> bool:
        """
        Load settings from a configuration file.
        
        Supports both Primer3 format (key=value) and JSON format,
        automatically detecting the format based on file extension.
        
        Args:
            filepath: Path to the settings file
            
        Returns:
            bool: True if settings were loaded successfully
            
        Raises:
            ConfigError: If file loading fails
            FileFormatError: If file format is invalid
            
        Example:
            >>> success = Config.load_from_file("my_config.json")
            >>> if success:
            ...     print("Settings loaded successfully")
        """
        logger.debug("=== CONFIG LOAD DEBUG ===")
        logger.debug(f"Loading configuration from {filepath}")
        
        try:
            # Initialize the singleton if not already done
            cls.get_instance()
            
            if not os.path.exists(filepath):
                error_msg = f"Configuration file not found: {filepath}"
                logger.error(error_msg)
                raise FileFormatError(error_msg)
            
            # JSON configuration
            if filepath.endswith('.json'):
                result = cls._load_from_json(filepath)
            # Primer3 format (key=value)
            else:
                result = cls._load_from_primer3_format(filepath)
            
            if result:
                logger.debug("Configuration loaded successfully")
            else:
                logger.warning("Configuration loading returned False")
                
            return result
                
        except Exception as e:
            error_msg = f"Failed to load settings from {filepath}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
        
    logger.debug("=== END CONFIG LOAD DEBUG ===")
    
    @classmethod
    def _load_from_json(cls, filepath: str) -> bool:
        """
        Load settings from a JSON file.
        
        Args:
            filepath: Path to the JSON file
            
        Returns:
            bool: True if settings were loaded successfully
            
        Raises:
            ConfigError: If JSON loading fails
        """
        try:
            logger.debug(f"Loading JSON configuration from {filepath}")
            
            with open(filepath, 'r') as f:
                settings = json.load(f)
            
            logger.debug(f"Loaded {len(settings)} settings from JSON")
            
            # Update class attributes based on JSON
            for key, value in settings.items():
                if hasattr(cls, key):
                    setattr(cls, key, value)
                    logger.debug(f"Updated {key} = {value}")
                elif key in cls.PRIMER3_SETTINGS:
                    cls.PRIMER3_SETTINGS[key] = value
                    logger.debug(f"Updated Primer3 setting {key} = {value}")
            
            # Handle special case for PRIMER_PRODUCT_SIZE_RANGE
            if "PRIMER_PRODUCT_SIZE_RANGE" in settings:
                value = settings["PRIMER_PRODUCT_SIZE_RANGE"]
                if isinstance(value, list):
                    cls.PRIMER_PRODUCT_SIZE_RANGE = value
                    # Also update MIN_SEGMENT_LENGTH
                    if value and value[0]:
                        cls.MIN_SEGMENT_LENGTH = value[0][0]
                        logger.debug(f"Updated MIN_SEGMENT_LENGTH to {cls.MIN_SEGMENT_LENGTH}")
            
            return True
        except Exception as e:
            error_msg = f"Failed to load JSON settings from {filepath}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def _load_from_primer3_format(cls, filepath: str) -> bool:
        """
        Load settings from a Primer3 settings file (key=value format).
        
        Args:
            filepath: Path to the Primer3 settings file
            
        Returns:
            bool: True if settings were loaded successfully
            
        Raises:
            ConfigError: If Primer3 format loading fails
        """
        try:
            logger.debug(f"Loading Primer3 format configuration from {filepath}")
            
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
            
            logger.debug(f"Parsed {len(settings)} settings from Primer3 format")
            
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
                    logger.debug(f"Updated {config_attr} = {settings[primer3_key]}")
            
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
                        logger.debug(f"Updated MIN_SEGMENT_LENGTH to {cls.MIN_SEGMENT_LENGTH}")
            
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
                        logger.debug("Enabled USE_CUSTOM_DB due to DB_FASTA setting")
            
            return True
        except Exception as e:
            error_msg = f"Failed to load Primer3 settings from {filepath}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def save_to_file(cls, filepath: str, format_type: str = "primer3") -> bool:
        """
        Save current settings to a file.
        
        Args:
            filepath: Path to save the settings
            format_type: Format to use ("primer3" or "json")
            
        Returns:
            bool: True if settings were saved successfully
            
        Raises:
            ConfigError: If saving fails
            
        Example:
            >>> Config.save_to_file("my_config.json", "json")
            True
        """
        logger.debug(f"Saving configuration to {filepath} in {format_type} format")
        
        try:
            if format_type.lower() == "json":
                result = cls._save_to_json(filepath)
            else:  # Default to primer3 format
                result = cls._save_to_primer3_format(filepath)
            
            if result:
                logger.debug("Configuration saved successfully")
            else:
                logger.warning("Configuration saving returned False")
                
            return result
        except Exception as e:
            error_msg = f"Failed to save settings to {filepath}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def _save_to_primer3_format(cls, filepath: str) -> bool:
        """
        Save settings in Primer3 format (key=value).
        
        Args:
            filepath: Path to save the settings
            
        Returns:
            bool: True if settings were saved successfully
            
        Raises:
            ConfigError: If Primer3 format saving fails
        """
        try:
            with open(filepath, 'w') as f:
                f.write(cls.format_settings_for_file())
            logger.debug(f"Saved Primer3 format settings to {filepath}")
            return True
        except Exception as e:
            error_msg = f"Failed to save Primer3 settings to {filepath}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e

    @classmethod
    def _initialize_db_path(cls):
        """
        Find a BLAST database path for use using unified config system.
        
        Returns:
            str or None: Path to BLAST database or None if not found
            
        Raises:
            ConfigError: If database initialization fails
        """
        logger.debug("=== DB PATH INITIALIZATION DEBUG ===")
        logger.debug("Searching for BLAST database path")
        
        try:
            # First check for a saved database path in unified config
            saved_db_path = cls.load_user_setting("blast_db_path", None)
            
            if saved_db_path:
                logger.debug(f"Found saved database path: {saved_db_path}")
                # Verify this path still has database files
                if os.path.exists(saved_db_path + ".nhr"):
                    logger.debug(f"Using configured BLAST database: {saved_db_path}")
                    return saved_db_path
                else:
                    logger.debug(f"Saved database path {saved_db_path} is not valid (missing .nhr file)")
            else:
                logger.debug("No saved database path found in unified config")
            
            # If no saved config or the saved path is invalid, 
            # fall back to looking in standard directories
            possible_locations = [
                # System directory
                "/usr/local/share/ddprimer/blast_db",
                # User directory
                os.path.join(Path.home(), ".ddprimer", "blast_db")
            ]
            
            for location in possible_locations:
                if os.path.exists(location):
                    logger.debug(f"Checking for BLAST databases in {location}")
                    # Look for any .nhr files which indicate BLAST databases
                    db_files = [f for f in os.listdir(location) if f.endswith(".nhr")]
                    if db_files:
                        # Use the first database found, without the extension
                        db_name = os.path.splitext(db_files[0])[0]
                        db_path = os.path.join(location, db_name)
                        logger.debug(f"Found existing BLAST database: {db_path}")
                        return db_path
                    else:
                        logger.debug(f"No database files (.nhr) found in {location}")
            
            # If no database found, return None
            logger.debug("No existing BLAST database found")
            return None
            
        except Exception as e:
            error_msg = f"Failed to initialize database path"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
        
        logger.debug("=== END DB PATH INITIALIZATION DEBUG ===")
    
    @classmethod
    def save_database_config(cls, db_path):
        """
        Save the database path to unified config system for future runs.
        
        Args:
            db_path (str): Path to the BLAST database
            
        Raises:
            ConfigError: If saving database config fails
            
        Example:
            >>> Config.save_database_config("/path/to/blast_db")
        """
        logger.debug(f"Saving database config: {db_path}")
        
        try:
            cls.save_user_setting("blast_db_path", db_path)
            logger.debug(f"Database config saved to unified storage: {db_path}")
            
        except Exception as e:
            error_msg = f"Failed to save database config for path {db_path}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def _save_to_json(cls, filepath: str) -> bool:
        """
        Save settings in JSON format.
        
        Args:
            filepath: Path to save the settings
            
        Returns:
            bool: True if settings were saved successfully
            
        Raises:
            ConfigError: If JSON format saving fails
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
            
            logger.debug(f"Saved JSON format settings to {filepath}")
            return True
        except Exception as e:
            error_msg = f"Failed to save JSON settings to {filepath}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def get_all_settings(cls) -> Dict[str, Any]:
        """
        Get all settings as a dictionary.
        
        Returns:
            dict: Dictionary of all configuration settings
            
        Example:
            >>> settings = Config.get_all_settings()
            >>> print(settings['PRIMER_MIN_SIZE'])
            18
        """
        settings = {}
        
        # Add all class variables that don't start with underscore
        for key in dir(cls):
            if not key.startswith('_') and not callable(getattr(cls, key)):
                settings[key] = getattr(cls, key)
        
        logger.debug(f"Retrieved {len(settings)} configuration settings")
        return settings
    
    @staticmethod
    def debug(message: str) -> None:
        """
        Print debug messages if debug mode is enabled.
        
        Args:
            message: The debug message to print
            
        Example:
            >>> Config.debug("This is a debug message")
        """
        if Config.DEBUG_MODE:
            print(f"[DEBUG] {message}")


class ConfigError(Exception):
    """Error with configuration parameters or operations."""
    pass


class FileFormatError(Exception):
    """Error with file formatting or parsing."""
    pass