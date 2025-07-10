#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration display module for ddPrimer pipeline.

Contains functionality for:
1. Structured display of configuration settings
2. Color-coded output for improved readability
3. Categorized setting organization
4. Primer3 settings visualization

This module provides comprehensive configuration display capabilities
for the ddPrimer pipeline, allowing users to view current settings
in an organized and readable format.
"""

import textwrap
import colorama
from colorama import Fore, Style
import logging

logger = logging.getLogger(__name__)


def display_config(config_cls):
    """
    Display all configuration settings in a structured, easy-to-read format.
    
    Args:
        config_cls: The Config class containing all settings
        
    Example:
        >>> from ddprimer.config import Config
        >>> display_config(Config)
    """
    logger.debug("Starting configuration display")
    
    try:
        # Initialize colorama for cross-platform colored output
        colorama.init()
        
        settings = config_cls.get_all_settings()
        logger.debug(f"Retrieved {len(settings)} configuration settings")
        
        # Print header
        print(f"\n{Fore.WHITE}{'='*80}{Style.RESET_ALL}")
        print(f"{Fore.WHITE}{'ddPrimer Configuration Settings':^80}{Style.RESET_ALL}")
        print(f"{Fore.WHITE}{'='*80}{Style.RESET_ALL}\n")
        
        # Group settings by category based on class structure
        categories = {
            "Performance Settings": [
                "NUM_PROCESSES", "SHOW_PROGRESS"
            ],
            "Basic Primer3 Parameters": [
                "PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE",
                "PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM",
                "PRIMER_MIN_GC", "PRIMER_MAX_GC", "PRIMER_PRODUCT_SIZE_RANGE",
            ],
            "Pipeline parameters": [                
                "MIN_SEGMENT_LENGTH", "RETAIN_TYPES", "GENE_OVERLAP_MARGIN", "RESTRICTION_SITE", 
                "PENALTY_MAX", "MAX_PRIMER_PAIRS_PER_SEGMENT", "PREFER_PROBE_MORE_C_THAN_G", 
                "SEQUENCE_MIN_GC", "SEQUENCE_MAX_GC", "MIN_CHROMOSOME_SIZE"
            ],
            "VCF/SNP Processing parameters": [    
                "VCF_ALLELE_FREQUENCY_THRESHOLD", "VCF_QUALITY_THRESHOLD", 
                "VCF_FLANKING_MASK_SIZE", "VCF_USE_SOFT_MASKING"
            ],
            "ViennaRNA Parameters": [
                "THERMO_TEMPERATURE", "THERMO_SODIUM", "THERMO_MAGNESIUM"
            ],
            "BLAST Database Options": [
                "BLAST_WORD_SIZE", "BLAST_EVALUE", 
                "BLAST_MAX_TARGET_SEQS", "BLAST_REWARD", "BLAST_PENALTY", 
                "BLAST_GAPOPEN", "BLAST_GAPEXTEND", "BLAST_FILTER_FACTOR"
            ]
        }
        
        # Create "Other" category for any settings not explicitly categorized
        categorized_keys = []
        for keys in categories.values():
            categorized_keys.extend(keys)
        
        other_keys = [key for key in settings.keys() if key not in categorized_keys 
                     and not key.startswith('_') and key != "PRIMER3_SETTINGS"]
        
        if other_keys:
            categories["Other"] = other_keys
        
        # Print settings by category
        for category, keys in categories.items():
            print(f"{Fore.WHITE}{category}{Style.RESET_ALL}")
            print(f"{Fore.WHITE}{'-' * len(category)}{Style.RESET_ALL}")
            
            for key in keys:
                if key in settings:
                    value = settings[key]
                    # Format value for display
                    if isinstance(value, list) and len(str(value)) > 60:
                        formatted_value = "\n" + textwrap.indent(str(value), " " * 4)
                    elif isinstance(value, str) and value is None:
                        formatted_value = "None"
                    else:
                        formatted_value = str(value)
                    
                    print(f"{Fore.CYAN}{key}{Style.RESET_ALL}: {formatted_value}")
            print()
        
        # Handle Primer3 settings separately - these are extensive
        print(f"{Fore.WHITE}Extensive Primer3 Settings{Style.RESET_ALL}")
        print(f"{Fore.WHITE}{'-' * 14}{Style.RESET_ALL}")
        print(f"{Fore.CYAN}Use \"--config primer3\" to display all Primer3 settings{Style.RESET_ALL}\n")
        
        # Print footer with usage instructions
        print(f"{Fore.WHITE}{'-'*80}{Style.RESET_ALL}")
        print(f"\n{Fore.WHITE}To modify the settings:{Style.RESET_ALL}")
        print(f"{Fore.WHITE}Use custom config: {Fore.CYAN}ddprimer --config your_config.json{Style.RESET_ALL}")
        print(f"{Fore.WHITE}Generate a template config file: {Fore.CYAN}ddprimer --config template{Style.RESET_ALL}")
        print(f"\n{Fore.WHITE}{'='*80}{Style.RESET_ALL}")
        print(f"{Style.RESET_ALL}\n")
        
        logger.debug("Configuration display completed successfully")
        
    except Exception as e:
        error_msg = f"Configuration display failed"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise ConfigDisplayError(error_msg) from e


def display_primer3_settings(config_cls):
    """
    Display all Primer3 settings in a structured format.
    
    Args:
        config_cls: The Config class containing Primer3 settings
        
    Raises:
        ConfigDisplayError: If display operation fails
        
    Example:
        >>> from ddprimer.config import Config
        >>> display_primer3_settings(Config)
    """
    logger.debug("Starting Primer3 settings display")
    
    try:
        # Initialize colorama for cross-platform colored output
        colorama.init()
        
        primer3_settings = config_cls.PRIMER3_SETTINGS
        logger.debug(f"Retrieved {len(primer3_settings)} Primer3 settings")
        
        # Print header
        print(f"\n{Fore.WHITE}{'='*80}{Style.RESET_ALL}")
        print(f"{Fore.WHITE}{'ddPrimer Primer3 Configuration Settings':^80}{Style.RESET_ALL}")
        print(f"{Fore.WHITE}{'='*80}{Style.RESET_ALL}\n")
        
        # Group Primer3 settings by categories
        categories = {
            "Basic Primer3 Parameters": [
                "PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE",
                "PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM",
                "PRIMER_MIN_GC", "PRIMER_MAX_GC", "PRIMER_PRODUCT_SIZE_RANGE",
            ],
            "General Settings": [key for key in primer3_settings.keys() if key.startswith("P3_")],
            "Primer Name Settings": [key for key in primer3_settings.keys() if key.startswith("P3P_PRIMER_NAME")],
            "Primer Conditions": [key for key in primer3_settings.keys() if key.startswith("PRIMER_") and 
                                 not key.startswith("PRIMER_INTERNAL_") and 
                                 not key.startswith("PRIMER_PAIR_") and
                                 not key.startswith("PRIMER_PICK_") and
                                 not key.startswith("PRIMER_PRODUCT_") and
                                 not key.startswith("PRIMER_WT_")],
            "Internal Oligo Parameters": [key for key in primer3_settings.keys() if key.startswith("PRIMER_INTERNAL_")],
            "Primer Pair Parameters": [key for key in primer3_settings.keys() if key.startswith("PRIMER_PAIR_")],
            "Selection and Product Parameters": [
                key for key in primer3_settings.keys() 
                if key.startswith("PRIMER_PICK_") or key.startswith("PRIMER_PRODUCT_")
            ],
            "Chemistry and Thermodynamic Parameters": [
                key for key in primer3_settings.keys() 
                if key.startswith("PRIMER_SALT_") or 
                key.startswith("PRIMER_SEQUENCING_") or
                key.startswith("PRIMER_TM_") or
                key.startswith("PRIMER_THERMODYNAMIC_")
            ],
            "Penalty Weights": [key for key in primer3_settings.keys() if key.startswith("PRIMER_WT_")]
        }
        
        # Print settings by category
        for category, keys in categories.items():
            if keys:  # Only print categories that have keys
                print(f"{Fore.WHITE}{category}{Style.RESET_ALL}")
                print(f"{Fore.WHITE}{'-' * len(category)}{Style.RESET_ALL}")
                
                for key in sorted(keys):
                    value = primer3_settings[key]
                    print(f"{Fore.CYAN}{key}{Style.RESET_ALL}: {value}")
                print()
        
        # Collect any uncategorized settings
        all_categorized = []
        for keys in categories.values():
            all_categorized.extend(keys)
        
        uncategorized = [key for key in primer3_settings.keys() if key not in all_categorized]
        
        if uncategorized:
            print(f"{Fore.WHITE}Other Primer3 Settings{Style.RESET_ALL}")
            print(f"{Fore.WHITE}{'-' * 20}{Style.RESET_ALL}")
            for key in sorted(uncategorized):
                value = primer3_settings[key]
                print(f"{Fore.CYAN}{key}{Style.RESET_ALL}: {value}")
            print()
        
        # Print footer with usage instructions
        print(f"{Fore.WHITE}{'-'*80}{Style.RESET_ALL}")
        print(f"\n{Fore.WHITE}To modify the settings:{Style.RESET_ALL}")
        print(f"{Fore.WHITE}Use custom config: {Fore.CYAN}ddprimer --config your_config.json{Style.RESET_ALL}")
        print(f"{Fore.WHITE}Generate a template config file: {Fore.CYAN}ddprimer --config template{Style.RESET_ALL}")
        print(f"\n{Fore.WHITE}{'='*80}{Style.RESET_ALL}")
        print(f"{Style.RESET_ALL}\n")
        
        logger.debug("Primer3 settings display completed successfully")
        
    except Exception as e:
        error_msg = f"Primer3 settings display failed"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise ConfigDisplayError(error_msg) from e


class ConfigDisplayError(Exception):
    """Error during configuration display operations."""
    pass