#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration display module for the ddPrimer pipeline.
"""

import textwrap
import colorama
from colorama import Fore, Style

def display_config(config_cls):
    """
    Display all configuration settings in a structured, easy-to-read format.
    
    Args:
        config_cls: The Config class
    """
    # Initialize colorama for cross-platform colored output
    colorama.init()
    
    settings = config_cls.get_all_settings()
    
    # Print header
    print(f"\n{Fore.CYAN}{'='*80}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{'ddPrimer Configuration Settings':^80}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{'='*80}{Style.RESET_ALL}\n")
    
    # Group settings by category based on class structure
    categories = {
        "Pipeline Mode Options": [
            "DEBUG_MODE", "DISABLE_INTERNAL_OLIGO"
        ],
        "Performance Settings": [
            "NUM_PROCESSES", "BATCH_SIZE", "MAF_CHUNK_SIZE", "SHOW_PROGRESS"
        ],
        "Design Parameters": [
            "PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE",
            "PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM",
            "PRIMER_MIN_GC", "PRIMER_MAX_GC", "PRIMER_PRODUCT_SIZE_RANGE",
            "MIN_SEGMENT_LENGTH", "RETAIN_TYPES", "FILTER_MEANINGFUL_NAMES",
            "COUNT_AMBIGUOUS_AS_MISMATCH", "GENE_OVERLAP_MARGIN",
            "RESTRICTION_SITE", "PENALTY_MAX", "MAX_PRIMER_PAIRS_PER_SEGMENT",
            "PREFER_PROBE_MORE_C_THAN_G", "SEQUENCE_MIN_GC", "SEQUENCE_MAX_GC"
        ],
        "BLAST Database Options": [
            "DB_FASTA", "DB_OUTPUT_DIR", "DB_NAME", "USE_CUSTOM_DB", 
            "DB_PATH", "BLAST_WORD_SIZE", "BLAST_EVALUE", "BLAST_MAX_TARGET_SEQS",
            "BLAST_REWARD", "BLAST_PENALTY", "BLAST_GAPOPEN", "BLAST_GAPEXTEND",
            "BLAST_FILTER_FACTOR"
        ],
        "ViennaRNA Parameters": [
            "THERMO_TEMPERATURE", "THERMO_SODIUM", "THERMO_MAGNESIUM"
        ],
        "Alignment Parameters": [
            "MIN_IDENTITY", "MIN_LENGTH", "LASTZ_OPTIONS"
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
        print(f"{Fore.GREEN}{category}{Style.RESET_ALL}")
        print(f"{Fore.GREEN}{'-' * len(category)}{Style.RESET_ALL}")
        
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
                
                print(f"{Fore.YELLOW}{key}{Style.RESET_ALL}: {formatted_value}")
        print()
    
    # Handle Primer3 settings separately - these are extensive
    print(f"{Fore.GREEN}Primer3 Settings{Style.RESET_ALL}")
    print(f"{Fore.GREEN}{'-' * 14}{Style.RESET_ALL}")
    print(f"{Fore.YELLOW}Use --config all to display all Primer3 settings{Style.RESET_ALL}\n")
    
    # Print footer with usage instructions
    print(f"{Fore.CYAN}{'='*80}{Style.RESET_ALL}")
    print(f"\n{Fore.WHITE}Configuration Options:{Style.RESET_ALL}")
    print(f"- View basic settings: {Fore.YELLOW}ddprimer --config{Style.RESET_ALL}")
    print(f"- View all settings (including Primer3): {Fore.YELLOW}ddprimer --config all{Style.RESET_ALL}")
    print(f"- Generate a template config file: {Fore.YELLOW}ddprimer --config template{Style.RESET_ALL}")
    print(f"- Generate a template in specific directory: {Fore.YELLOW}ddprimer --config template --output /path/to/dir{Style.RESET_ALL}")
    print(f"- Use custom config: {Fore.YELLOW}ddprimer --config your_config.json{Style.RESET_ALL}")
    print(f"\nExample config file format:")
    print(f"{Fore.BLUE}{{")
    print(f'    "PRIMER_MIN_SIZE": 18,')
    print(f'    "PRIMER_OPT_SIZE": 20,')
    print(f'    "PRIMER_MAX_SIZE": 25,')
    print(f'    "PRIMER_PRODUCT_SIZE_RANGE": [[90, 200]]')
    print(f"}}{Style.RESET_ALL}\n")

def display_primer3_settings(config_cls):
    """
    Display all Primer3 settings in a structured format.
    
    Args:
        config_cls: The Config class
    """
    # Initialize colorama for cross-platform colored output
    colorama.init()
    
    primer3_settings = config_cls.PRIMER3_SETTINGS
    
    # Print header
    print(f"\n{Fore.CYAN}{'='*80}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{'ddPrimer Primer3 Configuration Settings':^80}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{'='*80}{Style.RESET_ALL}\n")
    
    # Group Primer3 settings by categories
    categories = {
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
            print(f"{Fore.GREEN}{category}{Style.RESET_ALL}")
            print(f"{Fore.GREEN}{'-' * len(category)}{Style.RESET_ALL}")
            
            for key in sorted(keys):
                value = primer3_settings[key]
                print(f"{Fore.YELLOW}{key}{Style.RESET_ALL}: {value}")
            print()
    
    # Collect any uncategorized settings
    all_categorized = []
    for keys in categories.values():
        all_categorized.extend(keys)
    
    uncategorized = [key for key in primer3_settings.keys() if key not in all_categorized]
    
    if uncategorized:
        print(f"{Fore.GREEN}Other Primer3 Settings{Style.RESET_ALL}")
        print(f"{Fore.GREEN}{'-' * 20}{Style.RESET_ALL}")
        for key in sorted(uncategorized):
            value = primer3_settings[key]
            print(f"{Fore.YELLOW}{key}{Style.RESET_ALL}: {value}")
        print()
    
    # Print footer with usage instructions
    print(f"{Fore.CYAN}{'='*80}{Style.RESET_ALL}")
    print(f"\n{Fore.WHITE}To modify these Primer3 settings:{Style.RESET_ALL}")
    print(f"1. Create a Primer3 settings file or add a \"PRIMER3_SETTINGS\" section to your JSON config")
    print(f"2. Run: {Fore.YELLOW}ddprimer --config your_settings.txt{Style.RESET_ALL}")
    print(f"   or: {Fore.YELLOW}ddprimer --config your_config.json{Style.RESET_ALL}")
    print(f"\nTo view these settings again: {Fore.YELLOW}ddprimer --config all{Style.RESET_ALL}\n")