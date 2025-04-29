#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration template generator for the ddPrimer pipeline.
"""

import os
import json
import sys
from datetime import datetime
import colorama
from colorama import Fore, Style

def generate_config_template(config_cls, filename=None):
    """
    Generate a template configuration file based on current settings.
    
    Args:
        config_cls: The Config class
        filename (str, optional): Filename to save the template. If None, uses default name.
        
    Returns:
        str: Path to the generated template file
    """
    # Initialize colorama for cross-platform colored output
    colorama.init()
    
    # Get current working directory
    cwd = os.getcwd()
    
    # If no filename is provided, create a default one
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d")
        filename = f"ddprimer_config_template_{timestamp}.json"
    
    # Make sure filename has .json extension
    if not filename.lower().endswith('.json'):
        filename += '.json'
    
    # Create full path
    filepath = os.path.join(cwd, filename)
    
    # Helper function to safely get attributes
    def safe_get_attr(obj, attr, default=None):
        return getattr(obj, attr, default)
    
    # Create template dictionary with commonly modified settings
    template = {
        # Performance Settings
        "NUM_PROCESSES": safe_get_attr(config_cls, "NUM_PROCESSES", 4),
        "BATCH_SIZE": safe_get_attr(config_cls, "BATCH_SIZE", 100),
        "SHOW_PROGRESS": safe_get_attr(config_cls, "SHOW_PROGRESS", True),
        
        # Design Parameters
        "PRIMER_MIN_SIZE": safe_get_attr(config_cls, "PRIMER_MIN_SIZE", 18),
        "PRIMER_OPT_SIZE": safe_get_attr(config_cls, "PRIMER_OPT_SIZE", 20),
        "PRIMER_MAX_SIZE": safe_get_attr(config_cls, "PRIMER_MAX_SIZE", 23),
        "PRIMER_MIN_TM": safe_get_attr(config_cls, "PRIMER_MIN_TM", 50.0),
        "PRIMER_OPT_TM": safe_get_attr(config_cls, "PRIMER_OPT_TM", 57.5),
        "PRIMER_MAX_TM": safe_get_attr(config_cls, "PRIMER_MAX_TM", 65.0),
        "PRIMER_MIN_GC": safe_get_attr(config_cls, "PRIMER_MIN_GC", 50.0),
        "PRIMER_MAX_GC": safe_get_attr(config_cls, "PRIMER_MAX_GC", 60.0),
        "PRIMER_PRODUCT_SIZE_RANGE": safe_get_attr(config_cls, "PRIMER_PRODUCT_SIZE_RANGE", [[90, 200]]),
        "MAX_PRIMER_PAIRS_PER_SEGMENT": safe_get_attr(config_cls, "MAX_PRIMER_PAIRS_PER_SEGMENT", 3),
    }
    
    # Add validation options if they exist
    validation_mode = safe_get_attr(config_cls, "VALIDATION_MODE", None)
    if validation_mode is not None:
        template["VALIDATION_MODE"] = validation_mode
    
    amp_mismatches = safe_get_attr(config_cls, "ALLOW_AMP_MISMATCHES", None)
    if amp_mismatches is not None:
        template["ALLOW_AMP_MISMATCHES"] = amp_mismatches
    
    # Add BLAST settings if they exist
    db_path = safe_get_attr(config_cls, "DB_PATH", None)
    if db_path is not None:
        template["DB_PATH"] = db_path
    
    # Get Primer3 settings if they exist
    primer3_settings = safe_get_attr(config_cls, "PRIMER3_SETTINGS", {})
    if primer3_settings:
        # Include only a subset of Primer3 settings
        p3_subset = {}
        important_p3_settings = [
            "PRIMER_PICK_INTERNAL_OLIGO",
            "PRIMER_GC_CLAMP",
            "PRIMER_MAX_POLY_X",
            "PRIMER_PAIR_MAX_DIFF_TM"
        ]
        
        for key in important_p3_settings:
            if key in primer3_settings:
                p3_subset[key] = primer3_settings[key]
        
        if p3_subset:
            template["PRIMER3_SETTINGS"] = p3_subset
    
    # Add comments as string values at the top of the template
    template_with_comments = {
        "# ddPrimer Configuration Template": "Generated on " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "# Instructions": "Modify the values below and save this file. Use it with: ddprimer --config your_config.json",
        "# For full settings list": "Run 'ddprimer --config all'",
        "# Documentation": "Visit https://github.com/yourusername/ddprimer for full documentation",
    }
    
    # Merge comments with actual settings
    template_with_comments.update(template)
    
    # Write to file
    try:
        with open(filepath, 'w') as f:
            json.dump(template_with_comments, f, indent=4)
        
        print(f"\n{Fore.GREEN}Configuration template generated successfully!{Style.RESET_ALL}")
        print(f"Template saved to: {Fore.YELLOW}{filepath}{Style.RESET_ALL}")
        print(f"\nTo use this template:")
        print(f"1. Edit the file with your preferred settings")
        print(f"2. Run: {Fore.YELLOW}ddprimer --config {filename}{Style.RESET_ALL}")
        print()
        
        return filepath
    except Exception as e:
        print(f"\n{Fore.RED}Error generating template: {str(e)}{Style.RESET_ALL}")
        return None