#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration template generator for ddPrimer pipeline.

Contains functionality for:
1. JSON configuration template generation with default values
2. Automatic directory creation and file management
3. Color-coded console output for user guidance
4. Comprehensive configuration documentation

This module provides template generation capabilities for the ddPrimer
pipeline, allowing users to create customizable configuration files
with appropriate defaults and documentation.
"""

import os
import json
from datetime import datetime
import colorama
from colorama import Fore, Style
import logging

logger = logging.getLogger(__name__)


def generate_config_template(config_cls, filename=None, output_dir=None):
    """
    Generate a template configuration file based on current settings.
    
    Creates a JSON configuration template with commonly modified settings
    and helpful comments to guide users in customizing their ddPrimer
    configuration.
    
    Args:
        config_cls: The Config class containing default settings
        filename (str, optional): Filename to save template. Uses default if None.
        output_dir (str, optional): Directory to save template. Uses current if None.
        
    Returns:
        str: Path to the generated template file, or None if generation failed
        
    Raises:
        TemplateGenerationError: If template generation fails
        
    Example:
        >>> from ddprimer.config import Config
        >>> template_path = generate_config_template(Config, "my_config.json")
        >>> print(f"Template created: {template_path}")
    """
    logger.debug("=== TEMPLATE GENERATION DEBUG ===")
    logger.debug(f"Generating config template: filename={filename}, output_dir={output_dir}")
    
    try:
        # Initialize colorama for cross-platform colored output
        colorama.init()
        
        # Determine output directory
        if output_dir is None:
            output_dir = os.getcwd()
            logger.debug("Using current directory as output")
        else:
            # Create the directory if it doesn't exist
            if not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir)
                    print(f"{Fore.GREEN}Created output directory: {output_dir}{Style.RESET_ALL}")
                    logger.debug(f"Created output directory: {output_dir}")
                except Exception as e:
                    error_msg = f"Failed to create output directory: {output_dir}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    print(f"{Fore.RED}Error creating output directory: {str(e)}{Style.RESET_ALL}")
                    print(f"{Fore.YELLOW}Using current directory instead.{Style.RESET_ALL}")
                    output_dir = os.getcwd()
        
        # If no filename is provided, create a default one
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d")
            filename = f"ddprimer_config_template_{timestamp}.json"
            logger.debug(f"Generated default filename: {filename}")
        
        # Make sure filename has .json extension
        if not filename.lower().endswith('.json'):
            filename += '.json'
            logger.debug(f"Added .json extension: {filename}")
        
        # Create full path
        filepath = os.path.join(output_dir, filename)
        logger.debug(f"Full template path: {filepath}")
        
        # Create template dictionary with commonly modified settings
        template = _build_template_dict(config_cls)
        logger.debug(f"Built template with {len(template)} settings")
        
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
            print(f"\n{Fore.WHITE}{'='*80}{Style.RESET_ALL}")
            print(f"{Fore.WHITE}Configuration Template File Generator")
            print(f"{Fore.WHITE}{'='*80}{Style.RESET_ALL}")
            print(f"\nTemplate saved to: {Fore.CYAN}{filepath}{Style.RESET_ALL}")
            print(f"\nTo use this template:")
            print(f"1. Edit the file with your preferred settings")
            print(f"2. Run: {Fore.CYAN}ddprimer --config {filepath}{Style.RESET_ALL}")
            print(f"\n{Fore.WHITE}{'='*80}{Style.RESET_ALL}\n")
            
            logger.debug("Template generation completed successfully")
            return filepath
            
        except Exception as e:
            error_msg = f"Failed to write template file to {filepath}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            print(f"\n{Fore.RED}Error generating template: {str(e)}{Style.RESET_ALL}")
            raise TemplateGenerationError(error_msg) from e
            
    except Exception as e:
        error_msg = f"Template generation failed"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        if not isinstance(e, TemplateGenerationError):
            raise TemplateGenerationError(error_msg) from e
        raise
    
    logger.debug("=== END TEMPLATE GENERATION DEBUG ===")


def _build_template_dict(config_cls):
    """
    Build the template dictionary with commonly modified settings.
    
    Args:
        config_cls: The Config class containing default settings
        
    Returns:
        dict: Template dictionary with organized settings
        
    Raises:
        TemplateGenerationError: If template building fails
    """
    logger.debug("Building template dictionary from config class")
    
    try:
        # Helper function to safely get attributes
        def safe_get_attr(obj, attr, default=None):
            """Safely get attribute value with fallback to default."""
            try:
                return getattr(obj, attr, default)
            except Exception as e:
                logger.debug(f"Failed to get attribute {attr}: {str(e)}")
                return default
        
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
            logger.debug("Added VALIDATION_MODE to template")
        
        amp_mismatches = safe_get_attr(config_cls, "ALLOW_AMP_MISMATCHES", None)
        if amp_mismatches is not None:
            template["ALLOW_AMP_MISMATCHES"] = amp_mismatches
            logger.debug("Added ALLOW_AMP_MISMATCHES to template")
        
        # Add BLAST settings if they exist
        db_path = safe_get_attr(config_cls, "DB_PATH", None)
        if db_path is not None:
            template["DB_PATH"] = db_path
            logger.debug("Added DB_PATH to template")
        
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
                logger.debug(f"Added {len(p3_subset)} Primer3 settings to template")
        
        logger.debug(f"Template dictionary built with {len(template)} settings")
        return template
        
    except Exception as e:
        error_msg = f"Failed to build template dictionary"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise TemplateGenerationError(error_msg) from e


class TemplateGenerationError(Exception):
    """Error during configuration template generation."""
    pass