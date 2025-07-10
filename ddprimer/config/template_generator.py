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
        
        # Generate default filename if none provided
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d")
            filename = f"ddprimer_config_template_{timestamp}.json"
            logger.debug(f"Generated default filename: {filename}")
        
        # Ensure filename has .json extension
        if not filename.lower().endswith('.json'):
            filename += '.json'
            logger.debug(f"Added .json extension: {filename}")
        
        # Create full path
        filepath = os.path.join(output_dir, filename)
        logger.debug(f"Full template path: {filepath}")
        
        # Create template dictionary with commonly modified settings
        template = _build_template_dict(config_cls)
        logger.debug(f"Built template with {len(template)} settings")
        
        # Write to file
        try:
            with open(filepath, 'w') as f:
                json.dump(template, f, indent=4)
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
            error_msg = f"Template file creation failed at {filepath}: {str(e)}"
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
            
            # Pipeline parameters
            "MIN_SEGMENT_LENGTH": safe_get_attr(config_cls, "MIN_SEGMENT_LENGTH", 90),
            "RETAIN_TYPES": safe_get_attr(config_cls, "RETAIN_TYPES", ["gene"]),
            "GENE_OVERLAP_MARGIN": safe_get_attr(config_cls, "GENE_OVERLAP_MARGIN", 25),
            "RESTRICTION_SITE": safe_get_attr(config_cls, "RESTRICTION_SITE", "GGCC"),
            "PENALTY_MAX": safe_get_attr(config_cls, "PENALTY_MAX", 5.0),
            "MAX_PRIMER_PAIRS_PER_SEGMENT": safe_get_attr(config_cls, "MAX_PRIMER_PAIRS_PER_SEGMENT", 3),
            "PREFER_PROBE_MORE_C_THAN_G": safe_get_attr(config_cls, "PREFER_PROBE_MORE_C_THAN_G", True),
            "SEQUENCE_MIN_GC": safe_get_attr(config_cls, "SEQUENCE_MIN_GC", 50.0),
            "SEQUENCE_MAX_GC": safe_get_attr(config_cls, "SEQUENCE_MAX_GC", 60.0),
            "MIN_CHROMOSOME_SIZE": safe_get_attr(config_cls, "MIN_CHROMOSOME_SIZE", 2000000),
            
            # VCF/SNP Processing parameters
            "VCF_ALLELE_FREQUENCY_THRESHOLD": safe_get_attr(config_cls, "VCF_ALLELE_FREQUENCY_THRESHOLD", None),
            "VCF_QUALITY_THRESHOLD": safe_get_attr(config_cls, "VCF_QUALITY_THRESHOLD", None),
            "VCF_FLANKING_MASK_SIZE": safe_get_attr(config_cls, "VCF_FLANKING_MASK_SIZE", 0),
            "VCF_USE_SOFT_MASKING": safe_get_attr(config_cls, "VCF_USE_SOFT_MASKING", False),
            
            # ViennaRNA Parameters
            "THERMO_TEMPERATURE": safe_get_attr(config_cls, "THERMO_TEMPERATURE", 37),
            "THERMO_SODIUM": safe_get_attr(config_cls, "THERMO_SODIUM", 0.05),
            "THERMO_MAGNESIUM": safe_get_attr(config_cls, "THERMO_MAGNESIUM", 0.0),
            
            # BLAST Database Options
            "DB_PATH": safe_get_attr(config_cls, "DB_PATH", None),
            "USE_CUSTOM_DB": safe_get_attr(config_cls, "USE_CUSTOM_DB", False),
            "BLAST_WORD_SIZE": safe_get_attr(config_cls, "BLAST_WORD_SIZE", 7),
            "BLAST_EVALUE": safe_get_attr(config_cls, "BLAST_EVALUE", 10),
            "BLAST_MAX_TARGET_SEQS": safe_get_attr(config_cls, "BLAST_MAX_TARGET_SEQS", 100),
            "BLAST_REWARD": safe_get_attr(config_cls, "BLAST_REWARD", 2),
            "BLAST_PENALTY": safe_get_attr(config_cls, "BLAST_PENALTY", -3),
            "BLAST_GAPOPEN": safe_get_attr(config_cls, "BLAST_GAPOPEN", 5),
            "BLAST_GAPEXTEND": safe_get_attr(config_cls, "BLAST_GAPEXTEND", 2),
            "BLAST_FILTER_FACTOR": safe_get_attr(config_cls, "BLAST_FILTER_FACTOR", 100),
        }
        
        # Get Primer3 settings if they exist
        primer3_settings = safe_get_attr(config_cls, "PRIMER3_SETTINGS", {})
        if primer3_settings:
            # Include only a subset of Primer3 settings that users commonly modify
            p3_subset = {}
            important_p3_settings = [
                "PRIMER_PICK_INTERNAL_OLIGO",
                "PRIMER_GC_CLAMP",
                "PRIMER_MAX_POLY_X",
                "PRIMER_PAIR_MAX_DIFF_TM",
                "PRIMER_SALT_MONOVALENT",
                "PRIMER_SALT_DIVALENT",
                "PRIMER_ANNEALING_TEMP",
                "PRIMER_DNA_CONC"
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