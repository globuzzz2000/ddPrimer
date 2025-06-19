#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Logging configuration module for ddPrimer pipeline.

Contains functionality for:
1. Module-specific debug level control based on filenames
2. Enhanced debug filtering and formatting with colors  
3. Singleton pattern for consistent logging configuration
4. Automatic log file management and rotation

This module provides comprehensive logging configuration for the ddPrimer
pipeline, enabling granular debug control per module and enhanced
debugging capabilities with colored output.
"""

import os
import logging
from datetime import datetime
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class ModuleDebugConfig:
    """
    Configuration for module-specific debug settings.
    
    This class manages debug level configuration for individual modules,
    providing filename-based debug control for intuitive usage and
    comprehensive module coverage.
    
    Attributes:
        MODULE_DEBUG_LEVELS: Dictionary mapping module names to log levels
        FILENAME_TO_MODULE: Dictionary mapping filenames to full module paths
        
    Example:
        >>> config = ModuleDebugConfig()
        >>> config.MODULE_DEBUG_LEVELS['ddprimer.pipeline']
        20
    """
    
    # All modules default to INFO level
    MODULE_DEBUG_LEVELS = {
        'ddprimer.pipeline': logging.INFO,
        'ddprimer.modes.standard_mode': logging.INFO,
        'ddprimer.modes.direct_mode': logging.INFO,
        'ddprimer.modes.alignment_mode': logging.INFO,
        'ddprimer.modes.common': logging.INFO,
        'ddprimer.core.snp_processor': logging.INFO,
        'ddprimer.core.annotation_processor': logging.INFO,
        'ddprimer.core.blast_processor': logging.INFO,
        'ddprimer.core.vienna_processor': logging.INFO,
        'ddprimer.core.primer3_processor': logging.INFO,
        'ddprimer.core.primer_processor': logging.INFO,
        'ddprimer.core.sequence_processor': logging.INFO,
        'ddprimer.helpers.alignment_workflow': logging.INFO,
        'ddprimer.helpers.direct_masking': logging.INFO,
        'ddprimer.helpers.lastz_runner': logging.INFO,
        'ddprimer.helpers.maf_parser': logging.INFO,
        'ddprimer.utils.blast_db_creator': logging.INFO,
        'ddprimer.utils.blast_verification': logging.INFO,
        'ddprimer.utils.chromosome_mapper': logging.INFO,
        'ddprimer.utils.common_utils': logging.INFO,
        'ddprimer.utils.db_selector': logging.INFO,
        'ddprimer.utils.file_io': logging.INFO,
        'ddprimer.utils.model_organism_manager': logging.INFO,
        'ddprimer.utils.sequence_utils': logging.INFO,
        'ddprimer.config.config': logging.INFO,
        'ddprimer.config.config_display': logging.INFO,
        'ddprimer.config.template_generator': logging.INFO,
    }
    
    # Filename to module mapping for intuitive usage
    FILENAME_TO_MODULE = {
        # Main pipeline
        'pipeline': 'ddprimer.pipeline',
        
        # Mode files
        'standard_mode': 'ddprimer.modes.standard_mode',
        'direct_mode': 'ddprimer.modes.direct_mode', 
        'alignment_mode': 'ddprimer.modes.alignment_mode',
        'common': 'ddprimer.modes.common',
        
        # Core processing files
        'snp_processor': 'ddprimer.core.snp_processor',
        'annotation_processor': 'ddprimer.core.annotation_processor',
        'blast_processor': 'ddprimer.core.blast_processor',
        'vienna_processor': 'ddprimer.core.vienna_processor',
        'primer3_processor': 'ddprimer.core.primer3_processor',
        'primer_processor': 'ddprimer.core.primer_processor',
        'sequence_processor': 'ddprimer.core.sequence_processor',
        
        # Helper files
        'alignment_workflow': 'ddprimer.helpers.alignment_workflow',
        'direct_masking': 'ddprimer.helpers.direct_masking',
        'lastz_runner': 'ddprimer.helpers.lastz_runner',
        'maf_parser': 'ddprimer.helpers.maf_parser',
        
        # Utility files
        'blast_db_creator': 'ddprimer.utils.blast_db_creator',
        'blast_verification': 'ddprimer.utils.blast_verification',
        'chromosome_mapper': 'ddprimer.utils.chromosome_mapper',
        'common_utils': 'ddprimer.utils.common_utils',
        'db_selector': 'ddprimer.utils.db_selector',
        'file_io': 'ddprimer.utils.file_io',
        'model_organism_manager': 'ddprimer.utils.model_organism_manager',
        'sequence_utils': 'ddprimer.utils.sequence_utils',
        
        # Config files  
        'config': 'ddprimer.config.config',
        'config_display': 'ddprimer.config.config_display',
        'template_generator': 'ddprimer.config.template_generator',
    }


class SimpleDebugFormatter(logging.Formatter):
    """
    Simple formatter with colors for debug output.
    
    Provides color-coded log output for improved readability during
    debugging, with support for different log levels and optional
    color disabling for environments that don't support ANSI colors.
    
    Attributes:
        use_colors: Whether to use ANSI color codes in output
        COLORS: Dictionary mapping log levels to ANSI color codes
        
    Example:
        >>> formatter = SimpleDebugFormatter(use_colors=True)
        >>> handler.setFormatter(formatter)
    """
    
    def __init__(self, fmt=None, datefmt=None, use_colors=False):
        """
        Initialize formatter with optional color support.
        
        Args:
            fmt: Log message format string
            datefmt: Date format string
            use_colors: Whether to enable ANSI color codes
        """
        super().__init__(fmt, datefmt)
        self.use_colors = use_colors
        
        # ANSI color codes
        self.COLORS = {
            'DEBUG': '\033[37m',      # White
            'INFO': '\033[1;37m',     # Bold White
            'WARNING': '\033[33m',    # Yellow
            'ERROR': '\033[31m',      # Red
            'CRITICAL': '\033[35m',   # Magenta
        }
        self.RESET = '\033[0m'
    
    def format(self, record):
        """
        Format the log record with colors.
        
        Args:
            record: LogRecord instance to format
            
        Returns:
            str: Formatted log message with optional colors
        """
        formatted_message = super().format(record)
        
        # Add colors if enabled
        if self.use_colors and record.levelname in self.COLORS:
            color_code = self.COLORS[record.levelname]
            formatted_message = f"{color_code}{formatted_message}{self.RESET}"
        
        return formatted_message


class EnhancedDebugFilter(logging.Filter):
    """
    Filter that controls module-specific debug output.
    
    Provides granular control over debug output by filtering log records
    based on module-specific log levels, allowing users to enable debug
    mode for specific modules while keeping others at INFO level.
    
    Attributes:
        module_levels: Dictionary mapping module names to minimum log levels
        
    Example:
        >>> filter_obj = EnhancedDebugFilter({'ddprimer.pipeline': logging.DEBUG})
        >>> handler.addFilter(filter_obj)
    """
    
    def __init__(self, module_levels: Dict[str, int]):
        """
        Initialize filter with module-specific log levels.
        
        Args:
            module_levels: Dictionary mapping module names to minimum log levels
        """
        super().__init__()
        self.module_levels = module_levels
    
    def filter(self, record):
        """
        Filter records based on module-specific levels.
        
        Args:
            record: LogRecord instance to filter
            
        Returns:
            bool: True if record should be logged, False otherwise
        """
        module_name = record.name
        
        # Try exact match first
        if module_name in self.module_levels:
            return record.levelno >= self.module_levels[module_name]
        
        # Try parent module matches
        for module_pattern, level in self.module_levels.items():
            if module_name.startswith(module_pattern + '.') or module_name == module_pattern:
                return record.levelno >= level
        
        # For modules not in our config, default to INFO level
        return record.levelno >= logging.INFO


def setup_logging(debug: Union[bool, List[str], str] = False) -> str:
    """
    Configure logging with filename-based debug control.
    
    Sets up comprehensive logging configuration with support for module-specific
    debug levels, colored console output, and automatic log file management.
    
    Args:
        debug: Debug configuration options:
               - False: No debug logging
               - True: Universal debug for all modules
               - str: Single filename for debug (e.g., 'pipeline', 'SNP_processor')
               - List[str]: List of filenames for debug
        
    Returns:
        str: Path to the created log file
        
    Raises:
        LoggingConfigError: If logging setup fails
        
    Example:
        >>> log_file = setup_logging(debug=['pipeline', 'blast_processor'])
        >>> print(f"Logs saved to: {log_file}")
    """
    logger.debug("=== LOGGING SETUP DEBUG ===")
    logger.debug(f"Setting up logging with debug={debug}")
    
    try:
        from .config import Config
        
        # Normalize debug input
        debug_enabled, debug_modules = _normalize_debug_input(debug)
        logger.debug(f"Normalized debug: enabled={debug_enabled}, modules={debug_modules}")
        
        # Start with default module configuration
        module_config = ModuleDebugConfig.MODULE_DEBUG_LEVELS.copy()
        
        # Apply module-specific debug configuration
        if debug_modules:
            # Start with INFO level for all modules
            for module in module_config:
                module_config[module] = logging.INFO
            
            # Enable DEBUG for specified modules
            for module_name in debug_modules:
                full_module_name = _resolve_module_name(module_name)
                if full_module_name and full_module_name in module_config:
                    module_config[full_module_name] = logging.DEBUG
                    logger.debug(f"Enabled DEBUG for {full_module_name}")
        elif debug_enabled:
            # Universal debug - enable DEBUG for all modules
            for module in module_config:
                module_config[module] = logging.DEBUG
            logger.debug("Enabled DEBUG for all modules")
        
        # Set up log directory
        log_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer", "logs")
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, f"ddPrimer_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        
        # Configure root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)
        
        # Clear existing handlers
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)
        
        # File handler (detailed, no colors)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG if debug_enabled else logging.INFO)
        file_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
        )
        file_handler.setFormatter(file_formatter)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG if debug_enabled else logging.INFO)
        
        if debug_enabled:
            # Enhanced debug format with colors
            console_formatter = SimpleDebugFormatter(
                fmt='%(levelname)-8s [%(name)s] %(message)s',
                use_colors=True
            )
            
            # Add module-specific filter
            debug_filter = EnhancedDebugFilter(module_config)
            console_handler.addFilter(debug_filter)
        else:
            # Simple format for normal operation
            console_formatter = SimpleDebugFormatter(
                fmt='%(message)s',
                use_colors=False
            )
        
        console_handler.setFormatter(console_formatter)
        
        # Add handlers to root logger
        root_logger.addHandler(file_handler)
        root_logger.addHandler(console_handler)
        
        # Configure main logger
        main_logger = logging.getLogger("ddPrimer")
        
        if debug_enabled:
            main_logger.debug("Debug logging enabled")
            main_logger.debug(f"Log file: {log_file}")
            if debug_modules:
                resolved_modules = []
                unresolved_modules = []
                for module_name in debug_modules:
                    full_module_name = _resolve_module_name(module_name)
                    if full_module_name:
                        resolved_modules.append(f"{module_name} → {full_module_name}")
                    else:
                        unresolved_modules.append(module_name)
                
                if resolved_modules:
                    main_logger.debug(f"Debug enabled for modules: {', '.join(resolved_modules)}")
                if unresolved_modules:
                    main_logger.warning(f"Unknown module names: {', '.join(unresolved_modules)}")
            else:
                main_logger.debug("Debug enabled for all modules")
        
        # Update Config
        Config.DEBUG_MODE = debug_enabled
        logger.debug("Updated Config.DEBUG_MODE")
        
        logger.debug("Logging setup completed successfully")
        return log_file
        
    except Exception as e:
        error_msg = f"Failed to setup logging configuration"
        print(f"ERROR: {error_msg}: {str(e)}")  # Can't use logger here since setup failed
        raise LoggingConfigError(error_msg) from e
    
    logger.debug("=== END LOGGING SETUP DEBUG ===")


def _normalize_debug_input(debug: Union[bool, List[str], str]) -> tuple[bool, Optional[List[str]]]:
    """
    Normalize various debug input formats to (debug_enabled, debug_modules).
    
    Args:
        debug: Debug input in various formats
        
    Returns:
        Tuple of (debug_enabled: bool, debug_modules: Optional[List[str]])
        
    Example:
        >>> enabled, modules = _normalize_debug_input(['pipeline', 'blast'])
        >>> print(enabled, modules)
        True ['pipeline', 'blast']
    """
    try:
        from .config import Config
        
        if isinstance(debug, bool):
            debug_enabled = debug or Config.DEBUG_MODE
            debug_modules = None
        elif isinstance(debug, str):
            # Single module name as string
            debug_enabled = True
            debug_modules = [debug]
        elif isinstance(debug, list):
            # List of module names
            debug_enabled = True
            debug_modules = debug
        else:
            # Unknown format, default to False
            debug_enabled = False
            debug_modules = None
        
        return debug_enabled, debug_modules
        
    except Exception as e:
        error_msg = f"Failed to normalize debug input: {debug}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise LoggingConfigError(error_msg) from e


def _resolve_module_name(filename: str) -> Optional[str]:
    """
    Resolve a filename to its full module path.
    
    Args:
        filename: Filename like 'pipeline', 'snp_processor', etc.
        
    Returns:
        Full module path or None if not found
        
    Example:
        >>> module_name = _resolve_module_name('pipeline')
        >>> print(module_name)
        ddprimer.pipeline
    """
    # Try exact match first
    if filename in ModuleDebugConfig.FILENAME_TO_MODULE:
        return ModuleDebugConfig.FILENAME_TO_MODULE[filename]
    
    # If it's already a full module name, check if it exists
    if filename in ModuleDebugConfig.MODULE_DEBUG_LEVELS:
        return filename
        
    return None


class LoggingConfigError(Exception):
    """Error during logging configuration setup."""
    pass