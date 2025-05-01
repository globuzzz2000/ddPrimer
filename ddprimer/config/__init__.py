#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration package for the ddPrimer pipeline.
"""

# In ddprimer/config/__init__.py
from .config import Config
from .logging_config import setup_logging
from .config_display import display_config, display_primer3_settings
from .template_generator import generate_config_template

# Import all exceptions
from .exceptions import (
    DDPrimerError,
    FileError, 
    FileSelectionError, 
    FileFormatError,
    ConfigError,
    SequenceProcessingError,
    BlastError,
    BlastDBError,
    BlastExecutionError,
    NupackError,
    Primer3Error,
    SNPVerificationError,
    PrimerDesignError,
    ValidationError,
    AlignmentError,
    WorkflowError,
    ExternalToolError
)

__all__ = [
    # Configuration
    'Config', 
    'setup_logging', 
    'display_config', 
    'display_primer3_settings',
    'generate_config_template',
    
    # Exceptions
    'DDPrimerError',
    'FileError',
    'FileSelectionError',
    'FileFormatError',
    'ConfigError',
    'SequenceProcessingError',
    'BlastError',
    'BlastDBError',
    'BlastExecutionError',
    'NupackError',
    'Primer3Error',
    'SNPVerificationError',
    'PrimerDesignError',
    'ValidationError',
    'AlignmentError',
    'WorkflowError',
    'ExternalToolError'
]