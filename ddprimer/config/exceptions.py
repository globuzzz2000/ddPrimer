#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exception handling module for the ddPrimer pipeline.
"""

from typing import Optional


class PrimerDesignError(Exception):
    """Base exception class for all primer design errors."""
    
    def __init__(self, message: str = "An error occurred in the primer design pipeline"):
        self.message = message
        super().__init__(self.message)


class FileError(PrimerDesignError):
    """Exception raised for file-related errors."""
    
    def __init__(self, message: str = "File error", filename: Optional[str] = None):
        self.filename = filename
        msg = message
        if filename:
            msg = f"{message}: '{filename}'"
        super().__init__(msg)


class FileFormatError(FileError):
    """Exception raised for file format errors."""
    
    def __init__(self, message: str = "Invalid file format", filename: Optional[str] = None):
        super().__init__(message, filename)


class FileMissingError(FileError):
    """Exception raised when a required file is missing."""
    
    def __init__(self, message: str = "Required file is missing", filename: Optional[str] = None):
        super().__init__(message, filename)


class SequenceError(PrimerDesignError):
    """Exception raised for sequence-related errors."""
    
    def __init__(self, message: str = "Sequence error", sequence_id: Optional[str] = None):
        self.sequence_id = sequence_id
        msg = message
        if sequence_id:
            msg = f"{message} for sequence '{sequence_id}'"
        super().__init__(msg)


class PrimerError(PrimerDesignError):
    """Exception raised for primer-related errors."""
    
    def __init__(self, message: str = "Primer error", primer_id: Optional[str] = None):
        self.primer_id = primer_id
        msg = message
        if primer_id:
            msg = f"{message} for primer '{primer_id}'"
        super().__init__(msg)


class Primer3Error(PrimerDesignError):
    """Exception raised for Primer3 execution errors."""
    
    def __init__(self, message: str = "Primer3 error", details: Optional[str] = None):
        self.details = details
        msg = message
        if details:
            msg = f"{message}: {details}"
        super().__init__(msg)


class BlastError(PrimerDesignError):
    """Exception raised for BLAST-related errors."""
    
    def __init__(self, message: str = "BLAST error", details: Optional[str] = None):
        self.details = details
        msg = message
        if details:
            msg = f"{message}: {details}"
        super().__init__(msg)


class ConfigError(PrimerDesignError):
    """Exception raised for configuration errors."""
    
    def __init__(self, message: str = "Configuration error", setting: Optional[str] = None):
        self.setting = setting
        msg = message
        if setting:
            msg = f"{message}: '{setting}'"
        super().__init__(msg)


class ValidationError(PrimerDesignError):
    """Exception raised for validation errors."""
    
    def __init__(self, message: str = "Validation error", details: Optional[str] = None):
        self.details = details
        msg = message
        if details:
            msg = f"{message}: {details}"
        super().__init__(msg)