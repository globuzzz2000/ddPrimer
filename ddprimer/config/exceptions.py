#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Custom exceptions for ddPrimer pipeline.

This module defines exception classes used throughout the ddPrimer pipeline
to provide more specific error information and improve error handling.
"""


class DDPrimerError(Exception):
    """Base exception class for all ddPrimer-specific errors."""
    pass


class FileError(DDPrimerError):
    """Base class for file-related errors."""
    pass


class FileSelectionError(FileError):
    """Error during file selection via GUI or CLI."""
    pass


class FileFormatError(FileError):
    """Error with file formatting or parsing."""
    pass


class ConfigError(DDPrimerError):
    """Error with configuration parameters."""
    pass


class SequenceProcessingError(DDPrimerError):
    """Error during sequence processing."""
    pass


class BlastError(DDPrimerError):
    """Base class for BLAST-related errors."""
    pass


class BlastDBError(BlastError):
    """Error with BLAST database creation or access."""
    pass


class BlastExecutionError(BlastError):
    """Error when executing BLAST commands."""
    pass



class Primer3Error(SequenceProcessingError):
    """Error during Primer3 execution or parsing."""
    pass


class SNPVerificationError(DDPrimerError):
    """Error during SNP verification or checking."""
    pass


class PrimerDesignError(DDPrimerError):
    """Error during primer design process."""
    pass


class ValidationError(DDPrimerError):
    """Error during validation of primers or parameters."""
    pass


class AlignmentError(DDPrimerError):
    """Error during sequence alignment."""
    pass


class WorkflowError(DDPrimerError):
    """Error in workflow execution."""
    pass


class ExternalToolError(DDPrimerError):
    """Error related to external tools like Primer3, Vienna, etc."""
    
    def __init__(self, message, tool_name=None, command=None, return_code=None, stdout=None, stderr=None):
        """
        Initialize with extended information about the external tool error.
        
        Args:
            message (str): Error message
            tool_name (str, optional): Name of the external tool
            command (str, optional): Command that was executed
            return_code (int, optional): Return code from the command
            stdout (str, optional): Standard output from the command
            stderr (str, optional): Standard error from the command
        """
        self.tool_name = tool_name
        self.command = command
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr
        
        detailed_message = message
        if tool_name:
            detailed_message = f"{tool_name} error: {message}"
        if return_code is not None:
            detailed_message += f" (return code: {return_code})"
            
        super().__init__(detailed_message)