#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Custom exceptions for ddPrimer pipeline.

Contains functionality for:
1. Hierarchical exception structure for different error types
2. Enhanced error information with context and debugging data
3. External tool error handling with command details
4. Pipeline-specific error categorization

This module defines exception classes used throughout the ddPrimer pipeline
to provide more specific error information and improve error handling
across all components of the primer design workflow.
"""


class DDPrimerError(Exception):
    """
    Base exception class for all ddPrimer-specific errors.
    
    This is the root exception class from which all other ddPrimer
    exceptions inherit. It provides a common base for catching
    any pipeline-related error.
    
    Example:
        >>> try:
        ...     # Some ddPrimer operation
        ...     pass
        ... except DDPrimerError as e:
        ...     print(f"ddPrimer error occurred: {e}")
    """
    pass


class FileError(DDPrimerError):
    """
    Base class for file-related errors.
    
    Used for errors involving file operations, including reading,
    writing, parsing, and file access issues. This serves as a
    parent class for more specific file-related exceptions.
    
    Example:
        >>> raise FileError("Unable to read input file: permission denied")
    """
    pass


class FileSelectionError(FileError):
    """
    Error during file selection via GUI or CLI.
    
    Raised when there are issues with the file selection process,
    including GUI failures, invalid file paths, or user cancellation
    in interactive file selection modes.
    
    Example:
        >>> raise FileSelectionError("File selection dialog failed to open")
    """
    pass


class FileFormatError(FileError):
    """
    Error with file formatting or parsing.
    
    Raised when input files have invalid formats, corrupted content,
    or cannot be parsed according to their expected structure
    (FASTA, VCF, GFF, etc.).
    
    Example:
        >>> raise FileFormatError("Invalid FASTA format: missing sequence data")
    """
    pass


class ConfigError(DDPrimerError):
    """
    Error with configuration parameters.
    
    Raised when there are issues with configuration settings,
    invalid parameter values, or problems loading/saving
    configuration files.
    
    Example:
        >>> raise ConfigError("Invalid primer size range: min > max")
    """
    pass


class SequenceProcessingError(DDPrimerError):
    """
    Error during sequence processing.
    
    Raised when there are issues with DNA sequence operations,
    including sequence validation, masking, filtering, or
    any sequence-specific processing tasks.
    
    Example:
        >>> raise SequenceProcessingError("Invalid nucleotide characters in sequence")
    """
    pass


class BlastError(DDPrimerError):
    """
    Base class for BLAST-related errors.
    
    Parent class for all BLAST-specific errors, including
    database issues and execution problems. Use more specific
    subclasses when possible.
    
    Example:
        >>> raise BlastError("BLAST operation failed")
    """
    pass


class BlastDBError(BlastError):
    """
    Error with BLAST database creation or access.
    
    Raised when there are problems creating, accessing, or
    validating BLAST databases, including missing database
    files or corruption issues.
    
    Example:
        >>> raise BlastDBError("BLAST database files are missing or corrupted")
    """
    pass


class BlastExecutionError(BlastError):
    """
    Error when executing BLAST commands.
    
    Raised when BLAST commands fail to execute properly,
    including timeout issues, invalid parameters, or
    unexpected BLAST tool behavior.
    
    Example:
        >>> raise BlastExecutionError("BLAST search timed out after 300 seconds")
    """
    pass


class Primer3Error(SequenceProcessingError):
    """
    Error during Primer3 execution or parsing.
    
    Raised when there are issues with Primer3 primer design,
    including execution failures, parameter errors, or
    problems parsing Primer3 output.
    
    Example:
        >>> raise Primer3Error("Primer3 failed to design primers for sequence")
    """
    pass


class SNPVerificationError(DDPrimerError):
    """
    Error during SNP verification or checking.
    
    Raised when there are issues with SNP processing,
    including VCF file parsing errors, chromosome mapping
    problems, or SNP validation failures.
    
    Example:
        >>> raise SNPVerificationError("VCF chromosome names do not match FASTA sequences")
    """
    pass


class PrimerDesignError(DDPrimerError):
    """
    Error during primer design process.
    
    Raised when there are issues with the overall primer
    design workflow, including filtering failures, quality
    assessment problems, or design constraint violations.
    
    Example:
        >>> raise PrimerDesignError("No primers passed quality filtering criteria")
    """
    pass


class ValidationError(DDPrimerError):
    """
    Error during validation of primers or parameters.
    
    Raised when validation checks fail for primers, sequences,
    or pipeline parameters, including constraint violations
    or quality threshold failures.
    
    Example:
        >>> raise ValidationError("Primer melting temperature outside acceptable range")
    """
    pass


class AlignmentError(DDPrimerError):
    """
    Error during sequence alignment.
    
    Raised when there are issues with sequence alignment
    operations, including LastZ execution failures, MAF
    file parsing errors, or alignment processing problems.
    
    Example:
        >>> raise AlignmentError("LastZ alignment failed: insufficient memory")
    """
    pass


class WorkflowError(DDPrimerError):
    """
    Error in workflow execution.
    
    Raised when there are issues with overall pipeline
    workflow execution, including mode selection problems,
    workflow coordination failures, or step sequencing errors.
    
    Example:
        >>> raise WorkflowError("Invalid workflow mode specified")
    """
    pass


class ExternalToolError(DDPrimerError):
    """
    Error related to external tools like Primer3, Vienna, BLAST, etc.
    
    This exception provides detailed information about external tool
    failures, including command details, return codes, and output
    for comprehensive debugging and error reporting.
    
    Attributes:
        tool_name: Name of the external tool that failed
        command: Command that was executed
        return_code: Return code from the command
        stdout: Standard output from the command
        stderr: Standard error from the command
        
    Example:
        >>> raise ExternalToolError(
        ...     "ViennaRNA calculation failed",
        ...     tool_name="RNAfold",
        ...     command="RNAfold -T 37",
        ...     return_code=1,
        ...     stderr="Invalid temperature parameter"
        ... )
    """
    
    def __init__(self, message, tool_name=None, command=None, return_code=None, stdout=None, stderr=None):
        """
        Initialize with extended information about the external tool error.
        
        Creates an ExternalToolError with comprehensive debugging information
        about the failed external tool execution, including command details
        and output for troubleshooting.
        
        Args:
            message: Error message describing the failure
            tool_name: Name of the external tool (e.g., "primer3", "blastn")
            command: Command that was executed (for debugging)
            return_code: Return code from the command execution
            stdout: Standard output from the command
            stderr: Standard error from the command
            
        Example:
            >>> error = ExternalToolError(
            ...     "Tool execution failed",
            ...     tool_name="primer3",
            ...     return_code=1
            ... )
            >>> print(error.tool_name)  # "primer3"
        """
        self.tool_name = tool_name
        self.command = command
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr
        
        # Build detailed error message
        detailed_message = message
        if tool_name:
            detailed_message = f"{tool_name} error: {message}"
        if return_code is not None:
            detailed_message += f" (return code: {return_code})"
            
        super().__init__(detailed_message)
    
    def get_debug_info(self):
        """
        Get comprehensive debug information about the external tool error.
        
        Returns a dictionary containing all available debugging information
        about the failed external tool execution for logging and debugging.
        
        Returns:
            Dictionary with debug information including command, outputs, etc.
            
        Example:
            >>> error = ExternalToolError("Failed", tool_name="blast")
            >>> debug_info = error.get_debug_info()
            >>> print(debug_info['tool_name'])  # "blast"
        """
        return {
            'tool_name': self.tool_name,
            'command': self.command,
            'return_code': self.return_code,
            'stdout': self.stdout,
            'stderr': self.stderr,
            'message': str(self)
        }