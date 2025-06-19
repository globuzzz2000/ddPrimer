# ddPrimer Coding Standards

## What NOT to Change

- **Never modify existing `.info()` calls**
- **Don't add new `.info()` calls** 
- **Don't remove existing `.info()` calls**
- **Don't change existing progress reporting**

## Logging Standards

### Logger Initialization
- **Use module-level loggers only**: `logger = logging.getLogger(__name__)` at module top
- **Remove class-level loggers**: No `self.logger` or class attributes like `logger = logging.getLogger("ddPrimer.module_name")`
- **Remove hardcoded names**: Replace `"ddPrimer.module_name"` with `__name__`

### Error Logging Patterns
```python
# Standard error logging with context
logger.error(f"Operation failed during {operation_name}: {context}")
logger.debug(f"Error details: {str(e)}", exc_info=True)

# Before raising custom exceptions
logger.error(error_msg)
raise CustomException(error_msg) from e

# For warnings about recoverable issues
logger.warning(f"Issue detected but continuing: {details}")
```

### Debug Message Structure
```python
# Debug sections for major operations
logger.debug("=== SECTION NAME DEBUG ===")
# ... debug content ...
logger.debug("=== END SECTION NAME DEBUG ===")

# Consistent formatting with relevant data
logger.debug(f"Processing {count} items with param={value}")
logger.debug(f"Results: success={success_count}, failed={fail_count}")
```

### Reducing Debug Redundancy Between Modules
Avoid duplicating debug information across modules. The main pipeline module should focus on high-level workflow, while core modules provide detailed operation info:

```python
# MAIN PIPELINE - High-level workflow only
logger.debug("Starting primer design workflow")
logger.debug(f"CHECKPOINT 1: {len(restriction_fragments)} restriction fragments")

# CORE MODULES - Detailed operation info
logger.debug(f"=== RESTRICTION SITE CUTTING ===")
logger.debug(f"Using restriction site pattern: {Config.RESTRICTION_SITE}")
logger.debug(f"Generated {fragment_count} valid fragments from {seq_id}")

# AVOID - Redundant logging in both main and core modules
# Main: logger.debug("Processing sequences with VCF variants...")
# Core: logger.debug("Processing sequences with VCF variants...")  # REDUNDANT
```

### Performance Considerations
```python
# Guard expensive debug operations
if logger.isEnabledFor(logging.DEBUG):
    expensive_debug_info = generate_complex_debug_data()
    logger.debug(f"Complex analysis: {expensive_debug_info}")

# Avoid logging in tight loops unless critical
```

### Debug Message Cleanup
Remove redundant debug statements that merely announce function entry/exit:
```python
# REMOVE - These add no value
logger.debug("Starting function X")
logger.debug("Function X completed")
logger.debug("Entering method Y")
logger.debug("Exiting method Y")

# KEEP - These provide useful information
logger.debug(f"Processing {len(sequences)} sequences")
logger.debug(f"Found {variant_count} variants in chromosome {chr_name}")
logger.debug(f"Analysis results: {summary_stats}")
```

## Docstring Standards

### Function/Method Docstrings
```python
def process_data(sequences: List[str], threshold: float = 0.5) -> Dict[str, Any]:
    """
    Process DNA sequences and filter by quality threshold.
    
    Performs quality analysis on input sequences and filters out those
    below the specified threshold. Returns detailed analysis results.
    
    Args:
        sequences: List of DNA sequences to process
        threshold: Quality threshold for filtering (0.0-1.0)
        
    Returns:
        Dictionary containing processed sequences and analysis results
        
    Raises:
        SequenceProcessingError: If sequences contain invalid characters
        ValueError: If threshold is outside valid range [0.0, 1.0]
        FileError: If required reference files cannot be accessed
        
    Example:
        >>> sequences = ["ATCG", "GCTA"] 
        >>> results = process_data(sequences, threshold=0.8)
        >>> len(results['filtered_sequences'])
        2
    """
```

### Class Docstrings
```python
class SequenceProcessor:
    """
    Handles DNA sequence processing and quality analysis.
    
    This class provides methods for analyzing DNA sequences, applying
    quality filters, and generating summary statistics. It supports
    both single sequences and batch processing.
    
    Attributes:
        threshold: Default quality threshold for filtering
        processed_count: Number of sequences processed in current session
        
    Example:
        >>> processor = SequenceProcessor(threshold=0.8)
        >>> results = processor.process_batch(sequences)
        >>> processor.get_statistics()
    """
```

### Module Docstrings
```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence processing module for ddPrimer pipeline.

Contains functionality for:
1. DNA sequence quality analysis
2. Filtering and validation 
3. Batch processing operations
4. Statistical reporting

This module integrates with the broader ddPrimer pipeline to provide
robust sequence processing capabilities for primer design workflows.
"""
```

## Error Handling Standards

### Custom Exception Usage
```python
# Use specific exception types
raise SequenceProcessingError(f"Invalid sequence format in {filename}")
raise FileError(f"Cannot access required file: {filepath}")
raise ExternalToolError(f"BLAST execution failed", tool_name="blastn")

# Chain exceptions to preserve context
try:
    risky_operation()
except ValueError as e:
    error_msg = f"Data validation failed for {data_id}"
    logger.error(error_msg)
    raise SequenceProcessingError(error_msg) from e
```

### Error Context Requirements
- **Always include what was being attempted**: "Failed to process sequence X"
- **Include relevant identifiers**: sequence names, file paths, parameters
- **Provide actionable information**: what the user could check or fix
- **Log before raising**: ensure errors are captured even if exceptions are caught upstream

### Exception Documentation
```python
def critical_function(data: str) -> bool:
    """
    Process critical data with validation.
    
    Args:
        data: Input data string
        
    Returns:
        True if processing successful
        
    Raises:
        SequenceProcessingError: If data format is invalid or processing fails
        FileError: If required reference files cannot be accessed
        ExternalToolError: If external tools (BLAST, Primer3) fail
    """
```

## Code Comment Cleanup

### Remove Implementation History Comments
Remove comments that reference changes, fixes, or development history:

```python
# REMOVE these types of comments:
# Fixed issue with coordinate mapping
# Enhanced error handling for edge cases  
# New implementation using bcftools
# Improved performance in this section
# Temporary fix for issue #123
# TODO: Refactor this when time permits
# Updated to handle new file format
# Workaround for threading issue

# KEEP functional comments:
# Convert 1-based coordinates to 0-based for Python
# Handle edge case where sequence ends before primer
# Apply thermodynamic constraints from configuration
```

### Method and Variable Name Cleanup
Remove "fix", "new", "enhanced", etc. from names:

```python
# CHANGE from descriptive of change to descriptive of purpose:
def process_sequences_fixed() -> None:  # CHANGE TO: def process_sequences()
def validate_input_enhanced() -> bool:  # CHANGE TO: def validate_input()
def new_coordinate_mapper() -> Dict:    # CHANGE TO: def coordinate_mapper()

# Variables:
enhanced_results = {}     # CHANGE TO: results = {}
fixed_coordinates = []    # CHANGE TO: coordinates = []
new_algorithm_params = {} # CHANGE TO: algorithm_params = {}
```

### Comment Content Standards
Focus on explaining WHY and WHAT, not WHEN or HOW it changed:

```python
# GOOD - Explains purpose and context
# Use reverse complement for probe design to optimize C>G ratio
# BLAST word size of 7 recommended for short primer sequences
# Minimum segment length must accommodate primer product size range

# BAD - Documents change history instead of purpose  
# Fixed the reverse complement calculation
# Updated BLAST parameters after testing
# New validation logic added here
```

## Migration Checklist

### Logging Updates
- [ ] Replace hardcoded logger names with `logger = logging.getLogger(__name__)`
- [ ] Remove class-level logger attributes
- [ ] Add context to all error messages
- [ ] Add `exc_info=True` to debug logging for exceptions
- [ ] Standardize debug section formatting
- [ ] Guard expensive debug operations
- [ ] Remove redundant debug entry/exit statements

### Docstring Updates  
- [ ] Add missing function/method docstrings
- [ ] Include all parameters in Args section
- [ ] Document return values with types
- [ ] List all possible exceptions in Raises section
- [ ] Add usage examples where helpful
- [ ] Update class docstrings with attributes and examples

### Error Handling Updates
- [ ] Use specific custom exception types
- [ ] Chain exceptions with `from e` 
- [ ] Log errors before raising
- [ ] Include relevant context in error messages
- [ ] Document all exceptions in docstrings

### Code Comment Cleanup
- [ ] Remove implementation history references
- [ ] Remove "fix", "new", "enhanced" from method names
- [ ] Remove "TODO" and "FIXME" comments unless actionable
- [ ] Focus comments on explaining purpose, not change history
- [ ] Remove redundant comments that just restate the code

## Example Conversion

### Before
```python
class BlastProcessor:
    logger = logging.getLogger("ddPrimer.blast_processor")
    
    def blast_short_seq_enhanced(self, seq):
        # New improved implementation
        self.logger.debug("Starting BLAST operation")  # Redundant
        if not seq:
            self.logger.debug("Entering validation branch")  # Redundant
            return None, None
        try:
            # Fixed coordinate handling issue
            # process
            pass
        except Exception as e:
            self.logger.error("BLAST Error")  # No context
            return None, None
        finally:
            self.logger.debug("BLAST operation completed")  # Redundant
```

### After
```python
logger = logging.getLogger(__name__)

class BlastProcessor:
    """Handles BLAST operations for primer specificity checking."""
    
    def blast_short_seq(self, seq: str) -> tuple[Optional[float], Optional[float]]:
        """
        Run BLASTn for short sequences and return the two best e-values.
        
        Args:
            seq: DNA sequence to BLAST against database
            
        Returns:
            Tuple of (best_evalue, second_best_evalue), or (None, None) if failed
            
        Raises:
            SequenceProcessingError: If BLAST execution fails
        """
        if not seq or not isinstance(seq, str) or not seq.strip():
            logger.debug("Empty or invalid sequence provided to BLAST")
            return None, None
            
        try:
            # process
            pass
        except Exception as e:
            error_msg = f"BLAST execution failed for sequence {seq[:20]}..."
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise SequenceProcessingError(error_msg) from e
```