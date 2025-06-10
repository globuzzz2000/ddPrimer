# ddPrimer Coding Standards

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

### Performance Considerations
```python
# Guard expensive debug operations
if logger.isEnabledFor(logging.DEBUG):
    expensive_debug_info = generate_complex_debug_data()
    logger.debug(f"Complex analysis: {expensive_debug_info}")

# Avoid logging in tight loops unless critical
```

### What NOT to Change
- **Never modify existing `.info()` calls**
- **Don't add new `.info()` calls** 
- **Don't remove existing `.info()` calls**
- **Don't change existing progress reporting**

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

## Migration Checklist

### Logging Updates
- [ ] Replace hardcoded logger names with `logger = logging.getLogger(__name__)`
- [ ] Remove class-level logger attributes
- [ ] Add context to all error messages
- [ ] Add `exc_info=True` to debug logging for exceptions
- [ ] Standardize debug section formatting
- [ ] Guard expensive debug operations

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

## Example Conversion

### Before
```python
class BlastProcessor:
    logger = logging.getLogger("ddPrimer.blast_processor")
    
    def blast_short_seq(self, seq):
        if not seq:
            return None, None
        try:
            # process
            pass
        except Exception as e:
            self.logger.error("BLAST Error")
            return None, None
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