# ddPrimer Coding Standards

## What NOT to Change

- **Never modify existing `.info()` calls**
- **Don't add new `.info()` calls** 
- **Don't remove existing `.info()` calls**
- **Don't change existing progress reporting**

## Architecture: Workflow Wrapper Integration

### Workflow Wrapper Structure
Each core and utils module must include a workflow wrapper section at the top of the main class:

```python
class ModuleProcessor:
    """Class docstring."""
    
    #############################################################################
    #                           Workflow Wrappers
    #############################################################################
    
    @classmethod  # or @staticmethod as appropriate
    def method_name_workflow(cls, args) -> ReturnType:
        """
        Brief description for workflow integration.
        
        Longer description explaining what this wrapper does in the context
        of the overall pipeline workflow.
        
        Args:
            args: Method arguments
            
        Returns:
            Return description
            
        Raises:
            SpecificModuleError: If specific operation fails (e.g., SequenceProcessingError)
            PipelineError: If workflow coordination fails
        """
        logger.debug("=== WORKFLOW: OPERATION NAME ===")
        logger.debug(f"Brief setup description with key parameters")
        
        try:
            # Delegate to core implementation methods
            result = cls.core_implementation_method(args)
            
            logger.debug("Operation complete with summary statistics")
            logger.debug("=== END WORKFLOW: OPERATION NAME ===")
            
            return result
            
        except Exception as e:
            error_msg = f"Error in operation workflow: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END WORKFLOW: OPERATION NAME ===")
            
            # Re-raise specific module errors, wrap others as PipelineError
            if isinstance(e, (SequenceProcessingError, PrimerDesignError, FileError, ExternalToolError)):
                raise  # Let specific errors bubble up
            else:
                raise PipelineError(error_msg) from e
    
    #############################################################################
    
    def core_implementation_method(self, args):
        # Existing detailed implementation...
```

### Workflow Wrapper Guidelines

#### Naming Convention
- **Workflow wrappers must end with `_workflow`**
- Use descriptive names that indicate the workflow step: `process_sequences_with_vcf_batch`, `filter_primers_workflow`
- Use `@classmethod` for operations that don't need instance state
- Use `@staticmethod` for pure functional operations
- Use instance methods only when module configuration/state is required

#### Placement Requirements
- **Must be placed at the very top of the class**, immediately after class docstring
- **Must be wrapped in the comment block structure** with proper formatting
- **Must be separated from regular methods** by the closing comment block

#### Content Requirements
- **Debug boundaries**: Start and end with `=== WORKFLOW: OPERATION NAME ===`
- **Delegation focus**: Primary purpose is delegating to core implementation methods
- **Minimal logic**: Should contain orchestration logic only, not implementation details
- **Smart error handling**: Re-raise specific module errors, wrap generic errors as `PipelineError`
- **Progress logging**: Include brief setup and completion messages with key metrics

#### Error Handling Strategy
```python
# PREFERRED - Let specific errors bubble up for proper handling
except Exception as e:
    error_msg = f"Error in operation workflow: {str(e)}"
    logger.error(error_msg)
    logger.debug(f"Error details: {str(e)}", exc_info=True)
    logger.debug("=== END WORKFLOW: OPERATION NAME ===")
    
    # Re-raise specific module errors, wrap others as PipelineError
    if isinstance(e, (SequenceProcessingError, PrimerDesignError, FileError, ExternalToolError)):
        raise  # Let specific errors bubble up
    else:
        raise PipelineError(error_msg) from e
```

#### What Workflow Wrappers Should Do
- **Provide workflow-level entry points** for orchestration
- **Handle workflow-level logging** with debug boundaries
- **Delegate to core implementation methods** without duplicating logic
- **Collect and report summary statistics** (counts, success rates)
- **Maintain error context** for debugging while preserving specific error types

#### What Workflow Wrappers Should NOT Do
- **No complex business logic**: Complex processing should remain in core methods
- **No .info() logging**: Use only .debug() and .error() logging levels
- **No progress bars**: Leave progress reporting to core methods
- **No detailed statistics**: Summary only, details in core methods
- **No parameter validation**: Core methods should handle validation
- **No duplicate error handling**: Don't catch and re-raise the same error type

## Logging Standards

### Logger Initialization
- **Use module-level loggers only**: `logger = logging.getLogger(__name__)` at module top
- **Remove class-level loggers**: No `self.logger` or class attributes like `logger = logging.getLogger("ddPrimer.module_name")`
- **Remove hardcoded names**: Replace `"ddPrimer.module_name"` with `__name__`

### Logging Hierarchy and Responsibility

#### Main Pipeline (main.py) - User-Facing Workflow Logging
```python
# WORKFLOW PHASES - Major user-visible steps
logger.info("Processing sequences with VCF variants...")
logger.info("Calculating thermodynamic properties with ViennaRNA...")
logger.info("Pipeline completed successfully")

# WORKFLOW CHECKPOINTS - High-level progress only
logger.debug("MAIN: Pipeline mode determined as 'direct'")
logger.debug("MAIN: Output directory prepared: /path/to/output")

# USER-FRIENDLY ERRORS - Convert technical errors for users
logger.error(f"VCF sequence processing failed: {str(e)}")

# DELEGATION - What we're asking modules to do
logger.debug("MAIN: Delegating to workflow orchestrator")
```

#### Workflow Orchestrator (workflow.py) - Coordination Logging
```python
# WORKFLOW ORCHESTRATION - High-level coordination
logger.debug("WORKFLOW: Delegating BLAST database setup")
logger.debug("WORKFLOW: Delegating primer design")

# COORDINATION ERRORS - Workflow-level error handling
logger.error(f"Primer design failed: {str(e)}")

# WORKFLOW BOUNDARIES - Clear stage transitions
logger.debug("=== WORKFLOW ORCHESTRATOR: PIPELINE EXECUTION ===")
logger.debug("=== END WORKFLOW ORCHESTRATOR: PIPELINE EXECUTION ===")
```

#### Workflow Wrappers - Operation Boundaries
```python
# OPERATION BOUNDARIES - Clear start/end markers
logger.debug("=== WORKFLOW: RESTRICTION SITE PROCESSING ===")
logger.debug("=== END WORKFLOW: RESTRICTION SITE PROCESSING ===")

# OPERATION SETUP