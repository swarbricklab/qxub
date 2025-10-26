"""
Execution package for qxub.

This package provides job execution functionality including:
- Core job submission and monitoring (core module)
- Context-specific executors (executors module)
- Unified execution context handling (context module)

The execution system supports:
- Conda environments
- Module loading
- Singularity containers
- Default (bare) execution

Architecture:
- Single-threaded design for simplicity and reliability
- Unified job submission and monitoring via core.submit_and_monitor_job
- Context validation and variable expansion
- Resource tracking and history logging integration
"""

# Unified execution context system
from .context import (
    ExecutionContext,
    create_conda_context,
    create_default_context,
    create_module_context,
    create_singularity_context,
    execute_unified,
)

# Core execution functionality
from .core import (
    build_submission_variables,
    expand_submission_variables,
    expand_variables_preserving_quotes,
    submit_and_monitor_job,
    validate_execution_context,
)

# Context-specific executors
from .executors import (
    execute_conda,
    execute_default,
    execute_module,
    execute_singularity,
)

# Public API exports
__all__ = [
    # Core execution functions
    "submit_and_monitor_job",
    "validate_execution_context",
    "build_submission_variables",
    "expand_submission_variables",
    "expand_variables_preserving_quotes",
    # Context-specific executors
    "execute_conda",
    "execute_module",
    "execute_singularity",
    "execute_default",
    # Unified execution context
    "ExecutionContext",
    "execute_unified",
    "create_conda_context",
    "create_module_context",
    "create_singularity_context",
    "create_default_context",
]


# Convenience functions for common patterns
def get_executor_for_context(context_type: str):
    """Get the appropriate executor function for a context type."""
    executors = {
        "conda": execute_conda,
        "module": execute_module,
        "singularity": execute_singularity,
        "default": execute_default,
    }
    return executors.get(context_type, execute_default)


def create_context_for_type(context_type: str, context_value=None):
    """Create an ExecutionContext for the given type and value."""
    if context_type == "conda":
        return create_conda_context(context_value or "base")
    elif context_type == "module":
        modules = (
            context_value
            if isinstance(context_value, list)
            else [context_value] if context_value else []
        )
        return create_module_context(modules)
    elif context_type == "singularity":
        return create_singularity_context(context_value or "")
    else:
        return create_default_context()
