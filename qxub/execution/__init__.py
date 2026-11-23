"""
Execution package for qxub.

This package provides job execution functionality including:
- The one submission pipeline, and the programmatic API (submit module)
- Execution contexts and the ``qxub exec`` path (context module)
- Variable expansion and context validation (core module)
- Execution mode detection for local vs remote (mode module)

The execution system supports:
- Conda environments
- Module loading
- Singularity containers
- Default (bare) execution
- Local and remote execution modes

Architecture:
- Single-threaded design for simplicity and reliability
- One submission pipeline (submit.dispatch) through the platform's scheduler backend
- Context validation and variable expansion
- Resource tracking and history logging integration
- Platform-aware execution routing (local vs remote)
"""

# Re-export QsubError so callers can catch it
from ..core.scheduler import QsubError

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
    expand_submission_variables,
    expand_variables_preserving_quotes,
    validate_execution_context,
)

# Execution mode detection
from .mode import (
    ExecutionMode,
    get_execution_mode,
    get_remote_config,
    validate_remote_config,
)

# Programmatic submission API
from .submit import SubmitResult, dispatch, submit_job

# Public API exports
__all__ = [
    # Core execution functions
    "validate_execution_context",
    "expand_submission_variables",
    "expand_variables_preserving_quotes",
    # Unified execution context
    "ExecutionContext",
    "execute_unified",
    "create_conda_context",
    "create_module_context",
    "create_singularity_context",
    "create_default_context",
    # Execution mode detection
    "ExecutionMode",
    "get_execution_mode",
    "get_remote_config",
    "validate_remote_config",
    # Programmatic submission API
    "submit_job",
    "dispatch",
    "SubmitResult",
    "QsubError",
]


# Convenience functions for common patterns
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
