"""
Execution package for qxub.

This package provides job execution functionality including:
- The one submission pipeline, and the programmatic API (submit module)
- Execution contexts and the ``qxub exec`` path (context module)
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
- Context validation
- Resource tracking and history logging integration
- Platform-aware execution routing (local vs remote)
"""

# Re-export QsubError so callers can catch it
from ..core.scheduler import QsubError

# Unified execution context system
from .context import ExecutionContext, execute_unified

# Execution mode detection
from .mode import ExecutionMode, get_execution_mode

# Programmatic submission API
from .submit import SubmitResult, dispatch, submit_job

# Public API exports
__all__ = [
    # Unified execution context
    "ExecutionContext",
    "execute_unified",
    # Execution mode detection
    "ExecutionMode",
    "get_execution_mode",
    # Programmatic submission API
    "submit_job",
    "dispatch",
    "SubmitResult",
    "QsubError",
]
