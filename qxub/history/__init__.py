"""
qxub.history - Command History and Recipe Management Package

This package provides command history tracking, recipe management, and execution
record persistence for qxub jobs. It supports both computational recipes and
detailed execution logging.

Package Structure:
    base.py       - Core history data structures and constants
    manager.py    - History management operations and database interface

Public API:
    # Core History Management
    HistoryManager                        # Main history management interface

    # History Types
    RECIPE_HISTORY                        # Computational recipe tracking
    EXECUTION_HISTORY                     # Detailed execution records

    # History Operations
    record_recipe(cmd, resources)         # Record computational recipe
    record_execution(job_details)         # Record execution details
    get_recent_recipes(limit)             # Retrieve recent recipes
    get_execution_history(filters)        # Query execution history

Example Usage:
    from qxub.history import HistoryManager, RECIPE_HISTORY

    # Create history manager
    history = HistoryManager()

    # Record a computational recipe
    history.record_recipe(
        command="python train.py --epochs 100",
        resources={"mem": "8GB", "walltime": "2:00:00", "ncpus": 4}
    )

    # Get recent recipes
    recent = history.get_recent_recipes(limit=10)
"""

# Import all public APIs for convenient access
from .base import CommandHistoryLogger, history_logger
from .manager import HistoryManager

# Define what gets imported with "from qxub.history import *"
__all__ = [
    # Classes
    "HistoryManager",
    "CommandHistoryLogger",
    # Global instances
    "history_logger",
]
