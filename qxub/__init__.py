"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.

v2.2 adds remote execution capabilities for seamless job submission to remote HPC systems.
v2.3 adds parallel job execution with --terse output and 'qxub monitor' for pipeline workflows.
"""

# __init__.py

__version__ = "2.3.6"

# Import main CLI
from .cli import qxub

# Import remote execution components (v2.2)
try:
    from .remote_config import RemoteConfig
    from .remote_config_loader import get_remote_config, load_remote_configurations
    from .remote_executor import RemoteExecutorFactory

    __all__ = [
        "qxub",
        "RemoteConfig",
        "load_remote_configurations",
        "get_remote_config",
        "RemoteExecutorFactory",
    ]
except ImportError:
    # Remote execution dependencies not available
    __all__ = ["qxub"]

from . import config, scheduler
from .cli import qxub

# Backwards compatibility for resources and history packages (Phase 1 Migration)
# Backwards compatibility for config package (Phase 2 Migration)
# These imports allow existing code to continue working while we migrate
try:
    # Import key resource utilities for backwards compatibility
    # Import key config components for backwards compatibility
    from .config import (
        ConfigManager,
        ShortcutManager,
        config_manager,
        get_config,
        set_config,
    )

    # Import key history components for backwards compatibility
    from .history import HistoryManager, history_logger
    from .resources import (
        ResourceMapper,
        ResourceTracker,
        format_memory_size,
        format_walltime,
        parse_memory_size,
        parse_walltime,
    )

    # Make them available at package level for existing imports like:
    # from qxub import parse_memory_size, config_manager, history_logger
    __all__.extend(
        [
            "parse_memory_size",
            "parse_walltime",
            "format_memory_size",
            "format_walltime",
            "ResourceTracker",
            "ResourceMapper",
            "HistoryManager",
            "history_logger",
            "ConfigManager",
            "config_manager",
            "ShortcutManager",
            "get_config",
            "set_config",
        ]
    )
except ImportError:
    # Package import failed - migration may be in progress
    pass
