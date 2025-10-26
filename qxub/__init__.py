"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.

v2.2 adds remote execution capabilities for seamless job submission to remote HPC systems.
v2.3 adds parallel job execution with --terse output and 'qxub monitor' for pipeline workflows.
"""

# __init__.py

__version__ = "3.2.0"

# Import main CLI
from . import cli as cli_module

qxub = cli_module.qxub

# Import remote execution components (v2.2)
try:
    from .remote import (  # noqa: F401
        RemoteConfig,
        RemoteExecutorFactory,
        get_remote_config,
        load_remote_configurations,
    )

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

# Backwards compatibility imports (must come after conditional remote imports)
from . import config  # noqa: F401,E402
from .core import scheduler  # noqa: F401,E402

# CLI is imported at the top of the file

# Backwards compatibility for resources and history packages (Phase 1 Migration)
# Backwards compatibility for config package (Phase 2 Migration)
# These imports allow existing code to continue working while we migrate
try:
    # Import key resource utilities for backwards compatibility
    # Import key config components for backwards compatibility
    from .config import (  # noqa: F401
        ConfigManager,
        ShortcutManager,
        config_manager,
        get_config,
        set_config,
    )

    # Import key history components for backwards compatibility
    from .history import HistoryManager, history_logger  # noqa: F401
    from .resources import (  # noqa: F401
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
