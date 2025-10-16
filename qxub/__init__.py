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
