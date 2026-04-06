"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.

v2.2 adds remote execution capabilities for seamless job submission to remote HPC systems.
v2.3 adds parallel job execution with --terse output and 'qxub monitor' for pipeline workflows.
"""

# __init__.py

# Suppress system Python requests/urllib3 version warnings
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning, module="requests")
try:
    from requests.exceptions import RequestsDependencyWarning

    warnings.filterwarnings("ignore", category=RequestsDependencyWarning)
except ImportError:
    pass

__version__ = "3.5.2"

# Library-standard NullHandler: prevents "No handler found" warnings and
# avoids triggering basicConfig() when qxub is used as a library.
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())


# Import main CLI and heavy subpackages lazily to avoid pulling in the
# entire Click/omegaconf dependency chain when qxub is imported as a library.
def __getattr__(name):
    # CLI entry point
    if name == "qxub":
        from . import cli as cli_module

        globals()["qxub"] = cli_module.qxub
        return cli_module.qxub

    # Remote execution
    if name == "PlatformRemoteExecutor":
        from .remote import PlatformRemoteExecutor

        globals()["PlatformRemoteExecutor"] = PlatformRemoteExecutor
        return PlatformRemoteExecutor
    if name == "RemoteExecutionError":
        from .remote import RemoteExecutionError

        globals()["RemoteExecutionError"] = RemoteExecutionError
        return RemoteExecutionError
    if name == "build_remote_command":
        from .remote import build_remote_command

        globals()["build_remote_command"] = build_remote_command
        return build_remote_command

    # Backwards-compatibility imports (config, history, resources)
    if name in _COMPAT_IMPORTS:
        module_path, attr = _COMPAT_IMPORTS[name]
        import importlib

        mod = importlib.import_module(module_path, __name__)
        value = getattr(mod, attr)
        globals()[name] = value
        return value

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "qxub",
    "PlatformRemoteExecutor",
    "RemoteExecutionError",
    "build_remote_command",
    # Backwards compatibility
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

# Backwards-compatibility lazy imports: heavy subpackages are only loaded
# when the caller actually accesses one of these names.
_COMPAT_IMPORTS = {
    # config
    "config_manager": (".config", "config_manager"),
    "ConfigManager": (".config", "ConfigManager"),
    "ShortcutManager": (".config", "ShortcutManager"),
    "get_config": (".config", "get_config"),
    "set_config": (".config", "set_config"),
    # history
    "HistoryManager": (".history", "HistoryManager"),
    "history_logger": (".history", "history_logger"),
    # resources
    "ResourceTracker": (".resources", "ResourceTracker"),
    "ResourceMapper": (".resources", "ResourceMapper"),
    "parse_memory_size": (".resources", "parse_memory_size"),
    "parse_walltime": (".resources", "parse_walltime"),
    "format_memory_size": (".resources", "format_memory_size"),
    "format_walltime": (".resources", "format_walltime"),
}
