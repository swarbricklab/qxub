"""
Remote execution package for qxub v3.3+.

This package provides platform-aware remote execution capabilities:
- PlatformRemoteExecutor: SSH-based executor using platform config
- Command serialization for remote execution
- Remote execution error handling

Main classes:
- PlatformRemoteExecutor: Execute qxub commands on remote platforms via SSH
- RemoteExecutionError: Base exception for remote execution errors
"""

# Command building utilities
from .command_builder import build_remote_command

# Platform-aware remote executor
from .platform_executor import PlatformRemoteExecutor, RemoteExecutionError

__all__ = [
    "PlatformRemoteExecutor",
    "RemoteExecutionError",
    "build_remote_command",
]

# Core remote configuration
from .config import RemoteConfig, RemoteConfigError, UnsupportedProtocolError

# Remote execution
from .core import ConnectionError, RemoteExecutionError
from .core import UnsupportedProtocolError as CoreUnsupportedProtocolError

# Remote executor factory and backends
from .executor import RemoteExecutorFactory

# Configuration loading
from .loader import (
    ConfigLoadError,
    get_remote_config,
    get_user_config_path,
    load_remote_configurations,
    load_user_config,
)

__all__ = [
    # Configuration classes
    "RemoteConfig",
    "RemoteConfigError",
    "UnsupportedProtocolError",
    # Core execution classes
    "RemoteExecutionError",
    "ConnectionError",
    "CoreUnsupportedProtocolError",
    # Executor factory
    "RemoteExecutorFactory",
    # Configuration loading
    "ConfigLoadError",
    "get_remote_config",
    "get_user_config_path",
    "load_remote_configurations",
    "load_user_config",
]
