"""
Remote execution package for qxub.

This package provides remote execution capabilities for qxub v2.2+ including:
- URL-based remote configuration
- SSH-based remote execution (with future protocol extensibility)
- Remote configuration loading and validation
- Remote execution backends and error handling

Main classes and functions:
- RemoteConfig: Remote execution configuration
- RemoteExecutorFactory: Factory for creating remote executors
- load_remote_configurations: Load user remote configurations
- get_remote_config: Get specific remote configuration

Supported protocols:
- SSH: ssh://user@hostname:port
- Future: AWS, Kubernetes, etc.
"""

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
