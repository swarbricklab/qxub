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
