"""
Remote execution backend for qxub v2.2

This module handles URL-based remote execution with support for multiple protocols.
Currently implements SSH-only execution with future extensibility for AWS, Kubernetes, etc.
"""

import logging
import subprocess
from abc import ABC, abstractmethod
from typing import Any, List, Optional

logger = logging.getLogger(__name__)


class RemoteExecutionError(Exception):
    """Base class for remote execution errors."""

    pass


class UnsupportedProtocolError(RemoteExecutionError):
    """Raised when protocol is not supported."""

    pass


class ConnectionError(RemoteExecutionError):
    """Raised when connection fails."""

    pass


class ConfigError(RemoteExecutionError):
    """Raised when configuration is invalid."""

    pass
