"""
Remote execution backend for qxub v2.2

This module handles URL-based remote execution with support for multiple protocols.
Currently implements SSH-only execution with future extensibility for AWS, Kubernetes, etc.
"""

import logging
import os
import subprocess
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urlparse

import yaml

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


@dataclass
class RemoteConfig:
    """Configuration for remote execution endpoints."""

    url: str  # protocol://host:port format
    qxub_env: str
    platform_file: str
    config: Optional[str] = None  # Protocol-specific config file
    project_root_dir: Optional[str] = None
    _parsed_url: Optional[object] = field(default=None, init=False)

    def __post_init__(self):
        """Validate and set smart defaults."""
        self._parsed_url = urlparse(self.url)

        if not self._parsed_url.scheme:
            raise ConfigError(f"URL must include protocol: {self.url}")

        if self._parsed_url.scheme not in ["ssh"]:
            raise UnsupportedProtocolError(
                f"Unsupported protocol: {self._parsed_url.scheme}"
            )

        # Set smart defaults for config file
        if self.config is None:
            self.config = self._get_default_config_path()

    def _get_default_config_path(self) -> str:
        """Get default config file path based on protocol."""
        defaults = {
            "ssh": "~/.ssh/config",
            "k8s": "~/.kube/config",  # Future
            "aws": "~/.aws/config",  # Future
        }
        return defaults.get(self._parsed_url.scheme, "")

    @property
    def protocol(self) -> str:
        """Get connection protocol."""
        return self._parsed_url.scheme

    @property
    def hostname(self) -> str:
        """Get hostname from URL."""
        return self._parsed_url.hostname or ""

    @property
    def port(self) -> Optional[int]:
        """Get port from URL."""
        return self._parsed_url.port


def load_remote_configurations() -> Dict[str, RemoteConfig]:
    """Load remote configurations from user config file."""

    config_path = Path.home() / ".config" / "qxub" / "config.yaml"

    if not config_path.exists():
        return {}

    try:
        with open(config_path) as f:
            config_data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigError(f"Invalid YAML in config file: {e}")

    if not config_data:
        return {}

    remotes = {}
    for name, remote_data in config_data.get("remotes", {}).items():
        try:
            # Substitute environment variables
            for key, value in remote_data.items():
                if isinstance(value, str):
                    remote_data[key] = os.path.expandvars(value)

            remotes[name] = RemoteConfig(
                url=remote_data["url"],
                qxub_env=remote_data["qxub_env"],
                platform_file=remote_data["platform_file"],
                config=remote_data.get("config"),
                project_root_dir=remote_data.get("project_root_dir"),
            )
        except KeyError as e:
            raise ConfigError(f"Missing required field in remote '{name}': {e}")
        except Exception as e:
            raise ConfigError(f"Error loading remote '{name}': {e}")

    return remotes
