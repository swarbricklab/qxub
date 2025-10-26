"""
Remote execution configuration and management for qxub v2.2.

This module provides URL-based remote configuration with support for multiple
protocols (starting with SSH-only in v2.2).
"""

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from urllib.parse import ParseResult, urlparse

logger = logging.getLogger(__name__)


class RemoteConfigError(Exception):
    """Base exception for remote configuration errors."""

    pass


class UnsupportedProtocolError(RemoteConfigError):
    """Raised when an unsupported protocol is specified."""

    pass


@dataclass
class RemoteConfig:
    """Configuration for remote execution endpoints."""

    name: str
    url: str
    platform: Optional[str] = (
        None  # Platform name (optional - remote will auto-detect if not specified)
    )
    conda_env: Optional[str] = (
        None  # Conda environment (optional - remote will use default)
    )
    working_dir: Optional[str] = None  # Working directory (optional - uses current dir)
    config: Optional[str] = None  # Protocol-specific config file (e.g., ~/.ssh/config)
    force_tty: Optional[bool] = (
        None  # Force TTY allocation (None=auto-detect, True=force, False=disable)
    )

    # Parsed URL components (set in __post_init__)
    parsed_url: ParseResult = field(init=False)

    def __post_init__(self):
        """Validate configuration and set smart defaults."""
        self.parsed_url = urlparse(self.url)

        # Validate URL format
        if not self.parsed_url.scheme:
            raise RemoteConfigError(
                f"URL must include protocol (e.g., ssh://host): {self.url}"
            )

        if not self.parsed_url.hostname:
            raise RemoteConfigError(f"URL must include hostname: {self.url}")

        # Check supported protocols
        if self.parsed_url.scheme not in self._supported_protocols():
            supported = ", ".join(self._supported_protocols())
            raise UnsupportedProtocolError(
                f"Protocol '{self.parsed_url.scheme}' not supported. "
                f"Supported protocols: {supported}"
            )

        # Set smart defaults for config file
        if self.config is None:
            self.config = self._get_default_config_path()

        # Expand user paths
        if self.config and self.config.startswith("~"):
            self.config = str(Path(self.config).expanduser())

    @staticmethod
    def _supported_protocols() -> list[str]:
        """Get list of currently supported protocols."""
        return ["ssh"]

    def _get_default_config_path(self) -> str:
        """Get default config file path based on protocol."""
        defaults = {
            "ssh": "~/.ssh/config",
            "k8s": "~/.kube/config",  # Future
            "aws": "~/.aws/config",  # Future
        }
        default_path = defaults.get(self.parsed_url.scheme, "")
        if default_path:
            return str(Path(default_path).expanduser())
        return ""

    @property
    def protocol(self) -> str:
        """Get connection protocol."""
        return self.parsed_url.scheme

    @property
    def hostname(self) -> str:
        """Get hostname from URL."""
        return self.parsed_url.hostname or ""

    @property
    def port(self) -> Optional[int]:
        """Get port from URL."""
        return self.parsed_url.port

    @property
    def username(self) -> Optional[str]:
        """Get username from URL."""
        return self.parsed_url.username

    def substitute_environment_vars(self, value: str) -> str:
        """Substitute environment variables in configuration values."""
        if not value:
            return value

        # Substitute common environment variables
        substitutions = {
            "${USER}": os.environ.get("USER", ""),
            "${PROJECT}": os.environ.get("PROJECT", ""),
            "${HOME}": os.environ.get("HOME", ""),
        }

        result = value
        for var, replacement in substitutions.items():
            result = result.replace(var, replacement)

        return result

    def get_working_dir(self) -> str:
        """Get working directory with environment variable substitution."""
        if not self.working_dir:
            return ""
        return self.substitute_environment_vars(self.working_dir)

    def determine_remote_working_dir(
        self, local_cwd: Path, explicit_execdir: Optional[str] = None
    ) -> str:
        """
        Determine the remote working directory.

        Args:
            local_cwd: Current working directory on local machine
            explicit_execdir: Explicitly specified execution directory

        Returns:
            Remote working directory path
        """
        if explicit_execdir:
            return explicit_execdir

        # Use configured working directory if set
        if self.working_dir:
            return self.get_working_dir()

        # Smart default: preserve current directory structure in home
        return f"~/{local_cwd.name}" if local_cwd.name != "." else "~"

    def validate(self) -> list[str]:
        """
        Validate the remote configuration.

        Returns:
            List of validation errors (empty if valid)
        """
        errors = []

        # Check config file exists if specified
        if self.config and not Path(self.config).expanduser().exists():
            errors.append(f"Config file not found: {self.config}")

        return errors

    def __str__(self) -> str:
        """String representation for debugging."""
        return (
            f"RemoteConfig(name={self.name}, url={self.url}, protocol={self.protocol})"
        )
