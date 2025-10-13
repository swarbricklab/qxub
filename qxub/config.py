"""
Configuration management for qxub including platform definitions and user preferences.

This module handles:
- Platform configuration loading from system and user locations
- User preference management and validation
- Configuration file discovery and validation
- Default settings and environment-specific overrides
"""

import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

logger = logging.getLogger(__name__)


class QxubConfig:
    """Central configuration manager for qxub."""

    def __init__(self):
        self.platform_search_paths = [
            Path("/etc/qxub/platforms"),
            Path.home() / ".config" / "qxub" / "platforms",
            Path.home() / ".qxub" / "platforms",
        ]
        self.config_search_paths = [
            Path("/etc/qxub"),
            Path.home() / ".config" / "qxub",
            Path.home() / ".qxub",
        ]
        self._user_config = {}
        self._system_config = {}
        self._load_config()

    def _load_config(self):
        """Load configuration from all available sources."""
        # Load system config
        system_config_file = Path("/etc/qxub/config.yaml")
        if system_config_file.exists():
            try:
                with open(system_config_file) as f:
                    self._system_config = yaml.safe_load(f) or {}
                logger.debug(f"Loaded system config from {system_config_file}")
            except Exception as e:
                logger.warning(f"Failed to load system config: {e}")

        # Load user config
        user_config_paths = [
            Path.home() / ".config" / "qxub" / "config.yaml",
            Path.home() / ".qxub" / "config.yaml",
        ]

        for config_path in user_config_paths:
            if config_path.exists():
                try:
                    with open(config_path) as f:
                        self._user_config = yaml.safe_load(f) or {}
                    logger.debug(f"Loaded user config from {config_path}")
                    break
                except Exception as e:
                    logger.warning(
                        f"Failed to load user config from {config_path}: {e}"
                    )

    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value with user config taking precedence."""
        # Check user config first
        if key in self._user_config:
            return self._user_config[key]

        # Fall back to system config
        if key in self._system_config:
            return self._system_config[key]

        return default

    def get_platform_search_paths(self) -> List[Path]:
        """Get platform search paths from config, environment, or defaults."""
        import os

        # Check environment variable first
        env_paths = os.getenv("QXUB_PLATFORM_PATHS")
        if env_paths:
            # Support both single path and colon-separated paths
            if ":" in env_paths:
                return [Path(p.strip()) for p in env_paths.split(":")]
            else:
                return [Path(env_paths)]

        # Then check config
        configured_paths = self.get("platform_search_paths", [])
        if configured_paths:
            return [Path(p) for p in configured_paths]

        # Finally use defaults
        return self.platform_search_paths

    def get_default_platform(self) -> Optional[str]:
        """Get the default platform name from config."""
        return self.get("default_platform")

    def get_platform_preferences(self) -> Dict[str, Any]:
        """Get platform-specific preferences."""
        return self.get("platform_preferences", {})

    def get_queue_preferences(self) -> Dict[str, Any]:
        """Get queue selection preferences."""
        return self.get(
            "queue_preferences",
            {
                "optimization": "balanced",  # cost, speed, balanced
                "adjustment_policy": "suggest",  # auto, suggest, user, error
                "auto_select": True,
            },
        )

    def get_default_resources(self) -> Dict[str, Any]:
        """Get default resource specifications."""
        return self.get(
            "default_resources", {"cpus": 1, "memory": "4GB", "walltime": "1:00:00"}
        )

    def set_user_preference(self, key: str, value: Any):
        """Set a user preference (runtime only, not persisted)."""
        self._user_config[key] = value

    def create_user_config_dir(self) -> Path:
        """Create user config directory if it doesn't exist."""
        config_dir = Path.home() / ".config" / "qxub"
        config_dir.mkdir(parents=True, exist_ok=True)
        return config_dir

    def save_user_config(self):
        """Save current user configuration to file."""
        config_dir = self.create_user_config_dir()
        config_file = config_dir / "config.yaml"

        try:
            with open(config_file, "w") as f:
                yaml.dump(self._user_config, f, default_flow_style=False)
            logger.info(f"Saved user config to {config_file}")
        except Exception as e:
            logger.error(f"Failed to save user config: {e}")
            raise

    def init_user_config(self) -> Path:
        """Initialize user config directory with examples."""
        config_dir = self.create_user_config_dir()
        config_file = config_dir / "config.yaml"
        platforms_dir = config_dir / "platforms"

        # Create platforms directory
        platforms_dir.mkdir(exist_ok=True)

        # Create example config if it doesn't exist
        if not config_file.exists():
            example_config = {
                "default_platform": None,
                "queue_preferences": {
                    "optimization": "balanced",
                    "adjustment_policy": "suggest",
                    "auto_select": True,
                },
                "default_resources": {
                    "cpus": 1,
                    "memory": "4GB",
                    "walltime": "1:00:00",
                },
                "platform_preferences": {},
                "verbosity": 0,
            }

            with open(config_file, "w") as f:
                yaml.dump(example_config, f, default_flow_style=False)

            logger.info(f"Created example config at {config_file}")

        return config_dir


# Global config instance
_config = None


def get_config() -> QxubConfig:
    """Get global config instance."""
    global _config
    if _config is None:
        _config = QxubConfig()
    return _config


def setup_logging(verbosity: int = None):
    """
    Configures the logging level based on the verbosity provided by the user.

    Args:
        verbosity (int): The number of '-v' flags used, or from config.
                       - 0: ERROR level (default)
                       - 1: WARNING level
                       - 2: INFO level
                       - 3 or more: DEBUG level

    This function adjusts the logging output to provide more detailed information
    as verbosity increases, allowing users to control the granularity of log messages.
    """
    if verbosity is None:
        verbosity = get_config().get("verbosity", 0)

    # Get the root logger and configure it directly
    root_logger = logging.getLogger()

    # Remove existing handlers to avoid conflicts
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Set the level and format based on verbosity
    if verbosity == 1:
        level = logging.WARNING
        format_str = "%(levelname)s: %(message)s"
    elif verbosity == 2:
        level = logging.INFO
        format_str = "%(levelname)s: %(message)s"
    elif verbosity >= 3:
        level = logging.DEBUG
        format_str = "%(levelname)s:%(name)s: %(message)s"
    else:
        level = logging.ERROR
        format_str = "%(levelname)s: %(message)s"

    # Configure the root logger
    root_logger.setLevel(level)
    handler = logging.StreamHandler()
    handler.setLevel(level)
    formatter = logging.Formatter(format_str)
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)


def validate_config() -> List[str]:
    """Validate current configuration and return any issues."""
    issues = []
    config = get_config()

    # Check platform search paths exist
    for path in config.get_platform_search_paths():
        if not path.exists():
            issues.append(f"Platform search path does not exist: {path}")

    # Check default platform is available
    default_platform = config.get_default_platform()
    if default_platform:
        from .platform import get_platform

        platform = get_platform(default_platform)
        if not platform:
            issues.append(f"Default platform '{default_platform}' not found")

    # Validate queue preferences
    queue_prefs = config.get_queue_preferences()
    valid_optimizations = ["cost", "speed", "balanced"]
    valid_policies = ["auto", "suggest", "user", "error"]

    if queue_prefs.get("optimization") not in valid_optimizations:
        issues.append(
            f"Invalid optimization preference: {queue_prefs.get('optimization')}"
        )

    if queue_prefs.get("adjustment_policy") not in valid_policies:
        issues.append(
            f"Invalid adjustment policy: {queue_prefs.get('adjustment_policy')}"
        )

    return issues


def get_effective_config() -> Dict[str, Any]:
    """Get complete effective configuration with all sources merged."""
    config = get_config()

    return {
        "platform_search_paths": [str(p) for p in config.get_platform_search_paths()],
        "default_platform": config.get_default_platform(),
        "platform_preferences": config.get_platform_preferences(),
        "queue_preferences": config.get_queue_preferences(),
        "default_resources": config.get_default_resources(),
        "verbosity": config.get("verbosity", 0),
    }
