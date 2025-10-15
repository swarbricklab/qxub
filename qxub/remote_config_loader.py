"""
Configuration loading and management for qxub remote execution.

This module handles loading remote configurations from user config files
and provides validation and error handling.
"""

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

from .remote_config import RemoteConfig, RemoteConfigError

logger = logging.getLogger(__name__)


class ConfigLoadError(Exception):
    """Raised when configuration loading fails."""

    pass


def get_user_config_path() -> Path:
    """Get the path to the user's qxub configuration file."""
    return Path.home() / ".config" / "qxub" / "config.yaml"


def load_user_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """
    Load the user's qxub configuration file.

    Args:
        config_file: Optional path to alternative config file

    Returns:
        Configuration dictionary

    Raises:
        ConfigLoadError: If configuration file cannot be loaded
    """
    if config_file:
        config_path = Path(config_file).expanduser().resolve()
    else:
        config_path = get_user_config_path()

    if not config_path.exists():
        return {}

    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        if config is None:
            return {}

        if not isinstance(config, dict):
            raise ConfigLoadError(
                f"Configuration file must contain a YAML dictionary: {config_path}"
            )

        return config

    except yaml.YAMLError as e:
        raise ConfigLoadError(f"Invalid YAML in configuration file {config_path}: {e}")
    except IOError as e:
        raise ConfigLoadError(f"Cannot read configuration file {config_path}: {e}")


def load_remote_configurations(
    config_file: Optional[str] = None,
) -> Dict[str, RemoteConfig]:
    """
    Load remote configurations from user config file.

    Args:
        config_file: Optional path to alternative config file

    Returns:
        Dictionary mapping remote names to RemoteConfig objects

    Raises:
        ConfigLoadError: If configuration loading or validation fails
    """
    try:
        config = load_user_config(config_file)
        remotes_data = config.get("remotes", {})

        if not isinstance(remotes_data, dict):
            raise ConfigLoadError("'remotes' section must be a dictionary")

        remotes = {}
        errors = []

        for name, remote_data in remotes_data.items():
            try:
                # Validate required fields
                if not isinstance(remote_data, dict):
                    errors.append(f"Remote '{name}' configuration must be a dictionary")
                    continue

                # Extract and validate fields
                url = remote_data.get("url")
                if not url:
                    errors.append(f"Remote '{name}' missing required field: url")
                    continue

                # Create RemoteConfig object with simplified schema
                remote_config = RemoteConfig(
                    name=name,
                    url=url,
                    platform=remote_data.get("platform"),  # Optional
                    conda_env=remote_data.get("conda_env"),  # Optional
                    working_dir=remote_data.get("working_dir"),  # Optional
                    config=remote_data.get("config"),  # Optional
                )

                # Validate the configuration
                validation_errors = remote_config.validate()
                if validation_errors:
                    for error in validation_errors:
                        errors.append(f"Remote '{name}': {error}")
                    continue

                remotes[name] = remote_config

            except RemoteConfigError as e:
                errors.append(f"Remote '{name}': {e}")
            except Exception as e:
                errors.append(f"Remote '{name}': Unexpected error: {e}")

        if errors:
            error_msg = "Configuration errors found:\\n" + "\\n".join(
                f"  - {error}" for error in errors
            )
            raise ConfigLoadError(error_msg)

        logger.info(
            f"Loaded {len(remotes)} remote configurations: {', '.join(remotes.keys())}"
        )
        return remotes

    except ConfigLoadError:
        raise
    except Exception as e:
        raise ConfigLoadError(f"Unexpected error loading remote configurations: {e}")


def get_remote_config(
    remote_name: str, config_file: Optional[str] = None
) -> RemoteConfig:
    """
    Get configuration for a specific remote.

    Args:
        remote_name: Name of the remote to retrieve
        config_file: Optional path to alternative config file

    Returns:
        RemoteConfig object for the specified remote

    Raises:
        ConfigLoadError: If remote not found or configuration invalid
    """
    remotes = load_remote_configurations(config_file)

    if remote_name not in remotes:
        available = list(remotes.keys())
        if available:
            available_str = ", ".join(available)
            raise ConfigLoadError(
                f"Remote '{remote_name}' not found. "
                f"Available remotes: {available_str}"
            )
        else:
            config_path = get_user_config_path()
            raise ConfigLoadError(
                f"Remote '{remote_name}' not found. "
                f"No remotes configured in {config_path}"
            )

    return remotes[remote_name]


def validate_remote_config(remote_name: str) -> tuple[bool, list[str]]:
    """
    Validate a remote configuration without raising exceptions.

    Args:
        remote_name: Name of the remote to validate

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    try:
        remote_config = get_remote_config(remote_name)
        errors = remote_config.validate()
        return len(errors) == 0, errors
    except ConfigLoadError as e:
        return False, [str(e)]
    except Exception as e:
        return False, [f"Unexpected error: {e}"]


def create_example_config(config_path: Optional[Path] = None) -> None:
    """
    Create an example configuration file.

    Args:
        config_path: Path to create config file (defaults to user config path)
    """
    if config_path is None:
        config_path = get_user_config_path()

    # Create directory if it doesn't exist
    config_path.parent.mkdir(parents=True, exist_ok=True)

    example_config = {
        "remotes": {
            "example_cluster": {
                "url": "ssh://cluster.example.edu",
                "qxub_env": "qxub",
                "platform_file": "/shared/qxub/platforms/cluster.yaml",
                "project_root_dir": "/home/${USER}/projects",
            }
        }
    }

    with open(config_path, "w") as f:
        yaml.dump(example_config, f, default_flow_style=False, sort_keys=False)

    logger.info(f"Created example configuration file: {config_path}")


def list_remotes() -> Dict[str, str]:
    """
    List all configured remotes with their URLs.

    Returns:
        Dictionary mapping remote names to their URLs
    """
    try:
        remotes = load_remote_configurations()
        return {name: config.url for name, config in remotes.items()}
    except ConfigLoadError:
        return {}
