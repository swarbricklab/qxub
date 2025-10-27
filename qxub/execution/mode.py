"""
Execution mode detection for qxub.

Determines whether job execution should be local (direct PBS submission)
or remote (SSH to remote platform first).
"""

import logging
from enum import Enum
from typing import Optional

logger = logging.getLogger(__name__)


class ExecutionMode(Enum):
    """Execution mode for job submission."""

    LOCAL = "local"
    REMOTE = "remote"


def get_execution_mode(
    platform_config: Optional[dict],
) -> ExecutionMode:
    """
    Determine execution mode from platform configuration.

    The execution mode is determined by the presence of a 'remote' section
    in the platform configuration:
    - Has 'remote' section → REMOTE execution (SSH to platform)
    - No 'remote' section → LOCAL execution (direct PBS submission)

    Args:
        platform_config: Platform configuration dictionary from config manager

    Returns:
        ExecutionMode.LOCAL or ExecutionMode.REMOTE

    Examples:
        >>> # Local platform config (no remote section)
        >>> config = {
        ...     'name': 'nci_gadi',
        ...     'definition': 'file:///apps/qxub/platforms/nci_gadi.yaml'
        ... }
        >>> get_execution_mode(config)
        ExecutionMode.LOCAL

        >>> # Remote platform config (has remote section)
        >>> config = {
        ...     'name': 'nci_gadi',
        ...     'definition': 'https://github.com/.../nci_gadi.yaml',
        ...     'remote': {
        ...         'host': 'gadi',
        ...         'working_dir': '/scratch/a56/{user}'
        ...     }
        ... }
        >>> get_execution_mode(config)
        ExecutionMode.REMOTE
    """
    # If no platform config, default to local
    if not platform_config:
        logger.debug("No platform config provided, defaulting to local execution")
        return ExecutionMode.LOCAL

    # Check for remote section
    if "remote" in platform_config and platform_config["remote"]:
        logger.debug(
            f"Platform {platform_config.get('name', 'unknown')} has remote config, "
            "using remote execution"
        )
        return ExecutionMode.REMOTE

    logger.debug(
        f"Platform {platform_config.get('name', 'unknown')} has no remote config, "
        "using local execution"
    )
    return ExecutionMode.LOCAL


def validate_remote_config(remote_config: dict) -> list[str]:
    """
    Validate remote execution configuration.

    Args:
        remote_config: Remote section of platform config

    Returns:
        List of validation error messages (empty if valid)
    """
    errors = []

    if not remote_config:
        errors.append("Remote config is empty")
        return errors

    # Required fields
    if "host" not in remote_config:
        errors.append("Remote config missing required field: host")

    # Optional but recommended fields
    if "working_dir" not in remote_config:
        logger.warning("Remote config does not specify working_dir")

    return errors


def get_remote_config(platform_config: dict) -> Optional[dict]:
    """
    Extract remote configuration from platform config.

    Args:
        platform_config: Platform configuration dictionary

    Returns:
        Remote config dictionary or None if not present
    """
    if not platform_config:
        return None

    return platform_config.get("remote")
