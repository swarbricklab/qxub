"""
Integration of remote execution with existing platform system.

This shows how the remote execution capability integrates with
the existing Platform class and queue management.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..remote import RemoteConfig, execute_remote_platform
from .core import Platform

logger = logging.getLogger(__name__)


def load_platform_with_remote(platform_data: Dict[str, Any]) -> Platform:
    """
    Load platform configuration including remote execution settings.

    This extends the existing platform loading to handle remote configurations.
    """

    # Load basic platform (existing functionality)
    platform = Platform(
        name=platform_data["platform"]["name"],
        type=platform_data["platform"]["type"],
        host=platform_data["platform"]["host"],
        description=platform_data["platform"].get("description", ""),
        # ... queue loading etc. (existing code)
    )

    # Add remote configuration if present
    if "remote" in platform_data["platform"]:
        remote_data = platform_data["platform"]["remote"]
        platform.remote_config = RemoteConfig(
            host=remote_data["host"],
            port=remote_data.get("port", 22),
            user=remote_data.get("user"),
            ssh_key=remote_data.get("ssh_key"),
            connection_timeout=remote_data.get("connection_timeout", 30),
            qxub_path=remote_data.get("qxub_path", "qxub"),
            remote_cwd=remote_data.get("remote_cwd", "preserve_relative"),
            remote_base_path=remote_data.get("remote_base_path"),
            remote_setup_commands=remote_data.get("remote_setup_commands", []),
        )
    else:
        platform.remote_config = None

    return platform


def execute_with_platform(platform: Platform, args: List[str], cwd: Path) -> int:
    """
    Execute qxub command, handling both local and remote platforms.

    This is the main integration point that decides whether to execute
    locally or remotely based on platform configuration.
    """

    # Check if this platform requires remote execution
    if hasattr(platform, "remote_config") and platform.remote_config:
        logger.info(f"Executing on remote platform: {platform.name} @ {platform.host}")

        # Build platform config dict for remote execution
        platform_config = {
            "remote": {
                "host": platform.remote_config.host,
                "port": platform.remote_config.port,
                "user": platform.remote_config.user,
                "ssh_key": platform.remote_config.ssh_key,
                "connection_timeout": platform.remote_config.connection_timeout,
                "qxub_path": platform.remote_config.qxub_path,
                "remote_cwd": platform.remote_config.remote_cwd,
                "remote_base_path": platform.remote_config.remote_base_path,
                "remote_setup_commands": platform.remote_config.remote_setup_commands,
            }
        }

        # Execute remotely
        return execute_remote_platform(platform_config, args, cwd)

    else:
        # Local execution (existing qxub functionality)
        logger.info(f"Executing on local platform: {platform.name}")
        return execute_local_platform(platform, args, cwd)


def execute_local_platform(platform: Platform, args: List[str], cwd: Path) -> int:
    """Execute qxub command locally (existing functionality)."""
    # This would contain the existing qxub execution logic
    # Just a placeholder for now
    logger.info("Local execution - would call existing qxub logic")
    return 0


# Example usage and platform configuration


def create_example_remote_platform():
    """Example of how a remote platform would be configured."""

    platform_yaml = """
    platform:
      name: nci_gadi_remote
      type: pbs_pro
      host: gadi.nci.org.au
      description: "NCI Gadi via remote SSH execution"

      # Remote execution configuration
      remote:
        host: gadi.nci.org.au
        port: 22
        user: jr9959
        ssh_key: ~/.ssh/id_rsa
        connection_timeout: 30
        qxub_path: /g/data/a56/software/qsub_tools/venv/bin/qxub
        remote_cwd: preserve_relative
        remote_base_path: /g/data/a56
        remote_setup_commands:
          - "module load python3"
          - "source /g/data/a56/software/qsub_tools/venv/bin/activate"

      queues:
        normal:
          type: standard
          limits:
            max_cores: 48
            max_memory_gb: 192
            max_walltime: "48:00:00"
          priority: normal
          su_billing_rate: 1.0

        gpu:
          type: gpu
          limits:
            max_cores: 48
            max_memory_gb: 192
            max_gpus: 4
            max_walltime: "24:00:00"
          priority: high
          su_billing_rate: 3.0

      auto_selection:
        - condition: "gpus > 0"
          queue: gpu
        - condition: "cores <= 48"
          queue: normal
          is_default: true
    """

    return platform_yaml


# CLI Integration Example


def enhanced_cli_execution(args: List[str]):
    """
    Enhanced CLI execution that handles remote platforms.

    This shows how the CLI would be modified to support remote execution.
    """

    # Parse arguments to identify platform
    platform_name = extract_platform_from_args(args)

    if not platform_name:
        print("No platform specified, using local execution")
        return execute_local_command(args)

    # Load platform configuration
    platform = load_platform_configuration(platform_name)

    # Execute with platform (automatically handles local vs remote)
    current_dir = Path.cwd()
    return execute_with_platform(platform, args, current_dir)


def extract_platform_from_args(args: List[str]) -> Optional[str]:
    """Extract platform name from command line arguments."""
    for i, arg in enumerate(args):
        if arg == "--platform" and i + 1 < len(args):
            return args[i + 1]
        if arg.startswith("--platform="):
            return arg.split("=", 1)[1]
    return None


def load_platform_configuration(platform_name: str) -> Platform:
    """Load platform configuration from YAML files."""
    # This would load from docs/platforms/{platform_name}.yaml
    # and call load_platform_with_remote()
    pass


def execute_local_command(args: List[str]) -> int:
    """Execute command locally without platform."""
    # Existing qxub local execution
    pass
