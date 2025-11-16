"""
Platform-aware remote execution adapter for v3.3.0.

This module provides a simplified remote executor that works directly with
platform configuration dictionaries, replacing the v2.2 URL-based RemoteConfig
approach with SSH hostname delegation to ~/.ssh/config.
"""

import logging
import select
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


class RemoteExecutionError(Exception):
    """Base exception for remote execution errors."""

    pass


class PlatformRemoteExecutor:
    """
    SSH-based remote executor for platform-aware execution.

    Uses platform configuration with remote: section containing:
    - host: SSH hostname (from ~/.ssh/config)
    - working_dir: Remote working directory (with variable substitution)
    - conda_init: Optional custom conda initialization commands
    """

    def __init__(self, platform_name: str, remote_config: Dict[str, Any]):
        """
        Initialize remote executor from platform config.

        Args:
            platform_name: Name of the platform (for QXUB_PLATFORM env var)
            remote_config: The 'remote:' section from platform config containing:
                - host: SSH hostname
                - working_dir: Remote working directory
                - conda_init: Optional conda initialization (default: standard hook)

        Raises:
            RemoteExecutionError: If required config fields are missing
        """
        self.platform_name = platform_name
        self.hostname = remote_config.get("host")
        self.working_dir = remote_config.get("working_dir")
        self.conda_init = remote_config.get("conda_init")

        if not self.hostname:
            raise RemoteExecutionError(
                f"Platform '{platform_name}' remote config missing required 'host' field"
            )

        if not self.working_dir:
            raise RemoteExecutionError(
                f"Platform '{platform_name}' remote config missing required 'working_dir' field"
            )

        # Expand variables in working_dir (e.g., {user}, {project})
        self.working_dir = self._expand_variables(self.working_dir)

    def _expand_variables(self, path: str) -> str:
        """
        Expand template variables in path strings.

        Follows the same pattern as compute node variable resolution:
        - {var}: Resolved on config system (laptop/CI environment)
        - {{var}}: Resolved on remote system (Gadi execution environment)

        Supported variables:
        - {user}: Current username from config system (laptop/CI user)
        - {{user}}: Deferred to remote host as $USER (remote user)
        - {project}: $PROJECT environment variable (expanded locally)
        - {{project}}: Deferred to remote host as $PROJECT (remote env)

        Args:
            path: Path string with optional variables

        Returns:
            Path string with local variables expanded and remote variables deferred
        """
        import getpass
        import os

        # First handle double-brace variables (remote evaluation)
        # {{user}} -> $USER (will be resolved by remote shell)
        result = path.replace("{{user}}", "$USER")

        # {{project}} -> $PROJECT (will be resolved by remote shell)
        result = result.replace("{{project}}", "$PROJECT")

        # Then handle single-brace variables (local evaluation)
        # {user} -> actual username from config system
        local_user = getpass.getuser()
        result = result.replace("{user}", local_user)

        # {project} -> actual project from local environment
        local_project = os.environ.get("PROJECT", "")
        result = result.replace("{project}", local_project)

        return result

    def execute(
        self,
        remote_command: str,
        stream_output: bool = True,
        verbose: int = 0,
    ) -> int:
        """
        Execute qxub command on remote platform via SSH.

        Args:
            remote_command: The qxub command to execute (e.g., "qxub exec --env pytorch -- python train.py")
            stream_output: Whether to stream output in real-time (default: True)
            verbose: Verbosity level for execution details

        Returns:
            Exit code from remote execution

        Raises:
            RemoteExecutionError: If SSH execution fails
        """
        ssh_command = self._build_ssh_command(remote_command)

        # Log the command for debugging - show full command for troubleshooting
        logger.info(f"Executing SSH command to {self.hostname}: {remote_command}")
        logger.debug(f"Full SSH command: {ssh_command}")

        # Show SSH connection info in verbose mode
        if verbose >= 2:
            print(f"🌐 SSH connection: ssh {self.hostname}", file=sys.stderr)
            print(f"📁 Remote directory: {self.working_dir}", file=sys.stderr)
            print(f"🔧 Remote command: {remote_command}", file=sys.stderr)
        elif verbose >= 1:
            print(f"🌐 Executing on {self.platform_name} via SSH", file=sys.stderr)

        try:
            if stream_output:
                return self._execute_with_streaming(ssh_command)
            else:
                result = subprocess.run(ssh_command, capture_output=True, text=True)
                if result.stdout:
                    print(result.stdout, end="")
                if result.stderr:
                    print(result.stderr, end="", file=sys.stderr)
                return result.returncode

        except FileNotFoundError:
            raise RemoteExecutionError(
                "SSH command not found. Please install OpenSSH client."
            )
        except subprocess.SubprocessError as e:
            raise RemoteExecutionError(f"SSH execution failed: {e}")

    def _build_ssh_command(self, remote_command: str) -> list:
        """
        Build SSH command with proper options.

        Uses ~/.ssh/config for all SSH configuration (keys, ports, usernames, etc).
        Only specifies hostname and command to execute.

        Args:
            remote_command: The qxub command to execute remotely

        Returns:
            List of SSH command arguments
        """
        ssh_cmd = ["ssh"]

        # TTY allocation for interactive experience
        if self._should_allocate_tty():
            ssh_cmd.append("-t")

        # Connection options - minimal, let ~/.ssh/config control the rest
        ssh_cmd.extend(
            [
                "-o",
                "BatchMode=yes",  # Don't prompt for passwords
            ]
        )

        # Add hostname (SSH will look up all settings in ~/.ssh/config)
        ssh_cmd.append(self.hostname)

        # Build remote command with environment setup
        wrapped_command = self._wrap_remote_command(remote_command)
        ssh_cmd.append(wrapped_command)

        return ssh_cmd

    def _wrap_remote_command(self, command: str) -> str:
        """
        Wrap remote command with necessary environment setup.

        Args:
            command: The qxub command to execute

        Returns:
            Complete command string with cd, conda init, env vars, etc.
        """
        commands = []

        # Change to working directory
        commands.append(f"cd {self.working_dir}")

        # Initialize conda if specified
        if self.conda_init:
            # Use custom conda initialization - handle multi-line strings
            # Split by newlines and add each line as a separate command
            init_lines = self.conda_init.strip().split("\n")
            for line in init_lines:
                line = line.strip()
                if line:  # Skip empty lines
                    commands.append(line)
        else:
            # Use standard conda initialization
            commands.append('eval "$(conda shell.bash hook)"')

        # Set platform environment variable for remote qxub
        commands.append(f"export QXUB_PLATFORM={self.platform_name}")

        # Force unbuffered Python output for real-time streaming over SSH
        commands.append("export PYTHONUNBUFFERED=1")

        # Execute the actual qxub command
        commands.append(command)

        # Join with && to ensure all commands succeed
        return " && ".join(commands)

    def _should_allocate_tty(self) -> bool:
        """
        Determine if TTY allocation would be beneficial for remote execution.

        For qxub remote execution, we generally want clean output without
        interactive progress indicators that can cause staggered output.
        TTY allocation can cause the remote qxub to output progress indicators
        to /dev/tty which interferes with clean command output.

        Returns:
            False - Disable TTY allocation for cleaner remote execution output
        """
        # Always return False for now to avoid staggered output issues
        # TODO: Consider making this configurable if interactive features are needed
        return False

    def _execute_with_streaming(self, ssh_command: list) -> int:
        """
        Execute SSH command with real-time output streaming.

        Args:
            ssh_command: SSH command as list of arguments

        Returns:
            Exit code from SSH process
        """
        try:
            process = subprocess.Popen(
                ssh_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,  # Line buffered
                universal_newlines=True,
            )

            # Stream output in real-time using select
            while True:
                if hasattr(select, "select"):
                    # Unix-like systems - efficient multiplexing
                    ready, _, _ = select.select(
                        [process.stdout, process.stderr], [], [], 0.1
                    )

                    for stream in ready:
                        if stream == process.stdout:
                            line = process.stdout.readline()
                            if line:
                                print(line, end="")
                        elif stream == process.stderr:
                            line = process.stderr.readline()
                            if line:
                                print(line, end="", file=sys.stderr)
                else:
                    # Windows fallback - less efficient but functional
                    while True:
                        output = process.stdout.readline()
                        if output:
                            print(output, end="")
                        else:
                            break

                # Check if process has finished
                if process.poll() is not None:
                    break

            # Read any remaining output
            stdout, stderr = process.communicate()
            if stdout:
                print(stdout, end="")
            if stderr:
                print(stderr, end="", file=sys.stderr)

            return process.returncode

        except KeyboardInterrupt:
            logger.info("Received interrupt signal, terminating remote process")
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
            return 130  # Standard exit code for SIGINT

    def test_connection(self) -> tuple[bool, str]:
        """
        Test SSH connection to remote host.

        Returns:
            Tuple of (success, error_message)
        """
        test_cmd = [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=10",
            self.hostname,
            "echo connection_test",
        ]

        logger.debug(f"Testing SSH connection: {' '.join(test_cmd)}")

        try:
            result = subprocess.run(
                test_cmd, capture_output=True, timeout=15, text=True
            )

            if result.returncode == 0 and "connection_test" in result.stdout:
                logger.debug(f"SSH connection test successful to {self.hostname}")
                return True, ""
            else:
                error_msg = (
                    f"SSH connection failed (exit code {result.returncode}): "
                    f"{result.stderr.strip()}"
                )
                logger.warning(error_msg)
                return False, error_msg

        except subprocess.TimeoutExpired:
            error_msg = f"SSH connection timed out after 15 seconds"
            logger.warning(f"{error_msg} for {self.hostname}")
            return False, error_msg
        except FileNotFoundError:
            error_msg = "SSH command not found. Please install OpenSSH client."
            logger.error(error_msg)
            return False, error_msg
        except Exception as e:
            error_msg = f"SSH connection test failed: {e}"
            logger.warning(error_msg)
            return False, error_msg
