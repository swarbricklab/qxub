"""
Remote execution backends for qxub v2.2.

This module provides the execution backends for different protocol        # Show SSH command in verbose mode
        if verbose >= 2:
            print(f"🔗 SSH connection: ssh {self.config.hostname}", file=sys.stderr)
            # Show the command with proper quoting for display
            ssh_display = ssh_command[:-1] + [shlex.quote(ssh_command[-1])]
            print(f"🔧 SSH command: {' '.join(ssh_display)}", file=sys.stderr)tarting with SSH-based remote execution.
"""

import logging
import select
import shlex
import subprocess
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, List, Optional

from .remote_config import RemoteConfig

logger = logging.getLogger(__name__)


class RemoteExecutionError(Exception):
    """Base exception for remote execution errors."""

    pass


class ConnectionError(RemoteExecutionError):
    """Raised when connection to remote system fails."""

    def __init__(self, message: str, protocol: str = "ssh"):
        super().__init__(message)
        self.protocol = protocol
        self.suggestions = self._get_suggestions()

    def _get_suggestions(self) -> List[str]:
        """Get protocol-specific troubleshooting suggestions."""
        if self.protocol == "ssh":
            return [
                "Check SSH configuration in ~/.ssh/config",
                "Verify network connectivity and VPN if required",
                "Test connection manually: ssh hostname echo 'test'",
                "Check SSH key permissions: chmod 600 ~/.ssh/id_*",
            ]
        return []


class RemoteExecutor(ABC):
    """Abstract base class for remote execution backends."""

    def __init__(self, config: RemoteConfig):
        self.config = config

    @abstractmethod
    def execute(
        self,
        command: str,
        working_dir: str,
        stream_output: bool = True,
        verbose: int = 0,
    ) -> int:
        """
        Execute command on remote system.

        Args:
            command: Command to execute
            working_dir: Remote working directory
            stream_output: Whether to stream output in real-time
            verbose: Verbosity level for execution details

        Returns:
            Exit code from remote execution
        """
        pass

    @abstractmethod
    def test_connection(self) -> bool:
        """
        Test if connection to remote system is available.

        Returns:
            True if connection successful, False otherwise
        """
        pass


class SSHRemoteExecutor(RemoteExecutor):
    """SSH-based remote execution."""

    def __init__(self, config: RemoteConfig):
        super().__init__(config)
        if config.protocol != "ssh":
            raise ValueError(
                f"SSHRemoteExecutor requires SSH protocol, got: {config.protocol}"
            )

    def _should_allocate_tty(self) -> bool:
        """
        Auto-detect if TTY allocation would be beneficial.

        Allocates TTY if local session is interactive, which preserves
        progress indicators, colors, and other TTY-dependent features
        in remote qxub execution.

        Returns:
            True if TTY should be allocated, False otherwise
        """
        try:
            # Allocate TTY if local session has both stdout and stderr as TTYs
            # This preserves the interactive experience for remote execution
            return sys.stdout.isatty() and sys.stderr.isatty()
        except (AttributeError, OSError):
            # Fallback to False if TTY detection fails
            return False

    def execute(
        self,
        command: str,
        working_dir: str,
        stream_output: bool = True,
        verbose: int = 0,
    ) -> int:
        """Execute command via SSH."""
        ssh_command = self._build_ssh_command(command, working_dir)

        # Log the command for debugging
        logger.info(
            f"Executing SSH command: {' '.join(ssh_command[:3])} ... (command truncated)"
        )
        logger.debug(f"Full SSH command: {ssh_command}")

        # Show SSH command in verbose mode
        if verbose >= 2:
            print(f"� SSH connection: ssh {self.config.hostname}", file=sys.stderr)
            print(f"�🔧 SSH command: {' '.join(ssh_command)}", file=sys.stderr)

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

    def _build_ssh_command(self, command: str, working_dir: str) -> List[str]:
        """Build SSH command with proper options."""
        ssh_cmd = ["ssh"]

        # Add config file if specified and exists
        if self.config.config and Path(self.config.config).exists():
            ssh_cmd.extend(["-F", self.config.config])

        # Add port if specified in URL
        if self.config.port:
            ssh_cmd.extend(["-p", str(self.config.port)])

        # TTY allocation logic
        should_allocate_tty = False
        if self.config.force_tty is True:
            # Explicitly requested
            should_allocate_tty = True
        elif self.config.force_tty is False:
            # Explicitly disabled
            should_allocate_tty = False
        else:
            # Auto-detect (None)
            should_allocate_tty = self._should_allocate_tty()

        if should_allocate_tty:
            ssh_cmd.append("-t")

        # Add connection options for better reliability
        ssh_cmd.extend(
            [
                "-o",
                "BatchMode=yes",  # Don't prompt for passwords
                "-o",
                "StrictHostKeyChecking=yes",  # Security: require known hosts
            ]
        )

        # Add hostname (with username if specified in URL)
        if self.config.username:
            ssh_cmd.append(f"{self.config.username}@{self.config.hostname}")
        else:
            ssh_cmd.append(self.config.hostname)

        # Build remote command
        remote_cmd = self._build_remote_command(command, working_dir)
        ssh_cmd.append(remote_cmd)

        return ssh_cmd

    def _build_remote_command(self, command: str, working_dir: str) -> str:
        """Build the command to execute on the remote system."""
        commands = []

        # Change to working directory
        commands.append(f"cd {working_dir}")

        # Initialize conda environment if specified
        if self.config.conda_env:
            # Initialize conda properly for non-interactive shells
            commands.append('eval "$(conda shell.bash hook)"')
            # Activate conda environment
            commands.append(f"conda activate {self.config.conda_env}")

        # Set platform override if specified (otherwise let remote auto-detect)
        if self.config.platform:
            commands.append(f"export QXUB_PLATFORM_OVERRIDE={self.config.platform}")

        # Execute the actual command
        commands.append(command)

        # Join with && to ensure all commands succeed
        return " && ".join(commands)

    def _execute_with_streaming(self, ssh_command: List[str]) -> int:
        """Execute SSH command with real-time output streaming."""
        try:
            process = subprocess.Popen(
                ssh_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,  # Line buffered
                universal_newlines=True,
            )

            # Stream output in real-time
            while True:
                # Use select to check for available output
                if hasattr(select, "select"):
                    # Unix-like systems
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
                    # Windows fallback - less efficient but works
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

    def test_connection(self) -> bool:
        """Test SSH connection to remote host."""
        test_cmd = ["ssh"]

        # Add config file if specified
        if self.config.config and Path(self.config.config).exists():
            test_cmd.extend(["-F", self.config.config])

        # Add port if specified
        if self.config.port:
            test_cmd.extend(["-p", str(self.config.port)])

        # Add connection options
        test_cmd.extend(
            [
                "-o",
                "BatchMode=yes",
                "-o",
                "ConnectTimeout=10",
                "-o",
                "StrictHostKeyChecking=yes",
            ]
        )

        # Add hostname
        if self.config.username:
            test_cmd.append(f"{self.config.username}@{self.config.hostname}")
        else:
            test_cmd.append(self.config.hostname)

        # Simple test command
        test_cmd.append("echo connection_test")

        logger.debug(f"Testing SSH connection with command: {' '.join(test_cmd)}")

        try:
            result = subprocess.run(
                test_cmd, capture_output=True, timeout=15, text=True
            )

            if result.returncode == 0 and "connection_test" in result.stdout:
                logger.debug(
                    f"SSH connection test successful to {self.config.hostname}"
                )
                return True
            else:
                logger.warning(
                    f"SSH connection test failed to {self.config.hostname}: "
                    f"returncode={result.returncode}, stdout='{result.stdout.strip()}', "
                    f"stderr='{result.stderr.strip()}'"
                )
                # Store error details for better error messages
                self._last_connection_error = {
                    "returncode": result.returncode,
                    "stdout": result.stdout.strip(),
                    "stderr": result.stderr.strip(),
                    "command": " ".join(test_cmd),
                }
                return False

        except subprocess.TimeoutExpired:
            logger.warning(f"SSH connection test timed out for {self.config.hostname}")
            self._last_connection_error = {
                "error": "timeout",
                "message": "Connection timed out after 15 seconds",
            }
            return False
        except FileNotFoundError:
            logger.error("SSH command not found")
            self._last_connection_error = {
                "error": "ssh_not_found",
                "message": "SSH command not found. Please install OpenSSH client.",
            }
            return False
        except Exception as e:
            logger.warning(f"SSH connection test failed: {e}")
            self._last_connection_error = {"error": "exception", "message": str(e)}
            return False

    def get_connection_error_details(self) -> dict:
        """Get detailed information about the last connection failure."""
        return getattr(self, "_last_connection_error", {})

    def test_remote_qxub(self) -> tuple[bool, str]:
        """
        Test if qxub is available in the remote conda environment.

        Returns:
            Tuple of (success, error_message)
        """
        test_cmd = self._build_ssh_command("qxub --version", "~")

        try:
            result = subprocess.run(
                test_cmd, capture_output=True, timeout=30, text=True
            )

            if result.returncode == 0:
                return True, ""
            else:
                return False, f"qxub command failed: {result.stderr.strip()}"

        except subprocess.TimeoutExpired:
            return False, "Timeout testing remote qxub installation"
        except Exception as e:
            return False, f"Error testing remote qxub: {e}"


class RemoteExecutorFactory:
    """Factory for creating protocol-specific remote executors."""

    @staticmethod
    def create(config: RemoteConfig) -> RemoteExecutor:
        """
        Create appropriate executor for the given configuration.

        Args:
            config: Remote configuration

        Returns:
            Remote executor instance

        Raises:
            UnsupportedProtocolError: If protocol is not supported
        """
        executors = {
            "ssh": SSHRemoteExecutor,
            # Future protocols can be added here:
            # 'aws': AWSBatchExecutor,
            # 'k8s': KubernetesExecutor,
        }

        executor_class = executors.get(config.protocol)
        if not executor_class:
            supported = ", ".join(executors.keys())
            from .remote_config import UnsupportedProtocolError

            raise UnsupportedProtocolError(
                f"Protocol '{config.protocol}' not supported. "
                f"Supported protocols: {supported}"
            )

        return executor_class(config)
