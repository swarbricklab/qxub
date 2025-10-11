"""
Remote execution backends for qxub v2.2.

This module provides the execution backends for different protocols,
starting with SSH-based remote execution.
"""

import logging
import select
import subprocess
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import IO, Any, List, Optional

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
        self, command: str, working_dir: str, stream_output: bool = True
    ) -> int:
        """
        Execute command on remote system.

        Args:
            command: Command to execute
            working_dir: Remote working directory
            stream_output: Whether to stream output in real-time

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

    def execute(
        self, command: str, working_dir: str, stream_output: bool = True
    ) -> int:
        """Execute command via SSH."""
        ssh_command = self._build_ssh_command(command, working_dir)

        logger.info(
            f"Executing SSH command: {' '.join(ssh_command[:3])} ... (command truncated)"
        )
        logger.debug(f"Full SSH command: {ssh_command}")

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

        # Activate conda environment
        commands.append(f"conda activate {self.config.qxub_env}")

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

        try:
            result = subprocess.run(
                test_cmd, capture_output=True, timeout=15, text=True
            )
            return result.returncode == 0 and "connection_test" in result.stdout

        except subprocess.TimeoutExpired:
            logger.warning(f"SSH connection test timed out for {self.config.hostname}")
            return False
        except FileNotFoundError:
            logger.error("SSH command not found")
            return False
        except Exception as e:
            logger.warning(f"SSH connection test failed: {e}")
            return False

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
