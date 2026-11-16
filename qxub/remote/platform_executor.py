"""
Platform-aware remote execution adapter for v3.3.0.

This module provides a simplified remote executor that works directly with
platform configuration dictionaries, replacing the v2.2 URL-based RemoteConfig
approach with SSH hostname delegation to ~/.ssh/config.
"""

import logging
import os
import random
import re
import select
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# PBS job IDs look like "12345678.gadi-pbs" or just "12345678". If we ever see
# one in the remote stdout we know the submission succeeded and any retry
# would risk a duplicate submission.
_PBS_JOB_ID_RE = re.compile(r"\b(\d{4,})\.[A-Za-z][\w-]*\b")

# SSH transport-level failure: the ssh client itself returns 255 for things
# like "Connection closed by <host>", auth failure, network reset, etc. The
# remote command is guaranteed not to have started in this case (or to have
# been killed before completing), so retrying is generally safe — *unless*
# we already saw a PBS job id in stdout, in which case the job was submitted.
_SSH_TRANSPORT_EXIT = 255


# Defaults for SSH robustness and retry behaviour. All can be overridden via
# environment variables so CI runners can tune behaviour without code changes.
_DEFAULT_CONNECT_TIMEOUT = 10  # seconds, per SSH ConnectTimeout / probe timeout
_DEFAULT_MAX_RETRIES = 3  # total connection probe attempts
_DEFAULT_RETRY_DELAY = 5.0  # seconds, base for exponential backoff
_DEFAULT_RETRY_MAX_DELAY = 60.0  # seconds, cap on backoff sleep
_DEFAULT_SERVER_ALIVE_INTERVAL = 15  # seconds between SSH keepalive probes
_DEFAULT_SERVER_ALIVE_COUNT_MAX = 4  # missed keepalives before SSH disconnects

# Watchdog for the *real* remote command. Some CI runs see the SSH session go
# completely silent mid-stream even with ServerAlive* set (the remote qxub or
# its conda init can wedge). If no output arrives on stdout/stderr for
# IDLE_TIMEOUT seconds we kill SSH and retry up to EXEC_MAX_RETRIES times.
# Set QXUB_REMOTE_IDLE_TIMEOUT=0 to disable the watchdog entirely.
_DEFAULT_IDLE_TIMEOUT = 600  # seconds of silence before treating as hung
_DEFAULT_EXEC_MAX_RETRIES = 2  # extra attempts after the first try (total 3)


def _env_int(name: str, default: int) -> int:
    """Read a positive integer from an env var, falling back to default."""
    raw = os.environ.get(name)
    if not raw:
        return default
    try:
        value = int(raw)
        return value if value >= 0 else default
    except ValueError:
        logger.warning("Invalid value for %s=%r, using default %s", name, raw, default)
        return default


def _env_float(name: str, default: float) -> float:
    """Read a non-negative float from an env var, falling back to default."""
    raw = os.environ.get(name)
    if not raw:
        return default
    try:
        value = float(raw)
        return value if value >= 0 else default
    except ValueError:
        logger.warning("Invalid value for %s=%r, using default %s", name, raw, default)
        return default


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

        idle_timeout = _env_int("QXUB_REMOTE_IDLE_TIMEOUT", _DEFAULT_IDLE_TIMEOUT)
        max_retries = _env_int(
            "QXUB_REMOTE_EXEC_MAX_RETRIES", _DEFAULT_EXEC_MAX_RETRIES
        )
        base_delay = _env_float("QXUB_REMOTE_RETRY_DELAY", _DEFAULT_RETRY_DELAY)
        max_delay = _env_float("QXUB_REMOTE_RETRY_MAX_DELAY", _DEFAULT_RETRY_MAX_DELAY)

        total_attempts = max_retries + 1
        last_exit = 1

        for attempt in range(1, total_attempts + 1):
            try:
                if stream_output:
                    exit_code, hung, captured_stdout = self._execute_with_streaming(
                        ssh_command, idle_timeout=idle_timeout
                    )
                else:
                    # Non-streaming path: rely on subprocess timeout = idle_timeout
                    # acting as a wall-clock cap. Treat timeout as a hang.
                    hung = False
                    captured_stdout = ""
                    try:
                        result = subprocess.run(
                            ssh_command,
                            capture_output=True,
                            text=True,
                            timeout=idle_timeout if idle_timeout > 0 else None,
                        )
                        if result.stdout:
                            print(result.stdout, end="")
                            captured_stdout = result.stdout
                        if result.stderr:
                            print(result.stderr, end="", file=sys.stderr)
                        exit_code = result.returncode
                    except subprocess.TimeoutExpired:
                        hung = True
                        exit_code = 124  # conventional timeout exit code

            except FileNotFoundError:
                raise RemoteExecutionError(
                    "SSH command not found. Please install OpenSSH client."
                )
            except subprocess.SubprocessError as e:
                raise RemoteExecutionError(f"SSH execution failed: {e}")

            last_exit = exit_code

            # Did the remote get far enough to submit a PBS job?
            job_id_match = _PBS_JOB_ID_RE.search(captured_stdout or "")
            submitted_job_id = job_id_match.group(0) if job_id_match else None

            # Success path: exit 0 and not hung.
            if not hung and exit_code == 0:
                return exit_code

            # Classify the failure.
            ssh_transport_failure = not hung and exit_code == _SSH_TRANSPORT_EXIT
            retryable = hung or ssh_transport_failure

            if not retryable:
                # Real remote failure (non-zero exit from the user command).
                # Don't retry — the user's command itself failed.
                return exit_code

            # If a job was already submitted, retrying would duplicate it.
            if submitted_job_id is not None:
                print(
                    f"⚠️  SSH connection to {self.hostname} dropped "
                    f"(exit {exit_code})"
                    + (f" after {idle_timeout}s of silence" if hung else "")
                    + f", but PBS job {submitted_job_id} was already "
                    f"submitted. NOT retrying to avoid duplicate "
                    f"submission. Check job status with: "
                    f"`ssh {self.hostname} qstat {submitted_job_id}`",
                    file=sys.stderr,
                )
                return exit_code

            # No job submitted — safe to retry.
            if attempt >= total_attempts:
                reason = (
                    f"no output for {idle_timeout}s"
                    if hung
                    else f"SSH transport failure (exit {exit_code})"
                )
                print(
                    f"❌ Remote command on {self.hostname} failed: {reason}. "
                    f"Giving up after {attempt} attempt(s).",
                    file=sys.stderr,
                )
                return exit_code

            delay = min(base_delay * (2 ** (attempt - 1)), max_delay)
            # Small jitter to avoid thundering-herd on shared CI runners.
            delay += random.uniform(0, min(1.0, delay * 0.1))
            reason = (
                f"produced no output for {idle_timeout}s (assuming hang)"
                if hung
                else f"SSH dropped the connection (exit {exit_code})"
            )
            print(
                f"⚠️  Remote command on {self.hostname} {reason}. "
                f"No PBS job id seen in output, so retrying is safe. "
                f"Retrying ({attempt}/{max_retries}) in {delay:.1f}s.",
                file=sys.stderr,
            )
            time.sleep(delay)

        return last_exit

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

        # Connection options - minimal, let ~/.ssh/config control the rest.
        # The robustness options below mitigate silent network drops that have
        # been observed when running from CI runners against NCI: ConnectTimeout
        # bounds dial time, and the ServerAlive* pair causes SSH to give up
        # rather than hang indefinitely if the connection stalls mid-command.
        connect_timeout = _env_int(
            "QXUB_REMOTE_CONNECT_TIMEOUT", _DEFAULT_CONNECT_TIMEOUT
        )
        alive_interval = _env_int(
            "QXUB_REMOTE_SERVER_ALIVE_INTERVAL", _DEFAULT_SERVER_ALIVE_INTERVAL
        )
        alive_count = _env_int(
            "QXUB_REMOTE_SERVER_ALIVE_COUNT_MAX", _DEFAULT_SERVER_ALIVE_COUNT_MAX
        )
        ssh_cmd.extend(
            [
                "-o",
                "BatchMode=yes",  # Don't prompt for passwords
                "-o",
                f"ConnectTimeout={connect_timeout}",
                "-o",
                f"ServerAliveInterval={alive_interval}",
                "-o",
                f"ServerAliveCountMax={alive_count}",
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

    def _execute_with_streaming(
        self, ssh_command: list, idle_timeout: int = 0
    ) -> tuple[int, bool, str]:
        """
        Execute SSH command with real-time output streaming and idle watchdog.

        Args:
            ssh_command: SSH command as list of arguments
            idle_timeout: Seconds of no stdout/stderr output before the remote
                command is considered hung and the SSH process is terminated.
                Set to 0 to disable the watchdog.

        Returns:
            Tuple of (exit_code, hung, captured_stdout) where ``hung`` is True
            if the watchdog killed the process due to no output, and
            ``captured_stdout`` is the full remote stdout (used by the caller
            to detect whether a PBS job id was emitted before any failure).
        """
        process = None
        hung = False
        stdout_buf: List[str] = []
        try:
            process = subprocess.Popen(
                ssh_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,  # Line buffered
                universal_newlines=True,
            )

            last_output = time.monotonic()

            # Stream output in real-time using select
            while True:
                if hasattr(select, "select"):
                    # Unix-like systems - efficient multiplexing
                    ready, _, _ = select.select(
                        [process.stdout, process.stderr], [], [], 0.5
                    )

                    for stream in ready:
                        if stream == process.stdout:
                            line = process.stdout.readline()
                            if line:
                                print(line, end="")
                                stdout_buf.append(line)
                                last_output = time.monotonic()
                        elif stream == process.stderr:
                            line = process.stderr.readline()
                            if line:
                                print(line, end="", file=sys.stderr)
                                last_output = time.monotonic()
                else:
                    # Windows fallback - less efficient but functional
                    while True:
                        output = process.stdout.readline()
                        if output:
                            print(output, end="")
                            stdout_buf.append(output)
                            last_output = time.monotonic()
                        else:
                            break

                # Idle watchdog — only fires while the process is still alive.
                if (
                    idle_timeout > 0
                    and process.poll() is None
                    and (time.monotonic() - last_output) > idle_timeout
                ):
                    hung = True
                    logger.warning(
                        "No output from remote for %ss; terminating SSH",
                        idle_timeout,
                    )
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait(timeout=5)
                    break

                # Check if process has finished
                if process.poll() is not None:
                    break

            # Read any remaining output
            stdout, stderr = process.communicate()
            if stdout:
                print(stdout, end="")
                stdout_buf.append(stdout)
            if stderr:
                print(stderr, end="", file=sys.stderr)

            exit_code = process.returncode if process.returncode is not None else 124
            return exit_code, hung, "".join(stdout_buf)

        except KeyboardInterrupt:
            logger.info("Received interrupt signal, terminating remote process")
            if process is not None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
            return 130, False, "".join(stdout_buf)  # Standard exit code for SIGINT

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
