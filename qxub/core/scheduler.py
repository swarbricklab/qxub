"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

import logging
import shlex
import subprocess
import sys
import threading
import time

import click

from ..resources import (
    bytes_to_human,
    calculate_efficiency,
    parse_exec_host,
    parse_exec_vnode,
    parse_timestamp,
    size_to_bytes,
    time_to_seconds,
)


def print_status(message, final=False):
    """Print a status message that overwrites the previous one"""
    import os

    # Check if we're in a remote SSH context
    is_remote_ssh = any(
        [
            os.environ.get("SSH_CLIENT"),
            os.environ.get("SSH_CONNECTION"),
            os.environ.get("SSH_TTY"),
        ]
    )

    try:
        with open("/dev/tty", "w") as tty:
            if is_remote_ssh:
                # In remote SSH: force cursor to start of line, then print message
                print(f"\r{message}", file=tty, flush=True)
            else:
                # Local TTY: use carriage returns for overwriting
                if final:
                    # Final message - print normally and move to next line
                    print(f"{message}", file=tty)
                else:
                    # Temporary message - overwrite without newline
                    print(f"\r{message}", end="", flush=True, file=tty)
    except (OSError, IOError):
        # Fallback to stderr if /dev/tty is not available (non-interactive context like SSH)
        # Important messages like job IDs need to be visible even without a TTY
        print(f"{message}", file=sys.stderr, flush=True)


class SimpleSpinner:
    """Simple non-threaded spinner for single-thread job monitoring."""

    def __init__(self, quiet=False):
        self.quiet = quiet or self._is_remote_ssh()
        self.spinner_chars = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
        self.current_char_index = 0

    def _is_remote_ssh(self):
        """Detect if we're running in a remote SSH session where spinners don't work well."""

        # Check for SSH environment variables that indicate remote execution
        # ssh_indicators = [
        #     os.environ.get("SSH_CLIENT"),
        #     os.environ.get("SSH_CONNECTION"),
        #     os.environ.get("SSH_TTY"),
        # ]

        # For now, allow spinners in SSH sessions for better UX
        # TODO: Make this configurable if needed
        return False  # any(ssh_indicators)

    def show_next_frame(self):
        """Show the next frame of the spinner animation."""
        if self.quiet:
            return

        char = self.spinner_chars[self.current_char_index % len(self.spinner_chars)]
        try:
            with open("/dev/tty", "w") as tty:
                print(f"\r{char}", end="", flush=True, file=tty)
        except (OSError, IOError):
            # Ignore if /dev/tty not available
            pass

        self.current_char_index += 1

    def clear(self):
        """Clear the spinner line."""
        if self.quiet:
            return

        try:
            with open("/dev/tty", "w") as tty:
                print("\r ", end="", flush=True, file=tty)
        except (OSError, IOError):
            # Ignore if /dev/tty not available
            pass


class JobSpinner:  # pylint: disable=too-many-instance-attributes
    """Context manager for displaying a spinner during job operations."""

    def __init__(self, message="", quiet=False, show_message=False, coordinator=None):
        self.message = message
        self.quiet = (
            quiet or self._is_remote_ssh()
        )  # Disable spinner in remote SSH contexts
        self.show_message = show_message  # Default to False (no messages)
        self.coordinator = coordinator
        self.spinner_chars = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
        self.spinning = False
        self.thread = None
        self.original_line_len = 0

    def _is_remote_ssh(self):
        """Detect if we're running in a remote SSH session where spinners don't work well."""
        import os

        # Check for SSH environment variables that indicate remote execution
        ssh_indicators = [
            os.environ.get("SSH_CLIENT"),
            os.environ.get("SSH_CONNECTION"),
            os.environ.get("SSH_TTY"),
        ]

        # If any SSH indicator is present, we're likely in an SSH session
        return any(ssh_indicators)

    def _spin(self):
        """Run the spinner animation"""
        i = 0
        while self.spinning:
            # Check if any event requires stopping the spinner
            if self.coordinator and self.coordinator.should_stop_spinner():
                self._clear_spinner()
                self.coordinator.signal_spinner_cleared()
                break

            char = self.spinner_chars[i % len(self.spinner_chars)]
            try:
                with open("/dev/tty", "w") as tty:
                    if self.show_message:
                        line = f"{self.message} {char}"
                        print(f"\r{line}", end="", flush=True, file=tty)
                        self.original_line_len = len(line)
                    else:
                        print(f"\r{char}", end="", flush=True, file=tty)
                        self.original_line_len = 1
            except (OSError, IOError):
                # Fallback to /dev/null if /dev/tty not available
                with open("/dev/null", "w") as devnull:
                    if self.show_message:
                        line = f"{self.message} {char}"
                        print(f"\r{line}", end="", flush=True, file=devnull)
                    else:
                        print(f"\r{char}", end="", flush=True, file=devnull)
            time.sleep(0.1)
            i += 1

    def _clear_spinner(self):
        """Clear the spinner line"""
        clear_line = " " * (self.original_line_len + 5)
        try:
            with open("/dev/tty", "w") as tty:
                print(f"\r{clear_line}\r", end="", flush=True, file=tty)
        except (OSError, IOError):
            # Fallback to /dev/null if /dev/tty not available
            with open("/dev/null", "w") as devnull:
                print(f"\r{clear_line}\r", end="", flush=True, file=devnull)

    def __enter__(self):
        if not self.quiet:
            self.spinning = True
            self.thread = threading.Thread(target=self._spin)
            self.thread.daemon = True
            self.thread.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.quiet and self.spinning:
            self.spinning = False
            if self.thread:
                self.thread.join(timeout=0.2)
            # Clear the spinner line completely if not already cleared
            if self.coordinator and not self.coordinator.spinner_cleared.is_set():
                self._clear_spinner()
                self.coordinator.signal_spinner_cleared()


def qsub(cmd, quiet=False):
    """
    Submits a job via qsub based on the given command (cmd)

    Args:
        cmd (chr): qsub command, including all of the options generated by qt
        quiet (bool): whether to suppress spinner output

    Returns:
        The job id for the submitted job
    """
    # Make sure that this is a qsub cmd -- don't pass on random commands!
    if not cmd[:4] == "qsub":
        click.echo("Expected qsub comment. Exiting")
        sys.exit(1)

    # Submit job directly without spinner (spinner is handled in monitor)
    # pylint: disable=W1510
    result = subprocess.run(
        shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        logging.debug("Job submission failed")
        click.echo(result.stderr)

        # Map PBS validation errors to exit code 2 for consistency
        if result.returncode == 166:  # PBS validation error
            sys.exit(2)
        else:
            sys.exit(result.returncode)
    else:
        logging.debug("Job submitted successfully")
        return result.stdout.rstrip("\n")


def qdel(job_id, quiet=False):
    """
    Deletes a PBS job using qdel command

    Args:
        job_id (str): The PBS job id to delete
        quiet (bool): whether to suppress output

    Returns:
        bool: True if deletion was successful, False otherwise
    """
    logging.debug("Deleting job %s", job_id)

    try:
        # pylint: disable=W1510
        result = subprocess.run(
            ["qdel", job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.returncode == 0:
            if not quiet:
                click.echo(f"🗑️  Job {job_id} deleted successfully")
            logging.info("Job %s deleted successfully", job_id)
            return True

        if not quiet:
            click.echo(f"Failed to delete job {job_id}: {result.stderr.strip()}")
        logging.warning("Failed to delete job %s: %s", job_id, result.stderr.strip())
        return False
    except Exception as e:  # pylint: disable=broad-except
        if not quiet:
            click.echo(f"Error deleting job {job_id}: {e}")
        logging.error("Error deleting job %s: %s", job_id, e)
        return False


def job_status(job_id):
    """
    Check the current status of a job.

    Returns:
        H: held (not enough quota)
    """
    logging.debug("Job id: %s", job_id)

    # Use DSV format with ASCII Unit Separator for robust parsing
    delimiter = "\x1f"  # ASCII Unit Separator - designed for field separation
    qstat_result = subprocess.run(
        ["qstat", "-fx", "-F", "dsv", "-D", delimiter, job_id],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    # Check if qstat failed
    if qstat_result.returncode != 0:
        logging.debug("qstat failed for job %s", job_id)
        click.echo(qstat_result.stderr)
        click.echo(qstat_result.args)
        sys.exit(qstat_result.returncode)

    # Parse DSV output to find job_state
    output = qstat_result.stdout.strip()
    if not output:
        logging.error("Empty qstat output for job %s", job_id)
        return "C"  # Assume completed if no output


def wait_for_job_exit_status(job_id, initial_wait=5, poll_interval=5, max_attempts=12):
    """
    Wait for PBS to clean up the job and get its exit status.

    Args:
        job_id: The PBS job ID
        initial_wait: Initial wait time for PBS cleanup (seconds)
        poll_interval: Interval between polling attempts (seconds)
        max_attempts: Maximum number of polling attempts

    Returns:
        int: The job's exit status, or 1 if unavailable
    """
    logging.debug("Waiting %d seconds for PBS cleanup", initial_wait)
    time.sleep(initial_wait)

    for attempt in range(max_attempts):
        exit_status = job_exit_status(job_id)
        if exit_status is not None:
            logging.info("Job %s exit status: %d", job_id, exit_status)
            return exit_status

        if attempt < max_attempts - 1:  # Don't sleep on the last attempt
            logging.debug(
                "Exit status not available yet, waiting %d seconds (attempt %d/%d)",
                poll_interval,
                attempt + 1,
                max_attempts,
            )
            time.sleep(poll_interval)

    logging.warning(
        "Could not get exit status for job %s after %d attempts", job_id, max_attempts
    )
    return 1  # Return failure exit code if we can't determine the real status


def job_exit_status(job_id):
    """
    Get the exit status of a completed job.

    Returns:
        int: The job's exit status, or None if not available
    """
    job_data = get_job_resource_data(job_id)
    if job_data and "exit_status" in job_data:
        return job_data["exit_status"]
    return None


def get_job_resource_data(job_id):
    """
    Get comprehensive resource data for a completed job from qstat -fx.

    Args:
        job_id: The PBS job ID

    Returns:
        dict: Complete job resource information, or None if not available

    The returned dictionary contains:
        - exit_status: Job exit code
        - resources_requested: Dict of requested resources
        - resources_used: Dict of actually used resources
        - execution: Dict of execution environment info
        - timing: Dict of job timing information
        - metadata: Dict of job metadata (name, owner, etc.)
    """
    logging.debug("Getting resource data for job %s", job_id)

    # Use standard qstat -fx format for full resource information
    qstat_result = subprocess.run(
        ["qstat", "-fx", job_id],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if qstat_result.returncode != 0:
        logging.debug("qstat failed for job %s: %s", job_id, qstat_result.stderr)
        return None

    output = qstat_result.stdout.strip()
    if not output:
        logging.debug("Empty qstat output for job %s", job_id)
        return None

    return parse_qstat_fx_output(output)


def parse_qstat_fx_output(qstat_output):
    """
    Parse qstat -fx output to extract comprehensive job resource information.

    Args:
        qstat_output: Raw output from qstat -fx command

    Returns:
        dict: Parsed job resource data
    """
    if not qstat_output:
        return None

    data = {
        "exit_status": None,
        "resources_requested": {},
        "resources_used": {},
        "execution": {},
        "timing": {},
        "metadata": {},
    }

    # Parse line by line
    for line in qstat_output.split("\n"):
        line = line.strip()
        if not line or "=" not in line:
            continue

        # Split on first '=' to handle values with '=' in them
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()

        # Parse different categories of information
        if key == "Exit_status":
            try:
                data["exit_status"] = int(value)
            except ValueError:
                logging.warning("Could not parse exit status: %s", value)

        elif key.startswith("Resource_List."):
            # Requested resources
            resource_key = key[14:]  # Remove 'Resource_List.' prefix
            data["resources_requested"][resource_key] = value

        elif key.startswith("resources_used."):
            # Used resources
            resource_key = key[15:]  # Remove 'resources_used.' prefix
            data["resources_used"][resource_key] = value

        elif key in [
            "exec_host",
            "exec_vnode",
            "queue",
            "session_id",
            "project",
            "group_list",
        ]:
            # Execution environment
            data["execution"][key] = value

        elif key in ["qtime", "stime", "etime", "mtime", "ctime", "obittime"]:
            # Timing information
            data["timing"][key] = value

        elif key in ["Job_Id", "Job_Name", "Job_Owner", "job_state"]:
            # Job metadata
            data["metadata"][key] = value

    # Post-process resource data for easier consumption
    data = _enhance_resource_data(data)

    return data


def _enhance_resource_data(data):
    """
    Enhance parsed resource data with calculated fields and conversions.

    Args:
        data: Raw parsed data from parse_qstat_fx_output

    Returns:
        dict: Enhanced data with calculated efficiency metrics
    """
    if not data:
        return data

    # Convert resource values to standardized formats
    requested = data["resources_requested"]
    used = data["resources_used"]

    # Parse timing information
    timing = data["timing"]
    for time_key in list(timing.keys()):  # Convert to list to avoid iteration issues
        if timing[time_key]:
            parsed_time = parse_timestamp(timing[time_key])
            if parsed_time:
                timing[f"{time_key}_parsed"] = parsed_time.isoformat()

    # Calculate derived timing metrics
    if "qtime_parsed" in timing and "stime_parsed" in timing:
        try:
            from datetime import datetime

            qtime = datetime.fromisoformat(timing["qtime_parsed"])
            stime = datetime.fromisoformat(timing["stime_parsed"])
            timing["queue_wait_seconds"] = int((stime - qtime).total_seconds())
        except Exception:
            pass

    if "stime_parsed" in timing and "mtime_parsed" in timing:
        try:
            from datetime import datetime

            stime = datetime.fromisoformat(timing["stime_parsed"])
            mtime = datetime.fromisoformat(timing["mtime_parsed"])
            timing["execution_seconds"] = int((mtime - stime).total_seconds())
        except Exception:
            pass

    # Parse execution environment
    execution = data["execution"]
    if "exec_host" in execution:
        execution["exec_host_parsed"] = parse_exec_host(execution["exec_host"])
    if "exec_vnode" in execution:
        execution["exec_vnode_parsed"] = parse_exec_vnode(execution["exec_vnode"])

    # Calculate resource efficiency metrics
    efficiency = {}

    # Memory efficiency
    if "mem" in requested and "mem" in used:
        try:
            mem_requested_bytes = size_to_bytes(requested["mem"])
            mem_used_bytes = size_to_bytes(used["mem"])
            efficiency["memory_efficiency"] = calculate_efficiency(
                mem_used_bytes, mem_requested_bytes
            )
            efficiency["memory_requested_human"] = bytes_to_human(mem_requested_bytes)
            efficiency["memory_used_human"] = bytes_to_human(mem_used_bytes)
        except Exception as e:
            logging.warning("Could not calculate memory efficiency: %s", e)

    # Time efficiency
    if "walltime" in requested and "walltime" in used:
        try:
            time_requested_seconds = time_to_seconds(requested["walltime"])
            time_used_seconds = time_to_seconds(used["walltime"])
            efficiency["time_efficiency"] = calculate_efficiency(
                time_used_seconds, time_requested_seconds
            )
        except Exception as e:
            logging.warning("Could not calculate time efficiency: %s", e)

    # CPU efficiency (from cpupercent if available)
    if "cpupercent" in used:
        try:
            efficiency["cpu_efficiency"] = float(used["cpupercent"])
        except Exception:
            pass

    # Jobfs efficiency
    if "jobfs" in requested and "jobfs" in used:
        try:
            jobfs_requested_bytes = size_to_bytes(requested["jobfs"])
            jobfs_used_bytes = size_to_bytes(used["jobfs"])
            efficiency["jobfs_efficiency"] = calculate_efficiency(
                jobfs_used_bytes, jobfs_requested_bytes
            )
            efficiency["jobfs_requested_human"] = bytes_to_human(jobfs_requested_bytes)
            efficiency["jobfs_used_human"] = bytes_to_human(jobfs_used_bytes)
        except Exception as e:
            logging.warning("Could not calculate jobfs efficiency: %s", e)

    data["efficiency"] = efficiency

    return data


def _format_completion_message(job_id, exit_status, status_dir):
    """
    Format a completion message with appropriate emoji and hook status information.

    Args:
        job_id: The PBS job ID
        exit_status: The final exit status of the job
        status_dir: Path to the status directory containing hook exit codes

    Returns:
        str: Formatted completion message
    """
    import pathlib

    # Choose emoji based on exit status
    if exit_status == 0:
        emoji = "✅"
        status_text = "Job completed successfully"
    else:
        emoji = "❌"
        status_text = f"Job failed with exit code {exit_status}"

    # Base message
    message = f"{emoji} {status_text}"

    # Check for hook status information
    hook_info = []

    # Check pre-command status
    pre_exit_file = pathlib.Path(status_dir) / "pre_exit_code"
    if pre_exit_file.exists():
        try:
            with open(pre_exit_file, "r", encoding="utf-8") as f:
                pre_exit_code = int(f.readline().strip())
                if pre_exit_code == 0:
                    hook_info.append("pre ✅")
                else:
                    hook_info.append(f"pre ❌({pre_exit_code})")
        except (OSError, ValueError):
            hook_info.append("pre ❓")

    # Check main command status (always present)
    main_exit_file = pathlib.Path(status_dir) / "main_exit_code"
    if main_exit_file.exists():
        try:
            with open(main_exit_file, "r", encoding="utf-8") as f:
                main_exit_code = int(f.readline().strip())
                if main_exit_code == 0:
                    hook_info.append("main ✅")
                else:
                    hook_info.append(f"main ❌({main_exit_code})")
        except (OSError, ValueError):
            hook_info.append("main ❓")

    # Check post-command status
    post_exit_file = pathlib.Path(status_dir) / "post_exit_code"
    if post_exit_file.exists():
        try:
            with open(post_exit_file, "r", encoding="utf-8") as f:
                post_exit_code = int(f.readline().strip())
                if post_exit_code == 0:
                    hook_info.append("post ✅")
                else:
                    hook_info.append(f"post ❌({post_exit_code})")
        except (OSError, ValueError):
            hook_info.append("post ❓")

    # Add hook information if any hooks were specified
    if hook_info:
        message += f" ({', '.join(hook_info)})"

    return message


def collect_job_resources_after_completion(job_id, out_file):
    """
    Collect detailed resource usage after job completion by parsing the joblog.
    Waits for PBS to finalize the joblog with resource information.
    """

    try:
        # Derive joblog path from out_file path
        joblog_path = str(out_file).replace(".out", ".log")

        logging.debug("Waiting for PBS to finalize joblog for job %s", job_id)

        # Wait 60 seconds for PBS to write final resource information to joblog
        time.sleep(60)

        # Try to parse the joblog with retry logic
        max_retries = 3
        retry_delay = 60  # Wait 60 seconds between retries

        for attempt in range(max_retries):
            logging.debug(
                "Attempting to collect resource usage from joblog for job %s "
                "(attempt %d/%d)",
                job_id,
                attempt + 1,
                max_retries,
            )

            try:
                from ..resources import parse_joblog_resources, resource_tracker

                # Parse resource data from joblog
                resource_data = parse_joblog_resources(joblog_path)

                if resource_data:
                    # Update the existing job record with resource data
                    success = resource_tracker.update_job_resources(
                        job_id, resource_data
                    )
                    if success:
                        logging.debug(
                            f"Successfully updated resource data for job {job_id}"
                        )
                        return

                    logging.debug(f"Failed to update resource data for job {job_id}")
                else:
                    logging.debug(
                        f"Could not parse resource data from joblog: {joblog_path}"
                    )

            except Exception as e:
                logging.debug(
                    f"Error collecting resources for job {job_id} (attempt {attempt + 1}): {e}"
                )

            # If not the last attempt, wait before retrying
            if attempt < max_retries - 1:
                logging.debug(
                    f"Retrying resource collection in {retry_delay} seconds..."
                )
                time.sleep(retry_delay)

        logging.debug(
            f"Failed to collect resource data for job {job_id} after {max_retries} attempts"
        )

    except Exception as e:
        logging.debug("Resource collection failed for job %s: %s", job_id, e)


def start_background_resource_collection(
    job_id, joblog_path, foreground=False, verbose=0
):
    """
    Start resource collection in background subprocess using fork.

    This creates a truly independent subprocess that survives after the parent
    exits, avoiding threading complexities and allowing the parent to return
    to the shell immediately.

    Logs progress to ~/.config/qxub/resource_collection.log for monitoring.
    """
    import logging
    import os
    import time
    from datetime import datetime

    def log_message(message):
        """Log message to dedicated resource collection log file."""
        try:
            log_dir = os.path.expanduser("~/.config/qxub")
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(log_dir, "resource_collection.log")

            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            with open(log_file, "a") as f:
                f.write(f"[{timestamp}] {message}\n")
        except Exception as e:
            logging.debug(f"Error writing to resource collection log: {e}")

    def collect_resources_background():
        """Actual resource collection logic."""
        try:
            log_message(f"Starting resource collection for job {job_id}")
            log_message(f"Joblog path: {joblog_path}")

            from ..resources import parse_joblog_resources, resource_tracker

            # Mark job as pending resource collection
            pending_data = {"status": "pending_resources"}
            resource_tracker.update_job_resources(job_id, pending_data)
            log_message("Marked job as pending resource collection")

            # Wait for PBS to finalize the joblog
            log_message("Waiting 60 seconds for PBS to finalize joblog...")
            time.sleep(60)

            # Check if file exists and parse it
            if os.path.exists(joblog_path):
                log_message("Joblog file exists, parsing...")
                resource_data = parse_joblog_resources(joblog_path)

                if resource_data:
                    log_message(f"Successfully parsed resource data")
                    success = resource_tracker.update_job_resources(
                        job_id, resource_data
                    )
                    if success:
                        log_message("Successfully updated database with resource data")
                    else:
                        log_message("ERROR: Failed to update database")
                else:
                    log_message("ERROR: Failed to parse resource data from joblog")
                    # Show file info for debugging
                    try:
                        with open(joblog_path, "r") as f:
                            lines = f.readlines()
                            log_message(f"Joblog has {len(lines)} lines")
                            if len(lines) > 10:
                                log_message(f"Last 10 lines: {lines[-10:]}")
                    except Exception as e:
                        log_message(f"Error reading joblog for debugging: {e}")
            else:
                log_message(f"ERROR: Joblog file does not exist: {joblog_path}")
                # List available files for debugging
                try:
                    log_dir = os.path.dirname(joblog_path)
                    if os.path.exists(log_dir):
                        files = [f for f in os.listdir(log_dir) if f.endswith(".log")]
                        log_message(f"Available .log files in {log_dir}: {files}")
                except Exception as e:
                    log_message(f"Error listing directory: {e}")

        except Exception as e:
            log_message(f"FATAL ERROR in resource collection: {e}")
            import traceback

            log_message(f"Traceback: {traceback.format_exc()}")

    # If foreground requested, run collection inline for easier debugging
    if foreground:
        # Also emit to stdout so user sees verbose progress
        print(
            f"🔍 DEBUG: Running resource collection for job {job_id} in foreground (verbose debug)"
        )
        collect_resources_background()
        return None

    # Fork a child process for background collection
    try:
        pid = os.fork()
    except OSError as e:
        logging.warning(f"Failed to fork for background resource collection: {e}")
        # Fall back to no collection rather than failing the whole job
        return None

    if pid > 0:
        # Parent process - return immediately
        if verbose >= 2:
            print(f"🔍 DEBUG: Started background resource collection (PID {pid})")
            print(
                f"🔍 DEBUG: Monitor progress: tail -f ~/.config/qxub/resource_collection.log"
            )
        return None

    # Child process continues here
    # Close standard file descriptors to fully detach
    try:
        sys.stdin.close()
        sys.stdout.close()
        sys.stderr.close()

        # Redirect stdin/stdout/stderr to /dev/null
        devnull = os.open(os.devnull, os.O_RDWR)
        os.dup2(devnull, 0)  # stdin
        os.dup2(devnull, 1)  # stdout
        os.dup2(devnull, 2)  # stderr
        if devnull > 2:
            os.close(devnull)
    except Exception as e:
        # If descriptor closing fails, log but continue
        log_message(f"Warning: Failed to close descriptors: {e}")

    # Do the actual resource collection work
    collect_resources_background()

    # Child process exits when done (don't return to caller)
    os._exit(0)

    return thread


def monitor_job_single_thread(
    job_id, out_file, err_file, quiet=False, joblog_file=None, verbose=0
):
    """
    Single-thread job monitoring with proper spinner integration using status files.

    Flow:
    1. "Job constructed" progress message (handled by caller)
    2. Spinner while waiting for submission (handled by caller)
    3. "Job submitted" (handled by caller)
    4. Spinner while waiting for job to start
    5. "Job started"
    6. Job output streaming
    7. "Job completed" with status

    Returns:
        int: Job exit status
    """
    from pathlib import Path

    # Determine status directory based on output file location
    status_dir = Path(out_file).parent / "status"

    # Status files to monitor (using job ID directly for uniqueness)
    job_started_file = status_dir / f"job_started_{job_id}"
    final_exit_code_file = status_dir / f"final_exit_code_{job_id}"

    # Track job state
    job_started_notified = False

    # 4. Spinner while waiting for job to start - check status files every 0.5 seconds
    spinner = SimpleSpinner(quiet=quiet)
    try:
        while True:
            # Show spinner frame
            spinner.show_next_frame()

            # Check for job started (more frequent polling for responsiveness)
            if not job_started_notified and job_started_file.exists():
                spinner.clear()  # Clear spinner before showing message
                if not quiet:
                    print_status("🚀 Job started", final=True)
                    job_started_notified = True
                break

            # Check if job finished without running (rare but possible)
            if final_exit_code_file.exists():
                spinner.clear()  # Clear spinner before showing message
                logging.info("Job %s completed before starting", job_id)
                exit_status = read_exit_status_from_file(final_exit_code_file)

                # Completion message
                if not quiet:
                    completion_msg = _format_completion_message(
                        job_id, exit_status, status_dir
                    )
                    print_status(completion_msg, final=True)

                return exit_status

            time.sleep(0.5)  # Check every 500ms for good responsiveness
    finally:
        spinner.clear()  # Always clear spinner when exiting

    # 6. Show output streaming message
    if not quiet:
        print_status("📡 Streaming job output", final=True)

    # 7. Stream job output
    exit_status = stream_job_output_with_status_files(
        job_id, out_file, err_file, final_exit_code_file
    )

    # 8. Completion message
    if not quiet:
        completion_msg = _format_completion_message(job_id, exit_status, status_dir)
        print_status(completion_msg, final=True)

    # 9. Start background resource collection and return control immediately
    if not quiet and verbose < 2:
        print_status(
            "🎉 Job completed! Resource collection starting in background...",
            final=True,
        )

    # Use provided joblog path or derive from out_file path as fallback
    if joblog_file:
        joblog_path = str(joblog_file)
        if verbose >= 2:
            print(f"🔍 DEBUG: Using provided joblog path: {joblog_path}")
    else:
        # Fallback: derive joblog path from out_file path (change .out to .log)
        joblog_path = str(out_file).replace(".out", ".log")
        if verbose >= 2:
            print(f"🔍 DEBUG: Derived joblog path from out_file: {joblog_path}")

    # If user requested very verbose debug (e.g. -vvv) run collection in foreground
    run_foreground = bool(verbose and int(verbose) >= 3)
    if run_foreground and verbose >= 2:
        print(f"🔍 DEBUG: Foreground resource collection enabled (verbose={verbose})")

    thread = start_background_resource_collection(
        job_id, joblog_path, foreground=run_foreground, verbose=verbose
    )

    return exit_status


def read_exit_status_from_file(exit_code_file):
    """Read exit status from a status file."""
    try:
        with open(str(exit_code_file), "r", encoding="utf-8") as f:
            content = f.read().strip()
            # Status file format: first line is exit code, second line is timestamp
            first_line = content.split("\n")[0]
            return int(first_line)
    except (OSError, IOError, ValueError) as e:
        logging.warning(f"Could not read exit status from {exit_code_file}: {e}")
        return 1  # Default to failure


def stream_job_output_with_status_files(
    job_id, out_file, err_file, final_exit_code_file
):
    """
    Stream job output files until job completion using status files for detection.

    Returns:
        int: Job exit status
    """
    import os

    # Keep track of file positions to avoid re-reading
    out_pos = 0
    err_pos = 0

    # Monitor job completion via status file
    while True:
        # Stream new output from STDOUT
        if os.path.exists(out_file):
            try:
                with open(out_file, "r", encoding="utf-8") as f:
                    f.seek(out_pos)
                    for line in f:
                        print(line.rstrip(), file=sys.stdout, flush=True)
                    out_pos = f.tell()
            except (OSError, IOError):
                pass

        # Stream new output from STDERR
        if os.path.exists(err_file):
            try:
                with open(err_file, "r", encoding="utf-8") as f:
                    f.seek(err_pos)
                    for line in f:
                        print(line.rstrip(), file=sys.stderr, flush=True)
                    err_pos = f.tell()
            except (OSError, IOError):
                pass

        # Check if job finished using status file
        if final_exit_code_file.exists():
            # Read any final output
            if os.path.exists(out_file):
                try:
                    with open(out_file, "r", encoding="utf-8") as f:
                        f.seek(out_pos)
                        for line in f:
                            print(line.rstrip(), file=sys.stdout, flush=True)
                except (OSError, IOError):
                    pass

            if os.path.exists(err_file):
                try:
                    with open(err_file, "r", encoding="utf-8") as f:
                        f.seek(err_pos)
                        for line in f:
                            print(line.rstrip(), file=sys.stderr, flush=True)
                except (OSError, IOError):
                    pass

            # Get final exit status from status file
            exit_status = read_exit_status_from_file(final_exit_code_file)
            return exit_status

        time.sleep(0.2)  # Check every 200ms for good output responsiveness


def stream_job_output(job_id, out_file, err_file):
    """
    Stream job output files until job completion.

    Returns:
        int: Job exit status
    """
    import os

    # Keep track of file positions to avoid re-reading
    out_pos = 0
    err_pos = 0

    # Monitor job status and stream output
    while True:
        status = job_status(job_id)

        # Stream new output from STDOUT
        if os.path.exists(out_file):
            try:
                with open(out_file, "r", encoding="utf-8") as f:
                    f.seek(out_pos)
                    for line in f:
                        print(line.rstrip(), file=sys.stdout, flush=True)
                    out_pos = f.tell()
            except (OSError, IOError):
                pass

        # Stream new output from STDERR
        if os.path.exists(err_file):
            try:
                with open(err_file, "r", encoding="utf-8") as f:
                    f.seek(err_pos)
                    for line in f:
                        print(line.rstrip(), file=sys.stderr, flush=True)
                    err_pos = f.tell()
            except (OSError, IOError):
                pass

        # Check if job finished
        if status in ["F", "H", "C"]:
            # Read any final output
            if os.path.exists(out_file):
                try:
                    with open(out_file, "r", encoding="utf-8") as f:
                        f.seek(out_pos)
                        for line in f:
                            print(line.rstrip(), file=sys.stdout, flush=True)
                except (OSError, IOError):
                    pass

            if os.path.exists(err_file):
                try:
                    with open(err_file, "r", encoding="utf-8") as f:
                        f.seek(err_pos)
                        for line in f:
                            print(line.rstrip(), file=sys.stderr, flush=True)
                except (OSError, IOError):
                    pass

            # Get final exit status
            exit_status = wait_for_job_exit_status(job_id)
            return exit_status

        time.sleep(1)  # Check every 1 second
