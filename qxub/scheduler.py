"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

import logging
import os
import pathlib
import shlex
import signal
import subprocess
import sys
import threading
import time

import click
import tailer

from .resource_parser import (
    bytes_to_human,
    calculate_efficiency,
    parse_exec_host,
    parse_exec_vnode,
    parse_timestamp,
    size_to_bytes,
    time_to_seconds,
)


class OutputCoordinator:
    """Coordinates between spinner, output threads, and monitoring threads."""

    def __init__(self):
        self.output_started = threading.Event()
        self.spinner_cleared = threading.Event()
        self.job_completed = threading.Event()
        self.shutdown_requested = threading.Event()
        self.eof_detected = threading.Event()
        self.submission_complete = (
            threading.Event()
        )  # New event for submission completion
        # New events for job status changes
        self.job_running = threading.Event()  # Job status changed to Running
        self.job_finished = threading.Event()  # Job status changed to Finished/Held
        self.job_exit_status = None  # Store the job's exit status

    def signal_output_started(self):
        """Called by tail threads when they start producing output."""
        self.output_started.set()

    def wait_for_output_started(self, timeout=None):
        """Called by spinner to wait for output to start."""
        return self.output_started.wait(timeout)

    def signal_spinner_cleared(self):
        """Called by spinner when it has cleared itself."""
        self.spinner_cleared.set()

    def signal_job_completed(self):
        """Called by monitor when job completes."""
        self.job_completed.set()

    def signal_shutdown(self):
        """Called when shutdown is requested (e.g., Ctrl-C)."""
        self.shutdown_requested.set()

    def signal_eof(self):
        """Called by tail threads when they reach EOF."""
        self.eof_detected.set()

    def should_shutdown(self):
        """Check if threads should shutdown."""
        return (
            self.shutdown_requested.is_set()
            or self.job_completed.is_set()
            or self.eof_detected.is_set()
        )

    def wait_for_spinner_cleared(self, timeout=None):
        """Called by tail threads to wait for spinner to clear."""
        return self.spinner_cleared.wait(timeout)

    def signal_submission_complete(self):
        """Called when job submission messages are complete."""
        self.submission_complete.set()

    def wait_for_submission_complete(self, timeout=None):
        """Called by monitor thread to wait for submission to complete."""
        return self.submission_complete.wait(timeout)

    def signal_job_running(self):
        """Called when job status changes to Running."""
        self.job_running.set()

    def signal_job_finished(self):
        """Called when job status changes to Finished or Held."""
        self.job_finished.set()

    def should_stop_spinner(self):
        """Check if spinner should stop due to any relevant event."""
        return (
            self.output_started.is_set()
            or self.job_running.is_set()
            or self.job_finished.is_set()
            or self.shutdown_requested.is_set()
            or self.job_completed.is_set()
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
        # Fallback to /dev/null if /dev/tty is not available (non-interactive context)
        # This ensures progress messages don't interfere with stdout redirection
        with open("/dev/null", "w") as devnull:
            print(f"{message}", file=devnull, flush=True)


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
        if not line or not "=" in line:
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


def job_status(job_id):
    """
    Check the current status of a completed job_status function.

    Returns:
        string: Job state (R, Q, H, F, C, etc.)
    """
    logging.debug("Getting job status for job %s", job_id)

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
        return "C"  # Assume completed if qstat fails

    # Parse DSV output to find job_state
    output = qstat_result.stdout.strip()
    if not output:
        logging.error("Empty qstat output for job %s", job_id)
        return "C"  # Assume completed if no output

    # Split on delimiter and look for job_state field
    fields = output.split(delimiter)
    for field in fields:
        if field.startswith("job_state="):
            status = field.split("=", 1)[1]
            logging.debug("Job status: %s", status)
            return status

    # Fallback: if job_state not found, try simple qstat
    logging.warning("job_state field not found in DSV output, trying simple qstat")
    try:
        fallback_result = subprocess.run(
            ["qstat", job_id],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if fallback_result.returncode == 0 and fallback_result.stdout:
            # Parse simple qstat output
            lines = fallback_result.stdout.strip().split("\n")
            if len(lines) >= 2:  # Header + job line
                job_line = lines[-1].split()
                if len(job_line) >= 8:
                    status = job_line[-1]  # Last column is status
                    logging.debug("Job status from simple qstat: %s", status)
                    return status

        # Final fallback - assume completed
        logging.warning("All qstat methods failed, assuming job completed")
        return "C"
    except Exception as e:
        logging.error("Exception in job_status fallback: %s", e)
        return "C"


def monitor_qstat(job_id, quiet=False, coordinator=None, success_msg=None):
    """
    Monitors for job completion by checking job status periodically.
    Now uses event-driven spinner control.
    """
    logging.debug("Starting job monitoring with initial check")

    # Track job state
    job_started_notified = False
    last_status = None

    # Wait for job submission to complete before starting
    if coordinator:
        coordinator.wait_for_submission_complete()

    # Start event-driven spinner - it will stop automatically on status changes
    logging.debug("Starting event-driven spinner and monitoring loop")

    # Use proper context manager for spinner
    with JobSpinner("", show_message=False, quiet=quiet, coordinator=coordinator):
        while True:
            # Check if we should shutdown early
            if coordinator and coordinator.should_shutdown():
                logging.info("Monitor shutdown requested")
                return

            status = job_status(job_id)

            # Signal status changes to coordinator for spinner control
            if status != last_status:
                logging.debug(f"Job status changed from {last_status} to {status}")

                if status == "R":  # Job started running
                    if coordinator:
                        coordinator.signal_job_running()
                    if not job_started_notified and not quiet:
                        print_status("🚀 Job started running", final=True)
                        job_started_notified = True

                elif status in ["F", "H"]:  # Job finished or held
                    if coordinator:
                        coordinator.signal_job_finished()
                        coordinator.signal_job_completed()

                    # Wait for PBS cleanup and get exit status
                    logging.info("Job %s completed with status %s", job_id, status)
                    logging.info("Waiting for PBS cleanup and getting exit status...")
                    exit_status = wait_for_job_exit_status(job_id)
                    if coordinator:
                        coordinator.job_exit_status = exit_status
                    return exit_status

                last_status = status

            logging.debug("Job %s status: %s", job_id, status)

            # Sleep in smaller intervals so we can respond to shutdown requests
            for _ in range(6):  # 6 * 5 = 30 seconds total
                if coordinator and coordinator.should_shutdown():
                    logging.info("Monitor shutdown during sleep")
                    return 0  # Return success if interrupted
                time.sleep(5)


def monitor_files(
    job_id, out_file, err_file, quiet=False, coordinator=None, success_msg=None
):
    """
    Monitors for job completion by watching status files instead of polling qstat.
    This is much more efficient as it uses filesystem notifications rather than
    30-second qstat polling that creates load on the PBS scheduler.
    """
    import os
    import pathlib

    logging.debug("Starting file-based job monitoring")

    # Determine status directory based on output file location
    status_dir = pathlib.Path(out_file).parent / "status"

    # Status files to monitor (using job ID directly for uniqueness)
    job_started_file = status_dir / f"job_started_{job_id}"
    main_started_file = status_dir / f"main_started_{job_id}"
    final_exit_code_file = status_dir / f"final_exit_code_{job_id}"

    # Track job state
    job_started_notified = False
    job_running_notified = False

    # Wait for job submission to complete before starting
    if coordinator:
        coordinator.wait_for_submission_complete()

    logging.debug("Starting file-based monitoring loop")

    # Use proper context manager for spinner
    with JobSpinner("", show_message=False, quiet=quiet, coordinator=coordinator):
        while True:
            # Check if we should shutdown early
            if coordinator and coordinator.should_shutdown():
                logging.info("File monitor shutdown requested")
                return

            # Check for job started
            if not job_started_notified and job_started_file.exists():
                if coordinator:
                    coordinator.signal_job_running()
                if not quiet:
                    print_status("🚀 Job started running", final=True)
                job_started_notified = True
                logging.debug("Job started detected via status file")

            # Check for main command started (indicates job is actively running)
            if not job_running_notified and main_started_file.exists():
                job_running_notified = True
                logging.debug("Main command started detected via status file")

            # Check for job completion
            if final_exit_code_file.exists():
                if coordinator:
                    coordinator.signal_job_finished()
                    coordinator.signal_job_completed()

                # Start background resource collection
                start_background_resource_collection(job_id, final_exit_code_file)

                # Read exit status from file
                try:
                    with open(final_exit_code_file, "r") as f:
                        lines = f.read().strip().split("\n")
                        exit_status = (
                            int(lines[0]) if lines and lines[0].isdigit() else 0
                        )
                except (OSError, ValueError, IndexError) as e:
                    logging.warning(
                        f"Could not read exit status from {final_exit_code_file}: {e}"
                    )
                    exit_status = 0

                logging.info(
                    "Job %s completed with exit status %s (from status file)",
                    job_id,
                    exit_status,
                )

                if coordinator:
                    coordinator.job_exit_status = exit_status
                    # Signal shutdown to stop tail threads
                    coordinator.signal_shutdown()

                return exit_status

            # Sleep much shorter than qstat polling - we want responsive file detection
            # But not too aggressive to avoid excessive CPU usage
            for _ in range(2):  # 2 * 1 = 2 seconds total (vs 30s for qstat)
                if coordinator and coordinator.should_shutdown():
                    logging.info("File monitor shutdown during sleep")
                    return 0  # Return success if interrupted
                time.sleep(1)


def start_background_resource_collection(job_id, final_exit_code_file):
    """
    Start background resource collection that waits 60 seconds after job completion
    and then collects detailed resource usage via qstat -f. This runs independently
    and doesn't block the main monitoring loop.
    """

    def collect_resources():
        try:
            # Wait for job completion signal
            while not final_exit_code_file.exists():
                time.sleep(2)  # Check every 2 seconds for completion

            logging.debug(
                f"Job {job_id} completed, starting 60s delay for resource collection"
            )

            # Wait 60 seconds for PBS cleanup and detailed stats to be available
            time.sleep(60)

            logging.debug(f"Collecting detailed resource usage for job {job_id}")

            # Try to collect detailed job statistics
            try:
                # Import here to avoid circular imports
                from .resource_parser import parse_qstat_output
                from .resource_tracker import store_resource_usage

                # Get detailed job info via qstat -f
                result = subprocess.run(
                    ["qstat", "-f", job_id], capture_output=True, text=True, timeout=30
                )

                if result.returncode == 0:
                    # Parse and store resource usage
                    job_info = parse_qstat_output(result.stdout)
                    if job_info:
                        store_resource_usage(job_id, job_info)
                        logging.info(
                            f"Successfully collected and stored resource usage for job {job_id}"
                        )
                    else:
                        logging.warning(
                            f"Could not parse qstat output for job {job_id}"
                        )
                else:
                    logging.warning(
                        f"qstat -f failed for job {job_id}: {result.stderr}"
                    )

            except ImportError:
                logging.debug(
                    "Resource tracking modules not available, skipping collection"
                )
            except subprocess.TimeoutExpired:
                logging.warning(f"qstat -f timed out for job {job_id}")
            except Exception as e:
                logging.warning(f"Error collecting resources for job {job_id}: {e}")

        except Exception as e:
            logging.error(
                f"Background resource collection failed for job {job_id}: {e}"
            )

    # Start background thread for resource collection
    resource_thread = threading.Thread(target=collect_resources, daemon=True)
    resource_thread.start()
    logging.debug(f"Started background resource collection for job {job_id}")
    return resource_thread


def tail(log_file, destination, coordinator=None):
    """
    Tails the given log_file until either EOF or job completion.
    Output is directed to the specified destination (STDOUT or STDERR)
    """
    if not destination in ["STDOUT", "STDERR"]:
        click.echo("Unknown destination for redirection")
        sys.exit(2)
    is_err = destination == "STDERR"
    output_started = False

    try:
        with open(log_file, "r", encoding="utf-8") as f:
            # First, read any existing content in the file
            while True:
                line = f.readline()
                if not line:  # No more existing content
                    break

                # Check if we should shutdown
                if coordinator and coordinator.should_shutdown():
                    logging.debug(f"Tail {destination} shutdown requested")
                    return

                # Signal that output has started on first line
                if not output_started and coordinator:
                    coordinator.signal_output_started()
                    # Clear any leftover spinner characters before streaming output
                    import os

                    is_remote_ssh = any(
                        [
                            os.environ.get("SSH_CLIENT"),
                            os.environ.get("SSH_CONNECTION"),
                            os.environ.get("SSH_TTY"),
                        ]
                    )

                    if not is_remote_ssh:  # Only clear spinner in local contexts
                        try:
                            with open("/dev/tty", "w") as tty:
                                print("\r", end="", flush=True, file=tty)
                        except (OSError, IOError):
                            pass  # Ignore if /dev/tty is not available
                    coordinator.signal_spinner_cleared()
                    output_started = True

                print(
                    line.rstrip(), file=sys.stderr if is_err else sys.stdout, flush=True
                )

            # Now follow for any new content (for longer-running jobs)
            for line in tailer.follow(f):
                # Check if we should shutdown
                if coordinator and coordinator.should_shutdown():
                    logging.debug(f"Tail {destination} shutdown requested")
                    break

                # Signal that output has started on first line (if we didn't already)
                if not output_started and coordinator:
                    coordinator.signal_output_started()
                    # Clear any leftover spinner characters before streaming output
                    import os

                    is_remote_ssh = any(
                        [
                            os.environ.get("SSH_CLIENT"),
                            os.environ.get("SSH_CONNECTION"),
                            os.environ.get("SSH_TTY"),
                        ]
                    )

                    if not is_remote_ssh:  # Only clear spinner in local contexts
                        try:
                            with open("/dev/tty", "w") as tty:
                                print("\r", end="", flush=True, file=tty)
                        except (OSError, IOError):
                            pass  # Ignore if /dev/tty is not available
                    coordinator.signal_spinner_cleared()
                    output_started = True

                print(
                    line.rstrip(), file=sys.stderr if is_err else sys.stdout, flush=True
                )

            # If we exit the loop normally, we've reached EOF
            if coordinator:
                coordinator.signal_eof()
                logging.debug(f"Tail {destination} reached EOF")

    except Exception as e:
        logging.error(f"Error in tail {destination}: {e}")
    finally:
        logging.debug(f"Tail {destination} thread exiting")


def start_job_monitoring(job_id, out_file, err_file, quiet=False, success_msg=None):
    """
    Start job monitoring threads and return coordinator for external signaling.

    Returns:
        tuple: (coordinator, monitor_function) where monitor_function() waits for completion
    """
    # Create coordinator for thread synchronization
    coordinator = OutputCoordinator()

    # Check if we should use file-based monitoring (default) or fall back to qstat
    # File-based monitoring is preferred for performance, but qstat is available as fallback
    use_file_monitoring = (
        os.environ.get("QXUB_USE_QSTAT_MONITORING", "false").lower() != "true"
    )

    if use_file_monitoring:
        logging.debug("Using file-based monitoring (preferred)")
        monitor_target = monitor_files
        monitor_args = (job_id, out_file, err_file, quiet, coordinator, success_msg)
    else:
        logging.debug("Using qstat polling monitoring (fallback mode)")
        monitor_target = monitor_qstat
        monitor_args = (job_id, quiet, coordinator, success_msg)

    # Create threads for job monitoring and log tailing
    monitor_thread = threading.Thread(target=monitor_target, args=monitor_args)
    out_thread = threading.Thread(
        target=tail, args=(out_file, "STDOUT", coordinator), daemon=True
    )
    err_thread = threading.Thread(
        target=tail, args=(err_file, "STDERR", coordinator), daemon=True
    )

    # Set up signal handler for Ctrl-C
    def signal_handler(signum, frame):
        logging.info("Interrupt received, shutting down threads...")
        # Clean up job first
        print(" " * 100, end="", flush=True)  # Clear current line
        print("🛑 Interrupted! Cleaning up job...")
        success = qdel(job_id, quiet=False)
        if success:
            print("✅ Job cleanup completed")
        else:
            print(
                f"⚠️  Job cleanup failed - you may need to manually run: qdel {job_id}"
            )
        # Then signal coordinator to shutdown threads
        coordinator.signal_shutdown()

    original_handler = signal.signal(signal.SIGINT, signal_handler)

    # Start all threads
    monitor_thread.start()
    out_thread.start()
    err_thread.start()

    def wait_for_completion():
        """Wait for job monitoring to complete and return exit status"""
        try:
            # Wait for job monitoring to complete or shutdown signal
            monitor_thread.join()

            # Wait for tail threads to finish after shutdown signal
            # Give them a reasonable timeout since they should exit quickly once signaled
            out_thread.join(timeout=5)
            err_thread.join(timeout=5)

            return coordinator.job_exit_status or 0
        finally:
            # Restore original signal handler
            signal.signal(signal.SIGINT, original_handler)

    return coordinator, wait_for_completion


def monitor_and_tail(
    job_id, out_file, err_file, quiet=False, success_msg=None, coordinator=None
):
    """
    Monitors the job status using qstat and tails the output (STDOUT and STDERR) logs
    until the job is finished. Stops both tailing and monitoring upon job completion.

    Args:
        job_id: The PBS job id to monitor.
        out_file: Path to the STDOUT log file to tail.
        err_file: Path to the STDERR log file to tail.
        quiet: Whether to suppress spinner output.
        success_msg: The success message to display with spinner.
        coordinator: Optional existing coordinator, creates new one if None

    Returns:
        int: The job's exit status
    """
    # Create coordinator for thread synchronization if not provided
    if coordinator is None:
        coordinator = OutputCoordinator()

    # Create threads for job monitoring and log tailing
    qstat_thread = threading.Thread(
        target=monitor_qstat, args=(job_id, quiet, coordinator, success_msg)
    )
    out_thread = threading.Thread(
        target=tail, args=(out_file, "STDOUT", coordinator), daemon=True
    )
    err_thread = threading.Thread(
        target=tail, args=(err_file, "STDERR", coordinator), daemon=True
    )

    # Set up signal handler for Ctrl-C
    def signal_handler(signum, frame):
        logging.info("Interrupt received, shutting down threads...")
        coordinator.signal_shutdown()

    original_handler = signal.signal(signal.SIGINT, signal_handler)

    try:
        # Start tail threads immediately to capture output
        out_thread.start()
        err_thread.start()

        # Start monitoring thread immediately (it will wait for submission signal)
        qstat_thread.start()

        # Wait for job monitoring to complete or shutdown signal
        qstat_thread.join()

        # Signal shutdown to remaining threads
        coordinator.signal_shutdown()

        # Give threads a moment to cleanup gracefully
        time.sleep(0.5)

        logging.info("Job monitoring completed")

        # Return the job's exit status
        return (
            coordinator.job_exit_status
            if coordinator.job_exit_status is not None
            else 0
        )

    except KeyboardInterrupt:
        logging.info("KeyboardInterrupt received, shutting down...")
        coordinator.signal_shutdown()
        qstat_thread.join(timeout=2)  # Give monitor thread time to cleanup
        return 130  # Standard exit code for SIGINT

    finally:
        # Restore original signal handler
        signal.signal(signal.SIGINT, original_handler)
