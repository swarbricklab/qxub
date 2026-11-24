"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

import logging
import sys
import time

import click

from .. import backends
from ..backends import SubmitError
from ..resources import time_to_seconds

logger = logging.getLogger(__name__)


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


def same_job(job_id: str) -> str:
    """A key that matches ``12345`` with ``12345.gadi-pbs``."""
    return backends.backend_for_job(job_id).same_job(job_id)


# Submission failures. The old name stays for callers that catch it.
QsubError = SubmitError


def qdel(job_id, quiet=False):
    """Cancel a job through the platform's scheduler backend.

    Args:
        job_id (str): The job ID to cancel
        quiet (bool): whether to suppress output

    Returns:
        bool: True if the scheduler accepted the cancellation
    """
    logger.debug("Deleting job %s", job_id)
    ok, reason = backends.backend_for_job(job_id).cancel(job_id)
    if ok:
        if not quiet:
            click.echo(f"🗑️  Job {job_id} deleted successfully")
        logger.info("Job %s deleted successfully", job_id)
        return True
    if not quiet:
        click.echo(f"Failed to delete job {job_id}: {reason}")
    logger.warning("Failed to delete job %s: %s", job_id, reason)
    return False


def get_default_status_dir():
    """Get the default status directory path for qxub jobs.

    Returns:
        Path to the status directory: ``{scratch}/qxub/status`` from the
        platform definition, or ``$XDG_STATE_HOME/qxub/status``.
    """
    from ..config.site import default_status_dir

    return default_status_dir()


def status_file_key(job_id, status_dir):
    """The name a job's status files use.

    Jobs name them after the qxub run ID (``QXUB_RUN_ID``) and leave a
    ``native_<job ID>`` file holding it, so a scheduler job ID can be looked
    up without the database. Jobs from older qxub used the job ID itself.
    """
    from pathlib import Path

    from ..queue.db import is_virtual_id

    if is_virtual_id(job_id):
        return job_id
    try:
        run_id = (Path(status_dir) / f"native_{job_id}").read_text().strip()
    except OSError:
        return job_id
    return run_id or job_id


def job_status_from_files(job_id, status_dir=None, log_dir=None):
    """
    Check the current status of a job using status files instead of qstat.

    This is faster and more reliable than polling qstat, especially when the
    PBS scheduler is under heavy load.

    Status detection priority:
    1. final_exit_code_{key} file - definitive completion status
    2. PBS job log file with Exit Status - job was killed by PBS
    3. the run file (runs/{key}.run) saying state=finished - e.g. a signal
    4. job_started_{key} file - job is running
    5. No files - job is queued or unknown

    ``key`` is the run ID (see :func:`status_file_key`).

    Args:
        job_id: PBS job ID (e.g., "12345.gadi-pbs") or qxub run ID
        status_dir: Path to status directory. If None, uses default.
        log_dir: Path to log directory (for PBS job logs). If None, uses parent of status_dir.

    Returns:
        Tuple of (status, exit_code) where:
        - status: "Q" (queued/unknown), "R" (running), "C" (completed), or "F" (failed)
        - exit_code: int if completed, None otherwise
    """
    import os
    from pathlib import Path

    if status_dir is None:
        status_dir = get_default_status_dir()
    else:
        status_dir = Path(status_dir)

    # Log directory is typically the parent of status directory
    if log_dir is None:
        log_dir = status_dir.parent
    else:
        log_dir = Path(log_dir)

    # Check for status files
    key = status_file_key(job_id, status_dir)
    job_started_file = status_dir / f"job_started_{key}"
    final_exit_code_file = status_dir / f"final_exit_code_{key}"

    # Check if job has completed via our exit code file
    if final_exit_code_file.exists():
        try:
            with open(final_exit_code_file, "r", encoding="utf-8") as f:
                exit_code = int(f.readline().strip())
                if exit_code == 0:
                    return "C", exit_code  # Completed successfully
                else:
                    return "F", exit_code  # Failed with non-zero exit
        except (ValueError, IOError) as e:
            logger.warning(f"Could not read exit code from {final_exit_code_file}: {e}")
            return "C", None  # Assume completed if file exists but unreadable

    # Check for PBS job log with exit status (appears after job ends).
    # This covers two cases:
    # 1. job_started_ exists: PBS killed the job before our cleanup code ran.
    # 2. job_started_ missing: job failed before qxub's startup code could write
    #    the file (e.g. conda activate failed, module load error, OOM at startup).
    #    PBS still writes the resource-usage block, so we can detect completion.
    exit_code = _check_pbs_job_log_for_exit(job_id, log_dir)
    if exit_code is not None:
        logger.debug(f"Found PBS exit status {exit_code} for job {job_id} in log file")
        return ("C" if exit_code == 0 else "F"), exit_code

    # The run file records every exit, a signal included (state=finished,
    # exit_code), so a job stopped by the scheduler or `qxub cancel` isn't
    # reported running when its job log is elsewhere.
    from ..queue import runs as run_files

    record = run_files.read_run(status_dir / "runs" / f"{key}{run_files.SUFFIX}")
    if record and record.get("state") == "finished":
        code = run_files.exit_code(record)
        if code is not None:
            return ("C" if code == 0 else "F"), code

    # Check if job has started (and PBS log not yet written — still running)
    if job_started_file.exists():
        return "R", None

    # No status files found - job is either queued or unknown
    return "Q", None


def read_exit_from_joblog_file(joblog_path):
    """
    Read the PBS exit status from a specific job log file.

    PBS appends a resource usage block to the .log file when the job ends,
    including when the job is killed (e.g. walltime exceeded).  The block
    contains a line like:
        Exit Status:        -29 (Job failed due to exceeding walltime)

    This function is intentionally cheap — it opens one known file, so it
    can be called in a tight polling loop without directory scanning.

    Args:
        joblog_path: Path to the PBS .log file for this job

    Returns:
        int: Exit code if the resource-usage block has been written, else None
    """
    import re

    try:
        with open(str(joblog_path), "r", encoding="utf-8") as f:
            content = f.read()
        exit_match = re.search(r"Exit Status:\s*(-?\d+)", content)
        if exit_match:
            return int(exit_match.group(1))
    except (IOError, OSError):
        pass
    return None


def _check_pbs_job_log_for_exit(job_id, log_dir):
    """
    Check PBS job log files for exit status when job was killed by PBS.

    When PBS kills a job (e.g., walltime exceeded), our cleanup code doesn't run
    and no final_exit_code file is created. However, PBS appends a resource usage
    block to the job log file that includes the exit status.

    Args:
        job_id: PBS job ID (e.g., "12345.gadi-pbs")
        log_dir: Directory containing job log files

    Returns:
        int: Exit code if found in PBS log, None otherwise
    """
    from pathlib import Path

    log_dir = Path(log_dir)
    if not log_dir.exists():
        return None

    try:
        for log_file in log_dir.glob("*.log"):
            try:
                with open(log_file, "r", encoding="utf-8") as f:
                    content = f.read()

                # Only process log files that belong to this job
                if f"Job Id:             {job_id}" not in content:
                    continue

                import re

                exit_match = re.search(r"Exit Status:\s*(-?\d+)", content)
                if exit_match:
                    return int(exit_match.group(1))

            except (IOError, OSError) as e:
                logger.debug(f"Could not read log file {log_file}: {e}")
                continue

    except Exception as e:
        logger.debug(f"Error searching log files in {log_dir}: {e}")

    return None


def get_job_resource_data(job_id):
    """
    Everything the scheduler knows about a job, or ``None``.

    Asks the platform's scheduler backend (PBS: ``qstat -fx``).

    The returned dictionary contains:
        - exit_status: Job exit code
        - resources_requested: Dict of requested resources
        - resources_used: Dict of actually used resources
        - execution: Dict of execution environment info
        - timing: Dict of job timing information
        - metadata: Dict of job metadata (name, owner, etc.)
        - efficiency: Dict of calculated efficiencies
    """
    logger.debug("Getting resource data for job %s", job_id)
    backend = backends.backend_for_job(job_id)
    return backend.find_details(backend.job_details([job_id]), job_id)


def _format_completion_message(job_id, exit_status, status_dir, key=None):
    """
    Format a completion message with appropriate emoji and hook status information.

    Args:
        job_id: The PBS job ID
        exit_status: The final exit status of the job
        status_dir: Path to the status directory containing hook exit codes
        key: The name the job's status files use (its run ID); defaults
            to :func:`status_file_key`

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

    # Exit codes of --pre, the main command and --post, if they ran
    key = key or status_file_key(job_id, status_dir)
    hook_info = []
    for step in ("pre", "main", "post"):
        exit_file = pathlib.Path(status_dir) / f"{step}_exit_code_{key}"
        if not exit_file.exists():
            continue
        try:
            with open(exit_file, "r", encoding="utf-8") as f:
                step_exit_code = int(f.readline().strip())
        except (OSError, ValueError):
            hook_info.append(f"{step} ❓")
            continue
        if step_exit_code == 0:
            hook_info.append(f"{step} ✅")
        else:
            hook_info.append(f"{step} ❌({step_exit_code})")

    # Add hook information if any hooks were specified
    if hook_info:
        message += f" ({', '.join(hook_info)})"

    return message


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
            logger.debug(f"Error writing to resource collection log: {e}")

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
        logger.warning(f"Failed to fork for background resource collection: {e}")
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


def monitor_job_single_thread(
    job_id,
    out_file,
    err_file,
    quiet=False,
    joblog_file=None,
    verbose=0,
    walltime_str=None,
    joblog_check_interval=60,
    walltime_offset_sec=0,
    run_id=None,
    status_dir=None,
):
    """
    Single-thread job monitoring with proper spinner integration using status files.

    The job writes its status files to ``status_dir`` (default: the
    platform's status directory, passed to it as ``QXUB_STATUS_DIR``),
    named after ``run_id`` (``QXUB_RUN_ID``).

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

    status_dir = Path(status_dir) if status_dir else get_default_status_dir()
    key = run_id or job_id

    # Status files to monitor
    job_started_file = status_dir / f"job_started_{key}"
    final_exit_code_file = status_dir / f"final_exit_code_{key}"

    # Track job state
    job_started_notified = False

    # A job that dies before it writes any status file (its status directory
    # can't be created, the shell can't start, ...) would leave us waiting
    # forever: also look for the exit status PBS appends to the joblog. One
    # small file read, while the job is queued usually a missing file.
    joblog_wait_check = min(joblog_check_interval, 15)
    last_joblog_check = time.monotonic()

    # 4. Spinner while waiting for job to start - check status files every 0.5 seconds
    spinner = SimpleSpinner(quiet=quiet)
    try:
        while True:
            if (
                joblog_file
                and time.monotonic() - last_joblog_check >= joblog_wait_check
            ):
                last_joblog_check = time.monotonic()
                pbs_exit = read_exit_from_joblog_file(joblog_file)
                if pbs_exit is not None and not (
                    job_started_file.exists() or final_exit_code_file.exists()
                ):
                    # The job ended without writing status files here (it
                    # may have written them elsewhere, or none at all).
                    spinner.clear()
                    if not quiet:
                        if pbs_exit == 0:
                            message = "✅ Job completed successfully"
                        else:
                            message = f"❌ Job failed with exit code {pbs_exit}"
                        print_status(
                            f"{message} (from its log, {joblog_file}; "
                            f"no status files in {status_dir})",
                            final=True,
                        )
                    return pbs_exit

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
                logger.info("Job %s completed before starting", job_id)
                exit_status = read_exit_status_from_file(final_exit_code_file)

                # Completion message
                if not quiet:
                    completion_msg = _format_completion_message(
                        job_id, exit_status, status_dir, key
                    )
                    print_status(completion_msg, final=True)

                return exit_status

            time.sleep(0.5)  # Check every 500ms for good responsiveness
    finally:
        spinner.clear()  # Always clear spinner when exiting

    # Derive joblog path early — needed for both walltime-kill detection during
    # streaming and background resource collection afterwards.
    if joblog_file:
        joblog_path = str(joblog_file)
        if verbose >= 2:
            print(f"🔍 DEBUG: Using provided joblog path: {joblog_path}")
    else:
        # Fallback: derive joblog path from out_file path (change .out to .log)
        joblog_path = str(out_file).replace(".out", ".log")
        if verbose >= 2:
            print(f"🔍 DEBUG: Derived joblog path from out_file: {joblog_path}")

    # Parse walltime string to seconds for NFS-friendly joblog gating.
    # walltime_offset_sec shifts the gate: positive = extra grace period after
    # walltime before first check; negative = start checking earlier.
    walltime_sec = None
    if walltime_str:
        try:
            walltime_sec = max(0, time_to_seconds(walltime_str) + walltime_offset_sec)
            logger.debug("Monitoring with walltime gate: %s s", walltime_sec)
        except Exception:
            pass  # Unknown format — fall back to immediate 60 s polling

    # 6. Show output streaming message
    if not quiet:
        print_status("📡 Streaming job output", final=True)

    # 7. Stream job output (also watches joblog for PBS-killed / walltime jobs)
    exit_status = stream_job_output_with_status_files(
        job_id,
        out_file,
        err_file,
        final_exit_code_file,
        quiet=quiet,
        joblog_file=joblog_path,
        walltime_sec=walltime_sec,
        joblog_check_interval=joblog_check_interval,
    )

    # 8. Completion message
    if not quiet:
        completion_msg = _format_completion_message(
            job_id, exit_status, status_dir, key
        )
        print_status(completion_msg, final=True)

    # 9. Start background resource collection and return control immediately
    if not quiet and verbose < 2:
        print_status(
            "🎉 Job completed! Resource collection starting in background...",
            final=True,
        )

    # If user requested very verbose debug (e.g. -vvv) run collection in foreground
    run_foreground = bool(verbose and int(verbose) >= 3)
    if run_foreground and verbose >= 2:
        print(f"🔍 DEBUG: Foreground resource collection enabled (verbose={verbose})")

    start_background_resource_collection(
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
        logger.warning(f"Could not read exit status from {exit_code_file}: {e}")
        return 1  # Default to failure


def stream_job_output_with_status_files(
    job_id,
    out_file,
    err_file,
    final_exit_code_file,
    quiet=False,
    joblog_file=None,
    walltime_sec=None,
    joblog_check_interval=60,
):
    """
    Stream job output files until job completion using status files for detection.

    Args:
        job_id: The PBS job ID
        out_file: Path to stdout file
        err_file: Path to stderr file
        final_exit_code_file: Path to status file containing exit code
        quiet: If True, suppress output streaming (for terse/quiet modes)
        joblog_file: Optional path to the PBS .log file.  When provided, the
            loop also watches for PBS writing its resource-usage block (which
            contains an ``Exit Status:`` line).  This catches jobs killed by
            PBS (e.g. walltime exceeded) before our in-job cleanup code runs
            and writes final_exit_code_file.
        walltime_sec: Expected walltime in seconds (optional, already offset-
            adjusted by caller).  When set, joblog checks are deferred until
            this many seconds have elapsed.  When not set, joblog checks start
            immediately but are still rate-limited to catch unexpected PBS kills.
        joblog_check_interval: Seconds between joblog reads once eligible
            (default 60).  Increase on busy NFS systems to reduce load.

    Returns:
        int: Job exit status
    """
    import os

    # Keep track of file positions to avoid re-reading
    out_pos = 0
    err_pos = 0

    # Joblog polling state -- NFS-friendly rate limiting.
    # - If walltime_sec is known: don't even open the joblog until that many
    #   seconds have elapsed (zero NFS reads during normal operation).
    # - If walltime_sec is unknown: start checking immediately as a safety net
    #   for unexpected PBS kills (node failure, qdel, etc.).
    # Either way, cap checks at one per joblog_check_interval seconds.
    monitoring_start = time.time()
    last_joblog_check = None  # timestamp of last actual joblog open

    # Monitor job completion via status file
    while True:
        # Stream new output from STDOUT (skip if quiet mode)
        if not quiet and os.path.exists(out_file):
            try:
                with open(out_file, "r", encoding="utf-8") as f:
                    f.seek(out_pos)
                    for line in f:
                        print(line.rstrip(), file=sys.stdout, flush=True)
                    out_pos = f.tell()
            except (OSError, IOError):
                pass
        elif quiet and os.path.exists(out_file):
            # Still need to track position even in quiet mode for final read
            try:
                with open(out_file, "r", encoding="utf-8") as f:
                    f.seek(0, 2)  # Seek to end
                    out_pos = f.tell()
            except (OSError, IOError):
                pass

        # Stream new output from STDERR (skip if quiet mode)
        if not quiet and os.path.exists(err_file):
            try:
                with open(err_file, "r", encoding="utf-8") as f:
                    f.seek(err_pos)
                    for line in f:
                        print(line.rstrip(), file=sys.stderr, flush=True)
                    err_pos = f.tell()
            except (OSError, IOError):
                pass
        elif quiet and os.path.exists(err_file):
            # Still need to track position even in quiet mode for final read
            try:
                with open(err_file, "r", encoding="utf-8") as f:
                    f.seek(0, 2)  # Seek to end
                    err_pos = f.tell()
            except (OSError, IOError):
                pass

        # Check if job finished using status file (normal completion path)
        if final_exit_code_file.exists():
            # Read any final output (skip if quiet mode)
            if not quiet and os.path.exists(out_file):
                try:
                    with open(out_file, "r", encoding="utf-8") as f:
                        f.seek(out_pos)
                        for line in f:
                            print(line.rstrip(), file=sys.stdout, flush=True)
                except (OSError, IOError):
                    pass

            if not quiet and os.path.exists(err_file):
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

        # Check PBS job log for exit status (walltime-kill / PBS-kill path).
        # PBS writes an Exit Status line to the .log file even when it kills
        # a job before our cleanup code can write final_exit_code_file.
        # We rate-limit to joblog_check_interval and defer until walltime
        # has elapsed (when known) to minimise NFS pressure.
        if joblog_file:
            now = time.time()
            elapsed = now - monitoring_start
            past_walltime = walltime_sec is None or elapsed >= walltime_sec
            due_for_check = last_joblog_check is None or (
                now - last_joblog_check >= joblog_check_interval
            )
            if past_walltime and due_for_check:
                last_joblog_check = now
                pbs_exit = read_exit_from_joblog_file(joblog_file)
                if pbs_exit is not None:
                    # Flush whatever partial output the job managed to write
                    if not quiet and os.path.exists(out_file):
                        try:
                            with open(out_file, "r", encoding="utf-8") as f:
                                f.seek(out_pos)
                                for line in f:
                                    print(line.rstrip(), file=sys.stdout, flush=True)
                        except (OSError, IOError):
                            pass
                    if not quiet and os.path.exists(err_file):
                        try:
                            with open(err_file, "r", encoding="utf-8") as f:
                                f.seek(err_pos)
                                for line in f:
                                    print(line.rstrip(), file=sys.stderr, flush=True)
                        except (OSError, IOError):
                            pass
                    logger.debug(
                        "Job %s terminated by PBS (exit %s), detected via job log",
                        job_id,
                        pbs_exit,
                    )
                    return pbs_exit

        time.sleep(0.2)  # Check every 200ms for good output responsiveness
