"""
Monitor CLI for qxub - monitor multiple PBS jobs until completion.

This module provides the `qxub monitor` command for monitoring multiple
PBS jobs simultaneously. It accepts job IDs from stdin, arguments, or both,
and blocks until all jobs are complete.
"""

import logging
import signal
import sys
import threading
import time
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import click
from rich.console import Console
from rich.live import Live
from rich.table import Table
from rich.text import Text

from .config import config_manager
from .core.scheduler import get_default_status_dir, job_status_from_files

console = Console()


class MultiJobMonitor:
    """Monitor multiple PBS jobs simultaneously using file-based status checking."""

    def __init__(
        self,
        job_ids: List[str],
        quiet: bool = False,
        refresh_interval: int = 30,
        suffix_to_remove: Optional[str] = None,
        show_names_only: bool = False,
        show_job_ids_only: bool = False,
        summary_mode: bool = False,
        status_dir: Optional[str] = None,
    ):
        """
        Initialize the multi-job monitor.

        Args:
            job_ids: List of PBS job IDs to monitor
            quiet: If True, show minimal output
            refresh_interval: Seconds between status checks
            suffix_to_remove: Optional suffix to remove from job IDs in display
            show_names_only: If True, only show job names
            show_job_ids_only: If True, only show job IDs
            summary_mode: If True, show only summary counts
            status_dir: Path to status directory. If None, uses default.
        """
        self.job_ids = set(job_ids)
        self.quiet = quiet
        self.refresh_interval = refresh_interval
        self.suffix_to_remove = suffix_to_remove
        self.show_names_only = show_names_only
        self.show_job_ids_only = show_job_ids_only
        self.summary_mode = summary_mode
        self.shutdown_requested = threading.Event()

        # Status directory for file-based monitoring
        self.status_dir = status_dir if status_dir else get_default_status_dir()

        # Track job states and metadata
        self.job_states: Dict[str, str] = {}
        self.job_names: Dict[str, str] = {}
        self.job_exit_codes: Dict[str, Optional[int]] = {}
        self.completed_jobs: Set[str] = set()
        self.failed_jobs: Set[str] = set()  # Jobs with non-zero exit codes

        # Initialize all jobs as queued (waiting for status files)
        for job_id in self.job_ids:
            self.job_states[job_id] = "Q"
            self.job_names[job_id] = "Unknown"
            self.job_exit_codes[job_id] = None

        # Set up signal handler for graceful shutdown
        signal.signal(signal.SIGINT, self._signal_handler)

    def _signal_handler(self, signum, frame):
        """Handle SIGINT (Ctrl+C) gracefully."""
        logging.info("Interrupt received, shutting down monitor...")
        self.shutdown_requested.set()

    def _strip_suffix(self, job_id: str) -> str:
        """Remove suffix from job ID for display."""
        if self.suffix_to_remove and job_id.endswith(self.suffix_to_remove):
            return job_id[: -len(self.suffix_to_remove)]
        return job_id

    def _check_job_status(
        self, job_id: str
    ) -> Tuple[str, Optional[str], Optional[int]]:
        """
        Check status of a single job using file-based monitoring.

        This approach is faster and more reliable than polling qstat,
        especially when the PBS scheduler is under heavy load.

        Args:
            job_id: PBS job ID

        Returns:
            Tuple of (status, job_name, exit_code)
            Status can be Q (queued), R (running), C (completed), F (failed)
        """
        try:
            # Use file-based status checking (fast, no qstat dependency)
            status, exit_code = job_status_from_files(job_id, self.status_dir)

            # Extract job name from job_id if possible (e.g., "qx-20260128_184415")
            # The job name is typically the prefix before the timestamp
            job_name = None
            # We don't have the job name from files, leave as None
            # The display will use "Unknown" which is fine

            return status, job_name, exit_code

        except Exception as e:
            logging.debug(f"Failed to check status for job {job_id}: {e}")
            return "Q", None, None  # Assume queued if we can't check

    def _get_status_emoji(self, status: str, exit_code: Optional[int] = None) -> str:
        """Get emoji representation of job status."""
        if status == "Q":
            return "⏳"  # Queued/Waiting
        elif status == "R":
            return "🔄"  # Running
        elif status == "C":  # Completed successfully
            return "✅"
        elif status == "F":  # Failed
            return "❌"
        else:
            return "❓"  # Unknown

    def _update_all_statuses(self) -> None:
        """Update status for all monitored jobs."""
        for job_id in self.job_ids:
            if job_id in self.completed_jobs:
                continue

            status, job_name, exit_code = self._check_job_status(job_id)

            # Update job metadata
            if job_name:
                self.job_names[job_id] = job_name
            if exit_code is not None:
                self.job_exit_codes[job_id] = exit_code

            if status in ["C", "F"]:  # Completed or Failed
                self.completed_jobs.add(job_id)
                self.job_states[job_id] = status
                # Track failed jobs (non-zero exit codes)
                if status == "F" or (exit_code is not None and exit_code != 0):
                    self.failed_jobs.add(job_id)
            else:
                self.job_states[job_id] = status

    def _create_status_table(self) -> Table:
        """Create a rich table showing job statuses."""
        table = Table(title="Job Monitoring Status")

        # Add columns based on display mode
        if self.show_names_only:
            table.add_column("Name", style="cyan", no_wrap=True)
        elif self.show_job_ids_only:
            table.add_column("Job ID", style="cyan", no_wrap=True)
        else:
            table.add_column("Job ID", style="cyan", no_wrap=True)
            table.add_column("Name", style="blue", no_wrap=True)

        table.add_column("Status", style="magenta", justify="center")

        for job_id in sorted(self.job_ids):
            status = self.job_states.get(job_id, "?")
            job_name = self.job_names.get(job_id, "Unknown")
            exit_code = self.job_exit_codes.get(job_id)

            # Get emoji status
            emoji = self._get_status_emoji(status, exit_code)

            # Format display ID (remove suffix if requested)
            display_id = self._strip_suffix(job_id)

            # Add row based on display mode
            if self.show_names_only:
                table.add_row(job_name, emoji)
            elif self.show_job_ids_only:
                table.add_row(display_id, emoji)
            else:
                table.add_row(display_id, job_name, emoji)

        return table

    def _get_summary_stats(self) -> Dict[str, int]:
        """Get summary statistics of job states."""
        stats = defaultdict(int)
        for status in self.job_states.values():
            stats[status] += 1
        return dict(stats)

    def _create_summary_display(self) -> Table:
        """Create a compact summary table showing job state counts."""
        table = Table(title=f"Job Summary ({len(self.job_ids)} total jobs)")
        table.add_column("State", style="cyan", no_wrap=True)
        table.add_column("Count", style="green", justify="right")
        table.add_column("Status", style="magenta", justify="center")

        # Get statistics
        stats = self._get_summary_stats()

        # Define state display order and mapping for file-based monitoring
        state_display = [
            ("Q", "Queued", "⏳"),
            ("R", "Running", "🔄"),
            ("C", "Completed", "✅"),
            ("F", "Failed", "❌"),
        ]

        # Count completed vs failed jobs for more detailed stats
        completed_success = 0
        completed_failed = 0
        for job_id in self.job_ids:
            status = self.job_states.get(job_id, "Q")
            if status == "C":
                completed_success += 1
            elif status == "F":
                completed_failed += 1

        # Add rows for states that have jobs
        for state_code, state_name, emoji in state_display:
            count = stats.get(state_code, 0)
            if count > 0:
                table.add_row(state_name, str(count), emoji)

        return table

    def monitor_jobs(self) -> int:
        """
        Monitor all jobs until completion.

        Returns:
            Exit code: 0 if all jobs completed with exit code 0, 1 if any failed
        """
        if not self.job_ids:
            click.echo("No job IDs to monitor.", err=True)
            return 1

        logging.info(f"Monitoring {len(self.job_ids)} jobs...")

        # Initial status check
        self._update_all_statuses()

        if self.quiet:
            # Quiet mode: just show progress counts
            while True:
                if self.shutdown_requested.is_set():
                    click.echo("\nMonitoring interrupted")
                    return 130  # SIGINT exit code

                # Check if all jobs are done
                if len(self.completed_jobs) >= len(self.job_ids):
                    break

                click.echo(
                    f"{len(self.completed_jobs)}/{len(self.job_ids)} jobs complete"
                )

                # Wait with interrupt checking and countdown
                self._wait_with_countdown("Checking again in")
                if self.shutdown_requested.is_set():
                    return 130

                self._update_all_statuses()

            # Final status update to get exit codes for completed jobs
            self._update_all_statuses()

            # Final summary
            successful_jobs = len(self.completed_jobs) - len(self.failed_jobs)
            click.echo(
                f"All jobs complete: {successful_jobs} successful, {len(self.failed_jobs)} failed"
            )

        else:
            # Rich interactive mode
            try:
                # Choose display method based on summary mode
                display_method = (
                    self._create_summary_display
                    if self.summary_mode
                    else self._create_status_table
                )

                with Live(
                    display_method(), refresh_per_second=0.5, console=console
                ) as live:
                    while True:
                        if self.shutdown_requested.is_set():
                            console.print("\n[yellow]Monitoring interrupted[/yellow]")
                            return 130

                        # Check if all jobs are done
                        if len(self.completed_jobs) >= len(self.job_ids):
                            # Do a final update to get exit codes
                            self._update_all_statuses()
                            # Show final status with updated exit codes
                            live.update(display_method())
                            break

                        # Update display
                        live.update(display_method())

                        # Wait with interrupt checking and countdown
                        self._wait_with_countdown_rich(live, display_method)
                        if self.shutdown_requested.is_set():
                            return 130

                        self._update_all_statuses()

                # Final summary outside the Live context
                console.print("\n[bold green]Monitoring Complete![/bold green]")
                successful_jobs = len(self.completed_jobs) - len(self.failed_jobs)
                console.print(f"✅ Successful: {successful_jobs}")
                if self.failed_jobs:
                    console.print(f"❌ Failed: {len(self.failed_jobs)}")

            except Exception as e:
                console.print(f"[red]Error during monitoring: {e}[/red]")
                return 1

        # Return appropriate exit code - 0 only if ALL jobs succeeded
        if self.failed_jobs:
            return 1  # Some jobs failed
        return 0  # All jobs completed successfully with exit code 0

    def _wait_with_countdown(self, message: str) -> None:
        """Wait with countdown display for quiet mode."""
        for remaining in range(self.refresh_interval, 0, -1):
            if self.shutdown_requested.is_set():
                return
            click.echo(f"\r{message} {remaining} seconds...", nl=False)
            time.sleep(1)
        click.echo()  # New line after countdown

    def _wait_with_countdown_rich(self, live: Live, display_method) -> None:
        """Wait with countdown display for rich mode."""
        for remaining in range(self.refresh_interval, 0, -1):
            if self.shutdown_requested.is_set():
                return

            # Create table with countdown using the appropriate display method
            table = display_method()
            countdown_text = Text(f"Next update in {remaining} seconds...", style="dim")
            table.caption = countdown_text
            live.update(table)
            time.sleep(1)


@click.command(name="monitor")
@click.argument("job_ids", nargs=-1)
@click.option(
    "--quiet", "-q", is_flag=True, help="Quiet mode: show minimal progress output"
)
@click.option(
    "--interval", "-i", default=30, help="Refresh interval in seconds (default: 30)"
)
@click.option(
    "--suffix", help="Suffix to remove from job IDs in display (e.g., '.gadi-pbs')"
)
@click.option("--name-only", is_flag=True, help="Show only job names in display")
@click.option("--job-id-only", is_flag=True, help="Show only job IDs in display")
@click.option(
    "--summary",
    is_flag=True,
    help="Show only summary counts of jobs in each state (useful for many jobs)",
)
@click.pass_context
def monitor_cli(ctx, job_ids, quiet, interval, suffix, name_only, job_id_only, summary):
    """
    Monitor PBS jobs until completion.

    Monitor multiple PBS jobs and block until all are complete.
    Job IDs can be provided as arguments or read from stdin.

    Uses file-based status checking (fast, no qstat dependency) by monitoring
    status files created by qxub job scripts in the status directory.

    The exit code is 0 only if ALL monitored jobs complete with exit code 0.
    If any job fails (non-zero exit code), the exit code is 1.

    Examples:

        # Monitor specific jobs
        qxub monitor 12345.gadi-pbs 12346.gadi-pbs

        # Monitor jobs from pipeline, hiding .gadi-pbs suffix
        find -name "*.csv" -exec qxub --terse --env myenv -- process.py {} \\; | qxub monitor --suffix .gadi-pbs

        # Show only job names
        qxub monitor --name-only 12345.gadi-pbs 12346.gadi-pbs

        # Quiet mode for scripting
        echo "12345.gadi-pbs" | qxub monitor --quiet
    """
    # Validate mutually exclusive options
    if name_only and job_id_only:
        click.echo("Error: Cannot specify both --name-only and --job-id-only", err=True)
        return 1

    # Get defaults from config if not provided
    try:
        config = config_manager.load_config()
        monitor_config = config.get("monitor", {})

        if suffix is None:
            suffix = monitor_config.get("default_suffix", None)

        # Get default interval if not specified or use config default
        if interval == 30:  # Check if it's the default value
            config_interval = monitor_config.get("default_interval", None)
            if config_interval is not None:
                interval = config_interval

    except Exception:
        # If config loading fails, use command line defaults
        pass

    # Collect job IDs from arguments
    all_job_ids = list(job_ids)

    # Read additional job IDs from stdin if available
    if not sys.stdin.isatty():
        try:
            stdin_input = sys.stdin.read().strip()
            if stdin_input:
                stdin_job_ids = [
                    line.strip() for line in stdin_input.splitlines() if line.strip()
                ]
                all_job_ids.extend(stdin_job_ids)
        except Exception as e:
            click.echo(f"Error reading from stdin: {e}", err=True)
            return 1

    # Validate we have job IDs
    if not all_job_ids:
        click.echo(
            "Error: No job IDs provided. Provide job IDs as arguments or via stdin.",
            err=True,
        )
        click.echo("Example: qxub monitor 12345.gadi-pbs 12346.gadi-pbs", err=True)
        return 1

    # Remove duplicates while preserving order
    unique_job_ids = []
    seen = set()
    for job_id in all_job_ids:
        if job_id not in seen:
            unique_job_ids.append(job_id)
            seen.add(job_id)

    # Basic validation of job ID format
    valid_job_ids = []
    for job_id in unique_job_ids:
        if job_id and ("." in job_id or job_id.isdigit()):
            valid_job_ids.append(job_id)
        else:
            click.echo(f"Warning: Skipping invalid job ID format: {job_id}", err=True)

    if not valid_job_ids:
        click.echo("Error: No valid job IDs found.", err=True)
        return 1

    # Create monitor and start monitoring
    monitor = MultiJobMonitor(
        valid_job_ids,
        quiet=quiet,
        refresh_interval=interval,
        suffix_to_remove=suffix,
        show_names_only=name_only,
        show_job_ids_only=job_id_only,
        summary_mode=summary,
    )
    exit_code = monitor.monitor_jobs()
    ctx.exit(exit_code)
