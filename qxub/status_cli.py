"""
Status CLI for qxub - view job status and monitoring.

Provides the `qxub status` command for viewing job status without qstat polling.
"""

import json
import logging
import os

logger = logging.getLogger(__name__)
import sys
from datetime import datetime, timedelta
from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from .core.scheduler import job_status_from_files
from .queue import is_virtual_id, resolve_virtual_id
from .resources import resource_tracker
from .resources.tracker import DatabaseError

console = Console()


def _get_status_check_logger():
    """Return a logger that writes to ~/.local/share/qxub/status_check.log."""
    _logger = logging.getLogger("qxub.status_check_audit")
    if not _logger.handlers:
        log_dir = (
            Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
            / "qxub"
        )
        log_dir.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(log_dir / "status_check.log")
        handler.setFormatter(
            logging.Formatter("%(asctime)s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        )
        _logger.addHandler(handler)
        _logger.setLevel(logging.DEBUG)
    return _logger


def _log_status_check(job_id, result, source, detail=""):
    """Append a one-line audit record for a qxtat check invocation."""
    try:
        audit = _get_status_check_logger()
        msg = f"job={job_id}  result={result}  source={source}"
        if detail:
            msg += f"  {detail}"
        audit.info(msg)
    except Exception:
        pass  # Audit logging must never break status checks


@click.group(invoke_without_command=True, name="status")
@click.pass_context
def status_cli(ctx):
    """View job status and manage job tracking database."""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@status_cli.command()
@click.option(
    "--status", help="Filter by status (submitted, running, completed, failed)"
)
@click.option("--limit", default=20, help="Maximum number of jobs to show")
@click.option(
    "--all", "show_all", is_flag=True, help="Show all jobs (ignore 7-day limit)"
)
@click.option("--id-only", is_flag=True, help="Output only job IDs (useful for piping)")
def list(status, limit, show_all, id_only):
    """List recent jobs with their status."""

    if status:
        jobs = resource_tracker.get_jobs_by_status(status=status, limit=limit)
        title = f"Jobs with status: {status}"
    else:
        jobs = resource_tracker.get_jobs_by_status(limit=limit)
        title = "Recent Jobs"

    if not jobs:
        if not id_only:
            console.print(f"[yellow]No jobs found[/yellow]")
        return

    # If --id-only, just print job IDs for piping
    if id_only:
        for job in jobs:
            click.echo(job["job_id"])
        return

    # Create table
    table = Table(title=title, show_header=True, header_style="bold magenta")
    table.add_column("Job ID", style="cyan", no_wrap=True)
    table.add_column("Status", justify="center")
    table.add_column("Command", style="green", max_width=40)
    table.add_column("Submitted", style="blue")
    table.add_column("Duration", justify="right")

    for job in jobs:
        # Status with emoji
        status_display = _format_status(job["status"])

        # Command truncation
        command = job["command"] or "N/A"
        if len(command) > 40:
            command = command[:37] + "..."

        # Time formatting
        submitted = _format_time(job["submitted_at"])
        duration = _calculate_duration(job)

        # Clean job ID display
        job_id = job["job_id"]
        if job_id.endswith(".gadi-pbs"):
            job_id = job_id[:-9]  # Remove .gadi-pbs suffix

        table.add_row(job_id, status_display, command, submitted, duration)

    console.print(table)


@status_cli.command()
def summary():
    """Show summary of job counts by status."""

    status_counts = resource_tracker.get_status_summary()

    if not status_counts:
        console.print("[yellow]No recent jobs found[/yellow]")
        return

    # Create summary table
    table = Table(
        title="Job Status Summary (Last 7 Days)",
        show_header=True,
        header_style="bold magenta",
    )
    table.add_column("Status", style="cyan")
    table.add_column("Count", justify="right", style="bold")

    # Order statuses for better display
    status_order = ["submitted", "running", "completed", "failed"]

    for status in status_order:
        if status in status_counts:
            status_display = _format_status(status)
            table.add_row(status_display, str(status_counts[status]))

    # Add any other statuses not in the standard order
    for status, count in status_counts.items():
        if status not in status_order:
            status_display = _format_status(status)
            table.add_row(status_display, str(count))

    console.print(table)


@status_cli.command()
@click.argument("job_id")
def show(job_id):
    """Show detailed status for a specific job."""

    # Add .gadi-pbs suffix if not present
    if not job_id.endswith(".gadi-pbs"):
        job_id = f"{job_id}.gadi-pbs"

    job = resource_tracker.get_job_status(job_id)

    if not job:
        console.print(f"[red]Job {job_id} not found[/red]")
        return

    # Display job details
    console.print(f"\n[bold magenta]Job Details: {job_id}[/bold magenta]")
    console.print(f"Status: {_format_status(job['status'])}")
    console.print(f"Command: [green]{job['command'] or 'N/A'}[/green]")

    if job["submitted_at"]:
        console.print(f"Submitted: [blue]{_format_time(job['submitted_at'])}[/blue]")

    if job["started_at"]:
        console.print(f"Started: [blue]{_format_time(job['started_at'])}[/blue]")

    if job["completed_at"]:
        console.print(f"Completed: [blue]{_format_time(job['completed_at'])}[/blue]")

    if job["exit_code"] is not None:
        exit_style = "green" if job["exit_code"] == 0 else "red"
        console.print(f"Exit Code: [{exit_style}]{job['exit_code']}[/{exit_style}]")

    duration = _calculate_duration(job)
    if duration != "N/A":
        console.print(f"Duration: {duration}")


@status_cli.command()
@click.option("--days", default=30, help="Remove jobs older than this many days")
@click.option(
    "--dry-run", is_flag=True, help="Show what would be deleted without deleting"
)
def cleanup(days, dry_run):
    """Clean up old job records from the database."""

    if dry_run:
        console.print(
            f"[yellow]DRY RUN: Would delete jobs older than {days} days[/yellow]"
        )
        # TODO: Add actual dry-run logic to show what would be deleted
        return

    deleted_count = resource_tracker.cleanup_old_jobs(days_old=days)

    if deleted_count > 0:
        console.print(f"[green]✅ Cleaned up {deleted_count} old job records[/green]")
    else:
        console.print("[blue]No old job records to clean up[/blue]")


@status_cli.command()
@click.argument("job_id")
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["snakemake", "json", "exitcode"], case_sensitive=False),
    default="snakemake",
    help="Output format for workflow engine integration",
)
@click.option(
    "--snakemake",
    is_flag=True,
    help="Force snakemake output format (alias for --format=snakemake)",
)
@click.option(
    "--log-dir",
    "log_dir",
    default=None,
    help="Directory containing PBS job log files (avoids path guessing)",
)
def check(job_id, output_format, snakemake, log_dir):
    """Check job status for workflow engine integration.

    Returns machine-readable status information suitable for workflow engines
    like Snakemake, Nextflow, or CWL.

    Exit codes:
    - 0: Normal operation (status successfully determined)
    - 1: Error (invalid job ID, system error, etc.)

    For --format=snakemake (default):
    - Outputs: "success", "running", or "failed"

    For --format=json:
    - Outputs: {"status": "...", "job_id": "...", "timestamp": "..."}

    For --format=exitcode:
    - Outputs: (empty), exit code 0=running, 1=failed, 2=completed
    """
    if snakemake:
        output_format = "snakemake"

    # -----------------------------------------------------------------------
    # Dispatch-on-status-check hook (Phase 3 will fill this in; no-op now)
    # -----------------------------------------------------------------------
    _maybe_dispatch_pending()

    # -----------------------------------------------------------------------
    # Virtual job ID resolution (Phase 2+)
    # -----------------------------------------------------------------------
    if is_virtual_id(job_id):
        entry = resolve_virtual_id(job_id)
        if entry is None:
            if output_format == "json":
                print(
                    json.dumps({"error": "Virtual job ID not found", "job_id": job_id})
                )
            else:
                sys.stderr.write(f"qxub: Virtual job ID {job_id} not found\n")
            sys.exit(1)

        virtual_status = entry.get("status", "unknown")

        if virtual_status in ("initiated", "pending"):
            # Still waiting to be dispatched to PBS
            _output_status("running", None, job_id, output_format)
            return

        if virtual_status == "cancelled":
            _output_status("failed", 1, job_id, output_format)
            return

        if virtual_status in ("completed", "failed"):
            exit_code = entry.get("exit_code")
            _output_status(virtual_status, exit_code, job_id, output_format)
            return

        # status == 'dispatched' (or any other): resolve to real PBS job ID
        pbs_job_id = entry.get("pbs_job_id")
        if not pbs_job_id:
            # No PBS ID yet — treat as still running
            _output_status("running", None, job_id, output_format)
            return

        # Continue below with the real PBS job ID
        job_id = pbs_job_id

    # -----------------------------------------------------------------------
    # Standard PBS job ID lookup — files-first to minimise DB pressure
    # -----------------------------------------------------------------------
    # Add .gadi-pbs suffix if not present
    if not job_id.endswith(".gadi-pbs"):
        job_id = f"{job_id}.gadi-pbs"

    # -------------------------------------------------------------------
    # 1. Check status files on disk FIRST (no DB, no locking)
    # -------------------------------------------------------------------
    file_status, file_exit_code = job_status_from_files(job_id, log_dir=log_dir)

    if file_status in ("C", "F"):
        # Job has finished — persist to DB (best-effort) then report.
        status = "completed" if file_status == "C" else "failed"
        _persist_completion(job_id, status, file_exit_code)
        _log_status_check(
            job_id,
            status,
            "files",
            f"file_status={file_status} exit_code={file_exit_code} log_dir={log_dir}",
        )
        _output_status(status, file_exit_code, job_id, output_format)
        return

    if file_status == "R":
        # Job is still running according to status files — no DB needed.
        _log_status_check(job_id, "running", "files", f"log_dir={log_dir}")
        _output_status("running", None, job_id, output_format)
        return

    # -------------------------------------------------------------------
    # 2. No conclusive file evidence (status == "Q" / unknown).
    #    Fall back to the DB for jobs that haven't written files yet
    #    (e.g. just submitted, or pre-file-status era).
    # -------------------------------------------------------------------
    try:
        job_info = resource_tracker.get_job_status(job_id)
    except DatabaseError as exc:
        logger.debug("Database error during status check: %s", exc)
        job_info = None

    if job_info:
        db_status = job_info.get("status", "unknown")
        if db_status in ("completed", "failed"):
            _log_status_check(
                job_id,
                db_status,
                "db",
                f"exit_code={job_info.get('exit_code')}",
            )
            _output_status(db_status, job_info.get("exit_code"), job_id, output_format)
            return
        # DB says submitted/running but no files yet — still in flight.
        if db_status in ("submitted", "running"):
            _log_status_check(job_id, "running", "db", f"db_status={db_status}")
            _output_status("running", None, job_id, output_format)
            return

    # Genuinely unknown — report "running" so Snakemake retries.
    _log_status_check(
        job_id,
        "running",
        "unknown",
        f"file_status={file_status} job_info={'found' if job_info else 'none'}"
        f" log_dir={log_dir}",
    )
    _output_status("running", None, job_id, output_format)


def _persist_completion(job_id: str, status: str, exit_code) -> None:
    """Best-effort DB update after file-based completion detection.

    Persists the resolved status and exit code so that subsequent calls
    (and other tools like ``qxub status list``) see the final state.
    Also backfills resource efficiency data from the PBS joblog.

    All DB operations are wrapped in try/except — a locked or unavailable
    database must never prevent the caller from reporting the correct status.
    """
    try:
        resource_tracker.update_job_status(job_id, status)
        if exit_code is not None:
            resource_tracker.update_job_exit_code(job_id, exit_code)
    except Exception as exc:
        _log_status_check(job_id, "db-write-failed", "persist", f"error={exc}")

    # Backfill resource efficiency from the PBS joblog stored on disk.
    try:
        import os as _os

        from .history import history_manager
        from .resources import parse_joblog_resources

        files = history_manager.get_execution_file_paths(job_id)
        joblog = (files or {}).get("joblog")
        if joblog and _os.path.exists(joblog):
            resource_data = parse_joblog_resources(joblog)
            if resource_data:
                resource_tracker.update_job_resources(job_id, resource_data)
    except Exception:
        pass  # Non-fatal: resource data remains empty until next check


def _maybe_dispatch_pending() -> None:
    """Trigger dispatch of pending queue entries if headroom allows.

    This is a no-op in Phase 2 (all entries are dispatched immediately).
    Phase 3 will implement the actual headroom check and dispatch logic here.
    """
    # TODO Phase 3: check active_count vs job_limit, dispatch pending entries
    pass


def _output_status(
    status: str,
    exit_code,
    job_id: str,
    output_format: str,
) -> None:
    """Emit status in the requested output format and exit appropriately."""
    if output_format == "snakemake":
        if status == "completed":
            print("success" if (exit_code is None or exit_code == 0) else "failed")
        elif status in ("initiated", "dispatched", "submitted", "running"):
            print("running")
        else:
            print("failed")

    elif output_format == "json":
        result = {
            "status": status,
            "job_id": job_id,
            "timestamp": datetime.now().isoformat(),
            "exit_code": exit_code,
        }
        print(json.dumps(result))

    elif output_format == "exitcode":
        if status == "completed":
            sys.exit(2 if (exit_code is None or exit_code == 0) else 1)
        elif status in ("initiated", "dispatched", "submitted", "running"):
            sys.exit(0)
        else:
            sys.exit(1)


def _format_status(status):
    """Format status with emoji and color."""
    status_map = {
        "initiated": "⚪ Initiated",
        "dispatched": "🟠 Dispatched",
        "submitted": "🟡 Submitted",
        "running": "🔵 Running",
        "completed": "✅ Completed",
        "failed": "❌ Failed",
        "unknown": "❓ Unknown",
    }
    return status_map.get(status, f"❓ {status}")


def _format_time(timestamp_str):
    """Format timestamp for display."""
    if not timestamp_str:
        return "N/A"

    try:
        dt = datetime.fromisoformat(timestamp_str.replace("Z", "+00:00"))
        now = datetime.now(dt.tzinfo) if dt.tzinfo else datetime.now()

        # If within last 24 hours, show relative time
        diff = now - dt
        if diff < timedelta(days=1):
            hours = int(diff.total_seconds() / 3600)
            minutes = int((diff.total_seconds() % 3600) / 60)
            if hours > 0:
                return f"{hours}h {minutes}m ago"
            elif minutes > 0:
                return f"{minutes}m ago"
            else:
                return "just now"
        else:
            return dt.strftime("%m-%d %H:%M")
    except Exception:
        return timestamp_str[:16]  # Fallback to raw timestamp


def _calculate_duration(job):
    """Calculate job duration."""
    if not job["started_at"]:
        return "N/A"

    try:
        start = datetime.fromisoformat(job["started_at"].replace("Z", "+00:00"))

        if job["completed_at"]:
            end = datetime.fromisoformat(job["completed_at"].replace("Z", "+00:00"))
        else:
            end = datetime.now(start.tzinfo) if start.tzinfo else datetime.now()

        duration = end - start

        if duration.days > 0:
            return f"{duration.days}d {duration.seconds // 3600}h"
        elif duration.seconds >= 3600:
            hours = duration.seconds // 3600
            minutes = (duration.seconds % 3600) // 60
            return f"{hours}h {minutes}m"
        elif duration.seconds >= 60:
            minutes = duration.seconds // 60
            seconds = duration.seconds % 60
            return f"{minutes}m {seconds}s"
        else:
            return f"{duration.seconds}s"

    except Exception:
        return "N/A"
