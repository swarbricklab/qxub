"""
Status CLI for qxub - view job status and monitoring.

Provides the `qxub status` command for viewing job status without qstat polling.
"""

from datetime import datetime, timedelta

import click
from rich.console import Console
from rich.table import Table
from rich.text import Text

from .resource_tracker import resource_tracker

console = Console()


@click.group(name="status")
def status_cli():
    """View job status and manage job tracking database."""
    pass


@status_cli.command()
@click.option(
    "--status", help="Filter by status (submitted, running, completed, failed)"
)
@click.option("--limit", default=20, help="Maximum number of jobs to show")
@click.option(
    "--all", "show_all", is_flag=True, help="Show all jobs (ignore 7-day limit)"
)
def list(status, limit, show_all):
    """List recent jobs with their status."""

    if status:
        jobs = resource_tracker.get_jobs_by_status(status=status, limit=limit)
        title = f"Jobs with status: {status}"
    else:
        jobs = resource_tracker.get_jobs_by_status(limit=limit)
        title = "Recent Jobs"

    if not jobs:
        console.print(f"[yellow]No jobs found[/yellow]")
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


def _format_status(status):
    """Format status with emoji and color."""
    status_map = {
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
            return f"{duration.days}d {duration.seconds//3600}h"
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
