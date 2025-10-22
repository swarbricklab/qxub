"""
Cancel CLI for qxub - Cancel PBS jobs with qxub database integration.

Provides job cancellation functionality for snakemake compatibility and
general job management. Supports cancellation by job ID, name, or pattern
matching with automatic database status updates.
"""

import logging
import re
import subprocess
import time
from typing import List, Optional

import click

from .resource_tracker import ResourceTracker
from .scheduler import get_job_resource_data, qdel


@click.group(name="cancel")
def cancel_cli():
    """Cancel PBS jobs with database tracking."""
    pass


@cancel_cli.command(name="id")
@click.argument("job_ids", nargs=-1, required=False)
@click.option(
    "--dry-run",
    "-n",
    is_flag=True,
    help="Show what would be cancelled without doing it",
)
@click.option("--force", "-f", is_flag=True, help="Skip confirmation prompt")
@click.option("--quiet", "-q", is_flag=True, help="Minimize output messages")
@click.option(
    "--wait", "-w", is_flag=True, help="Wait for job deletion to be confirmed"
)
@click.option(
    "--timeout", default=30, type=int, help="Maximum wait time in seconds (default: 30)"
)
def cancel_by_id(job_ids, dry_run, force, quiet, wait, timeout):
    """Cancel jobs by PBS job ID(s)."""
    import sys

    # Collect job IDs from arguments or stdin
    all_job_ids = list(job_ids) if job_ids else []

    # Auto-detect stdin if no job IDs provided and stdin is available
    if not all_job_ids and not sys.stdin.isatty():
        # Read job IDs from stdin
        if not quiet:
            click.echo("Reading job IDs from stdin...")
        try:
            for line in sys.stdin:
                job_id = line.strip()
                if job_id:  # Skip empty lines
                    all_job_ids.append(job_id)
        except KeyboardInterrupt:
            click.echo("\n❌ Cancelled reading from stdin")
            return 1

    if not all_job_ids:
        click.echo("❌ No job IDs provided", err=True)
        click.echo("Examples:")
        click.echo("  qxub cancel id 123456 123457")
        click.echo("  echo '123456' | qxub cancel id")
        click.echo("  qxub status list --status running --id-only | qxub cancel id")
        return 1

    # Dry run - just show what would be cancelled
    if dry_run:
        click.echo(f"🔍 Dry run: Would cancel {len(all_job_ids)} job(s):")
        for job_id in all_job_ids:
            click.echo(f"  • {job_id}")
        return 0

    # Confirmation prompt
    if not force:
        job_list = ", ".join(all_job_ids)
        if not click.confirm(f"Cancel {len(all_job_ids)} job(s): {job_list}?"):
            click.echo("❌ Cancellation aborted")
            return 0

    # Cancel each job
    success_count = 0
    failed_jobs = []
    tracker = ResourceTracker()

    for job_id in all_job_ids:
        if not quiet:
            click.echo(f"🗑️  Cancelling job {job_id}...")

        # Send qdel command
        if qdel(job_id, quiet=True):
            success_count += 1
            # Update database status
            tracker.update_job_status(job_id, "cancelled")

            if wait:
                # Wait for deletion confirmation
                if _wait_for_deletion(job_id, timeout, quiet):
                    if not quiet:
                        click.echo(f"✅ Job {job_id} successfully cancelled")
                else:
                    if not quiet:
                        click.echo(
                            f"⚠️  Job {job_id} deletion timeout (may still be cancelling)"
                        )
            elif not quiet:
                click.echo(f"✅ Job {job_id} cancellation initiated")
        else:
            failed_jobs.append(job_id)

    # Summary
    if not quiet:
        if success_count > 0:
            click.echo(f"✅ Successfully cancelled {success_count} job(s)")
        if failed_jobs:
            click.echo(
                f"❌ Failed to cancel {len(failed_jobs)} job(s): {', '.join(failed_jobs)}"
            )

    return 0 if len(failed_jobs) == 0 else 1


@cancel_cli.command(name="name")
@click.argument("job_name")
@click.option("--pattern", "-p", is_flag=True, help="Treat job_name as a regex pattern")
@click.option("--force", "-f", is_flag=True, help="Skip confirmation prompt")
@click.option("--quiet", "-q", is_flag=True, help="Minimize output messages")
@click.option(
    "--wait", "-w", is_flag=True, help="Wait for job deletion to be confirmed"
)
@click.option(
    "--timeout", default=30, type=int, help="Maximum wait time in seconds (default: 30)"
)
def cancel_by_name(job_name, pattern, force, quiet, wait, timeout):
    """Cancel jobs by PBS job name or pattern."""
    # Find jobs matching the name/pattern
    matching_jobs = _find_jobs_by_name(job_name, pattern, quiet)

    if not matching_jobs:
        if not quiet:
            click.echo(f"❌ No jobs found matching '{job_name}'")
        return 1

    # Show matching jobs
    if not quiet:
        if pattern:
            click.echo(
                f"Found {len(matching_jobs)} job(s) matching pattern '{job_name}':"
            )
        else:
            click.echo(f"Found {len(matching_jobs)} job(s) with name '{job_name}':")
        for job_id, name in matching_jobs:
            click.echo(f"  • {job_id}: {name}")

    # Confirmation prompt
    if not force:
        job_list = ", ".join(job_id for job_id, _ in matching_jobs)
        if not click.confirm(f"Cancel {len(matching_jobs)} job(s): {job_list}?"):
            click.echo("❌ Cancellation aborted")
            return 0

    # Cancel jobs using the ID cancellation logic
    job_ids = [job_id for job_id, _ in matching_jobs]

    # Use click.Context to call the other command
    ctx = click.get_current_context()
    return ctx.invoke(
        cancel_by_id,
        job_ids=job_ids,
        force=True,
        quiet=quiet,
        wait=wait,
        timeout=timeout,
    )


@cancel_cli.command(name="status")
@click.option(
    "--status",
    "-s",
    type=click.Choice(["submitted", "running", "completed", "failed"]),
    help="Cancel jobs with specific status",
)
@click.option(
    "--dry-run",
    "-n",
    is_flag=True,
    help="Show what would be cancelled without doing it",
)
@click.option("--force", "-f", is_flag=True, help="Skip confirmation prompt")
@click.option("--quiet", "-q", is_flag=True, help="Minimize output messages")
@click.option("--limit", type=int, help="Maximum number of jobs to cancel")
def cancel_by_status(status, dry_run, force, quiet, limit):
    """Cancel qxub-tracked jobs by status."""
    from .resource_tracker import ResourceTracker

    if not status:
        click.echo("❌ Please specify --status", err=True)
        click.echo("Available statuses: submitted, running, completed, failed")
        return 1

    # Get jobs from qxub database
    tracker = ResourceTracker()
    jobs = tracker.get_jobs_by_status(status, limit)

    if not jobs:
        if not quiet:
            click.echo(f"❌ No jobs found with status '{status}'")
        return 1

    # Extract job IDs
    job_ids = [job["job_id"] for job in jobs]

    # Show matching jobs
    if not quiet:
        click.echo(f"Found {len(job_ids)} job(s) with status '{status}':")
        for job in jobs[:10]:  # Show first 10
            cmd_preview = (
                (job.get("command", "")[:50] + "...")
                if job.get("command", "")
                else "N/A"
            )
            click.echo(f"  • {job['job_id']}: {cmd_preview}")
        if len(jobs) > 10:
            click.echo(f"  ... and {len(jobs) - 10} more")

    # Dry run - just show what would be cancelled
    if dry_run:
        if quiet:
            # Just output job IDs for piping
            for job_id in job_ids:
                click.echo(job_id)
        else:
            click.echo(
                f"🔍 Dry run: Would cancel {len(job_ids)} job(s) with status '{status}'"
            )
        return 0

    # Use the existing cancel by ID functionality
    ctx = click.get_current_context()
    return ctx.invoke(
        cancel_by_id,
        job_ids=job_ids,
        dry_run=False,
        force=force,
        quiet=quiet,
        wait=False,
        timeout=30,
    )


@cancel_cli.command(name="filter")
@click.option(
    "--status",
    "-s",
    type=click.Choice(["submitted", "running", "completed", "failed"]),
    help="Filter by job status",
)
@click.option("--command", "-c", help="Filter by command pattern (substring match)")
@click.option("--job-name", help="Filter by job name pattern (regex)")
@click.option("--job-id", help="Filter by job ID pattern (regex)")
@click.option(
    "--before",
    help="Filter jobs submitted before time (e.g., '2h', '30m', '1d', '2023-10-22 14:30')",
)
@click.option(
    "--after",
    help="Filter jobs submitted after time (e.g., '2h', '30m', '1d', '2023-10-22 14:30')",
)
@click.option(
    "--dry-run",
    "-n",
    is_flag=True,
    help="Show what would be cancelled without doing it",
)
@click.option("--force", "-f", is_flag=True, help="Skip confirmation prompt")
@click.option("--quiet", "-q", is_flag=True, help="Minimize output messages")
@click.option("--limit", type=int, help="Maximum number of jobs to cancel")
def cancel_by_filter(
    status, command, job_name, job_id, before, after, dry_run, force, quiet, limit
):
    """Cancel qxub-tracked jobs using advanced filters."""
    import re
    from datetime import datetime, timedelta

    from .resource_tracker import ResourceTracker
    from .scheduler import get_job_resource_data

    # Get all jobs and apply filters
    tracker = ResourceTracker()
    jobs = tracker.get_jobs_by_status(limit=None)  # Get all jobs first

    if not jobs:
        if not quiet:
            click.echo("❌ No jobs found in qxub database")
        return 1

    # Apply filters
    filtered_jobs = jobs
    filter_descriptions = []

    # Status filter
    if status:
        filtered_jobs = [job for job in filtered_jobs if job.get("status") == status]
        filter_descriptions.append(f"Status: {status}")

    # Command filter
    if command:
        pattern = re.compile(re.escape(command), re.IGNORECASE)
        filtered_jobs = [
            job
            for job in filtered_jobs
            if job.get("command") and pattern.search(job["command"])
        ]
        filter_descriptions.append(f"Command contains: '{command}'")

    # Job name filter (need to get from PBS)
    if job_name:
        try:
            name_pattern = re.compile(job_name, re.IGNORECASE)
            jobs_with_names = []
            for job in filtered_jobs:
                job_data = get_job_resource_data(job["job_id"])
                if job_data and "metadata" in job_data:
                    pbs_job_name = job_data["metadata"].get("Job_Name", "")
                    if pbs_job_name and name_pattern.search(pbs_job_name):
                        jobs_with_names.append(job)
            filtered_jobs = jobs_with_names
            filter_descriptions.append(f"Job name matches: '{job_name}'")
        except re.error as e:
            if not quiet:
                click.echo(
                    f"❌ Invalid regex pattern for job name '{job_name}': {e}", err=True
                )
            return 1

    # Job ID filter
    if job_id:
        try:
            id_pattern = re.compile(job_id, re.IGNORECASE)
            filtered_jobs = [
                job for job in filtered_jobs if id_pattern.search(job["job_id"])
            ]
            filter_descriptions.append(f"Job ID matches: '{job_id}'")
        except re.error as e:
            if not quiet:
                click.echo(
                    f"❌ Invalid regex pattern for job ID '{job_id}': {e}", err=True
                )
            return 1

    # Enhanced time filters with better parsing
    if before or after:
        now = datetime.now()

        def parse_time(time_str):
            """Parse time string supporting hours, minutes, days, and datetime formats."""
            time_str = time_str.strip()

            # Try relative time formats first
            if time_str.endswith("m"):
                minutes = int(time_str[:-1])
                return now - timedelta(minutes=minutes)
            elif time_str.endswith("h"):
                hours = int(time_str[:-1])
                return now - timedelta(hours=hours)
            elif time_str.endswith("d"):
                days = int(time_str[:-1])
                return now - timedelta(days=days)
            else:
                # Try parsing as datetime
                # Support formats: "2023-10-22 14:30", "2023-10-22", "14:30"
                if ":" in time_str and "-" not in time_str:
                    # Just time, assume today
                    time_part = datetime.strptime(time_str, "%H:%M").time()
                    return datetime.combine(now.date(), time_part)
                elif " " in time_str:
                    # Full datetime
                    return datetime.strptime(time_str, "%Y-%m-%d %H:%M")
                else:
                    # Just date
                    return datetime.strptime(time_str, "%Y-%m-%d")

        try:
            if before:
                before_time = parse_time(before)
                filtered_jobs = [
                    job
                    for job in filtered_jobs
                    if job.get("submitted_at")
                    and datetime.fromisoformat(job["submitted_at"]) < before_time
                ]
                filter_descriptions.append(f"Submitted before: {before}")

            if after:
                after_time = parse_time(after)
                filtered_jobs = [
                    job
                    for job in filtered_jobs
                    if job.get("submitted_at")
                    and datetime.fromisoformat(job["submitted_at"]) > after_time
                ]
                filter_descriptions.append(f"Submitted after: {after}")
        except ValueError as e:
            if not quiet:
                click.echo(f"❌ Invalid time format: {e}", err=True)
                click.echo(
                    "Supported formats: '2h', '30m', '1d', '2023-10-22', '2023-10-22 14:30', '14:30'"
                )
            return 1

    # Apply limit
    if limit:
        filtered_jobs = filtered_jobs[:limit]

    if not filtered_jobs:
        if not quiet:
            click.echo("❌ No jobs match the specified filters")
        return 1

    # Extract job IDs
    job_ids = [job["job_id"] for job in filtered_jobs]

    # Show filtering results
    if not quiet:
        click.echo(f"Found {len(job_ids)} job(s) matching filters:")
        for desc in filter_descriptions:
            click.echo(f"  {desc}")
        click.echo("")

        for job in filtered_jobs[:10]:  # Show first 10
            cmd_preview = (
                (job.get("command", "")[:50] + "...")
                if job.get("command", "")
                else "N/A"
            )
            click.echo(f"  • {job['job_id']}: {cmd_preview}")
        if len(filtered_jobs) > 10:
            click.echo(f"  ... and {len(filtered_jobs) - 10} more")

    # Dry run - just show what would be cancelled
    if dry_run:
        if quiet:
            # Just output job IDs for piping
            for job_id in job_ids:
                click.echo(job_id)
        else:
            click.echo(f"🔍 Dry run: Would cancel {len(job_ids)} job(s)")
        return 0

    # Use the existing cancel by ID functionality
    ctx = click.get_current_context()
    return ctx.invoke(
        cancel_by_id,
        job_ids=job_ids,
        dry_run=False,
        force=force,
        quiet=quiet,
        wait=False,
        timeout=30,
    )


@cancel_cli.command(name="all")
@click.option(
    "--user", "-u", help="Cancel only jobs for specific user (default: current user)"
)
@click.option("--queue", "-q", help="Cancel only jobs in specific queue")
@click.option(
    "--status",
    "-s",
    type=click.Choice(["Q", "R", "H"]),
    help="Cancel only jobs with specific status (Q=queued, R=running, H=held)",
)
@click.option("--force", "-f", is_flag=True, help="Skip confirmation prompt")
@click.option("--quiet", is_flag=True, help="Minimize output messages")
def cancel_all(user, queue, status, force, quiet):
    """Cancel all jobs matching criteria (DANGEROUS - use with caution)."""
    import os

    # Default to current user for safety
    if not user:
        user = os.environ.get("USER", os.environ.get("USERNAME", ""))

    # Build qstat command
    qstat_cmd = ["qstat", "-u", user]
    if queue:
        qstat_cmd.extend(["-q", queue])

    # Get job list
    try:
        result = subprocess.run(qstat_cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split("\n")[2:]  # Skip header lines

        matching_jobs = []
        for line in lines:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 5:
                job_id = parts[0]
                job_status = parts[4]
                if status is None or job_status == status:
                    matching_jobs.append(job_id)

    except subprocess.CalledProcessError:
        click.echo(f"❌ Failed to get job list for user '{user}'", err=True)
        return 1

    if not matching_jobs:
        if not quiet:
            click.echo(
                f"❌ No jobs found for user '{user}'"
                + (f" in queue '{queue}'" if queue else "")
                + (f" with status '{status}'" if status else "")
            )
        return 1

    # Safety check and confirmation
    click.echo(
        f"⚠️  WARNING: This will cancel {len(matching_jobs)} job(s) for user '{user}'"
    )
    if queue:
        click.echo(f"   Queue filter: {queue}")
    if status:
        click.echo(f"   Status filter: {status}")
    click.echo(
        f"   Jobs: {', '.join(matching_jobs[:10])}"
        + (f" and {len(matching_jobs)-10} more..." if len(matching_jobs) > 10 else "")
    )

    if not force:
        if not click.confirm("Are you absolutely sure you want to cancel these jobs?"):
            click.echo("❌ Cancellation aborted")
            return 0

    # Cancel jobs using the ID cancellation logic
    ctx = click.get_current_context()
    return ctx.invoke(
        cancel_by_id,
        job_ids=matching_jobs,
        force=True,
        quiet=quiet,
        wait=False,
        timeout=30,
    )


def _find_jobs_by_name(
    job_name: str, pattern: bool = False, quiet: bool = False
) -> List[tuple]:
    """
    Find PBS jobs by name or pattern.

    Returns:
        List of (job_id, job_name) tuples
    """
    try:
        # Get all jobs for current user
        result = subprocess.run(
            ["qstat", "-f"], capture_output=True, text=True, check=True
        )

        jobs = []
        current_job_id = None
        current_job_name = None

        for line in result.stdout.split("\n"):
            line = line.strip()

            # Look for job ID line
            if line.startswith("Job Id:"):
                current_job_id = line.split(":", 1)[1].strip()
                current_job_name = None

            # Look for job name line
            elif line.startswith("Job_Name"):
                if "=" in line:
                    current_job_name = line.split("=", 1)[1].strip()

                    # Check if this job matches our criteria
                    if current_job_id and current_job_name:
                        if pattern:
                            try:
                                if re.search(job_name, current_job_name):
                                    jobs.append((current_job_id, current_job_name))
                            except re.error as e:
                                if not quiet:
                                    click.echo(
                                        f"❌ Invalid regex pattern '{job_name}': {e}",
                                        err=True,
                                    )
                                return []
                        else:
                            if current_job_name == job_name:
                                jobs.append((current_job_id, current_job_name))

        return jobs

    except subprocess.CalledProcessError as e:
        if not quiet:
            click.echo(f"❌ Failed to query jobs: {e}", err=True)
        return []


def _wait_for_deletion(job_id: str, timeout: int, quiet: bool = False) -> bool:
    """
    Wait for job deletion to be confirmed by checking job status.

    Returns:
        True if job is confirmed deleted, False if timeout
    """
    start_time = time.time()

    while time.time() - start_time < timeout:
        try:
            # Check if job still exists
            result = subprocess.run(
                ["qstat", job_id], capture_output=True, text=True, check=False
            )

            if result.returncode != 0:
                # Job not found - deletion confirmed
                return True

            # Job still exists, wait a bit
            time.sleep(1)

        except Exception as e:
            if not quiet:
                logging.debug(f"Error checking job status during deletion wait: {e}")
            break

    return False
