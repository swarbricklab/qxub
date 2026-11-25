"""
Cancel CLI for qxub - Cancel PBS jobs with qxub database integration.

Provides job cancellation functionality for snakemake compatibility and
general job management. Supports cancellation by job ID, name, or pattern
matching with automatic database status updates.
"""

import logging

logger = logging.getLogger(__name__)
import json
import re
import time
from collections import Counter
from datetime import datetime, timedelta
from typing import List

import click

from . import backends
from .backends import SchedulerQueryError
from .resources import ResourceTracker


@click.group(name="cancel")
def cancel_cli():
    """Cancel PBS jobs with database tracking."""


@cancel_cli.result_callback()
def _exit_with_subcommand_status(result, **_kwargs):
    """Use a subcommand's integer return value as the exit status.

    Click ignores the values that commands return, so without this a
    subcommand's ``return 1`` (for example after a failed ``qdel``) would
    still exit 0.
    """
    if isinstance(result, int) and result:
        click.get_current_context().exit(result)


UNFINISHED = ("initiated", "dispatched", "submitted", "running")
# Leave rows this young alone even if qstat doesn't list them yet.
STALE_AFTER = timedelta(minutes=10)


def _still_queued_or_running(
    tracker, jobs: List[dict], quiet: bool, dry_run: bool = False
):
    """The database rows whose job the scheduler still lists.

    Only those can be cancelled, so ``cancel status``/``filter`` never pick a
    row left ``submitted`` long ago. One scheduler query
    (``list_jobs`` on the scheduler backend) for all rows. A row whose job name wasn't
    stored gets the (possibly truncated) name from that listing.

    Unfinished rows that the scheduler no longer lists, and that are older
    than :data:`STALE_AFTER`, are finished: from the statistics snapshot if
    it knows how the job ended, else as ``unknown``. No run file said how
    they ended, and PBS has forgotten them. A dry run only reports them.

    Returns ``None`` (after saying why) if the scheduler can't be asked.
    """
    try:
        live = {
            (owner.name, owner.same_job(j.job_id)): j
            for owner, j in backends.live_jobs()
        }
    except SchedulerQueryError as exc:
        click.echo(
            f"❌ Can't tell which jobs are still queued or running: {exc}", err=True
        )
        return None

    cutoff = (datetime.now() - STALE_AFTER).isoformat()
    candidates, gone = [], []
    current = backends.get_backend()
    for job in jobs:
        backend = backends.backend_for_job(job["job_id"], job.get("backend"))
        key = (backend.name, backend.same_job(job["job_id"]))
        if key in live:
            if not job.get("job_name"):
                job = {**job, "job_name": live[key].name}
            candidates.append(job)
        elif (
            job.get("status") in UNFINISHED
            and backend.looks_like_job_id(job["job_id"])
            # Another backend's jobs count as gone only if it can be asked
            and (backend.name == current.name or backend.available())
            and (job.get("submitted_at") or "") < cutoff
        ):
            gone.append(job["job_id"])
    if gone and dry_run:
        if not quiet:
            click.echo(
                f"ℹ️  {len(gone)} old job(s) are still 'submitted' or 'running' in"
                " the database but no longer queued or running (left as they are"
                " in a dry run)"
            )
    elif gone:
        counts = tracker.settle_gone_jobs(gone)
        if not quiet and any(counts.values()):
            click.echo(
                f"ℹ️  Finished {sum(counts.values())} old job(s) that were still"
                " 'submitted' or 'running' in the database but no longer queued"
                f" or running: {counts['from_snapshot']} from the stats snapshot,"
                f" {counts['unknown']} as unknown"
            )
    return candidates


def _tags_match(job: dict, wanted: List[str]) -> bool:
    """True if the job has every ``KEY=VALUE`` tag in ``wanted``."""
    try:
        tags = set(json.loads(job.get("tags") or "[]"))
    except (TypeError, ValueError):
        return False
    return all(tag in tags for tag in wanted)


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

        # The job's own backend: the one recorded at submission if known
        backend = backends.backend_for_job(job_id, tracker.recorded_backend(job_id))
        ok, reason = backend.cancel(job_id)
        if ok:
            success_count += 1
            for db_id in _db_job_ids(tracker, backend, job_id):
                tracker.update_job_status(db_id, "cancelled")

            if wait:
                # Wait for deletion confirmation
                if _wait_for_deletion(job_id, timeout, quiet, backend):
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
            if not quiet:
                click.echo(f"❌ Could not cancel {job_id}: {reason}", err=True)

    # Summary
    if not quiet:
        if success_count > 0:
            click.echo(f"✅ Successfully cancelled {success_count} job(s)")
        if failed_jobs:
            click.echo(
                f"❌ Failed to cancel {len(failed_jobs)} job(s): {', '.join(failed_jobs)}"
            )

    return 0 if len(failed_jobs) == 0 else 1


def _db_job_ids(tracker, backend, job_id: str) -> List[str]:
    """The DB rows for ``job_id``, under any spelling of it.

    ``cancel id 12345`` must also mark the stored ``12345.gadi-pbs`` row.
    Falls back to ``job_id`` itself if the DB has no matching row.
    """
    key = backend.same_job(job_id)
    found = [
        db_id
        for db_id in tracker.job_ids_starting_with(key)
        if backend.same_job(db_id) == key
    ]
    return found or [job_id]


@cancel_cli.command(name="name")
@click.argument("job_name")
@click.option("--pattern", "-p", is_flag=True, help="Treat job_name as a regex pattern")
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
def cancel_by_name(job_name, pattern, dry_run, force, quiet, wait, timeout):
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

    if dry_run:
        if quiet:
            for job_id, _ in matching_jobs:
                click.echo(job_id)
        else:
            click.echo(f"🔍 Dry run: Would cancel {len(matching_jobs)} job(s)")
        return 0

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
    from .resources import ResourceTracker

    if not status:
        click.echo("❌ Please specify --status", err=True)
        click.echo("Available statuses: submitted, running, completed, failed")
        return 1

    # Get jobs from qxub database
    tracker = ResourceTracker()
    tracker.sync_from_run_files()
    jobs = _still_queued_or_running(
        tracker, tracker.get_jobs_by_status(), quiet, dry_run
    )
    if jobs is None:
        return 1
    jobs = [job for job in jobs if job.get("status") == status]
    if limit:
        jobs = jobs[:limit]

    if not jobs:
        if not quiet:
            click.echo(f"❌ No queued or running jobs found with status '{status}'")
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
    "--tag",
    "tags",
    multiple=True,
    metavar="KEY=VALUE",
    help="Only jobs with this tag (from exec --tag/--tags); repeat to require several",
)
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
    status,
    command,
    job_name,
    job_id,
    tags,
    before,
    after,
    dry_run,
    force,
    quiet,
    limit,
):
    """Cancel qxub-tracked jobs that are still queued or running, by filters.

    Every filter works on the qxub database and one listing of your jobs
    from the scheduler; no per-job scheduler queries.
    """
    tracker = ResourceTracker()
    tracker.sync_from_run_files()
    jobs = _still_queued_or_running(
        tracker, tracker.get_jobs_by_status(), quiet, dry_run
    )
    if jobs is None:
        return 1

    if not jobs:
        if not quiet:
            click.echo("❌ No qxub jobs are queued or running")
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

    # Tag filter
    if tags:
        bad = [t for t in tags if "=" not in t]
        if bad:
            click.echo(f"❌ Tags must be KEY=VALUE: {', '.join(bad)}", err=True)
            return 1
        filtered_jobs = [job for job in filtered_jobs if _tags_match(job, list(tags))]
        filter_descriptions.append(f"Tags: {', '.join(tags)}")

    # Job name filter: the name stored at submission, else the scheduler's
    if job_name:
        try:
            name_pattern = re.compile(job_name, re.IGNORECASE)
        except re.error as e:
            if not quiet:
                click.echo(
                    f"❌ Invalid regex pattern for job name '{job_name}': {e}", err=True
                )
            return 1
        filtered_jobs = [
            job
            for job in filtered_jobs
            if job.get("job_name") and name_pattern.search(job["job_name"])
        ]
        filter_descriptions.append(f"Job name matches: '{job_name}'")

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
            name = f" [{job['job_name']}]" if job.get("job_name") else ""
            click.echo(f"  • {job['job_id']}{name}: {cmd_preview}")
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
@click.option(
    "--all-backends",
    is_flag=True,
    help="Also cancel jobs of other schedulers on this machine (local jobs "
    "under a cluster platform, cluster jobs under the local one)",
)
@click.option("--force", "-f", is_flag=True, help="Skip confirmation prompt")
@click.option("--quiet", is_flag=True, help="Minimize output messages")
def cancel_all(
    user, queue, status, all_backends, force, quiet
):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    """Cancel all jobs matching criteria (DANGEROUS - use with caution)."""
    import os

    # Default to current user for safety
    if not user:
        user = os.environ.get("USER", os.environ.get("USERNAME", ""))

    # The scheduler lists the queue a job is *in*, which on sites with
    # routing queues differs from the one it was submitted to (Gadi: `copyq`
    # -> `copyq-exec`). Accept either so `--queue copyq` works as users expect.
    accepted_queues = {queue, f"{queue}-exec"} if queue else None

    # Only the current platform's jobs unless asked: trying jobs out on the
    # local platform, `cancel all -f` shouldn't kill your cluster jobs.
    try:
        if all_backends:
            live = backends.live_jobs(user)
        else:
            current = backends.get_backend()
            live = [(current, job) for job in current.list_jobs(user)]
        matching = [
            (backend, job.job_id)
            for backend, job in live
            if (accepted_queues is None or job.queue in accepted_queues)
            and (status is None or job.state == status)
        ]
    except SchedulerQueryError as exc:
        click.echo(f"❌ Failed to get job list for user '{user}': {exc}", err=True)
        return 1

    matching_jobs = [job_id for _backend, job_id in matching]
    by_backend = Counter(backend.name for backend, _job_id in matching)

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
        # Which scheduler's, when there is more than one
        + (
            " (" + ", ".join(f"{n} {name}" for name, n in by_backend.items()) + ")"
            if len(by_backend) > 1
            else ""
        )
    )
    if queue:
        click.echo(f"   Queue filter: {queue}")
    if status:
        click.echo(f"   Status filter: {status}")
    click.echo(
        f"   Jobs: {', '.join(matching_jobs[:10])}"
        + (f" and {len(matching_jobs) - 10} more..." if len(matching_jobs) > 10 else "")
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
    Find your queued or running jobs by name or pattern.

    Asks the scheduler for whole names (PBS: ``qstat -f``), because the
    short listing truncates them.

    Returns:
        List of (job_id, job_name) tuples
    """
    if pattern:
        try:
            regex = re.compile(job_name)
        except re.error as e:
            if not quiet:
                click.echo(f"❌ Invalid regex pattern '{job_name}': {e}", err=True)
            return []

    try:
        live = [job for _backend, job in backends.live_jobs(full_names=True)]
    except SchedulerQueryError as e:
        if not quiet:
            click.echo(f"❌ Failed to query jobs: {e}", err=True)
        return []

    return [
        (job.job_id, job.name)
        for job in live
        if job.name and (regex.search(job.name) if pattern else job.name == job_name)
    ]


def _wait_for_deletion(
    job_id: str, timeout: int, quiet: bool = False, backend=None
) -> bool:
    """
    Wait until the scheduler no longer lists the job.

    Returns:
        True if job is confirmed deleted, False if timeout
    """
    backend = backend or backends.backend_for_job(job_id)
    key = backend.same_job(job_id)
    start_time = time.time()

    while time.time() - start_time < timeout:
        try:
            if all(backend.same_job(j.job_id) != key for j in backend.list_jobs()):
                return True
        except SchedulerQueryError as e:
            if not quiet:
                logger.debug(f"Error checking job status during deletion wait: {e}")
            break
        time.sleep(1)

    return False
