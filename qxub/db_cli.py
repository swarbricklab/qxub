"""``qxub db``: the resource statistics snapshot.

See ``qxub.queue.consolidate`` and ``docs/dev/sqlite-robustness.md``.
"""

from datetime import datetime
from pathlib import Path

import click


@click.group(name="db")
def db_cli():
    """Resource statistics snapshot built from run files."""


@db_cli.command()
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show what would be added without publishing or removing anything.",
)
def consolidate(dry_run):
    """Add finished runs to the statistics snapshot and tidy up run files.

    Safe to run at any time, from any login node: only one round runs at a
    time, and a round that doesn't get its turn does nothing. Useful in a
    Snakemake ``onsuccess:`` hook.
    """
    from .queue.consolidate import consolidate as run_round

    result = run_round(dry_run=dry_run)
    if dry_run and not result.skipped:
        click.echo("Dry run: nothing was published or removed.")
    click.echo(result.summary())
    for error in result.errors:
        click.echo(f"Error: {error}", err=True)
    if result.errors:
        raise SystemExit(1)


@db_cli.command()
@click.argument(
    "source",
    required=False,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show what would be added without publishing anything.",
)
def recover(source, dry_run):
    """Add the jobs in the old team job DB to the snapshot.

    SOURCE defaults to the old ``shared_db.path``. It is only read: a
    damaged file is copied and rebuilt with ``sqlite3 .recover``. Jobs the
    snapshot already has are kept as they are, so it is safe to run again.
    """
    from .queue.db import legacy_shared_db_path
    from .queue.recover import recover as run_recover
    from .queue.recover import source_problem

    source = source or legacy_shared_db_path()
    if source is None:
        raise click.UsageError(
            "No SOURCE given and shared_db.path isn't set; name the old DB file."
        )
    problem = source_problem(source)
    if problem:
        click.echo(f"Warning: {problem}", err=True)
    result = run_recover(source, dry_run=dry_run)
    if dry_run and not result.skipped:
        click.echo("Dry run: nothing was published.")
    click.echo(result.summary())
    for error in result.errors:
        click.echo(f"Error: {error}", err=True)
    if result.errors or result.skipped:
        raise SystemExit(1)


@db_cli.command()
def status():
    """Show the snapshot's age and size, and run files waiting for it."""
    from .queue import runs
    from .queue.consolidate import (
        ATTEMPT_MARKER,
        Lease,
        lease_path,
        open_snapshot,
        stats_db_path,
        writable_problem,
    )
    from .queue.db import legacy_shared_db_path
    from .queue.recover import recovered_from

    snapshot = stats_db_path()
    click.echo(f"Snapshot:   {snapshot}")
    problem = writable_problem(snapshot)
    if problem:
        click.echo(f"Warning:    cannot publish: {problem}")
    if snapshot.exists():
        conn = open_snapshot(snapshot)
        try:
            row = conn.execute(
                "SELECT value FROM snapshot_meta WHERE key='published_at'"
            ).fetchone()
            jobs = conn.execute("SELECT COUNT(*) FROM job_resources").fetchone()[0]
            with_usage = conn.execute(
                "SELECT COUNT(*) FROM job_resources WHERE mem_used_mb IS NOT NULL"
            ).fetchone()[0]
        finally:
            conn.close()
        published = row[0] if row else None
        click.echo(f"Published:  {published or 'unknown'}{_age(published)}")
        click.echo(f"Jobs:       {jobs} ({with_usage} with resource usage)")
    else:
        click.echo("Published:  never")

    legacy = legacy_shared_db_path()
    if legacy and legacy.exists():
        done = recovered_from(snapshot)
        click.echo(
            f"Old job DB: {legacy} (no longer written; "
            + (f"recovered: {done})" if done else "run `qxub db recover` once)")
        )

    runs_dir = runs.runs_dir()
    count = len(list(runs_dir.glob(f"*{runs.SUFFIX}"))) if runs_dir.is_dir() else 0
    click.echo(f"Run files:  {count} in {runs_dir}")
    marker = runs_dir / ATTEMPT_MARKER
    if marker.exists():
        last = datetime.fromtimestamp(marker.stat().st_mtime)
        click.echo(
            f"Last round: {last.isoformat(timespec='seconds')}"
            f"{_age(last.isoformat())}"
        )
    holder = Lease(lease_path(snapshot)).holder()
    if holder:
        click.echo(
            "Running:    {}@{} (pid {}), {}s".format(
                holder.get("user", "?"),
                holder.get("host", "?"),
                holder.get("pid", "?"),
                holder.get("age_seconds", "?"),
            )
        )


def _age(iso):
    try:
        delta = datetime.now() - datetime.fromisoformat(iso)
    except (TypeError, ValueError):
        return ""
    hours = delta.total_seconds() / 3600
    if hours < 48:
        return f" ({hours:.1f} h ago)"
    return f" ({hours / 24:.1f} days ago)"
