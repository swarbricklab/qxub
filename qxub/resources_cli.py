"""
CLI commands for resource efficiency tracking and analysis.
"""

import csv as csv_mod
import json
import sys
from datetime import datetime
from pathlib import Path

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .resources.tracker import SnapshotReader, parse_date_filter

console = Console()


@click.group()
def resources():
    """View resource usage and efficiency statistics.

    Reports read the published statistics snapshot (see ``qxub db status``),
    which holds finished jobs. Run ``qxub db consolidate`` to add recent ones.
    """
    pass


def _snapshot() -> SnapshotReader:
    """The snapshot to report on; exits with a hint if there isn't one yet."""
    from .queue.consolidate import stats_db_path

    path = stats_db_path()
    if not path.exists():
        console.print(
            f"📊 No statistics snapshot yet at {path}.\n"
            "Run `qxub db consolidate` to build it from finished jobs.",
            style="yellow",
        )
        raise SystemExit(0)
    return SnapshotReader(path)


@resources.command()
@click.option("--limit", default=20, help="Number of recent jobs to show")
@click.option(
    "--tag",
    multiple=True,
    metavar="TAG",
    help="Filter by tag (can be used multiple times, e.g. --tag rule=align)",
)
@click.option("--user", default=None, help="Filter by username")
@click.option(
    "--since",
    default=None,
    metavar="DATE",
    help="Show jobs submitted after DATE (e.g. '7d', '2h', '1w', '2026-02-20')",
)
@click.option(
    "--before",
    default=None,
    metavar="DATE",
    help="Show jobs submitted before DATE (e.g. '2026-02-25')",
)
@click.option(
    "--queue",
    "queue_name",
    default=None,
    help="Filter by queue name (e.g. normal, hugemem)",
)
@click.option(
    "--status",
    default=None,
    type=click.Choice(
        ["submitted", "running", "completed", "failed", "cancelled"],
        case_sensitive=False,
    ),
    help="Filter by job status",
)
@click.option(
    "--csv",
    "output_csv",
    is_flag=True,
    help="Output as CSV to stdout (redirect to file with > output.csv)",
)
def list(limit, tag, user, since, before, queue_name, status, output_csv):
    """List recent finished jobs with resource efficiency."""
    # Parse date filter args, emitting a clear error for bad formats.
    since_dt = before_dt = None
    for label, raw in (("--since", since), ("--before", before)):
        if raw:
            try:
                parsed = parse_date_filter(raw)
            except ValueError as exc:
                raise click.BadParameter(str(exc), param_hint=label) from exc
            if label == "--since":
                since_dt = parsed
            else:
                before_dt = parsed

    jobs = _snapshot().get_recent_jobs(
        limit,
        tags=tag if tag else None,
        user=user,
        since=since_dt,
        before=before_dt,
        queue_name=queue_name,
        status=status,
    )

    if not jobs:
        if not output_csv:
            console.print("📊 No resource data found", style="yellow")
        return

    if output_csv:
        writer = csv_mod.writer(sys.stdout)
        writer.writerow(
            [
                "job_id",
                "username",
                "status",
                "submitted_at",
                "command",
                "exit_code",
                "queue",
                "mem_requested_mb",
                "mem_used_mb",
                "mem_efficiency",
                "time_requested_sec",
                "time_used_sec",
                "time_efficiency",
                "cpus_requested",
                "cpu_efficiency",
                "tags",
            ]
        )
        for job in jobs:
            raw_tags = job.get("tags") or "[]"
            try:
                tags_display = ", ".join(json.loads(raw_tags))
            except Exception:
                tags_display = ""
            writer.writerow(
                [
                    job.get("job_id", ""),
                    job.get("username", ""),
                    job.get("status", ""),
                    job.get("submitted_at") or job.get("timestamp", ""),
                    job.get("command", ""),
                    job.get("exit_code", ""),
                    job.get("queue", ""),
                    job.get("mem_requested_mb", ""),
                    job.get("mem_used_mb", ""),
                    job.get("mem_efficiency", ""),
                    job.get("time_requested_sec", ""),
                    job.get("time_used_sec", ""),
                    job.get("time_efficiency", ""),
                    job.get("cpus_requested", ""),
                    job.get("cpu_efficiency", ""),
                    tags_display,
                ]
            )
        return

    table = Table(title="Recent Jobs - Resource Efficiency")
    table.add_column("Submitted", style="dim", no_wrap=True, width=12)
    table.add_column("User", style="dim", no_wrap=True, width=10)
    table.add_column("Job ID", style="cyan", no_wrap=True, width=10)
    table.add_column("Status", no_wrap=True, width=5)
    table.add_column("Command", style="white", width=26)
    table.add_column("Exit", style="blue", no_wrap=True, width=6)
    table.add_column("Mem", style="green", no_wrap=True, width=5)
    table.add_column("Time", style="yellow", no_wrap=True, width=6)
    table.add_column("CPU", style="magenta", no_wrap=True, width=5)
    table.add_column("Queue", style="dim", no_wrap=True, width=12)
    table.add_column("Tags", style="cyan", width=22)

    for job in jobs:
        # Format job ID (show first 8 chars)
        job_id = job["job_id"][:8] if job["job_id"] else "unknown"

        # Format submitted_at timestamp
        raw_ts = job.get("submitted_at") or job.get("timestamp") or ""
        try:
            submitted_str = datetime.fromisoformat(raw_ts).strftime("%m-%d %H:%M")
        except Exception:
            submitted_str = raw_ts[:10] if raw_ts else "?"

        # Format tags
        raw_tags = job.get("tags") or "[]"
        try:
            tags_display = " ".join(json.loads(raw_tags))
        except Exception:
            tags_display = ""
        # Format command (smart truncation)
        command = job["command"] or ""
        if len(command) > 35:
            # Try to preserve the important parts - keep the subcommand and first part of actual command
            if " -- " in command:
                # Split on -- separator
                parts = command.split(" -- ", 1)
                cmd_part = parts[1] if len(parts) > 1 else parts[0]
                if len(cmd_part) > 25:
                    command = cmd_part[:22] + "..."
                else:
                    command = cmd_part
            else:
                # Find the actual command part (after subcommand like 'conda --env base')
                words = command.split()
                if len(words) > 4:  # qxub <subcommand> --option value <actual_command>
                    actual_cmd_start = 4
                    for i, word in enumerate(words[1:], 1):  # Skip 'qxub'
                        if not word.startswith("-") and i > 2:  # Found actual command
                            actual_cmd_start = i
                            break
                    if actual_cmd_start < len(words):
                        actual_cmd = " ".join(words[actual_cmd_start:])
                        if len(actual_cmd) > 25:
                            command = actual_cmd[:22] + "..."
                        else:
                            command = actual_cmd
                    else:
                        command = command[:32] + "..."
                else:
                    command = command[:32] + "..."

        # Format status
        status = job.get("status") or "unknown"
        status_map = {
            "submitted": "[yellow]sub[/yellow]",
            "running": "[cyan]run[/cyan]",
            "completed": "[green]done[/green]",
            "failed": "[red]fail[/red]",
            "cancelled": "[dim]canc[/dim]",
        }
        status_str = status_map.get(status, f"[dim]{status[:4]}[/dim]")

        # Format exit code
        exit_code = job["exit_code"]
        if exit_code == 0:
            exit_style = "[green]✓[/green]"
        elif exit_code is None:
            exit_style = "[yellow]?[/yellow]"
        else:
            exit_style = f"[red]{exit_code}[/red]"

        # Format efficiency with color coding
        def format_efficiency(eff, unit=""):
            if eff is None:
                return "[dim]?[/dim]"
            color = "green" if eff >= 70 else "yellow" if eff >= 50 else "red"
            return f"[{color}]{eff:.0f}%{unit}[/{color}]"

        mem_eff = format_efficiency(job["mem_efficiency"])
        time_eff = format_efficiency(job["time_efficiency"])
        cpu_eff = format_efficiency(job["cpu_efficiency"])

        queue = job["queue"] or "?"

        username = job.get("username") or "?"

        table.add_row(
            submitted_str,
            username,
            job_id,
            status_str,
            command,
            exit_style,
            mem_eff,
            time_eff,
            cpu_eff,
            queue,
            tags_display,
        )

    console.print(table)


@resources.command()
def stats():
    """Show overall resource efficiency statistics."""
    stats = _snapshot().get_efficiency_stats()

    if not stats or stats.get("total_jobs", 0) == 0:
        console.print("📊 No resource data available for statistics", style="yellow")
        return

    # Create efficiency panel
    efficiency_text = f"""
[bold]Overall Efficiency Averages:[/bold]
  Memory:  {stats['avg_mem_efficiency']:.1f}%
  Time:    {stats['avg_time_efficiency']:.1f}%
  CPU:     {stats['avg_cpu_efficiency']:.1f}%
  Jobfs:   {stats['avg_jobfs_efficiency']:.1f}%

[bold]Low Efficiency Jobs (< 50%):[/bold]
  Memory:  {stats['low_mem_efficiency_jobs']} jobs
  Time:    {stats['low_time_efficiency_jobs']} jobs
  CPU:     {stats['low_cpu_efficiency_jobs']} jobs
"""

    console.print(
        Panel(
            efficiency_text,
            title=f"📊 Resource Statistics ({stats['total_jobs']} jobs)",
            border_style="blue",
        )
    )


@resources.command()
@click.option("--threshold", default=50.0, help="Efficiency threshold (default: 50%)")
@click.option("--limit", default=10, help="Number of jobs to show")
def inefficient(threshold, limit):
    """Show jobs with low resource efficiency."""
    jobs = _snapshot().get_inefficient_jobs(threshold, limit)

    if not jobs:
        console.print(
            f"✅ No jobs found with efficiency below {threshold}%", style="green"
        )
        return

    table = Table(title=f"Inefficient Jobs (< {threshold}% efficiency)")
    table.add_column("Job ID", style="cyan", no_wrap=True)
    table.add_column("Command", style="white")
    table.add_column("Memory", style="red")
    table.add_column("Time", style="yellow")
    table.add_column("CPU", style="magenta")
    table.add_column("Date", style="dim")

    for job in jobs:
        job_id = job["job_id"][:8] if job["job_id"] else "unknown"

        command = job["command"] or ""
        if len(command) > 30:
            command = command[:27] + "..."

        # Format resource usage with requested vs used
        def format_resource_usage(eff, used, requested, unit):
            if eff is None or used is None or requested is None:
                return "[dim]?[/dim]"
            color = "red" if eff < 25 else "yellow" if eff < 50 else "green"
            return f"[{color}]{eff:.0f}%[/{color}] ({used}{unit}/{requested}{unit})"

        mem_usage = format_resource_usage(
            job["mem_efficiency"], job["mem_used_mb"], job["mem_requested_mb"], "MB"
        )

        time_usage = format_resource_usage(
            job["time_efficiency"],
            job["time_used_sec"] // 60 if job["time_used_sec"] else None,
            job["time_requested_sec"] // 60 if job["time_requested_sec"] else None,
            "min",
        )

        cpu_eff = job["cpu_efficiency"]
        cpu_usage = (
            f"[red]{cpu_eff:.0f}%[/red]"
            if cpu_eff and cpu_eff < threshold
            else f"{cpu_eff:.0f}%" if cpu_eff else "[dim]?[/dim]"
        )

        # Format date
        try:
            from datetime import datetime

            dt = datetime.fromisoformat(job["timestamp"])
            date_str = dt.strftime("%m-%d %H:%M")
        except Exception:
            date_str = "?"

        table.add_row(job_id, command, mem_usage, time_usage, cpu_usage, date_str)

    console.print(table)

    # Show recommendations
    console.print("\n💡 [bold]Optimization Tips:[/bold]")
    console.print(
        "• Low memory efficiency: Request less memory or use memory-optimized algorithms"
    )
    console.print(
        "• Low time efficiency: Request shorter walltime or optimize your code"
    )
    console.print(
        "• Low CPU efficiency: Use more CPU cores or optimize for parallelization"
    )


@resources.command()
@click.option("--days", default=7, help="Number of days to show trends for")
def trends(days):
    """Show resource usage trends over time."""
    trend_data = _snapshot().get_resource_trends(days)

    if not trend_data:
        console.print(
            f"📊 No resource data available for the last {days} days", style="yellow"
        )
        return

    table = Table(title=f"Resource Usage Trends (Last {days} days)")
    table.add_column("Date", style="cyan")
    table.add_column("Jobs", style="blue", no_wrap=True)
    table.add_column("Avg Mem Eff", style="green", no_wrap=True)
    table.add_column("Avg Time Eff", style="yellow", no_wrap=True)
    table.add_column("Avg CPU Eff", style="magenta", no_wrap=True)
    table.add_column("Total Mem", style="dim", no_wrap=True)

    for day in trend_data:
        date = day["date"]
        job_count = day["job_count"]

        # Format efficiencies
        mem_eff = f"{day['avg_mem_eff']:.0f}%" if day["avg_mem_eff"] else "?"
        time_eff = f"{day['avg_time_eff']:.0f}%" if day["avg_time_eff"] else "?"
        cpu_eff = f"{day['avg_cpu_eff']:.0f}%" if day["avg_cpu_eff"] else "?"

        # Format total memory
        total_mem_mb = day["total_mem_req_mb"] or 0
        total_mem_gb = total_mem_mb / 1024
        mem_str = f"{total_mem_gb:.1f}GB" if total_mem_gb > 0 else "?"

        table.add_row(date, str(job_count), mem_eff, time_eff, cpu_eff, mem_str)

    console.print(table)


@resources.command()
@click.argument("output_path")
@click.option("--limit", help="Maximum number of records to export")
def export(output_path, limit):
    """Export resource data to CSV file."""
    reader = _snapshot()
    try:
        output_file = Path(output_path)
        reader.export_csv(output_file, limit)
        console.print(f"✅ Exported resource data to {output_path}", style="green")

        # Show summary
        job_count = len(reader.get_recent_jobs(limit or 1000000))
        console.print(f"📊 Exported {job_count} job records")

    except Exception as e:
        console.print(f"❌ Failed to export data: {e}", style="red")


@resources.command()
@click.argument("job_id")
def show(job_id):
    """Show detailed resource information for a specific job."""
    from .core.scheduler import get_job_info
    from .resources import bytes_to_human
    from .resources.spec import format_duration

    def size(value):
        return bytes_to_human(value) if value is not None else "?"

    def duration(value):
        return format_duration(value) if value is not None else "?"

    def percent(used, requested):
        if used is None or not requested:
            return "?"
        return f"{used / requested * 100:.1f}%"

    try:
        info = get_job_info(job_id)

        if not info:
            console.print(f"❌ No resource data found for job {job_id}", style="red")
            return

        console.print(f"📊 [bold]Resource Details for Job {job_id}[/bold]\n")

        ended = f" ({info.ended_by})" if info.ended_by else ""
        console.print(f"State: {info.state.value}{ended}")
        console.print(
            "Exit Status: "
            + ("?" if info.exit_status is None else str(info.exit_status))
        )
        console.print(f"Queue: {info.queue or '?'}")
        console.print(f"Host: {info.host or '?'}")

        usage = info.usage
        table = Table(title="Resource Usage Comparison")
        table.add_column("Resource", style="cyan")
        table.add_column("Requested", style="blue")
        table.add_column("Used", style="green")
        table.add_column("Efficiency", style="magenta")
        table.add_row(
            "Memory",
            size(usage.memory_requested),
            size(usage.memory_used),
            percent(usage.memory_used, usage.memory_requested),
        )
        table.add_row(
            "Walltime",
            duration(usage.walltime_requested),
            duration(usage.walltime_used),
            percent(usage.walltime_used, usage.walltime_requested),
        )
        cpus = usage.cpus_used or usage.cpus_requested
        table.add_row(
            "CPU",
            f"{usage.cpus_requested or '?'} cores",
            f"{duration(usage.cpu_time)} CPU time",
            percent(usage.cpu_time, (usage.walltime_used or 0) * (cpus or 0)),
        )
        table.add_row(
            "Job Filesystem",
            size(usage.disk_requested),
            size(usage.disk_used),
            percent(usage.disk_used, usage.disk_requested),
        )
        console.print(table)

        if info.queue_wait is not None:
            console.print("\n⏱️  [bold]Timing:[/bold]")
            console.print(f"Queue wait: {info.queue_wait} seconds")
            if info.run_time is not None:
                console.print(f"Execution: {info.run_time} seconds")

    except Exception as e:
        console.print(f"❌ Error retrieving resource data: {e}", style="red")
