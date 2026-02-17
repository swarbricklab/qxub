"""
Replay command for qxub interactive session recordings.

This module provides the 'replay' subcommand that allows users to:
- List available session recordings
- View transcript contents
- Replay sessions with timing data
"""

import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import click


def _get_config_manager():
    """Get the current config manager instance."""
    from .config import manager

    return manager.config_manager


def _get_default_record_dir() -> Path | None:
    """Get default transcript recording directory from config."""
    try:
        config_mgr = _get_config_manager()
        record_dir = config_mgr.get_config_value("defaults.interactive.record_dir")
        if record_dir:
            # Resolve template variables
            template_vars = {
                "user": os.environ.get("USER", "unknown"),
                "project": os.environ.get("PROJECT", "unknown"),
            }
            resolved = config_mgr.resolve_templates(record_dir, template_vars)
            return Path(resolved)
        return None
    except Exception:
        return None


def _find_recordings(directory: Path | None = None) -> list[dict]:
    """
    Find all session recordings in the given directory.

    Returns a list of dicts with:
        - transcript: Path to transcript file
        - timing: Path to timing file (or None)
        - timestamp: datetime parsed from filename
        - duration: Approximate duration in seconds (from timing file)
    """
    recordings = []

    if directory is None:
        directory = _get_default_record_dir()

    if directory is None or not directory.exists():
        return recordings

    # Find transcript files matching the pattern
    for transcript in sorted(directory.glob("qxub-transcript-*.txt"), reverse=True):
        # Extract timestamp from filename: qxub-transcript-YYYYMMDD-HHMMSS.txt
        try:
            name = transcript.stem  # qxub-transcript-YYYYMMDD-HHMMSS
            date_part = name.replace("qxub-transcript-", "")
            timestamp = datetime.strptime(date_part, "%Y%m%d-%H%M%S")
        except ValueError:
            timestamp = None

        # Find corresponding timing file
        timing_name = transcript.name.replace("qxub-transcript-", "qxub-timing-")
        timing = transcript.parent / timing_name

        if not timing.exists():
            timing = None

        # Estimate duration from timing file if available
        duration = None
        if timing and timing.exists():
            try:
                with open(timing) as f:
                    total = 0.0
                    for line in f:
                        parts = line.strip().split()
                        if parts:
                            try:
                                total += float(parts[0])
                            except ValueError:
                                pass
                    duration = int(total)
            except Exception:
                pass

        recordings.append(
            {
                "transcript": transcript,
                "timing": timing,
                "timestamp": timestamp,
                "duration": duration,
            }
        )

    return recordings


def _format_duration(seconds: int | None) -> str:
    """Format duration in human-readable form."""
    if seconds is None:
        return "unknown"
    if seconds < 60:
        return f"{seconds}s"
    elif seconds < 3600:
        mins = seconds // 60
        secs = seconds % 60
        return f"{mins}m{secs}s"
    else:
        hours = seconds // 3600
        mins = (seconds % 3600) // 60
        return f"{hours}h{mins}m"


@click.command("replay")
@click.argument("session", required=False)
@click.option(
    "--list",
    "-l",
    "list_recordings",
    is_flag=True,
    help="List available session recordings",
)
@click.option(
    "--view",
    "-v",
    is_flag=True,
    help="View transcript content instead of replaying",
)
@click.option(
    "--speed",
    "-s",
    type=float,
    default=1.0,
    help="Playback speed multiplier (e.g., 2 for 2x speed)",
)
@click.option(
    "--max-delay",
    "-m",
    type=float,
    default=2.0,
    help="Maximum delay between events in seconds",
)
@click.option(
    "--dir",
    "-d",
    "record_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Directory containing recordings (defaults to config record_dir)",
)
def replay_cli(
    session: str | None,
    list_recordings: bool,
    view: bool,
    speed: float,
    max_delay: float,
    record_dir: Path | None,
):
    """
    Replay interactive session recordings.

    SESSION can be:
    - A number (1 = most recent, 2 = second most recent, etc.)
    - A timestamp (YYYYMMDD-HHMMSS)
    - A path to a transcript file

    Examples:
        qxub replay --list              # List available recordings
        qxub replay 1                   # Replay most recent session
        qxub replay --view 1            # View transcript of most recent
        qxub replay --speed 2 1         # Replay at 2x speed
        qxub replay 20260217-143000     # Replay specific session
    """
    # Resolve recording directory
    search_dir = record_dir or _get_default_record_dir()

    if list_recordings:
        # List all available recordings
        recordings = _find_recordings(search_dir)

        if not recordings:
            if search_dir:
                click.echo(f"No recordings found in {search_dir}")
            else:
                click.echo("No recordings found.")
                click.echo("")
                click.echo("To record sessions, use: qxi --record ...")
                click.echo("Or set defaults.interactive.record_dir in config")
            return

        click.echo(f"📼 Session recordings in {search_dir}:")
        click.echo("")

        for i, rec in enumerate(recordings, 1):
            ts = rec["timestamp"]
            ts_str = ts.strftime("%Y-%m-%d %H:%M") if ts else "unknown"
            duration = _format_duration(rec["duration"])
            has_timing = "✓" if rec["timing"] else "✗"

            click.echo(f"  {i:3d}  {ts_str}  ({duration})  timing:{has_timing}")

        click.echo("")
        click.echo("Use 'qxub replay N' to replay session N")
        click.echo("Use 'qxub replay --view N' to view transcript")
        return

    if session is None:
        # No session specified - show help
        click.echo("Usage: qxub replay [OPTIONS] SESSION")
        click.echo("")
        click.echo("Use 'qxub replay --list' to see available recordings")
        click.echo("Use 'qxub replay --help' for full options")
        return

    # Resolve session to a recording
    recording = None
    recordings = _find_recordings(search_dir)

    # Try parsing as a number
    try:
        idx = int(session)
        if 1 <= idx <= len(recordings):
            recording = recordings[idx - 1]
        else:
            click.echo(f"❌ Session {idx} not found. Use --list to see available.")
            sys.exit(1)
    except ValueError:
        pass

    # Try parsing as a timestamp
    if recording is None:
        for rec in recordings:
            if rec["timestamp"]:
                ts_str = rec["timestamp"].strftime("%Y%m%d-%H%M%S")
                if session in ts_str:
                    recording = rec
                    break

    # Try as a file path
    if recording is None:
        path = Path(session)
        if path.exists():
            recording = {
                "transcript": path,
                "timing": path.parent / path.name.replace("transcript", "timing"),
                "timestamp": None,
                "duration": None,
            }
            if not recording["timing"].exists():
                recording["timing"] = None

    if recording is None:
        click.echo(f"❌ Could not find session: {session}")
        click.echo("Use 'qxub replay --list' to see available recordings")
        sys.exit(1)

    transcript = recording["transcript"]
    timing = recording["timing"]

    if view:
        # View mode - just cat the transcript
        click.echo(f"📄 Transcript: {transcript}")
        click.echo("=" * 60)
        try:
            with open(transcript) as f:
                click.echo(f.read())
        except Exception as e:
            click.echo(f"❌ Error reading transcript: {e}")
            sys.exit(1)
        return

    # Replay mode
    if timing is None or not timing.exists():
        click.echo(f"⚠️  No timing file found for this session")
        click.echo(f"   Cannot replay without timing data")
        click.echo(f"   Use --view to see the transcript instead")
        sys.exit(1)

    click.echo(f"▶️  Replaying session: {transcript.name}")
    if recording["timestamp"]:
        click.echo(
            f"   Recorded: {recording['timestamp'].strftime('%Y-%m-%d %H:%M:%S')}"
        )
    click.echo(f"   Speed: {speed}x, max delay: {max_delay}s")
    click.echo("   Press Ctrl+C to stop")
    click.echo("")

    # Build scriptreplay command
    # scriptreplay [-m maxdelay] [-d divisor] timingfile [typescript]
    cmd = ["scriptreplay"]

    if max_delay > 0:
        cmd.extend(["-m", str(max_delay)])

    if speed != 1.0:
        # -d is divisor, so speed 2 means divisor 0.5
        divisor = 1.0 / speed
        cmd.extend(["-d", str(divisor)])

    cmd.extend([str(timing), str(transcript)])

    try:
        subprocess.run(cmd)
    except KeyboardInterrupt:
        click.echo("")
        click.echo("⏹️  Replay stopped")
    except FileNotFoundError:
        click.echo("❌ scriptreplay command not found")
        click.echo("   This is usually part of the 'util-linux' package")
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error during replay: {e}")
        sys.exit(1)
