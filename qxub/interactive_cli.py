"""
Interactive session subcommand for qxub.

This module provides the 'interactive' subcommand that creates interactive
PBS sessions with environment activation, resource allocation, and optional
tmux integration for session persistence.
"""

import os
import sys
from pathlib import Path

import click


def _get_config_manager():
    """Get the current config manager instance."""
    from .config import manager

    return manager.config_manager


def _get_default_shell() -> str:
    """Get the user's preferred shell."""
    return os.environ.get("SHELL", "/bin/bash")


def _get_default_runtime() -> str:
    """Get default runtime from config or fallback."""
    try:
        config_mgr = _get_config_manager()
        interactive_defaults = config_mgr.get("defaults.interactive", {})
        return interactive_defaults.get("runtime", "4:00:00")
    except Exception:
        return "4:00:00"


def _get_default_queue() -> str:
    """Get default queue from config or fallback."""
    try:
        config_mgr = _get_config_manager()
        interactive_defaults = config_mgr.get("defaults.interactive", {})
        return interactive_defaults.get("queue", "normal")
    except Exception:
        return "normal"


def _encode_command(cmd: str) -> str:
    """Base64 encode a command for safe shell passing."""
    return base64.b64encode(cmd.encode("utf-8")).decode("ascii")


def _build_interactive_script(
    env: str,
    shell: str,
    working_dir: str,
    pre_cmd: str | None,
    post_cmd: str | None,
    tmux_session: str | None,
) -> str:
    """
    Build the shell script that runs inside the interactive PBS job.

    This script:
    1. Activates the conda environment
    2. Runs pre-commands
    3. Optionally starts tmux
    4. Drops into interactive shell
    5. Runs post-commands on exit (best effort)
    """
    lines = [
        "#!/bin/bash",
        "",
        "# qxub interactive session",
        f'echo "🖥️  Interactive session starting on $(hostname)"',
        f'echo "📁 Working directory: {working_dir}"',
        f'echo "🐍 Conda environment: {env}"',
        "",
    ]

    # Change to working directory
    lines.extend(
        [
            f'cd "{working_dir}" || {{ echo "❌ Cannot cd to {working_dir}"; exit 1; }}',
            "",
        ]
    )

    # Conda activation
    lines.extend(
        [
            "# Activate conda environment",
            "if ! command -v conda > /dev/null 2>&1; then",
            '    echo "❌ ERROR: conda command not found"',
            "    exit 1",
            "fi",
            "",
            'eval "$(conda shell.bash hook)"',
            f"if ! conda activate {env}; then",
            f'    echo "❌ ERROR: Failed to activate conda environment: {env}"',
            "    exit 1",
            "fi",
            f'echo "✅ Activated conda environment: {env}"',
            "",
        ]
    )

    # Pre-command
    if pre_cmd:
        lines.extend(
            [
                "# Run pre-command",
                f'echo "⚡ Running pre-command..."',
                f"{pre_cmd}",
                "PRE_EXIT=$?",
                "if [ $PRE_EXIT -ne 0 ]; then",
                '    echo "❌ Pre-command failed with exit code $PRE_EXIT"',
                "    exit $PRE_EXIT",
                "fi",
                'echo "✅ Pre-command completed"',
                "",
            ]
        )

    # Post-command trap (best effort)
    if post_cmd:
        lines.extend(
            [
                "# Set up post-command trap",
                "_qxub_post_cmd() {",
                f'    echo "🏁 Running post-command..."',
                f"    {post_cmd}",
                '    echo "✅ Post-command completed"',
                "}",
                "",
                "# Trap EXIT to run post-command (best effort - may not run if killed)",
                "trap _qxub_post_cmd EXIT",
                "",
                "# Provide shortcut function to cleanly exit",
                "qxub_exit() {",
                '    echo "👋 Exiting interactive session..."',
                "    exit 0",
                "}",
                "echo \"💡 Tip: Run 'qxub_exit' to cleanly exit and run post-command\"",
                "",
            ]
        )

    # Main interactive part
    if tmux_session:
        # Tmux mode: create or attach to session
        lines.extend(
            [
                "# Tmux session management",
                f'TMUX_SESSION="{tmux_session}"',
                "",
                'if tmux has-session -t "$TMUX_SESSION" 2>/dev/null; then',
                '    echo "🔗 Attaching to existing tmux session: $TMUX_SESSION"',
                '    exec tmux attach-session -t "$TMUX_SESSION"',
                "else",
                '    echo "🆕 Creating new tmux session: $TMUX_SESSION"',
                f'    exec tmux new-session -s "$TMUX_SESSION" "{shell}"',
                "fi",
            ]
        )
    else:
        # Direct shell mode
        lines.extend(
            [
                "# Start interactive shell",
                f'echo "🚀 Starting interactive shell ({shell})"',
                'echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"',
                'echo ""',
                "",
                f"exec {shell}",
            ]
        )

    return "\n".join(lines)


def _build_qsub_command(
    resources: list[str],
    queue: str,
    project: str | None,
    name: str,
    working_dir: str,
    storage: str | None,
) -> list[str]:
    """Build the qsub -I command with all options."""
    cmd = ["qsub", "-I"]

    # Job name
    cmd.extend(["-N", name])

    # Queue
    cmd.extend(["-q", queue])

    # Project
    if project:
        cmd.extend(["-P", project])
    else:
        # Try to get from environment
        env_project = os.environ.get("PROJECT")
        if env_project:
            cmd.extend(["-P", env_project])

    # Resources
    for res in resources:
        cmd.extend(["-l", res])

    # Storage
    if storage:
        cmd.extend(["-l", f"storage={storage}"])

    return cmd


@click.command(name="interactive")
@click.option(
    "--env",
    required=True,  # For now, require conda env (will expand later)
    help="Conda environment name for the interactive session",
)
@click.option(
    "-l",
    "--resources",
    multiple=True,
    help="PBS resource specification (e.g., 'walltime=4:00:00,mem=16GB')",
)
@click.option(
    "-q",
    "--queue",
    default=None,
    help=f"PBS queue name (default: from config or 'normal')",
)
@click.option(
    "-P",
    "--project",
    default=None,
    help="PBS project code (default: $PROJECT)",
)
@click.option(
    "-N",
    "--name",
    default=None,
    help="PBS job name (default: qxub-interactive)",
)
@click.option(
    "--mem",
    "--memory",
    help="Memory requirement (e.g., '16GB', '4000MB')",
)
@click.option(
    "--runtime",
    "--time",
    default=None,
    help="Session duration (default: from config or '4:00:00')",
)
@click.option(
    "--cpus",
    "--threads",
    type=int,
    help="Number of CPU cores",
)
@click.option(
    "--disk",
    "--jobfs",
    help="Local disk/jobfs requirement",
)
@click.option(
    "--volumes",
    "--storage",
    help="Storage volumes (e.g., 'gdata/a56+scratch/a56')",
)
@click.option(
    "--shell",
    default=None,
    help="Shell to use (default: $SHELL or /bin/bash)",
)
@click.option(
    "--working-dir",
    "--execdir",
    default=None,
    help="Starting directory (default: current directory)",
)
@click.option(
    "--pre",
    help="Commands to run before starting the shell",
)
@click.option(
    "--post",
    help="Commands to run on exit (best effort via trap)",
)
@click.option(
    "--tmux",
    "tmux_session",
    default=None,
    help="Start/attach to a tmux session with this name",
)
@click.option(
    "--dry",
    "--dry-run",
    is_flag=True,
    help="Show what would be executed without running",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity",
)
@click.pass_context
def interactive_cli(
    ctx,
    env,
    resources,
    queue,
    project,
    name,
    mem,
    runtime,
    cpus,
    disk,
    volumes,
    shell,
    working_dir,
    pre,
    post,
    tmux_session,
    dry,
    verbose,
):
    """
    Start an interactive PBS session with environment setup.

    This command creates an interactive PBS job (`qsub -I`) with your
    conda environment pre-activated, making it easy to work interactively
    on compute nodes.

    \b
    Basic usage:
        qxub interactive --env myenv

    \b
    With resources:
        qxub interactive --env myenv --mem 32GB --cpus 4 --runtime 2h

    \b
    With tmux for session persistence:
        qxub interactive --env myenv --tmux mysession
        # If disconnected, resubmit to reattach

    \b
    With pre/post commands:
        qxub interactive --env myenv --pre "cd /project && git pull"
        qxub interactive --env myenv --post "echo 'Session ended at $(date)'"

    \b
    Tips:
        - Use Ctrl+D or 'exit' to end the session normally
        - If --post is set, use 'qxub_exit' function to ensure post runs
        - Use --tmux for long sessions to enable reconnection
    """
    from .resources import ResourceMapper

    # Resolve defaults
    shell = shell or _get_default_shell()
    working_dir = working_dir or os.getcwd()
    queue = queue or _get_default_queue()
    runtime = runtime or _get_default_runtime()
    name = name or "qxub-interactive"

    # Build resources using ResourceMapper
    mapper = ResourceMapper()

    # Add walltime
    try:
        mapper.add_runtime(runtime)
    except Exception:
        # Fall back to raw format if parsing fails
        mapper.pbs_resources.append(f"walltime={runtime}")

    # Add memory if specified
    if mem:
        try:
            mapper.add_memory(mem)
        except Exception:
            mapper.pbs_resources.append(f"mem={mem}")

    # Add CPUs if specified
    if cpus:
        mapper.add_cpus(cpus)

    # Add jobfs if specified
    if disk:
        try:
            mapper.add_disk(disk)
        except Exception:
            mapper.pbs_resources.append(f"jobfs={disk}")

    # Get the combined resource string from mapper
    mapped_resources = mapper.get_pbs_resources()

    # Combine with any -l resources passed directly
    all_resources = list(resources) + mapped_resources

    # Build the interactive script
    script = _build_interactive_script(
        env=env,
        shell=shell,
        working_dir=working_dir,
        pre_cmd=pre,
        post_cmd=post,
        tmux_session=tmux_session,
    )

    # Build qsub command
    qsub_cmd = _build_qsub_command(
        resources=all_resources,
        queue=queue,
        project=project,
        name=name,
        working_dir=working_dir,
        storage=volumes,
    )

    if verbose or dry:
        click.echo("=" * 60)
        click.echo("qxub interactive session")
        click.echo("=" * 60)
        click.echo(f"Environment: {env}")
        click.echo(f"Shell: {shell}")
        click.echo(f"Working dir: {working_dir}")
        click.echo(f"Queue: {queue}")
        click.echo(f"Runtime: {runtime}")
        if mem:
            click.echo(f"Memory: {mem}")
        if cpus:
            click.echo(f"CPUs: {cpus}")
        if tmux_session:
            click.echo(f"Tmux session: {tmux_session}")
        if pre:
            click.echo(f"Pre-command: {pre}")
        if post:
            click.echo(f"Post-command: {post}")
        click.echo("")
        click.echo("qsub command:")
        click.echo(f"  {' '.join(qsub_cmd)}")
        click.echo("")

    if dry:
        click.echo("Interactive script that would run:")
        click.echo("-" * 40)
        click.echo(script)
        click.echo("-" * 40)
        return

    # Log to history using the existing API
    try:
        from .history import history_manager

        # Log execution with context - the history manager will extract
        # recipe info from sys.argv
        history_manager.log_execution(
            ctx=ctx,
            success=True,  # We assume success on start
            job_id=None,  # We don't have job_id yet (qsub -I hasn't run)
        )
    except Exception:
        pass  # Don't fail if history logging fails

    # Execute interactive session
    click.echo(f"🚀 Starting interactive session...")
    click.echo(f"   Environment: {env}")
    click.echo(f"   Queue: {queue}")
    click.echo(f"   Runtime: {runtime}")
    if tmux_session:
        click.echo(f"   Tmux: {tmux_session}")
    click.echo("")

    # Write the initialization script to a temporary file
    # This will be sourced when the interactive job starts
    import tempfile

    script_dir = Path(tempfile.gettempdir()) / "qxub_interactive"
    script_dir.mkdir(parents=True, exist_ok=True)

    # Use a unique filename based on timestamp
    import time

    script_file = script_dir / f"init_{int(time.time())}_{os.getpid()}.sh"

    with open(script_file, "w") as f:
        f.write(script)
    script_file.chmod(0o755)

    # Build the final qsub command that runs our script
    # qsub -I executes interactively, and we pass our script as the command
    qsub_cmd.extend(["--", str(script_file)])

    if verbose:
        click.echo(f"Script file: {script_file}")
        click.echo(f"Full command: {' '.join(qsub_cmd)}")
        click.echo("")

    # Use os.execvp to replace this process with qsub -I
    # This ensures proper TTY handling and signal propagation
    try:
        os.execvp("qsub", qsub_cmd)
    except Exception as e:
        # Clean up script file on error
        try:
            script_file.unlink()
        except Exception:
            pass
        click.echo(f"❌ Error starting interactive session: {e}")
        sys.exit(1)


def qxi_main():
    """Entry point for the qxi command (alias for qxub interactive)."""
    interactive_cli(standalone_mode=True)
