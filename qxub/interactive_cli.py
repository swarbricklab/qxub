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
        # Check interactive defaults first
        interactive_runtime = config_mgr.get_config_value(
            "defaults.interactive.runtime"
        )
        if interactive_runtime:
            return interactive_runtime
        # Fall back to general defaults
        general_runtime = config_mgr.get_config_value("defaults.runtime")
        if general_runtime:
            return general_runtime
        return "4:00:00"
    except Exception:
        return "4:00:00"


def _get_default_queue() -> str:
    """Get default queue from config or fallback."""
    try:
        config_mgr = _get_config_manager()
        # Check interactive defaults first
        interactive_queue = config_mgr.get_config_value("defaults.interactive.queue")
        if interactive_queue:
            return interactive_queue
        # Fall back to general defaults
        general_queue = config_mgr.get_config_value("defaults.queue")
        if general_queue:
            return general_queue
        return "normal"
    except Exception:
        return "normal"


def _get_default_volumes() -> str | None:
    """Get default storage volumes from config."""
    try:
        config_mgr = _get_config_manager()
        # Check interactive defaults first
        interactive_volumes = config_mgr.get_config_value(
            "defaults.interactive.volumes"
        )
        if interactive_volumes:
            return interactive_volumes
        # Fall back to general defaults
        return config_mgr.get_config_value("defaults.volumes")
    except Exception:
        return None


def _get_default_modules() -> list[str]:
    """Get default modules from config."""
    try:
        config_mgr = _get_config_manager()
        modules = config_mgr.get_config_value("defaults.interactive.modules")
        if modules:
            # Handle both string and list formats
            if isinstance(modules, str):
                import re

                return [m.strip() for m in re.split(r"[,\s]+", modules) if m.strip()]
            elif isinstance(modules, list):
                return [str(m).strip() for m in modules if m]
        return []
    except Exception:
        return []


def _get_default_env() -> str | None:
    """Get default conda environment from config."""
    try:
        config_mgr = _get_config_manager()
        return config_mgr.get_config_value("defaults.interactive.env")
    except Exception:
        return None


def _get_default_tmux_session() -> str | None:
    """Get default tmux session name from config."""
    try:
        config_mgr = _get_config_manager()
        return config_mgr.get_config_value("defaults.interactive.tmux_session")
    except Exception:
        return None


def _get_default_record_dir() -> str | None:
    """Get default transcript recording directory from config.

    Supports template variables like {user}, {project}, etc.
    """
    try:
        config_mgr = _get_config_manager()
        record_dir = config_mgr.get_config_value("defaults.interactive.record_dir")
        if record_dir:
            # Resolve template variables
            template_vars = {
                "user": os.environ.get("USER", "unknown"),
                "project": os.environ.get("PROJECT", "unknown"),
            }
            return config_mgr.resolve_templates(record_dir, template_vars)
        return None
    except Exception:
        return None


def _get_container_prompt(container_name: str) -> str:
    """Get container prompt from config or use default.

    The prompt can include {container} as a placeholder for the container name.
    Default: "\\[\\033[1;35m\\]📦 {container}\\[\\033[0m\\] \\[\\033[1;34m\\]\\W\\[\\033[0m\\] \\$ "
    """
    default_prompt = "\\[\\033[1;35m\\]📦 {container}\\[\\033[0m\\] \\[\\033[1;34m\\]\\W\\[\\033[0m\\] \\$ "
    try:
        config_mgr = _get_config_manager()
        custom_prompt = config_mgr.get_config_value(
            "defaults.interactive.container_prompt"
        )
        if custom_prompt:
            return custom_prompt.replace("{container}", container_name)
        return default_prompt.replace("{container}", container_name)
    except Exception:
        return default_prompt.replace("{container}", container_name)


def _encode_command(cmd: str) -> str:
    """Base64 encode a command for safe shell passing."""
    return base64.b64encode(cmd.encode("utf-8")).decode("ascii")


def _build_interactive_script(
    shell: str,
    working_dir: str,
    pre_cmd: str | None,
    post_cmd: str | None,
    conda_env: str | None = None,
    modules: list[str] | None = None,
    container: str | None = None,
    bind: str | None = None,
    verbose: bool = False,
    run_cmd: str | None = None,
    record: bool = False,
    record_dir: str | None = None,
) -> str:
    """
    Build the shell script that runs inside the interactive PBS job.

    This script:
    1. Loads modules (if any) - with environment diff tracking
    2. Activates conda environment (if any)
    3. Runs pre-commands
    4. Drops into interactive shell
    5. Runs post-commands on exit (best effort)
    """
    # Determine context type for display
    if conda_env and modules:
        context_type = "conda+modules"
        context_display = f"🐍 Conda: {conda_env} + 📦 Modules: {', '.join(modules)}"
    elif conda_env:
        context_type = "conda"
        context_display = f"🐍 Conda environment: {conda_env}"
    elif modules:
        context_type = "module"
        context_display = f"📦 Modules: {', '.join(modules)}"
    elif container:
        context_type = "container"
        context_display = f"📦 Container: {container}"
    else:
        context_type = "default"
        context_display = "🔧 Default environment (no conda/modules/container)"

    lines = [
        "#!/bin/bash",
        "",
        "# qxub interactive session",
        f'export QXUB_VERBOSE={"1" if verbose else "0"}',
        f'echo "🖥️  Interactive session starting on $(hostname)"',
        f'echo "📁 Working directory: {working_dir}"',
        f'echo "{context_display}"',
        "",
    ]

    # Change to working directory
    lines.extend(
        [
            f'cd "{working_dir}" || {{ echo "❌ Cannot cd to {working_dir}"; exit 1; }}',
            "",
        ]
    )

    # Source bash initialization
    # PBS jobs don't source .bashrc by default for non-login shells
    lines.extend(
        [
            "# Initialize shell environment (PBS doesn't source .bashrc by default)",
            "if [ -f ~/.bashrc ]; then",
            '    echo "🔧 Sourcing ~/.bashrc..."',
            "    source ~/.bashrc",
            "fi",
            "",
        ]
    )

    # Environment-specific activation
    # Order: modules first (if any), then conda (if any)
    # This allows utility modules (singularity, dvc, etc.) to coexist with conda envs

    if modules:
        # Module loading with environment diff tracking
        module_list_str = " ".join(modules)
        lines.extend(
            [
                "# Load environment modules",
                "if ! command -v module > /dev/null 2>&1; then",
                '    echo "⚠️  module command not found, trying to initialize..."',
                "    # Try common module initialization paths",
                "    if [ -f /etc/profile.d/modules.sh ]; then",
                "        source /etc/profile.d/modules.sh",
                "    elif [ -f /opt/modules/init/bash ]; then",
                "        source /opt/modules/init/bash",
                "    fi",
                "fi",
                "",
                "if ! command -v module > /dev/null 2>&1; then",
                '    echo "❌ ERROR: module command not found"',
                "    exit 1",
                "fi",
                "",
                "# Capture environment before loading modules (for diff tracking)",
                "QXUB_ENV_BEFORE=$(mktemp)",
                'printenv | sort > "$QXUB_ENV_BEFORE"',
                "",
            ]
        )
        # Load each module
        for mod in modules:
            lines.extend(
                [
                    f'echo "📦 Loading module: {mod}"',
                    f"if ! module load {mod}; then",
                    f'    echo "❌ ERROR: Failed to load module: {mod}"',
                    '    rm -f "$QXUB_ENV_BEFORE"',
                    "    exit 1",
                    "fi",
                ]
            )
        lines.extend(
            [
                f'echo "✅ Loaded modules: {module_list_str}"',
                "",
                "# Check for environment changes that might conflict with conda",
                "QXUB_ENV_AFTER=$(mktemp)",
                'printenv | sort > "$QXUB_ENV_AFTER"',
                "",
                "# Warn about potentially problematic environment changes",
                "for VAR in PYTHONPATH PYTHONHOME LD_LIBRARY_PATH R_HOME R_LIBS; do",
                '    if diff "$QXUB_ENV_BEFORE" "$QXUB_ENV_AFTER" | grep -q "^[<>].*$VAR="; then',
                '        echo "⚠️  Warning: Module modified $VAR - may conflict with conda"',
                "    fi",
                "done",
                "",
                "# Show full environment diff in verbose mode",
                'if [ "${QXUB_VERBOSE:-0}" = "1" ]; then',
                '    echo ""',
                '    echo "📊 Environment changes from modules:"',
                '    diff "$QXUB_ENV_BEFORE" "$QXUB_ENV_AFTER" | grep -E "^[<>]" || echo "  (no changes)"',
                '    echo ""',
                "fi",
                "",
                'rm -f "$QXUB_ENV_BEFORE" "$QXUB_ENV_AFTER"',
                "",
            ]
        )

    if conda_env:
        # Conda activation (after modules if both specified)
        lines.extend(
            [
                "# Activate conda environment",
                "if ! command -v conda > /dev/null 2>&1; then",
                '    echo "⚠️  conda not found in PATH, trying common locations..."',
                "    # Try common conda installation paths",
                "    for CONDA_PATH in \\",
                '        "$HOME/miniconda3" \\',
                '        "$HOME/anaconda3" \\',
                '        "$HOME/miniforge3" \\',
                '        "$HOME/mambaforge" \\',
                '        "/g/data/a56/conda/miniconda3" \\',
                '        "/apps/conda"; do',
                '        if [ -f "$CONDA_PATH/etc/profile.d/conda.sh" ]; then',
                '            echo "📦 Found conda at: $CONDA_PATH"',
                '            source "$CONDA_PATH/etc/profile.d/conda.sh"',
                "            break",
                "        fi",
                "    done",
                "fi",
                "",
                "# Final check for conda",
                "if ! command -v conda > /dev/null 2>&1; then",
                '    echo "❌ ERROR: conda command not found"',
                '    echo "💡 Tip: Ensure conda is initialized in ~/.bashrc"',
                "    exit 1",
                "fi",
                "",
                'eval "$(conda shell.bash hook)"',
                f"if ! conda activate {conda_env}; then",
                f'    echo "❌ ERROR: Failed to activate conda environment: {conda_env}"',
                "    exit 1",
                "fi",
                f'echo "✅ Activated conda environment: {conda_env}"',
                "",
            ]
        )

    if container:
        # Singularity container - load singularity module first
        # Note: container is mutually exclusive with conda/modules
        lines.extend(
            [
                "# Load singularity module",
                "if ! command -v singularity > /dev/null 2>&1; then",
                '    echo "📦 Loading singularity module..."',
                "    if ! module load singularity 2>/dev/null; then",
                '        echo "❌ ERROR: Failed to load singularity module"',
                "        exit 1",
                "    fi",
                "fi",
                "",
                f'echo "📦 Using container: {container}"',
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
    if container:
        # Container mode: exec into singularity shell directly
        # Extract container basename for prompt
        container_name = container.split("/")[-1].replace(".sif", "")
        bind_opts = f"--bind {bind} " if bind else ""

        # Get prompt from config or use default
        container_prompt = _get_container_prompt(container_name)

        # Set a custom prompt that shows we're in a container
        # SINGULARITYENV_* variables are passed through to the container
        lines.extend(
            [
                f'echo "🚀 Entering container: {container}"',
                'echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"',
                'echo ""',
                "",
                "# Set custom prompt for container shell",
                f'export SINGULARITYENV_PS1="{container_prompt}"',
                "",
                f'exec singularity shell {bind_opts}"{container}"',
            ]
        )
    else:
        # Direct shell mode - create a custom rcfile that re-activates environment
        # This is necessary because exec replaces the process, losing activation
        lines.extend(
            [
                "# Start interactive shell with environment",
                f'echo "🚀 Starting interactive shell ({shell})"',
                'echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"',
                'echo ""',
                "",
                "# Create a temporary rcfile that sources bashrc and activates environment",
                "QXUB_RCFILE=$(mktemp /tmp/qxub_rc.XXXXXX)",
                "cat > \"$QXUB_RCFILE\" << 'RCEOF'",
                "# Source standard bashrc first",
                "[ -f ~/.bashrc ] && source ~/.bashrc",
            ]
        )

        # Add environment-specific activation to rcfile
        # Order: modules first, then conda (same as main script)
        if modules:
            for mod in modules:
                lines.append(f"module load {mod} 2>/dev/null")

        if conda_env:
            lines.extend(
                [
                    "# Ensure conda is available and activate the environment",
                    'eval "$(conda shell.bash hook)" 2>/dev/null',
                    f"conda activate {conda_env}",
                ]
            )

        lines.extend(
            [
                "# Clean up the temp rcfile",
                'rm -f "$QXUB_RCFILE" 2>/dev/null',
                "RCEOF",
                "",
            ]
        )

        if run_cmd:
            # Run the specified command instead of dropping to shell
            #
            # Signal handling challenge in PBS:
            # PBS runs scripts in a non-interactive context. When Ctrl+C is pressed,
            # SIGINT goes to the process group, killing everything including R.
            #
            # Solution: Use the `script` command to create a pseudo-terminal (pty).
            # Programs like R detect they're running on a pty and set up proper
            # signal handlers that interrupt the current command instead of exiting.
            #
            # The `script` command:
            # - Creates a new pseudo-terminal
            # - Runs the command with that pty as stdin/stdout
            # - R sees a proper terminal and handles Ctrl+C correctly
            # - `-q`: quiet mode (no "Script started" messages)
            # - `-c`: run this command
            # - `-e`: return exit code of the child process
            # - With --record: save transcript and timing for playback
            lines.extend(
                [
                    f'echo "🚀 Running: {run_cmd}"',
                    'echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"',
                    'echo ""',
                    "",
                    "# Source the rcfile to set up environment (conda/modules)",
                    'source "$QXUB_RCFILE"',
                    "",
                ]
            )

            if record:
                # Record transcript and timing data
                # Use specified directory or working directory
                rec_dir = record_dir if record_dir else "$PWD"
                lines.extend(
                    [
                        "# Generate transcript filename with timestamp",
                        f'RECORD_DIR="{rec_dir}"',
                        'mkdir -p "$RECORD_DIR" 2>/dev/null',
                        'TRANSCRIPT_FILE="$RECORD_DIR/qxub-transcript-$(date +%Y%m%d-%H%M%S).txt"',
                        'TIMING_FILE="$RECORD_DIR/qxub-timing-$(date +%Y%m%d-%H%M%S).txt"',
                        'echo "📝 Recording session to: $TRANSCRIPT_FILE"',
                        'echo "⏱️  Timing data: $TIMING_FILE"',
                        'echo "   Replay with: scriptreplay "$TIMING_FILE" "$TRANSCRIPT_FILE""',
                        'echo ""',
                        "",
                        "# Use 'script' with recording enabled",
                        "# -e: return child's exit code",
                        "# -q: quiet mode (no 'Script started' messages)",
                        "# -t: output timing to file for replay",
                        f'exec script -eqc "{run_cmd}" -t"$TIMING_FILE" "$TRANSCRIPT_FILE"',
                    ]
                )
            else:
                lines.extend(
                    [
                        "# Use 'script' to create a pseudo-terminal for the command",
                        "# -e: return child's exit code, -q: quiet mode",
                        f'exec script -eqc "{run_cmd}" /dev/null',
                    ]
                )
        else:
            # Start interactive shell
            lines.extend(
                [
                    "# Start bash with our custom rcfile",
                    'exec /bin/bash --rcfile "$QXUB_RCFILE"',
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
    default=None,
    help="Conda environment name for the interactive session",
)
@click.option(
    "--mod",
    "--module",
    "mod",
    multiple=True,
    help="Module to load (can be used multiple times)",
)
@click.option(
    "--mods",
    "--modules",
    "mods",
    default=None,
    help="Comma or space-separated list of modules to load",
)
@click.option(
    "--sif",
    "--container",
    default=None,
    help="Singularity container (.sif file) for the interactive session",
)
@click.option(
    "--bind",
    default=None,
    help="Singularity bind mounts (e.g., '/data:/mnt', only used with --sif)",
)
@click.option(
    "--shortcut",
    default=None,
    help="Use a predefined shortcut for execution settings",
)
@click.option(
    "--alias",
    default=None,
    help="Use a predefined alias for execution settings",
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
@click.option(
    "--cmd",
    default=None,
    help="Command to run (alternative to positional arguments)",
)
@click.option(
    "--no-defaults",
    is_flag=True,
    help="Skip default modules and conda env from config",
)
@click.option(
    "--record",
    is_flag=True,
    help="Record session transcript",
)
@click.option(
    "--record-dir",
    default=None,
    help="Directory for transcript files (default: config or working directory)",
)
@click.argument("command", nargs=-1, required=False)
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
    cmd,
    no_defaults,
    record,
    record_dir,
    mod,
    mods,
    sif,
    bind,
    shortcut,
    alias,
    command,
):
    """
    Start an interactive PBS session with environment setup.

    This command creates an interactive PBS job (`qsub -I`) with your
    environment pre-activated, making it easy to work interactively
    on compute nodes.

    \b
    Basic usage (drops into shell):
        qxub interactive --env myenv
        qxi --env myenv

    \b
    Run a command directly (e.g., start Python/R interpreter):
        qxi --env myenv -- python
        qxi --env myenv -- ipython
        qxi --env myenv -- R
        qxi --cmd "python -i script.py"  # alternative for complex args

    \b
    With environment modules:
        qxub interactive --mod python3 --mod numpy
        qxub interactive --mods "python3 numpy scipy"

    \b
    With conda + modules (modules loaded first, then conda):
        qxub interactive --env myenv --mod singularity --mod dvc
        # Warnings shown if modules modify PYTHONPATH, LD_LIBRARY_PATH, etc.
        # Use -v to see full environment diff from module loading

    \b
    With Singularity container:
        qxub interactive --sif /path/to/container.sif
        qxub interactive --sif container.sif --bind /data:/mnt

    \b
    With resources:
        qxub interactive --env myenv --mem 32GB --cpus 4 --runtime 2h

    \b
    With tmux for session persistence (on login node):
        qxub interactive --env myenv --tmux mysession
        # Tmux runs on login node - survives PBS job ending
        # If disconnected: tmux attach -t mysession
        # To resubmit PBS job from inside tmux, run qxi again

    \b
    With pre/post commands:
        qxub interactive --env myenv --pre "cd /project && git pull"
        qxub interactive --mods python3 --post "echo 'Session ended'"

    \b
    With shortcuts (auto-detected or explicit):
        qxi -- python                  # Auto-detects 'python' shortcut
        qxi --shortcut pytorch         # Explicit shortcut usage
        qxi --shortcut gpu -- python   # Shortcut + command

    \b
    With aliases:
        qxi --alias train              # Use predefined alias settings

    \b
    Tips:
        - Use Ctrl+D or 'exit' to end the PBS session normally
        - If --post is set, use 'qxub_exit' function to ensure post runs
        - With --tmux: session persists on login node after PBS job ends
        - Reconnect to tmux with: tmux attach -t <session>
    """
    from .resources import ResourceMapper

    # Check for conflicting shortcut/alias usage
    if shortcut and alias:
        raise click.ClickException(
            "Cannot specify both --shortcut and --alias. Use only one."
        )

    # Handle shortcut processing
    if shortcut:
        from .shortcut_manager import ShortcutManager

        shortcut_manager = ShortcutManager()
        shortcut_def = shortcut_manager.get_shortcut(shortcut)

        if not shortcut_def:
            available_shortcuts = shortcut_manager.list_shortcuts()
            click.echo(f"❌ Shortcut '{shortcut}' not found")
            if available_shortcuts:
                click.echo("💡 Available shortcuts:")
                for name in sorted(available_shortcuts.keys()):
                    click.echo(f"  • {name}")
            raise click.ClickException(f"Shortcut '{shortcut}' not found")

        # Apply shortcut settings (CLI options override shortcut settings)
        if not env and shortcut_def.get("env"):
            env = shortcut_def["env"]
        if not mod and shortcut_def.get("mod"):
            mod_val = shortcut_def["mod"]
            mod = tuple(mod_val) if isinstance(mod_val, list) else (mod_val,)
        if not mods and shortcut_def.get("mods"):
            mods = shortcut_def["mods"]
        if not queue and shortcut_def.get("queue"):
            queue = shortcut_def["queue"]
        if not resources and shortcut_def.get("resources"):
            res_val = shortcut_def["resources"]
            resources = tuple(res_val) if isinstance(res_val, list) else (res_val,)
        if not sif and shortcut_def.get("sif"):
            sif = shortcut_def["sif"]
        if not bind and shortcut_def.get("bind"):
            bind = shortcut_def["bind"]

        click.echo(f"🎯 Using shortcut '{shortcut}'")

    # Handle alias processing
    elif alias:
        config_mgr = _get_config_manager()
        alias_def = config_mgr.get_alias(alias)

        if not alias_def:
            available_aliases = config_mgr.list_aliases()
            click.echo(f"❌ Alias '{alias}' not found")
            if available_aliases:
                click.echo("💡 Available aliases:")
                for alias_name in sorted(available_aliases):
                    click.echo(f"  • {alias_name}")
            raise click.ClickException(f"Alias '{alias}' not found")

        # Flatten hierarchical alias structure if present
        flat_alias = {}
        if "main" in alias_def and isinstance(alias_def["main"], dict):
            flat_alias.update(alias_def["main"])
        if "subcommand" in alias_def and isinstance(alias_def["subcommand"], dict):
            subcommand_section = alias_def["subcommand"]
            subcommand_type = subcommand_section.get("type")
            if subcommand_type == "conda" and "env" in subcommand_section:
                flat_alias["env"] = subcommand_section["env"]
            elif subcommand_type == "module":
                if "mod" in subcommand_section:
                    flat_alias["mod"] = subcommand_section["mod"]
                elif "mods" in subcommand_section:
                    flat_alias["mods"] = subcommand_section["mods"]
            elif subcommand_type == "sing" and "sif" in subcommand_section:
                flat_alias["sif"] = subcommand_section["sif"]
                if "bind" in subcommand_section:
                    flat_alias["bind"] = subcommand_section["bind"]
        # Also handle flat structure
        for key, value in alias_def.items():
            if key not in ["main", "subcommand", "target"]:
                flat_alias[key] = value

        # Apply alias settings (CLI options override alias settings)
        if not env and flat_alias.get("env"):
            env = flat_alias["env"]
        if not mod and flat_alias.get("mod"):
            mod_val = flat_alias["mod"]
            mod = tuple(mod_val) if isinstance(mod_val, list) else (mod_val,)
        if not mods and flat_alias.get("mods"):
            mods = flat_alias["mods"]
        if not queue and flat_alias.get("queue"):
            queue = flat_alias["queue"]
        if not resources and flat_alias.get("resources"):
            res_val = flat_alias["resources"]
            resources = tuple(res_val) if isinstance(res_val, list) else (res_val,)
        if not sif and flat_alias.get("sif"):
            sif = flat_alias["sif"]
        if not bind and flat_alias.get("bind"):
            bind = flat_alias["bind"]

        click.echo(f"🎯 Using alias '{alias}'")

    # Auto-detect shortcut from command (if no explicit shortcut/alias and command provided)
    elif command and not env and not mod and not mods and not sif:
        from .shortcut_manager import ShortcutManager

        shortcut_manager = ShortcutManager()
        shortcut_match = shortcut_manager.find_shortcut(list(command))

        if shortcut_match:
            shortcut_def = shortcut_match["definition"]
            shortcut_name = shortcut_match["name"]

            # Apply shortcut settings
            if shortcut_def.get("env"):
                env = shortcut_def["env"]
            if shortcut_def.get("mod"):
                mod_val = shortcut_def["mod"]
                mod = tuple(mod_val) if isinstance(mod_val, list) else (mod_val,)
            if shortcut_def.get("mods"):
                mods = shortcut_def["mods"]
            if not queue and shortcut_def.get("queue"):
                queue = shortcut_def["queue"]
            if not resources and shortcut_def.get("resources"):
                res_val = shortcut_def["resources"]
                resources = tuple(res_val) if isinstance(res_val, list) else (res_val,)

            click.echo(f"🎯 Auto-detected shortcut '{shortcut_name}'")

    # Resolve command from --cmd or positional arguments
    if cmd and command:
        raise click.ClickException(
            "Cannot specify both --cmd and positional command arguments"
        )

    run_cmd = None
    if cmd:
        run_cmd = cmd
    elif command:
        run_cmd = " ".join(command)

    # Process module options - combine --mod, --mods, and config defaults
    # Skip defaults when using container mode (--sif)
    module_list = list(mod) if mod else []
    if mods:
        # Split on comma or space
        import re

        module_list.extend(re.split(r"[,\s]+", mods.strip()))
    module_list = [m.strip() for m in module_list if m.strip()]

    # Add default modules from config (only if not using container and not --no-defaults)
    if not sif and not no_defaults:
        default_modules = _get_default_modules()
        for default_mod in default_modules:
            if default_mod not in module_list:
                module_list.append(default_mod)

    # Apply default conda environment from config (if --env not specified and not container)
    if not env and not sif and not no_defaults:
        env = _get_default_env()

    # Validate execution contexts:
    # - Container (--sif) is mutually exclusive with conda/modules
    # - Conda (--env) and modules (--mod/--mods) CAN be combined
    if sif and (env or module_list):
        raise click.ClickException(
            "Cannot use --sif (container) with --env or --mod/--mods. "
            "Container mode is mutually exclusive with other execution contexts."
        )

    # Warn if --bind is used without --sif
    if bind and not sif:
        raise click.ClickException("--bind can only be used with --sif (container).")

    # Resolve defaults
    shell = shell or _get_default_shell()
    working_dir = working_dir or os.getcwd()
    queue = queue or _get_default_queue()
    runtime = runtime or _get_default_runtime()
    name = name or "qxub-interactive"

    # Handle tmux: detect if already inside a tmux session
    already_in_tmux = os.environ.get("TMUX") is not None

    if already_in_tmux and tmux_session:
        # User explicitly requested --tmux but we're already in tmux
        raise click.ClickException(
            "Already inside a tmux session. Cannot create nested tmux sessions.\n"
            "   Run without --tmux to start the PBS job directly in this session."
        )

    if already_in_tmux:
        # Already in tmux - skip default session name, run qsub directly
        tmux_session = None
    elif not tmux_session and not no_defaults:
        # Not in tmux - apply default session from config
        tmux_session = _get_default_tmux_session()

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
    # Resolve record_dir from config if not specified
    resolved_record_dir = record_dir or _get_default_record_dir()
    script = _build_interactive_script(
        shell=shell,
        working_dir=working_dir,
        pre_cmd=pre,
        post_cmd=post,
        conda_env=env,
        modules=module_list if module_list else None,
        container=sif,
        bind=bind,
        verbose=verbose > 0,
        run_cmd=run_cmd,
        record=record,
        record_dir=resolved_record_dir,
    )

    # Resolve storage volumes - use default if not specified
    storage = volumes or _get_default_volumes()

    # Build qsub command
    qsub_cmd = _build_qsub_command(
        resources=all_resources,
        queue=queue,
        project=project,
        name=name,
        working_dir=working_dir,
        storage=storage,
    )

    # Determine context description for output
    if env and module_list:
        context_desc = f"Conda: {env} + Modules: {', '.join(module_list)}"
    elif env:
        context_desc = f"Conda env: {env}"
    elif module_list:
        context_desc = f"Modules: {', '.join(module_list)}"
    elif sif:
        context_desc = f"Container: {sif}"
        if bind:
            context_desc += f" (bind: {bind})"
    else:
        context_desc = "Default (no conda/modules/container)"

    if verbose or dry:
        click.echo("=" * 60)
        click.echo("qxub interactive session")
        click.echo("=" * 60)
        click.echo(f"Context: {context_desc}")
        click.echo(f"Shell: {shell}")
        click.echo(f"Working dir: {working_dir}")
        click.echo(f"Queue: {queue}")
        click.echo(f"Runtime: {runtime}")
        if mem:
            click.echo(f"Memory: {mem}")
        if cpus:
            click.echo(f"CPUs: {cpus}")
        if storage:
            click.echo(f"Storage: {storage}")
        if tmux_session:
            click.echo(f"Tmux session: {tmux_session}")
        if pre:
            click.echo(f"Pre-command: {pre}")
        if post:
            click.echo(f"Post-command: {post}")
        if run_cmd:
            click.echo(f"Run command: {run_cmd}")
        if record:
            record_location = resolved_record_dir or working_dir
            click.echo(f"Recording to: {record_location}")
        if no_defaults:
            click.echo("Config defaults: SKIPPED (--no-defaults)")
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
    click.echo("🚀 Starting interactive session...")
    click.echo(f"   {context_desc}")
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

    # Handle tmux: start tmux on LOGIN node, run qsub INSIDE tmux
    # This ensures the tmux session persists when the PBS job ends
    if tmux_session:
        import shlex
        import subprocess

        qsub_cmd_str = " ".join(shlex.quote(arg) for arg in qsub_cmd)

        # Check if tmux session already exists
        result = subprocess.run(
            ["tmux", "has-session", "-t", tmux_session],
            capture_output=True,
        )
        session_exists = result.returncode == 0

        if session_exists:
            # Session exists - get its current directory and send info message
            pane_path_result = subprocess.run(
                [
                    "tmux",
                    "display-message",
                    "-t",
                    tmux_session,
                    "-p",
                    "#{pane_current_path}",
                ],
                capture_output=True,
                text=True,
            )
            session_pwd = (
                pane_path_result.stdout.strip()
                if pane_path_result.returncode == 0
                else "unknown"
            )

            # Build info message to send to the session
            info_lines = [
                "",
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
                f"🔗 Attached to existing tmux session: {tmux_session}",
            ]
            if session_pwd != working_dir:
                info_lines.append(f"📁 Session PWD: {session_pwd}")
                info_lines.append(f"   (You requested: {working_dir})")
                info_lines.append("   Use 'cd' to change directory if needed")
            info_lines.append("   If PBS job ended, run 'qxi' again to start a new job")
            info_lines.append(
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            )
            info_lines.append("")

            # Send the info message to the tmux session (will appear after attach)
            info_message = "\\n".join(info_lines)
            subprocess.run(
                [
                    "tmux",
                    "send-keys",
                    "-t",
                    tmux_session,
                    f"echo -e '{info_message}'",
                    "Enter",
                ],
                capture_output=True,
            )

            click.echo(f"🔗 Attaching to existing tmux session: {tmux_session}")
            os.execvp("tmux", ["tmux", "attach-session", "-t", tmux_session])
        else:
            # Create new session and run qsub inside it
            click.echo(
                f"🆕 Creating tmux session '{tmux_session}' and starting PBS job..."
            )
            os.execvp(
                "tmux",
                ["tmux", "new-session", "-s", tmux_session, qsub_cmd_str],
            )
    else:
        # No tmux - use os.execvp to replace this process with qsub -I
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
