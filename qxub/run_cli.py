"""
qxub run command for unified execution: shortcuts, aliases, and direct execution.

The run command handles all qxub execution modes:
1. Shortcut execution with --shortcut flag
2. Direct execution with execution context options
3. Optional alias overlays for both modes
"""

import os

import click

from .config_manager import config_manager


@click.command(name="run")
@click.argument("command_parts", nargs=-1, required=True)
# Control flags
@click.option(
    "--shortcut", help="Use specific named shortcut (otherwise direct execution)"
)
@click.option("--alias", help="Apply alias settings to shortcut")
# Execution directory and basic options
@click.option(
    "--execdir",
    default=os.getcwd(),
    help="Execution directory (default: current directory)",
)
@click.option(
    "--out",
    help="STDOUT log file (default: configured or /scratch/$PROJECT/$USER/qt/timestamp/out)",
)
@click.option(
    "--err",
    help="STDERR log file (default: configured or /scratch/$PROJECT/$USER/qt/timestamp/err)",
)
@click.option("--joblog", help="PBS Pro job log (default: configured or {name}.log)")
@click.option(
    "--dry",
    "--dry-run",
    is_flag=True,
    default=False,
    help="Generate job submission command but don't submit",
)
@click.option("--quiet", is_flag=True, default=False, help="Display no output")
# PBS job options
@click.option(
    "-l", "--resources", multiple=True, help="Job resource (default: configured)"
)
@click.option(
    "-q",
    "--queue",
    help="Job queue (default: configured or normal, use 'auto' for intelligent selection)",
)
@click.option("-N", "--name", help="Job name (default: configured or qt)")
@click.option(
    "-P", "--project", help="PBS project code (default: configured or $PROJECT)"
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (use -v, -vv, -vvv for more detail)",
)
# Execution context options (mutually exclusive)
@click.option("--env", "--conda", help="Conda environment for execution")
@click.option("--mod", multiple=True, help="Environment module to load (repeatable)")
@click.option("--mods", "--modules", help="Comma-separated list of environment modules")
@click.option("--sif", "--sing", "--singularity", help="Singularity container image")
# Additional execution options
@click.option("--bind", help="Singularity bind mounts")
@click.option(
    "--template", help="Jobscript template (optional - for further customization)"
)
@click.option("--pre", help="Command to run before the main command")
@click.option("--post", help="Command to run after the main command")
@click.option(
    "--cmd",
    help="Command to execute (supports ${var} for submission-time and ${{var}} for execution-time variables)",
)
@click.pass_context
def run_cli(
    ctx,
    command_parts: tuple,
    shortcut: str,
    alias: str,
    execdir: str,
    out: str,
    err: str,
    joblog: str,
    dry: bool,
    quiet: bool,
    resources: tuple,
    queue: str,
    name: str,
    project: str,
    verbose: int,
    env: str,
    mod: tuple,
    mods: str,
    sif: str,
    bind: str,
    template: str,
    pre: str,
    post: str,
    cmd: str,
):
    """Execute commands using shortcuts or direct execution contexts.

    Two execution modes:
    1. Shortcut mode (--shortcut NAME): Use predefined shortcut configuration
    2. Direct mode (no --shortcut): Execute with provided options directly

    Examples:
        # Shortcut execution
        qxub run --shortcut dvc_doctor -- dvc doctor --help
        qxub run --shortcut jupyter --alias quick -- jupyter lab

        # Direct execution
        qxub run --env myenv -- python script.py
        qxub run --mod python3 --mod gcc -- make
        qxub run --queue express --sif container.sif -- ./app

        # Complex commands with variables
        qxub run --env myenv --cmd "python script.py --input \${HOME}/data.txt"

    Bash alias suggestion: alias qr="qxub run"
    """
    if not command_parts and not cmd:
        click.echo("❌ No command provided")
        click.echo("💡 Use: qxub run [options] -- <command> [args...]")
        click.echo('💡 Or:  qxub run [options] --cmd "<command>"')
        ctx.exit(2)

    # Handle command specification - --cmd vs -- syntax
    if cmd and command_parts:
        click.echo(
            "Error: Cannot specify both --cmd and command arguments after --\n"
            "Use either:\n"
            '  qxub run --env base --cmd "command with ${vars}"\n'
            "  qxub run --env base -- command args",
            err=True,
        )
        ctx.exit(2)
    elif cmd:
        # Convert --cmd string to tuple for compatibility
        command = (cmd,)
    else:
        # Get remaining arguments (these would be the command to execute)
        command = tuple(command_parts) if command_parts else tuple()

    if shortcut:
        # Shortcut execution mode
        _execute_with_shortcut(ctx, shortcut, alias, command, locals())
    else:
        # Direct execution mode
        _execute_directly(ctx, command, locals())


def _show_no_shortcut_found(command_parts: tuple) -> None:
    """Show helpful message when no shortcut matches."""
    from .shortcut_manager import ShortcutManager

    command_str = " ".join(command_parts)
    click.echo(f"❌ No shortcut found for command: '{command_str}'")

    # Show available shortcuts
    shortcut_manager = ShortcutManager()
    shortcuts = shortcut_manager.list_shortcuts()
    if shortcuts:
        click.echo("💡 Available shortcuts:")
        for name, definition in sorted(shortcuts.items()):
            cmd = definition.get("cmd", name)
            context = _get_context_description(definition)
            click.echo(f"  • {name}: {cmd} ({context})")
    else:
        click.echo("💡 No shortcuts defined yet.")

    click.echo(
        f'💡 Create one with: qxub config shortcut set "{command_str}" --env <env> -- {command_str}'
    )


def _get_context_description(definition: dict) -> str:
    """Get human-readable description of execution context."""
    if definition.get("env"):
        return f"conda: {definition['env']}"
    elif definition.get("mod"):
        return f"module: {definition['mod']}"
    elif definition.get("mods"):
        return f"modules: {definition['mods']}"
    elif definition.get("sif"):
        return f"singularity: {definition['sif']}"
    else:
        return "default"


def _execute_shortcut(ctx, shortcut_match: dict, overrides: dict) -> None:
    """Execute shortcut with CLI overrides only."""
    shortcut_name = shortcut_match["name"]
    shortcut_def = shortcut_match["definition"]
    remaining_args = shortcut_match["remaining_args"]

    # Report shortcut usage
    click.echo(f"🎯 Using shortcut '{shortcut_name}'")

    # Merge shortcut settings with CLI overrides
    merged_settings = _merge_shortcut_and_overrides(shortcut_def, overrides)

    # Show effective settings
    _show_effective_settings(merged_settings)

    # Execute command
    _execute_with_settings(ctx, merged_settings, remaining_args)


def _execute_shortcut_with_alias(
    ctx, shortcut_match: dict, alias_name: str, overrides: dict
) -> None:
    """Execute shortcut combined with alias settings."""
    shortcut_name = shortcut_match["name"]
    shortcut_def = shortcut_match["definition"]
    remaining_args = shortcut_match["remaining_args"]

    # Get alias definition
    alias_def = config_manager.get_alias(alias_name)
    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        click.echo("💡 Available aliases:")
        for alias in config_manager.list_aliases():
            click.echo(f"  • {alias}")
        ctx.exit(2)

    # Report combination usage
    click.echo(f"🎯 Using shortcut '{shortcut_name}' + alias '{alias_name}'")

    # Merge settings with precedence: CLI > Alias > Shortcut
    merged_settings = _merge_shortcut_alias_and_overrides(
        shortcut_def, alias_def, overrides
    )

    # Show effective settings
    _show_effective_settings(merged_settings)

    # Execute command
    _execute_with_settings(ctx, merged_settings, remaining_args)


def _merge_shortcut_and_overrides(shortcut_def: dict, overrides: dict) -> dict:
    """Merge shortcut definition with CLI overrides."""
    merged = {}

    # Start with shortcut settings (exclude cmd)
    for key, value in shortcut_def.items():
        if key != "cmd":
            merged[key] = value

    # Apply CLI overrides (higher precedence)
    for key, value in overrides.items():
        if value is not None and not (isinstance(value, tuple) and len(value) == 0):
            merged[key] = value

    return merged


def _merge_shortcut_alias_and_overrides(
    shortcut_def: dict, alias_def: dict, overrides: dict
) -> dict:
    """Merge shortcut, alias, and CLI overrides with proper precedence."""
    merged = {}

    # Start with shortcut settings (lowest precedence, exclude cmd)
    for key, value in shortcut_def.items():
        if key != "cmd":
            merged[key] = value

    # Apply alias settings (medium precedence, exclude cmd to favor shortcut command)
    for key, value in alias_def.items():
        if key != "cmd":
            merged[key] = value

    # Apply CLI overrides (highest precedence)
    for key, value in overrides.items():
        if value is not None and not (isinstance(value, tuple) and len(value) == 0):
            merged[key] = value

    return merged


def _show_effective_settings(settings: dict) -> None:
    """Display effective settings to user."""
    setting_parts = []

    # Execution context
    if settings.get("env"):
        setting_parts.append(f"env={settings['env']}")
    elif settings.get("mod"):
        setting_parts.append(f"mod={settings['mod']}")
    elif settings.get("mods"):
        setting_parts.append(f"mods={settings['mods']}")
    elif settings.get("sif"):
        setting_parts.append(f"sif={settings['sif']}")

    # PBS options
    if settings.get("queue"):
        setting_parts.append(f"queue={settings['queue']}")
    if settings.get("resources"):
        resources = settings["resources"]
        if isinstance(resources, (list, tuple)):
            setting_parts.append(f"resources={','.join(resources)}")
        else:
            setting_parts.append(f"resources={resources}")

    if setting_parts:
        click.echo(f"🔧 Settings: {', '.join(setting_parts)}")


def _execute_with_settings(ctx, settings: dict, remaining_args: list) -> None:
    """Execute command with merged settings."""
    # Import here to avoid circular imports
    from .cli import qxub as main_cli

    # Build command arguments for main CLI
    cmd_args = ["qxub"]

    # Add execution context
    if settings.get("env"):
        cmd_args.extend(["--env", settings["env"]])
    elif settings.get("mod"):
        cmd_args.extend(["--mod", settings["mod"]])
    elif settings.get("mods"):
        mods = settings["mods"]
        if isinstance(mods, str):
            mods = mods.split(",")
        for mod in mods:
            cmd_args.extend(["--mod", mod.strip()])
    elif settings.get("sif"):
        cmd_args.extend(["--sif", settings["sif"]])
    else:
        # Use --default for explicit default execution
        cmd_args.append("--default")

    # Add PBS options
    pbs_options = ["name", "queue", "project", "joblog", "out", "err", "pre", "post"]
    for option in pbs_options:
        value = settings.get(option)
        if value:
            cmd_args.extend([f"--{option}", str(value)])

    # Add resources
    resources = settings.get("resources")
    if resources:
        if isinstance(resources, (list, tuple)):
            for resource in resources:
                cmd_args.extend(["--resources", resource])
        else:
            cmd_args.extend(["--resources", resources])

    # Add bind paths for singularity
    bind_paths = settings.get("bind")
    if bind_paths:
        if isinstance(bind_paths, (list, tuple)):
            for bind_path in bind_paths:
                cmd_args.extend(["--bind", bind_path])
        else:
            cmd_args.extend(["--bind", bind_paths])

    # Add environment variables
    env_vars = settings.get("env-var")
    if env_vars:
        if isinstance(env_vars, (list, tuple)):
            for env_var in env_vars:
                cmd_args.extend(["--env-var", env_var])
        else:
            cmd_args.extend(["--env-var", env_vars])

    # Add separator and remaining command arguments
    cmd_args.append("--")
    cmd_args.extend(remaining_args)

    # Execute via main CLI
    try:
        # Remove 'qxub' from beginning since we're calling directly
        main_cli(cmd_args[1:], standalone_mode=False)
    except SystemExit as e:
        # Re-raise the exit code
        ctx.exit(e.code)
    except Exception as e:
        click.echo(f"❌ Error executing shortcut: {e}")
        ctx.exit(1)


def _execute_with_shortcut(
    ctx, shortcut_name: str, alias_name: str, command: tuple, options: dict
):
    """Execute using a named shortcut with optional alias overlay."""
    from .shortcut_manager import ShortcutManager

    # Get shortcut definition
    shortcut_manager = ShortcutManager()
    shortcuts = shortcut_manager.list_shortcuts()

    if shortcut_name not in shortcuts:
        click.echo(f"❌ Shortcut '{shortcut_name}' not found")
        click.echo("💡 Available shortcuts:")
        for name, definition in sorted(shortcuts.items()):
            cmd = definition.get("cmd", name)
            context = _get_context_description(definition)
            click.echo(f"  • {name}: {cmd} ({context})")
        ctx.exit(2)

    shortcut_def = shortcuts[shortcut_name]

    # Report shortcut usage
    if alias_name:
        click.echo(f"🎯 Using shortcut '{shortcut_name}' + alias '{alias_name}'")

        # Get alias definition
        alias_def = config_manager.get_alias(alias_name)
        if not alias_def:
            click.echo(f"❌ Alias '{alias_name}' not found")
            click.echo("💡 Available aliases:")
            for alias in config_manager.list_aliases():
                click.echo(f"  • {alias}")
            ctx.exit(2)

        # Merge settings with precedence: CLI > Alias > Shortcut
        merged_settings = _merge_shortcut_alias_and_overrides(
            shortcut_def, alias_def, options
        )
    else:
        click.echo(f"🎯 Using shortcut '{shortcut_name}'")
        # Merge shortcut settings with CLI overrides
        merged_settings = _merge_shortcut_and_overrides(shortcut_def, options)

    # Show effective settings
    _show_effective_settings(merged_settings)

    # Execute command
    _execute_with_settings(ctx, merged_settings, list(command))


def _execute_directly(ctx, command: tuple, options: dict):
    """Execute command directly using provided options without shortcuts."""
    # Import execution functions from main CLI
    from .execution import validate_execution_context

    # Extract execution options
    conda_env = options.get("env")
    mod = options.get("mod")
    mods = options.get("mods")
    sif = options.get("sif")

    # Handle module options: --mod (multiple) vs --mods (comma-separated)
    module_list = None
    if mod:
        module_list = list(mod)  # --mod can be used multiple times
    elif mods:
        module_list = [m.strip() for m in mods.split(",")]

    container = sif

    # Check if any execution context is specified
    execution_contexts = [conda_env, module_list, container]
    has_execution_context = any(execution_contexts)

    if has_execution_context:
        # Validate mutual exclusivity
        if sum(bool(x) for x in execution_contexts) > 1:
            click.echo("Error: Cannot specify multiple execution contexts", err=True)
            ctx.exit(2)
    else:
        # No execution context - use default
        has_execution_context, context_type = validate_execution_context(
            True, conda_env, module_list, container  # default=True
        )

    click.echo("🚀 Direct execution mode")

    # Build execution options for the unified execution system
    execution_options = {
        "execdir": options.get("execdir"),
        "out": options.get("out"),
        "err": options.get("err"),
        "joblog": options.get("joblog"),
        "dry": options.get("dry"),
        "quiet": options.get("quiet"),
        "resources": options.get("resources"),
        "queue": options.get("queue"),
        "name": options.get("name"),
        "project": options.get("project"),
        "verbose": options.get("verbose"),
        "template": options.get("template"),
        "pre": options.get("pre"),
        "post": options.get("post"),
    }

    # Import and use execution functions from main CLI
    # For now, delegate to main CLI via command construction
    _execute_via_main_cli(
        ctx, command, conda_env, module_list, container, execution_options
    )


def _execute_via_main_cli(
    ctx,
    command: tuple,
    conda_env: str,
    module_list: list,
    container: str,
    options: dict,
):
    """Execute by delegating to main CLI (temporary approach)."""
    from .cli import qxub as main_cli

    # Build command line arguments to pass to main CLI
    cmd_args = ["qxub"]

    # Add execution context
    if conda_env:
        cmd_args.extend(["--env", conda_env])
    elif module_list:
        for mod in module_list:
            cmd_args.extend(["--mod", mod])
    elif container:
        cmd_args.extend(["--sif", container])

    # Add PBS and other options
    option_mapping = {
        "execdir": "--execdir",
        "out": "--out",
        "err": "--err",
        "joblog": "--joblog",
        "dry": "--dry",
        "quiet": "--quiet",
        "queue": "--queue",
        "name": "--name",
        "project": "--project",
        "template": "--template",
        "pre": "--pre",
        "post": "--post",
    }

    for key, flag in option_mapping.items():
        value = options.get(key)
        if value is not None and value is not False:
            if isinstance(value, bool) and value:
                cmd_args.append(flag)
            elif not isinstance(value, bool):
                cmd_args.extend([flag, str(value)])

    # Add verbose flags
    verbose = options.get("verbose", 0)
    if verbose:
        cmd_args.extend(["-v"] * verbose)

    # Add resources
    resources = options.get("resources")
    if resources:
        for resource in resources:
            cmd_args.extend(["-l", resource])

    # Add separator and command arguments
    cmd_args.append("--")
    cmd_args.extend(command)

    # Execute via main CLI
    try:
        # Remove 'qxub' from beginning since we're calling directly
        main_cli(cmd_args[1:], standalone_mode=False)
    except SystemExit as e:
        # Re-raise the exit code
        ctx.exit(e.code)
    except Exception as e:
        click.echo(f"❌ Error executing command: {e}")
        ctx.exit(1)
