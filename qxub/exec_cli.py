"""
Execution subcommand for qxub - handles all job execution contexts.

This module provides the 'exec' subcommand that consolidates all execution
options and contexts (conda, modules, singularity, default) into a single
comprehensive interface.
"""

import click

from .config_handler import process_job_options
from .config_manager import config_manager
from .execution_context import ExecutionContext, execute_unified


def _get_shortcut_context_description(definition: dict) -> str:
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


@click.command(name="exec")
@click.option(
    "-l",
    "--resources",
    multiple=True,
    help="PBS resource specification (e.g., 'walltime=1:00:00,mem=4GB') (default: configured)",
)
@click.option("-q", "--queue", help="PBS queue name (default: configured or normal)")
@click.option(
    "-N", "--name", help="PBS job name (default: configured or qx-{timestamp})"
)
@click.option(
    "-P", "--project", help="PBS project code (default: configured or $PROJECT)"
)
@click.option("-v", "--variable", multiple=True, help="PBS environment variables")
@click.option(
    "--dry", is_flag=True, help="Show what would be executed without submitting"
)
@click.option("--quiet", is_flag=True, help="Suppress output")
@click.option(
    "--out",
    help="Output file path (default: configured or {log_dir}/{name}_{timestamp}.out)",
)
@click.option(
    "--err",
    help="Error file path (default: configured or {log_dir}/{name}_{timestamp}.err)",
)
@click.option("--execdir", help="Execution directory (default: current directory)")
@click.option(
    "--email", help="Email address for PBS notifications (default: configured)"
)
@click.option(
    "--email-opts", help="PBS email options (e.g., 'abe') (default: configured)"
)
@click.option("--array", help="Job array specification (e.g., '1-10' or '1-100:2')")
@click.option("--env", help="Conda environment name for execution")
@click.option(
    "--mod",
    "--module",
    help="Single module to load (can be used multiple times)",
    multiple=True,
)
@click.option("--mods", "--modules", help="Space-separated list of modules to load")
@click.option("--sif", "--container", help="Singularity container (.sif file) to use")
@click.option(
    "--default", is_flag=True, help="Use default execution (no special environment)"
)
@click.option("--bind", help="Singularity bind mounts (only used with --sif)")
@click.option("--template", help="Custom job script template file")
@click.option("--pre", help="Command to run before the main command")
@click.option("--post", help="Command to run after the main command")
@click.option("--cmd", help="Command to execute (alternative to positional arguments)")
@click.option("--shortcut", help="Use a predefined shortcut for execution settings")
@click.argument("command", nargs=-1, required=False)
@click.pass_context
def exec_cli(ctx, command, cmd, shortcut, **options):
    """
    Execute commands in various environments using PBS.

    This command supports multiple execution contexts:

    \b
    Conda environments:
        qxub exec --env myenv -- python script.py

    \b
    Environment modules:
        qxub exec --mod python3 --mod numpy -- python script.py
        qxub exec --mods "python3 numpy scipy" -- python script.py

    \b
    Singularity containers:
        qxub exec --sif container.sif -- python script.py
        qxub exec --sif container.sif --bind /data:/mnt -- ls /mnt

    \b
    Default execution:
        qxub exec --default -- echo "hello world"
        qxub exec -- echo "hello world"  # default is implicit

    \b
    Shortcuts:
        qxub exec --shortcut myshortcut -- additional args
        qxub exec -- myshortcut-command args  # automatic shortcut matching

    Commands can be specified either after -- or using --cmd:

    \b
        qxub exec --env myenv -- python script.py arg1 arg2
        qxub exec --env myenv --cmd "python script.py arg1 arg2"
    """
    # Handle command specification
    if cmd and command:
        raise click.ClickException(
            "Cannot specify both --cmd and positional command arguments"
        )

    if cmd:
        # Split the command string into components
        command = tuple(cmd.split())
    elif not command:
        raise click.ClickException(
            "Must specify a command either after -- or using --cmd"
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
                for name, definition in sorted(available_shortcuts.items()):
                    context = _get_shortcut_context_description(definition)
                    click.echo(f"  • {name}: {context}")
            else:
                click.echo("💡 No shortcuts defined yet.")
                click.echo(
                    "💡 Create one with: qxub shortcut set <name> --env <env> [options]"
                )
            ctx.exit(2)

        # Apply shortcut settings to options (CLI options override shortcut settings)
        for key, value in shortcut_def.items():
            # Map shortcut keys to option keys
            if key == "env" and not options["env"]:
                options["env"] = value
            elif key == "mod" and not options["mod"]:
                options["mod"] = (value,) if isinstance(value, str) else tuple(value)
            elif key == "mods" and not options["mods"]:
                options["mods"] = value
            elif key == "sif" and not options["sif"]:
                options["sif"] = value
            elif key == "queue" and not options["queue"]:
                options["queue"] = value
            elif key == "resources" and not options["resources"]:
                # Convert shortcut resources to the expected format
                if isinstance(value, (list, tuple)):
                    options["resources"] = tuple(value)
                else:
                    options["resources"] = (value,)
            elif key == "project" and not options["project"]:
                options["project"] = value
            elif key == "template" and not options["template"]:
                options["template"] = value
            elif key == "pre" and not options["pre"]:
                options["pre"] = value
            elif key == "post" and not options["post"]:
                options["post"] = value
            elif key == "bind" and not options["bind"]:
                options["bind"] = value

        click.echo(f"🎯 Using shortcut '{shortcut}'")
    elif not any(
        [
            options["env"],
            options["mod"],
            options["mods"],
            options["sif"],
            options["default"],
        ]
    ):
        # No shortcut and no explicit execution context - check for command-based shortcut matching
        from .shortcut_manager import ShortcutManager

        shortcut_manager = ShortcutManager()
        shortcut_match = shortcut_manager.find_shortcut(list(command))

        if shortcut_match:
            shortcut_name = shortcut_match["name"]
            shortcut_def = shortcut_match["definition"]
            remaining_args = shortcut_match["remaining_args"]

            # Apply shortcut settings (CLI options still override)
            for key, value in shortcut_def.items():
                if key == "env" and not options["env"]:
                    options["env"] = value
                elif key == "mod" and not options["mod"]:
                    options["mod"] = (
                        (value,) if isinstance(value, str) else tuple(value)
                    )
                elif key == "mods" and not options["mods"]:
                    options["mods"] = value
                elif key == "sif" and not options["sif"]:
                    options["sif"] = value
                elif key == "queue" and not options["queue"]:
                    options["queue"] = value
                elif key == "resources" and not options["resources"]:
                    if isinstance(value, (list, tuple)):
                        options["resources"] = tuple(value)
                    else:
                        options["resources"] = (value,)
                elif key == "project" and not options["project"]:
                    options["project"] = value

            # For shortcuts, keep the full original command (don't strip the prefix)
            # The shortcut defines the execution context, not the command itself

            click.echo(f"🎯 Found shortcut '{shortcut_name}' for command")
        # If no shortcut found, continue with default execution

    # Validate execution contexts
    execution_contexts = [
        options["env"],
        options["mod"] or options["mods"],
        options["sif"],
        options["default"],
    ]

    active_contexts = sum(bool(x) for x in execution_contexts)
    if active_contexts > 1:
        raise click.ClickException(
            "Cannot specify multiple execution contexts. "
            "Use only one of: --env, --mod/--mods, --sif, --default"
        )

    # Determine execution context
    if options["env"]:
        execution_context = ExecutionContext("conda", options["env"], "conda")
    elif options["mod"] or options["mods"]:
        # Combine --mod and --mods options
        modules = list(options["mod"]) if options["mod"] else []
        if options["mods"]:
            modules.extend(options["mods"].split())
        execution_context = ExecutionContext("module", modules, "module")
    elif options["sif"]:
        execution_context = ExecutionContext(
            "singularity", options["sif"], "singularity"
        )
    else:
        # Default execution (explicit --default or implicit)
        execution_context = ExecutionContext("default", None, "default")

    # Extract PBS-specific options for processing
    params = {
        "resources": options["resources"],
        "queue": options["queue"],
        "name": options["name"],
        "project": options["project"],
        "variable": options["variable"],
        "out": options["out"],
        "err": options["err"],
        "joblog": None,  # Will be set by config system
        "execdir": options["execdir"],
        "email": options["email"],
        "email_opts": options["email_opts"],
        "array": options["array"],
        "dry": options["dry"],
        "quiet": options["quiet"],
    }

    # Process configuration using the existing config system
    processed_params, qsub_options = process_job_options(params, config_manager)

    # Update ctx.obj with processed options
    if not hasattr(ctx, "obj") or ctx.obj is None:
        ctx.obj = {}
    ctx.obj.update(processed_params)
    ctx.obj["options"] = qsub_options

    # Execute the job using unified execution
    execute_unified(
        ctx,
        command,
        execution_context,
        template=options["template"],
        pre=options["pre"],
        post=options["post"],
        bind=options["bind"],
    )
