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


@click.command(name="exec")
@click.option(
    "-l",
    "--resource",
    help="PBS resource specification (e.g., 'walltime=1:00:00,mem=4GB')",
)
@click.option("-q", "--queue", help="PBS queue name")
@click.option("-N", "--job-name", help="PBS job name")
@click.option("-P", "--project", help="PBS project code")
@click.option("-v", "--variable", multiple=True, help="PBS environment variables")
@click.option(
    "--dry", is_flag=True, help="Show what would be executed without submitting"
)
@click.option("--quiet", is_flag=True, help="Suppress output")
@click.option("--out", help="Output file path")
@click.option("--err", help="Error file path")
@click.option("--execdir", help="Execution directory")
@click.option("--email", help="Email address for PBS notifications")
@click.option("--email-opts", help="PBS email options (e.g., 'abe')")
@click.option("--after", help="Job dependency (run after specified job)")
@click.option("--afterok", help="Job dependency (run after successful completion)")
@click.option("--afternotok", help="Job dependency (run after failed completion)")
@click.option("--afterany", help="Job dependency (run after any completion)")
@click.option("--hold", is_flag=True, help="Submit job in held state")
@click.option("--array", help="Job array specification")
@click.option("--priority", help="Job priority")
@click.option("--checkpoint", help="Checkpoint options")
@click.option("--interactive", is_flag=True, help="Interactive job")
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
@click.argument("command", nargs=-1, required=False)
@click.pass_context
def exec_cli(ctx, command, cmd, **options):
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
        "resources": (options["resource"],) if options["resource"] else (),
        "queue": options["queue"],
        "name": options["job_name"],
        "project": options["project"],
        "variable": options["variable"],
        "out": options["out"],
        "err": options["err"],
        "joblog": None,  # Will be set by defaults
        "execdir": options["execdir"],
        "email": options["email"],
        "email_opts": options["email_opts"],
        "after": options["after"],
        "afterok": options["afterok"],
        "afternotok": options["afternotok"],
        "afterany": options["afterany"],
        "hold": options["hold"],
        "array": options["array"],
        "priority": options["priority"],
        "checkpoint": options["checkpoint"],
        "interactive": options["interactive"],
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
