"""
Execution subcommand for qxub - handles all job execution contexts.

This module provides the 'exec' subcommand that consolidates all execution
options and contexts (conda, modules, singularity, default) into a single
comprehensive interface.
"""

import click

from .config import config_manager
from .config.handler import process_job_options
from .execution import (
    ExecutionContext,
    ExecutionMode,
    execute_unified,
    get_execution_mode,
)


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
# Workflow-friendly resource options
@click.option(
    "--mem",
    "--memory",
    help="Memory requirement (workflow-friendly). Examples: '4GB', '2000MB', '16g'. (default: configured). Converted to PBS format automatically.",
)
@click.option(
    "--runtime",
    "--time",
    help="Runtime/walltime limit (workflow-friendly). Examples: '2h', '30m', '1h30m', '02:30:00'. (default: configured). Converted to PBS format automatically.",
)
@click.option(
    "--cpus",
    "--threads",
    type=int,
    help="Number of CPU cores/threads (workflow-friendly). (default: configured). Converted to PBS ncpus automatically.",
)
@click.option(
    "--disk",
    "--jobfs",
    help="Local disk/jobfs requirement (workflow-friendly). Examples: '10GB', '500MB'. (default: configured). Converted to PBS format automatically.",
)
@click.option(
    "--volumes",
    "--storage",
    help="Storage volumes to mount (NCI format). Examples: 'gdata/a56', 'gdata/a56+gdata/px14'. (default: configured). Converted to PBS storage= automatically.",
)
@click.option(
    "-P", "--project", help="PBS project code (default: configured or $PROJECT)"
)
@click.option(
    "--dry",
    "--dry-run",
    is_flag=True,
    help="Show what would be executed without submitting",
)
@click.option("--quiet", is_flag=True, help="Suppress output")
@click.option(
    "--terse",
    is_flag=True,
    help="Terse output: emit only job ID and return immediately (for use in pipelines)",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (use -v, -vv, -vvv for more detail)",
)
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
@click.option("--platform", help="Target platform for execution (default: from config)")
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
@click.option("--alias", help="Use a predefined alias for execution settings")
@click.argument("command", nargs=-1, required=False)
@click.pass_context
def exec_cli(ctx, command, cmd, shortcut, alias, verbose, **options):
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

    \b
    Aliases:
        qxub exec --alias myalias -- additional args
        qxub exec --alias myalias  # if alias has default command

    Commands can be specified either after -- or using --cmd:

    \b
        qxub exec --env myenv -- python script.py arg1 arg2
        qxub exec --env myenv --cmd "python script.py arg1 arg2"
    """
    # Set up logging verbosity using the count-based system
    from .config.manager import setup_logging

    if options.get("quiet"):
        # Override verbose count when quiet is specified
        setup_logging(verbosity=0)  # ERROR level
    else:
        setup_logging(verbosity=verbose)

    # Handle command specification
    if cmd and command:
        raise click.ClickException(
            "Cannot specify both --cmd and positional command arguments"
        )

    if cmd:
        # Split the command string into components
        command = tuple(cmd.split())
    elif not command:
        # Commands are required unless using an alias with a default command
        if not alias:
            raise click.ClickException(
                "Must specify a command either after -- or using --cmd"
            )

    # Check for conflicting execution modes
    if shortcut and alias:
        raise click.ClickException(
            "Cannot specify both --shortcut and --alias. Use only one."
        )

    # Handle alias processing
    if alias:
        alias_def = config_manager.get_alias(alias)
        if not alias_def:
            available_aliases = config_manager.list_aliases()
            click.echo(f"❌ Alias '{alias}' not found")
            if available_aliases:
                click.echo("💡 Available aliases:")
                for alias_name in sorted(available_aliases):
                    click.echo(f"  • {alias_name}")
            else:
                click.echo("💡 No aliases defined yet.")
                click.echo("💡 Create one with: qxub config alias set <name> [options]")
            ctx.exit(2)

        # Apply alias settings to options (CLI options override alias settings)
        # Handle both hierarchical and flat alias structures

        # Flatten hierarchical alias structure if present
        flat_alias = {}

        if "main" in alias_def and isinstance(alias_def["main"], dict):
            # Hierarchical structure: extract main options
            main_section = alias_def["main"]
            flat_alias.update(main_section)

        if "subcommand" in alias_def and isinstance(alias_def["subcommand"], dict):
            # Hierarchical structure: extract subcommand options
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

        if "target" in alias_def and isinstance(alias_def["target"], dict):
            # Hierarchical structure: extract target command
            target_section = alias_def["target"]
            if "cmd" in target_section:
                flat_alias["cmd"] = target_section["cmd"]

        # Also handle flat structure for backward compatibility
        for key, value in alias_def.items():
            if key not in ["main", "subcommand", "target"]:
                flat_alias[key] = value

        # Apply flattened alias settings to options
        for key, value in flat_alias.items():
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

        # Set command from alias if no command was provided
        if not command and not cmd:
            alias_cmd = flat_alias.get("cmd")
            if alias_cmd:
                command = tuple(alias_cmd.split())

        click.echo(f"🎯 Using alias '{alias}'")

    # Handle shortcut processing
    elif shortcut:
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

    # Process workflow-friendly resource options with config defaults
    workflow_resources = []

    # Parse existing --resources to see what's already specified
    existing_resource_keys = set()
    if options["resources"]:
        for resource in options["resources"]:
            # Split by comma in case of comma-separated resources (e.g., "ncpus=48,mem=128GB")
            for part in resource.split(","):
                if "=" in part:
                    key = part.split("=", 1)[0].strip()
                    existing_resource_keys.add(key)

    # Map resource key aliases (e.g., ncpus -> cpus, walltime -> runtime)
    resource_key_map = {
        "ncpus": "cpus",
        "walltime": "runtime",
        "mem": "mem",
        "jobfs": "disk",
        "storage": "volumes",
    }

    # Check which resources are already specified (considering aliases)
    specified_resources = set()
    for key in existing_resource_keys:
        canonical_key = resource_key_map.get(key, key)
        specified_resources.add(canonical_key)

    # Helper function to get workflow-friendly config values with backward compatibility
    def get_workflow_resource_config(key):
        """
        Get workflow resource from config, checking new location first (defaults.KEY),
        then falling back to old root-level location with deprecation warning.
        """
        # Try new location: defaults.mem, defaults.cpus, etc.
        value = config_manager.get_config_value(f"defaults.{key}")
        if value is not None:
            return value

        # Fall back to old root-level location for backward compatibility
        value = config_manager.get_config_value(key)
        if value is not None:
            import logging

            logging.warning(
                f"Config key '{key}' at root level is deprecated. "
                f"Please move to 'defaults.{key}' in your config file."
            )
            return value

        return None

    # Collect all resource values (CLI args take precedence over config)
    # Only use config defaults if the resource wasn't already specified in --resources
    resource_values = {
        "mem": options.get("mem")
        or options.get("memory")
        or (
            get_workflow_resource_config("mem")
            if "mem" not in specified_resources
            else None
        ),
        "runtime": options.get("runtime")
        or options.get("time")
        or (
            get_workflow_resource_config("runtime")
            if "runtime" not in specified_resources
            else None
        ),
        "cpus": options.get("cpus")
        or options.get("threads")
        or (
            get_workflow_resource_config("cpus")
            if "cpus" not in specified_resources
            else None
        ),
        "disk": options.get("disk")
        or options.get("jobfs")
        or (
            get_workflow_resource_config("disk")
            if "disk" not in specified_resources
            else None
        ),
        "volumes": options.get("volumes")
        or options.get("storage")
        or (
            get_workflow_resource_config("volumes")
            if "volumes" not in specified_resources
            else None
        ),
    }

    # Only create mapper if we have any resource values to process
    if any(resource_values.values()):
        from .resources import ResourceMapper

        mapper = ResourceMapper()

        if resource_values["mem"]:
            mapper.add_memory(resource_values["mem"])

        if resource_values["runtime"]:
            mapper.add_runtime(resource_values["runtime"])

        if resource_values["cpus"]:
            mapper.add_cpus(resource_values["cpus"])

        if resource_values["disk"]:
            mapper.add_disk(resource_values["disk"])

        if resource_values["volumes"]:
            mapper.add_storage(resource_values["volumes"])

        workflow_resources.extend(mapper.get_pbs_resources())

    # Merge workflow resources with existing PBS resources
    all_resources = list(options["resources"]) if options["resources"] else []
    all_resources.extend(workflow_resources)

    # Track whether CPUs were explicitly specified (for graceful queue adjustment)
    cpus_explicit = (
        options.get("cpus") is not None
        or options.get("threads") is not None
        or any(r.startswith("ncpus=") for r in (options["resources"] or []))
    )

    # Extract PBS-specific options for processing
    params = {
        "resources": tuple(all_resources),  # Use merged resources
        "queue": options["queue"],
        "name": options["name"],
        "project": options["project"],
        "out": options["out"],
        "err": options["err"],
        "joblog": None,  # Will be set by config system
        "execdir": options["execdir"],
        "email": options["email"],
        "email_opts": options["email_opts"],
        "array": options["array"],
        "dry": options["dry"],
        "quiet": options["quiet"],
        "terse": options["terse"],
        "verbose": verbose,
        "cpus_explicit": cpus_explicit,  # Track for graceful adjustment
    }

    # Process configuration using the existing config system
    processed_params, qsub_options = process_job_options(params, config_manager)

    # Update ctx.obj with processed options
    if not hasattr(ctx, "obj") or ctx.obj is None:
        ctx.obj = {}
    ctx.obj.update(processed_params)
    ctx.obj["options"] = qsub_options

    # Platform selection and execution mode detection
    platform_name = (
        options.get("platform") or config_manager.get_default_platform_name()
    )
    platform_config = None
    execution_mode = ExecutionMode.LOCAL  # Default to local

    if platform_name:
        platform_config = config_manager.get_platform_config(platform_name)
        if not platform_config:
            # List available platforms for helpful error message
            available_platforms = config_manager.list_platforms()
            if available_platforms:
                platforms_str = ", ".join(available_platforms)
                raise click.ClickException(
                    f"Platform '{platform_name}' not found in configuration. "
                    f"Available platforms: {platforms_str}"
                )
            else:
                raise click.ClickException(
                    f"Platform '{platform_name}' not found. "
                    f"No platforms configured in {config_manager._get_user_config_dir() / 'config.yaml'}"
                )

        # Determine execution mode based on platform config
        execution_mode = get_execution_mode(platform_config)

        if execution_mode == ExecutionMode.REMOTE:
            # Remote execution path
            from .remote.command_builder import build_remote_command
            from .remote.platform_executor import PlatformRemoteExecutor

            remote_config_dict = platform_config.get("remote")
            if not remote_config_dict:
                raise click.ClickException(
                    f"Platform '{platform_name}' is configured for remote execution "
                    "but missing 'remote' configuration section"
                )

            # Build remote command that recreates the qxub invocation
            # Combine options for serialization (don't include --platform since we're already on the platform)
            remote_options = dict(options)
            remote_options.update(
                {
                    "verbose": verbose,
                    "dry": options.get("dry", False),
                    "quiet": options.get("quiet", False),
                    "terse": options.get("terse", False),
                }
            )

            # Serialize the execution context and options back to CLI arguments
            remote_command = build_remote_command(
                execution_context, remote_options, list(command)
            )

            if verbose >= 1:
                click.echo(
                    f"🌐 Remote execution on platform: {platform_name}", err=True
                )

            # Create and execute via SSH
            try:
                executor = PlatformRemoteExecutor(platform_name, remote_config_dict)
                exit_code = executor.execute(
                    remote_command,
                    stream_output=not options.get("quiet", False),
                    verbose=verbose,
                )
                ctx.exit(exit_code)
            except Exception as e:
                raise click.ClickException(f"Remote execution failed: {e}")

        # Local execution continues below
        if verbose >= 1:
            click.echo(f"📍 Executing on platform: {platform_name} (local)", err=True)

    # Execute the job using unified execution (local path)
    execute_unified(
        ctx,
        command,
        execution_context,
        template=options["template"],
        pre=options["pre"],
        post=options["post"],
        bind=options["bind"],
    )
