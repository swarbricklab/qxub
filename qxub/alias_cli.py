"""
Alias CLI commands for qxub.
"""

import click

from .config import config_manager


@click.command(name="alias")
@click.argument("alias_name")
@click.argument("command_args", nargs=-1)
@click.option("--cmd", help="Override the alias command")
@click.option("--name", help="Override job name")
@click.option("--queue", help="Override queue")
@click.option("--project", help="Override project")
@click.option("--resources", multiple=True, help="Override resource requests")
@click.option("--joblog", help="Override joblog path")
@click.option("--out", help="Override stdout path")
@click.option("--err", help="Override stderr path")
@click.option("--pre", help="Override pre-command")
@click.option("--post", help="Override post-command")
# Conda-specific overrides
@click.option("--env", help="Override conda environment")
# Module-specific overrides
@click.option("--mod", help="Override single module to load")
@click.option("--mods", help="Override multiple modules to load (comma-separated)")
# Singularity-specific overrides
@click.option("--sif", help="Override singularity container")
@click.option("--bind", multiple=True, help="Override bind paths")
@click.option("--env-var", multiple=True, help="Override environment variables")
@click.pass_context
def alias_cli(ctx, alias_name: str, command_args: tuple, **overrides):
    """Execute an alias with optional overrides and command arguments.

    ⚠️  DEPRECATION WARNING: This command is deprecated.
    Use 'qxub exec --alias <alias_name>' instead.

    This command will be removed in a future version.

    Special management commands:
        qxub alias list - List all available aliases
        qxub alias test <alias_name> - Test an alias without executing

    Examples:
        qxub alias dvc_push
        qxub alias train_model --queue gpuvolta
        qxub alias analysis -- input.bam output.bam
        qxub alias quick_task --cmd "python script.py"
    """
    # Show deprecation warning
    import warnings

    warnings.warn(
        "qxub alias is deprecated. Use 'qxub exec --alias <alias_name>' instead. "
        "This command will be removed in a future version.",
        DeprecationWarning,
        stacklevel=2,
    )
    click.echo("⚠️  DEPRECATION WARNING: 'qxub alias' is deprecated.")
    click.echo("💡 Use 'qxub exec --alias <alias_name>' instead.")
    click.echo()

    # Handle special management commands
    if alias_name == "list":
        _handle_alias_list()
        return
    elif alias_name == "test":
        if not command_args:
            click.echo("❌ 'test' command requires an alias name")
            ctx.exit(2)
        _handle_alias_test(command_args[0])
        return

    # Handle normal alias execution
    _handle_alias_execution(ctx, alias_name, command_args, overrides)


def _handle_alias_list():
    """List all available aliases."""
    aliases = config_manager.list_aliases()
    if not aliases:
        click.echo("No aliases configured.")
        return

    click.echo("Available aliases:")
    for alias_name in aliases:
        alias_def = config_manager.get_alias(alias_name)
        if alias_def:
            # Determine execution context
            if alias_def.get("env"):
                context = f"conda (env: {alias_def['env']})"
            elif alias_def.get("mod"):
                context = f"module (mod: {alias_def['mod']})"
            elif alias_def.get("mods"):
                context = f"modules (mods: {alias_def['mods']})"
            elif alias_def.get("sif"):
                context = f"singularity (sif: {alias_def['sif']})"
            else:
                context = "default"

            cmd = alias_def.get("cmd", "(requires command args)")
            click.echo(f"  • {alias_name}: {context} - {cmd}")


def _handle_alias_test(alias_name: str):
    """Test an alias without executing it (dry run)."""
    # Get the alias definition
    alias_def = config_manager.get_alias(alias_name)
    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        click.echo("Available aliases:")
        for alias in config_manager.list_aliases():
            click.echo(f"  • {alias}")
        return

    click.echo(f"🧪 Testing alias: {alias_name}")

    # Show alias definition
    click.echo("📋 Alias definition:")
    for key, value in alias_def.items():
        if isinstance(value, (list, tuple)):
            click.echo(f"  {key}: {', '.join(str(v) for v in value)}")
        else:
            click.echo(f"  {key}: {value}")

    # Determine execution context
    execution_context = "default"
    if alias_def.get("env"):
        execution_context = f"conda (env: {alias_def['env']})"
    elif alias_def.get("mod"):
        execution_context = f"module (mod: {alias_def['mod']})"
    elif alias_def.get("mods"):
        execution_context = f"modules (mods: {alias_def['mods']})"
    elif alias_def.get("sif"):
        execution_context = f"singularity (sif: {alias_def['sif']})"

    click.echo("✅ Alias validation:")
    click.echo(f"  • Execution context: {execution_context}")
    click.echo(f"  • Command: {alias_def.get('cmd', '(none - requires command args)')}")

    if alias_def.get("sif") and not alias_def.get("sif"):
        click.echo("  ⚠️  Warning: No Singularity container specified")

    click.echo("🎉 Alias test completed")


def _handle_alias_execution(ctx, alias_name: str, command_args: tuple, overrides: dict):
    """Execute an alias with optional additional arguments."""
    # Get the alias definition
    alias_def = config_manager.get_alias(alias_name)
    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        click.echo("Available aliases:")
        for alias in config_manager.list_aliases():
            click.echo(f"  • {alias}")
        ctx.exit(2)

    # Handle legacy format migration
    if "subcommand" in alias_def:
        click.echo("⚠️  Migrating legacy alias format...")
        # Convert old subcommand format to new unified CLI format
        subcommand = alias_def.get("subcommand")
        if subcommand == "conda":
            alias_def["env"] = alias_def.get("env")
        elif subcommand == "mod":
            alias_def["mod"] = alias_def.get("mod")
        elif subcommand == "sing":
            alias_def["sif"] = alias_def.get("sif")
        # Remove legacy fields
        alias_def.pop("subcommand", None)

    # Remove None values from overrides, and empty tuples from Click multiple options
    clean_overrides = {
        k: v
        for k, v in overrides.items()
        if v is not None and not (isinstance(v, tuple) and len(v) == 0)
    }

    # Build the command line for the unified CLI
    cmd_args = ["qxub"]

    # Add execution context options (from alias, can be overridden)
    if clean_overrides.get("env") or alias_def.get("env"):
        cmd_args.extend(["--env", clean_overrides.get("env") or alias_def["env"]])
    elif clean_overrides.get("mod") or alias_def.get("mod"):
        cmd_args.extend(["--mod", clean_overrides.get("mod") or alias_def["mod"]])
    elif clean_overrides.get("mods") or alias_def.get("mods"):
        mods = clean_overrides.get("mods") or alias_def["mods"]
        if isinstance(mods, str):
            mods = mods.split(",")
        for mod in mods:
            cmd_args.extend(["--mod", mod.strip()])
    elif clean_overrides.get("sif") or alias_def.get("sif"):
        cmd_args.extend(["--sif", clean_overrides.get("sif") or alias_def["sif"]])

    # Add all PBS-related options from alias (can be overridden)
    pbs_options = ["walltime", "mem", "ncpus", "jobfs", "queue", "project", "name"]
    for option in pbs_options:
        value = clean_overrides.get(option) or alias_def.get(option)
        if value:
            cmd_args.extend([f"--{option}", str(value)])

    # Add resources (can be multiple, and can be overridden)
    resources = clean_overrides.get("resources") or alias_def.get("resources")
    if resources:
        if isinstance(resources, (list, tuple)):
            for resource in resources:
                cmd_args.extend(["--resources", resource])
        else:
            cmd_args.extend(["--resources", resources])

    # Add separator before command
    cmd_args.append("--")

    # Add the command from alias or override
    cmd = clean_overrides.get("cmd") or alias_def.get("cmd")
    if cmd:
        cmd_args.extend(cmd.split())

    # Add additional arguments provided by user
    cmd_args.extend(command_args)

    click.echo(f"🚀 Executing alias: {alias_name}")
    click.echo(f"📝 Command: {' '.join(cmd_args)}")

    # Import here to avoid circular imports
    from .cli import qxub as main_cli

    # Execute the command by invoking the main CLI
    try:
        # Remove 'qxub' from the beginning since we're calling it directly
        main_cli.main(cmd_args[1:], standalone_mode=False)
    except SystemExit as e:
        # Re-raise the exit code
        ctx.exit(e.code)
    except Exception as e:
        click.echo(f"❌ Error executing alias: {e}")
        ctx.exit(1)
