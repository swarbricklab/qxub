"""
Alias execution CLI for qxub.

Provides commands for executing aliases with optional overrides and command arguments.
"""

# pylint: disable=import-outside-toplevel
# pylint: disable=import-outside-toplevel,cyclic-import
import click
from rich.console import Console

from .config_manager import config_manager


console = Console()


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
def alias_cli(
    ctx, alias_name: str, command_args: tuple, **overrides
):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """Execute an alias with optional overrides and command arguments.

    Examples:
        qxub alias dvc_push
        qxub alias train_model --queue gpuvolta
        qxub alias analysis -- input.bam output.bam
        qxub alias quick_task --cmd "python script.py"
        qxub alias module_task --mod python3
        qxub alias multi_task --mods "python3,samtools,gcc"
    """
    # Get the alias definition
    alias_def = config_manager.get_alias(alias_name)
    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        click.echo("Available aliases:")
        for alias in config_manager.list_aliases():
            click.echo(f"  • {alias}")
        ctx.exit(2)

    # Parse the new structured alias format
    main_options = alias_def.get("main", {})
    subcommand_def = alias_def.get("subcommand", {})
    target_def = alias_def.get("target", {})

    # Check for legacy format and migrate
    if "subcommand" in alias_def and isinstance(alias_def["subcommand"], str):
        # Legacy format - migrate to new structure
        subcommand_type = alias_def["subcommand"]
        main_options = {
            k: v
            for k, v in alias_def.items()
            if k not in ["subcommand", "env", "mod", "sif", "bind", "env_var"]
        }
        subcommand_def = {
            "type": subcommand_type,
            **{
                k: v
                for k, v in alias_def.items()
                if k in ["env", "mod", "sif", "bind", "env_var"]
            },
        }
        target_def = {"cmd": alias_def.get("cmd")}

    # Get subcommand type
    subcommand_type = subcommand_def.get("type")
    if not subcommand_type:
        click.echo(f"❌ Alias '{alias_name}' missing subcommand type")
        ctx.exit(2)

    if subcommand_type not in ["conda", "module", "sing"]:
        click.echo(
            f"❌ Invalid subcommand type in alias '{alias_name}': {subcommand_type}"
        )
        ctx.exit(2)

    # Prepare the command
    cmd = overrides.get("cmd") or target_def.get("cmd")
    if not cmd and not command_args:
        click.echo(f"❌ No command specified for alias '{alias_name}'")
        click.echo(
            "Either define 'cmd' in the target section or provide command arguments"
        )
        ctx.exit(2)

    # Append command arguments if provided
    if command_args:
        if cmd:
            cmd = f"{cmd} {' '.join(command_args)}"
        else:
            cmd = " ".join(command_args)

    # Remove None values from overrides, and empty tuples from Click multiple options
    clean_overrides = {
        k: v
        for k, v in overrides.items()
        if v is not None and not (isinstance(v, tuple) and len(v) == 0)
    }

    # Merge alias options with CLI overrides
    # Main options (for qxub command level)
    main_opts = main_options.copy()
    for key in [
        "name",
        "queue",
        "project",
        "joblog",
        "out",
        "err",
        "pre",
        "post",
        "resources",
    ]:
        if key in clean_overrides:
            main_opts[key] = clean_overrides[key]

    # Subcommand options
    sub_opts = subcommand_def.copy()
    for key in ["env", "mod", "mods", "sif", "bind", "env_var"]:
        if key in clean_overrides:
            sub_opts[key] = clean_overrides[key]

    # Build the qxub command with proper option ordering
    qxub_cmd = ["qxub"]

    # Add main options first
    if main_opts.get("name"):
        qxub_cmd.extend(["--name", main_opts["name"]])
    if main_opts.get("queue"):
        qxub_cmd.extend(["--queue", main_opts["queue"]])
    if main_opts.get("project"):
        qxub_cmd.extend(["--project", main_opts["project"]])
    if main_opts.get("joblog"):
        qxub_cmd.extend(["--joblog", main_opts["joblog"]])
    if main_opts.get("out"):
        qxub_cmd.extend(["--out", main_opts["out"]])
    if main_opts.get("err"):
        qxub_cmd.extend(["--err", main_opts["err"]])
    if main_opts.get("pre"):
        qxub_cmd.extend(["--pre", main_opts["pre"]])
    if main_opts.get("post"):
        qxub_cmd.extend(["--post", main_opts["post"]])

    # Add resources (can be multiple)
    if main_opts.get("resources"):
        resources = main_opts["resources"]
        # Handle both tuples (from Click multiple=True) and lists
        if isinstance(resources, (list, tuple)):
            for resource in resources:
                qxub_cmd.extend(["--resources", resource])
        else:
            qxub_cmd.extend(["--resources", resources])

    # Add the subcommand type
    qxub_cmd.append(subcommand_type)

    # Add subcommand-specific options
    if subcommand_type == "conda":
        if sub_opts.get("env"):
            qxub_cmd.extend(["--env", sub_opts["env"]])
        if sub_opts.get("pre"):
            qxub_cmd.extend(["--pre", sub_opts["pre"]])
        if sub_opts.get("post"):
            qxub_cmd.extend(["--post", sub_opts["post"]])
    elif subcommand_type == "module":
        if sub_opts.get("mod"):
            qxub_cmd.extend(["--mod", sub_opts["mod"]])
        if sub_opts.get("mods"):
            for mod in sub_opts["mods"]:
                qxub_cmd.extend(["--mod", mod])
        if sub_opts.get("pre"):
            qxub_cmd.extend(["--pre", sub_opts["pre"]])
        if sub_opts.get("post"):
            qxub_cmd.extend(["--post", sub_opts["post"]])
    elif subcommand_type == "sing":
        if sub_opts.get("sif"):
            qxub_cmd.extend(["--sif", sub_opts["sif"]])
        if sub_opts.get("bind"):
            bind_paths = sub_opts["bind"]
            # Handle both tuples (from Click multiple=True) and lists
            if isinstance(bind_paths, (list, tuple)):
                for bind in bind_paths:
                    qxub_cmd.extend(["--bind", bind])
            else:
                qxub_cmd.extend(["--bind", bind_paths])
        if sub_opts.get("env_var"):
            env_vars = sub_opts["env_var"]
            # Handle both tuples (from Click multiple=True) and lists
            if isinstance(env_vars, (list, tuple)):
                for env_var in env_vars:
                    qxub_cmd.extend(["--env-var", env_var])
            else:
                qxub_cmd.extend(["--env-var", env_vars])

    # Add the command to execute
    qxub_cmd.append(cmd)

    # Show what would be executed
    click.echo(f"🏃 Executing alias '{alias_name}' with command:")
    # Ensure all elements are strings before joining
    qxub_cmd_str = []
    for item in qxub_cmd:
        if isinstance(item, (list, tuple)):
            qxub_cmd_str.extend(str(x) for x in item)
        else:
            qxub_cmd_str.append(str(item))
    click.echo(f"   {' '.join(qxub_cmd_str)}")
    click.echo()

    # Execute the command
    import subprocess

    result = subprocess.run(qxub_cmd_str)

    # Exit with the same code as the executed command
    exit(result.returncode)


# Test command for aliases
@click.command(name="alias-test")
@click.argument("alias_name")
@click.pass_context
def alias_test_cli(ctx, alias_name: str):
    """Test an alias without executing it (dry run)."""
    # pylint: disable=duplicate-code
    # Get the alias definition
    alias_def = config_manager.get_alias(alias_name)
    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        ctx.exit(2)

    click.echo(f"🧪 Testing alias: {alias_name}")

    # Show resolved configuration
    resolved_options = config_manager.resolve_options({}, alias_name)

    click.echo("📋 Resolved configuration:")
    for key, value in resolved_options.items():
        if isinstance(value, dict):
            click.echo(f"  {key}:")
            for subkey, subvalue in value.items():
                click.echo(f"    {subkey}: {subvalue}")
        else:
            click.echo(f"  {key}: {value}")

    # Validate the alias
    subcommand = alias_def.get("subcommand")
    cmd = alias_def.get("cmd")

    click.echo("✅ Alias validation:")
    click.echo(f"  • Subcommand: {subcommand}")
    click.echo(f"  • Command: {cmd or '(none - requires command args)'}")

    if subcommand == "sing":
        sif = resolved_options.get("sing", {}).get("sif")
        if not sif:
            click.echo("  ⚠️  Warning: No Singularity container specified")
        else:
            click.echo(f"  • Container: {sif}")

    click.echo("🎉 Alias test completed")
