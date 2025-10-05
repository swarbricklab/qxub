"""
Configuration management CLI for qxub.

Provides commands for viewing, editing, and managing qxub configuration
including defaults and aliases.
"""

# pylint: disable=broad-exception-caught,protected-access,raise-missing-from,unnecessary-pass,unused-variable
import os
import subprocess
from typing import Optional
import click
from rich.console import Console
from rich.table import Table
from rich.syntax import Syntax
from omegaconf import OmegaConf

from .config_manager import config_manager


console = Console()


@click.group(name="config")
def config_cli():
    """Manage qxub configuration and aliases."""
    pass


@config_cli.command()
@click.argument("key_path", required=False)
def get(key_path: Optional[str]):
    """Get configuration value by key path (e.g., 'defaults.name')."""
    if key_path:
        value = config_manager.get_config_value(key_path)
        if value is None:
            click.echo(f"Configuration key '{key_path}' not found")
            raise click.Abort()
        click.echo(value)
    else:
        # Show all config
        if config_manager.merged_config:
            yaml_str = OmegaConf.to_yaml(config_manager.merged_config)
            syntax = Syntax(yaml_str, "yaml", theme="monokai", line_numbers=True)
            console.print(syntax)
        else:
            click.echo("No configuration found")


@config_cli.command()
@click.argument("key_path")
@click.argument("value")
def set_config(key_path: str, value: str):
    """Set configuration value by key path (e.g., 'defaults.name' 'myjob')."""
    try:
        # Try to parse as YAML to handle different types
        try:
            # Create a temporary YAML document to parse the value
            yaml_doc = f"temp: {value}"
            parsed_config = OmegaConf.create(yaml_doc)
            parsed_value = parsed_config.temp
        except Exception:
            # If parsing fails, treat as string
            parsed_value = value

        config_manager.set_user_config_value(key_path, parsed_value)
        click.echo(f"✅ Set {key_path} = {parsed_value}")
    except Exception as e:
        click.echo(f"❌ Error setting configuration: {e}")
        raise click.Abort()


@config_cli.command()
@click.argument("section", required=False)
def list_config(section: Optional[str]):
    """List configuration values, optionally filtered by section."""
    if not config_manager.merged_config:
        click.echo("No configuration found")
        return

    if section:
        section_config = config_manager.get_config_value(section)
        if section_config is None:
            click.echo(f"Section '{section}' not found")
            raise click.Abort()

        if isinstance(section_config, dict):
            config_dict = {section: section_config}
        else:
            config_dict = {section: section_config}
    else:
        config_dict = OmegaConf.to_container(config_manager.merged_config, resolve=True)

    yaml_str = OmegaConf.to_yaml(OmegaConf.create(config_dict))
    syntax = Syntax(yaml_str, "yaml", theme="monokai", line_numbers=True)
    console.print(syntax)


@config_cli.command()
def edit():
    """Open user configuration file in $EDITOR."""
    user_config_file = config_manager._get_user_config_dir() / "config.yaml"

    # Create file if it doesn't exist
    if not user_config_file.exists():
        user_config_file.parent.mkdir(parents=True, exist_ok=True)
        config_manager.create_user_config_template()
        click.echo(f"📄 Created config template: {user_config_file}")

    editor = os.environ.get("EDITOR", "nano")
    try:
        subprocess.run([editor, str(user_config_file)], check=True)
        config_manager.reload_configs()
        click.echo("✅ Configuration reloaded")
    except subprocess.CalledProcessError:
        click.echo(f"❌ Error opening editor: {editor}")
        raise click.Abort()
    except FileNotFoundError:
        click.echo(f"❌ Editor not found: {editor}")
        click.echo("Set the EDITOR environment variable or install a text editor")
        raise click.Abort()


@config_cli.command()
def validate():
    """Validate configuration files."""
    files = config_manager.get_config_files()
    valid = True

    for name, config_file in files.items():
        if not config_file.exists():
            continue

        try:
            OmegaConf.load(config_file)
            click.echo(f"✅ {name}: {config_file}")
        except Exception as e:
            click.echo(f"❌ {name}: {config_file} - {e}")
            valid = False

    if valid:
        click.echo("🎉 All configuration files are valid")
    else:
        click.echo("⚠️  Some configuration files have errors")
        raise click.Abort()


@config_cli.command()
def reset():
    """Reset user configuration to defaults."""
    user_config_file = config_manager._get_user_config_dir() / "config.yaml"

    if user_config_file.exists():
        if click.confirm(f"Delete user config file {user_config_file}?"):
            user_config_file.unlink()
            config_manager.reload_configs()
            click.echo("✅ User configuration reset")
        else:
            click.echo("Reset cancelled")
    else:
        click.echo("No user configuration file to reset")


@config_cli.command(name="show-files")
def show_files():
    """Show configuration file locations and status."""
    files = config_manager.get_config_files()

    table = Table(title="Configuration Files")
    table.add_column("Type", style="cyan")
    table.add_column("Path", style="yellow")
    table.add_column("Status", style="green")

    for name, config_file in files.items():
        if config_file.exists():
            status = "✅ Found"
            # Check if it's readable
            try:
                OmegaConf.load(config_file)
            except Exception:
                status = "❌ Invalid"
        else:
            status = "➖ Not found"

        table.add_row(name.replace("_", " ").title(), str(config_file), status)

    console.print(table)


@config_cli.command()
def init():
    """Create user configuration file template."""
    try:
        config_file = config_manager.create_user_config_template()
        click.echo(f"✅ Created configuration template: {config_file}")
        click.echo("Edit the file to customize your defaults and aliases")
    except click.ClickException:
        click.echo("❌ User configuration file already exists")
        click.echo(
            "Use 'qxub config edit' to modify or 'qxub config reset' to recreate"
        )


# Alias management subcommands
@config_cli.group(name="alias")
def alias_config():
    """Manage configuration aliases."""
    pass


@alias_config.command()
def list_aliases():
    """List all available aliases."""
    aliases = config_manager.list_aliases()

    if not aliases:
        click.echo("No aliases defined")
        return

    table = Table(title="Available Aliases")
    table.add_column("Name", style="cyan")
    table.add_column("Subcommand", style="yellow")
    table.add_column("Command", style="green")
    table.add_column("Description", style="white")

    for alias_name in sorted(aliases):
        alias_def = config_manager.get_alias(alias_name)

        # Handle both new hierarchical and legacy flat formats
        if "subcommand" in alias_def and isinstance(alias_def["subcommand"], dict):
            # New hierarchical format
            subcommand_def = alias_def.get("subcommand", {})
            subcommand = subcommand_def.get("type", "—")
            target_def = alias_def.get("target", {})
            cmd = target_def.get("cmd", "—")
            main_def = alias_def.get("main", {})
            job_name = main_def.get("name", alias_def.get("name", ""))
        else:
            # Legacy flat format
            subcommand = alias_def.get("subcommand", "—")
            cmd = alias_def.get("cmd", "—")
            job_name = alias_def.get("name", "")

        # Truncate long commands
        if len(cmd) > 50:
            cmd = cmd[:47] + "..."

        # Look for description or create one
        description = alias_def.get("description", "")
        if not description and job_name:
            description = f"Job: {job_name}"

        table.add_row(alias_name, subcommand, cmd, description)

    console.print(table)


@alias_config.command()
@click.argument("alias_name")
def show(alias_name: str):
    """Show detailed information about a specific alias."""
    alias_def = config_manager.get_alias(alias_name)

    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        raise click.Abort()

    click.echo(f"📋 Alias: {alias_name}")
    yaml_str = OmegaConf.to_yaml(OmegaConf.create(alias_def))
    syntax = Syntax(yaml_str, "yaml", theme="monokai", line_numbers=True)
    console.print(syntax)


@alias_config.command()
@click.argument("alias_name")
def test(alias_name: str):  # pylint: disable=too-many-branches
    """Test and validate an alias definition."""
    alias_def = config_manager.get_alias(alias_name)

    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        raise click.Abort()

    click.echo(f"🧪 Testing alias: {alias_name}")

    # Validate required fields
    errors = []
    warnings = []

    # Handle both new hierarchical and legacy flat formats
    if "subcommand" in alias_def and isinstance(alias_def["subcommand"], dict):
        # New hierarchical format
        subcommand_def = alias_def.get("subcommand", {})
        subcommand_type = subcommand_def.get("type")
        target_def = alias_def.get("target", {})
        cmd = target_def.get("cmd")
    else:
        # Legacy flat format
        subcommand_type = alias_def.get("subcommand")
        cmd = alias_def.get("cmd")

    if not subcommand_type:
        errors.append("Missing required field: subcommand.type")
    elif subcommand_type not in ["conda", "module", "sing"]:
        errors.append(f"Invalid subcommand type: {subcommand_type}")

    if not cmd:
        warnings.append(
            "No 'cmd' specified - alias can only be used with command override"
        )

    # Check subcommand-specific requirements
    if subcommand_type == "sing":
        if "sif" not in subcommand_def and "sif" not in alias_def:
            sing_config = alias_def[subcommand]
            if not sing_config.get("sif"):
                warnings.append("Singularity container (sif) not specified")

    # Test template resolution
    try:
        template_vars = config_manager.get_template_variables(
            name=alias_def.get("name", "test"),
            project=alias_def.get("project", "test"),
            queue=alias_def.get("queue", "test"),
        )
        resolved = config_manager.resolve_templates(alias_def, template_vars)
        click.echo("✅ Template resolution successful")
    except Exception as e:
        errors.append(f"Template resolution failed: {e}")

    # Report results
    if errors:
        click.echo("❌ Errors found:")
        for error in errors:
            click.echo(f"  • {error}")

    if warnings:
        click.echo("⚠️  Warnings:")
        for warning in warnings:
            click.echo(f"  • {warning}")

    if not errors and not warnings:
        click.echo("✅ Alias definition is valid")
    elif not errors:
        click.echo("✅ Alias definition is valid (with warnings)")
    else:
        raise click.Abort()


@alias_config.command()
@click.argument("alias_name")
@click.option(
    "--subcommand",
    type=click.Choice(["conda", "module", "sing"]),
    help="Which qxub subcommand to use",
)
@click.option("--cmd", help="Command to execute")
@click.option("--name", help="Job name")
@click.option("--queue", help="Queue name")
@click.option(
    "--resources", multiple=True, help="Resource requests (can specify multiple)"
)
@click.option("--env", help="Conda environment (for conda subcommand)")
@click.option("--mod", help="Single module to load (for module subcommand)")
@click.option(
    "--mods", help="Multiple modules to load, comma-separated (for module subcommand)"
)
@click.option("--sif", help="Singularity container (for sing subcommand)")
def set_alias(alias_name: str, **kwargs):
    """Create or update an alias."""
    # Remove None values
    alias_def = {k: v for k, v in kwargs.items() if v is not None}

    if not alias_def:
        click.echo("❌ No options provided")
        raise click.Abort()

    # Validate mutually exclusive options
    if "mod" in alias_def and "mods" in alias_def:
        click.echo("❌ Cannot specify both --mod and --mods")
        raise click.Abort()

    # Convert multiple values to lists
    if "resources" in alias_def and alias_def["resources"]:
        alias_def["resources"] = list(alias_def["resources"])

    # Handle module options
    if "mod" in alias_def:
        # Single module - convert to list
        alias_def["mod"] = [alias_def["mod"]]
    elif "mods" in alias_def:
        # Multiple modules - split by comma and clean up
        modules = [m.strip() for m in alias_def["mods"].split(",") if m.strip()]
        alias_def["mod"] = modules
        del alias_def["mods"]  # Remove the mods key, use mod for storage

    # Organize subcommand-specific options
    subcommand = alias_def.get("subcommand")
    if subcommand:
        subcommand_opts = {}
        if subcommand == "conda" and "env" in alias_def:
            subcommand_opts["env"] = alias_def.pop("env")
        elif subcommand == "module" and "mod" in alias_def:
            subcommand_opts["mod"] = alias_def.pop("mod")
        elif subcommand == "sing" and "sif" in alias_def:
            subcommand_opts["sif"] = alias_def.pop("sif")

        if subcommand_opts:
            alias_def[subcommand] = subcommand_opts

    try:
        # Set the alias in user config
        config_manager.set_user_config_value(f"aliases.{alias_name}", alias_def)
        click.echo(f"✅ Set alias '{alias_name}'")

        # Show the created alias
        click.echo("📋 Alias definition:")
        yaml_str = OmegaConf.to_yaml(OmegaConf.create(alias_def))
        syntax = Syntax(yaml_str, "yaml", theme="monokai")
        console.print(syntax)

    except Exception as e:
        click.echo(f"❌ Error creating alias: {e}")
        raise click.Abort()


@alias_config.command()
@click.argument("alias_name")
def delete(alias_name: str):
    """Delete an alias."""
    if not config_manager.get_alias(alias_name):
        click.echo(f"❌ Alias '{alias_name}' not found")
        raise click.Abort()

    if click.confirm(f"Delete alias '{alias_name}'?"):
        try:
            # Load user config
            user_config_file = config_manager._get_user_config_dir() / "config.yaml"
            if user_config_file.exists():
                user_config = OmegaConf.load(user_config_file)
                if "aliases" in user_config and alias_name in user_config.aliases:
                    del user_config.aliases[alias_name]
                    OmegaConf.save(user_config, user_config_file)
                    config_manager.reload_configs()
                    click.echo(f"✅ Deleted alias '{alias_name}'")
                else:
                    click.echo(f"❌ Alias '{alias_name}' not found in user config")
            else:
                click.echo("❌ No user config file found")
        except Exception as e:
            click.echo(f"❌ Error deleting alias: {e}")
            raise click.Abort()
    else:
        click.echo("Delete cancelled")
