"""
Configuration management CLI for qxub.

Provides commands for viewing, editing, and managing qxub configuration
including defaults and aliases.
"""

# pylint: disable=broad-exception-caught,protected-access,raise-missing-from,unnecessary-pass,unused-variable
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

import click
from omegaconf import OmegaConf
from rich.console import Console
from rich.syntax import Syntax
from rich.table import Table

from .config import config_manager
from .history import history_logger

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


@config_cli.command(name="set")
@click.argument("key_path")
@click.argument("value")
@click.option(
    "--global", "global_config", is_flag=True, help="Set value in user config (default)"
)
@click.option(
    "--system",
    is_flag=True,
    help="Set value in system config (requires appropriate permissions)",
)
@click.option(
    "--project", is_flag=True, help="Set value in project config (.qx/project.yaml)"
)
@click.option(
    "--local", is_flag=True, help="Set value in local config (.qx/local.yaml)"
)
@click.option("--test", is_flag=True, help="Set value in test config (.qx/test.yaml)")
def set_config(
    key_path: str,
    value: str,
    global_config: bool,
    system: bool,
    project: bool,
    local: bool,
    test: bool,
):
    """Set configuration value by key path (e.g., 'defaults.name' 'myjob')."""
    try:
        # Check for mutually exclusive options
        config_targets = [global_config, system, project, local, test]
        if sum(config_targets) > 1:
            click.echo(
                "Error: Cannot specify multiple config targets (--global, --system, --project, --local, --test)"
            )
            raise click.Abort()

        # Try to parse as YAML to handle different types
        try:
            # Create a temporary YAML document to parse the value
            yaml_doc = f"temp: {value}"
            parsed_config = OmegaConf.create(yaml_doc)
            parsed_value = parsed_config.temp
        except Exception:
            # If parsing fails, treat as string
            parsed_value = value

        if system:
            config_manager.set_system_config_value(key_path, parsed_value)
            click.echo(f"✅ Set {key_path} = {parsed_value} (system-wide)")
        elif global_config:
            config_manager.set_user_config_value(key_path, parsed_value)
            click.echo(f"✅ Set {key_path} = {parsed_value} (global/user config)")
        elif project:
            config_manager.set_project_config_value(key_path, parsed_value)
            click.echo(f"✅ Set {key_path} = {parsed_value} (project config)")
        elif local:
            config_manager.set_local_config_value(key_path, parsed_value)
            click.echo(f"✅ Set {key_path} = {parsed_value} (local config)")
        elif test:
            config_manager.set_test_config_value(key_path, parsed_value)
            click.echo(f"✅ Set {key_path} = {parsed_value} (test config)")
        else:
            # Default to user config when no specific target is specified
            config_manager.set_user_config_value(key_path, parsed_value)
            click.echo(f"✅ Set {key_path} = {parsed_value} (user config)")
    except Exception as e:
        click.echo(f"❌ Error setting configuration: {e}")
        raise click.Abort()


@config_cli.command(name="list")
@click.argument("section", required=False)
@click.option("--user-only", is_flag=True, help="Show only user configuration")
@click.option("--system-only", is_flag=True, help="Show only system configuration")
@click.option(
    "--show-origin", is_flag=True, help="Show source file path for each setting"
)
def list_config(
    section: Optional[str], user_only: bool, system_only: bool, show_origin: bool
):
    """List configuration values, optionally filtered by section."""
    if user_only and system_only:
        click.echo("Error: Cannot specify both --user-only and --system-only")
        raise click.Abort()

    if show_origin and (user_only or system_only):
        click.echo(
            "Error: --show-origin cannot be used with --user-only or --system-only"
        )
        raise click.Abort()

    if show_origin:
        _show_config_with_origin(section)
        return

    if user_only:
        config_to_show = config_manager.user_config
    elif system_only:
        config_to_show = config_manager.system_config
    else:
        config_to_show = config_manager.merged_config

    if not config_to_show:
        if user_only:
            click.echo("No user configuration found")
        elif system_only:
            click.echo("No system configuration found")
        else:
            click.echo("No configuration found")
        return

    if section:
        if user_only and config_manager.user_config:
            section_config = config_manager.user_config.get(section)
        elif system_only and config_manager.system_config:
            section_config = config_manager.system_config.get(section)
        else:
            section_config = config_manager.get_config_value(section)

        if section_config is None:
            click.echo(f"Section '{section}' not found")
            raise click.Abort()

        if isinstance(section_config, dict):
            config_dict = {section: section_config}
        else:
            config_dict = {section: section_config}
    else:
        config_dict = OmegaConf.to_container(config_to_show, resolve=True)

    yaml_str = OmegaConf.to_yaml(OmegaConf.create(config_dict))
    syntax = Syntax(yaml_str, "yaml", theme="monokai", line_numbers=True)
    console.print(syntax)


@config_cli.command()
def files():
    """List all configuration files and their locations."""
    files = config_manager.get_config_files()

    table = Table(title="Configuration Files")
    table.add_column("Type", style="cyan", no_wrap=True)
    table.add_column("File Path", style="green")
    table.add_column("Status", style="yellow")

    for name, config_file in files.items():
        if config_file.exists():
            status = "✅ Exists"
        else:
            status = "❌ Missing"

        table.add_row(name, str(config_file), status)

    console.print(table)


def _show_config_with_origin(section: Optional[str]):
    """Show configuration with origin information."""
    from rich.tree import Tree

    # Get config file paths
    config_files = config_manager.get_config_files()

    # Build merged config with origin tracking
    merged_with_origin = {}

    # Process configs in order of precedence (system -> user -> project -> local -> test)
    config_order = [
        ("system", config_manager.system_config),
        ("user", config_manager.user_config),
        ("project", config_manager.project_config),
        ("local", config_manager.local_config),
        ("test", config_manager.test_config),
    ]

    for config_type, config_obj in config_order:
        if config_obj:
            # Find matching config file
            config_file = None
            for name, path in config_files.items():
                if (
                    config_type == "system"
                    and name.startswith("system")
                    and path.exists()
                ) or (config_type == name and path.exists()):
                    config_file = path
                    break

            if config_file:
                config_dict = OmegaConf.to_container(config_obj, resolve=True)
                _add_origin_to_dict(merged_with_origin, config_dict, str(config_file))

    if not merged_with_origin:
        click.echo("No configuration found")
        return

    # Filter by section if specified
    if section:
        if section in merged_with_origin:
            display_dict = {section: merged_with_origin[section]}
        else:
            click.echo(f"Section '{section}' not found")
            raise click.Abort()
    else:
        display_dict = merged_with_origin

    # Display with Rich Tree
    tree = Tree("📋 Configuration (with origin)")
    _build_config_tree(tree, display_dict, "")
    console.print(tree)


def _add_origin_to_dict(target_dict, source_dict, origin_file):
    """Add values from source_dict to target_dict with origin tracking."""
    for key, value in source_dict.items():
        if isinstance(value, dict):
            if key not in target_dict:
                target_dict[key] = {}
            _add_origin_to_dict(target_dict[key], value, origin_file)
        else:
            # Store value with origin info
            target_dict[key] = {"_value": value, "_origin": origin_file}


def _build_config_tree(tree, config_dict, prefix):
    """Recursively build Rich tree from config dictionary with origin info."""
    for key, value in config_dict.items():
        if isinstance(value, dict) and "_value" in value and "_origin" in value:
            # Leaf node with origin info - show full path
            origin_path = value["_origin"]
            tree.add(
                f"[cyan]{key}[/cyan]: [green]{value['_value']}[/green] [dim]({origin_path})[/dim]"
            )
        elif isinstance(value, dict):
            # Branch node
            branch = tree.add(f"[cyan]{key}[/cyan]:")
            _build_config_tree(branch, value, f"{prefix}{key}.")


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


@config_cli.command(name="init-project")
@click.option(
    "--project-root",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Project root directory (defaults to current directory)",
)
def init_project(project_root: Optional[Path]):
    """Initialize project-level configuration directory (.qx)."""
    if project_root is None:
        project_root = Path.cwd()

    try:
        success = config_manager.init_project_config(project_root)
        if success:
            qx_dir = project_root / ".qx"
            console.print(
                f"✅ Initialized project configuration in {qx_dir}", style="green"
            )
            console.print("\nCreated files:", style="bold")
            console.print("  📁 .qx/")
            console.print("  ├── 📄 project.yaml   (git-tracked team settings)")
            console.print("  ├── 📄 test.yaml      (git-tracked CI settings)")
            console.print("  └── 🚫 .gitignore     (excludes local.yaml)")
            console.print(
                "\n💡 Create local.yaml for personal overrides (git-ignored)",
                style="dim",
            )
        else:
            console.print(
                "❌ Project already initialized (.qx directory exists)", style="red"
            )
    except Exception as e:
        console.print(f"❌ Failed to initialize project: {e}", style="red")


# Alias management subcommands
@config_cli.group(name="alias")
def alias_config():
    """Manage configuration aliases."""
    pass


@alias_config.command(name="list")
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
            sing_config = alias_def[subcommand_type]
            if not sing_config.get("sif"):
                warnings.append("Singularity container (sif) not specified")

    # Test template resolution
    try:
        template_vars = config_manager.get_template_variables(
            name=alias_def.get("name", "test"),
            project=alias_def.get("project", "test"),
            queue=alias_def.get("queue", "test"),
        )
        config_manager.resolve_templates(alias_def, template_vars)
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


@alias_config.command(name="set")
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
@click.option("--env", help="Conda environment to use")
@click.option("--mod", help="Single module to load (for module subcommand)")
@click.option(
    "--mods", help="Multiple modules to load, comma-separated (for module subcommand)"
)
@click.option("--sif", help="Singularity container (for sing subcommand)")
@click.option(
    "--user", "user_config", is_flag=True, help="Set alias in user config (default)"
)
@click.option(
    "--system",
    is_flag=True,
    help="Set alias in system config (requires appropriate permissions)",
)
def set_alias(alias_name: str, user_config: bool, system: bool, **kwargs):
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

    # Check for mutually exclusive scope options
    if user_config and system:
        raise click.ClickException(
            "Cannot specify both --user and --system. Choose one scope."
        )

    try:
        # Set the alias with appropriate scope
        if system:
            config_manager.set_system_config_value(f"aliases.{alias_name}", alias_def)
            scope_desc = "(system-wide)"
        else:
            # Default to user config
            config_manager.set_user_config_value(f"aliases.{alias_name}", alias_def)
            scope_desc = "(user config)"

        click.echo(f"✅ Set alias '{alias_name}' {scope_desc}")

        # Show the created alias
        click.echo("📋 Alias definition:")
        yaml_str = OmegaConf.to_yaml(OmegaConf.create(alias_def))
        syntax = Syntax(yaml_str, "yaml", theme="monokai")
        console.print(syntax)

    except (PermissionError, OSError) as e:
        if "Permission denied" in str(e) or "Read-only file system" in str(e):
            click.echo(f"❌ Permission denied writing to config file: {e}")
            click.echo("💡 Hint: Use 'sudo' for system-wide configuration")
        else:
            click.echo(f"❌ Error writing config file: {e}")
        raise click.Abort()
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


@config_cli.group(name="history")
def history_cli():
    """Manage qxub command history."""
    pass


@history_cli.command(name="list")
@click.option("--limit", "-n", default=10, help="Number of recent commands to show")
@click.option("--all", "show_all", is_flag=True, help="Show all available history")
def list_history(limit, show_all):
    """List recent qxub commands from history."""
    try:
        if show_all:
            limit = 1000  # Get maximum available

        commands = history_logger.get_recent_commands(limit)

        if not commands:
            console.print("📝 No command history found", style="yellow")
            return

        title_text = "All" if show_all else f"Last {min(limit, len(commands))}"
        table = Table(title=f"Recent qxub Commands ({title_text})")
        table.add_column("Time", style="cyan", no_wrap=True)
        table.add_column("Command", style="green")
        table.add_column("Type", style="yellow")
        table.add_column("Status", style="magenta")

        for cmd in commands:
            timestamp = cmd.get("timestamp", "Unknown")
            # Format timestamp to be more readable
            if timestamp != "Unknown":
                try:
                    from datetime import datetime

                    dt = datetime.fromisoformat(timestamp)
                    timestamp = dt.strftime("%Y-%m-%d %H:%M:%S")
                except:
                    pass

            command_line = cmd.get("command_line", "")
            # Truncate long commands
            if len(command_line) > 60:
                command_line = command_line[:57] + "..."

            cmd_type = ""
            if "subcommand" in cmd:
                cmd_type = cmd["subcommand"].get("type", "")
            elif "config" in command_line:
                cmd_type = "config"

            status = "✅ Success" if cmd.get("success", True) else "❌ Failed"

            table.add_row(timestamp, command_line, cmd_type, status)

        console.print(table)

    except Exception as e:
        console.print(f"❌ Error reading command history: {e}", style="red")


@history_cli.command(name="show")
@click.argument("index", type=int)
def show_history_item(index):
    """Show detailed information about a specific history item."""
    try:
        commands = history_logger.get_recent_commands(1000)  # Get all available

        if not commands:
            console.print("📝 No command history found", style="yellow")
            return

        # Convert to 1-based indexing (most recent is 1)
        if index < 1 or index > len(commands):
            console.print(
                f"❌ Invalid index. Available range: 1-{len(commands)}", style="red"
            )
            return

        # Get command (index 1 = most recent, so we need to reverse)
        cmd = commands[-(index)]

        console.print(f"📋 Command History Item #{index}", style="bold blue")
        console.print("")

        # Show basic info
        console.print(f"⏰ Time: {cmd.get('timestamp', 'Unknown')}")
        console.print(f"📂 Directory: {cmd.get('working_directory', 'Unknown')}")
        console.print(f"✅ Success: {cmd.get('success', True)}")
        if cmd.get("error"):
            console.print(f"❌ Error: {cmd['error']}", style="red")
        console.print("")

        # Show command line
        console.print("💻 Command Line:", style="bold")
        console.print(cmd.get("command_line", ""), style="green")
        console.print("")

        # Show alias-like structure
        alias_structure = {}
        if "main" in cmd:
            alias_structure["main"] = cmd["main"]
        if "subcommand" in cmd:
            alias_structure["subcommand"] = cmd["subcommand"]
        if "target" in cmd:
            alias_structure["target"] = cmd["target"]

        if alias_structure:
            console.print("🏗️ Alias-like Structure:", style="bold")
            yaml_content = OmegaConf.to_yaml(OmegaConf.create(alias_structure))
            syntax = Syntax(
                yaml_content, "yaml", theme="github-dark", line_numbers=True
            )
            console.print(syntax)

    except Exception as e:
        console.print(f"❌ Error showing history item: {e}", style="red")


@history_cli.command(name="clear")
@click.confirmation_option(prompt="Are you sure you want to clear all command history?")
def clear_history():
    """Clear all command history."""
    try:
        history_logger.clear_history()
        console.print("🗑️ Command history cleared", style="green")
    except Exception as e:
        console.print(f"❌ Error clearing history: {e}", style="red")


@history_cli.command(name="to-alias")
@click.argument("index", type=int)
@click.argument("alias_name")
@click.option("--overwrite", is_flag=True, help="Overwrite existing alias")
def history_to_alias(index, alias_name, overwrite):
    """Create an alias from a history item."""
    try:
        commands = history_logger.get_recent_commands(1000)  # Get all available

        if not commands:
            console.print("📝 No command history found", style="yellow")
            return

        # Convert to 1-based indexing (most recent is 1)
        if index < 1 or index > len(commands):
            console.print(
                f"❌ Invalid index. Available range: 1-{len(commands)}", style="red"
            )
            return

        # Get command (index 1 = most recent, so we need to reverse)
        cmd = commands[-(index)]

        # Only allow creation from execution commands, not config commands
        if "subcommand" not in cmd or cmd["subcommand"].get("type") not in [
            "conda",
            "module",
            "sing",
        ]:
            console.print(
                "❌ Can only create aliases from execution commands (conda, module, sing)",
                style="red",
            )
            return

        # Check if alias already exists
        existing_aliases = config_manager.list_aliases()
        if alias_name in existing_aliases and not overwrite:
            console.print(
                f"❌ Alias '{alias_name}' already exists. Use --overwrite to replace it.",
                style="red",
            )
            return

        # Create alias structure
        alias_def = {}
        if "main" in cmd:
            alias_def["main"] = cmd["main"]
        if "subcommand" in cmd:
            alias_def["subcommand"] = cmd["subcommand"]
        if "target" in cmd and cmd["target"]:
            # Convert target list to cmd string for alias format
            alias_def["target"] = {"cmd": " ".join(cmd["target"])}

        # Save the alias
        config_manager.save_alias(alias_name, alias_def)

        action = (
            "Updated" if (alias_name in existing_aliases and overwrite) else "Created"
        )
        console.print(
            f"✅ {action} alias '{alias_name}' from history item #{index}",
            style="green",
        )

        # Show the created alias
        console.print("\n📋 Created alias structure:", style="bold")
        yaml_content = OmegaConf.to_yaml(OmegaConf.create(alias_def))
        syntax = Syntax(yaml_content, "yaml", theme="github-dark", line_numbers=True)
        console.print(syntax)

    except Exception as e:
        console.print(f"❌ Error creating alias from history: {e}", style="red")


# Shortcut management subcommands
@config_cli.group(name="shortcut")
def shortcut_config():
    """Manage configuration shortcuts."""
    pass


@shortcut_config.command(name="set")
@click.argument("command_prefix")
@click.option("--env", help="Conda environment for execution")
@click.option("--mod", multiple=True, help="Environment module to load (repeatable)")
@click.option("--mods", help="Comma-separated list of environment modules")
@click.option("--sif", help="Singularity container image")
@click.option("--bind", help="Singularity bind mounts")
@click.option("-l", "--resources", help="PBS resource requirements")
@click.option(
    "--mem",
    "--memory",
    help="Memory requirement (workflow-friendly). Examples: '4GB', '2000MB', '16g'",
)
@click.option(
    "--runtime",
    "--time",
    help="Runtime/walltime limit (workflow-friendly). Examples: '2h', '30m', '1h30m', '02:30:00'",
)
@click.option(
    "--cpus",
    "--threads",
    type=int,
    help="Number of CPU cores/threads (workflow-friendly)",
)
@click.option(
    "--disk",
    "--jobfs",
    help="Local disk/jobfs requirement (workflow-friendly). Examples: '10GB', '500MB'",
)
@click.option(
    "--volumes",
    "--storage",
    help="Storage volumes to mount (NCI format). Examples: 'gdata/a56', 'gdata/a56+gdata/px14'",
)
@click.option("-q", "--queue", help="PBS queue name")
@click.option("-N", "--name", help="PBS job name")
@click.option("-P", "--project", help="PBS project code")
@click.option("--template", help="Job script template")
@click.option("--pre", help="Command to run before main command")
@click.option("--post", help="Command to run after main command")
@click.option("--cmd", help="Default command to execute (can be overridden)")
@click.option(
    "--user", "user_config", is_flag=True, help="Set shortcut in user config (default)"
)
@click.option(
    "--system",
    is_flag=True,
    help="Set shortcut in system config (requires appropriate permissions)",
)
def set_shortcut(command_prefix: str, user_config: bool, system: bool, **options):
    """
    Create or update a shortcut.

    Creates a shortcut for the specified command or command prefix.
    The shortcut will store execution context (conda env, modules, etc.)
    and PBS settings that will be automatically applied when you run
    matching commands.

    Examples:
        # Create a Python shortcut with workflow-friendly options
        qxub config shortcut set "python" --env datascience --mem 8GB --time 2h

        # Create a Jupyter shortcut with traditional PBS format
        qxub config shortcut set "jupyter" --env jupyter -l "mem=16GB,walltime=4:00:00"

        # Create a DVC shortcut with multiple modules
        qxub config shortcut set "dvc" --mod python3 --mod dvc

        # Create a high-memory analysis shortcut
        qxub config shortcut set "bigmem" --env scipy --mem 128GB --cpus 48 --time 8h

        # Create a Singularity ML shortcut
        qxub config shortcut set "tensorflow" --sif /path/to/tf.sif --bind /data:/mnt
    """
    from typing import Any, Dict

    from .shortcut_manager import shortcut_manager

    # Build shortcut definition
    definition: Dict[str, Any] = {}

    # Execution context (mutually exclusive)
    execution_contexts = [options["env"], options["mod"], options["sif"]]
    active_contexts = sum(bool(x) for x in execution_contexts)

    if active_contexts > 1:
        raise click.ClickException(
            "Cannot specify multiple execution contexts. "
            "Use only one of: --env, --mod/--mods, --sif"
        )

    # Set execution context
    if options["env"]:
        definition["env"] = options["env"]
    elif options["mod"]:
        definition["mod"] = list(options["mod"])
    elif options["mods"]:
        definition["mods"] = options["mods"]
    elif options["sif"]:
        definition["sif"] = options["sif"]
        if options["bind"]:
            definition["bind"] = options["bind"]

    # PBS options
    resources_list = []

    # Handle workflow-friendly resource options - convert to PBS format
    if options.get("mem") or options.get("memory"):
        from qxub.resources import ResourceMapper

        mapper = ResourceMapper()
        mapper.add_memory(options.get("mem") or options.get("memory"))
        resources_list.extend(mapper.get_pbs_resources())

    if options.get("runtime") or options.get("time"):
        from qxub.resources import ResourceMapper

        mapper = ResourceMapper()
        mapper.add_runtime(options.get("runtime") or options.get("time"))
        resources_list.extend(mapper.get_pbs_resources())

    if options.get("cpus") or options.get("threads"):
        from qxub.resources import ResourceMapper

        mapper = ResourceMapper()
        mapper.add_cpus(options.get("cpus") or options.get("threads"))
        resources_list.extend(mapper.get_pbs_resources())

    if options.get("disk") or options.get("jobfs"):
        from qxub.resources import ResourceMapper

        mapper = ResourceMapper()
        mapper.add_disk(options.get("disk") or options.get("jobfs"))
        resources_list.extend(mapper.get_pbs_resources())

    if options.get("volumes") or options.get("storage"):
        from qxub.resources import ResourceMapper

        mapper = ResourceMapper()
        mapper.add_storage(options.get("volumes") or options.get("storage"))
        resources_list.extend(mapper.get_pbs_resources())

    # Handle traditional PBS resources format
    if options["resources"]:
        resources_list.append(options["resources"])

    # Store resources if any were provided
    if resources_list:
        definition["resources"] = resources_list

    if options["queue"]:
        definition["queue"] = options["queue"]
    if options["name"]:
        definition["name"] = options["name"]
    if options["project"]:
        definition["project"] = options["project"]

    # Job script options
    if options["template"]:
        definition["template"] = options["template"]
    if options["pre"]:
        definition["pre"] = options["pre"]
    if options["post"]:
        definition["post"] = options["post"]

    # Command and metadata
    if options["cmd"]:
        definition["cmd"] = options["cmd"]

    # Validate we have at least some configuration
    if not definition:
        raise click.ClickException(
            "Shortcut must specify at least one option (execution context, PBS settings, etc.)"
        )

    # Check for mutually exclusive scope options
    if user_config and system:
        raise click.ClickException(
            "Cannot specify both --user and --system. Choose one scope."
        )

    # Save shortcut with appropriate scope
    try:
        if system:
            shortcut_manager.add_system_shortcut(command_prefix, definition)
            scope_desc = "(system-wide)"
        else:
            # Default to user config
            shortcut_manager.add_shortcut(command_prefix, definition)
            scope_desc = "(user config)"
    except (PermissionError, OSError) as e:
        if "Permission denied" in str(e) or "Read-only file system" in str(e):
            raise click.ClickException(
                f"❌ Permission denied writing to system config: {e}\n"
                f"💡 Hint: Use 'sudo' for system-wide configuration"
            )
        else:
            raise click.ClickException(f"❌ Error writing system config: {e}")

    # Show confirmation
    def _get_execution_context_description(definition: dict) -> str:
        """Get human-readable description of execution context."""
        if definition.get("env"):
            return f"conda: {definition['env']}"
        elif definition.get("mod"):
            modules = definition["mod"]
            if isinstance(modules, list):
                return f"modules: {', '.join(modules)}"
            return f"module: {modules}"
        elif definition.get("mods"):
            return f"modules: {definition['mods']}"
        elif definition.get("sif"):
            sif_desc = f"singularity: {definition['sif']}"
            if definition.get("bind"):
                sif_desc += f" (bind: {definition['bind']})"
            return sif_desc
        else:
            return "default"

    context_desc = _get_execution_context_description(definition)
    click.echo(f"✅ Shortcut '{command_prefix}' saved successfully {scope_desc}")
    click.echo(f"   Context: {context_desc}")
    if definition.get("cmd"):
        click.echo(f"   Command: {definition['cmd']}")

    # Show usage example
    click.echo(f"\n💡 Usage: qxub exec -- {command_prefix} [additional args]")
    click.echo(f"💡 Or explicit: qxub exec --shortcut '{command_prefix}' -- [command]")


@shortcut_config.command()
@click.option(
    "--show-origin", is_flag=True, help="Show origin (user/system) for each shortcut"
)
def list(show_origin: bool):
    """List all available shortcuts."""
    from .shortcut_manager import shortcut_manager

    if show_origin:
        shortcuts_data = shortcut_manager.list_shortcuts_with_origin()
        if not shortcuts_data:
            click.echo("No shortcuts defined yet.")
            click.echo(
                "\n💡 Create one with: qxub config shortcut set <command_prefix> --env <environment>"
            )
            return

        # Create table with origin information
        table = Table(title="Available Shortcuts")
        table.add_column("Command Prefix", style="cyan", no_wrap=True)
        table.add_column("Context", style="green")
        table.add_column("Command", style="yellow")
        table.add_column("Origin", style="dim")
    else:
        shortcuts = shortcut_manager.list_shortcuts()
        if not shortcuts:
            click.echo("No shortcuts defined yet.")
            click.echo(
                "\n💡 Create one with: qxub config shortcut set <command_prefix> --env <environment>"
            )
            return

        # Create table without origin information
        table = Table(title="Available Shortcuts")
        table.add_column("Command Prefix", style="cyan", no_wrap=True)
        table.add_column("Context", style="green")
        table.add_column("Command", style="yellow")

    def _get_execution_context_description(definition: dict) -> str:
        """Get human-readable description of execution context."""
        if definition.get("env"):
            return f"conda: {definition['env']}"
        elif definition.get("mod"):
            modules = definition["mod"]
            if isinstance(modules, list):
                return f"modules: {', '.join(modules)}"
            return f"module: {modules}"
        elif definition.get("mods"):
            return f"modules: {definition['mods']}"
        elif definition.get("sif"):
            sif_desc = f"singularity: {definition['sif']}"
            if definition.get("bind"):
                sif_desc += f" (bind: {definition['bind']})"
            return sif_desc
        else:
            return "default"

    if show_origin:
        for name, data in sorted(shortcuts_data.items()):
            definition = data["definition"]
            origin = data["origin"]
            context = _get_execution_context_description(definition)
            cmd = definition.get("cmd", "(dynamic)")
            table.add_row(name, context, cmd, origin)
    else:
        for name, definition in sorted(shortcuts.items()):
            context = _get_execution_context_description(definition)
            cmd = definition.get("cmd", "(dynamic)")
            table.add_row(name, context, cmd)

    console.print(table)

    click.echo(
        f"\n💡 Usage: qxub exec -- <command_prefix> [args] (automatic detection)"
    )
    click.echo(f"💡 Or explicit: qxub exec --shortcut <command_prefix> -- [command]")


@shortcut_config.command()
@click.argument("name")
def show(name: str):
    """Show detailed information about a specific shortcut."""
    import json

    from rich.syntax import Syntax

    from .shortcut_manager import shortcut_manager

    shortcut_def = shortcut_manager.get_shortcut(name)

    if not shortcut_def:
        click.echo(f"❌ Shortcut '{name}' not found")

        # Show similar shortcuts
        shortcuts = shortcut_manager.list_shortcuts()
        if shortcuts:
            click.echo("\n💡 Available shortcuts:")
            for shortcut_name in sorted(shortcuts.keys()):
                click.echo(f"  • {shortcut_name}")
        sys.exit(2)

    # Get source information
    shortcuts_with_origin = shortcut_manager.list_shortcuts_with_origin()
    source = shortcuts_with_origin.get(name, {}).get("origin", "unknown")
    source_display = "🌐 system config" if source == "system" else "👤 user config"

    # Display shortcut details
    click.echo(f"🎯 Shortcut: {name}")
    click.echo(f"   Source: {source_display}")

    def _get_execution_context_description(definition: dict) -> str:
        """Get human-readable description of execution context."""
        if definition.get("env"):
            return f"conda: {definition['env']}"
        elif definition.get("mod"):
            modules = definition["mod"]
            if isinstance(modules, list):
                return f"modules: {', '.join(modules)}"
            return f"module: {modules}"
        elif definition.get("mods"):
            return f"modules: {definition['mods']}"
        elif definition.get("sif"):
            sif_desc = f"singularity: {definition['sif']}"
            if definition.get("bind"):
                sif_desc += f" (bind: {definition['bind']})"
            return sif_desc
        else:
            return "default"

    # Display shortcut details
    click.echo(f"🎯 Shortcut: {name}")

    # Execution context
    context = _get_execution_context_description(shortcut_def)
    click.echo(f"   Context: {context}")

    # PBS settings
    pbs_settings = []
    if shortcut_def.get("queue"):
        pbs_settings.append(f"queue={shortcut_def['queue']}")
    if shortcut_def.get("resources"):
        pbs_settings.append(f"resources={shortcut_def['resources']}")
    if shortcut_def.get("name"):
        pbs_settings.append(f"job_name={shortcut_def['name']}")
    if shortcut_def.get("project"):
        pbs_settings.append(f"project={shortcut_def['project']}")

    if pbs_settings:
        click.echo(f"   PBS: {', '.join(pbs_settings)}")

    # Commands
    if shortcut_def.get("cmd"):
        click.echo(f"   Command: {shortcut_def['cmd']}")
    if shortcut_def.get("pre"):
        click.echo(f"   Pre-command: {shortcut_def['pre']}")
    if shortcut_def.get("post"):
        click.echo(f"   Post-command: {shortcut_def['post']}")

    # Description
    if shortcut_def.get("description"):
        click.echo(f"   Description: {shortcut_def['description']}")

    # Show raw JSON for advanced users
    if click.confirm("\n🔍 Show raw configuration?", default=False):
        json_content = json.dumps(shortcut_def, indent=2, sort_keys=True)
        syntax = Syntax(json_content, "json", theme="monokai", line_numbers=False)
        console.print(syntax)

    # Usage example
    click.echo(f"\n💡 Usage: qxub exec -- {name} [args] (automatic detection)")
    click.echo(f"💡 Or explicit: qxub exec --shortcut '{name}' -- [command]")


@shortcut_config.command()
@click.argument("name")
@click.option("--user", is_flag=True, help="Delete from user config (default)")
@click.option(
    "--system",
    is_flag=True,
    help="Delete from system config (requires write permissions)",
)
@click.option("--yes", is_flag=True, help="Confirm deletion without prompting")
def delete(name: str, user: bool, system: bool, yes: bool):
    """Delete a shortcut."""
    from .shortcut_manager import shortcut_manager

    # Validate mutually exclusive scope flags
    if user and system:
        raise click.ClickException(
            "Cannot specify both --user and --system flags. Choose one or omit both for default (user)."
        )

    # Manual confirmation if --yes not provided
    if not yes:
        if not click.confirm("Are you sure you want to delete this shortcut?"):
            click.echo("❌ Deletion cancelled")
            return

    try:
        if system:
            success = shortcut_manager.remove_system_shortcut(name)
            location = "system config"
        else:
            # Default to user if neither flag specified, or if --user explicitly set
            success = shortcut_manager.remove_shortcut(name)
            location = "user config"

        if success:
            click.echo(f"✅ Shortcut '{name}' deleted successfully from {location}")
        else:
            click.echo(f"❌ Shortcut '{name}' not found in {location}")

            # Show available shortcuts
            shortcuts = shortcut_manager.list_shortcuts()
            if shortcuts:
                click.echo("\n💡 Available shortcuts:")
                for shortcut_name in sorted(shortcuts.keys()):
                    click.echo(f"  • {shortcut_name}")

    except PermissionError:
        raise click.ClickException(
            f"Permission denied: Cannot write to system shortcuts file. "
            f"Check write permissions for {shortcut_manager._system_shortcuts_file}"
        )


@shortcut_config.command()
@click.argument("old_name")
@click.argument("new_name")
@click.option("--user", is_flag=True, help="Rename in user config (default)")
@click.option(
    "--system",
    is_flag=True,
    help="Rename in system config (requires write permissions)",
)
def rename(old_name: str, new_name: str, user: bool, system: bool):
    """Rename a shortcut."""
    from .shortcut_manager import shortcut_manager

    # Validate mutually exclusive scope flags
    if user and system:
        raise click.ClickException(
            "Cannot specify both --user and --system flags. Choose one or omit both for default (user)."
        )

    # Get existing shortcut - need to check the appropriate config
    shortcuts_with_origin = shortcut_manager.list_shortcuts_with_origin()
    shortcut_info = shortcuts_with_origin.get(old_name)

    if not shortcut_info:
        click.echo(f"❌ Shortcut '{old_name}' not found")
        return

    shortcut_def = shortcut_info["definition"]
    shortcut_source = shortcut_info["origin"]

    # Determine target scope
    if system:
        target_scope = "system"
    else:
        target_scope = "user"  # Default behavior

    # Warn if trying to rename across scopes
    if shortcut_source != target_scope:
        click.echo(
            f"⚠️  Warning: '{old_name}' exists in {shortcut_source} config but you're targeting {target_scope} config"
        )
        if not click.confirm(f"Continue with rename in {target_scope} config?"):
            click.echo("❌ Rename cancelled")
            return

    # Check if new name already exists in target scope
    existing_new = shortcut_manager.get_shortcut(new_name)
    if existing_new:
        new_info = shortcuts_with_origin.get(new_name)
        if new_info and new_info["origin"] == target_scope:
            if not click.confirm(
                f"Shortcut '{new_name}' already exists in {target_scope} config. Overwrite?"
            ):
                click.echo("❌ Rename cancelled")
                return

    try:
        # Add with new name in target scope
        if target_scope == "system":
            shortcut_manager.add_system_shortcut(new_name, shortcut_def)
        else:
            shortcut_manager.add_shortcut(new_name, shortcut_def)

        # Remove old name from original scope
        if shortcut_source == "system":
            shortcut_manager.remove_system_shortcut(old_name)
        else:
            shortcut_manager.remove_shortcut(old_name)

        click.echo(
            f"✅ Shortcut renamed from '{old_name}' to '{new_name}' in {target_scope} config"
        )

    except PermissionError:
        raise click.ClickException(
            f"Permission denied: Cannot write to system shortcuts file. "
            f"Check write permissions for {shortcut_manager._system_shortcuts_file}"
        )


@shortcut_config.command()
def files():
    """Show shortcuts configuration file locations."""
    from .shortcut_manager import shortcut_manager

    config_files = shortcut_manager.get_config_files()

    table = Table(title="Shortcuts Configuration Files")
    table.add_column("Type", style="cyan")
    table.add_column("Path", style="green")
    table.add_column("Exists", style="yellow")

    for file_type, (path, exists) in config_files.items():
        exists_str = "✅ Yes" if exists else "❌ No"
        table.add_row(file_type.title(), str(path), exists_str)

    console.print(table)


@shortcut_config.command()
def refresh():
    """Refresh shortcuts cache (reload from files)."""
    from .shortcut_manager import shortcut_manager

    shortcut_manager.refresh_cache()
    shortcuts = shortcut_manager.list_shortcuts()
    count = len(shortcuts)
    click.echo(f"✅ Shortcuts cache refreshed ({count} shortcuts loaded)")
