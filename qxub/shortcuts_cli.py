"""
Shortcuts management CLI commands for qxub.

Provides commands for creating, listing, editing, and managing shortcuts
for fast command execution with pre-configured environments.
"""

import json
from pathlib import Path
from typing import Any, Dict

import click
from rich.console import Console
from rich.syntax import Syntax
from rich.table import Table

from .shortcut_manager import shortcut_manager

console = Console()


@click.group(name="shortcut")
def shortcuts_cli():
    """Manage qxub shortcuts for fast command execution."""
    pass


@shortcuts_cli.command()
@click.argument("name")
@click.option("--env", help="Conda environment for execution")
@click.option("--mod", multiple=True, help="Environment module to load (repeatable)")
@click.option("--mods", help="Comma-separated list of environment modules")
@click.option("--sif", help="Singularity container image")
@click.option("--bind", help="Singularity bind mounts")
@click.option("-l", "--resources", multiple=True, help="PBS resource requirements")
@click.option("-q", "--queue", help="PBS queue name")
@click.option("-N", "--name", "job_name", help="PBS job name")
@click.option("-P", "--project", help="PBS project code")
@click.option("--template", help="Job script template")
@click.option("--pre", help="Command to run before main command")
@click.option("--post", help="Command to run after main command")
@click.option("--cmd", help="Default command to execute (can be overridden)")
@click.option("--description", help="Human-readable description of the shortcut")
def set(name: str, **options):
    """
    Create or update a shortcut.

    Creates a shortcut with the specified name and configuration. The shortcut
    will store execution context (conda env, modules, etc.) and PBS settings
    for quick reuse.

    Examples:
        # Create a Python data science shortcut
        qxub shortcut set "python-ds" --env datascience --queue express

        # Create a Jupyter shortcut with resource requirements
        qxub shortcut set "jupyter" --env jupyter -l "mem=16GB,walltime=4:00:00"

        # Create a DVC shortcut with modules
        qxub shortcut set "dvc doctor" --mod python3 --mod dvc --cmd "dvc doctor"

        # Create a Singularity ML shortcut
        qxub shortcut set "tensorflow" --sif /path/to/tf.sif --bind /data:/mnt
    """
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
    if options["resources"]:
        definition["resources"] = list(options["resources"])
    if options["queue"]:
        definition["queue"] = options["queue"]
    if options["job_name"]:
        definition["name"] = options["job_name"]
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
    if options["description"]:
        definition["description"] = options["description"]

    # Validate we have at least some configuration
    if not definition:
        raise click.ClickException(
            "Shortcut must specify at least one option (execution context, PBS settings, etc.)"
        )

    # Save shortcut
    shortcut_manager.add_shortcut(name, definition)

    # Show confirmation
    context_desc = _get_execution_context_description(definition)
    click.echo(f"✅ Shortcut '{name}' saved successfully")
    click.echo(f"   Context: {context_desc}")
    if definition.get("cmd"):
        click.echo(f"   Command: {definition['cmd']}")

    # Show usage example
    click.echo(f"\n💡 Usage: qxub exec --shortcut '{name}' -- [command]")


@shortcuts_cli.command()
def list():
    """List all available shortcuts."""
    shortcuts = shortcut_manager.list_shortcuts()

    if not shortcuts:
        click.echo("No shortcuts defined yet.")
        click.echo("\n💡 Create one with: qxub shortcut set <name> --env <environment>")
        return

    # Create table
    table = Table(title="Available Shortcuts")
    table.add_column("Name", style="cyan", no_wrap=True)
    table.add_column("Context", style="green")
    table.add_column("Command", style="yellow")
    table.add_column("Description", style="dim")

    for name, definition in sorted(shortcuts.items()):
        context = _get_execution_context_description(definition)
        cmd = definition.get("cmd", "(dynamic)")
        description = definition.get("description", "")
        table.add_row(name, context, cmd, description)

    console.print(table)

    click.echo(
        f"\n💡 Usage: qxub exec --shortcut <name> -- [command] or qxub exec -- <name> [args]"
    )


@shortcuts_cli.command()
@click.argument("name")
def show(name: str):
    """Show detailed information about a specific shortcut."""
    shortcut_def = shortcut_manager.get_shortcut(name)

    if not shortcut_def:
        click.echo(f"❌ Shortcut '{name}' not found")

        # Show similar shortcuts
        shortcuts = shortcut_manager.list_shortcuts()
        if shortcuts:
            click.echo("\n💡 Available shortcuts:")
            for shortcut_name in sorted(shortcuts.keys()):
                click.echo(f"  • {shortcut_name}")
        return

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
    click.echo(
        f"\n💡 Usage: qxub exec --shortcut '{name}' -- [command] or qxub exec -- {name} [args]"
    )


@shortcuts_cli.command()
@click.argument("name")
@click.confirmation_option(prompt="Are you sure you want to delete this shortcut?")
def delete(name: str):
    """Delete a shortcut."""
    if shortcut_manager.remove_shortcut(name):
        click.echo(f"✅ Shortcut '{name}' deleted successfully")
    else:
        click.echo(f"❌ Shortcut '{name}' not found")

        # Show available shortcuts
        shortcuts = shortcut_manager.list_shortcuts()
        if shortcuts:
            click.echo("\n💡 Available shortcuts:")
            for shortcut_name in sorted(shortcuts.keys()):
                click.echo(f"  • {shortcut_name}")


@shortcuts_cli.command()
@click.argument("old_name")
@click.argument("new_name")
def rename(old_name: str, new_name: str):
    """Rename a shortcut."""
    # Get existing shortcut
    shortcut_def = shortcut_manager.get_shortcut(old_name)
    if not shortcut_def:
        click.echo(f"❌ Shortcut '{old_name}' not found")
        return

    # Check if new name already exists
    if shortcut_manager.get_shortcut(new_name):
        if not click.confirm(f"Shortcut '{new_name}' already exists. Overwrite?"):
            click.echo("❌ Rename cancelled")
            return

    # Add with new name and remove old
    shortcut_manager.add_shortcut(new_name, shortcut_def)
    shortcut_manager.remove_shortcut(old_name)

    click.echo(f"✅ Shortcut renamed from '{old_name}' to '{new_name}'")


@shortcuts_cli.command()
def files():
    """Show shortcuts configuration file locations."""
    config_files = shortcut_manager.get_config_files()

    table = Table(title="Shortcuts Configuration Files")
    table.add_column("Type", style="cyan")
    table.add_column("Path", style="green")
    table.add_column("Exists", style="yellow")

    for file_type, (path, exists) in config_files.items():
        exists_str = "✅ Yes" if exists else "❌ No"
        table.add_row(file_type.title(), str(path), exists_str)

    console.print(table)


@shortcuts_cli.command()
def refresh():
    """Refresh shortcuts cache (reload from files)."""
    shortcut_manager.refresh_cache()
    shortcuts = shortcut_manager.list_shortcuts()
    count = len(shortcuts)
    click.echo(f"✅ Shortcuts cache refreshed ({count} shortcuts loaded)")


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
