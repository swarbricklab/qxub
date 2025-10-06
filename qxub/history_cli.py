"""
Top-level history CLI commands for qxub.

Provides commands for managing the dual-log history system:
- recipes: Manage computational recipes
- executions: Manage execution records
- conversion utilities
"""

import click
from rich.console import Console
from rich.syntax import Syntax
from rich.table import Table
from omegaconf import OmegaConf

from .history_manager import history_manager

console = Console()


@click.group()
def history():
    """Manage qxub execution history and recipes."""
    pass


@history.command()
@click.option("--limit", default=10, help="Number of recipes to show")
def recipes(limit):
    """List all computational recipes."""
    recipes_dict = history_manager.get_recipes()

    if not recipes_dict:
        console.print("📝 No recipes found", style="yellow")
        return

    table = Table(title="Computational Recipes")
    table.add_column("Hash", style="cyan", no_wrap=True)
    table.add_column("Executor", style="magenta")
    table.add_column("Target", style="green")
    table.add_column("Runs", justify="right", style="blue")
    table.add_column("Last Used", style="dim")

    # Sort by run count (most used first), then by last_seen
    sorted_recipes = sorted(
        recipes_dict.items(),
        key=lambda x: (
            x[1].get("metadata", {}).get("run_count", 0),
            x[1].get("metadata", {}).get("last_seen", ""),
        ),
        reverse=True,
    )

    for recipe_hash, recipe in sorted_recipes[:limit]:
        executor_info = recipe.get("executor", {})
        executor_type = executor_info.get("type", "unknown")

        # Build executor description
        if executor_type == "conda":
            executor_desc = f"conda --env {executor_info.get('env', '?')}"
        elif executor_type == "module":
            mods = executor_info.get("mod") or executor_info.get("mods", "?")
            if isinstance(mods, list):
                mods = ",".join(mods)
            executor_desc = f"module --mod {mods}"
        elif executor_type == "sing":
            sif = executor_info.get("sif", "?")
            executor_desc = f"sing --sif {sif}"
        else:
            executor_desc = executor_type

        target_cmd = recipe.get("target", {}).get("cmd", "?")

        metadata = recipe.get("metadata", {})
        run_count = metadata.get("run_count", 0)
        last_seen = metadata.get("last_seen", "never")
        if last_seen != "never":
            # Format timestamp
            try:
                from datetime import datetime

                dt = datetime.fromisoformat(last_seen.replace("Z", "+00:00"))
                last_seen = dt.strftime("%Y-%m-%d %H:%M")
            except Exception:
                pass

        table.add_row(
            recipe_hash,
            executor_desc,
            target_cmd[:50] + "..." if len(target_cmd) > 50 else target_cmd,
            str(run_count),
            last_seen,
        )

    console.print(table)


@history.command()
@click.option("--limit", default=10, help="Number of executions to show")
def executions(limit):
    """List recent executions."""
    executions_list = history_manager.get_executions(limit)

    if not executions_list:
        console.print("📝 No executions found", style="yellow")
        return

    table = Table(title="Recent Executions")
    table.add_column("Time", style="cyan")
    table.add_column("Recipe", style="magenta", no_wrap=True)
    table.add_column("Command", style="green")
    table.add_column("Status", style="blue")
    table.add_column("Directory", style="dim")

    for execution in executions_list:
        # Format timestamp
        try:
            timestamp_micro = int(execution["timestamp"])
            timestamp_sec = timestamp_micro / 1000000
            from datetime import datetime

            dt = datetime.fromtimestamp(timestamp_sec)
            time_str = dt.strftime("%m-%d %H:%M:%S")
        except Exception:
            time_str = execution["timestamp"][:10]

        recipe_hash = execution.get("recipe_hash", "unknown")[:8]
        command_line = execution.get("context", {}).get("command_line", "")

        # Extract just the command part for display
        if "conda" in command_line:
            cmd_start = command_line.find("conda")
        elif "module" in command_line:
            cmd_start = command_line.find("module")
        elif "sing" in command_line:
            cmd_start = command_line.find("sing")
        else:
            cmd_start = 0

        display_cmd = command_line[cmd_start:] if cmd_start > 0 else command_line
        display_cmd = display_cmd[:60] + "..." if len(display_cmd) > 60 else display_cmd

        status = execution.get("execution", {}).get("status", "unknown")
        status_style = (
            "green" if status == "completed" else "red" if status == "failed" else "yellow"
        )

        working_dir = execution.get("context", {}).get("working_directory", "")
        working_dir = working_dir.replace(str(click.get_app_dir("qxub", force_posix=True)), "~")
        working_dir = working_dir[-30:] if len(working_dir) > 30 else working_dir

        table.add_row(
            time_str,
            recipe_hash,
            display_cmd,
            f"[{status_style}]{status}[/{status_style}]",
            working_dir,
        )

    console.print(table)


@history.command()
@click.argument("recipe_hash")
def show(recipe_hash):
    """Show detailed information about a recipe."""
    recipe = history_manager.get_recipe_by_hash(recipe_hash)

    if not recipe:
        console.print(f"❌ Recipe not found: {recipe_hash}", style="red")
        return

    console.print(f"📋 Recipe Details: {recipe_hash}", style="bold")
    console.print()

    # Show recipe structure
    console.print("🏗️ Recipe Structure:", style="bold")
    recipe_copy = recipe.copy()
    if "metadata" in recipe_copy:
        del recipe_copy["metadata"]  # Don't show metadata in structure

    yaml_content = OmegaConf.to_yaml(OmegaConf.create(recipe_copy))
    syntax = Syntax(yaml_content, "yaml", theme="github-dark", line_numbers=True)
    console.print(syntax)

    # Show metadata
    metadata = recipe.get("metadata", {})
    if metadata:
        console.print("\n📊 Usage Statistics:", style="bold")
        console.print(f"First seen: {metadata.get('first_seen', 'unknown')}")
        console.print(f"Last seen: {metadata.get('last_seen', 'unknown')}")
        console.print(f"Run count: {metadata.get('run_count', 0)}")


@history.command()
@click.argument("recipe_hash")
@click.option("--limit", default=5, help="Number of executions to show")
def runs(recipe_hash, limit):
    """Show executions for a specific recipe."""
    recipe = history_manager.get_recipe_by_hash(recipe_hash)
    if not recipe:
        console.print(f"❌ Recipe not found: {recipe_hash}", style="red")
        return

    executions_list = history_manager.get_executions_for_recipe(recipe_hash, limit)

    if not executions_list:
        console.print(f"📝 No executions found for recipe {recipe_hash}", style="yellow")
        return

    console.print(f"📋 Executions for Recipe: {recipe_hash}", style="bold")
    console.print()

    table = Table()
    table.add_column("Time", style="cyan")
    table.add_column("Status", style="blue")
    table.add_column("Directory", style="green")
    table.add_column("Job ID", style="magenta")

    for execution in executions_list:
        # Format timestamp
        try:
            timestamp_micro = int(execution["timestamp"])
            timestamp_sec = timestamp_micro / 1000000
            from datetime import datetime

            dt = datetime.fromtimestamp(timestamp_sec)
            time_str = dt.strftime("%Y-%m-%d %H:%M:%S")
        except Exception:
            time_str = execution["timestamp"]

        status = execution.get("execution", {}).get("status", "unknown")
        status_style = (
            "green" if status == "completed" else "red" if status == "failed" else "yellow"
        )

        working_dir = execution.get("context", {}).get("working_directory", "")
        job_id = execution.get("execution", {}).get("job_id") or "—"

        table.add_row(time_str, f"[{status_style}]{status}[/{status_style}]", working_dir, job_id)

    console.print(table)


@history.command()
@click.confirmation_option(
    prompt="Are you sure you want to clear all history? This cannot be undone."
)
def clear():
    """Clear all history (recipes and executions)."""
    try:
        history_manager.clear_all_history()
        console.print("✅ History cleared successfully", style="green")
    except Exception as e:
        console.print(f"❌ Error clearing history: {e}", style="red")


@history.command(name="to-alias")
@click.argument("recipe_hash")
@click.argument("alias_name")
@click.option("--overwrite", is_flag=True, help="Overwrite existing alias")
def recipe_to_alias(recipe_hash, alias_name, overwrite):
    """Convert a recipe to an alias."""
    from .config_manager import config_manager

    try:
        # Get recipe
        recipe = history_manager.get_recipe_by_hash(recipe_hash)
        if not recipe:
            console.print(f"❌ Recipe not found: {recipe_hash}", style="red")
            return

        # Check if alias already exists
        existing_aliases = config_manager.list_aliases()
        if alias_name in existing_aliases and not overwrite:
            console.print(
                f"❌ Alias '{alias_name}' already exists. Use --overwrite to replace it.",
                style="red",
            )
            return

        # Convert to alias format
        alias_def = history_manager.convert_recipe_to_alias(recipe_hash)
        if not alias_def:
            console.print(f"❌ Could not convert recipe {recipe_hash} to alias format", style="red")
            return

        # Save the alias
        config_manager.save_alias(alias_name, alias_def)

        action = "Updated" if (alias_name in existing_aliases and overwrite) else "Created"
        console.print(f"✅ {action} alias '{alias_name}' from recipe {recipe_hash}", style="green")

        # Show the created alias
        console.print("\n📋 Created alias structure:", style="bold")
        yaml_content = OmegaConf.to_yaml(OmegaConf.create(alias_def))
        syntax = Syntax(yaml_content, "yaml", theme="github-dark", line_numbers=True)
        console.print(syntax)

    except Exception as e:
        console.print(f"❌ Error creating alias from recipe: {e}", style="red")


@history.command()
def latest():
    """Show the most recent execution."""
    executions_list = history_manager.get_executions(1)

    if not executions_list:
        console.print("📝 No executions found", style="yellow")
        return

    execution = executions_list[0]
    recipe_hash = execution.get("recipe_hash")

    console.print("📋 Latest Execution", style="bold")
    console.print()

    # Show execution details
    context = execution.get("context", {})
    exec_details = execution.get("execution", {})

    console.print(f"Recipe: {recipe_hash}")
    console.print(f"Command: {context.get('command_line', 'unknown')}")
    console.print(f"Directory: {context.get('working_directory', 'unknown')}")
    console.print(f"Status: {exec_details.get('status', 'unknown')}")
    console.print(f"Time: {exec_details.get('timestamp', 'unknown')}")

    # Show recipe details if available
    recipe = history_manager.get_recipe_by_hash(recipe_hash)
    if recipe:
        console.print("\n🏗️ Recipe Structure:", style="bold")
        recipe_copy = recipe.copy()
        if "metadata" in recipe_copy:
            del recipe_copy["metadata"]

        yaml_content = OmegaConf.to_yaml(OmegaConf.create(recipe_copy))
        syntax = Syntax(yaml_content, "yaml", theme="github-dark", line_numbers=True)
        console.print(syntax)
