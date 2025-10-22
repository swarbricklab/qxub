"""
New dual-log history system for qxub.

This system separates "what" (recipes) from "how/when" (executions):
- recipes.yaml: Unique computational recipes indexed by hash
- executions.yaml: Individual execution records indexed by microsecond timestamp
"""

import hashlib
import json
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from omegaconf import OmegaConf


class HistoryManager:
    """Manages dual-log history system with recipes and executions."""

    def __init__(self, config_dir: Optional[Path] = None):
        """Initialize history manager with config directory."""
        if config_dir is None:
            # Use XDG config directory
            config_dir = Path.home() / ".config" / "qxub"

        self.config_dir = Path(config_dir)
        self.recipes_file = self.config_dir / "recipes.yaml"
        self.executions_file = self.config_dir / "executions.yaml"

        # Ensure config directory exists
        self.config_dir.mkdir(parents=True, exist_ok=True)

    def compute_recipe_hash(self, recipe_data: Dict[str, Any]) -> str:
        """Compute stable hash for a recipe based on executor and target."""
        # Only hash the executor and target - ignore metadata and context
        hashable = {
            "executor": recipe_data.get("executor", {}),
            "target": recipe_data.get("target", {}),
        }

        # Canonical JSON for stable hashing
        canonical = json.dumps(hashable, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode()).hexdigest()[:12]

    def _extract_recipe_from_command(self, ctx, command_line: str) -> Dict[str, Any]:
        """Extract recipe definition from Click context and command line."""
        import sys

        recipe = {"executor": {}, "target": {}}

        if len(sys.argv) > 1:
            args = sys.argv[1:]  # Skip script name

            # Handle new v3 CLI structure: qxub exec --env test "command"
            if len(args) > 0 and args[0] == "exec":
                args = args[1:]  # Remove 'exec' subcommand

            # For unified CLI (qxub 2.0/3.0), determine executor type from execution context options
            executor_type = None
            executor_params = {}

            # Check for execution context options
            i = 0
            while i < len(args):
                arg = args[i]

                if arg == "--env" and i + 1 < len(args):
                    # Conda environment specified
                    executor_type = "conda"
                    executor_params["env"] = args[i + 1]
                    i += 2
                elif arg in ["--mod", "--mods"] and i + 1 < len(args):
                    # Module specified (handle both --mod and --mods)
                    executor_type = "module"
                    executor_params["name"] = args[i + 1]
                    i += 2
                elif arg == "--sif" and i + 1 < len(args):
                    # Singularity container specified
                    executor_type = "sing"
                    executor_params["image"] = args[i + 1]
                    i += 2
                elif arg == "--":
                    # Found command separator - everything after is the target command
                    i += 1
                    if i < len(args):
                        recipe["target"] = {"cmd": " ".join(args[i:])}
                    break
                elif (
                    arg.startswith("--")
                    and i + 1 < len(args)
                    and not args[i + 1].startswith("-")
                ):
                    # Other executor options
                    option_name = arg[2:]  # Remove --
                    executor_params[option_name] = args[i + 1]
                    i += 2
                elif arg.startswith("--"):
                    # Flag without value
                    option_name = arg[2:]  # Remove --
                    executor_params[option_name] = True
                    i += 1
                elif not arg.startswith("-"):
                    # This is the start of the target command (no -- separator found)
                    recipe["target"] = {"cmd": " ".join(args[i:])}
                    break
                else:
                    i += 1

            # Set executor type and params
            if executor_type:
                executor_params["type"] = executor_type
                recipe["executor"] = executor_params

            # If no explicit executor type found, check for legacy subcommand format for backward compatibility
            if not executor_type:
                executor_types = ["conda", "module", "sing"]
                for j, arg in enumerate(args):
                    if arg in executor_types:
                        executor_type = arg
                        recipe["executor"] = {"type": executor_type}
                        # Parse remaining arguments for this legacy format
                        k = j + 1
                        while k < len(args):
                            if args[k] == "--":
                                k += 1
                                if k < len(args):
                                    recipe["target"] = {"cmd": " ".join(args[k:])}
                                break
                            elif not args[k].startswith("-"):
                                recipe["target"] = {"cmd": " ".join(args[k:])}
                                break
                            k += 1
                        break

        return recipe

    def log_execution(
        self,
        ctx,
        success: bool = True,
        error: Optional[str] = None,
        job_id: Optional[str] = None,
        resource_data: Optional[Dict] = None,
    ) -> str:
        """Log a command execution and return the execution timestamp."""
        try:
            # Generate microsecond timestamp
            execution_timestamp = str(int(time.time() * 1000000))

            # Extract recipe from command
            import sys

            command_line = " ".join(sys.argv)
            recipe_data = self._extract_recipe_from_command(ctx, command_line)

            # Only proceed if we have a valid executor (not config commands, etc.)
            if not recipe_data.get("executor") or recipe_data["executor"].get(
                "type"
            ) not in [
                "conda",
                "module",
                "sing",
            ]:
                return execution_timestamp  # Skip logging for non-execution commands

            # Compute recipe hash
            recipe_hash = self.compute_recipe_hash(recipe_data)

            # Load or create recipes file
            if self.recipes_file.exists():
                try:
                    recipes_config = OmegaConf.load(self.recipes_file)
                    if not recipes_config:
                        recipes_config = OmegaConf.create({"recipes": {}})
                except Exception:
                    recipes_config = OmegaConf.create({"recipes": {}})
            else:
                recipes_config = OmegaConf.create({"recipes": {}})

            # Ensure recipes key exists
            if "recipes" not in recipes_config:
                recipes_config.recipes = {}

            # Add or update recipe
            if recipe_hash not in recipes_config.recipes:
                # New recipe
                recipes_config.recipes[recipe_hash] = {
                    **recipe_data,
                    "metadata": {
                        "first_seen": datetime.now().isoformat(),
                        "last_seen": datetime.now().isoformat(),
                        "run_count": 1,
                    },
                }
            else:
                # Update existing recipe metadata
                recipes_config.recipes[recipe_hash].metadata.last_seen = (
                    datetime.now().isoformat()
                )
                recipes_config.recipes[recipe_hash].metadata.run_count += 1

            # Save recipes file
            with open(self.recipes_file, "w", encoding="utf-8") as f:
                OmegaConf.save(recipes_config, f, resolve=True)

            # Load or create executions file
            if self.executions_file.exists():
                try:
                    executions_config = OmegaConf.load(self.executions_file)
                    if not executions_config:
                        executions_config = OmegaConf.create({"executions": {}})
                except Exception:
                    executions_config = OmegaConf.create({"executions": {}})
            else:
                executions_config = OmegaConf.create({"executions": {}})

            # Ensure executions key exists
            if "executions" not in executions_config:
                executions_config.executions = {}

            # Create execution record
            execution_record = {
                "recipe_hash": recipe_hash,
                "context": {
                    "command_line": command_line,
                    "working_directory": os.getcwd(),
                    "user": os.getenv("USER", "unknown"),
                },
                "requested": {"main_options": ctx.params.copy() if ctx.params else {}},
                "execution": {
                    "job_id": job_id,
                    "joblog": None,  # To be filled in later if available
                    "status": "completed" if success else "failed",
                    "timestamp": datetime.now().isoformat(),
                    "exit_code": (
                        resource_data.get("exit_status") if resource_data else None
                    ),
                },
            }

            # Add resource information if available
            if resource_data:
                execution_record["resources"] = {
                    "requested": resource_data.get("resources_requested", {}),
                    "used": resource_data.get("resources_used", {}),
                    "efficiency": resource_data.get("efficiency", {}),
                    "execution_env": resource_data.get("execution", {}),
                    "timing": resource_data.get("timing", {}),
                }

            if error:
                execution_record["execution"]["error"] = error

            # Add execution record
            executions_config.executions[execution_timestamp] = execution_record

            # Keep only last 1000 executions to prevent file from growing too large
            if len(executions_config.executions) > 1000:
                # Sort by timestamp and keep the most recent 1000
                sorted_timestamps = sorted(executions_config.executions.keys(), key=int)
                old_timestamps = sorted_timestamps[:-1000]
                for ts in old_timestamps:
                    del executions_config.executions[ts]

            # Save executions file
            with open(self.executions_file, "w", encoding="utf-8") as f:
                OmegaConf.save(executions_config, f, resolve=True)

            return execution_timestamp

        except Exception as e:
            # Don't let history logging break the main command
            import logging

            logging.debug("Failed to log execution to history: %s", str(e))
            return str(int(time.time() * 1000000))

    def update_execution_with_resources(
        self, execution_timestamp: str, job_id: str, resource_data: Dict[str, Any]
    ) -> bool:
        """Update an existing execution record with resource data after job completion."""
        if not self.executions_file.exists():
            return False

        try:
            # Load executions
            config = OmegaConf.load(self.executions_file)
            executions = config.get("executions", {})

            if execution_timestamp not in executions:
                return False

            # Update the execution record
            execution = executions[execution_timestamp]
            execution.execution.job_id = job_id
            execution.execution.exit_code = resource_data.get("exit_status")

            # Add resource information
            if resource_data:
                execution.resources = {
                    "requested": resource_data.get("resources_requested", {}),
                    "used": resource_data.get("resources_used", {}),
                    "efficiency": resource_data.get("efficiency", {}),
                    "execution_env": resource_data.get("execution", {}),
                    "timing": resource_data.get("timing", {}),
                }

            # Save updated executions
            with open(self.executions_file, "w", encoding="utf-8") as f:
                OmegaConf.save(config, f, resolve=True)

            return True

        except Exception as e:
            import logging

            logging.debug("Failed to update execution with resources: %s", str(e))
            return False

    def get_recipes(self) -> Dict[str, Any]:
        """Get all recipes."""
        if not self.recipes_file.exists():
            return {}

        try:
            config = OmegaConf.load(self.recipes_file)
            return config.get("recipes", {})
        except Exception:
            return {}

    def get_executions(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent executions."""
        if not self.executions_file.exists():
            return []

        try:
            config = OmegaConf.load(self.executions_file)
            executions = config.get("executions", {})

            # Sort by timestamp (most recent first) and limit
            sorted_timestamps = sorted(executions.keys(), key=int, reverse=True)
            recent_timestamps = sorted_timestamps[:limit]

            result = []
            for ts in recent_timestamps:
                execution = executions[ts].copy()
                execution["timestamp"] = ts
                result.append(execution)

            return result
        except Exception:
            return []

    def get_recipe_by_hash(self, recipe_hash: str) -> Optional[Dict[str, Any]]:
        """Get a specific recipe by hash."""
        recipes = self.get_recipes()
        return recipes.get(recipe_hash)

    def get_executions_for_recipe(
        self, recipe_hash: str, limit: int = 10
    ) -> List[Dict[str, Any]]:
        """Get executions for a specific recipe."""
        if not self.executions_file.exists():
            return []

        try:
            config = OmegaConf.load(self.executions_file)
            executions = config.get("executions", {})

            # Filter executions for this recipe
            matching_executions = []
            for ts, execution in executions.items():
                if execution.get("recipe_hash") == recipe_hash:
                    exec_copy = execution.copy()
                    exec_copy["timestamp"] = ts
                    matching_executions.append(exec_copy)

            # Sort by timestamp (most recent first) and limit
            matching_executions.sort(key=lambda x: int(x["timestamp"]), reverse=True)
            return matching_executions[:limit]
        except Exception:
            return []

    def clear_all_history(self):
        """Clear all history (recipes and executions)."""
        if self.recipes_file.exists():
            self.recipes_file.unlink()
        if self.executions_file.exists():
            self.executions_file.unlink()

    def convert_recipe_to_alias(self, recipe_hash: str) -> Optional[Dict[str, Any]]:
        """Convert a recipe to alias format for backward compatibility."""
        recipe = self.get_recipe_by_hash(recipe_hash)
        if not recipe:
            return None

        # Convert to alias format
        alias_def = {}

        # Map executor to subcommand format for compatibility
        if "executor" in recipe:
            alias_def["subcommand"] = recipe["executor"]

        # Map target format
        if "target" in recipe:
            alias_def["target"] = recipe["target"]

        return alias_def


# Global history manager instance
history_manager = HistoryManager()
