"""
Command history logging system for qxub.

Logs every qxub command in YAML format alongside the config file,
using the same hierarchical structure as aliases for easy conversion.
"""

import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import click
from omegaconf import DictConfig, OmegaConf


class CommandHistoryLogger:
    """Logs qxub commands in YAML format using alias-like structure."""

    def __init__(self):
        self.history_file = self._get_history_file_path()

    def _get_history_file_path(self) -> Path:
        """Get the history file path in user config directory."""
        xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
        if xdg_config_home:
            config_dir = Path(xdg_config_home) / "qxub"
        else:
            config_dir = Path.home() / ".config" / "qxub"

        config_dir.mkdir(parents=True, exist_ok=True)
        return config_dir / "history.yaml"

    def _extract_command_info(self, ctx: click.Context) -> Dict[str, Any]:
        """Extract command information from Click context."""
        # Get the full command line as originally typed
        command_line = " ".join(sys.argv)

        # Parse the command structure
        command_parts = []
        current_ctx = ctx

        # Walk up the context chain to get the full command path
        while current_ctx:
            if current_ctx.info_name != "qxub":  # Skip the root command
                command_parts.insert(0, current_ctx.info_name)
            current_ctx = current_ctx.parent

        # Determine the subcommand type
        subcommand_type = None
        if len(command_parts) > 0:
            first_part = command_parts[0]
            if first_part in ["conda", "module", "sing"]:
                subcommand_type = first_part
            elif first_part == "alias":
                subcommand_type = "alias"
            elif first_part in ["config", "history"]:
                subcommand_type = "config"

        return {
            "command_line": command_line,
            "subcommand_type": subcommand_type,
            "command_parts": command_parts,
        }

    def _create_history_entry(self, ctx, success=True):
        """Create a history entry from Click context and command line"""
        import sys

        entry = {
            "timestamp": datetime.now().isoformat(),
            "command_line": " ".join(sys.argv),
            "working_directory": os.getcwd(),
            "main": {},
            "success": success,
        }

        # Extract main command parameters from context
        if ctx.params:
            entry["main"] = ctx.params.copy()

        # Parse command line to find subcommand and target
        if len(sys.argv) > 1:
            args = sys.argv[1:]  # Skip script name

            # Find subcommand (conda, module, sing)
            subcommand_types = ["conda", "module", "sing"]
            subcommand_idx = None
            subcommand_type = None

            for i, arg in enumerate(args):
                if arg in subcommand_types:
                    subcommand_idx = i
                    subcommand_type = arg
                    break

            if subcommand_idx is not None:
                # Parse subcommand arguments and find target
                subcommand_params = {}
                i = subcommand_idx + 1

                while i < len(args):
                    arg = args[i]
                    if arg == "--":
                        # Skip the separator and start target from next argument
                        i += 1
                        if i < len(args):
                            entry["target"] = args[i:]
                        break
                    elif arg.startswith("--"):
                        # This is a subcommand option
                        option_name = arg[2:]  # Remove --
                        if (
                            i + 1 < len(args)
                            and not args[i + 1].startswith("-")
                            and args[i + 1] != "--"
                        ):
                            # Has a value
                            subcommand_params[option_name] = args[i + 1]
                            i += 2
                        else:
                            # Flag without value
                            subcommand_params[option_name] = True
                            i += 1
                    elif not arg.startswith("-"):
                        # This is the start of the target command
                        entry["target"] = args[i:]
                        break
                    else:
                        i += 1

                # Set subcommand section
                entry["subcommand"] = {"type": subcommand_type, **subcommand_params}

        return entry

    def log_command(
        self, ctx: click.Context, success: bool = True, error: Optional[str] = None
    ):
        """Log a command execution to the history file."""
        try:
            history_entry = self._create_history_entry(ctx)
            history_entry["success"] = success
            if error:
                history_entry["error"] = error

            # Load existing history or create new
            if self.history_file.exists():
                try:
                    history_config = OmegaConf.load(self.history_file)
                    if not history_config:
                        history_config = OmegaConf.create({})
                except Exception:
                    # If file is corrupted, start fresh
                    history_config = OmegaConf.create({})
            else:
                history_config = OmegaConf.create({})

            # Ensure history key exists
            if "history" not in history_config:
                history_config.history = []

            # Add new entry
            history_config.history.append(history_entry)

            # Keep only last 1000 entries to prevent file from growing too large
            if len(history_config.history) > 1000:
                history_config.history = history_config.history[-1000:]

            # Save to file
            with open(self.history_file, "w", encoding="utf-8") as f:
                OmegaConf.save(history_config, f, resolve=True)

        except Exception as e:
            # Don't let history logging break the main command
            # Just log the error and continue
            import logging

            logging.debug("Failed to log command to history: %s", str(e))

    def get_recent_commands(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent commands from history."""
        if not self.history_file.exists():
            return []

        try:
            history_config = OmegaConf.load(self.history_file)
            if "history" not in history_config:
                return []

            # Return most recent commands
            history_list = history_config.history
            return history_list[-limit:] if len(history_list) > limit else history_list

        except Exception:
            return []

    def clear_history(self):
        """Clear the command history."""
        if self.history_file.exists():
            self.history_file.unlink()


# Global instance
history_logger = CommandHistoryLogger()
