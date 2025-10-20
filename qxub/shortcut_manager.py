"""
Shortcut management system for qxub.

Shortcuts provide fast command execution with pre-configured environments and resources.
Unlike aliases (which are resource profiles), shortcuts are command-specific with
optimized JSON storage for instant loading.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import click


class ShortcutManager:
    """Manages command shortcuts with XDG config resolution."""

    def __init__(self):
        self._shortcuts_cache: Optional[Dict[str, Dict[str, Any]]] = None
        self._system_shortcuts_file = self._get_system_shortcuts_file()
        self._user_shortcuts_file = self._get_user_shortcuts_file()

    def _get_system_shortcuts_file(self) -> Path:
        """Get system shortcuts file location using XDG Base Directory spec."""
        xdg_config_dirs = os.environ.get("XDG_CONFIG_DIRS", "/etc/xdg").split(":")
        return Path(xdg_config_dirs[0]) / "qxub" / "shortcuts.json"

    def _get_user_shortcuts_file(self) -> Path:
        """Get user shortcuts file location using XDG Base Directory spec."""
        xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
        if xdg_config_home:
            return Path(xdg_config_home) / "qxub" / "shortcuts.json"
        return Path.home() / ".config" / "qxub" / "shortcuts.json"

    def _load_shortcuts(self) -> Dict[str, Dict[str, Any]]:
        """Load shortcuts with XDG precedence: user overrides system."""
        shortcuts = {}

        # Load system shortcuts first (lower precedence)
        if self._system_shortcuts_file.exists():
            try:
                with open(self._system_shortcuts_file, "r", encoding="utf-8") as f:
                    system_shortcuts = json.load(f)
                    if isinstance(system_shortcuts, dict):
                        shortcuts.update(system_shortcuts)
            except (json.JSONDecodeError, OSError) as e:
                click.echo(
                    f"⚠️  Warning: Could not load system shortcuts: {e}", err=True
                )

        # Load user shortcuts (higher precedence - overwrites system)
        if self._user_shortcuts_file.exists():
            try:
                with open(self._user_shortcuts_file, "r", encoding="utf-8") as f:
                    user_shortcuts = json.load(f)
                    if isinstance(user_shortcuts, dict):
                        shortcuts.update(user_shortcuts)
            except (json.JSONDecodeError, OSError) as e:
                click.echo(f"⚠️  Warning: Could not load user shortcuts: {e}", err=True)

        return shortcuts

    def _ensure_cache_loaded(self) -> None:
        """Ensure shortcuts cache is loaded."""
        if self._shortcuts_cache is None:
            self._shortcuts_cache = self._load_shortcuts()

    def find_shortcut(self, command_parts: List[str]) -> Optional[Dict[str, Any]]:
        """
        Find longest matching shortcut for command parts.

        Args:
            command_parts: List of command words (e.g., ["dvc", "doctor", "--help"])

        Returns:
            Dict with 'name', 'definition', and 'remaining_args' if found, None otherwise
        """
        self._ensure_cache_loaded()

        if not command_parts or not self._shortcuts_cache:
            return None

        # Try progressively shorter prefixes (longest match wins)
        for length in range(len(command_parts), 0, -1):
            candidate = " ".join(command_parts[:length])

            if candidate in self._shortcuts_cache:
                shortcut_def = self._shortcuts_cache[candidate]
                remaining_args = command_parts[length:]

                return {
                    "name": candidate,
                    "definition": shortcut_def,
                    "remaining_args": remaining_args,
                    "match_length": length,
                }

        return None

    def list_shortcuts(self) -> Dict[str, Dict[str, Any]]:
        """List all available shortcuts."""
        self._ensure_cache_loaded()
        return self._shortcuts_cache.copy() if self._shortcuts_cache else {}

    def get_shortcut(self, name: str) -> Optional[Dict[str, Any]]:
        """Get specific shortcut by name."""
        self._ensure_cache_loaded()
        return self._shortcuts_cache.get(name) if self._shortcuts_cache else None

    def add_shortcut(self, name: str, definition: Dict[str, Any]) -> None:
        """
        Add or update a shortcut in user config.

        Args:
            name: Shortcut name (command to match)
            definition: Shortcut settings (env, resources, etc.)
        """
        # Ensure user config directory exists
        self._user_shortcuts_file.parent.mkdir(parents=True, exist_ok=True)

        # Load existing user shortcuts
        user_shortcuts = {}
        if self._user_shortcuts_file.exists():
            try:
                with open(self._user_shortcuts_file, "r", encoding="utf-8") as f:
                    user_shortcuts = json.load(f)
            except (json.JSONDecodeError, OSError):
                # File exists but is invalid - start fresh
                pass

        # Add/update shortcut
        user_shortcuts[name] = definition

        # Save back to file
        with open(self._user_shortcuts_file, "w", encoding="utf-8") as f:
            json.dump(user_shortcuts, f, indent=2, sort_keys=True)

        # Invalidate cache
        self._shortcuts_cache = None

    def add_system_shortcut(self, name: str, definition: Dict[str, Any]) -> None:
        """
        Add or update a shortcut in system config.

        Args:
            name: Shortcut name (command to match)
            definition: Shortcut settings (env, resources, etc.)

        Raises:
            PermissionError: If no write permissions to system config
        """
        # Ensure system config directory exists
        self._system_shortcuts_file.parent.mkdir(parents=True, exist_ok=True)

        # Load existing system shortcuts
        system_shortcuts = {}
        if self._system_shortcuts_file.exists():
            try:
                with open(self._system_shortcuts_file, "r", encoding="utf-8") as f:
                    system_shortcuts = json.load(f)
            except (json.JSONDecodeError, OSError):
                # File exists but is invalid - start fresh
                pass

        # Add/update shortcut
        system_shortcuts[name] = definition

        # Save back to file
        with open(self._system_shortcuts_file, "w", encoding="utf-8") as f:
            json.dump(system_shortcuts, f, indent=2, sort_keys=True)

        # Invalidate cache
        self._shortcuts_cache = None

    def remove_shortcut(self, name: str) -> bool:
        """
        Remove a shortcut from user config.

        Args:
            name: Shortcut name to remove

        Returns:
            True if shortcut was removed, False if it didn't exist
        """
        if not self._user_shortcuts_file.exists():
            return False

        try:
            with open(self._user_shortcuts_file, "r", encoding="utf-8") as f:
                user_shortcuts = json.load(f)
        except (json.JSONDecodeError, OSError):
            return False

        if name not in user_shortcuts:
            return False

        # Remove shortcut
        del user_shortcuts[name]

        # Save back to file
        with open(self._user_shortcuts_file, "w", encoding="utf-8") as f:
            json.dump(user_shortcuts, f, indent=2, sort_keys=True)

        # Invalidate cache
        self._shortcuts_cache = None
        return True

    def refresh_cache(self) -> None:
        """Force refresh of shortcuts cache."""
        self._shortcuts_cache = None

    def get_config_files(self) -> Dict[str, Tuple[Path, bool]]:
        """Get shortcuts config file locations and existence status."""
        return {
            "system": (
                self._system_shortcuts_file,
                self._system_shortcuts_file.exists(),
            ),
            "user": (self._user_shortcuts_file, self._user_shortcuts_file.exists()),
        }


# Global shortcut manager instance
shortcut_manager = ShortcutManager()
