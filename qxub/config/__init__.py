"""
qxub.config - Configuration Management Package

This package provides comprehensive configuration management for qxub, including
hierarchical configuration loading, user preferences, shortcuts, and aliases.

Package Structure:
    base.py         - Logging set-up for the CLI
    manager.py      - Configuration loading, merging, and template resolution
    handler.py      - Configuration file operations and validation
    shortcuts.py    - Shortcut management and resolution
    aliases.py      - Alias system management (optional standalone)

Public API:
    # Core Configuration Management
    ConfigManager                         # Main configuration interface
    config_manager                        # Global configuration manager instance

    # Configuration Operations
    config_manager.get_config_value(key_path) -> Any

    # Shortcuts and Aliases
    ShortcutManager                       # Shortcut management interface
    get_shortcut(name) -> Dict           # Resolve shortcut by name
    list_shortcuts() -> List[str]         # List available shortcuts

Example Usage:
    from qxub.config import config_manager, ShortcutManager

    # Access configuration
    project = config_manager.get_config_value("defaults.project")

    # Manage shortcuts
    shortcuts = ShortcutManager()
    python_shortcut = shortcuts.get_shortcut("python")

    # Configuration validation
    issues = config_manager.validate()
    if issues:
        print("Configuration issues:", issues)
"""

# Import all public APIs for convenient access
from .base import setup_logging
from .manager import config_manager  # Global instance
from .manager import ConfigManager
from .shortcuts import ShortcutManager

# Optional aliases support (may be moved later)
try:
    from .aliases import StandaloneAliases  # noqa: F401

    ALIASES_AVAILABLE = True
except ImportError:
    ALIASES_AVAILABLE = False

# Define what gets imported with "from qxub.config import *"
__all__ = [
    # Core configuration
    "ConfigManager",
    "config_manager",
    # Shortcuts management
    "ShortcutManager",
    "setup_logging",
]

# Add aliases if available
if ALIASES_AVAILABLE:
    __all__.append("StandaloneAliases")
