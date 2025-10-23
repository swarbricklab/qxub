"""
qxub.config - Configuration Management Package

This package provides comprehensive configuration management for qxub, including
hierarchical configuration loading, user preferences, shortcuts, and aliases.

Package Structure:
    base.py         - Core configuration data structures and constants
    manager.py      - Configuration loading, merging, and template resolution
    handler.py      - Configuration file operations and validation
    shortcuts.py    - Shortcut management and resolution
    aliases.py      - Alias system management (optional standalone)

Public API:
    # Core Configuration Management
    ConfigManager                         # Main configuration interface
    config_manager                        # Global configuration manager instance

    # Configuration Operations
    get_config(key_path) -> Any          # Get configuration value
    set_config(key_path, value)          # Set configuration value
    validate_config() -> List[str]        # Validate configuration files

    # User Preferences
    get_user_preference(key) -> Any       # Get user-specific setting
    set_user_preference(key, value)       # Set user-specific setting

    # Shortcuts and Aliases
    ShortcutManager                       # Shortcut management interface
    get_shortcut(name) -> Dict           # Resolve shortcut by name
    list_shortcuts() -> List[str]         # List available shortcuts

Example Usage:
    from qxub.config import config_manager, get_config, ShortcutManager

    # Access configuration
    project = get_config("defaults.project")
    walltime = get_config("defaults.walltime")

    # Manage shortcuts
    shortcuts = ShortcutManager()
    python_shortcut = shortcuts.get_shortcut("python")

    # Configuration validation
    issues = config_manager.validate_config()
    if issues:
        print("Configuration issues:", issues)
"""

# Import all public APIs for convenient access
from .base import (  # Legacy configuration (deprecated)
    QxubConfig,
    get_effective_config,
    setup_logging,
)
from .manager import config_manager  # Global instance
from .manager import ConfigManager
from .shortcuts import ShortcutManager

# Optional aliases support (may be moved later)
try:
    from .aliases import StandaloneAliases

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
    # Legacy (deprecated)
    "QxubConfig",
    "get_effective_config",
    "setup_logging",
]

# Add aliases if available
if ALIASES_AVAILABLE:
    __all__.append("StandaloneAliases")


# Convenience functions that delegate to config_manager
def get_config(key_path: str, default=None):
    """Get configuration value by key path."""
    return config_manager.get(key_path, default)


def set_config(key_path: str, value, scope: str = "user"):
    """Set configuration value by key path."""
    return config_manager.set(key_path, value, scope)


def validate_config():
    """Validate configuration files and return list of issues."""
    return config_manager.validate()


def get_user_preference(key: str, default=None):
    """Get user-specific preference."""
    return config_manager.get_user_preference(key, default)


def set_user_preference(key: str, value):
    """Set user-specific preference."""
    return config_manager.set_user_preference(key, value)


# Add convenience functions to __all__
__all__.extend(
    [
        "get_config",
        "set_config",
        "validate_config",
        "get_effective_config",
        "setup_logging",
        "get_user_preference",
        "set_user_preference",
    ]
)
