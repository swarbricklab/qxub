"""
Configuration manager for qxub with XDG-compliant paths and template resolution.

Handles loading and merging configuration from system, user, and environment sources
with hierarchical precedence and template variable substitution.
"""

# pylint: disable=no-else-return,no-else-continue,broad-exception-caught,too-many-return-statements


import os
import pwd
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import click
from omegaconf import DictConfig, OmegaConf

# pylint: disable=broad-exception-caught,no-else-return,too-many-return-statements,no-else-continue,duplicate-code


class ConfigManager:
    """Manages qxub configuration loading, merging, and template resolution."""

    def __init__(self):
        self.system_config: Optional[DictConfig] = None
        self.user_config: Optional[DictConfig] = None
        self.project_config: Optional[DictConfig] = None
        self.local_config: Optional[DictConfig] = None
        self.test_config: Optional[DictConfig] = None
        self.merged_config: Optional[DictConfig] = None
        self._load_configs()

    def _get_xdg_config_dirs(self) -> List[Path]:
        """Get XDG config directories in precedence order."""
        xdg_config_dirs = os.environ.get("XDG_CONFIG_DIRS", "/etc/xdg")
        dirs = [Path(d) / "qxub" for d in xdg_config_dirs.split(":")]
        return dirs

    def _get_user_config_dir(self) -> Path:
        """Get user config directory following XDG spec."""
        xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
        if xdg_config_home:
            return Path(xdg_config_home) / "qxub"
        return Path.home() / ".config" / "qxub"

    def _find_project_root(self, start_path: Optional[Path] = None) -> Optional[Path]:
        """Find project root by searching for .qx directory or git/dvc markers."""
        if start_path is None:
            start_path = Path.cwd()

        current = start_path.resolve()

        # Search up the directory tree
        while current != current.parent:
            # Look for .qx directory first (most specific)
            qx_dir = current / ".qx"
            if qx_dir.exists() and qx_dir.is_dir():
                return current

            # Look for project markers (.git, .dvc, etc.)
            project_markers = [".git", ".dvc", "pyproject.toml", "setup.py"]
            if any((current / marker).exists() for marker in project_markers):
                return current

            current = current.parent

        return None

    def _get_project_config_dir(
        self, start_path: Optional[Path] = None
    ) -> Optional[Path]:
        """Get project config directory (.qx) if it exists."""
        project_root = self._find_project_root(start_path)
        if project_root:
            qx_dir = project_root / ".qx"
            if qx_dir.exists():
                return qx_dir
        return None

    def _load_project_configs(
        self,
    ) -> tuple[Optional[DictConfig], Optional[DictConfig], Optional[DictConfig]]:
        """Load project, local, and test configurations."""
        project_config = None
        local_config = None
        test_config = None

        project_dir = self._get_project_config_dir()
        if project_dir:
            # Load project config (.qx/project.yaml)
            project_file = project_dir / "project.yaml"
            if project_file.exists():
                try:
                    project_config = OmegaConf.load(project_file)
                except Exception as e:
                    click.echo(
                        f"Warning: Failed to load project config {project_file}: {e}",
                        err=True,
                    )

            # Load local config (.qx/local.yaml)
            local_file = project_dir / "local.yaml"
            if local_file.exists():
                try:
                    local_config = OmegaConf.load(local_file)
                except Exception as e:
                    click.echo(
                        f"Warning: Failed to load local config {local_file}: {e}",
                        err=True,
                    )

            # Load test config (.qx/test.yaml)
            test_file = project_dir / "test.yaml"
            if test_file.exists():
                try:
                    test_config = OmegaConf.load(test_file)
                except Exception as e:
                    click.echo(
                        f"Warning: Failed to load test config {test_file}: {e}",
                        err=True,
                    )

        return project_config, local_config, test_config

    def _load_system_config(
        self,
    ) -> Optional[DictConfig]:  # pylint: disable=broad-exception-caught
        """Load system-wide configuration."""
        for config_dir in self._get_xdg_config_dirs():
            config_file = config_dir / "config.yaml"
            if config_file.exists():
                try:
                    return OmegaConf.load(config_file)
                except Exception as e:  # pylint: disable=broad-exception-caught
                    click.echo(
                        f"Warning: Failed to load system config {config_file}: {e}",
                        err=True,
                    )
        return None

    def _load_user_config(self) -> Optional[DictConfig]:
        """Load user configuration."""
        config_file = self._get_user_config_dir() / "config.yaml"
        if config_file.exists():
            try:
                return OmegaConf.load(config_file)
            except Exception as e:  # pylint: disable=broad-exception-caught
                click.echo(
                    f"Warning: Failed to load user config {config_file}: {e}", err=True
                )
        return None

    def _load_configs(self):
        """Load and merge all configuration files."""
        self.system_config = self._load_system_config()
        self.user_config = self._load_user_config()
        self.project_config, self.local_config, self.test_config = (
            self._load_project_configs()
        )

        # Merge configurations in precedence order:
        # system < user < project < local < test
        configs = []
        if self.system_config:
            configs.append(self.system_config)
        if self.user_config:
            configs.append(self.user_config)
        if self.project_config:
            configs.append(self.project_config)
        if self.local_config:
            configs.append(self.local_config)
        if self.test_config:
            configs.append(self.test_config)

        if configs:
            self.merged_config = OmegaConf.merge(*configs)
        else:
            self.merged_config = OmegaConf.create({})

    def reload_configs(self):
        """Reload configuration files."""
        self._load_configs()

    def get_config_files(self) -> Dict[str, Path]:
        """Get paths to all relevant config files."""
        files = {}

        # System config files
        for i, config_dir in enumerate(self._get_xdg_config_dirs()):
            files[f"system_{i}"] = config_dir / "config.yaml"

        # User config file
        files["user"] = self._get_user_config_dir() / "config.yaml"

        # Project config files
        project_dir = self._get_project_config_dir()
        if project_dir:
            files["project"] = project_dir / "project.yaml"
            files["local"] = project_dir / "local.yaml"
            files["test"] = project_dir / "test.yaml"

        return files

    def get_template_variables(
        self,
        name: Optional[str] = None,
        project: Optional[str] = None,
        queue: Optional[str] = None,
    ) -> Dict[str, str]:
        """Get template variables for string substitution."""
        now = datetime.now()
        user = pwd.getpwuid(os.getuid()).pw_name

        # Start with built-in template variables
        template_vars = {
            "user": user,
            "project": project or "",
            "name": name or "",
            "queue": queue or "",
            "timestamp": now.strftime("%Y%m%d_%H%M%S"),
            "date": now.strftime("%Y%m%d"),
            "time": now.strftime("%H%M%S"),
        }

        # Add custom template variables from config
        if self.merged_config and "templates" in self.merged_config:
            config_templates = OmegaConf.to_container(
                self.merged_config.templates, resolve=True
            )
            if isinstance(config_templates, dict):
                template_vars.update(config_templates)

        return template_vars

    def resolve_templates(
        self, value: Any, template_vars: Dict[str, str], max_iterations: int = 10
    ) -> Any:  # pylint: disable=too-many-return-statements,no-else-return
        """Recursively resolve template variables in configuration values.

        Args:
            value: The value to resolve templates in
            template_vars: Dictionary of template variables
            max_iterations: Maximum number of recursive resolution attempts to prevent infinite loops
        """
        if isinstance(value, str):
            # Only resolve if there are template placeholders
            if "{" in value and "}" in value:
                resolved_value = value
                for iteration in range(max_iterations):
                    try:
                        new_value = resolved_value.format(**template_vars)
                        # If no more template variables to resolve, we're done
                        if new_value == resolved_value or (
                            "{" not in new_value or "}" not in new_value
                        ):
                            return new_value
                        resolved_value = new_value
                    except KeyError as e:
                        click.echo(
                            f"Warning: Unknown template variable {e} in '{resolved_value}'",
                            err=True,
                        )
                        return resolved_value
                # If we hit max iterations, warn about possible infinite loop
                click.echo(
                    f"Warning: Template resolution hit maximum iterations ({max_iterations}) for '{value}'. "
                    f"Possible circular reference.",
                    err=True,
                )
                return resolved_value
            return value
        if isinstance(value, (list, tuple)):
            return [
                self.resolve_templates(item, template_vars, max_iterations)
                for item in value
            ]
        if isinstance(value, dict):
            return {
                k: self.resolve_templates(v, template_vars, max_iterations)
                for k, v in value.items()
            }
        if OmegaConf.is_config(value):
            # Handle OmegaConf DictConfig
            resolved = {}
            for k, v in value.items():
                resolved[k] = self.resolve_templates(v, template_vars, max_iterations)
            return OmegaConf.create(resolved)
        return value

    def get_defaults(self) -> Dict[str, Any]:
        """Get default configuration values."""
        if not self.merged_config or "defaults" not in self.merged_config:
            return {}
        return OmegaConf.to_container(self.merged_config.defaults, resolve=True)

    def get_alias(self, alias_name: str) -> Optional[Dict[str, Any]]:
        """Get alias definition."""
        if (
            not self.merged_config
            or "aliases" not in self.merged_config
            or alias_name not in self.merged_config.aliases
        ):
            return None
        return OmegaConf.to_container(
            self.merged_config.aliases[alias_name], resolve=True
        )

    def list_aliases(self) -> List[str]:
        """List all available aliases."""
        if not self.merged_config or "aliases" not in self.merged_config:
            return []
        return list(self.merged_config.aliases.keys())

    def get_config_value(self, key_path: str) -> Any:
        """Get configuration value by dot-separated path (e.g., 'defaults.name')."""
        if not self.merged_config:
            return None
        try:
            return OmegaConf.select(self.merged_config, key_path)
        except Exception:  # pylint: disable=broad-exception-caught
            return None

    def set_user_config_value(self, key_path: str, value: Any):
        """Set a configuration value in user config file."""
        user_config_file = self._get_user_config_dir() / "config.yaml"

        # Create directory if it doesn't exist
        user_config_file.parent.mkdir(parents=True, exist_ok=True)

        # Load existing user config or create new
        if user_config_file.exists():
            user_config = OmegaConf.load(user_config_file)
        else:
            user_config = OmegaConf.create({})

        # Set the value using OmegaConf.update
        keys = key_path.split(".")
        current = user_config
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        current[keys[-1]] = value

        # Save the config
        OmegaConf.save(user_config, user_config_file)

        # Reload configs to pick up changes
        self.reload_configs()

    def set_system_config_value(self, key_path: str, value: Any):
        """Set a configuration value in system config file."""
        # Get the first (primary) system config directory
        system_config_dirs = self._get_xdg_config_dirs()
        if not system_config_dirs:
            raise click.ClickException("No system config directories found")

        system_config_file = system_config_dirs[0] / "config.yaml"

        # Create directory if it doesn't exist
        system_config_file.parent.mkdir(parents=True, exist_ok=True)

        # Load existing system config or create new
        if system_config_file.exists():
            system_config = OmegaConf.load(system_config_file)
        else:
            system_config = OmegaConf.create({})

        # Set the value using OmegaConf.update
        keys = key_path.split(".")
        current = system_config
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        current[keys[-1]] = value

        # Save the config
        OmegaConf.save(system_config, system_config_file)

        # Reload configs to pick up changes
        self.reload_configs()

    def set_project_config_value(self, key_path: str, value: Any):
        """Set a configuration value in project config file."""
        project_dir = self._get_project_config_dir()
        if not project_dir:
            raise click.ClickException(
                "No project directory found. Use 'qxub init' to initialize a project config."
            )

        project_config_file = project_dir / "project.yaml"

        # Load existing project config or create new
        if project_config_file.exists():
            project_config = OmegaConf.load(project_config_file)
        else:
            project_config = OmegaConf.create({})

        # Set the value
        keys = key_path.split(".")
        current = project_config
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        current[keys[-1]] = value

        # Save the config
        OmegaConf.save(project_config, project_config_file)

        # Reload configs to pick up changes
        self.reload_configs()

    def set_local_config_value(self, key_path: str, value: Any):
        """Set a configuration value in local config file."""
        project_dir = self._get_project_config_dir()
        if not project_dir:
            raise click.ClickException(
                "No project directory found. Use 'qxub init' to initialize a project config."
            )

        local_config_file = project_dir / "local.yaml"

        # Load existing local config or create new
        if local_config_file.exists():
            local_config = OmegaConf.load(local_config_file)
        else:
            local_config = OmegaConf.create({})

        # Set the value
        keys = key_path.split(".")
        current = local_config
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        current[keys[-1]] = value

        # Save the config
        OmegaConf.save(local_config, local_config_file)

        # Reload configs to pick up changes
        self.reload_configs()

    def set_test_config_value(self, key_path: str, value: Any):
        """Set a configuration value in test config file."""
        project_dir = self._get_project_config_dir()
        if not project_dir:
            raise click.ClickException(
                "No project directory found. Use 'qxub init' to initialize a project config."
            )

        test_config_file = project_dir / "test.yaml"

        # Load existing test config or create new
        if test_config_file.exists():
            test_config = OmegaConf.load(test_config_file)
        else:
            test_config = OmegaConf.create({})

        # Set the value
        keys = key_path.split(".")
        current = test_config
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        current[keys[-1]] = value

        # Save the config
        OmegaConf.save(test_config, test_config_file)

        # Reload configs to pick up changes
        self.reload_configs()

    def init_project_config(self, project_root: Path) -> bool:
        """Initialize project configuration directory and files.

        Args:
            project_root: Root directory of the project

        Returns:
            True if initialization was successful, False if already exists
        """
        qx_dir = project_root / ".qx"

        if qx_dir.exists():
            return False  # Already initialized

        # Create .qx directory
        qx_dir.mkdir(parents=True, exist_ok=True)

        # Create default project.yaml (git-tracked, team-shared settings)
        project_config_file = qx_dir / "project.yaml"
        project_defaults = {
            "qxub": {"defaults": {"walltime": "1:00:00", "queue": "auto"}}
        }
        OmegaConf.save(project_defaults, project_config_file)

        # Create default test.yaml (git-tracked, CI/testing settings)
        test_config_file = qx_dir / "test.yaml"
        test_defaults = {
            "qxub": {"defaults": {"walltime": "0:10:00", "queue": "copyq"}}
        }
        OmegaConf.save(test_defaults, test_config_file)

        # Create .gitignore for local.yaml (git-ignored, user-specific settings)
        gitignore_file = qx_dir / ".gitignore"
        gitignore_file.write_text("local.yaml\n")

        return True

    def create_user_config_template(self):
        """Create a template user config file."""
        user_config_file = self._get_user_config_dir() / "config.yaml"

        if user_config_file.exists():
            raise click.ClickException(
                f"User config file already exists: {user_config_file}"
            )

        # Create directory
        user_config_file.parent.mkdir(parents=True, exist_ok=True)

        # Create template config
        template = {
            "defaults": {
                "name": "qt",
                "queue": "normal",
                "project": "a56",
                "joblog": "{name}.log",
                "resources": ["mem=4GB", "ncpus=1"],
                "out": "/scratch/{project}/{user}/qt/{timestamp}/out",
                "err": "/scratch/{project}/{user}/qt/{timestamp}/err",
                "conda": {"env": "base", "pre": None, "post": None},
                "module": {"mod": ["python3"], "pre": None, "post": None},
                "sing": {
                    "sif": None,
                    "bind": ["/scratch", "/g/data"],
                    "env": [],
                    "pre": None,
                    "post": None,
                },
            },
            "aliases": {
                "example": {
                    "cmd": 'echo "Hello from qxub alias!"',
                    "name": "example",
                    "env": "base",  # Direct execution context option
                },
                "module_example": {
                    "cmd": "samtools --version",
                    "name": "samtools_test",
                    "mod": "samtools",  # Single module
                },
                "multi_module_example": {
                    "cmd": "python analysis.py",
                    "name": "analysis",
                    "mods": "python3,gcc",  # Multiple modules
                },
                "container_example": {
                    "cmd": "python script.py",
                    "name": "container_job",
                    "sif": "/containers/analysis.sif",
                    "bind": ["/data:/data"],
                },
                "default_example": {
                    "cmd": "echo 'Direct PBS execution'",
                    "name": "simple_job",
                    # No execution context = default execution
                },
            },
        }

        config = OmegaConf.create(template)
        OmegaConf.save(config, user_config_file)
        return user_config_file

    def resolve_options(
        self, cli_args: Dict[str, Any], alias_name: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Resolve final options from config hierarchy + alias + CLI args.

        Precedence: CLI Args > Alias > User Config > System Config > Built-in Defaults
        """
        # Start with config defaults
        resolved = self.get_defaults().copy()

        # Apply alias if specified
        if alias_name:
            alias_def = self.get_alias(alias_name)
            if not alias_def:
                raise click.ClickException(f"Alias '{alias_name}' not found")

            # Merge alias into resolved options
            for key, value in alias_def.items():
                if key == "subcommand":
                    # Don't merge subcommand into general options
                    continue
                elif key in ["conda", "module", "sing"]:
                    # Merge subcommand-specific options
                    if key not in resolved:
                        resolved[key] = {}
                    resolved[key].update(value)
                else:
                    # General options
                    resolved[key] = value

        # Apply CLI arguments (highest precedence)
        for key, value in cli_args.items():
            if value is not None:  # Only override if explicitly set
                resolved[key] = value

        # Resolve templates
        template_vars = self.get_template_variables(
            name=resolved.get("name"),
            project=resolved.get("project"),
            queue=resolved.get("queue"),
        )
        resolved = self.resolve_templates(resolved, template_vars)

        return resolved

    def save_alias(self, alias_name: str, alias_definition: Dict[str, Any]):
        """Save an alias to the user configuration."""
        user_config_file = self._get_user_config_dir() / "config.yaml"

        # Ensure user config directory exists
        user_config_file.parent.mkdir(parents=True, exist_ok=True)

        # Load existing user config or create new
        if user_config_file.exists():
            try:
                user_config = OmegaConf.load(user_config_file)
            except Exception:
                user_config = OmegaConf.create({})
        else:
            user_config = OmegaConf.create({})

        # Ensure aliases section exists
        if "aliases" not in user_config:
            user_config.aliases = {}

        # Add the new alias
        user_config.aliases[alias_name] = alias_definition

        # Save back to file
        with open(user_config_file, "w", encoding="utf-8") as f:
            OmegaConf.save(user_config, f, resolve=True)

        # Reload configs to pick up the new alias
        self._load_configs()

    def get_platform_search_paths(self) -> List[Path]:
        """Get platform search paths from config, environment, or defaults."""
        import os

        # Check environment variable first
        env_paths = os.getenv("QXUB_PLATFORM_PATHS")
        if env_paths:
            # Support both single path and colon-separated paths
            if ":" in env_paths:
                return [Path(p.strip()) for p in env_paths.split(":")]
            else:
                return [Path(env_paths)]

        # Then check config with template resolution
        configured_paths = self.get_config_value("platform_search_paths")
        if configured_paths is not None and configured_paths:
            # Get template variables (including project from defaults)
            defaults = self.get_defaults()
            project = defaults.get("project", "")
            template_vars = self.get_template_variables(project=project)

            # Resolve templates in each path
            resolved_paths = []
            for path in configured_paths:
                resolved_path = self.resolve_templates(path, template_vars)
                resolved_paths.append(Path(resolved_path))
            return resolved_paths

        # Finally use defaults
        return [
            Path("/etc/qxub/platforms"),
            self._get_user_config_dir() / "platforms",
            Path.home() / ".qxub" / "platforms",
        ]

    def get_default_platform(self) -> Optional[str]:
        """Get the default platform name from config."""
        return self.get_config_value("default_platform")

    def get_platform_preferences(self) -> Dict[str, Any]:
        """Get platform-specific preferences."""
        platform_prefs = self.get_config_value("platform_preferences")
        return platform_prefs if platform_prefs is not None else {}

    def get_queue_preferences(self) -> Dict[str, Any]:
        """Get queue selection preferences."""
        queue_prefs = self.get_config_value("queue_preferences")
        if queue_prefs is not None:
            return queue_prefs
        else:
            return {
                "optimization": "balanced",  # cost, speed, balanced
                "adjustment_policy": "suggest",  # auto, suggest, user, error
                "auto_select": True,
            }


def setup_logging(verbosity: int = None):
    """
    Configures the logging level based on the verbosity provided by the user.

    Args:
        verbosity (int): The number of '-v' flags used, or from config.
                       - 0: ERROR level (default)
                       - 1: WARNING level
                       - 2: INFO level
                       - 3 or more: DEBUG level

    This function adjusts the logging output to provide more detailed information
    as verbosity increases, allowing users to control the granularity of log messages.
    """
    import logging

    if verbosity is None:
        verbosity = config_manager.get_config_value("verbosity") or 0

    # Get the root logger and configure it directly
    root_logger = logging.getLogger()

    # Remove existing handlers to avoid conflicts
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Set the level and format based on verbosity
    if verbosity == 1:
        level = logging.WARNING
        format_str = "%(levelname)s: %(message)s"
    elif verbosity == 2:
        level = logging.INFO
        format_str = "%(levelname)s: %(message)s"
    elif verbosity >= 3:
        level = logging.DEBUG
        format_str = "%(levelname)s:%(name)s: %(message)s"
    else:
        level = logging.ERROR
        format_str = "%(levelname)s: %(message)s"

    # Configure the root logger
    root_logger.setLevel(level)
    handler = logging.StreamHandler()
    handler.setLevel(level)
    formatter = logging.Formatter(format_str)
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)


# Global config manager instance
config_manager = ConfigManager()
