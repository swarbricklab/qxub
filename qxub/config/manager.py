"""
Configuration manager for qxub with XDG-compliant paths and template resolution.

Handles loading and merging configuration from system, user, and environment sources
with hierarchical precedence and template variable substitution.
"""

# pylint: disable=no-else-return,no-else-continue,broad-exception-caught
# pylint: disable=too-many-return-statements


import os
import pwd
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import click
from omegaconf import DictConfig, OmegaConf

from .files import (
    file_trust_problem,
    open_trusted_dir,
    open_trusted_file,
    write_private_file,
    xdg_config_dirs,
)

# pylint: disable=broad-exception-caught,no-else-return,too-many-return-statements
# pylint: disable=no-else-continue,duplicate-code

# A project's .qx files come with the directory the user happens to be in: a
# repository they cloned, or a shared area someone else can write. Unless the
# user lists the project under security.trusted_projects, these keys are
# dropped from them, because they decide where and how commands run (ssh
# destination, remote bootstrap commands, platform files) or who may see
# a job's settings.
PROJECT_RESTRICTED_KEYS = (
    "platforms",
    "platform_search_paths",
    "defaults.platform",
    "default_platform",
    "defaults.allow_qsub_v_settings",
)
# Likewise an untrusted project may add aliases, but not change one of the
# same name from the system, user or --config file (its command, pre, post).
ALIASES_KEY = "aliases"
# Never taken from a project, trusted or not: a project can't trust itself.
PROJECT_FORBIDDEN_KEYS = ("security",)
TRUSTED_PROJECTS_KEY = "security.trusted_projects"
# Notification credentials and endpoints come only from the user's own
# config (or the system config), never from a project's files.
PROJECT_SECRET_KEYS = ("notifications",)
# Keys whose values are credentials; shown masked by `qxub config get/list`
SECRET_KEY_NAMES = frozenset(
    {"bot_token", "webhook_url", "token", "password", "secret", "api_key"}
)
MASKED = "********"


def _interpolated_paths(value: Any, path: str = "") -> List[str]:
    """Dotted paths of values in a plain container that use ``${...}``."""
    if isinstance(value, dict):
        found = []
        for key, item in value.items():
            found += _interpolated_paths(item, f"{path}.{key}" if path else str(key))
        return found
    if isinstance(value, list):
        return [path] if any(_interpolated_paths(item, "-") for item in value) else []
    if isinstance(value, str) and "${" in value:
        return [path]
    return []


def _drop_interpolations(value: Any) -> Any:
    """``value`` without any entry that uses ``${...}``."""
    if isinstance(value, dict):
        return {
            key: _drop_interpolations(item)
            for key, item in value.items()
            if not _interpolated_paths(item, "-") or isinstance(item, dict)
        }
    return value


def mask_secrets(value: Any) -> Any:
    """A copy of a plain container with credential values masked."""
    if isinstance(value, dict):
        return {
            key: (
                MASKED
                if str(key).rsplit(".", 1)[-1] in SECRET_KEY_NAMES
                and item not in (None, "")
                and not isinstance(item, (dict, list))
                else mask_secrets(item)
            )
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [mask_secrets(item) for item in value]
    return value


# Modes for the .qx directory and files qxub creates, whatever the umask:
# group- or world-writable project config is ignored as untrusted.
PROJECT_DIR_MODE = 0o755
PROJECT_FILE_MODE = 0o644


def _save_project_file(conf: Any, path: Path) -> None:
    """Save project config ``conf`` to ``path``; a new file is made 0644."""
    existed = path.exists()
    OmegaConf.save(conf, path)
    if not existed:
        os.chmod(path, PROJECT_FILE_MODE)


class ConfigManager:
    """Manages qxub configuration loading, merging, and template resolution."""

    def __init__(self, override_config_path: Optional[str] = None):
        self.system_config: Optional[DictConfig] = None
        self.user_config: Optional[DictConfig] = None
        self.project_config: Optional[DictConfig] = None
        self.local_config: Optional[DictConfig] = None
        self.test_config: Optional[DictConfig] = None
        self.override_config: Optional[DictConfig] = None
        self.merged_config: Optional[DictConfig] = None
        self.override_config_path = override_config_path
        self._load_configs()

    def _get_xdg_config_dirs(self) -> List[Path]:
        """Get XDG config directories in precedence order (absolute only)."""
        return [d / "qxub" for d in xdg_config_dirs()]

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

    def _trusted_project_roots(self) -> List[Path]:
        """Projects the user trusts, from the system, user and --config files."""
        roots = []
        for layer in (self.system_config, self.user_config, self.override_config):
            if layer is None:
                continue
            try:
                value = OmegaConf.select(layer, TRUSTED_PROJECTS_KEY, default=None)
                if OmegaConf.is_config(value):
                    value = OmegaConf.to_container(value, resolve=True)
            except Exception as e:
                click.echo(f"Warning: Bad {TRUSTED_PROJECTS_KEY}: {e}", err=True)
                continue
            if isinstance(value, str):
                value = [value]
            for entry in value or []:
                if entry:
                    roots.append(Path(os.path.expanduser(str(entry))).resolve())
        return roots

    def _is_trusted_project(self, project_root: Path) -> bool:
        """Whether the user listed ``project_root`` (or a parent) as trusted."""
        project_root = project_root.resolve()
        return any(
            root == project_root or root in project_root.parents
            for root in self._trusted_project_roots()
        )

    def _own_alias_names(self) -> set:
        """Alias names defined in the system, user and --config files."""
        names = set()
        for layer in (self.system_config, self.user_config, self.override_config):
            if layer is None or ALIASES_KEY not in layer:
                continue
            aliases = OmegaConf.to_container(layer, resolve=False)[ALIASES_KEY]
            if isinstance(aliases, dict):
                names.update(str(name) for name in aliases)
        return names

    def _restrict_project_layer(
        self, config: DictConfig, config_path: Path, trusted: bool
    ) -> DictConfig:
        """Drop the keys a project's config may not set (see module top)."""
        keys = PROJECT_FORBIDDEN_KEYS
        if not trusted:
            keys = keys + PROJECT_RESTRICTED_KEYS
        removed = []
        for key in keys:
            parent_key, _, leaf = key.rpartition(".")
            parent = OmegaConf.select(config, parent_key) if parent_key else config
            if isinstance(parent, DictConfig) and leaf in parent:
                del parent[leaf]
                removed.append(key)
        if not trusted and ALIASES_KEY in config:
            aliases = OmegaConf.to_container(config, resolve=False)[ALIASES_KEY]
            if not isinstance(aliases, dict):
                # Not a mapping of new aliases: it could only replace the user's
                del config[ALIASES_KEY]
                removed.append(ALIASES_KEY)
            else:
                own = self._own_alias_names()
                for name in aliases:
                    if str(name) in own:
                        del config[ALIASES_KEY][name]
                        removed.append(f"{ALIASES_KEY}.{name}")
        if removed:
            click.echo(
                f"Warning: Ignoring {', '.join(removed)} in {config_path}: "
                "project config can't choose where or how jobs run. If you "
                "trust this project, add its directory to "
                f"{TRUSTED_PROJECTS_KEY} in your user config.",
                err=True,
            )

        secrets = [key for key in PROJECT_SECRET_KEYS if key in config]
        for key in secrets:
            del config[key]
        if secrets:
            click.echo(
                f"Warning: Ignoring {', '.join(secrets)} in {config_path}: "
                "set notification settings in your user config.",
                err=True,
            )

        # ${...} would be resolved against the user's environment and their
        # private config (oc.env, references to notifications.*), so an
        # untrusted project's values are taken as plain data or not at all.
        if not trusted:
            plain = OmegaConf.to_container(config, resolve=False)
            interpolated = _interpolated_paths(plain)
            if interpolated:
                click.echo(
                    f"Warning: Ignoring {', '.join(interpolated)} in "
                    f"{config_path}: ${{...}} interpolation is only allowed "
                    "in your own config or a trusted project's. If you trust "
                    f"this project, add its directory to {TRUSTED_PROJECTS_KEY} "
                    "in your user config.",
                    err=True,
                )
                config = OmegaConf.create(_drop_interpolations(plain))
        return config

    def _load_project_configs(
        self,
    ) -> tuple[Optional[DictConfig], Optional[DictConfig], Optional[DictConfig]]:
        """Load project, local, and test configurations.

        A .qx directory or file that someone other than the user (or root)
        owns or can write, or that is a symlink, is ignored, unless the user
        trusts the project.
        """
        project_dir = self._get_project_config_dir()
        if not project_dir:
            return None, None, None

        trusted = self._is_trusted_project(project_dir.parent)
        config_files = ["project.yaml", "local.yaml", "test.yaml"]
        if trusted:
            configs = [
                self._load_config_file(project_dir / name, name)
                for name in config_files
            ]
        else:
            configs = self._load_untrusted_project_files(project_dir, config_files)
            if configs is None:
                return None, None, None

        return tuple(
            (
                self._restrict_project_layer(config, project_dir / name, trusted)
                if config is not None
                else None
            )
            for config, name in zip(configs, config_files)
        )

    def _load_untrusted_project_files(
        self, project_dir: Path, names: List[str]
    ) -> Optional[List[Optional[DictConfig]]]:
        """Load ``names`` from ``project_dir`` if neither can be changed by others.

        The directory and each file are opened without following symlinks,
        checked through the open descriptors and read from them, so nothing
        can be swapped in between the check and the read. None if the
        directory itself is not to be trusted.
        """
        dir_fd, problem = open_trusted_dir(project_dir)
        if problem:
            click.echo(
                f"Warning: Ignoring project config in {project_dir}: "
                f"{problem}. If you trust this project, add its directory "
                f"to {TRUSTED_PROJECTS_KEY} in your user config.",
                err=True,
            )
            return None
        configs = []
        try:
            for name in names:
                config_path = project_dir / name
                config = None
                f, problem = open_trusted_file(dir_fd, name, config_path)
                if problem:
                    click.echo(f"Warning: Ignoring {config_path}: {problem}", err=True)
                elif f is not None:
                    with f:
                        try:
                            config = OmegaConf.load(f)
                        except Exception as e:
                            click.echo(
                                f"Warning: Failed to load {name} config "
                                f"{config_path}: {e}",
                                err=True,
                            )
                configs.append(config)
        finally:
            os.close(dir_fd)
        return configs

    def _load_config_file(
        self, config_path: Path, config_name: str
    ) -> Optional[DictConfig]:
        """Load a single config file with error handling."""
        if not config_path.exists():
            return None

        try:
            return OmegaConf.load(config_path)
        except Exception as e:
            click.echo(
                f"Warning: Failed to load {config_name} config {config_path}: {e}",
                err=True,
            )
            return None

    def _load_override_config(self) -> Optional[DictConfig]:
        """Load override configuration from CLI --config option."""
        if not self.override_config_path:
            return None

        config_path = Path(self.override_config_path)
        if not config_path.exists():
            click.echo(
                f"Error: Override config file not found: {config_path}", err=True
            )
            return None

        try:
            return OmegaConf.load(config_path)
        except Exception as e:
            click.echo(
                f"Error: Failed to load override config {config_path}: {e}", err=True
            )
            return None

    def _load_system_config(
        self,
    ) -> Optional[DictConfig]:  # pylint: disable=broad-exception-caught
        """Load system-wide configuration."""
        for config_dir in self._get_xdg_config_dirs():
            config_file = config_dir / "config.yaml"
            if config_file.exists():
                # Anyone who can change a shared file could steer every
                # user's jobs, so a writable one is not used.
                problem = file_trust_problem(config_file, require_owner=False)
                if problem:
                    click.echo(
                        f"Warning: Ignoring system config {config_file}: {problem}",
                        err=True,
                    )
                    continue
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
        # Before the project files, which it can mark as trusted
        self.override_config = self._load_override_config()
        self.project_config, self.local_config, self.test_config = (
            self._load_project_configs()
        )

        # Merge configurations in precedence order:
        # system < user < project < local < test < override
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
        if self.override_config:
            configs.append(self.override_config)

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
            max_iterations: Maximum number of recursive resolution attempts to
                prevent infinite loops
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
                    f"Warning: Template resolution hit maximum iterations "
                    f"({max_iterations}) for '{value}'. Possible circular reference.",
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
        self.save_user_config(user_config)

        # Reload configs to pick up changes
        self.reload_configs()

    def save_user_config(self, user_config: DictConfig) -> Path:
        """Write the user config file, readable by the user only.

        It can hold notification tokens and webhook URLs. Interpolations
        such as ``${oc.env:TOKEN}`` are kept as written, not resolved.
        """
        user_config_file = self._get_user_config_dir() / "config.yaml"
        write_private_file(user_config_file, OmegaConf.to_yaml(user_config))
        return user_config_file

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
        _save_project_file(project_config, project_config_file)

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
        _save_project_file(local_config, local_config_file)

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
        _save_project_file(test_config, test_config_file)

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

        # Create .qx directory. 0755 whatever the umask: a group-writable
        # .qx would be ignored as untrusted.
        qx_dir.mkdir(parents=True, exist_ok=True)
        os.chmod(qx_dir, PROJECT_DIR_MODE)

        # Create default project.yaml (git-tracked, team-shared settings)
        project_config_file = qx_dir / "project.yaml"
        project_defaults = {
            "qxub": {"defaults": {"walltime": "1:00:00", "queue": "auto"}}
        }
        _save_project_file(OmegaConf.create(project_defaults), project_config_file)

        # Create default test.yaml (git-tracked, CI/testing settings)
        test_config_file = qx_dir / "test.yaml"
        test_defaults = {
            "qxub": {"defaults": {"walltime": "0:10:00", "queue": "copyq"}}
        }
        _save_project_file(OmegaConf.create(test_defaults), test_config_file)

        # Create .gitignore for local.yaml (git-ignored, user-specific settings)
        gitignore_file = qx_dir / ".gitignore"
        gitignore_file.write_text("local.yaml\n")
        os.chmod(gitignore_file, PROJECT_FILE_MODE)

        return True

    def create_user_config_template(self):
        """Create a template user config file."""
        user_config_file = self._get_user_config_dir() / "config.yaml"

        if user_config_file.exists():
            raise click.ClickException(
                f"User config file already exists: {user_config_file}"
            )

        # Create template config. Site-specific values (project, storage
        # volumes, bind paths) are left for the user to fill in; the platform
        # definition supplies scratch paths and the default project.
        defaults = {
            "name": "qt",
            "queue": "normal",
            "joblog": "{name}.log",
            "resources": ["mem=4GB", "ncpus=1"],
            # Workflow-friendly resource defaults
            "mem": "4GB",
            "cpus": 2,
            "runtime": "1h",
            "disk": "10GB",
            "conda": {"env": "base", "pre": None, "post": None},
            "module": {"mod": ["python3"], "pre": None, "post": None},
            "sing": {
                "sif": None,
                "bind": [],
                "env": [],
                "pre": None,
                "post": None,
            },
        }
        from .site import default_project

        project = default_project()
        if project:
            defaults["project"] = project

        template = {
            "defaults": defaults,
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
        return self.save_user_config(config)

    def save_alias(self, alias_name: str, alias_definition: Dict[str, Any]):
        """Save an alias to the user configuration."""
        user_config_file = self._get_user_config_dir() / "config.yaml"

        # Load existing user config or create new. A file that can't be
        # read is left alone rather than replaced by one holding the alias.
        if user_config_file.exists():
            user_config = OmegaConf.load(user_config_file)
        else:
            user_config = OmegaConf.create({})

        # Ensure aliases section exists
        if "aliases" not in user_config:
            user_config.aliases = {}

        # Add the new alias
        user_config.aliases[alias_name] = alias_definition

        # Save back to file, leaving the user's other values unresolved
        self.save_user_config(user_config)

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

    def get_platform_config(self, platform_name: str) -> Optional[Dict[str, Any]]:
        """
        Get platform configuration from the platforms registry.

        Args:
            platform_name: Name of the platform

        Returns:
            Platform configuration dict or None if not found

        Platform config structure:
            definition: Path/URL to platform definition file (optional for delegation)
            remote: Optional remote execution config
                host: SSH hostname (from ~/.ssh/config)
                working_dir: Remote working directory
                conda_init: Optional conda initialization commands

        Note: Platform name is the YAML key, no separate 'name' field needed.
        """
        if not self.merged_config or "platforms" not in self.merged_config:
            return None

        platforms = self.merged_config.platforms
        if platform_name not in platforms:
            return None

        platform_config = OmegaConf.to_container(platforms[platform_name], resolve=True)
        return platform_config

    def list_platforms(self) -> List[str]:
        """
        List all configured platforms.

        Returns:
            List of platform names
        """
        if not self.merged_config or "platforms" not in self.merged_config:
            return []
        return list(self.merged_config.platforms.keys())

    def get_default_platform_name(self) -> Optional[str]:
        """
        Get the default platform name from config.

        Returns:
            Platform name from defaults.platform, or None
        """
        defaults = self.get_defaults()
        return defaults.get("platform")

    def validate(self) -> List[str]:
        """Validate configuration and return list of issues."""
        issues = []
        if not self.merged_config:
            issues.append("No configuration loaded")
        return issues


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


def set_config_override(override_config_path: Optional[str] = None):
    """
    Reinitialize the global config manager with an override config file.

    This is used by the --config CLI option to load a custom configuration
    for testing or alternative environments.

    Args:
        override_config_path: Path to override config file (highest precedence)
    """
    global config_manager
    config_manager = ConfigManager(override_config_path=override_config_path)


# Global config manager instance
config_manager = ConfigManager()
