"""
Configuration management for qxub.

Handles loading and merging configuration from:
1. System config: ${XDG_CONFIG_DIRS}/qxub/config.yaml
2. User config: ~/.config/qxub/config.yaml
3. CLI arguments (highest precedence)

Supports template variable resolution and alias definitions.
# pylint: disable=broad-exception-caught,no-else-return,too-many-return-statements,no-else-continue
"""

import os
import pwd
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
from omegaconf import OmegaConf, DictConfig
import click


class ConfigManager:
    """Manages qxub configuration loading, merging, and template resolution."""

    def __init__(self):
        self.system_config: Optional[DictConfig] = None
        self.user_config: Optional[DictConfig] = None
        self.merged_config: Optional[DictConfig] = None
        self._load_configs()

    def _get_xdg_config_dirs(self) -> List[Path]:
        """Get XDG config directories in precedence order."""
        xdg_config_dirs = os.environ.get('XDG_CONFIG_DIRS', '/etc/xdg')
        dirs = [Path(d) / 'qxub' for d in xdg_config_dirs.split(':')]
        return dirs

    def _get_user_config_dir(self) -> Path:
        """Get user config directory following XDG spec."""
        xdg_config_home = os.environ.get('XDG_CONFIG_HOME')
        if xdg_config_home:
            return Path(xdg_config_home) / 'qxub'
        return Path.home() / '.config' / 'qxub'

    def _load_system_config(self) -> Optional[DictConfig]:
        """Load system-wide configuration."""
        for config_dir in self._get_xdg_config_dirs():
            config_file = config_dir / 'config.yaml'
            if config_file.exists():
                try:
                    return OmegaConf.load(config_file)
                except Exception as e:
                    click.echo(
                        f"Failed to load system config {config_file}: {e}",
                        err=True
                    )
        return None

    def _load_user_config(self) -> Optional[DictConfig]:
        """Load user configuration."""
        config_file = self._get_user_config_dir() / 'config.yaml'
        if config_file.exists():
            try:
                return OmegaConf.load(config_file)
            except Exception as e:  # pylint: disable=broad-exception-caught
                click.echo(f"Warning: Failed to load user config {config_file}: {e}", err=True)
        return None

    def _load_configs(self):
        """Load and merge all configuration files."""
        self.system_config = self._load_system_config()
        self.user_config = self._load_user_config()

        # Merge configurations: system < user
        configs = []
        if self.system_config:
            configs.append(self.system_config)
        if self.user_config:
            configs.append(self.user_config)

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
            files[f'system_{i}'] = config_dir / 'config.yaml'

        # User config file
        files['user'] = self._get_user_config_dir() / 'config.yaml'

        return files

    def get_template_variables(self, name: Optional[str] = None,
                             project: Optional[str] = None,
                             queue: Optional[str] = None) -> Dict[str, str]:
        """Get template variables for string substitution."""
        now = datetime.now()
        user = pwd.getpwuid(os.getuid()).pw_name

        return {
            'user': user,
            'project': project or '',
            'name': name or '',
            'queue': queue or '',
            'timestamp': now.strftime('%Y%m%d_%H%M%S'),
            'date': now.strftime('%Y%m%d'),
            'time': now.strftime('%H%M%S'),
        }

    def resolve_templates(self, value: Any, template_vars: Dict[str, str]) -> Any:
        """Recursively resolve template variables in configuration values."""
        if isinstance(value, str):
            # Only resolve if there are template placeholders
            if '{' in value and '}' in value:
                try:
                    return value.format(**template_vars)
                except KeyError as e:
                    click.echo(f"Warning: Unknown template variable {e} in '{value}'", err=True)
                    return value
            return value
        elif isinstance(value, (list, tuple)):
            return [self.resolve_templates(item, template_vars) for item in value]
        elif isinstance(value, dict):
            return {k: self.resolve_templates(v, template_vars) for k, v in value.items()}
        elif OmegaConf.is_config(value):
            # Handle OmegaConf DictConfig
            resolved = {}
            for k, v in value.items():
                resolved[k] = self.resolve_templates(v, template_vars)
            return OmegaConf.create(resolved)
        else:
            return value

    def get_defaults(self) -> Dict[str, Any]:
        """Get default configuration values."""
        if not self.merged_config or 'defaults' not in self.merged_config:
            return {}
        return OmegaConf.to_container(self.merged_config.defaults, resolve=True)

    def get_alias(self, alias_name: str) -> Optional[Dict[str, Any]]:
        """Get alias definition."""
        if (not self.merged_config or
            'aliases' not in self.merged_config or
            alias_name not in self.merged_config.aliases):
            return None
        return OmegaConf.to_container(self.merged_config.aliases[alias_name], resolve=True)

    def list_aliases(self) -> List[str]:
        """List all available aliases."""
        if not self.merged_config or 'aliases' not in self.merged_config:
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
        user_config_file = self._get_user_config_dir() / 'config.yaml'

        # Create directory if it doesn't exist
        user_config_file.parent.mkdir(parents=True, exist_ok=True)

        # Load existing user config or create new
        if user_config_file.exists():
            user_config = OmegaConf.load(user_config_file)
        else:
            user_config = OmegaConf.create({})

        # Set the value using OmegaConf.update
        keys = key_path.split('.')
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

    def create_user_config_template(self):
        """Create a template user config file."""
        user_config_file = self._get_user_config_dir() / 'config.yaml'

        if user_config_file.exists():
            raise click.ClickException(f"User config file already exists: {user_config_file}")

        # Create directory
        user_config_file.parent.mkdir(parents=True, exist_ok=True)

        # Create template config
        template = {
            'defaults': {
                'name': 'qt',
                'queue': 'normal',
                'project': 'a56',
                'joblog': '{name}.log',
                'resources': ['mem=4GB', 'ncpus=1'],
                'out': '/scratch/{project}/{user}/qt/{timestamp}/out',
                'err': '/scratch/{project}/{user}/qt/{timestamp}/err',
                'conda': {
                    'env': 'base',
                    'pre': None,
                    'post': None
                },
                'module': {
                    'mod': ['python3'],
                    'pre': None,
                    'post': None
                },
                'sing': {
                    'sif': None,
                    'bind': ['/scratch', '/g/data'],
                    'env': [],
                    'pre': None,
                    'post': None
                }
            },
            'aliases': {
                'example': {
                    'subcommand': 'conda',
                    'cmd': 'echo "Hello from qxub alias!"',
                    'name': 'example',
                    'conda': {
                        'env': 'base'
                    }
                }
            }
        }

        config = OmegaConf.create(template)
        OmegaConf.save(config, user_config_file)
        return user_config_file

    def resolve_options(self, cli_args: Dict[str, Any],
                       alias_name: Optional[str] = None) -> Dict[str, Any]:
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
                if key == 'subcommand':
                    # Don't merge subcommand into general options
                    continue
                elif key in ['conda', 'module', 'sing']:
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
            name=resolved.get('name'),
            project=resolved.get('project'),
            queue=resolved.get('queue')
        )
        resolved = self.resolve_templates(resolved, template_vars)

        return resolved


# Global config manager instance
config_manager = ConfigManager()
