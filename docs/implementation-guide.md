# URL-Based Remote Configuration Implementation Guide

## Overview

This document outlines the implementation approach for URL-based remote configuration in qxub v2.2, designed to be future-proof for multiple execution protocols while starting with SSH-only support.

## Configuration Schema

### Current Implementation (SSH-only)

```python
from dataclasses import dataclass
from typing import Optional
from urllib.parse import urlparse

@dataclass
class RemoteConfig:
    """Configuration for remote execution endpoints."""

    name: str
    url: str  # protocol://host:port format
    config: Optional[str] = None  # Protocol-specific config file
    qxub_env: str = "qxub"
    platform_file: str = ""
    project_root_dir: Optional[str] = None

    def __post_init__(self):
        """Validate and set smart defaults."""
        self.parsed_url = urlparse(self.url)

        if not self.parsed_url.scheme:
            raise ConfigError(f"URL must include protocol: {self.url}")

        if self.parsed_url.scheme not in ['ssh']:
            raise ConfigError(f"Unsupported protocol: {self.parsed_url.scheme}")

        # Set smart defaults for config file
        if self.config is None:
            self.config = self._get_default_config_path()

    def _get_default_config_path(self) -> str:
        """Get default config file path based on protocol."""
        defaults = {
            'ssh': '~/.ssh/config',
            'k8s': '~/.kube/config',  # Future
            'aws': '~/.aws/config',   # Future
        }
        return defaults.get(self.parsed_url.scheme, '')

    @property
    def protocol(self) -> str:
        """Get connection protocol."""
        return self.parsed_url.scheme

    @property
    def hostname(self) -> str:
        """Get hostname from URL."""
        return self.parsed_url.hostname or ''

    @property
    def port(self) -> Optional[int]:
        """Get port from URL."""
        return self.parsed_url.port
```

## Protocol-Specific Executors

### Executor Factory Pattern

```python
from abc import ABC, abstractmethod

class RemoteExecutor(ABC):
    """Abstract base class for remote execution backends."""

    def __init__(self, config: RemoteConfig):
        self.config = config

    @abstractmethod
    def execute(self, command: str, working_dir: str) -> int:
        """Execute command on remote system."""
        pass

    @abstractmethod
    def test_connection(self) -> bool:
        """Test if connection is available."""
        pass

class SSHRemoteExecutor(RemoteExecutor):
    """SSH-based remote execution."""

    def execute(self, command: str, working_dir: str) -> int:
        """Execute command via SSH."""
        ssh_command = self._build_ssh_command(command, working_dir)
        return subprocess.run(ssh_command).returncode

    def _build_ssh_command(self, command: str, working_dir: str) -> list[str]:
        """Build SSH command with proper options."""
        ssh_cmd = ['ssh']

        # Add config file if specified
        if self.config.config:
            ssh_cmd.extend(['-F', self.config.config])

        # Add port if specified in URL
        if self.config.port:
            ssh_cmd.extend(['-p', str(self.config.port)])

        # Add hostname
        ssh_cmd.append(self.config.hostname)

        # Add remote command
        remote_cmd = f"cd {working_dir} && conda activate {self.config.qxub_env} && {command}"
        ssh_cmd.append(remote_cmd)

        return ssh_cmd

    def test_connection(self) -> bool:
        """Test SSH connection."""
        test_cmd = ['ssh', '-o', 'ConnectTimeout=10']

        if self.config.config:
            test_cmd.extend(['-F', self.config.config])

        test_cmd.extend([self.config.hostname, 'echo', 'connection_test'])

        try:
            result = subprocess.run(test_cmd, capture_output=True, timeout=15)
            return result.returncode == 0
        except subprocess.TimeoutExpired:
            return False

class RemoteExecutorFactory:
    """Factory for creating protocol-specific executors."""

    @staticmethod
    def create(config: RemoteConfig) -> RemoteExecutor:
        """Create appropriate executor for protocol."""
        executors = {
            'ssh': SSHRemoteExecutor,
            # Future protocols:
            # 'aws': AWSBatchExecutor,
            # 'k8s': KubernetesExecutor,
        }

        executor_class = executors.get(config.protocol)
        if not executor_class:
            supported = ', '.join(executors.keys())
            raise UnsupportedProtocolError(
                f"Protocol '{config.protocol}' not supported. "
                f"Supported protocols: {supported}"
            )

        return executor_class(config)
```

## Configuration Loading

### YAML Configuration Parser

```python
import yaml
from pathlib import Path
from typing import Dict

def load_remote_configurations() -> Dict[str, RemoteConfig]:
    """Load remote configurations from user config file."""

    config_path = Path.home() / '.config' / 'qxub' / 'config.yaml'

    if not config_path.exists():
        return {}

    with open(config_path) as f:
        config_data = yaml.safe_load(f)

    remotes = {}
    for name, remote_data in config_data.get('remotes', {}).items():
        try:
            remotes[name] = RemoteConfig(
                name=name,
                url=remote_data['url'],
                config=remote_data.get('config'),
                qxub_env=remote_data.get('qxub_env', 'qxub'),
                platform_file=remote_data['platform_file'],
                project_root_dir=remote_data.get('project_root_dir')
            )
        except KeyError as e:
            raise ConfigError(f"Missing required field in remote '{name}': {e}")
        except Exception as e:
            raise ConfigError(f"Error loading remote '{name}': {e}")

    return remotes
```

## CLI Integration

### Command Line Processing

```python
def execute_remote_command(remote_name: str, args: list[str]) -> int:
    """Execute command on remote system."""

    # Load remote configurations
    remotes = load_remote_configurations()

    if remote_name not in remotes:
        available = ', '.join(remotes.keys())
        raise RemoteNotFoundError(
            f"Remote '{remote_name}' not found. "
            f"Available remotes: {available}"
        )

    remote_config = remotes[remote_name]

    # Create appropriate executor
    executor = RemoteExecutorFactory.create(remote_config)

    # Test connection
    if not executor.test_connection():
        raise ConnectionError(
            f"Cannot connect to {remote_config.hostname}. "
            f"Check your {remote_config.protocol} configuration."
        )

    # Build remote command
    remote_args = [arg for arg in args if not arg.startswith('--remote')]
    remote_args = ['--platform-file', remote_config.platform_file] + remote_args
    remote_command = f"qxub {' '.join(remote_args)}"

    # Determine working directory
    working_dir = determine_remote_working_dir(remote_config)

    # Execute
    return executor.execute(remote_command, working_dir)
```

## Future Protocol Extensions

### AWS Batch Example (Future v2.3+)

```python
class AWSBatchExecutor(RemoteExecutor):
    """AWS Batch execution backend."""

    def execute(self, command: str, working_dir: str) -> int:
        """Submit job to AWS Batch."""
        # Implementation would use boto3 to submit batch jobs
        # with qxub container image
        pass

    def test_connection(self) -> bool:
        """Test AWS credentials and batch access."""
        # Implementation would test AWS API access
        pass
```

### Kubernetes Example (Future v2.3+)

```python
class KubernetesExecutor(RemoteExecutor):
    """Kubernetes Job execution backend."""

    def execute(self, command: str, working_dir: str) -> int:
        """Create Kubernetes Job."""
        # Implementation would use kubectl or kubernetes python client
        # to create Job resources with qxub container
        pass

    def test_connection(self) -> bool:
        """Test cluster access."""
        # Implementation would test cluster connectivity
        pass
```

## Error Handling

### User-Friendly Error Messages

```python
class RemoteExecutionError(Exception):
    """Base class for remote execution errors."""
    pass

class UnsupportedProtocolError(RemoteExecutionError):
    """Raised when protocol is not supported."""

    def __init__(self, message: str):
        super().__init__(message)
        self.suggestion = "Currently only 'ssh://' URLs are supported."

class ConnectionError(RemoteExecutionError):
    """Raised when connection fails."""

    def __init__(self, message: str, protocol: str = 'ssh'):
        super().__init__(message)
        self.suggestions = {
            'ssh': [
                "Check SSH configuration in ~/.ssh/config",
                "Verify network connectivity",
                "Test connection: ssh hostname echo 'test'"
            ]
        }

    def get_suggestions(self) -> list[str]:
        return self.suggestions.get(self.protocol, [])
```

## Benefits of This Architecture

1. **🔮 Future-Proof**: Easy to add new protocols without breaking changes
2. **🧹 Clean Separation**: Protocol logic isolated in dedicated executors
3. **⚙️ Smart Defaults**: Reasonable config file defaults for each protocol
4. **🛠️ Extensible**: Factory pattern makes adding protocols straightforward
5. **🔒 Delegated Auth**: Each protocol uses its standard authentication
6. **📝 Consistent Config**: Same configuration schema for all protocols

This design provides the foundation for multi-protocol remote execution while maintaining simplicity for the current SSH-only implementation.
