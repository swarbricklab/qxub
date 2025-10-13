"""
Tests for remote execution functionality in qxub v2.2.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest
import yaml

from qxub.remote_config import RemoteConfig, RemoteConfigError, UnsupportedProtocolError
from qxub.remote_config_loader import (
    ConfigLoadError,
    get_remote_config,
    load_remote_configurations,
)
from qxub.remote_executor import (
    ConnectionError,
    RemoteExecutorFactory,
    SSHRemoteExecutor,
)


class TestRemoteConfig:
    """Test RemoteConfig class functionality."""

    def test_valid_ssh_config(self):
        """Test creating a valid SSH remote configuration."""
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
        )

        assert config.name == "test_remote"
        assert config.protocol == "ssh"
        assert config.hostname == "test.example.com"
        assert config.port is None
        assert config.qxub_env == "test-env"
        assert config.platform_file == "/path/to/platform.yaml"

    def test_ssh_config_with_port(self):
        """Test SSH configuration with custom port."""
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com:2222",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
        )

        assert config.hostname == "test.example.com"
        assert config.port == 2222

    def test_ssh_config_with_username(self):
        """Test SSH configuration with username in URL."""
        config = RemoteConfig(
            name="test_remote",
            url="ssh://user@test.example.com",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
        )

        assert config.hostname == "test.example.com"
        assert config.username == "user"

    def test_invalid_url_no_protocol(self):
        """Test that missing protocol raises error."""
        with pytest.raises(RemoteConfigError, match="URL must include protocol"):
            RemoteConfig(
                name="test_remote",
                url="test.example.com",
                qxub_env="test-env",
                platform_file="/path/to/platform.yaml",
            )

    def test_invalid_url_no_hostname(self):
        """Test that missing hostname raises error."""
        with pytest.raises(RemoteConfigError, match="URL must include hostname"):
            RemoteConfig(
                name="test_remote",
                url="ssh://",
                qxub_env="test-env",
                platform_file="/path/to/platform.yaml",
            )

    def test_unsupported_protocol(self):
        """Test that unsupported protocol raises error."""
        with pytest.raises(
            UnsupportedProtocolError, match="Protocol 'ftp' not supported"
        ):
            RemoteConfig(
                name="test_remote",
                url="ftp://test.example.com",
                qxub_env="test-env",
                platform_file="/path/to/platform.yaml",
            )

    def test_environment_variable_substitution(self):
        """Test environment variable substitution in project_root_dir."""
        with patch.dict(os.environ, {"USER": "testuser", "PROJECT": "testproject"}):
            config = RemoteConfig(
                name="test_remote",
                url="ssh://test.example.com",
                qxub_env="test-env",
                platform_file="/path/to/platform.yaml",
                project_root_dir="/scratch/${PROJECT}/${USER}/projects",
            )

            result = config.get_project_root_dir()
            assert result == "/scratch/testproject/testuser/projects"

    def test_working_directory_resolution(self):
        """Test remote working directory resolution logic."""
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
            project_root_dir="/remote/projects",
        )

        local_cwd = Path("/local/my-project")

        # Test smart default
        result = config.determine_remote_working_dir(local_cwd)
        assert result == "/remote/projects/my-project"

        # Test explicit execdir
        result = config.determine_remote_working_dir(local_cwd, "/explicit/path")
        assert result == "/explicit/path"

        # Test relative execdir
        result = config.determine_remote_working_dir(local_cwd, "subdir/experiment")
        assert result == "/remote/projects/subdir/experiment"

    def test_validation(self):
        """Test configuration validation."""
        # Valid config
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
        )
        errors = config.validate()
        assert len(errors) == 0

        # Invalid config - missing required fields
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com",
            qxub_env="",
            platform_file="",
        )
        errors = config.validate()
        assert "qxub_env is required" in errors
        assert "platform_file is required" in errors


class TestRemoteConfigLoader:
    """Test remote configuration loading functionality."""

    def test_load_valid_config(self):
        """Test loading a valid configuration file."""
        config_data = {
            "remotes": {
                "test_remote": {
                    "url": "ssh://test.example.com",
                    "qxub_env": "test-env",
                    "platform_file": "/path/to/platform.yaml",
                    "project_root_dir": "/scratch/${USER}/projects",
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name

        try:
            with patch(
                "qxub.remote_config_loader.get_user_config_path",
                return_value=Path(config_path),
            ):
                remotes = load_remote_configurations()

            assert "test_remote" in remotes
            config = remotes["test_remote"]
            assert config.name == "test_remote"
            assert config.url == "ssh://test.example.com"
            assert config.qxub_env == "test-env"

        finally:
            os.unlink(config_path)

    def test_load_missing_config_file(self):
        """Test handling of missing configuration file."""
        with patch(
            "qxub.remote_config_loader.get_user_config_path",
            return_value=Path("/nonexistent/config.yaml"),
        ):
            remotes = load_remote_configurations()
            assert remotes == {}

    def test_load_invalid_yaml(self):
        """Test handling of invalid YAML."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("invalid: yaml: content: [")
            config_path = f.name

        try:
            with patch(
                "qxub.remote_config_loader.get_user_config_path",
                return_value=Path(config_path),
            ):
                with pytest.raises(ConfigLoadError, match="Invalid YAML"):
                    load_remote_configurations()
        finally:
            os.unlink(config_path)

    def test_missing_required_fields(self):
        """Test handling of missing required fields."""
        config_data = {
            "remotes": {
                "test_remote": {
                    "url": "ssh://test.example.com",
                    # Missing qxub_env and platform_file
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name

        try:
            with patch(
                "qxub.remote_config_loader.get_user_config_path",
                return_value=Path(config_path),
            ):
                with pytest.raises(ConfigLoadError, match="missing required field"):
                    load_remote_configurations()
        finally:
            os.unlink(config_path)

    def test_get_remote_config(self):
        """Test getting a specific remote configuration."""
        config_data = {
            "remotes": {
                "test_remote": {
                    "url": "ssh://test.example.com",
                    "qxub_env": "test-env",
                    "platform_file": "/path/to/platform.yaml",
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name

        try:
            with patch(
                "qxub.remote_config_loader.get_user_config_path",
                return_value=Path(config_path),
            ):
                config = get_remote_config("test_remote")
                assert config.name == "test_remote"

                # Test non-existent remote
                with pytest.raises(
                    ConfigLoadError, match="Remote 'nonexistent' not found"
                ):
                    get_remote_config("nonexistent")
        finally:
            os.unlink(config_path)


class TestSSHRemoteExecutor:
    """Test SSH remote executor functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
        )
        self.executor = SSHRemoteExecutor(self.config)

    def test_build_ssh_command(self):
        """Test SSH command building."""
        command = "qxub --env pytorch -- python script.py"
        working_dir = "/remote/working/dir"

        ssh_cmd = self.executor._build_ssh_command(command, working_dir)

        assert ssh_cmd[0] == "ssh"
        assert "test.example.com" in ssh_cmd

        # Check that the remote command includes our components
        remote_command = ssh_cmd[-1]
        assert f"cd {working_dir}" in remote_command
        assert f"conda activate {self.config.qxub_env}" in remote_command
        assert command in remote_command

    def test_build_ssh_command_with_port(self):
        """Test SSH command building with custom port."""
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com:2222",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
        )
        executor = SSHRemoteExecutor(config)

        ssh_cmd = executor._build_ssh_command("test command", "/tmp")

        assert "-p" in ssh_cmd
        port_index = ssh_cmd.index("-p")
        assert ssh_cmd[port_index + 1] == "2222"

    def test_build_ssh_command_with_config_file(self):
        """Test SSH command building with config file."""
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
            config="/custom/ssh/config",
        )

        with patch("pathlib.Path.exists", return_value=True):
            executor = SSHRemoteExecutor(config)
            ssh_cmd = executor._build_ssh_command("test command", "/tmp")

            assert "-F" in ssh_cmd
            config_index = ssh_cmd.index("-F")
            assert ssh_cmd[config_index + 1] == "/custom/ssh/config"

    @patch("subprocess.run")
    def test_test_connection_success(self, mock_run):
        """Test successful connection test."""
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = "connection_test"
        mock_run.return_value = mock_result

        result = self.executor.test_connection()
        assert result is True

        # Verify the SSH command was called correctly
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert call_args[0] == "ssh"
        assert "test.example.com" in call_args
        assert "echo connection_test" in call_args

    @patch("subprocess.run")
    def test_test_connection_failure(self, mock_run):
        """Test failed connection test."""
        mock_result = Mock()
        mock_result.returncode = 1
        mock_result.stdout = ""
        mock_run.return_value = mock_result

        result = self.executor.test_connection()
        assert result is False

    @patch("subprocess.run")
    def test_test_connection_timeout(self, mock_run):
        """Test connection test timeout."""
        mock_run.side_effect = subprocess.TimeoutExpired("ssh", 15)

        result = self.executor.test_connection()
        assert result is False


class TestRemoteExecutorFactory:
    """Test remote executor factory."""

    def test_create_ssh_executor(self):
        """Test creating SSH executor."""
        config = RemoteConfig(
            name="test_remote",
            url="ssh://test.example.com",
            qxub_env="test-env",
            platform_file="/path/to/platform.yaml",
        )

        executor = RemoteExecutorFactory.create(config)
        assert isinstance(executor, SSHRemoteExecutor)

    def test_create_unsupported_executor(self):
        """Test creating executor for unsupported protocol."""
        # We need to bypass the RemoteConfig validation for this test
        config = Mock()
        config.protocol = "unsupported"

        with pytest.raises(
            UnsupportedProtocolError, match="Protocol 'unsupported' not supported"
        ):
            RemoteExecutorFactory.create(config)


if __name__ == "__main__":
    pytest.main([__file__])
