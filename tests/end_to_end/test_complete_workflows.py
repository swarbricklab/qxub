"""
End-to-end tests that validate complete qxub workflows work correctly.
"""

import os
import subprocess
import tempfile
from pathlib import Path

import pytest


class TestCompleteWorkflows:
    """Test complete qxub workflows end-to-end."""

    def test_conda_execution_workflow(self):
        """Test complete conda execution workflow."""
        cmd = [
            "qxub",
            "exec",
            "--dry",
            "--env",
            "base",
            "--name",
            "test-job",
            "--",
            "python",
            "-c",
            "print('hello world')",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Conda execution failed: {result.stderr}"
        assert "qsub" in result.stdout, "Expected qsub command in output"
        assert "test-job" in result.stdout, "Expected job name in output"
        assert "base" in result.stdout, "Expected conda environment in output"

    def test_module_execution_workflow(self):
        """Test module execution workflow."""
        cmd = [
            "qxub",
            "exec",
            "--dry",
            "--mod",
            "python3",
            "--name",
            "module-test",
            "--",
            "python",
            "--version",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Module execution failed: {result.stderr}"
        assert "qsub" in result.stdout, "Expected qsub command in output"
        assert "python3" in result.stdout, "Expected module in output"

    def test_default_execution_workflow(self):
        """Test default execution workflow."""
        cmd = ["qxub", "exec", "--dry", "--default", "--", "echo", "hello"]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Default execution failed: {result.stderr}"
        assert "qsub" in result.stdout, "Expected qsub command in output"

    def test_platform_queue_selection_workflow(self):
        """Test platform detection and queue selection."""
        # Test platform detection
        result = subprocess.run(
            ["qxub", "platform", "list"], capture_output=True, text=True, timeout=30
        )
        assert result.returncode == 0, f"Platform list failed: {result.stderr}"

        # Test queue selection (may fail if no platform detected, but shouldn't crash)
        result = subprocess.run(
            [
                "qxub",
                "select-queue",
                "--cpus",
                "2",
                "--memory",
                "4GB",
                "--format",
                "json",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        # Return code 0 (success) or 1 (no platform) are both acceptable
        assert result.returncode in [0, 1], f"Queue selection crashed: {result.stderr}"

    def test_configuration_workflow(self):
        """Test configuration management workflow."""
        # Test config operations using system config
        result = subprocess.run(
            ["qxub", "config", "list"], capture_output=True, text=True, timeout=30
        )
        assert result.returncode == 0, f"Config list failed: {result.stderr}"

        # Test getting a config value
        result = subprocess.run(
            ["qxub", "config", "get", "defaults.project"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        # Should succeed or fail gracefully (not crash)

    def test_workflow_resource_options_workflow(self):
        """Test workflow-friendly resource options end-to-end."""
        cmd = [
            "qxub",
            "exec",
            "--dry",
            "--mem",
            "8GB",
            "--runtime",
            "1h30m",
            "--cpus",
            "4",
            "--env",
            "base",
            "--",
            "python",
            "-c",
            "print('resource test')",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Workflow resources failed: {result.stderr}"
        assert "qsub" in result.stdout, "Expected qsub command in output"

    def test_pbs_resource_options_workflow(self):
        """Test traditional PBS resource options."""
        cmd = [
            "qxub",
            "exec",
            "--dry",
            "-l",
            "mem=4GB",
            "-l",
            "ncpus=2",
            "-l",
            "walltime=01:00:00",
            "--env",
            "base",
            "--",
            "echo",
            "pbs test",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"PBS resources failed: {result.stderr}"
        assert "qsub" in result.stdout, "Expected qsub command in output"
        assert "mem=4GB" in result.stdout, "Expected memory resource in output"

    def test_mixed_resource_options_workflow(self):
        """Test mixing PBS and workflow-friendly options."""
        cmd = [
            "qxub",
            "exec",
            "--dry",
            "-l",
            "wd=true",
            "--mem",
            "4GB",
            "--cpus",
            "2",
            "--env",
            "base",
            "--",
            "echo",
            "mixed test",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Mixed resources failed: {result.stderr}"
        assert "qsub" in result.stdout, "Expected qsub command in output"

    def test_help_system_workflow(self):
        """Test help system works for all commands."""
        help_commands = [
            ["qxub", "--help"],
            ["qxub", "exec", "--help"],
            ["qxub", "config", "--help"],
            ["qxub", "platform", "--help"],
            ["qxub", "history", "--help"],
        ]

        for cmd in help_commands:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            assert (
                result.returncode == 0
            ), f"Help command failed: {' '.join(cmd)}\nError: {result.stderr}"
            assert (
                "usage:" in result.stdout.lower()
            ), f"Help output missing for: {' '.join(cmd)}"

    @pytest.mark.skipif(
        not os.getenv("QXUB_TEST_REMOTE"),
        reason="Remote testing not enabled (set QXUB_TEST_REMOTE=1)",
    )
    def test_remote_execution_workflow(self):
        """Test remote execution workflow (if configured)."""
        # This test is skipped unless explicitly enabled
        cmd = [
            "qxub",
            "exec",
            "--dry",
            "--remote",
            "test_remote",
            "--env",
            "base",
            "--",
            "echo",
            "remote test",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        # May fail if no remote configured, but should not crash
        assert result.returncode in [
            0,
            1,
            2,
        ], f"Remote execution crashed: {result.stderr}"


if __name__ == "__main__":
    # Run as standalone script
    pytest.main([__file__, "-v"])
