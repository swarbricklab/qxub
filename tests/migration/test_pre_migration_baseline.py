"""
Comprehensive baseline test suite to run before migration begins.
Creates a reference point for all functionality.
"""

import subprocess
import sys
import time
from pathlib import Path

import pytest


class TestPreMigrationBaseline:
    """Baseline tests to establish reference behavior."""

    def test_all_cli_commands_exist(self):
        """Verify all CLI commands are accessible."""
        commands = [
            "qxub --help",
            "qxub config --help",
            "qxub exec --help",
            "qxub platform --help",
            "qxub history --help",
            "qxub monitor --help",
            "qxub cancel --help",
        ]
        for cmd in commands:
            result = subprocess.run(
                cmd.split(), capture_output=True, text=True, timeout=30
            )
            assert (
                result.returncode == 0
            ), f"Command failed: {cmd}\nError: {result.stderr}"

    def test_import_performance_baseline(self):
        """Measure current import performance as baseline."""
        start_time = time.time()
        import qxub

        import_time = time.time() - start_time

        # Store baseline (should be < 1.0 seconds on NCI)
        assert import_time < 1.0, f"Import too slow: {import_time}s"

        # Save baseline for comparison
        baseline_file = Path("tests/migration/.baseline_import_time")
        baseline_file.parent.mkdir(exist_ok=True)
        with open(baseline_file, "w") as f:
            f.write(str(import_time))

        print(f"📊 Baseline import time: {import_time:.3f}s")

    def test_all_imports_work(self):
        """Verify all current imports work."""
        import_tests = [
            "from qxub.config_manager import config_manager",
            "from qxub.platform import get_platform",
            "from qxub.core.scheduler import qsub, qdel",
            "from qxub.resources import parse_memory_size",
            "from qxub.execution_context import execute_unified",
        ]

        for import_test in import_tests:
            try:
                exec(import_test)
            except ImportError as e:
                pytest.fail(f"Import failed: {import_test}\nError: {e}")

    def test_dry_run_execution(self):
        """Test basic execution functionality."""
        cmd = ["qxub", "exec", "--dry-run", "--", "echo", "test"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Dry run failed: {result.stderr}"
        assert "qsub" in result.stdout, "Expected qsub command in dry run output"

    def test_configuration_access(self):
        """Test configuration system."""
        cmd = ["qxub", "config", "get", "defaults.project"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        # Should not crash, regardless of output
        assert result.returncode in [0, 1], f"Config command failed: {result.stderr}"

    def test_platform_detection(self):
        """Test platform detection."""
        cmd = ["qxub", "platform", "list"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Platform list failed: {result.stderr}"

    def test_conda_execution_context(self):
        """Test conda execution context works."""
        cmd = [
            "qxub",
            "exec",
            "--dry-run",
            "--env",
            "base",
            "--",
            "python",
            "--version",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Conda execution failed: {result.stderr}"
        assert "base" in result.stdout, "Expected conda environment in output"

    def test_module_execution_context(self):
        """Test module execution context works."""
        cmd = [
            "qxub",
            "exec",
            "--dry-run",
            "--mod",
            "python3",
            "--",
            "python",
            "--version",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Module execution failed: {result.stderr}"
        assert "python3" in result.stdout, "Expected module in output"

    def test_workflow_resource_options(self):
        """Test workflow-friendly resource options work."""
        cmd = [
            "qxub",
            "exec",
            "--dry",
            "--mem",
            "4GB",
            "--runtime",
            "2h",
            "--cpus",
            "2",
            "--",
            "echo",
            "test",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, f"Workflow resources failed: {result.stderr}"

    def test_resource_parsing(self):
        """Test resource parsing utilities work."""
        from qxub.resources import parse_memory_size, parse_walltime

        # Test memory parsing
        assert parse_memory_size("4GB") > 0
        assert parse_memory_size("1000MB") > 0

        # Test walltime parsing (returns seconds as int or None)
        walltime_seconds = parse_walltime("2:00:00")
        assert walltime_seconds == 7200  # 2 hours in seconds

    def test_platform_queue_selection(self):
        """Test platform queue selection works."""
        cmd = [
            "qxub",
            "select-queue",
            "--cpus",
            "2",
            "--memory",
            "4GB",
            "--format",
            "json",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        # May succeed or fail depending on platform, but should not crash
        assert result.returncode in [0, 1], f"Queue selection crashed: {result.stderr}"


if __name__ == "__main__":
    # Run as standalone script
    pytest.main([__file__, "-v"])
