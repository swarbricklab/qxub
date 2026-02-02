"""
Performance regression tests for migration.
"""

import subprocess
import time
from pathlib import Path

import pytest


class TestPerformanceRegression:
    """Test performance hasn't degraded after migration."""

    def test_import_performance(self):
        """Test import performance hasn't degraded."""
        # Load baseline if available
        baseline_file = Path("tests/migration/.baseline_import_time")
        baseline = None
        if baseline_file.exists():
            try:
                with open(baseline_file) as f:
                    baseline = float(f.read().strip())
            except (FileNotFoundError, ValueError):
                pass

        # Measure current import time
        start_time = time.time()
        import qxub

        import_time = time.time() - start_time

        print(f"📊 Current import time: {import_time:.3f}s")

        if baseline:
            print(f"📊 Baseline import time: {baseline:.3f}s")
            # Allow 50% degradation during migration (temporary)
            max_allowed = baseline * 1.5
            assert (
                import_time <= max_allowed
            ), f"Import time degraded: {import_time:.3f}s > {max_allowed:.3f}s (baseline: {baseline:.3f}s)"
        else:
            # Without baseline, just ensure reasonable performance
            assert import_time < 1.0, f"Import too slow: {import_time:.3f}s"

    def test_cli_startup_time(self):
        """Test CLI startup time hasn't degraded."""
        start_time = time.time()
        result = subprocess.run(
            ["qxub", "--help"], capture_output=True, text=True, timeout=10
        )
        startup_time = time.time() - start_time

        print(f"📊 CLI startup time: {startup_time:.3f}s")

        assert result.returncode == 0, "CLI help command failed"
        assert startup_time < 3.0, f"CLI startup too slow: {startup_time:.3f}s"

    def test_dry_run_performance(self):
        """Test dry run performance."""
        start_time = time.time()
        result = subprocess.run(
            ["qxub", "exec", "--dry", "--env", "base", "--", "echo", "test"],
            capture_output=True,
            text=True,
            timeout=15,
        )
        dry_run_time = time.time() - start_time

        print(f"📊 Dry run time: {dry_run_time:.3f}s")

        assert result.returncode == 0, "Dry run failed"
        assert dry_run_time < 15.0, f"Dry run too slow: {dry_run_time:.3f}s"

    def test_config_access_performance(self):
        """Test configuration access performance."""
        start_time = time.time()
        result = subprocess.run(
            ["qxub", "config", "get", "defaults.project"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        config_time = time.time() - start_time

        print(f"📊 Config access time: {config_time:.3f}s")

        # Should not crash, performance should be reasonable
        assert result.returncode in [0, 1], "Config command crashed"
        assert config_time < 3.0, f"Config access too slow: {config_time:.3f}s"

    def test_platform_list_performance(self):
        """Test platform listing performance."""
        start_time = time.time()
        result = subprocess.run(
            ["qxub", "platform", "list"], capture_output=True, text=True, timeout=10
        )
        platform_time = time.time() - start_time

        print(f"📊 Platform list time: {platform_time:.3f}s")

        assert result.returncode == 0, "Platform list failed"
        assert platform_time < 3.0, f"Platform list too slow: {platform_time:.3f}s"


if __name__ == "__main__":
    # Run as standalone script
    pytest.main([__file__, "-v"])
