"""
Performance regression tests for migration.

Each limit is checked against the fastest of a few runs (``timed_run``), so
a busy shared login node doesn't fail them; a real slowdown still does.
"""

import pytest
from conftest import timed_run


class TestPerformanceRegression:
    """Test performance hasn't degraded after migration."""

    def test_import_performance(self, cold_import_time):
        """A cold import should take under a second."""
        print(f"📊 Cold import time: {cold_import_time:.3f}s")
        assert cold_import_time < 1.0, f"Import too slow: {cold_import_time:.3f}s"

    def test_cli_startup_time(self):
        """Test CLI startup time hasn't degraded."""
        result, startup_time = timed_run(["qxub", "--help"], limit=3.0, timeout=10)

        print(f"📊 CLI startup time: {startup_time:.3f}s")

        assert result.returncode == 0, "CLI help command failed"
        assert startup_time < 3.0, f"CLI startup too slow: {startup_time:.3f}s"

    def test_dry_run_performance(self):
        """Test dry run performance."""
        result, dry_run_time = timed_run(
            ["qxub", "exec", "--dry", "--env", "base", "--", "echo", "test"],
            limit=15.0,
            timeout=15,
        )

        print(f"📊 Dry run time: {dry_run_time:.3f}s")

        assert result.returncode == 0, "Dry run failed"
        assert dry_run_time < 15.0, f"Dry run too slow: {dry_run_time:.3f}s"

    def test_config_access_performance(self):
        """Test configuration access performance."""
        result, config_time = timed_run(
            ["qxub", "config", "get", "defaults.project"], limit=3.0, timeout=10
        )

        print(f"📊 Config access time: {config_time:.3f}s")

        # Should not crash, performance should be reasonable
        assert result.returncode in [0, 1], "Config command crashed"
        assert config_time < 3.0, f"Config access too slow: {config_time:.3f}s"

    def test_platform_list_performance(self):
        """Test platform listing performance."""
        result, platform_time = timed_run(
            ["qxub", "platform", "list"], limit=3.0, timeout=10
        )

        print(f"📊 Platform list time: {platform_time:.3f}s")

        assert result.returncode == 0, "Platform list failed"
        assert platform_time < 3.0, f"Platform list too slow: {platform_time:.3f}s"


if __name__ == "__main__":
    # Run as standalone script
    pytest.main([__file__, "-v"])
