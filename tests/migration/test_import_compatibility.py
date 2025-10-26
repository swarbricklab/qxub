"""
Test import compatibility during migration.
This test should pass at every phase of migration.
"""

import importlib

import pytest


class TestImportCompatibility:
    """Test that all public APIs remain accessible."""

    def test_core_module_imports(self):
        """Test core module imports work."""
        # These should always work
        import qxub
        from qxub import __version__

        assert qxub is not None
        assert __version__ is not None

    def test_config_imports(self):
        """Test configuration imports work."""
        # Current imports that should always work
        from qxub.config_manager import config_manager

        assert config_manager is not None

    def test_platform_imports(self):
        """Test platform imports work."""
        # Current imports that should always work
        from qxub.platform import get_platform, list_platforms

        assert callable(get_platform)
        assert callable(list_platforms)

    def test_scheduling_imports(self):
        """Test scheduling imports work."""
        # Current imports that should always work
        from qxub.scheduler import job_status, qdel, qsub

        assert callable(qsub)
        assert callable(qdel)
        assert callable(job_status)

    def test_resource_imports(self):
        """Test resource imports work."""
        # Current imports that should always work
        from qxub.resource_utils import (
            format_walltime,
            parse_memory_size,
            parse_walltime,
        )

        assert callable(parse_memory_size)
        assert callable(format_walltime)
        assert callable(parse_walltime)

    def test_execution_imports(self):
        """Test execution imports work."""
        # Current imports that should always work
        from qxub.execution_context import ExecutionContext, execute_unified

        assert callable(execute_unified)
        assert ExecutionContext is not None

    def test_cli_imports(self):
        """Test CLI imports work."""
        # These imports should work for CLI functionality
        from qxub.cli import qxub as main_cli
        from qxub.exec_cli import exec_cli

        assert main_cli is not None
        assert exec_cli is not None

    def test_optional_imports(self):
        """Test optional imports work or fail gracefully."""
        try:
            from qxub.remote_config import RemoteConfig

            assert RemoteConfig is not None
        except ImportError:
            # Remote functionality may not be available
            pass

    def test_import_performance(self):
        """Test that imports are reasonably fast."""
        import time

        start_time = time.time()
        importlib.import_module("qxub")
        import_time = time.time() - start_time

        # Should import in under 0.5 seconds
        assert import_time < 0.5, f"Import too slow: {import_time:.3f}s"


class TestCircularDependencies:
    """Test for circular dependencies."""

    def test_no_obvious_circular_imports(self):
        """Test that basic imports don't have circular dependencies."""
        # These imports should not cause circular import errors
        import qxub.config_manager
        import qxub.execution
        import qxub.platform
        import qxub.resource_utils
        import qxub.scheduler

        # If we get here without ImportError, basic imports work


if __name__ == "__main__":
    # Run as standalone script
    pytest.main([__file__, "-v"])
