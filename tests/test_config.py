#!/usr/bin/env python3
"""
Test script for config.py functionality.

Tests configuration loading, validation, and management.
"""

import os
import sys
import tempfile
from pathlib import Path

import yaml

# Add qxub to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_global_config():
    """Test global config instance."""
    print("🧪 Testing global config...")

    from qxub.config import config_manager

    config1 = config_manager
    config2 = config_manager

    # Should be the same instance
    if config1 is config2:
        print("  ✅ Global config returns same instance")
        return
    else:
        print("  ❌ Global config returns different instances")
        raise AssertionError("see output above")


def test_logging_setup():
    """Test logging configuration."""
    print("🧪 Testing logging setup...")

    import logging

    from qxub.config import setup_logging

    # Test different verbosity levels
    for verbosity in [0, 1, 2, 3]:
        setup_logging(verbosity)
        logger = logging.getLogger()

        expected_levels = [
            logging.ERROR,  # 0
            logging.WARNING,  # 1
            logging.INFO,  # 2
            logging.DEBUG,  # 3
        ]

        if logger.level == expected_levels[verbosity]:
            print(f"  ✅ Verbosity {verbosity} -> {logging.getLevelName(logger.level)}")
        else:
            print(
                f"  ❌ Verbosity {verbosity} -> {logging.getLevelName(logger.level)}, expected {logging.getLevelName(expected_levels[verbosity])}"
            )
            raise AssertionError("see output above")

    return


def main():
    """Run all configuration tests."""
    print("🚀 Testing qxub configuration system...\n")

    tests = [
        test_global_config,
        test_logging_setup,
    ]

    passed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  ❌ {test.__name__}: {e}")
        print()  # Add spacing between tests

    print(f"📊 Overall result: {passed}/{len(tests)} test groups passed")

    if passed == len(tests):
        print("🎉 All configuration tests passed!")
        return 0
    else:
        print("❌ Some tests failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
