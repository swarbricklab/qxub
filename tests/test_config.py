#!/usr/bin/env python3
"""
Test script for config.py functionality.

Tests configuration loading, validation, and management.
"""

import sys
import os
import tempfile
import yaml
from pathlib import Path

# Add qxub to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_config_creation():
    """Test QxubConfig creation and basic functionality."""
    print("🧪 Testing config creation...")

    from qxub.config import QxubConfig

    config = QxubConfig()

    # Test default values
    default_resources = config.get_default_resources()
    queue_prefs = config.get_queue_preferences()

    if (
        default_resources.get("cpus") == 1
        and queue_prefs.get("optimization") == "balanced"
    ):
        print("  ✅ Config created with correct defaults")
        return True
    else:
        print("  ❌ Config defaults are incorrect")
        return False


def test_config_search_paths():
    """Test platform search path management."""
    print("🧪 Testing config search paths...")

    from qxub.config import QxubConfig

    config = QxubConfig()
    search_paths = config.get_platform_search_paths()

    if len(search_paths) >= 2 and all(isinstance(p, Path) for p in search_paths):
        print(f"  ✅ Found {len(search_paths)} search paths")
        for path in search_paths:
            print(f"      {path}")
        return True
    else:
        print("  ❌ Search paths not configured correctly")
        return False


def test_config_file_loading():
    """Test configuration file loading."""
    print("🧪 Testing config file loading...")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Create a test config file
        config_data = {
            "default_platform": "test_platform",
            "queue_preferences": {"optimization": "cost", "adjustment_policy": "auto"},
            "default_resources": {"cpus": 2, "memory": "8GB", "walltime": "2:00:00"},
            "verbosity": 2,
        }

        config_file = temp_path / "config.yaml"
        with open(config_file, "w") as f:
            yaml.dump(config_data, f)

        # Mock the config paths
        from qxub.config import QxubConfig

        # Create config with temporary home directory
        original_home = os.environ.get("HOME")
        try:
            os.environ["HOME"] = str(temp_path)

            # Create the expected config directory structure
            config_dir = temp_path / ".config" / "qxub"
            config_dir.mkdir(parents=True)
            (config_dir / "config.yaml").write_text(yaml.dump(config_data))

            config = QxubConfig()

            # Test loaded values
            if (
                config.get_default_platform() == "test_platform"
                and config.get("verbosity") == 2
                and config.get_default_resources()["cpus"] == 2
            ):
                print("  ✅ Config file loaded successfully")
                return True
            else:
                print("  ❌ Config file not loaded correctly")
                return False
        finally:
            if original_home:
                os.environ["HOME"] = original_home
            elif "HOME" in os.environ:
                del os.environ["HOME"]


def test_config_validation():
    """Test configuration validation."""
    print("🧪 Testing config validation...")

    from qxub.config import validate_config, QxubConfig

    # Test with default config
    config = QxubConfig()
    issues = validate_config()

    # We expect some issues since we don't have real platform files
    if isinstance(issues, list):
        print(f"  ✅ Validation returned {len(issues)} issues")
        for issue in issues[:3]:  # Show first 3 issues
            print(f"      {issue}")
        return True
    else:
        print("  ❌ Validation did not return a list")
        return False


def test_config_user_preferences():
    """Test setting and getting user preferences."""
    print("🧪 Testing user preferences...")

    from qxub.config import QxubConfig

    config = QxubConfig()

    # Set a preference
    config.set_user_preference("test_key", "test_value")

    # Get the preference
    value = config.get("test_key")

    if value == "test_value":
        print("  ✅ User preference set and retrieved correctly")
        return True
    else:
        print(f"  ❌ User preference not set correctly: {value}")
        return False


def test_effective_config():
    """Test getting effective configuration."""
    print("🧪 Testing effective config...")

    from qxub.config import get_effective_config

    effective_config = get_effective_config()

    required_keys = [
        "platform_search_paths",
        "default_platform",
        "platform_preferences",
        "queue_preferences",
        "default_resources",
        "verbosity",
    ]

    missing_keys = [key for key in required_keys if key not in effective_config]

    if not missing_keys:
        print("  ✅ Effective config has all required keys")
        return True
    else:
        print(f"  ❌ Effective config missing keys: {missing_keys}")
        return False


def test_global_config():
    """Test global config instance."""
    print("🧪 Testing global config...")

    from qxub.config import get_config

    config1 = get_config()
    config2 = get_config()

    # Should be the same instance
    if config1 is config2:
        print("  ✅ Global config returns same instance")
        return True
    else:
        print("  ❌ Global config returns different instances")
        return False


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
            return False

    return True


def main():
    """Run all configuration tests."""
    print("🚀 Testing qxub configuration system...\n")

    tests = [
        test_config_creation,
        test_config_search_paths,
        test_config_file_loading,
        test_config_validation,
        test_config_user_preferences,
        test_effective_config,
        test_global_config,
        test_logging_setup,
    ]

    passed = 0
    for test in tests:
        if test():
            passed += 1
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
