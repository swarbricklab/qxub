#!/usr/bin/env python3
"""
Test script for platform.py functionality.

Tests platform loading, queue validation, and queue selection.
"""

import os
import sys
import tempfile
from pathlib import Path

import yaml

# Add qxub to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from qxub.platform import (
    Platform,
    PlatformLoader,
    Queue,
    QueueLimits,
    QueueSelector,
    ResourceValidationResult,
    WalltimeRule,
    select_best_queue,
)


def create_test_platform_file(temp_dir: Path) -> Path:
    """Create a test platform definition file."""
    platform_data = {
        "name": "test_platform",
        "type": "pbs",
        "host": "test.example.com",
        "description": "Test platform for unit testing",
        "queues": [
            {
                "name": "normal",
                "type": "compute",
                "priority": "normal",
                "su_billing_rate": 1.0,
                "limits": {
                    "max_cpus": 48,
                    "min_cpus": 1,
                    "max_memory": "192GB",
                    "min_memory": "1GB",
                },
                "walltime_rules": [{"cores": "1-48", "max_walltime": "48:00:00"}],
                "default_walltime": "1:00:00",
            },
            {
                "name": "express",
                "type": "compute",
                "priority": "high",
                "su_billing_rate": 2.0,
                "limits": {"max_cpus": 12, "min_cpus": 1, "max_memory": "48GB"},
                "walltime_rules": [{"cores": "1-12", "max_walltime": "2:00:00"}],
            },
            {
                "name": "gpu",
                "type": "gpu",
                "priority": "normal",
                "su_billing_rate": 3.0,
                "limits": {
                    "max_cpus": 12,
                    "min_cpus": 1,
                    "max_gpus": 4,
                    "min_gpus": 1,
                    "max_memory": "96GB",
                },
                "walltime_rules": [{"cores": "1-12", "max_walltime": "24:00:00"}],
            },
        ],
        "auto_selection_rules": [
            {"condition": "gpu_requested > 0", "queue": "gpu"},
            {"condition": "walltime <= 2:00:00", "queue": "express"},
            {"condition": "", "queue": "normal", "is_default": True},
        ],
    }

    platform_file = temp_dir / "test_platform.yaml"
    with open(platform_file, "w") as f:
        yaml.dump(platform_data, f, default_flow_style=False)

    return platform_file


def test_queue_limits():
    """Test QueueLimits functionality."""
    print("🧪 Testing QueueLimits...")

    limits = QueueLimits(max_cpus=48, min_cpus=1, max_memory="192GB", min_memory="1GB")

    if limits.max_cpus == 48 and limits.max_memory == "192GB":
        print("  ✅ QueueLimits creation successful")
        return True
    else:
        print("  ❌ QueueLimits creation failed")
        return False


def test_walltime_rule():
    """Test WalltimeRule functionality."""
    print("🧪 Testing WalltimeRule...")

    rule = WalltimeRule(cores="1-48", max_walltime="48:00:00")

    tests = [
        (24, True),  # Within range
        (48, True),  # At upper bound
        (1, True),  # At lower bound
        (96, False),  # Above range
    ]

    passed = 0
    for core_count, expected in tests:
        result = rule.matches_core_count(core_count)
        if result == expected:
            print(f"  ✅ {core_count} cores -> {result}")
            passed += 1
        else:
            print(f"  ❌ {core_count} cores -> {result}, expected {expected}")

    print(f"WalltimeRule: {passed}/{len(tests)} tests passed\n")
    return passed == len(tests)


def test_queue_validation():
    """Test Queue resource validation."""
    print("🧪 Testing Queue validation...")

    limits = QueueLimits(max_cpus=48, min_cpus=1, max_memory="192GB")
    walltime_rules = [WalltimeRule(cores="1-48", max_walltime="48:00:00")]

    queue = Queue(
        name="test_queue", type="compute", limits=limits, walltime_rules=walltime_rules
    )

    test_cases = [
        # Valid resources
        ({"cpus": 4, "memory": "8GB", "walltime": "2:00:00"}, True),
        ({"cpus": 1, "memory": "1GB", "walltime": "1:00:00"}, True),
        # Invalid resources
        ({"cpus": 96, "memory": "8GB", "walltime": "2:00:00"}, False),  # Too many CPUs
        (
            {"cpus": 4, "memory": "256GB", "walltime": "2:00:00"},
            False,
        ),  # Too much memory
        (
            {"cpus": 4, "memory": "8GB", "walltime": "72:00:00"},
            False,
        ),  # Too much walltime
    ]

    passed = 0
    for resources, should_be_valid in test_cases:
        result = queue.validate_resources(resources)
        if result.is_valid == should_be_valid:
            print(f"  ✅ {resources} -> valid={result.is_valid}")
            if not result.is_valid and result.errors:
                print(f"      Errors: {result.errors}")
            passed += 1
        else:
            print(
                f"  ❌ {resources} -> valid={result.is_valid}, expected {should_be_valid}"
            )

    print(f"Queue validation: {passed}/{len(test_cases)} tests passed\n")
    return passed == len(test_cases)


def test_platform_loading():
    """Test platform loading from YAML."""
    print("🧪 Testing platform loading...")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        platform_file = create_test_platform_file(temp_path)

        # Create loader with test directory
        loader = PlatformLoader(search_paths=[temp_path])

        # Check if platform was loaded
        platforms = loader.list_platforms()
        if "test_platform" in platforms:
            print("  ✅ Platform loaded successfully")

            platform = loader.get_platform("test_platform")
            if platform and len(platform.queues) == 3:
                print(f"  ✅ Platform has {len(platform.queues)} queues")
                return True
            else:
                print("  ❌ Platform queues not loaded correctly")
                return False
        else:
            print("  ❌ Platform not found in loader")
            return False


def test_queue_selection():
    """Test intelligent queue selection."""
    print("🧪 Testing queue selection...")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        platform_file = create_test_platform_file(temp_path)

        loader = PlatformLoader(search_paths=[temp_path])
        platform = loader.get_platform("test_platform")

        if not platform:
            print("  ❌ Could not load test platform")
            return False

        selector = QueueSelector(platform)

        test_cases = [
            # GPU job should select gpu queue
            ({"cpus": 2, "memory": "8GB", "gpus": 1}, "gpu"),
            # Short job should select express queue
            ({"cpus": 4, "memory": "8GB", "walltime": "1:00:00"}, "express"),
            # Regular job should select normal queue
            ({"cpus": 8, "memory": "16GB", "walltime": "12:00:00"}, "normal"),
        ]

        passed = 0
        for resources, expected_queue in test_cases:
            result = selector.select_queue(resources)
            if result.best_queue == expected_queue:
                print(f"  ✅ {resources} -> {result.best_queue}")
                passed += 1
            else:
                print(
                    f"  ❌ {resources} -> {result.best_queue}, expected {expected_queue}"
                )
                if result.warnings:
                    print(f"      Warnings: {result.warnings}")

        print(f"Queue selection: {passed}/{len(test_cases)} tests passed\n")
        return passed == len(test_cases)


def test_cost_estimation():
    """Test SU cost estimation."""
    print("🧪 Testing cost estimation...")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        platform_file = create_test_platform_file(temp_path)

        loader = PlatformLoader(search_paths=[temp_path])
        platform = loader.get_platform("test_platform")

        if not platform:
            print("  ❌ Could not load test platform")
            return False

        selector = QueueSelector(platform)

        # Test cost estimation for normal queue (1.0 SU/CPU·hour)
        resources = {
            "cpus": 4,
            "memory": "8GB",
            "walltime": "4:00:00",
        }  # Longer walltime to avoid express queue
        result = selector.select_queue(resources)

        expected_cost = 4 * 4.0 * 1.0  # 4 CPUs * 4 hours * 1.0 rate

        if result.estimated_cost and abs(result.estimated_cost - expected_cost) < 0.001:
            print(
                f"  ✅ Cost estimation: {result.estimated_cost} SU (expected {expected_cost})"
            )
            return True
        else:
            print(
                f"  ❌ Cost estimation: {result.estimated_cost}, expected {expected_cost}"
            )
            return False


def test_high_level_function():
    """Test the high-level select_best_queue function."""
    print("🧪 Testing high-level queue selection...")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        platform_file = create_test_platform_file(temp_path)

        # Mock the global platform loader to use our test directory
        import qxub.platform as platform_module

        original_loader = platform_module._platform_loader
        platform_module._platform_loader = PlatformLoader(search_paths=[temp_path])

        try:
            resources = {"cpus": 4, "memory": "8GB", "walltime": "4:00:00"}
            result = select_best_queue(resources, platform_name="test_platform")

            if result.best_queue:
                print(f"  ✅ Selected queue: {result.best_queue}")
                return True
            else:
                print(f"  ❌ No queue selected")
                if result.warnings:
                    print(f"      Warnings: {result.warnings}")
                return False
        finally:
            # Restore original loader
            platform_module._platform_loader = original_loader


def main():
    """Run all platform tests."""
    print("🚀 Testing qxub platform system...\n")

    tests = [
        test_queue_limits,
        test_walltime_rule,
        test_queue_validation,
        test_platform_loading,
        test_queue_selection,
        test_cost_estimation,
        test_high_level_function,
    ]

    passed = 0
    for test in tests:
        if test():
            passed += 1

    print(f"📊 Overall result: {passed}/{len(tests)} test groups passed")

    if passed == len(tests):
        print("🎉 All platform tests passed!")
        return 0
    else:
        print("❌ Some tests failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
