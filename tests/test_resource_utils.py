#!/usr/bin/env python3
"""
Test script for resource_utils.py functionality.

Tests memory parsing, walltime parsing, condition evaluation, and resource comparisons.
"""

import os
import sys
from pathlib import Path

# Add qxub to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from qxub.resources.utils import (
    compare_memory,
    compare_walltime,
    evaluate_condition,
    format_memory_size,
    format_walltime,
    parse_memory_size,
    parse_walltime,
    suggest_resource_adjustment,
)


def test_memory_parsing():
    """Test memory size parsing functionality."""
    print("🧪 Testing memory parsing...")

    test_cases = [
        ("4GB", 4 * 1024 * 1024 * 1024),
        ("512MB", 512 * 1024 * 1024),
        ("2TB", 2 * 1024 * 1024 * 1024 * 1024),
        ("1024", 1024 * 1024 * 1024),  # Default to MB
        ("8G", 8 * 1024 * 1024 * 1024),
        ("256M", 256 * 1024 * 1024),
        ("1.5GB", int(1.5 * 1024 * 1024 * 1024)),
    ]

    passed = 0
    for memory_str, expected in test_cases:
        result = parse_memory_size(memory_str)
        if result == expected:
            print(f"  ✅ {memory_str} -> {result} bytes")
            passed += 1
        else:
            print(f"  ❌ {memory_str} -> {result}, expected {expected}")

    print(f"Memory parsing: {passed}/{len(test_cases)} tests passed\n")
    assert passed == len(test_cases), f"{passed}/{len(test_cases)} cases passed"


def test_memory_formatting():
    """Test memory size formatting."""
    print("🧪 Testing memory formatting...")

    test_cases = [
        (4 * 1024 * 1024 * 1024, "4GB"),
        (512 * 1024 * 1024, "512MB"),
        (1024, "1KB"),
        (512, "512B"),
    ]

    passed = 0
    for bytes_size, expected in test_cases:
        result = format_memory_size(bytes_size)
        if result == expected:
            print(f"  ✅ {bytes_size} bytes -> {result}")
            passed += 1
        else:
            print(f"  ❌ {bytes_size} bytes -> {result}, expected {expected}")

    print(f"Memory formatting: {passed}/{len(test_cases)} tests passed\n")
    assert passed == len(test_cases), f"{passed}/{len(test_cases)} cases passed"


def test_walltime_parsing():
    """Test walltime parsing functionality."""
    print("🧪 Testing walltime parsing...")

    test_cases = [
        ("48:00:00", 48.0),
        ("2:30:00", 2.5),
        ("1:15:30", 1.25833333),  # Approximate
        ("1h", 1.0),
        ("30m", 0.5),
        ("45s", 45 / 3600),
        ("2", 2.0),  # Default to hours
    ]

    passed = 0
    for walltime_str, expected in test_cases:
        result = parse_walltime(walltime_str)
        if result is not None and abs(result - expected) < 0.001:
            print(f"  ✅ {walltime_str} -> {result:.3f} hours")
            passed += 1
        else:
            print(f"  ❌ {walltime_str} -> {result}, expected {expected}")

    print(f"Walltime parsing: {passed}/{len(test_cases)} tests passed\n")
    assert passed == len(test_cases), f"{passed}/{len(test_cases)} cases passed"


def test_walltime_formatting():
    """Test walltime formatting."""
    print("🧪 Testing walltime formatting...")

    test_cases = [
        (48.0, "48:00:00"),
        (2.5, "02:30:00"),
        (1.25, "01:15:00"),
        (0.5, "00:30:00"),
    ]

    passed = 0
    for hours, expected in test_cases:
        result = format_walltime(hours)
        if result == expected:
            print(f"  ✅ {hours} hours -> {result}")
            passed += 1
        else:
            print(f"  ❌ {hours} hours -> {result}, expected {expected}")

    print(f"Walltime formatting: {passed}/{len(test_cases)} tests passed\n")
    assert passed == len(test_cases), f"{passed}/{len(test_cases)} cases passed"


def test_comparisons():
    """Test memory and walltime comparisons."""
    print("🧪 Testing comparisons...")

    memory_tests = [
        ("4GB", "2GB", 1),
        ("512MB", "1GB", -1),
        ("1GB", "1024MB", 0),
    ]

    walltime_tests = [
        ("48:00:00", "24:00:00", 1),
        ("1:30:00", "2:00:00", -1),
        ("1h", "1:00:00", 0),
    ]

    passed = 0
    total = len(memory_tests) + len(walltime_tests)

    for mem1, mem2, expected in memory_tests:
        result = compare_memory(mem1, mem2)
        if result == expected:
            print(f"  ✅ {mem1} vs {mem2} -> {result}")
            passed += 1
        else:
            print(f"  ❌ {mem1} vs {mem2} -> {result}, expected {expected}")

    for time1, time2, expected in walltime_tests:
        result = compare_walltime(time1, time2)
        if result == expected:
            print(f"  ✅ {time1} vs {time2} -> {result}")
            passed += 1
        else:
            print(f"  ❌ {time1} vs {time2} -> {result}, expected {expected}")

    print(f"Comparisons: {passed}/{total} tests passed\n")
    assert passed == total, f"{passed}/{total} cases passed"


def test_condition_evaluation():
    """Test resource condition evaluation."""
    print("🧪 Testing condition evaluation...")

    resources = {"cpus": 4, "memory": "8GB", "walltime": "2:00:00", "gpus": 1}

    test_cases = [
        ("cpus > 2", True),
        ("cpus <= 4", True),
        ("memory > 4GB", True),
        ("memory < 16GB", True),
        ("walltime > 1:00:00", True),
        ("walltime <= 2:00:00", True),
        ("gpu_requested > 0", True),
        ("gpu_requested == 1", True),
        ("cpus > 10", False),
        ("memory > 16GB", False),
    ]

    passed = 0
    for condition, expected in test_cases:
        result = evaluate_condition(condition, resources)
        if result == expected:
            print(f"  ✅ '{condition}' -> {result}")
            passed += 1
        else:
            print(f"  ❌ '{condition}' -> {result}, expected {expected}")

    print(f"Condition evaluation: {passed}/{len(test_cases)} tests passed\n")
    assert passed == len(test_cases), f"{passed}/{len(test_cases)} cases passed"


def test_resource_adjustment():
    """Test resource adjustment suggestions."""
    print("🧪 Testing resource adjustment suggestions...")

    resources = {"cpus": 2, "memory": "2GB", "walltime": "1:00:00"}
    queue_limits = {
        "min_cpus": 4,
        "max_cpus": 48,
        "min_memory": "4GB",
        "max_memory": "192GB",
    }

    suggestions = suggest_resource_adjustment(resources, queue_limits, "suggest")

    if suggestions and "cpus" in suggestions and "memory" in suggestions:
        print(f"  ✅ Generated suggestions: {suggestions}")
        return
    else:
        print(f"  ❌ No suggestions generated or missing expected keys")
        raise AssertionError("see output above")


def main():
    """Run all resource utility tests."""
    print("🚀 Testing qxub resource utilities...\n")

    tests = [
        test_memory_parsing,
        test_memory_formatting,
        test_walltime_parsing,
        test_walltime_formatting,
        test_comparisons,
        test_condition_evaluation,
        test_resource_adjustment,
    ]

    passed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  ❌ {test.__name__}: {e}")

    print(f"📊 Overall result: {passed}/{len(tests)} test groups passed")

    if passed == len(tests):
        print("🎉 All resource utility tests passed!")
        return 0
    else:
        print("❌ Some tests failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
