#!/usr/bin/env python3
"""
Test script to validate workflow-friendly resource options work correctly.
"""

import subprocess
import sys
import tempfile
from pathlib import Path


def test_workflow_resources():
    """Test that workflow resource options are converted to PBS format."""
    print("Testing workflow-friendly resource options...")

    # Change to the qxub project directory
    qxub_dir = Path(__file__).parent

    # Test 1: Memory conversion
    print("\n1. Testing memory conversion...")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "qxub.cli",
            "exec",
            "--mem",
            "4GB",
            "--dry",
            "echo 'test'",
        ],
        cwd=qxub_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        print("✓ Memory option accepted")
        if "mem=4GB" in result.stdout or "mem=4gb" in result.stdout:
            print("✓ Memory converted to PBS format")
        else:
            print("? Memory conversion not visible in dry run output")
    else:
        print("✗ Memory option failed")
        print("STDERR:", result.stderr)

    # Test 2: Runtime conversion
    print("\n2. Testing runtime conversion...")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "qxub.cli",
            "exec",
            "--runtime",
            "2h",
            "--dry",
            "echo 'test'",
        ],
        cwd=qxub_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        print("✓ Runtime option accepted")
        if "walltime=02:00:00" in result.stdout or "walltime=2:00:00" in result.stdout:
            print("✓ Runtime converted to PBS format")
        else:
            print("? Runtime conversion not visible in dry run output")
    else:
        print("✗ Runtime option failed")
        print("STDERR:", result.stderr)

    # Test 3: CPU conversion
    print("\n3. Testing CPU conversion...")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "qxub.cli",
            "exec",
            "--cpus",
            "4",
            "--dry",
            "echo 'test'",
        ],
        cwd=qxub_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        print("✓ CPU option accepted")
        if "ncpus=4" in result.stdout:
            print("✓ CPU converted to PBS format")
        else:
            print("? CPU conversion not visible in dry run output")
    else:
        print("✗ CPU option failed")
        print("STDERR:", result.stderr)

    # Test 4: Combined options
    print("\n4. Testing combined workflow options...")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "qxub.cli",
            "exec",
            "--mem",
            "8GB",
            "--runtime",
            "1h30m",
            "--cpus",
            "2",
            "--dry",
            "echo 'combined test'",
        ],
        cwd=qxub_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        print("✓ Combined workflow options accepted")
        print("Dry run output snippet:")
        # Show relevant lines from output
        lines = result.stdout.split("\n")
        for line in lines:
            if any(
                keyword in line.lower()
                for keyword in ["mem=", "walltime=", "ncpus=", "qsub"]
            ):
                print(f"  {line}")
    else:
        print("✗ Combined options failed")
        print("STDERR:", result.stderr)

    # Test 5: Mixed PBS and workflow options
    print("\n5. Testing mixed PBS and workflow options...")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "qxub.cli",
            "exec",
            "-l",
            "wd=true",
            "--mem",
            "4GB",
            "--dry",
            "echo 'mixed test'",
        ],
        cwd=qxub_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        print("✓ Mixed PBS and workflow options accepted")
        # Check that both wd=true and mem=4GB appear
        if "wd=true" in result.stdout and (
            "mem=4GB" in result.stdout or "mem=4gb" in result.stdout
        ):
            print("✓ Both PBS and workflow resources preserved")
        else:
            print("? Resource mixing not clearly visible")
    else:
        print("✗ Mixed options failed")
        print("STDERR:", result.stderr)


if __name__ == "__main__":
    test_workflow_resources()
