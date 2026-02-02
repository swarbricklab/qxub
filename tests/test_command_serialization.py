#!/usr/bin/env python3
"""
Quick test script for command serialization.
"""
import sys

sys.path.insert(0, "/g/data/a56/software/qsub_tools")

from qxub.execution_context import ExecutionContext
from qxub.remote.command_builder import build_remote_command


def test_command_serialization():
    """Test that command serialization works correctly."""
    print("Testing command serialization...")

    # Test 1: Conda environment with basic options
    print("\n1. Conda environment with basic options:")
    ctx = ExecutionContext("conda", "pytorch", "conda")
    options = {
        "queue": "auto",
        "mem": "8GB",
        "cpus": 4,
        "walltime": "2:00:00",
        "project": "a56",
    }
    command = ["python", "train.py", "--epochs", "100"]
    result = build_remote_command(ctx, options, command)
    print(f"   Command: {result}")
    assert "qxub exec" in result
    assert "--env pytorch" in result
    assert "--mem 8GB" in result
    assert "python train.py" in result

    # Test 2: Module environment
    print("\n2. Module environment:")
    ctx = ExecutionContext("module", ["python3", "numpy"], "module")
    options = {"queue": "normal", "runtime": "1h30m"}
    command = ["python", "script.py"]
    result = build_remote_command(ctx, options, command)
    print(f"   Command: {result}")
    assert "--mod python3" in result
    assert "--mod numpy" in result
    assert "--runtime 1h30m" in result

    # Test 3: Singularity container
    print("\n3. Singularity container:")
    ctx = ExecutionContext("singularity", "/path/to/container.sif", "singularity")
    options = {"bind": "/data:/mnt", "cpus": 2}
    command = ["ls", "/mnt"]
    result = build_remote_command(ctx, options, command)
    print(f"   Command: {result}")
    assert "--sif" in result
    assert "--bind" in result
    assert "ls /mnt" in result

    # Test 4: Default execution
    print("\n4. Default execution:")
    ctx = ExecutionContext("default", None, "default")
    options = {"dry": True, "verbose": 2}
    command = ["echo", "hello"]
    result = build_remote_command(ctx, options, command)
    print(f"   Command: {result}")
    assert "--default" in result
    assert "--dry" in result
    assert "-vv" in result
    assert "echo hello" in result

    print("\n✅ All serialization tests passed!")


if __name__ == "__main__":
    test_command_serialization()
