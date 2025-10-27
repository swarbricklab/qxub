#!/usr/bin/env python3
"""
Integration test for remote execution - checks the full flow without SSH.
"""
import os
import sys

sys.path.insert(0, "/g/data/a56/software/qsub_tools")

# Mock SSH execution by intercepting subprocess
import subprocess
from unittest.mock import MagicMock, patch

from qxub.execution_context import ExecutionContext
from qxub.remote.command_builder import build_remote_command
from qxub.remote.platform_executor import PlatformRemoteExecutor


def test_full_remote_execution_flow():
    """Test the complete remote execution flow with mocked SSH."""
    print("Testing full remote execution flow...")

    # Setup execution context and options
    execution_context = ExecutionContext("conda", "pytorch", "conda")
    options = {
        "queue": "normal",
        "mem": "8GB",
        "cpus": 4,
        "walltime": "1:00:00",
        "project": "a56",
        "dry": True,
        "verbose": 2,
    }
    command = ["python", "train.py", "--epochs", "100"]

    # Test 1: Command serialization
    print("\n1. Command serialization:")
    remote_command = build_remote_command(execution_context, options, command)
    print(f"   Serialized: {remote_command}")
    assert "qxub exec" in remote_command
    assert "--env pytorch" in remote_command
    assert "--dry" in remote_command

    # Test 2: Platform executor setup
    print("\n2. Platform executor setup:")
    remote_config = {
        "host": "gadi",
        "working_dir": "/scratch/a56/{user}",
        "conda_init": 'eval "$(conda shell.bash hook)"',
    }
    executor = PlatformRemoteExecutor("nci_gadi", remote_config)
    print(f"   Platform: {executor.platform_name}")
    print(f"   Host: {executor.hostname}")
    print(f"   Working dir: {executor.working_dir}")

    # Test 3: SSH command building
    print("\n3. SSH command building:")
    ssh_command = executor._build_ssh_command(remote_command)
    print(f"   SSH command parts:")
    for i, part in enumerate(ssh_command):
        print(f"     [{i}] {part}")
    assert ssh_command[0] == "ssh"
    assert "gadi" in ssh_command
    assert "cd /scratch/a56/" in ssh_command[-1]
    assert "QXUB_PLATFORM=nci_gadi" in ssh_command[-1]
    assert remote_command in ssh_command[-1]

    # Test 4: Verify command structure (skip actual execution)
    print("\n4. Command structure verification:")
    wrapped_command = ssh_command[-1]  # Last element is the remote command
    print(f"   Wrapped command:")
    print(f"     {wrapped_command}")

    # Verify all required components
    assert "cd /scratch/a56/" in wrapped_command, "Missing cd to working dir"
    assert 'eval "$(conda shell.bash hook)"' in wrapped_command, "Missing conda init"
    assert (
        "export QXUB_PLATFORM=nci_gadi" in wrapped_command
    ), "Missing platform env var"
    assert "qxub exec" in wrapped_command, "Missing qxub exec command"
    assert "--env pytorch" in wrapped_command, "Missing conda env"
    assert "--dry" in wrapped_command, "Missing dry-run flag"
    assert "python train.py" in wrapped_command, "Missing user command"

    print("   ✓ All components present in wrapped command")

    print("\n✅ Full remote execution flow test passed!")


if __name__ == "__main__":
    test_full_remote_execution_flow()
