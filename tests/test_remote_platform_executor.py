#!/usr/bin/env python3
"""
Test remote execution mode detection and platform executor.
"""
import sys

sys.path.insert(0, "/g/data/a56/software/qsub_tools")

from qxub.execution.mode import ExecutionMode, get_execution_mode
from qxub.remote.platform_executor import PlatformRemoteExecutor


def test_execution_mode_detection():
    """Test that execution mode is correctly detected."""
    print("Testing execution mode detection...")

    # Test 1: Local platform (no remote section)
    print("\n1. Local platform (no remote section):")
    local_config = {
        "name": "nci_gadi",
        "definition": "file:///path/to/nci_gadi.yaml",
    }
    mode = get_execution_mode(local_config)
    print(f"   Mode: {mode}")
    assert mode == ExecutionMode.LOCAL, f"Expected LOCAL, got {mode}"

    # Test 2: Remote platform (has remote section)
    print("\n2. Remote platform (has remote section):")
    remote_config = {
        "name": "pawsey_setonix",
        "definition": "https://example.com/pawsey_setonix.yaml",
        "remote": {
            "host": "setonix",
            "working_dir": "/scratch/{project}/{user}",
        },
    }
    mode = get_execution_mode(remote_config)
    print(f"   Mode: {mode}")
    assert mode == ExecutionMode.REMOTE, f"Expected REMOTE, got {mode}"

    # Test 3: Remote platform with conda_init
    print("\n3. Remote platform with custom conda_init:")
    remote_config_custom = {
        "name": "custom_platform",
        "definition": "file:///path/to/custom.yaml",
        "remote": {
            "host": "custom-hpc",
            "working_dir": "/work/{user}",
            "conda_init": "source /opt/conda/etc/profile.d/conda.sh",
        },
    }
    mode = get_execution_mode(remote_config_custom)
    print(f"   Mode: {mode}")
    assert mode == ExecutionMode.REMOTE, f"Expected REMOTE, got {mode}"

    print("\n✅ All execution mode detection tests passed!")


def test_platform_executor():
    """Test PlatformRemoteExecutor initialization and command building."""
    print("\n\nTesting PlatformRemoteExecutor...")

    # Test 1: Basic initialization
    print("\n1. Basic initialization:")
    remote_config = {
        "host": "gadi",
        "working_dir": "/scratch/a56/{user}",
    }
    executor = PlatformRemoteExecutor("nci_gadi", remote_config)
    print(f"   Platform: {executor.platform_name}")
    print(f"   Hostname: {executor.hostname}")
    print(f"   Working dir: {executor.working_dir}")
    assert executor.platform_name == "nci_gadi"
    assert executor.hostname == "gadi"
    assert "/scratch/a56/" in executor.working_dir

    # Test 2: Custom conda_init
    print("\n2. Custom conda_init:")
    remote_config_custom = {
        "host": "setonix",
        "working_dir": "/scratch/{project}/{user}",
        "conda_init": "module load conda",
    }
    executor = PlatformRemoteExecutor("pawsey_setonix", remote_config_custom)
    print(f"   Conda init: {executor.conda_init}")
    assert executor.conda_init == "module load conda"

    # Test 3: Command wrapping
    print("\n3. Command wrapping:")
    remote_command = "qxub exec --env pytorch -- python train.py"
    wrapped = executor._wrap_remote_command(remote_command)
    print(f"   Wrapped: {wrapped}")
    assert "cd /scratch/" in wrapped
    assert "QXUB_PLATFORM=pawsey_setonix" in wrapped
    assert remote_command in wrapped
    assert "module load conda" in wrapped

    print("\n✅ All platform executor tests passed!")


if __name__ == "__main__":
    test_execution_mode_detection()
    test_platform_executor()
