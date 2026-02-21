#!/usr/bin/env python3
"""
Test script for interactive defaults configuration.
"""

import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

# Add qxub to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from qxub.interactive_cli import interactive_cli


def test_interactive_cli_defaults():
    print("🧪 Testing interactive CLI defaults (memory & cpus)...")
    runner = CliRunner()

    # Mock _get_config_manager to return a config with default memory/cpus
    with patch("qxub.interactive_cli._get_config_manager") as mock_get_config_manager:
        mock_config_mgr = MagicMock()

        # Mock getting config values
        def get_config_value_side_effect(key):
            if key == "defaults.interactive.mem":
                return "8GB"  # Interactive specific default
            if key == "defaults.mem":
                return "4GB"  # General default
            # CPU defaults
            if key == "defaults.interactive.cpus":
                return "4"
            if key == "defaults.cpus":
                return "2"
            # Return None for others to use coded defaults
            return None

        mock_config_mgr.get_config_value.side_effect = get_config_value_side_effect
        mock_get_config_manager.return_value = mock_config_mgr

        # Don't need to patch return value of get_config_value generally, just side_effect

        # Run with --dry-run so we see output
        result = runner.invoke(interactive_cli, ["--dry-run"])

        success = True

        # Check for memory setting
        if "mem=8GB" in result.output:
            print("  ✅ Memory default correctly set to 8GB")
        else:
            print("  ❌ Memory default NOT found or incorrect")
            success = False

        # Check for CPU setting (it's passed as -l ncpus=4 or similar depending on implementation)
        # interactive_cli adds CPUs via mapper.add_cpus(cpus)
        # resulting in -l ncpus=4 in qsub command
        # Let's check output for ncpus=4
        if "ncpus=4" in result.output:
            print("  ✅ CPU default correctly set to 4")
        else:
            print("  ❌ CPU default NOT found or incorrect")
            # Might be cpus=4 depending on mapper implementation? Let's check mapper.
            # Usually PBS uses ncpus
            success = False

        if success:
            return True
        else:
            print("Output was:")
            print(result.output)
            return False


if __name__ == "__main__":
    if test_interactive_cli_defaults():
        sys.exit(0)
    else:
        sys.exit(1)
