#!/usr/bin/env python3
"""Tests for enforce_queue_minimums() in config/handler.py."""

import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent))

from qxub.config.handler import enforce_queue_minimums, select_auto_queue


def _write_platform(temp_dir: Path) -> Path:
    """Write a minimal test platform YAML with gpuvolta-like constraints."""
    data = {
        "platform": {
            "name": "test_gadi",
            "type": "pbs_pro",
            "host": "gadi",
            "queues": [
                {
                    "name": "normal",
                    "type": "standard",
                    "priority": "normal",
                    "su_billing_rate": 2.0,
                    "limits": {"max_cpus": 48},
                },
                {
                    "name": "gpuvolta",
                    "type": "gpu",
                    "gpu_type": "v100",
                    "cpus_per_gpu": 12,
                    "priority": "normal",
                    "su_billing_rate": 3.0,
                    "limits": {
                        "max_cpus": 960,
                        "max_gpus": 4,
                        "min_gpus": 1,
                    },
                },
                {
                    "name": "diskqueue",
                    "type": "standard",
                    "priority": "normal",
                    "su_billing_rate": 1.0,
                    "limits": {
                        "max_cpus": 48,
                        "min_local_storage": "100GB",
                        "max_local_storage": "400GB",
                    },
                },
            ],
            "auto_selection_rules": [
                {"condition": "gpu_requested > 0", "queue": "gpuvolta"},
                {"condition": "cpus <= 48", "queue": "normal", "is_default": True},
            ],
        }
    }
    path = temp_dir / "test_gadi.yaml"
    path.write_text(yaml.dump(data))
    return path


def _params(queue, resources, cpus_explicit=False, queue_was_auto=False):
    return {
        "queue": queue,
        "queue_was_auto": queue_was_auto,
        "resources": list(resources),
        "cpus_explicit": cpus_explicit,
    }


# ---------------------------------------------------------------------------
# CPU enforcement — cpus_per_gpu multiple constraint
# ---------------------------------------------------------------------------


def test_cpu_below_minimum_is_raised(tmp_path):
    """CPUs below cpus_per_gpu minimum should be raised to cpus_per_gpu."""
    _write_platform(tmp_path)
    params = _params("gpuvolta", ["ncpus=1", "ngpus=1"], cpus_explicit=True)

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    ncpus = next(r for r in result["resources"] if r.startswith("ncpus="))
    assert ncpus == "ncpus=12", f"Expected ncpus=12, got {ncpus}"
    mock_echo.assert_called_once()
    warning_text = mock_echo.call_args[0][0]
    assert "1" in warning_text and "12" in warning_text
    print("  PASS: ncpus=1 raised to ncpus=12 for gpuvolta")


def test_cpu_non_multiple_rounded_up(tmp_path):
    """CPUs that aren't a multiple of cpus_per_gpu should be rounded up."""
    _write_platform(tmp_path)
    params = _params("gpuvolta", ["ncpus=13", "ngpus=1"], cpus_explicit=True)

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    ncpus = next(r for r in result["resources"] if r.startswith("ncpus="))
    assert ncpus == "ncpus=24", f"Expected ncpus=24, got {ncpus}"
    mock_echo.assert_called_once()
    print("  PASS: ncpus=13 rounded up to ncpus=24 for gpuvolta")


def test_cpu_exact_multiple_unchanged(tmp_path):
    """CPUs that are an exact multiple of cpus_per_gpu should not be changed."""
    _write_platform(tmp_path)
    params = _params("gpuvolta", ["ncpus=24", "ngpus=2"], cpus_explicit=True)

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    ncpus = next(r for r in result["resources"] if r.startswith("ncpus="))
    assert ncpus == "ncpus=24"
    mock_echo.assert_not_called()
    print("  PASS: ncpus=24 unchanged for gpuvolta")


def test_cpu_larger_compliant_multiple_unchanged(tmp_path):
    """User specifying 36 CPUs (multiple of 12) should not be overridden."""
    _write_platform(tmp_path)
    params = _params("gpuvolta", ["ncpus=36", "ngpus=3"], cpus_explicit=True)

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    ncpus = next(r for r in result["resources"] if r.startswith("ncpus="))
    assert ncpus == "ncpus=36"
    mock_echo.assert_not_called()
    print("  PASS: ncpus=36 unchanged for gpuvolta")


def test_no_cpus_in_resources_not_touched(tmp_path):
    """If ncpus isn't in resources, enforcement should not inject it."""
    _write_platform(tmp_path)
    params = _params("gpuvolta", ["ngpus=1"], cpus_explicit=False)

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    assert not any(r.startswith("ncpus=") for r in result["resources"])
    mock_echo.assert_not_called()
    print(
        "  PASS: no ncpus injected when not in resources (adjust_cpus_for_gpus handles that)"
    )


def test_queue_without_cpu_constraint_not_touched(tmp_path):
    """A queue with no cpus_per_gpu and no min_cpus should leave CPUs alone."""
    _write_platform(tmp_path)
    params = _params("normal", ["ncpus=3"])

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    ncpus = next(r for r in result["resources"] if r.startswith("ncpus="))
    assert ncpus == "ncpus=3"
    mock_echo.assert_not_called()
    print("  PASS: ncpus=3 unchanged on normal queue")


# ---------------------------------------------------------------------------
# CPU enforcement — auto-selected vs. explicit queue
# ---------------------------------------------------------------------------


def test_warning_mentions_auto_selected_context(tmp_path):
    """Warning text should say 'auto-selected queue' when queue_was_auto=True."""
    _write_platform(tmp_path)
    params = _params(
        "gpuvolta", ["ncpus=1", "ngpus=1"], cpus_explicit=True, queue_was_auto=True
    )

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            enforce_queue_minimums(params)

    warning_text = mock_echo.call_args[0][0]
    assert "auto-selected" in warning_text
    print("  PASS: warning mentions 'auto-selected' for auto queue path")


def test_warning_mentions_explicit_queue_context(tmp_path):
    """Warning text should say 'queue gpuvolta' when queue was explicitly named."""
    _write_platform(tmp_path)
    params = _params(
        "gpuvolta", ["ncpus=5", "ngpus=1"], cpus_explicit=True, queue_was_auto=False
    )

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            enforce_queue_minimums(params)

    warning_text = mock_echo.call_args[0][0]
    assert "auto-selected" not in warning_text
    assert "gpuvolta" in warning_text
    print("  PASS: warning mentions queue name for explicit queue path")


# ---------------------------------------------------------------------------
# Disk (jobfs) enforcement
# ---------------------------------------------------------------------------


def test_disk_below_minimum_is_raised(tmp_path):
    """jobfs below min_local_storage should be raised with a warning."""
    _write_platform(tmp_path)
    params = _params("diskqueue", ["ncpus=4", "jobfs=50GB"])

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    jobfs = next(r for r in result["resources"] if r.startswith("jobfs="))
    assert jobfs == "jobfs=100GB", f"Expected jobfs=100GB, got {jobfs}"
    mock_echo.assert_called_once()
    print("  PASS: jobfs=50GB raised to jobfs=100GB for diskqueue")


def test_disk_above_minimum_unchanged(tmp_path):
    """jobfs above min_local_storage should not be changed."""
    _write_platform(tmp_path)
    params = _params("diskqueue", ["ncpus=4", "jobfs=200GB"])

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    jobfs = next(r for r in result["resources"] if r.startswith("jobfs="))
    assert jobfs == "jobfs=200GB"
    mock_echo.assert_not_called()
    print("  PASS: jobfs=200GB unchanged for diskqueue")


def test_disk_not_specified_gets_minimum_injected_silently(tmp_path):
    """When jobfs not specified and queue has min_local_storage, inject it silently."""
    _write_platform(tmp_path)
    params = _params("diskqueue", ["ncpus=4"])

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo") as mock_echo:
            result = enforce_queue_minimums(params)

    jobfs = next((r for r in result["resources"] if r.startswith("jobfs=")), None)
    assert jobfs == "jobfs=100GB", f"Expected jobfs=100GB injected, got {jobfs}"
    mock_echo.assert_not_called()
    print("  PASS: jobfs=100GB injected silently when not specified")


def test_disk_not_specified_no_minimum_not_injected(tmp_path):
    """When queue has no min_local_storage, no jobfs should be injected."""
    _write_platform(tmp_path)
    params = _params("gpuvolta", ["ncpus=12", "ngpus=1"])

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        result = enforce_queue_minimums(params)

    assert not any(r.startswith("jobfs=") for r in result["resources"])
    print("  PASS: no jobfs injected when queue has no min_local_storage")


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


def test_unknown_queue_returns_params_unchanged(tmp_path):
    """Unknown queue name should return params unchanged without raising."""
    _write_platform(tmp_path)
    params = _params("nonexistent", ["ncpus=4"])

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        result = enforce_queue_minimums(params)

    assert result["resources"] == ["ncpus=4"]
    print("  PASS: unknown queue leaves params unchanged")


def test_auto_queue_sentinel_skipped():
    """queue='auto' should be skipped (not yet resolved)."""
    params = _params("auto", ["ncpus=1"])
    result = enforce_queue_minimums(params)
    assert result["resources"] == ["ncpus=1"]
    print("  PASS: queue='auto' sentinel skipped")


def test_other_resources_preserved(tmp_path):
    """Non-CPU/disk resources should be preserved through enforcement."""
    _write_platform(tmp_path)
    params = _params(
        "gpuvolta",
        ["mem=32GB", "ncpus=1", "ngpus=1", "walltime=2:00:00"],
        cpus_explicit=True,
    )

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        with patch("click.echo"):
            result = enforce_queue_minimums(params)

    assert "mem=32GB" in result["resources"]
    assert "walltime=2:00:00" in result["resources"]
    assert "ngpus=1" in result["resources"]
    print("  PASS: other resources preserved through CPU enforcement")


# ---------------------------------------------------------------------------
# Integration: select_auto_queue sets queue_was_auto
# ---------------------------------------------------------------------------


def test_select_auto_queue_sets_flag(tmp_path):
    """select_auto_queue should set queue_was_auto=True."""
    _write_platform(tmp_path)
    params = {
        "queue": "auto",
        "resources": ["ncpus=4", "ngpus=1"],
        "internet": False,
        "gpu_type": None,
    }

    with patch.dict(os.environ, {"QXUB_PLATFORM_PATHS": str(tmp_path)}):
        result = select_auto_queue(params)

    assert result.get("queue_was_auto") is True
    print("  PASS: select_auto_queue sets queue_was_auto=True")


def test_select_explicit_queue_does_not_set_flag():
    """select_auto_queue should not set queue_was_auto when queue is not 'auto'."""
    params = {"queue": "normal", "resources": []}
    result = select_auto_queue(params)
    assert result.get("queue_was_auto") is not True
    print("  PASS: explicit queue does not set queue_was_auto")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def main():
    tests = [
        test_cpu_below_minimum_is_raised,
        test_cpu_non_multiple_rounded_up,
        test_cpu_exact_multiple_unchanged,
        test_cpu_larger_compliant_multiple_unchanged,
        test_no_cpus_in_resources_not_touched,
        test_queue_without_cpu_constraint_not_touched,
        test_warning_mentions_auto_selected_context,
        test_warning_mentions_explicit_queue_context,
        test_disk_below_minimum_is_raised,
        test_disk_above_minimum_unchanged,
        test_disk_not_specified_gets_minimum_injected_silently,
        test_disk_not_specified_no_minimum_not_injected,
        test_unknown_queue_returns_params_unchanged,
        test_auto_queue_sentinel_skipped,
        test_other_resources_preserved,
        test_select_auto_queue_sets_flag,
        test_select_explicit_queue_does_not_set_flag,
    ]

    print(f"Running {len(tests)} tests for enforce_queue_minimums...\n")
    passed = 0
    for test in tests:
        name = test.__name__
        try:
            with tempfile.TemporaryDirectory() as tmp:
                try:
                    test(Path(tmp))
                except TypeError:
                    test()
            passed += 1
        except Exception as e:
            print(f"  FAIL: {name}: {e}")

    print(f"\n{passed}/{len(tests)} tests passed")
    return 0 if passed == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())
