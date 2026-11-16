"""Shared pytest fixtures for the qxub test suite."""

import subprocess
import sys

import pytest

# Upper bound for a cold `import qxub`, in seconds.
IMPORT_TIME_LIMIT = 1.0


def measure_cold_import_time(module: str = "qxub") -> float:
    """Time importing ``module`` in a fresh interpreter.

    Timing ``import qxub`` inside the test process only measures a
    ``sys.modules`` cache hit, so the import runs in a subprocess instead.
    """
    code = (
        "import time; t = time.perf_counter(); "
        f"import {module}; print(time.perf_counter() - t)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    )
    return float(result.stdout.strip().splitlines()[-1])


@pytest.fixture
def cold_import_time():
    """Seconds taken by a cold ``import qxub``."""
    return measure_cold_import_time()
