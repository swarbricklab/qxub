"""Helpers for tests against a real scheduler in a container (ci/README.md).

These tests run only inside a test cluster container, where
``ci/run-cluster-tests.sh <scheduler>`` sets ``QXUB_TEST_CLUSTER``. Each
module marks itself with :func:`needs` and is skipped anywhere else, and
imports the :func:`cluster` fixture from here. (Not a ``conftest.py``: a
second one would clash with ``tests/conftest.py``, which these import.)
"""

import os
import subprocess
import time
from pathlib import Path

import pytest
from conftest import FAKE_PBS_BIN, FAKE_SLURM_BIN, QxubSandbox

HERE = Path(__file__).resolve().parent
CLUSTER = os.environ.get("QXUB_TEST_CLUSTER", "")


def needs(scheduler: str):
    """Skip a module unless it runs in ``scheduler``'s test cluster."""
    return pytest.mark.skipif(
        CLUSTER != scheduler,
        reason=f"needs the {scheduler} test cluster (ci/run-cluster-tests.sh)",
    )


@pytest.fixture
def cluster(tmp_path):
    """A :class:`QxubSandbox` that reaches the container's real scheduler:
    no fake commands on ``PATH``, and the cluster's platform selected.

    ``monitoring.joblog_check_interval_sec`` is 5, so qx asks the scheduler
    whether a job ended every 5 s rather than every minute.
    """
    box = QxubSandbox(tmp_path)
    fakes = {str(FAKE_PBS_BIN), str(FAKE_SLURM_BIN)}
    box.env["PATH"] = os.pathsep.join(
        p for p in box.env["PATH"].split(os.pathsep) if p not in fakes
    )
    box.env["QXUB_PLATFORM_FILE"] = str(HERE / f"{CLUSTER}.yaml")
    box.write_config("monitoring:\n  joblog_check_interval_sec: 5\n")
    return box


def start_qx(box, *command, args=()):
    """``qxub exec [args] -- command`` in the background."""
    return subprocess.Popen(
        ["qxub", "exec", *args, "--", *command],
        cwd=box.cwd,
        env=box.env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def wait_for(condition, timeout=60, what="the condition"):
    """Poll ``condition()`` until it returns something true; return that."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        found = condition()
        if found:
            return found
        time.sleep(1)
    pytest.fail(f"timed out after {timeout}s waiting for {what}")
