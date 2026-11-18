"""Queue ``auto`` respects each queue's walltime limit.

Gadi's longest walltime is 48 hours, and it gets shorter as the CPU count
grows (see ``walltime_rules`` in ``docs/platforms/nci_gadi.yaml``).
Selection used to ignore this and pick ``normal`` for any walltime, so
``qsub`` rejected the job.
"""

import pytest
from conftest import GADI_PLATFORM_DIR

from qxub.platforms import CorePlatformLoader


@pytest.fixture(scope="module")
def gadi():
    return CorePlatformLoader(search_paths=[GADI_PLATFORM_DIR]).get_platform("nci_gadi")


@pytest.mark.parametrize(
    "resources, expected",
    [
        ({"cpus": 1, "walltime": "48:00:00"}, "normal"),
        ({"cpus": 1}, "normal"),
        ({"cpus": 1500, "walltime": "8:00:00"}, "normalbw"),
        ({"cpus": 4, "memory": "600GB", "walltime": "24:00:00"}, "hugemem"),
        ({"internet": True, "walltime": "2:00:00"}, "copyq"),
    ],
)
def test_queues_within_their_walltime_limit_are_selected(gadi, resources, expected):
    assert gadi.select_queue(resources) == expected


@pytest.mark.parametrize(
    "resources",
    [
        {"cpus": 1, "walltime": "72:00:00"},
        {"cpus": 1500, "walltime": "20:00:00"},
        {"internet": True, "walltime": "12:00:00"},
    ],
    ids=["one_cpu_72h", "1500_cpus_20h", "copyq_12h"],
)
def test_nothing_is_selected_when_no_queue_allows_the_walltime(gadi, resources):
    assert gadi.select_queue(resources) is None


@pytest.mark.parametrize("command", ["exec", "interactive"])
def test_auto_rejects_a_walltime_no_queue_allows(sandbox, command):
    sandbox.use_gadi_platform()
    args = ["--dry-run", "--queue", "auto", "--runtime", "72h", "--cpus", "1"]
    if command == "exec":
        args += ["--", "true"]

    result = sandbox.run(command, *args)

    assert result.returncode == 1, result.stdout
    assert "No queue allows a walltime of" in result.stderr
    assert "--runtime" in result.stderr
    assert sandbox.pbs_calls() == []


def test_auto_within_the_limit_still_selects_normal(sandbox):
    sandbox.use_gadi_platform()

    result = sandbox.run(
        "exec", "--dry-run", "--queue", "auto", "--runtime", "48h", "--", "true"
    )

    assert result.returncode == 0, result.stderr
    assert " -q normal " in result.stdout
