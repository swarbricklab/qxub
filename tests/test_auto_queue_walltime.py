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


def test_a_later_rule_stands_in_when_it_can_take_the_whole_request(gadi):
    # normalsl allows 10h at 1000 CPUs; the default, normal, allows 24h.
    assert gadi.select_queue({"cpus": 1000}) == "normalsl"
    assert gadi.select_queue({"cpus": 1000, "walltime": "20:00:00"}) == "normal"


GPU = {"gpu_requested": 1, "gpus": 1}


@pytest.mark.parametrize(
    "resources",
    [
        {**GPU, "cpus": 192, "memory": "96GB", "walltime": "30:00:00"},
        {"cpus": 96, "memory": "600GB", "walltime": "30:00:00"},
        {"cpus": 192, "memory": "600GB", "walltime": "11:00:00"},
    ],
    ids=["gpu_job_to_cpu_queue", "memory_too_big", "too_many_cpus"],
)
def test_no_stand_in_that_cannot_take_the_rest_of_the_request(gadi, resources):
    # Each of these once fell through to a queue that drops the GPUs or
    # fails the queue's own memory or CPU limits.
    assert gadi.select_queue(resources) is None


def _grid():
    for cpus in (1, 12, 48, 96, 192, 1000, 1500, 3500):
        for memory in (None, "96GB", "300GB", "600GB", "1200GB"):
            for walltime in (None, "4:00:00", "11:00:00", "30:00:00", "72:00:00"):
                for extra in ({}, GPU, {"internet": True}):
                    request = {"cpus": cpus, **extra}
                    if memory:
                        request["memory"] = memory
                    if walltime:
                        request["walltime"] = walltime
                    yield request


def test_any_stand_in_queue_can_take_the_whole_request(gadi):
    """A queue other than the rules' own choice must fit the whole request."""
    for request in _grid():
        picked = gadi.select_queue(request)
        rules_choice = gadi.select_queue(
            {k: v for k, v in request.items() if k != "walltime"}
        )
        if picked is None or picked == rules_choice:
            continue
        queue = gadi.get_queue(picked)
        assert not queue.validate_resources(request).errors, (request, picked)
        assert bool(queue.gpu_type) == ("gpus" in request), (request, picked)


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


def test_gpu_job_over_the_gpu_walltime_limit_is_rejected(sandbox):
    sandbox.use_gadi_platform()

    result = sandbox.run(
        "exec", "--dry-run", "--queue", "auto", "-l", "ngpus=1", "--cpus", "192",
        "--mem", "96GB", "--runtime", "30h", "--", "true",
    )  # fmt: skip

    assert result.returncode == 1, result.stdout
    assert "No queue allows a walltime of" in result.stderr


def test_auto_within_the_limit_still_selects_normal(sandbox):
    sandbox.use_gadi_platform()

    result = sandbox.run(
        "exec", "--dry-run", "--queue", "auto", "--runtime", "48h", "--", "true"
    )

    assert result.returncode == 0, result.stderr
    assert " -q normal " in result.stdout
