"""Submit, cancel and accounting paths exercised against the fake PBS.

The fake ``qsub``/``qstat``/``qdel`` live in ``tests/fakes/pbs`` and are put
first on ``PATH`` by the ``sandbox`` fixture (see ``conftest.py``).
"""

import json
from pathlib import Path

import pytest
from conftest import FAKE_PBS_BIN

from qxub.core import scheduler

FIXTURES = Path(__file__).parent / "fixtures"


# -- submission --------------------------------------------------------------


def test_terse_submit_prints_the_job_id(sandbox):
    result = sandbox.run("exec", "--terse", "-P", "ab12", "--", "echo", "hi")

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1000.gadi-pbs"
    (call,) = sandbox.pbs_calls("qsub")
    assert call["argv"][0] == "-v"
    assert call["argv"][-1].endswith("qdefault.pbs")
    assert Path(call["cwd"]) == sandbox.cwd


def test_job_ids_are_sequential(sandbox):
    ids = [
        sandbox.run("exec", "--terse", "-P", "ab12", "--", "true").stdout.strip()
        for _ in range(2)
    ]
    assert ids == ["1000.gadi-pbs", "1001.gadi-pbs"]


@pytest.mark.parametrize("mode", ["--terse", "--quiet"])
def test_qsub_failure_exit_code_and_message_are_passed_through(sandbox, mode):
    sandbox.fail_qsub(38, "qsub: Job violates queue and/or server resource limits")

    result = sandbox.run("exec", mode, "-P", "ab12", "--", "echo", "hi")

    assert result.returncode == 38
    assert "violates queue" in result.stdout + result.stderr
    assert len(sandbox.pbs_calls("qsub")) == 1


def test_dry_run_does_not_call_qsub(sandbox):
    result = sandbox.run("exec", "--dry-run", "-P", "ab12", "--", "echo", "hi")

    assert result.returncode == 0
    assert sandbox.pbs_calls() == []


# -- cancel ------------------------------------------------------------------


def test_cancel_by_id_calls_qdel(sandbox):
    result = sandbox.run("cancel", "id", "--force", "1000.gadi-pbs", "1001.gadi-pbs")

    assert result.returncode == 0, result.stderr
    assert [c["argv"] for c in sandbox.pbs_calls("qdel")] == [
        ["1000.gadi-pbs"],
        ["1001.gadi-pbs"],
    ]


def test_cancel_reports_a_qdel_failure(sandbox):
    sandbox.fail_qdel(35, "qdel: Unknown Job Id 999.gadi-pbs")

    result = sandbox.run("cancel", "id", "--force", "999.gadi-pbs")

    assert len(sandbox.pbs_calls("qdel")) == 1
    assert "Failed to cancel 1 job(s): 999.gadi-pbs" in result.stdout


@pytest.mark.xfail(
    strict=True,
    reason="known bug: cancel returns 1 from a click command, which click "
    "ignores, so a failed qdel still exits 0",
)
def test_cancel_exits_non_zero_when_qdel_fails(sandbox):
    sandbox.fail_qdel(35, "qdel: Unknown Job Id 999.gadi-pbs")

    assert sandbox.run("cancel", "id", "--force", "999.gadi-pbs").returncode != 0


def test_cancel_dry_run_does_not_call_qdel(sandbox):
    result = sandbox.run("cancel", "id", "--dry-run", "1000.gadi-pbs")

    assert result.returncode == 0
    assert sandbox.pbs_calls() == []


# -- qstat -fx accounting ----------------------------------------------------


@pytest.fixture
def fake_pbs_in_process(sandbox, monkeypatch):
    """Put the fake PBS on PATH for code called directly from the test."""
    monkeypatch.setenv("PATH", sandbox.env["PATH"])
    monkeypatch.setenv("FAKE_PBS_STATE", str(sandbox.pbs_state))
    return sandbox


def test_job_resource_data_from_qstat_fx(fake_pbs_in_process, golden):
    sandbox = fake_pbs_in_process
    sandbox.set_qstat("1000.gadi-pbs", (FIXTURES / "qstat_fx_finished.txt").read_text())

    data = scheduler.get_job_resource_data("1000.gadi-pbs")

    assert sandbox.pbs_calls("qstat")[0]["argv"] == ["-fx", "1000.gadi-pbs"]
    assert data["exit_status"] == 0
    golden("qstat/finished_job.json", json.dumps(data, indent=2, sort_keys=True) + "\n")


@pytest.mark.usefixtures("fake_pbs_in_process")
def test_job_resource_data_for_unknown_job_is_none():
    assert scheduler.get_job_resource_data("424242.gadi-pbs") is None


# -- known bugs pinned by the golden files -----------------------------------


@pytest.mark.xfail(
    strict=True, reason="known bug: --array is accepted but never sent to qsub"
)
def test_array_reaches_qsub(sandbox):
    sandbox.run("exec", "--terse", "-P", "ab12", "--array", "1-10", "--", "true")

    (call,) = sandbox.pbs_calls("qsub")
    assert "-J" in call["argv"]


@pytest.mark.xfail(
    strict=True,
    reason="known bug: the quotes around a comma-separated --bind are removed "
    "by shlex.split before qsub, so PBS splits the value at the comma",
)
def test_multi_path_bind_survives_qsub_variable_parsing(sandbox):
    sandbox.run(
        "exec", "--terse", "-P", "ab12", "--sif", "/c/t.sif",
        "--bind", "/a,/b", "--", "true",
    )  # fmt: skip

    (call,) = sandbox.pbs_calls("qsub")
    variables = call["argv"][call["argv"].index("-v") + 1]
    assert '"/a,/b"' in variables


def test_fake_pbs_is_first_on_path(sandbox):
    assert sandbox.env["PATH"].split(":")[0] == str(FAKE_PBS_BIN)
