"""Submit, cancel and accounting paths exercised against the fake PBS.

The fake ``qsub``/``qstat``/``qdel`` live in ``tests/fakes/pbs`` and are put
first on ``PATH`` by the ``sandbox`` fixture (see ``conftest.py``).
"""

import json
from pathlib import Path

import pytest
from conftest import FAKE_PBS_BIN, job_variables

from qxub.core import scheduler

FIXTURES = Path(__file__).parent / "fixtures"


# -- submission --------------------------------------------------------------


def test_terse_submit_prints_the_job_id(sandbox):
    result = sandbox.run("exec", "--terse", "-P", "ab12", "--", "echo", "hi")

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1000.gadi-pbs"
    (call,) = sandbox.pbs_calls("qsub")
    assert call["argv"][0] == "-v"
    assert call["argv"][-1].endswith("job.sh")
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
    assert "violates queue" in result.stderr
    # `id=$(qx --terse ...)` must not capture the error as a job ID.
    assert "violates queue" not in result.stdout
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


def test_cancel_exits_non_zero_when_qdel_fails(sandbox):
    sandbox.fail_qdel(35, "qdel: Unknown Job Id 999.gadi-pbs")

    assert sandbox.run("cancel", "id", "--force", "999.gadi-pbs").returncode == 1


def test_cancel_without_job_ids_exits_non_zero(sandbox):
    assert sandbox.run("cancel", "id", input="").returncode == 1


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


# -- --array and --bind --------------------------------------------------------


def test_array_is_rejected_without_submitting(sandbox):
    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--array", "1-10", "--", "true"
    )

    assert result.returncode == 1
    assert "Job arrays (--array) are not supported yet" in result.stderr
    assert sandbox.pbs_calls() == []


def test_array_from_a_shortcut_is_rejected(sandbox):
    sandbox.config_dir.mkdir(parents=True)
    (sandbox.config_dir / "shortcuts.json").write_text('{"sweep": {"array": "1-4"}}')

    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--shortcut", "sweep", "--", "true"
    )

    assert result.returncode == 1, result.stdout + result.stderr
    assert sandbox.pbs_calls() == []


def test_multi_path_bind_reaches_the_job_as_one_variable(sandbox):
    sandbox.run(
        "exec", "--terse", "-P", "ab12", "--sif", "/c/t.sif",
        "--bind", "/a,/b", "--", "true",
    )  # fmt: skip

    (call,) = sandbox.pbs_calls("qsub")
    variables = job_variables(call["argv"])
    assert (variables["sif"], variables["bind"]) == ("/c/t.sif", "/a,/b")


@pytest.mark.parametrize(
    "context",
    [[], ["--mod", "python3"], ["--env", "myenv"], ["--sif", "/c/my image.sif"]],
)
def test_paths_with_spaces_reach_qsub_whole(sandbox, context):
    """#140 review: an unquoted cwd was split by shlex and qsub refused it."""
    execdir = sandbox.cwd / "exec dir with spaces"
    execdir.mkdir()
    out = sandbox.cwd / "log dir" / "job.out"
    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--execdir", str(execdir),
        "--out", str(out), *context, "--", "pwd",
    )  # fmt: skip

    assert result.returncode == 0, result.stdout + result.stderr
    (call,) = sandbox.pbs_calls("qsub")
    assert call["argv"][-1].endswith("job.sh")
    # The default joblog is inside the execdir and goes to -o (#140 retest).
    joblog = call["argv"][call["argv"].index("-o") + 1]
    assert joblog.startswith(str(execdir) + "/"), call["argv"]
    variables = job_variables(call["argv"])
    assert variables["cwd"] == str(execdir)
    assert variables["out"] == str(out)
    if "--sif" in context:
        assert variables["sif"] == "/c/my image.sif"


def test_fake_pbs_is_first_on_path(sandbox):
    assert sandbox.env["PATH"].split(":")[0] == str(FAKE_PBS_BIN)
