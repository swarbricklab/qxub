"""Jobs record their state in a run file, not the job database (Phase 2).

See ``docs/dev/sqlite-robustness.md``. On Gadi's ``localflock`` Lustre
mounts, SQLite writes from many compute nodes corrupted the shared database.
Now the submitting process writes ``<status dir>/runs/<run_id>.run`` before
``qsub``, the job adds to it with plain shell, and login-side readers apply
it to the database.
"""

import base64
import os
import signal
import subprocess
import time
from pathlib import Path

import pytest
from conftest import REPO_ROOT, job_variables

from qxub.config import site
from qxub.queue import connection, db, runs
from qxub.resources.tracker import ResourceTracker

TEMPLATE = REPO_ROOT / "qxub" / "jobscripts" / "job.sh"


@pytest.fixture
def state_home(tmp_path, monkeypatch):
    """No platform, so run files go under ``$XDG_STATE_HOME``."""
    monkeypatch.setattr(site, "current_platform", lambda: None)
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "state"))
    connection.forget_schema_cache()
    return tmp_path / "state" / "qxub" / "status" / "runs"


# -- the run file format ----------------------------------------------------------


def test_write_and_read_round_trip(state_home):
    path = runs.write_submitted(
        "qx-1",
        {"state": "submitted", "tags": ["a=1", "b"], "note": "two\nlines"},
    )

    assert path == state_home / "qx-1.run"
    record = runs.read_run(path)
    assert record["run_id"] == "qx-1"
    assert record["tags"] == '["a=1","b"]'
    assert record["note"] == "two lines"
    assert not [p for p in state_home.iterdir() if p.name.startswith(".")]


def test_last_value_wins(tmp_path):
    path = tmp_path / "r.run"
    path.write_text("state=submitted\nx=a=b\nstate=finished\nexit_code=3\n")
    record = runs.read_run(path)
    assert record["state"] == "finished"
    assert record["x"] == "a=b"
    assert runs.exit_code(record) == 3


def test_unwritable_run_dir_warns_once_and_returns_none(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(runs, "runs_dir", lambda: tmp_path / "file" / "runs")
    (tmp_path / "file").write_text("not a directory")
    monkeypatch.setattr(runs, "_warned", False)

    assert runs.write_submitted("qx-1", {}) is None
    assert runs.write_submitted("qx-2", {}) is None
    assert sum(r.levelname == "WARNING" for r in caplog.records) == 1


# -- the job script updates it ------------------------------------------------------


def template_env(tmp_path, command: str, run_file: Path) -> dict:
    return {
        "PATH": "/usr/bin:/bin",
        "HOME": str(tmp_path),
        "cmd_b64": base64.b64encode(command.encode()).decode(),
        "cwd": str(tmp_path),
        "out": str(tmp_path / "logs" / "job.out"),
        "err": str(tmp_path / "logs" / "job.err"),
        "quiet": "false",
        "PBS_JOBID": "123.gadi-pbs",
        "QXUB_RUN_ID": "qx-1",
        "QXUB_RUN_FILE": str(run_file),
    }


def run_template(tmp_path, command: str, run_file: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", str(TEMPLATE)],
        env=template_env(tmp_path, command, run_file),
        capture_output=True,
        text=True,
        timeout=30,
    )


@pytest.mark.parametrize("command, code", [("true", 0), ("exit 3", 3), ("false", 1)])
def test_job_records_running_and_finished(tmp_path, command, code):
    run_file = tmp_path / "runs" / "qx-1.run"
    run_file.parent.mkdir()
    run_file.write_text("run_id=qx-1\nstate=submitted\nqueue=normal\n")

    result = run_template(tmp_path, command, run_file)

    assert result.returncode == code, result.stdout + result.stderr
    record = runs.read_run(run_file)
    assert record["state"] == "finished"
    assert runs.exit_code(record) == code
    assert record["signal"] == ""
    assert record["native_job_id"] == "123.gadi-pbs"
    assert record["started_at"] and record["finished_at"] and record["host"]
    assert record["queue"] == "normal", "submission fields are kept"
    assert [p.name for p in run_file.parent.iterdir()] == ["qx-1.run"]


@pytest.mark.parametrize(
    "sig, code, name",
    [
        (signal.SIGTERM, 143, "TERM"),
        (signal.SIGINT, 130, "INT"),
        (signal.SIGHUP, 129, "HUP"),
    ],
)
def test_job_killed_by_a_signal_is_not_recorded_as_success(tmp_path, sig, code, name):
    """Review of #133: walltime and qdel kills were recorded as exit_code=0.

    PBS sends SIGTERM to the job's processes, as here to the process group.
    """
    run_file = tmp_path / "runs" / "qx-1.run"
    proc = subprocess.Popen(
        ["bash", str(TEMPLATE)],
        env=template_env(tmp_path, "sleep 30", run_file),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    deadline = time.time() + 20
    while (runs.read_run(run_file) or {}).get("state") != "running":
        assert time.time() < deadline, "job never recorded state=running"
        time.sleep(0.05)
    time.sleep(0.3)  # let it reach the main command

    os.killpg(proc.pid, sig)
    proc.wait(timeout=20)

    record = runs.read_run(run_file)
    assert record["state"] == "finished"
    assert runs.exit_code(record) == code
    assert record["signal"] == name


def test_job_creates_the_run_file_if_submission_could_not(tmp_path):
    run_file = tmp_path / "runs" / "qx-1.run"

    run_template(tmp_path, "true", run_file)

    assert runs.read_run(run_file)["state"] == "finished"


def test_job_without_run_file_variable_still_runs(tmp_path):
    result = subprocess.run(
        ["bash", str(TEMPLATE)],
        env={
            "PATH": "/usr/bin:/bin",
            "HOME": str(tmp_path),
            "cmd_b64": base64.b64encode(b"true").decode(),
            "cwd": str(tmp_path),
            "out": str(tmp_path / "job.out"),
            "err": str(tmp_path / "job.err"),
        },
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr


# -- submission writes it and tells the job where it is -------------------------------


def qsub_var(call, name):
    return job_variables(call["argv"]).get(name)


def test_exec_writes_run_file_and_passes_its_path(sandbox):
    result = sandbox.run(
        "exec", "--terse", "-q", "normal", "--tags", "w=1", "--", "echo", "hi"
    )
    assert result.returncode == 0, result.stderr

    (call,) = sandbox.pbs_calls("qsub")
    run_id = qsub_var(call, "QXUB_RUN_ID")
    run_file = Path(qsub_var(call, "QXUB_RUN_FILE"))
    assert run_file.name == f"{run_id}.run"
    assert qsub_var(call, "QXUB_SHARED_DB") is None

    record = runs.read_run(run_file)
    assert record["state"] == "submitted"
    assert base64.b64decode(record["command_b64"]).decode() == "echo hi"
    assert record["queue"] == "normal"
    assert "w=1" in record["tags"]


def test_submit_job_writes_run_file(sandbox, monkeypatch):
    for key, value in sandbox.env.items():
        monkeypatch.setenv(key, value)
    monkeypatch.chdir(sandbox.cwd)
    monkeypatch.setattr(site, "current_platform", lambda: None)
    from qxub.execution.submit import submit_job

    result = submit_job(["echo", "hi"], queue="normal", tags=["w=1"])

    (call,) = sandbox.pbs_calls("qsub")
    run_file = Path(qsub_var(call, "QXUB_RUN_FILE"))
    assert run_file.name == f"{result.virtual_id}.run"
    assert runs.read_run(run_file)["state"] == "submitted"


# -- readers apply it to the database ---------------------------------------------------


@pytest.fixture
def tracked_job(state_home, tmp_path):
    """A job submitted to a fresh DB, with its run file."""
    db_path = tmp_path / "qxub.db"
    run_id = db.create_queue_entry(
        "echo hi", db_path=db_path, entry_id="qx-abc", pbs_job_id="77.gadi-pbs"
    )
    tracker = ResourceTracker(db_path=db_path)
    tracker.log_job_submitted("77.gadi-pbs", command="echo hi")
    runs.write_submitted(run_id, {"state": "submitted"})
    return tracker, runs.run_file(run_id)


def append(run_file, text):
    run_file.write_text(run_file.read_text() + text)


def test_sync_applies_finished_state(tracked_job):
    tracker, run_file = tracked_job
    append(run_file, "state=running\nstarted_at=2026-11-20T10:00:00\n")
    append(run_file, "state=finished\nexit_code=4\nfinished_at=2026-11-20T10:05:00\n")

    assert tracker.sync_from_run_files() == 1

    job = tracker.get_job_status("77.gadi-pbs")
    assert job["status"] == "failed"
    assert job["exit_code"] == 4
    assert job["completed_at"] == "2026-11-20T10:05:00"
    (legacy,) = tracker.get_jobs_by_status("failed")
    assert legacy["job_id"] == "77.gadi-pbs"
    assert tracker.sync_from_run_files() == 0, "finished rows are not revisited"


def test_sync_applies_running_state(tracked_job):
    tracker, run_file = tracked_job
    append(run_file, "state=running\nstarted_at=2026-11-20T10:00:00\n")

    assert tracker.sync_from_run_files(["77.gadi-pbs"]) == 1
    assert tracker.get_job_status("77.gadi-pbs")["status"] == "running"
    assert tracker.sync_from_run_files(["77.gadi-pbs"]) == 0


def test_sync_ignores_jobs_without_run_files(tracked_job):
    tracker, run_file = tracked_job
    run_file.unlink()
    assert tracker.sync_from_run_files() == 0
    assert tracker.sync_from_run_files([]) == 0


def test_qxtat_check_uses_the_run_file(tracked_job, monkeypatch):
    from click.testing import CliRunner

    from qxub import status_cli

    tracker, run_file = tracked_job
    append(run_file, "state=finished\nexit_code=0\n")
    monkeypatch.setattr(status_cli, "resource_tracker", tracker)

    result = CliRunner().invoke(status_cli.check, ["77.gadi-pbs"])

    assert result.output.strip() == "success"
