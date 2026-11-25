"""Why a job ended is recorded and shown, not just its exit code (#65).

``qxub cancel`` notes it in the run file before stopping a running job;
the PBS backend reads a walltime kill from the exit status (-29); ``qx``
says which it was, and the database and statistics keep it.
"""

import subprocess
import sys
import time

from conftest import _CLI_SHIM, REPO_ROOT
from test_consolidate import JOBLOG, NOW, area, finished, pbs_info, rows, run

from qxub.backends import (
    ENDED_BY_CANCEL,
    ENDED_BY_WALLTIME,
    LocalBackend,
    PbsBackend,
    describe_ended_by,
)
from qxub.core import scheduler
from qxub.queue import runs
from qxub.resources.parser import parse_joblog_resources

LOCAL_PLATFORM = REPO_ROOT / "docs" / "platforms" / "local" / "local.yaml"

__all__ = ["area"]  # the fixture, used by name below


def test_backends_read_the_reason_from_an_exit_status():
    assert PbsBackend().ended_by(-29) == ENDED_BY_WALLTIME
    assert PbsBackend().ended_by(271) == ENDED_BY_CANCEL  # qdel (#170 review)
    assert PbsBackend().ended_by(143) is None
    assert LocalBackend().ended_by(-29) is None
    assert describe_ended_by(ENDED_BY_WALLTIME) == "stopped at its walltime limit"


def test_joblog_parser_keeps_a_negative_exit_status(tmp_path):
    log = tmp_path / "j.log"
    log.write_text(JOBLOG.format(when=NOW, job_id="1.gadi-pbs", code=-29))
    assert parse_joblog_resources(str(log))["exit_status"] == -29


# -- cancelling ---------------------------------------------------------------


class _Backend(LocalBackend):
    def __init__(self, ok):
        self.ok = ok

    def cancel(self, job_id):
        return self.ok, "" if self.ok else "no such job"


def _running_job(tmp_path, monkeypatch, state="running"):
    status = tmp_path / "status"
    (status / "runs").mkdir(parents=True)
    (status / "native_77.gadi-pbs").write_text("qx-r\n")
    path = status / "runs" / "qx-r.run"
    path.write_text(runs.format_fields({"run_id": "qx-r", "state": state}))
    monkeypatch.setattr(scheduler, "get_default_status_dir", lambda: status)
    return path


def test_cancel_marks_a_running_job(tmp_path, monkeypatch):
    path = _running_job(tmp_path, monkeypatch)

    assert scheduler.cancel_job("77.gadi-pbs", _Backend(True)) == (True, "")
    assert runs.read_run(path)["ended_by"] == ENDED_BY_CANCEL


def test_failed_cancel_takes_the_mark_back(tmp_path, monkeypatch):
    path = _running_job(tmp_path, monkeypatch)

    ok, _ = scheduler.cancel_job("77.gadi-pbs", _Backend(False))
    assert not ok and not runs.read_run(path)["ended_by"]


def test_finished_job_is_not_marked(tmp_path, monkeypatch):
    path = _running_job(tmp_path, monkeypatch, state="finished")
    scheduler.cancel_job("77.gadi-pbs", _Backend(True))
    assert "ended_by" not in runs.read_run(path)


# -- qx -----------------------------------------------------------------------


def test_completion_message_says_why():
    message = scheduler._format_completion_message
    assert message("1", 143, "/nowhere", "k", ENDED_BY_CANCEL).startswith(
        "🛑 Job was cancelled (exit code 143)"
    )
    assert message("1", -29, "/nowhere", "k", ENDED_BY_WALLTIME).startswith(
        "❌ Job was stopped at its walltime limit (exit code -29)"
    )
    assert message("1", 3, "/nowhere", "k").startswith("❌ Job failed with exit code 3")


def test_walltime_from_the_joblog_is_added_to_the_run_file(tmp_path):
    path = tmp_path / "qx-w.run"
    path.write_text(runs.format_fields({"state": "finished", "exit_code": 143}))

    assert scheduler.job_end_reason("5.gadi-pbs", -29, path) == ENDED_BY_WALLTIME
    assert runs.read_run(path)["ended_by"] == ENDED_BY_WALLTIME
    assert scheduler.job_end_reason("5.gadi-pbs", 1, None) is None


def test_qx_stops_when_its_job_is_cancelled(sandbox):
    """Before, a job stopped by `qxub cancel` wrote no final_exit_code
    file, so qx noticed only at its next joblog check (about 50 s on Gadi)."""
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    qx = subprocess.Popen(
        [sys.executable, "-c", _CLI_SHIM, "exec", "--", "echo begun; sleep 60"],
        cwd=sandbox.cwd,
        env=sandbox.env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        local = sandbox.home / ".local" / "state" / "qxub" / "local"
        deadline = time.time() + 20
        job_id = None
        while time.time() < deadline and job_id is None:
            started = [
                d.name
                for d in (local.iterdir() if local.exists() else [])
                if (d / "pid").exists()
            ]
            job_id = started[0] if started else None
            time.sleep(0.2)
        assert job_id, "the job didn't start"
        time.sleep(1)  # let it record that it's running

        assert sandbox.run("cancel", "id", "-f", job_id).returncode == 0
        out, _ = qx.communicate(timeout=20)
    finally:
        qx.kill()

    assert "🛑 Job was cancelled (exit code 143)" in out, out
    assert qx.returncode == 143


# -- storing it ---------------------------------------------------------------


def test_consolidation_stores_why(area):
    runs_dir, snapshot = area
    joblog = runs_dir.parent / "w.log"
    joblog.write_text(JOBLOG.format(when=NOW, job_id="1.gadi-pbs", code=-29))
    # The job's TERM trap recorded 143; PBS's log says it was the walltime
    finished(runs_dir, "qx-w", "1.gadi-pbs", code=143, joblog=joblog)
    finished(runs_dir, "qx-c", "2.gadi-pbs", code=143, ended_by=ENDED_BY_CANCEL)
    finished(runs_dir, "qx-q", "3.gadi-pbs", code=1)

    run(
        area,
        lookup=lambda ids: {
            j: pbs_info(j, state="F", code=143 if j == "2.gadi-pbs" else 1) for j in ids
        },
    )

    got = {r["job_id"]: r["ended_by"] for r in rows(snapshot)}
    assert got == {
        "1.gadi-pbs": ENDED_BY_WALLTIME,
        "2.gadi-pbs": ENDED_BY_CANCEL,
        "3.gadi-pbs": None,
    }


def test_status_show_says_why(sandbox):
    """The reason reaches the database from the run file, as qx records
    a walltime kill it read from the job log."""
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    job_id = sandbox.run("exec", "--terse", "--", "exit 3").stdout.strip()
    sandbox.run("monitor", "--quiet", job_id)
    status = sandbox.home / ".local" / "state" / "qxub" / "status"
    (path,) = (status / "runs").glob("*.run")
    runs.update(path, {"ended_by": ENDED_BY_WALLTIME})

    out = sandbox.run("status", "show", job_id).stdout

    assert "Exit Code: 3" in out and "Ended: stopped at its walltime limit" in out, out


def test_local_job_cancelled_under_another_platform_is_marked(sandbox):
    """#170 review: the run file is in the submitting platform's status
    directory, which a cancel under Gadi wouldn't look in."""
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    job_id = sandbox.run("exec", "--terse", "--", "sleep 30").stdout.strip()
    status = sandbox.home / ".local" / "state" / "qxub" / "status"
    deadline = time.time() + 10
    while time.time() < deadline and not list(status.glob("job_started_*")):
        time.sleep(0.2)
    del sandbox.env["QXUB_PLATFORM_FILE"]
    sandbox.use_gadi_platform()
    sandbox.env["PROJECT"] = "a56"

    assert sandbox.run("cancel", "id", "-f", job_id).returncode == 0

    (path,) = (status / "runs").glob("*.run")
    assert runs.read_run(path).get("ended_by") == ENDED_BY_CANCEL
