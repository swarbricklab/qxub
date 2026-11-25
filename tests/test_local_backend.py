"""The local backend runs qxub jobs as background processes, no scheduler.

These run real jobs (``job.sh`` in a background process) through the CLI,
with ``docs/platforms/local.yaml`` as the platform.
"""

import time

import pytest
from conftest import REPO_ROOT

from qxub.backends import BACKENDS, LocalBackend, backend_for_type
from qxub.core.scheduler import job_status_from_files

LOCAL_PLATFORM = REPO_ROOT / "docs" / "platforms" / "local" / "local.yaml"


@pytest.fixture
def local(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    return sandbox


def _wait_for(check, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        value = check()
        if value:
            return value
        time.sleep(0.2)
    raise AssertionError("timed out")


def _status(sandbox, job_id):
    result = sandbox.run("status", "check", job_id)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_platform_type_local_selects_it():
    assert BACKENDS["local"] is LocalBackend
    assert isinstance(backend_for_type("local"), LocalBackend)


def test_exec_runs_the_job_and_reports_it(local):
    result = local.run(
        "exec", "--pre", "export G=hi", "--",
        'echo "$G from local"; echo oops >&2; exit 3',
    )  # fmt: skip

    assert result.returncode == 3
    assert result.stdout.strip() == "hi from local"
    assert "oops" in result.stderr
    assert "(pre ✅, main ❌(3))" in result.stderr
    assert "local-" in result.stderr  # the job ID
    assert local.pbs_calls() == [], "no scheduler involved"


def test_dry_run_shows_the_local_command(local):
    result = local.run("exec", "--dry-run", "--", "true")

    assert result.returncode == 0, result.stderr
    (line,) = [l for l in result.stdout.splitlines() if l.startswith("Dry run")]
    assert "QXUB_JOB_FILE=" in line and line.endswith("job.sh")
    assert "export QXUB_JOB_ID_VAR=QXUB_LOCAL_JOBID" in result.stdout


def test_terse_job_runs_in_the_background(local):
    result = local.run("exec", "--terse", "--", "sleep 1; echo done")
    assert result.returncode == 0, result.stderr
    job_id = result.stdout.strip()
    assert LocalBackend().looks_like_job_id(job_id)

    assert _wait_for(lambda: _status(local, job_id) == "success")
    monitor = local.run("monitor", "--quiet", job_id)
    assert monitor.returncode == 0, monitor.stdout + monitor.stderr


def test_cancel_stops_a_local_job(local, monkeypatch):
    result = local.run("exec", "--terse", "--", "sleep 30")
    job_id = result.stdout.strip()
    monkeypatch.setenv("XDG_STATE_HOME", local.env["XDG_STATE_HOME"])
    monkeypatch.setenv("HOME", local.env["HOME"])
    backend = LocalBackend()
    _wait_for(lambda: job_id in [j.job_id for j in backend.list_jobs()])

    cancel = local.run("cancel", "id", "-f", job_id)
    assert cancel.returncode == 0, cancel.stdout + cancel.stderr

    info = _wait_for(
        lambda: (
            backend.find_info(job_id) if backend.find_info(job_id).is_finished else None
        )
    )
    assert info.exit_status == 143  # the job script's TERM trap
    assert info.ended_by == "cancelled"
    assert job_id not in [j.job_id for j in backend.list_jobs()]
    assert _status(local, job_id) == "failed"


def test_interactive_runs_the_session_script_here(local):
    result = local.run("interactive", "--dry-run", "--runtime", "15m")

    assert result.returncode == 0, result.stderr
    assert "QXUB_LOCAL_JOBID=local-" in result.stdout
    assert 'echo "${QXUB_LOCAL_JOBID}"' in result.stdout


def test_status_check_reads_a_finished_run_file(tmp_path):
    """A job stopped by a signal leaves no final_exit_code file, and its
    job log may be elsewhere; its run file still says how it ended."""
    status = tmp_path / "status"
    (status / "runs").mkdir(parents=True)
    (status / "native_77.gadi-pbs").write_text("qx-7\n")
    (status / "job_started_qx-7").write_text("now\n")
    run = status / "runs" / "qx-7.run"
    run.write_text("run_id=qx-7\nstate=running\n")

    assert job_status_from_files("77.gadi-pbs", status) == ("R", None)

    run.write_text(run.read_text() + "state=finished\nexit_code=143\nsignal=TERM\n")
    assert job_status_from_files("77.gadi-pbs", status) == ("F", 143)


def test_a_job_on_another_host_is_not_judged_by_its_pid(tmp_path, monkeypatch):
    """#166 review: with a shared home, another machine checked the PID
    locally and reported a running job as killed."""
    import json

    import qxub.backends.local as local_mod

    monkeypatch.setattr(local_mod, "jobs_dir", lambda: tmp_path)
    job = tmp_path / "local-elsewhere"
    job.mkdir()
    (job / "meta.json").write_text(
        json.dumps(
            {
                "name": "j",
                "queue": "local",
                "job_log": "/x.log",
                "submitted": "2026-11-24T00:00:00",
                "host": "some-other-host",
            }
        )
    )
    (job / "pid").write_text("999999\n")  # not running here
    backend = LocalBackend()

    info = backend.find_info("local-elsewhere")
    assert info.exit_status is None and not info.is_finished
    assert info.host == "some-other-host"
    assert [j.job_id for j in backend.list_jobs()] == ["local-elsewhere"]
    ok, reason = backend.cancel("local-elsewhere")
    assert not ok and "some-other-host" in reason

    (job / "exit_status").write_text("0\n")  # it finished over there
    assert backend.find_info("local-elsewhere").is_finished

    # On its own host, a missing process still means it was killed
    meta = json.loads((job / "meta.json").read_text())
    meta["host"] = local_mod.socket.gethostname()
    (job / "meta.json").write_text(json.dumps(meta))
    (job / "exit_status").unlink()
    assert (
        backend.job_details(["local-elsewhere"])["local-elsewhere"]["exit_status"] == -1
    )
