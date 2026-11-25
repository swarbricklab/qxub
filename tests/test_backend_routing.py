"""Each job goes to its own backend, whatever the current platform (#65, 7b).

Listing live jobs covers every backend this machine can ask; consolidation
asks the backend a run file recorded; a remote qxub announces submissions
in a backend-neutral way; and jobs get their own name as ``QXUB_JOB_NAME``.
"""

import subprocess
import sys
import time

import pytest
from conftest import REPO_ROOT

from qxub import backends
from qxub.backends import LocalBackend
from qxub.queue import consolidate as cons
from qxub.remote import handshake
from qxub.remote.platform_executor import PlatformRemoteExecutor

LOCAL_PLATFORM = REPO_ROOT / "docs" / "platforms" / "local" / "local.yaml"


def _local_job(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    result = sandbox.run("exec", "--terse", "--", "sleep 30")
    del sandbox.env["QXUB_PLATFORM_FILE"]
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def _local_job_finished(sandbox, job_id, timeout=10):
    done = sandbox.home / ".local" / "state" / "qxub" / "local" / job_id
    deadline = time.time() + timeout
    while time.time() < deadline:
        if (done / "exit_status").exists():
            return True
        time.sleep(0.2)
    return False


# -- cancelling ---------------------------------------------------------------


def test_cancel_all_cancels_only_this_platforms_jobs_unless_asked(sandbox):
    """#168 review: under the local platform, `cancel all -f` shouldn't
    also kill your cluster jobs; other backends' jobs are opt-in."""
    local_id = _local_job(sandbox)
    sandbox.use_gadi_platform()
    sandbox.env["PROJECT"] = "a56"
    pbs = sandbox.run("exec", "--terse", "--", "true")
    pbs_id = pbs.stdout.strip()
    assert pbs.returncode == 0 and pbs_id in sandbox.live_jobs(), pbs.stderr

    only_pbs = sandbox.run("cancel", "all", input="n\n")
    assert pbs_id in only_pbs.stdout and local_id not in only_pbs.stdout
    assert "(1 " not in only_pbs.stdout, "no breakdown for one backend"

    asked = sandbox.run("cancel", "all", "--all-backends", input="n\n")
    assert "2 job(s)" in asked.stdout and "(1 pbs, 1 local)" in asked.stdout
    assert sandbox.pbs_calls("qdel") == [], "answered no"

    result = sandbox.run("cancel", "all", "--all-backends", "-f")
    assert result.returncode == 0, result.stdout + result.stderr
    assert [c["argv"][-1] for c in sandbox.pbs_calls("qdel")] == [pbs_id]
    assert _local_job_finished(sandbox, local_id)


def test_cancel_name_finds_a_local_job_under_gadi(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    result = sandbox.run("exec", "--terse", "--name", "mine", "--", "sleep 30")
    job_id = result.stdout.strip()
    del sandbox.env["QXUB_PLATFORM_FILE"]
    sandbox.use_gadi_platform()
    sandbox.env["PROJECT"] = "a56"

    result = sandbox.run("cancel", "name", "-f", "mine")

    assert result.returncode == 0, result.stdout + result.stderr
    assert sandbox.pbs_calls("qdel") == []
    assert _local_job_finished(sandbox, job_id)


def test_local_jobs_are_only_yours(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path))
    (tmp_path / "qxub" / "local").mkdir(parents=True)
    monkeypatch.setenv("USER", "me")
    assert LocalBackend().available()
    assert LocalBackend().list_jobs("someone-else") == []


def test_another_backend_that_fails_is_skipped(monkeypatch):
    class Broken(LocalBackend):
        def list_jobs(self, user=None, full_names=False):
            raise backends.SchedulerQueryError("no")

    monkeypatch.setitem(backends.BACKENDS, "broken", Broken)
    monkeypatch.setattr(backends, "get_backend", lambda platform=None: LocalBackend())
    monkeypatch.setattr(LocalBackend, "list_jobs", lambda *a, **k: [])
    monkeypatch.setattr(backends.PbsBackend, "available", lambda self: False)

    assert backends.live_jobs() == []


# -- consolidation ------------------------------------------------------------


def test_consolidation_asks_the_backend_the_run_file_recorded(tmp_path, monkeypatch):
    runs_dir = tmp_path / "runs"
    runs_dir.mkdir()
    from qxub.queue import runs

    (runs_dir / "qx-a.run").write_text(
        runs.format_fields(
            {
                "run_id": "qx-a",
                "state": "finished",
                "exit_code": 0,
                "native_job_id": "12345",
                "backend": "local",
                "submitted_at": "2026-11-20T09:00:00",
                "finished_at": "2026-11-20T10:00:00",
            }
        )
    )
    asked = []
    cons.consolidate(
        runs_dir,
        tmp_path / "stats.db",
        qstat=lambda ids: asked.append(dict(ids)) or {},
        now=cons.datetime(2026, 11, 20, 12),
    )
    assert asked == [{"12345": "local"}]

    chosen = []
    real = backends.backend_for_job
    monkeypatch.setattr(
        backends,
        "backend_for_job",
        lambda job_id, name=None: chosen.append((job_id, name)) or real(job_id, name),
    )
    monkeypatch.setattr(LocalBackend, "job_details", lambda self, ids: {})
    cons.query_qstat({"12345": "local"})
    assert chosen == [("12345", "local")]


# -- remote submission handshake ----------------------------------------------


def _executor(monkeypatch, script):
    """An executor whose "ssh" runs ``script`` locally."""
    executor = PlatformRemoteExecutor("p", {"host": "h", "working_dir": "/w"})
    monkeypatch.setattr(
        executor, "_build_ssh_command", lambda command: ["bash", "-c", script]
    )
    monkeypatch.setenv("QXUB_REMOTE_RETRY_DELAY", "0")
    return executor


@pytest.mark.parametrize("stream", [True, False])
def test_dropped_connection_after_an_announced_job_is_not_retried(
    monkeypatch, capfd, tmp_path, stream
):
    attempts = tmp_path / "attempts"
    script = (
        f"echo x >> {attempts}; echo working; "
        f"echo '{handshake.MARKER} job_id=987 backend=slurm' >&2; exit 255"
    )
    executor = _executor(monkeypatch, script)

    assert executor.execute("qxub exec -- true", stream_output=stream) == 255

    assert attempts.read_text().count("x") == 1, "not retried"
    out, err = capfd.readouterr()
    assert "job 987 was already submitted" in err
    assert handshake.MARKER not in out + err.replace("already submitted", "")


def test_dropped_connection_before_any_job_is_retried(monkeypatch, capfd, tmp_path):
    attempts = tmp_path / "attempts"
    executor = _executor(monkeypatch, f"echo x >> {attempts}; exit 255")
    monkeypatch.setenv("QXUB_REMOTE_EXEC_MAX_RETRIES", "1")

    assert executor.execute("qxub exec -- true") == 255
    assert attempts.read_text().count("x") == 2
    assert "No submission was reported" in capfd.readouterr().err


def test_remote_command_asks_for_the_announcement():
    executor = PlatformRemoteExecutor("p", {"host": "h", "working_dir": "/w"})
    assert f"export {handshake.ENV}=1" in executor._wrap_remote_command("qxub x")


def test_qxub_announces_its_submission_when_asked(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    result = sandbox.run(
        "exec", "--terse", "--", "true", env={handshake.ENV: "1"}
    )  # fmt: skip

    assert result.returncode == 0, result.stderr
    job_id = result.stdout.strip()
    assert handshake.find(result.stderr) == (job_id, "local")

    quiet = sandbox.run("exec", "--terse", "--", "true")
    assert handshake.MARKER not in quiet.stderr


def test_announcement_is_cleared_for_jobs_started_from_here(monkeypatch):
    monkeypatch.setenv(handshake.ENV, "1")
    assert handshake.requested()
    assert not handshake.requested()


# -- job name -----------------------------------------------------------------


def test_jobs_get_their_name(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    result = sandbox.run("exec", "--name", "named", "--", 'echo "[$QXUB_JOB_NAME]"')

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "[named]"
