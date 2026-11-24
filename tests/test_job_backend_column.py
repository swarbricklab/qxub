"""Step 7: the job database records each job's scheduler backend.

With two backends (PBS and local), a job must be asked about and cancelled
through the backend that ran it, whatever the current platform is.
"""

import sqlite3

import pytest
from conftest import REPO_ROOT

from qxub import backends
from qxub.backends import LocalBackend, PbsBackend, backend_for_job
from qxub.queue import connection, runs
from qxub.resources.tracker import ResourceTracker

LOCAL_PLATFORM = REPO_ROOT / "docs" / "platforms" / "local" / "local.yaml"


@pytest.fixture(autouse=True)
def fresh_schema_cache():
    connection.forget_schema_cache()
    yield
    connection.forget_schema_cache()


def test_version_2_database_gets_the_column_and_a_backfill(tmp_path):
    path = tmp_path / "qxub.db"
    old = sqlite3.connect(path)
    old.execute(
        "CREATE TABLE queue (entry_id TEXT PRIMARY KEY, command_b64 TEXT,"
        " pbs_job_id TEXT, status TEXT, user TEXT)"
    )
    old.executemany(
        "INSERT INTO queue (entry_id, command_b64, pbs_job_id, status, user)"
        " VALUES (?, '', ?, 'done', 'me')",
        [("qx-1", "123.gadi-pbs"), ("qx-2", "local-abc"), ("qx-3", None)],
    )
    old.execute("CREATE TABLE job_resources (job_id TEXT PRIMARY KEY, timestamp TEXT)")
    old.executemany(
        "INSERT INTO job_resources VALUES (?, '2026-01-01')",
        [("123.gadi-pbs",), ("local-abc",)],
    )
    old.execute("PRAGMA user_version=2")
    old.commit()
    old.close()

    connection.connect(path).close()

    conn = sqlite3.connect(path)
    assert conn.execute("PRAGMA user_version").fetchone()[0] == 3
    assert dict(conn.execute("SELECT entry_id, backend FROM queue")) == {
        "qx-1": "pbs",
        "qx-2": "local",
        "qx-3": None,  # never submitted
    }
    assert dict(conn.execute("SELECT job_id, backend FROM job_resources")) == {
        "123.gadi-pbs": "pbs",
        "local-abc": "local",
    }


def test_backend_for_job_prefers_the_recorded_backend(monkeypatch):
    monkeypatch.setattr(backends, "get_backend", lambda platform=None: PbsBackend())

    assert isinstance(backend_for_job("123.gadi-pbs"), PbsBackend)
    assert isinstance(backend_for_job("local-abc"), LocalBackend)  # by its shape
    assert isinstance(backend_for_job("123", "local"), LocalBackend)  # recorded
    assert isinstance(backend_for_job("qx-1"), PbsBackend)  # unknown: current


def test_submission_records_the_backend_everywhere(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    result = sandbox.run("exec", "--terse", "--", "true")
    assert result.returncode == 0, result.stderr
    job_id = result.stdout.strip()

    state = sandbox.home / ".local" / "state" / "qxub" / "status"
    (run_file,) = (state / "runs").glob("*.run")
    assert runs.read_run(run_file)["backend"] == "local"
    db = next(sandbox.home.rglob("qxub.db"))
    conn = sqlite3.connect(db)
    assert conn.execute(
        "SELECT backend FROM queue WHERE pbs_job_id=?", (job_id,)
    ).fetchone() == ("local",)
    assert conn.execute(
        "SELECT backend FROM job_resources WHERE job_id=?", (job_id,)
    ).fetchone() == ("local",)


@pytest.mark.parametrize("other", ["none", "gadi"])
def test_a_local_job_is_cancelled_locally_under_another_platform(sandbox, other):
    """Submitted with the local platform, cancelled under a PBS one.

    #167 review: the Gadi platform has a scratch path, so its status
    directory differs from the local platform's; the local backend must
    still find its jobs.
    """
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    job_id = sandbox.run("exec", "--terse", "--", "sleep 30").stdout.strip()

    del sandbox.env["QXUB_PLATFORM_FILE"]  # the current backend is now PBS
    if other == "gadi":
        sandbox.use_gadi_platform()
        sandbox.env["PROJECT"] = "a56"
    result = sandbox.run("cancel", "id", "-f", job_id)

    assert result.returncode == 0, result.stdout + result.stderr
    assert sandbox.pbs_calls("qdel") == [], "not sent to PBS"


def test_recorded_backend_reads_the_db(tmp_path, monkeypatch):
    monkeypatch.setenv("QXUB_DB_PATH", str(tmp_path / "qxub.db"))
    tracker = ResourceTracker()
    tracker.log_job_submitted("local-xyz", command="true", backend="local")
    tracker.log_job_submitted("9.gadi-pbs", command="true")  # older caller

    assert tracker.recorded_backend("local-xyz") == "local"
    assert tracker.recorded_backend("9.gadi-pbs") is None
    assert tracker.recorded_backend("nope") is None
