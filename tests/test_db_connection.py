"""Job database connections on shared filesystems (sqlite-robustness Phase 1).

See ``docs/dev/sqlite-robustness.md``. Each job used to open the shared DB
about four times, and every open ran schema DDL that needs a write lock.
"""

import sqlite3
import subprocess
import sys
import textwrap
import threading
import time

import pytest

from qxub.queue import connection, db
from qxub.resources.tracker import ResourceTracker


@pytest.fixture(autouse=True)
def fresh_schema_cache():
    connection.forget_schema_cache()
    yield
    connection.forget_schema_cache()


@pytest.fixture
def traced(monkeypatch):
    """Record every statement run on connections opened by ``connect``."""
    statements = []
    real_open = connection._open

    def tracing_open(path, timeout):
        conn = real_open(path, timeout)
        conn.set_trace_callback(statements.append)
        return conn

    monkeypatch.setattr(connection, "_open", tracing_open)
    return statements


def ddl(statements):
    return [s for s in statements if s.lstrip().upper().startswith(("CREATE", "ALTER"))]


# -- schema ---------------------------------------------------------------------


def test_new_database_gets_schema_and_version(tmp_path):
    path = tmp_path / "qxub.db"
    connection.connect(path).close()

    conn = sqlite3.connect(path)
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master")}
    assert {"queue", "project_jobs", "meta", "job_resources"} <= tables
    assert conn.execute("PRAGMA user_version").fetchone()[0] == (
        connection.SCHEMA_VERSION
    )
    assert conn.execute("PRAGMA journal_mode").fetchone()[0] == "delete"


def test_up_to_date_database_runs_no_ddl(tmp_path, traced):
    path = tmp_path / "qxub.db"
    connection.connect(path).close()
    assert ddl(traced), "first connection creates the schema"

    # Another process: nothing cached, but the version is current.
    connection.forget_schema_cache()
    traced.clear()
    db.create_queue_entry("echo hi", db_path=path)
    ResourceTracker(db_path=path).update_job_status("1.pbs", "running")

    assert ddl(traced) == []
    assert not [s for s in traced if "journal_mode=" in s.replace(" ", "")]


def test_schema_is_checked_once_per_process(tmp_path, traced):
    path = tmp_path / "qxub.db"
    connection.connect(path).close()
    traced.clear()

    for _ in range(3):
        connection.connect(path).close()

    assert not [s for s in traced if "user_version" in s]


def test_old_database_is_migrated_and_versioned(tmp_path):
    path = tmp_path / "qxub.db"
    old = sqlite3.connect(path)
    old.execute("PRAGMA journal_mode=WAL")
    old.execute("CREATE TABLE job_resources (job_id TEXT PRIMARY KEY, timestamp TEXT)")
    old.execute("INSERT INTO job_resources VALUES ('1.pbs', '2026-01-01')")
    old.commit()
    old.close()

    connection.connect(path).close()

    conn = sqlite3.connect(path)
    columns = {r[1] for r in conn.execute("PRAGMA table_info(job_resources)")}
    assert {"status", "tags", "username", "joblog_path"} <= columns
    assert conn.execute("SELECT status FROM job_resources").fetchone() == ("completed",)
    assert conn.execute("PRAGMA journal_mode").fetchone()[0] == "delete"
    assert conn.execute("PRAGMA user_version").fetchone()[0] == (
        connection.SCHEMA_VERSION
    )


# -- writes and retries -----------------------------------------------------------


def test_writes_start_with_begin_immediate(tmp_path, traced):
    path = tmp_path / "qxub.db"
    db.create_queue_entry("echo hi", db_path=path, pbs_job_id="1.pbs")
    traced.clear()

    db.mark_running("1.pbs", db_path=path)

    assert "BEGIN IMMEDIATE" in traced


def hold_lock(path, seconds):
    """Hold an exclusive lock on ``path`` from another connection."""
    holder = sqlite3.connect(path, check_same_thread=False)
    holder.execute("BEGIN EXCLUSIVE")

    def release():
        time.sleep(seconds)
        holder.rollback()
        holder.close()

    thread = threading.Thread(target=release)
    thread.start()
    return thread


def test_locked_write_is_retried(tmp_path, monkeypatch):
    path = tmp_path / "qxub.db"
    connection.connect(path).close()
    monkeypatch.setattr(connection, "RETRY_ATTEMPTS", 20)
    monkeypatch.setattr(connection, "RETRY_BASE_DELAY", 0.05)
    monkeypatch.setattr(connection, "RETRY_MAX_DELAY", 0.1)
    thread = hold_lock(path, 0.5)

    # A busy timeout far shorter than the lock: only the retries get it through.
    conn = connection.connect(path, timeout=0.01)
    with conn:
        conn.execute("INSERT INTO meta (key, value) VALUES ('k', 'v')")
    thread.join()

    assert sqlite3.connect(path).execute("SELECT value FROM meta").fetchone() == ("v",)


def test_retries_give_up_eventually(tmp_path, monkeypatch):
    path = tmp_path / "qxub.db"
    connection.connect(path).close()
    monkeypatch.setattr(connection, "RETRY_ATTEMPTS", 2)
    monkeypatch.setattr(connection, "RETRY_BASE_DELAY", 0.01)
    thread = hold_lock(path, 1.0)

    conn = connection.connect(path, timeout=0.01)
    with pytest.raises(sqlite3.OperationalError, match="locked"):
        conn.execute("INSERT INTO meta (key, value) VALUES ('k', 'v')")
    conn.close()
    thread.join()


@pytest.mark.parametrize(
    "message, in_transaction, retry",
    [
        ("database is locked", False, True),
        ("database is locked", True, True),  # the transaction is intact
        ("disk I/O error", False, True),
        ("disk I/O error", True, False),  # may have rolled back: don't resume
        ("no such table: queue", False, False),
    ],
)
def test_which_errors_are_retried(message, in_transaction, retry):
    exc = sqlite3.OperationalError(message)
    assert connection._should_retry(exc, in_transaction) is retry


# -- corruption -------------------------------------------------------------------


def test_corrupt_database_is_left_in_place_with_one_warning(tmp_path, caplog):
    """Review of #132: moving it aside cascaded under cross-node writers."""
    path = tmp_path / "qxub.db"
    junk = b"this is not a database" * 100
    path.write_bytes(junk)

    for _ in range(3):
        with pytest.raises(sqlite3.DatabaseError):
            connection.connect(path)

    assert [p.name for p in tmp_path.iterdir()] == ["qxub.db"]
    assert path.read_bytes() == junk
    warnings = [r for r in caplog.records if r.levelname == "WARNING"]
    assert len(warnings) == 1 and "is corrupt" in warnings[0].getMessage()


def test_failed_migration_does_not_stamp_the_version(tmp_path, monkeypatch):
    """Review of #132: a swallowed ALTER error marked the DB current."""
    path = tmp_path / "qxub.db"
    old = sqlite3.connect(path)
    old.execute(
        "CREATE TABLE queue (entry_id TEXT PRIMARY KEY, command_b64 TEXT,"
        " pbs_job_id TEXT, status TEXT, user TEXT)"
    )
    old.commit()
    old.close()
    monkeypatch.setattr(db, "_QUEUE_NEW_COLUMNS", [("not a valid name", "TEXT")])

    with pytest.raises(sqlite3.OperationalError):
        connection.connect(path)

    assert sqlite3.connect(path).execute("PRAGMA user_version").fetchone() == (0,)


# -- many processes ---------------------------------------------------------------

WORKER = textwrap.dedent(
    """
    import sys
    from qxub.queue import db
    from qxub.resources.tracker import ResourceTracker

    path, worker, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    tracker = ResourceTracker(db_path=path)
    for i in range(n):
        job = f"{worker}-{i}.pbs"
        db.create_queue_entry("echo", db_path=path, pbs_job_id=job)
        db.mark_running(job, db_path=path)
        tracker.update_job_status(job, "running")
        db.mark_complete(job, 0, db_path=path)
    """
)


def test_parallel_processes_record_every_job(tmp_path):
    path = tmp_path / "qxub.db"
    workers, per_worker = 8, 15

    procs = [
        subprocess.Popen(
            [sys.executable, "-c", WORKER, str(path), f"w{w}", str(per_worker)],
            stderr=subprocess.PIPE,
            text=True,
        )
        for w in range(workers)
    ]
    errors = [p.communicate(timeout=120)[1] for p in procs]
    assert [p.returncode for p in procs] == [0] * workers
    assert not [e for e in errors if "Failed" in e], errors

    conn = sqlite3.connect(path)
    rows = conn.execute("SELECT status, COUNT(*) FROM queue GROUP BY status").fetchall()
    assert rows == [("completed", workers * per_worker)]
