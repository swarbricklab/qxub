"""
Regression tests for the qxub queue database schema handling.

Guards against the silent-tracking-failure bug where a shared DB created under
the Phase-2 schema (``queue.virtual_id`` PK) survived the unified-DB merge
(``queue.entry_id`` PK).  ``init_db`` could not migrate the legacy table, so
every ``create_queue_entry`` insert failed with ``no such column`` — silently,
because the write path swallows exceptions.  See ``qxub/queue/db.py``.
"""

import sqlite3

import pytest

from qxub.queue import db


# Minimal reproduction of the Phase-2 (#55) ``queue`` schema.
_LEGACY_QUEUE_DDL = """
    CREATE TABLE queue (
        virtual_id TEXT PRIMARY KEY,
        pbs_job_id TEXT,
        job_name   TEXT,
        status     TEXT DEFAULT 'initiated',
        exit_code  INTEGER,
        created_at TEXT,
        updated_at TEXT,
        ncpus      INTEGER,
        mem        TEXT,
        walltime   TEXT,
        queue      TEXT,
        project    TEXT,
        tags       TEXT,
        jobfs      TEXT
    );
    CREATE INDEX idx_queue_pbs_job_id ON queue(pbs_job_id);
    CREATE INDEX idx_queue_status ON queue(status);
"""


def _make_legacy_db(path, rows=None):
    conn = sqlite3.connect(str(path))
    conn.executescript(_LEGACY_QUEUE_DDL)
    for r in rows or []:
        conn.execute(
            "INSERT INTO queue (virtual_id, pbs_job_id, job_name, status, created_at)"
            " VALUES (?, ?, ?, ?, ?)",
            r,
        )
    conn.commit()
    conn.close()


def _columns(path, table="queue"):
    conn = sqlite3.connect(str(path))
    try:
        return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
    finally:
        conn.close()


def test_fresh_db_gets_unified_schema(tmp_path):
    """A brand-new DB should be created with the unified ``entry_id`` schema."""
    dbp = tmp_path / "fresh.db"
    db.init_db(db_path=dbp)
    cols = _columns(dbp)
    assert "entry_id" in cols
    assert "submitted_at" in cols
    assert "command_b64" in cols


def test_empty_legacy_queue_is_rebuilt(tmp_path):
    """An empty legacy queue table is rebuilt to the unified schema and usable."""
    dbp = tmp_path / "legacy_empty.db"
    _make_legacy_db(dbp)

    db.init_db(db_path=dbp)

    cols = _columns(dbp)
    assert "entry_id" in cols and "virtual_id" not in cols

    # The empty legacy table is dropped, not backed up.
    conn = sqlite3.connect(str(dbp))
    names = {r[0] for r in conn.execute("SELECT name FROM sqlite_master")}
    conn.close()
    assert "queue_legacy_backup" not in names

    # And a real submission now records a row (the original bug: 0 rows).
    vid = db.create_queue_entry(
        "echo hi", pbs_job_id="1.gadi-pbs", tags=["t"], db_path=dbp
    )
    entry = db.get_queue_entry("1.gadi-pbs", db_path=dbp)
    assert entry is not None
    assert entry["entry_id"] == vid
    assert entry["command"] == "echo hi"


def test_nonempty_legacy_queue_is_backed_up(tmp_path):
    """A legacy queue with rows is preserved as ``queue_legacy_backup``."""
    dbp = tmp_path / "legacy_rows.db"
    _make_legacy_db(dbp, rows=[("qx-old", "111.gadi", "oldjob", "running", "2026-01-01")])

    db.init_db(db_path=dbp)

    assert "entry_id" in _columns(dbp)
    conn = sqlite3.connect(str(dbp))
    try:
        backup = conn.execute("SELECT COUNT(*) FROM queue_legacy_backup").fetchone()[0]
    finally:
        conn.close()
    assert backup == 1


def test_init_db_is_idempotent_on_unified_db(tmp_path):
    """Re-running init_db on an already-unified DB must not rebuild or warn."""
    dbp = tmp_path / "unified.db"
    db.init_db(db_path=dbp)
    db.create_queue_entry("echo hi", pbs_job_id="2.gadi-pbs", db_path=dbp)

    db.init_db(db_path=dbp)  # second call — should be a no-op

    conn = sqlite3.connect(str(dbp))
    try:
        n = conn.execute("SELECT COUNT(*) FROM queue").fetchone()[0]
        names = {r[0] for r in conn.execute("SELECT name FROM sqlite_master")}
    finally:
        conn.close()
    assert n == 1
    assert "queue_legacy_backup" not in names


def test_full_lifecycle_on_healed_legacy_db(tmp_path):
    """create -> mark_running -> mark_complete works after healing."""
    dbp = tmp_path / "lifecycle.db"
    _make_legacy_db(dbp)

    vid = db.create_queue_entry(
        "run job", pbs_job_id="3.gadi-pbs", tags=["x"], db_path=dbp
    )
    db.mark_running("3.gadi-pbs", db_path=dbp)
    db.mark_complete("3.gadi-pbs", 0, db_path=dbp)

    entry = db.resolve_virtual_id(vid, db_path=dbp)
    assert entry is not None
    assert entry["status"] == "completed"
    assert entry["exit_code"] == 0
