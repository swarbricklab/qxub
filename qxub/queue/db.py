"""
Unified database for qxub: queue, project job tracking, and metadata.

DB path resolution (priority order):
1. QXUB_DB_PATH environment variable (for testing/override)
2. shared_db.path from config  (system → user → project precedence)
3. ~/.config/qxub/qxub.db  (per-user fallback)

All tables live in a single SQLite file.  WAL mode is enabled for safe
concurrent access from multiple workflow nodes.
"""

import base64
import json
import logging

logger = logging.getLogger(__name__)
import os
import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------------
# DB path resolution
# ---------------------------------------------------------------------------


def get_db_path() -> Path:
    """Resolve the active qxub database path.

    Priority:
    1. ``QXUB_DB_PATH`` environment variable (testing / manual override)
    2. ``QXUB_SHARED_DB`` environment variable (set by job scripts on compute nodes)
    3. ``shared_db.path`` from the loaded config hierarchy
    4. ``~/.config/qxub/qxub.db`` per-user fallback
    """
    env_override = os.environ.get("QXUB_DB_PATH")
    if env_override:
        return Path(env_override)

    shared_db_env = os.environ.get("QXUB_SHARED_DB")
    if shared_db_env:
        return Path(shared_db_env)

    try:
        from ..config import config_manager  # pylint: disable=import-outside-toplevel

        shared_path = config_manager.get_config_value("shared_db.path")
        if shared_path:
            return Path(shared_path)
    except Exception:  # pylint: disable=broad-except
        pass  # Config not loaded yet or key absent — use fallback

    config_dir = Path.home() / ".config" / "qxub"
    config_dir.mkdir(parents=True, exist_ok=True)
    return config_dir / "qxub.db"


# ---------------------------------------------------------------------------
# Connection helper
# ---------------------------------------------------------------------------


@contextmanager
def get_connection(db_path: Optional[Path] = None):
    """Context manager yielding a SQLite connection in WAL mode.

    Commits on clean exit, rolls back on exception.
    """
    if db_path is None:
        db_path = get_db_path()

    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(str(db_path), timeout=5.0)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    conn.execute("PRAGMA mmap_size=0")
    conn.row_factory = sqlite3.Row

    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Schema initialisation / migration
# ---------------------------------------------------------------------------


def init_db(db_path: Optional[Path] = None) -> None:
    """Create or migrate the unified qxub database schema.

    This is idempotent — safe to call on every startup.

    Tables created:
    - ``queue``         — virtual job IDs + PBS mapping (Phase 2+)
    - ``project_jobs``  — PBS-level job tracking across users
    - ``meta``          — key/value store for DB-level state (active_count, etc.)
    """
    with get_connection(db_path) as conn:
        # ------------------------------------------------------------------
        # queue: unified job table — virtual ID as PK, PBS ID nullable
        # ------------------------------------------------------------------
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS queue (
                entry_id         TEXT PRIMARY KEY,
                submitted_at     TEXT NOT NULL,
                user             TEXT NOT NULL,
                tags             TEXT DEFAULT '[]',
                squashable       INTEGER DEFAULT 0,
                squash_condition TEXT,
                squash_batch_id  TEXT,
                command_b64      TEXT NOT NULL,
                command          TEXT,
                working_dir      TEXT,
                exec_context     TEXT,
                resources        TEXT,
                status           TEXT DEFAULT 'initiated',
                pbs_job_id       TEXT,
                initiated_at     TEXT,
                dispatched_at    TEXT,
                started_at       TEXT,
                completed_at     TEXT,
                last_status_update TEXT,
                exit_code        INTEGER,

                -- Identity / paths
                username         TEXT,
                joblog_path      TEXT,

                -- PBS execution details
                queue_name       TEXT,
                hostname         TEXT,

                -- Requested resources (standardised units)
                mem_requested_mb    INTEGER,
                time_requested_sec  INTEGER,
                cpus_requested      INTEGER,
                jobfs_requested_mb  INTEGER,

                -- Used resources
                mem_used_mb         INTEGER,
                time_used_sec       INTEGER,
                cpus_used           INTEGER,
                jobfs_used_mb       INTEGER,
                cpu_percent         REAL,

                -- Efficiency metrics (%%)
                mem_efficiency      REAL,
                time_efficiency     REAL,
                cpu_efficiency      REAL,
                jobfs_efficiency    REAL,

                -- Timing
                queue_wait_sec      INTEGER,
                execution_sec       INTEGER
            )
        """
        )

        # ------------------------------------------------------------------
        # project_jobs: PBS-level job tracking (all users on project)
        # ------------------------------------------------------------------
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS project_jobs (
                pbs_job_id          TEXT PRIMARY KEY,
                entry_id            TEXT,
                user                TEXT NOT NULL,
                submitted_at        TEXT NOT NULL,
                started_at          TEXT,
                completed_at        TEXT,
                status              TEXT DEFAULT 'submitted',
                exit_code           INTEGER,
                tags                TEXT DEFAULT '[]',
                resources_requested TEXT,
                resources_used      TEXT
            )
        """
        )

        # ------------------------------------------------------------------
        # meta: key-value store (e.g. active_count)
        # ------------------------------------------------------------------
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS meta (
                key        TEXT PRIMARY KEY,
                value      TEXT,
                updated_at TEXT
            )
        """
        )

        # Useful indexes
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_queue_pbs_job_id ON queue(pbs_job_id)"
        )
        conn.execute("CREATE INDEX IF NOT EXISTS idx_queue_status ON queue(status)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_queue_user ON queue(user)")
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_project_jobs_status ON project_jobs(status)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_project_jobs_user ON project_jobs(user)"
        )

        # Schema migration for existing queue tables
        _migrate_queue_schema(conn)


# ---------------------------------------------------------------------------
# Queue schema migration
# ---------------------------------------------------------------------------

# Columns added to the ``queue`` table for the unified-DB merge (issue #63).
_QUEUE_NEW_COLUMNS = [
    ("command", "TEXT"),
    ("initiated_at", "TEXT"),
    ("started_at", "TEXT"),
    ("last_status_update", "TEXT"),
    ("username", "TEXT"),
    ("joblog_path", "TEXT"),
    ("queue_name", "TEXT"),
    ("hostname", "TEXT"),
    ("mem_requested_mb", "INTEGER"),
    ("time_requested_sec", "INTEGER"),
    ("cpus_requested", "INTEGER"),
    ("jobfs_requested_mb", "INTEGER"),
    ("mem_used_mb", "INTEGER"),
    ("time_used_sec", "INTEGER"),
    ("cpus_used", "INTEGER"),
    ("jobfs_used_mb", "INTEGER"),
    ("cpu_percent", "REAL"),
    ("mem_efficiency", "REAL"),
    ("time_efficiency", "REAL"),
    ("cpu_efficiency", "REAL"),
    ("jobfs_efficiency", "REAL"),
    ("queue_wait_sec", "INTEGER"),
    ("execution_sec", "INTEGER"),
]


def _migrate_queue_schema(conn) -> None:
    """Add new resource-tracking columns to an existing ``queue`` table.

    Idempotent — safe to run on every startup.
    """
    try:
        cursor = conn.execute("PRAGMA table_info(queue)")
        existing = {row[1] for row in cursor.fetchall()}

        for col_name, col_type in _QUEUE_NEW_COLUMNS:
            if col_name not in existing:
                conn.execute(f"ALTER TABLE queue ADD COLUMN {col_name} {col_type}")
                logger.debug("queue: added column %s %s", col_name, col_type)

        # Backfill: change old default status 'dispatched' to keep working
        # (new default is 'initiated' but existing rows already have a status)
        conn.commit()
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("Queue schema migration failed (may be normal): %s", exc)


# ---------------------------------------------------------------------------
# Virtual job ID helpers
# ---------------------------------------------------------------------------


def new_virtual_id() -> str:
    """Generate a new qxub virtual job ID (``qx-{uuid4}``)."""
    return f"qx-{uuid.uuid4()}"


def is_virtual_id(job_id: str) -> bool:
    """Return True if *job_id* is a qxub virtual ID (starts with ``qx-``)."""
    return isinstance(job_id, str) and job_id.startswith("qx-")


# ---------------------------------------------------------------------------
# Queue entry lifecycle
# ---------------------------------------------------------------------------


def create_queue_entry(
    command: str,
    pbs_job_id: Optional[str] = None,
    tags: Optional[List[str]] = None,
    working_dir: Optional[str] = None,
    exec_context: Optional[Dict[str, Any]] = None,
    resources: Optional[Dict[str, Any]] = None,
    joblog_path: Optional[str] = None,
    queue_name: Optional[str] = None,
    cpus_requested: Optional[int] = None,
    db_path: Optional[Path] = None,
) -> str:
    """Register a job in the unified queue table.

    Must be called **before** ``qsub()`` so that a stable virtual ID exists
    even while PBS is blocking.  The ``pbs_job_id`` field is populated later
    via :func:`update_queue_entry`.

    Args:
        command:          Human-readable command string (stored base64-encoded).
        pbs_job_id:       Real PBS job ID (``None`` when called pre-qsub).
        tags:             List of tag strings.
        working_dir:      Working directory at submission time.
        exec_context:     Dict describing execution context (type + env/mod/sif).
        resources:        Dict of requested resources (mem, cpus, walltime, etc.).
        joblog_path:      Expected path to the PBS joblog file.
        queue_name:       PBS queue name (e.g. ``normal``, ``express``).
        cpus_requested:   Number of CPUs requested.
        db_path:          Override DB path (defaults to ``get_db_path()``).

    Returns:
        The virtual job ID string (``qx-{uuid4}``).
    """
    entry_id = new_virtual_id()
    user = os.environ.get("USER", "unknown")
    now = datetime.now().isoformat()
    command_b64 = base64.b64encode(command.encode("utf-8")).decode("ascii")

    # Status depends on whether we already have a PBS job ID
    status = "dispatched" if pbs_job_id else "initiated"

    try:
        # Ensure tables exist (no-op if already created)
        init_db(db_path)

        with get_connection(db_path) as conn:
            conn.execute(
                """
                INSERT INTO queue (
                    entry_id, submitted_at, user, tags, command_b64, command,
                    working_dir, exec_context, resources,
                    status, pbs_job_id, initiated_at, dispatched_at,
                    last_status_update, username, joblog_path,
                    queue_name, cpus_requested
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    entry_id,
                    now,
                    user,
                    json.dumps(tags or []),
                    command_b64,
                    command,
                    working_dir,
                    json.dumps(exec_context) if exec_context else None,
                    json.dumps(resources) if resources else None,
                    status,
                    pbs_job_id,
                    now,  # initiated_at
                    now if pbs_job_id else None,  # dispatched_at
                    now,  # last_status_update
                    user,  # username
                    joblog_path,
                    queue_name,
                    cpus_requested,
                ),
            )

            # Also insert into project_jobs if we have a PBS ID
            if pbs_job_id:
                conn.execute(
                    """
                    INSERT OR IGNORE INTO project_jobs (
                        pbs_job_id, entry_id, user, submitted_at, status, tags
                    ) VALUES (?, ?, ?, ?, 'submitted', ?)
                    """,
                    (pbs_job_id, entry_id, user, now, json.dumps(tags or [])),
                )
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("Failed to create queue entry: %s", exc)

    return entry_id


def get_queue_entry(pbs_job_id: str, db_path: Optional[Path] = None) -> Optional[dict]:
    """Look up a queue entry by PBS job ID.

    Returns a dict with the queue row columns, or *None* if not found.
    Failures are silently swallowed (returns *None*).
    """
    try:
        with get_connection(db_path) as conn:
            conn.row_factory = sqlite3.Row
            row = conn.execute(
                "SELECT * FROM queue WHERE pbs_job_id=?", (pbs_job_id,)
            ).fetchone()
            return dict(row) if row else None
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("get_queue_entry failed for %s: %s", pbs_job_id, exc)
        return None


def update_queue_entry(
    virtual_id: str,
    pbs_job_id: Optional[str] = None,
    status: Optional[str] = None,
    joblog_path: Optional[str] = None,
    db_path: Optional[Path] = None,
) -> bool:
    """Update an existing queue entry after ``qsub()`` returns.

    Typically called to fill in the real PBS job ID and transition the
    status from ``initiated`` → ``dispatched``.  Also used to record
    ``failed`` if ``qsub()`` raises an exception.

    Args:
        virtual_id:   The ``qx-{uuid4}`` returned by :func:`create_queue_entry`.
        pbs_job_id:   Real PBS job ID (or ``None`` on failure).
        status:       New status string (``dispatched`` or ``failed``).
        joblog_path:  Path to the PBS joblog (may only be known after qsub).
        db_path:      Override DB path.

    Returns:
        ``True`` if the row was updated, ``False`` otherwise.
    """
    now = datetime.now().isoformat()
    try:
        with get_connection(db_path) as conn:
            set_parts = ["last_status_update=?"]
            params: list = [now]

            if pbs_job_id is not None:
                set_parts.append("pbs_job_id=?")
                params.append(pbs_job_id)
                set_parts.append("dispatched_at=?")
                params.append(now)

            if status is not None:
                set_parts.append("status=?")
                params.append(status)

            if joblog_path is not None:
                set_parts.append("joblog_path=?")
                params.append(joblog_path)

            params.append(virtual_id)
            result = conn.execute(
                f"UPDATE queue SET {', '.join(set_parts)} WHERE entry_id=?",
                params,
            )

            # If we now have a PBS ID, also insert into project_jobs
            if pbs_job_id and result.rowcount > 0:
                row = conn.execute(
                    "SELECT user, submitted_at, tags FROM queue WHERE entry_id=?",
                    (virtual_id,),
                ).fetchone()
                if row:
                    conn.execute(
                        """
                        INSERT OR IGNORE INTO project_jobs (
                            pbs_job_id, entry_id, user, submitted_at, status, tags
                        ) VALUES (?, ?, ?, ?, 'submitted', ?)
                        """,
                        (pbs_job_id, virtual_id, row[0], row[1], row[2]),
                    )

            return result.rowcount > 0
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("update_queue_entry failed for %s: %s", virtual_id, exc)
        return False


def get_entries_bulk(
    virtual_ids: List[str],
    db_path: Optional[Path] = None,
) -> Dict[str, Dict[str, Any]]:
    """Get status for multiple queue entries in a single SQL query.

    Args:
        virtual_ids:  List of virtual job IDs (``qx-{uuid4}``).
        db_path:      Override DB path.

    Returns:
        Dict mapping ``virtual_id`` → entry dict.  IDs not found are omitted.
    """
    if not virtual_ids:
        return {}

    try:
        with get_connection(db_path) as conn:
            conn.row_factory = sqlite3.Row
            placeholders = ",".join("?" for _ in virtual_ids)
            rows = conn.execute(
                f"""
                SELECT entry_id, status, pbs_job_id, exit_code,
                       submitted_at, initiated_at, dispatched_at,
                       started_at, completed_at, tags, command,
                       working_dir, queue_name
                FROM queue
                WHERE entry_id IN ({placeholders})
                """,
                virtual_ids,
            ).fetchall()
            return {row["entry_id"]: dict(row) for row in rows}
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("get_entries_bulk failed: %s", exc)
        return {}


def mark_running(pbs_job_id: str, db_path: Optional[Path] = None) -> None:
    """Mark a PBS job as running in ``project_jobs`` and ``queue``.

    Called from within the running job script (via ``QXUB_SHARED_DB``).
    Failures are silently swallowed — status files remain the authoritative
    fallback.
    """
    now = datetime.now().isoformat()
    try:
        with get_connection(db_path) as conn:
            conn.execute(
                "UPDATE project_jobs SET status='running', started_at=? WHERE pbs_job_id=?",
                (now, pbs_job_id),
            )
            conn.execute(
                """
                UPDATE queue
                SET status='running', started_at=?, last_status_update=?
                WHERE pbs_job_id=?
                """,
                (now, now, pbs_job_id),
            )
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("mark_running failed for %s: %s", pbs_job_id, exc)


def mark_complete(
    pbs_job_id: str,
    exit_code: int,
    db_path: Optional[Path] = None,
) -> None:
    """Mark a PBS job as completed or failed in ``project_jobs`` and ``queue``.

    Called from within the running job script (via ``QXUB_SHARED_DB``).
    Failures are silently swallowed.
    """
    now = datetime.now().isoformat()
    status = "completed" if exit_code == 0 else "failed"
    try:
        with get_connection(db_path) as conn:
            conn.execute(
                """
                UPDATE project_jobs
                SET status=?, completed_at=?, exit_code=?
                WHERE pbs_job_id=?
                """,
                (status, now, exit_code, pbs_job_id),
            )
            conn.execute(
                """
                UPDATE queue
                SET status=?, completed_at=?, exit_code=?, last_status_update=?
                WHERE pbs_job_id=?
                """,
                (status, now, exit_code, now, pbs_job_id),
            )
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("mark_complete failed for %s: %s", pbs_job_id, exc)


# ---------------------------------------------------------------------------
# Virtual ID resolution
# ---------------------------------------------------------------------------


def resolve_virtual_id(
    virtual_id: str,
    db_path: Optional[Path] = None,
) -> Optional[Dict[str, Any]]:
    """Resolve a virtual job ID to its queue entry dict.

    Returns ``None`` if the ID is not found.

    Typical usage::

        entry = resolve_virtual_id("qx-550e8400-...")
        if entry and entry["pbs_job_id"]:
            real_id = entry["pbs_job_id"]
    """
    try:
        with get_connection(db_path) as conn:
            row = conn.execute(
                """
                SELECT entry_id, status, pbs_job_id, exit_code,
                       submitted_at, completed_at, tags
                FROM queue
                WHERE entry_id = ?
                """,
                (virtual_id,),
            ).fetchone()
            return dict(row) if row else None
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("resolve_virtual_id failed for %s: %s", virtual_id, exc)
        return None


def get_project_job_status(
    pbs_job_id: str,
    db_path: Optional[Path] = None,
) -> Optional[Dict[str, Any]]:
    """Return the ``project_jobs`` row for a given PBS job ID, or ``None``."""
    try:
        with get_connection(db_path) as conn:
            row = conn.execute(
                """
                SELECT pbs_job_id, entry_id, status, exit_code,
                       submitted_at, started_at, completed_at, tags
                FROM project_jobs
                WHERE pbs_job_id = ?
                """,
                (pbs_job_id,),
            ).fetchone()
            return dict(row) if row else None
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("get_project_job_status failed for %s: %s", pbs_job_id, exc)
        return None
