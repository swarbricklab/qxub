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
        # queue: tracks all qxub submissions with virtual IDs
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
                working_dir      TEXT,
                exec_context     TEXT,
                resources        TEXT,
                status           TEXT DEFAULT 'dispatched',
                pbs_job_id       TEXT,
                dispatched_at    TEXT,
                completed_at     TEXT,
                exit_code        INTEGER
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
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_project_jobs_status ON project_jobs(status)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_project_jobs_user ON project_jobs(user)"
        )


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
    pbs_job_id: str,
    command: str,
    tags: Optional[List[str]] = None,
    working_dir: Optional[str] = None,
    exec_context: Optional[Dict[str, Any]] = None,
    resources: Optional[Dict[str, Any]] = None,
    db_path: Optional[Path] = None,
) -> str:
    """Register a newly dispatched PBS job in the queue and project_jobs tables.

    Called by ``qxub exec`` immediately after a successful ``qsub``.

    Args:
        pbs_job_id:   Real PBS job ID returned by qsub.
        command:      Human-readable command string (stored base64-encoded).
        tags:         List of tag strings (e.g. ``["rule=align", "workflow=brca"]``).
        working_dir:  Working directory at submission time.
        exec_context: Dict describing execution context (type + env/mod/sif).
        resources:    Dict of requested resources (mem, cpus, walltime, etc.).
        db_path:      Override DB path (defaults to ``get_db_path()``).

    Returns:
        The virtual job ID string (``qx-{uuid4}``).
    """
    entry_id = new_virtual_id()
    user = os.environ.get("USER", "unknown")
    now = datetime.now().isoformat()
    command_b64 = base64.b64encode(command.encode("utf-8")).decode("ascii")

    try:
        # Ensure tables exist (no-op if already created)
        init_db(db_path)

        with get_connection(db_path) as conn:
            conn.execute(
                """
                INSERT INTO queue (
                    entry_id, submitted_at, user, tags, command_b64,
                    working_dir, exec_context, resources,
                    status, pbs_job_id, dispatched_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 'dispatched', ?, ?)
                """,
                (
                    entry_id,
                    now,
                    user,
                    json.dumps(tags or []),
                    command_b64,
                    working_dir,
                    json.dumps(exec_context) if exec_context else None,
                    json.dumps(resources) if resources else None,
                    pbs_job_id,
                    now,
                ),
            )
            conn.execute(
                """
                INSERT OR IGNORE INTO project_jobs (
                    pbs_job_id, entry_id, user, submitted_at, status, tags
                ) VALUES (?, ?, ?, ?, 'submitted', ?)
                """,
                (pbs_job_id, entry_id, user, now, json.dumps(tags or [])),
            )
    except Exception as exc:  # pylint: disable=broad-except
        logging.debug("Failed to create queue entry: %s", exc)

    return entry_id


def mark_running(pbs_job_id: str, db_path: Optional[Path] = None) -> None:
    """Mark a PBS job as running in ``project_jobs``.

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
    except Exception as exc:  # pylint: disable=broad-except
        logging.debug("mark_running failed for %s: %s", pbs_job_id, exc)


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
                SET status=?, completed_at=?, exit_code=?
                WHERE pbs_job_id=?
                """,
                (status, now, exit_code, pbs_job_id),
            )
    except Exception as exc:  # pylint: disable=broad-except
        logging.debug("mark_complete failed for %s: %s", pbs_job_id, exc)


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
        logging.debug("resolve_virtual_id failed for %s: %s", virtual_id, exc)
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
        logging.debug("get_project_job_status failed for %s: %s", pbs_job_id, exc)
        return None
