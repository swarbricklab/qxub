"""SQLite connections that behave on shared filesystems (Lustre, GPFS, NFS).

Phase 1 of ``docs/dev/sqlite-robustness.md``. Every qxub connection to the
job database comes from :func:`connect`, which:

- disables mmap (``PRAGMA mmap_size=0``, #69) and sets a busy timeout;
- starts write transactions with ``BEGIN IMMEDIATE``, so a writer takes its
  lock up front. A deferred transaction that later upgrades from read to
  write can fail with ``SQLITE_BUSY`` straight away, which the busy timeout
  does not retry;
- retries a statement or commit that fails with a transient error, with
  jittered backoff (:class:`RetryingConnection`);
- creates or migrates the schema only when ``PRAGMA user_version`` is behind
  :data:`SCHEMA_VERSION`, once per process (:func:`ensure_schema`). Before,
  every connection ran DDL, and each job opened the DB about four times;
- moves a corrupt database aside so the next connection starts a fresh one
  (:func:`quarantine_if_corrupt`).
"""

import logging
import os
import random
import sqlite3
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import Optional, Set

logger = logging.getLogger(__name__)

# Bump when a table or column is added; ensure_schema() then runs the
# (idempotent) DDL and migrations once on each existing database.
SCHEMA_VERSION = 1

BUSY_TIMEOUT_SECONDS = 15.0
RETRY_ATTEMPTS = 5
RETRY_BASE_DELAY = 0.2
RETRY_MAX_DELAY = 3.0

# Messages of sqlite3.OperationalError worth retrying.
_BUSY_MESSAGES = ("database is locked", "database is busy", "database table is locked")
_IO_MESSAGES = ("disk i/o error", "unable to open database file")
_CORRUPT_MESSAGES = ("database disk image is malformed", "file is not a database")

_schema_ready: Set[str] = set()
_schema_lock = threading.Lock()


def _message(exc: BaseException) -> str:
    return str(exc).lower()


def is_busy(exc: BaseException) -> bool:
    return isinstance(exc, sqlite3.OperationalError) and any(
        m in _message(exc) for m in _BUSY_MESSAGES
    )


def is_io_error(exc: BaseException) -> bool:
    return isinstance(exc, sqlite3.OperationalError) and any(
        m in _message(exc) for m in _IO_MESSAGES
    )


def is_corrupt(exc: BaseException) -> bool:
    return isinstance(exc, sqlite3.DatabaseError) and any(
        m in _message(exc) for m in _CORRUPT_MESSAGES
    )


def _should_retry(exc: BaseException, was_in_transaction: bool) -> bool:
    """Whether a failed statement can safely be run again.

    A busy/locked error leaves the transaction as it was, so retrying the
    statement (or the commit) is safe. An I/O error may have rolled the
    transaction back, so it is retried only when no transaction was open
    before the statement, i.e. when nothing earlier would be lost.
    """
    if is_busy(exc):
        return True
    return is_io_error(exc) and not was_in_transaction


def _backoff(attempt: int) -> float:
    """Full-jitter exponential backoff, so waiting writers spread out."""
    return random.uniform(0, min(RETRY_MAX_DELAY, RETRY_BASE_DELAY * 2**attempt))


class RetryingConnection(sqlite3.Connection):
    """A connection whose statements and commits retry transient errors."""

    def _retry(self, op, *args):
        for attempt in range(RETRY_ATTEMPTS):
            was_in_transaction = self.in_transaction
            try:
                return op(*args)
            except sqlite3.DatabaseError as exc:
                if is_corrupt(exc) and getattr(self, "db_path", None):
                    quarantine_if_corrupt(self.db_path)
                    raise
                if attempt == RETRY_ATTEMPTS - 1 or not _should_retry(
                    exc, was_in_transaction
                ):
                    raise
                if self.in_transaction and not was_in_transaction:
                    # The implicit BEGIN succeeded but the statement failed.
                    super().rollback()
                delay = _backoff(attempt)
                logger.debug(
                    "SQLite %s (attempt %d), retrying in %.2fs", exc, attempt + 1, delay
                )
                time.sleep(delay)
        raise AssertionError("unreachable")

    def execute(self, *args):
        return self._retry(super().execute, *args)

    def executemany(self, *args):
        return self._retry(super().executemany, *args)

    def commit(self):
        return self._retry(super().commit)

    def __exit__(self, exc_type, exc, tb):
        # sqlite3's own __exit__ calls the C commit, bypassing the retry.
        if exc_type is None:
            self.commit()
        else:
            self.rollback()
        return False


def connect(db_path, timeout: float = BUSY_TIMEOUT_SECONDS) -> RetryingConnection:
    """Open the job database with the settings above; the schema is ready."""
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = _open(db_path, timeout)
    try:
        ensure_schema(conn, db_path)
    except Exception:
        conn.close()
        raise
    return conn


def _open(db_path: Path, timeout: float) -> RetryingConnection:
    conn = sqlite3.connect(
        str(db_path),
        timeout=timeout,
        isolation_level="IMMEDIATE",  # writes start with BEGIN IMMEDIATE
        factory=RetryingConnection,
    )
    conn.db_path = db_path
    conn.execute(f"PRAGMA busy_timeout={int(timeout * 1000)}")
    conn.execute("PRAGMA mmap_size=0")
    return conn


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


def ensure_schema(conn: sqlite3.Connection, db_path: Path) -> None:
    """Create or migrate the schema if this database's version is behind.

    The version check is a single read; the DDL (which needs a write lock)
    runs only on a new or older database, and at most once per process.
    """
    key = str(Path(db_path).resolve())
    if key in _schema_ready:
        return
    with _schema_lock:
        if key in _schema_ready:
            return
        try:
            if _user_version(conn) < SCHEMA_VERSION:
                _upgrade(conn)
        except sqlite3.DatabaseError as exc:
            if is_corrupt(exc):
                quarantine_if_corrupt(db_path)
            raise
        _schema_ready.add(key)


def _user_version(conn: sqlite3.Connection) -> int:
    return conn.execute("PRAGMA user_version").fetchone()[0]


def _upgrade(conn: sqlite3.Connection) -> None:
    # WAL is persistent in the file and its -shm is always mmap'd (#75).
    # The journal mode cannot change inside a transaction, so do it first.
    mode = conn.execute("PRAGMA journal_mode").fetchone()[0]
    if str(mode).lower() != "delete":
        conn.execute("PRAGMA journal_mode=DELETE")

    from ..resources.tracker import create_job_resources_schema
    from .db import create_queue_schema

    conn.execute("BEGIN IMMEDIATE")
    try:
        # Another process may have upgraded while we waited for the lock.
        if _user_version(conn) < SCHEMA_VERSION:
            create_queue_schema(conn)
            create_job_resources_schema(conn)
            conn.execute(f"PRAGMA user_version={SCHEMA_VERSION}")
        conn.commit()
    except BaseException:
        conn.rollback()
        raise


def forget_schema_cache() -> None:
    """For tests: make the next connection check the schema version again."""
    _schema_ready.clear()


# ---------------------------------------------------------------------------
# Corruption
# ---------------------------------------------------------------------------


def quarantine_if_corrupt(db_path) -> Optional[Path]:
    """Move a corrupt database aside so the next connection starts afresh.

    Called after an operation fails with "malformed" or "not a database". A
    torn read on a shared filesystem can raise the same error on a healthy
    file, so the file is checked first with ``PRAGMA quick_check`` and only
    moved if that fails too. The rollback journal is moved with it: a hot
    journal left behind would be applied to the new database.

    Returns the new path of the corrupt file, or None if it was not moved.
    """
    db_path = Path(db_path)
    try:
        before = db_path.stat()
    except FileNotFoundError:
        return None
    if _passes_quick_check(db_path):
        return None
    try:
        after = db_path.stat()
    except FileNotFoundError:
        return None
    if (after.st_ino, after.st_mtime_ns) != (before.st_ino, before.st_mtime_ns):
        return None  # replaced or written meanwhile; let the next failure decide

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    target = db_path.with_name(f"{db_path.name}.corrupt-{stamp}-{os.getpid()}")
    try:
        os.replace(db_path, target)
    except FileNotFoundError:
        return None  # another process moved it first
    for suffix in ("-journal", "-wal", "-shm"):
        side = db_path.with_name(db_path.name + suffix)
        try:
            os.replace(side, target.with_name(target.name + suffix))
        except FileNotFoundError:
            pass
    _schema_ready.discard(str(db_path.resolve()))
    logger.warning(
        "qxub job database %s was corrupt; moved it to %s and started a new one",
        db_path,
        target,
    )
    return target


def _passes_quick_check(db_path: Path) -> bool:
    try:
        conn = sqlite3.connect(str(db_path), timeout=BUSY_TIMEOUT_SECONDS)
        try:
            conn.execute("PRAGMA mmap_size=0")
            rows = conn.execute("PRAGMA quick_check").fetchall()
        finally:
            conn.close()
    except sqlite3.DatabaseError as exc:
        if is_corrupt(exc):
            return False
        return True  # busy or I/O trouble: not evidence of corruption
    return rows == [("ok",)]
