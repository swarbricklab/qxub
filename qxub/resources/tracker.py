"""
Simple resource tracking system for qxub.

Tracks only key metrics: job_id, requested vs used resources, efficiency.
Uses SQLite for simple querying and analysis.
"""

import json
import logging

logger = logging.getLogger(__name__)
import os
import re
import sqlite3
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from ..queue.connection import BUSY_TIMEOUT_SECONDS, connect
from .parser import size_to_bytes, time_to_seconds


def parse_date_filter(value: str) -> str:
    """Parse a date filter string into an ISO datetime string.

    Accepts relative offsets (``7d``, ``2h``, ``1w``, ``30m``) or ISO date /
    datetime strings (``2026-02-20``, ``2026-02-20T10:00:00``).
    """
    value = value.strip()
    m = re.fullmatch(r"(\d+)([dhwm])", value, re.IGNORECASE)
    if m:
        n, unit = int(m.group(1)), m.group(2).lower()
        delta = {
            "d": timedelta(days=n),
            "h": timedelta(hours=n),
            "w": timedelta(weeks=n),
            "m": timedelta(minutes=n),
        }[unit]
        return (datetime.now() - delta).isoformat(sep=" ", timespec="seconds")
    try:
        dt = datetime.fromisoformat(value)
        return dt.isoformat(sep=" ", timespec="seconds")
    except ValueError:
        raise ValueError(
            f"Cannot parse date '{value}'. Use e.g. '7d', '2h', '1w', '2026-02-20'."
        )


def _resolve_tracker_db_path() -> Path:
    """Resolve the default database path for ResourceTracker.

    Delegates to ``qxub.queue.db.get_db_path()`` so the resource tracker
    always writes to the same database as the queue system.  Also performs
    a one-time migration of the legacy ``resources.db`` → ``qxub.db`` for
    existing user installs.
    """
    from ..queue.db import get_db_path  # pylint: disable=import-outside-toplevel

    resolved = get_db_path()

    # One-time migration: copy resources.db → qxub.db for existing installs
    config_dir = Path.home() / ".config" / "qxub"
    old_path = config_dir / "resources.db"
    if old_path.exists() and not resolved.exists():
        import shutil  # pylint: disable=import-outside-toplevel

        try:
            shutil.copy2(str(old_path), str(resolved))
            logger.debug("Migrated resources.db → %s", resolved)
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("Could not migrate resources.db: %s", exc)

    return resolved


# Busy timeout for each attempt; RetryingConnection retries on top of it.
_SQLITE_TIMEOUT = BUSY_TIMEOUT_SECONDS
# Status reads (qxtat check, polled by Snakemake) retry here three times on
# top of the connection's own retries, so each attempt waits less.
_STATUS_READ_TIMEOUT = 5.0


def create_job_resources_schema(conn) -> None:
    """Create or migrate the ``job_resources`` table. Idempotent; no commit.

    Called by ``qxub.queue.connection.ensure_schema`` inside its transaction.
    """
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS job_resources (
            job_id TEXT PRIMARY KEY,
            timestamp TEXT NOT NULL,
            command TEXT,
            exit_code INTEGER,

            -- Job status tracking (new for v3)
            status TEXT DEFAULT 'submitted',
            submitted_at TEXT,
            started_at TEXT,
            completed_at TEXT,
            last_status_update TEXT,

            -- Requested resources
            mem_requested_mb INTEGER,
            time_requested_sec INTEGER,
            cpus_requested INTEGER,
            jobfs_requested_mb INTEGER,

            -- Used resources
            mem_used_mb INTEGER,
            time_used_sec INTEGER,
            cpus_used INTEGER,
            jobfs_used_mb INTEGER,
            cpu_percent REAL,

            -- Efficiency metrics (%)
            mem_efficiency REAL,
            time_efficiency REAL,
            cpu_efficiency REAL,
            jobfs_efficiency REAL,

            -- Execution details
            queue TEXT,
            hostname TEXT,

            -- Timing
            queue_wait_sec INTEGER,
            execution_sec INTEGER,

            -- Tags
            tags TEXT DEFAULT '[]',

            -- Identity
            username TEXT,

            -- Path to PBS joblog for deferred resource backfill
            joblog_path TEXT,

            -- Scheduler job name (-N), so cancel can match it without PBS
            job_name TEXT
        )
        """
    )
    columns = {row[1] for row in conn.execute("PRAGMA table_info(job_resources)")}
    if "status" not in columns:
        logger.debug("Migrating job_resources to add status tracking columns")
        conn.execute(
            "ALTER TABLE job_resources ADD COLUMN status TEXT DEFAULT 'completed'"
        )
        for column in (
            "submitted_at",
            "started_at",
            "completed_at",
            "last_status_update",
        ):
            conn.execute(f"ALTER TABLE job_resources ADD COLUMN {column} TEXT")
        conn.execute(
            """
            UPDATE job_resources
            SET status='completed',
                completed_at=timestamp,
                last_status_update=timestamp
            WHERE status IS NULL OR status = 'completed'
            """
        )
    for column, column_type in (
        ("tags", "TEXT DEFAULT '[]'"),
        ("username", "TEXT"),
        ("joblog_path", "TEXT"),
        ("job_name", "TEXT"),
    ):
        if column not in columns:
            conn.execute(f"ALTER TABLE job_resources ADD COLUMN {column} {column_type}")
    # An old writer stored one status JSON-quoted ('"running"').
    conn.execute(
        "UPDATE job_resources SET status=TRIM(status, '\"') WHERE status LIKE '\"%'"
    )


def _snapshot_outcomes(job_ids: Sequence[str]) -> Dict[str, tuple]:
    """``{job_id: (status, exit_code, completed_at)}`` for finished jobs in
    the statistics snapshot. Empty if there is no snapshot or it can't be read.
    """
    from ..queue.consolidate import (  # pylint: disable=import-outside-toplevel
        open_snapshot,
    )

    try:
        conn = open_snapshot()
    except (OSError, sqlite3.Error):
        return {}
    found: Dict[str, tuple] = {}
    try:
        ids = list(job_ids)
        for i in range(0, len(ids), 500):
            chunk = ids[i : i + 500]
            rows = conn.execute(
                "SELECT job_id, status, exit_code, completed_at FROM job_resources"
                f" WHERE job_id IN ({','.join('?' * len(chunk))})"
                " AND status IN ('completed', 'failed', 'cancelled')",
                chunk,
            )
            for job_id, status, code, finished in rows:
                found[job_id] = (status, code, finished)
    except sqlite3.Error as exc:
        logger.debug("Reading the snapshot failed: %s", exc)
    finally:
        conn.close()
    return found


class DatabaseError(Exception):
    """Raised when a sqlite query fails due to a transient or I/O error.

    Distinct from a *missing row* (which returns ``None``).  Callers can
    catch this to distinguish "DB unreachable" from "job not found".
    """


class ResourceTracker:
    """Simple resource tracking focused on efficiency metrics."""

    def __init__(self, db_path: Optional[Path] = None):
        """Initialize resource tracker with SQLite database."""
        if db_path is None:
            db_path = _resolve_tracker_db_path()

        self.db_path = Path(db_path)
        self._init_database()

    def _connect(self, timeout: float = _SQLITE_TIMEOUT):
        """A connection to the job database (see ``qxub.queue.connection``).

        mmap is disabled, writes start with ``BEGIN IMMEDIATE``, transient
        errors are retried, and the schema is created if out of date.
        """
        return connect(self.db_path, timeout=timeout)

    def _init_database(self):
        """Make sure the schema exists (a single version read once it does)."""
        self._connect().close()

    def log_job_resources(
        self, job_id: str, resource_data: Dict[str, Any], command: Optional[str] = None
    ) -> bool:
        """
        Log resource usage for a completed job.

        Args:
            job_id: PBS job ID
            resource_data: Resource data from get_job_resource_data()
            command: Command that was executed (optional)

        Returns:
            bool: True if logged successfully
        """
        try:
            # Extract key metrics only
            requested = resource_data.get("resources_requested", {})
            used = resource_data.get("resources_used", {})
            efficiency = resource_data.get("efficiency", {})
            execution = resource_data.get("execution", {})
            timing = resource_data.get("timing", {})

            # Clean up command for better display
            clean_command = self._clean_command(command)

            # Convert to standardized units
            record = {
                "job_id": job_id,
                "timestamp": datetime.now().isoformat(),
                "command": clean_command,
                "exit_code": resource_data.get("exit_status"),
                # Requested resources (convert to MB/seconds)
                "mem_requested_mb": self._parse_size_to_mb(requested.get("mem")),
                "time_requested_sec": self._parse_time_to_sec(
                    requested.get("walltime")
                ),
                "cpus_requested": self._parse_int(requested.get("ncpus")),
                "jobfs_requested_mb": self._parse_size_to_mb(requested.get("jobfs")),
                # Used resources
                "mem_used_mb": self._parse_size_to_mb(used.get("mem")),
                "time_used_sec": self._parse_time_to_sec(used.get("walltime")),
                "cpus_used": self._parse_int(used.get("ncpus")),
                "jobfs_used_mb": self._parse_size_to_mb(used.get("jobfs")),
                "cpu_percent": self._parse_float(used.get("cpupercent")),
                # Efficiency metrics
                "mem_efficiency": efficiency.get("memory_efficiency"),
                "time_efficiency": efficiency.get("time_efficiency"),
                "cpu_efficiency": efficiency.get("cpu_efficiency"),
                "jobfs_efficiency": efficiency.get("jobfs_efficiency"),
                # Execution details
                "queue": execution.get("queue"),
                "hostname": execution.get("exec_host_parsed", {}).get("hostname"),
                # Timing
                "queue_wait_sec": timing.get("queue_wait_seconds"),
                "execution_sec": timing.get("execution_seconds"),
            }

            # Insert into database
            with self._connect() as conn:
                placeholders = ", ".join(["?" for _ in record])
                columns = ", ".join(record.keys())

                conn.execute(
                    f"INSERT OR REPLACE INTO job_resources ({columns}) VALUES ({placeholders})",
                    list(record.values()),
                )
                conn.commit()

            logger.debug("Logged resource data for job %s", job_id)
            return True

        except Exception as e:
            logger.debug("Failed to log resource data for job %s: %s", job_id, e)
            return False

    def update_job_resources(self, job_id: str, resource_data: Dict[str, Any]) -> bool:
        """
        Update resource usage for a job using parsed joblog data.

        Args:
            job_id: PBS job ID
            resource_data: Parsed resource data from parse_joblog_resources()

        Returns:
            bool: True if updated successfully
        """
        try:
            # Calculate efficiency metrics
            memory_efficiency = None
            time_efficiency = None
            jobfs_efficiency = None

            if resource_data.get("memory_requested_bytes") and resource_data.get(
                "memory_used_bytes"
            ):
                memory_efficiency = round(
                    (
                        resource_data["memory_used_bytes"]
                        / resource_data["memory_requested_bytes"]
                    )
                    * 100,
                    1,
                )

            if resource_data.get("walltime_requested_seconds") and resource_data.get(
                "walltime_used_seconds"
            ):
                time_efficiency = round(
                    (
                        resource_data["walltime_used_seconds"]
                        / resource_data["walltime_requested_seconds"]
                    )
                    * 100,
                    1,
                )

            if resource_data.get("jobfs_requested_bytes") and resource_data.get(
                "jobfs_used_bytes"
            ):
                jobfs_efficiency = round(
                    (
                        resource_data["jobfs_used_bytes"]
                        / resource_data["jobfs_requested_bytes"]
                    )
                    * 100,
                    1,
                )

            # Update the existing record with resource data
            update_fields = {
                "exit_code": resource_data.get("exit_status"),
                # Convert bytes to MB for storage
                "mem_requested_mb": (
                    round(
                        resource_data.get("memory_requested_bytes", 0) / (1024 * 1024),
                        2,
                    )
                    if resource_data.get("memory_requested_bytes")
                    else None
                ),
                "mem_used_mb": (
                    round(resource_data.get("memory_used_bytes", 0) / (1024 * 1024), 2)
                    if resource_data.get("memory_used_bytes")
                    else None
                ),
                "time_requested_sec": resource_data.get("walltime_requested_seconds"),
                "time_used_sec": resource_data.get("walltime_used_seconds"),
                "cpus_requested": resource_data.get("ncpus_requested"),
                "cpus_used": resource_data.get("ncpus_used"),
                "jobfs_requested_mb": (
                    round(
                        resource_data.get("jobfs_requested_bytes", 0) / (1024 * 1024), 2
                    )
                    if resource_data.get("jobfs_requested_bytes")
                    else None
                ),
                "jobfs_used_mb": (
                    round(resource_data.get("jobfs_used_bytes", 0) / (1024 * 1024), 2)
                    if resource_data.get("jobfs_used_bytes")
                    else None
                ),
                # Efficiency metrics
                "mem_efficiency": memory_efficiency,
                "time_efficiency": time_efficiency,
                "jobfs_efficiency": jobfs_efficiency,
            }

            # Build UPDATE SQL
            set_clauses = []
            values = []
            for field, value in update_fields.items():
                if value is not None:
                    set_clauses.append(f"{field} = ?")
                    values.append(value)

            if not set_clauses:
                logger.debug("No valid resource data to update for job %s", job_id)
                return False

            values.append(job_id)  # For WHERE clause

            with self._connect() as conn:
                # Legacy table
                sql = f"UPDATE job_resources SET {', '.join(set_clauses)} WHERE job_id = ?"
                result = conn.execute(sql, values)

                if result.rowcount == 0:
                    logger.debug(
                        "No existing record found to update for job %s", job_id
                    )

                # Unified queue table — same columns, different WHERE
                queue_values = values[:-1] + [job_id]
                queue_sql = (
                    f"UPDATE queue SET {', '.join(set_clauses)} WHERE pbs_job_id = ?"
                )
                conn.execute(queue_sql, queue_values)

                conn.commit()

            logger.debug("Updated resource data for job %s", job_id)
            return True

        except Exception as e:
            logger.debug("Failed to update resource data for job %s: %s", job_id, e)
            return False

    def _clean_command(self, command: Optional[str]) -> Optional[str]:
        """Clean up command string for better display."""
        if not command:
            return command

        # Remove the full path to qxub executable and replace with just 'qxub'
        import re

        # Pattern to match the full path to qxub executable
        # Examples: /g/data/ab12/software/venv/bin/qxub
        #           /some/other/path/to/qxub
        pattern = r"^[^\s]*qxub"
        cleaned = re.sub(pattern, "qxub", command)

        return cleaned

    def _parse_size_to_mb(self, size_str: Optional[str]) -> Optional[int]:
        """Parse size string to MB."""
        if not size_str:
            return None
        try:
            bytes_val = size_to_bytes(size_str)
            return int(bytes_val / (1024 * 1024))  # Convert to MB
        except Exception:
            return None

    def _parse_time_to_sec(self, time_str: Optional[str]) -> Optional[int]:
        """Parse time string to seconds."""
        if not time_str:
            return None
        try:
            return time_to_seconds(time_str)
        except Exception:
            return None

    def _parse_int(self, value) -> Optional[int]:
        """Parse integer value."""
        if value is None:
            return None
        try:
            return int(value)
        except Exception:
            return None

    def _parse_float(self, value) -> Optional[float]:
        """Parse float value."""
        if value is None:
            return None
        try:
            return float(value)
        except Exception:
            return None

    def get_recent_jobs(
        self,
        limit: int = 20,
        tags: Optional[Sequence[str]] = None,
        user: Optional[str] = None,
        since: Optional[str] = None,
        before: Optional[str] = None,
        queue_name: Optional[str] = None,
        status: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Get recent jobs with resource data, with optional filters.

        Args:
            limit:      Maximum rows to return.
            tags:       All listed tag values must be present on the job.
            user:       Exact username match.
            since:      ISO datetime string; only jobs submitted at or after this time.
            before:     ISO datetime string; only jobs submitted at or before this time.
            queue_name: Exact queue name match.
            status:     Exact status match (submitted/running/completed/failed/cancelled).
        """
        conditions: List[str] = []
        params: List[Any] = []

        if tags:
            for t in tags:
                conditions.append(
                    "EXISTS (SELECT 1 FROM json_each(tags) WHERE value = ?)"
                )
                params.append(t)
        if user:
            conditions.append("username = ?")
            params.append(user)
        if since:
            conditions.append("COALESCE(submitted_at, timestamp) >= ?")
            params.append(since)
        if before:
            conditions.append("COALESCE(submitted_at, timestamp) <= ?")
            params.append(before)
        if queue_name:
            conditions.append("queue = ?")
            params.append(queue_name)
        if status:
            conditions.append("status = ?")
            params.append(status)

        where_clause = ("WHERE " + " AND ".join(conditions)) if conditions else ""
        query = f"""
            SELECT * FROM job_resources
            {where_clause}
            ORDER BY COALESCE(submitted_at, timestamp) DESC
            LIMIT ?
        """
        params.append(limit)

        with self._connect() as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]

    def get_efficiency_stats(self) -> Dict[str, Any]:
        """Get overall efficiency statistics."""
        with self._connect() as conn:
            cursor = conn.execute(
                """
                SELECT
                    COUNT(*) as total_jobs,
                    AVG(mem_efficiency) as avg_mem_efficiency,
                    AVG(time_efficiency) as avg_time_efficiency,
                    AVG(cpu_efficiency) as avg_cpu_efficiency,
                    AVG(jobfs_efficiency) as avg_jobfs_efficiency,
                    COUNT(CASE WHEN mem_efficiency < 50 THEN 1 END) as low_mem_efficiency_jobs,
                    COUNT(CASE WHEN time_efficiency < 50 THEN 1 END) as low_time_efficiency_jobs,
                    COUNT(CASE WHEN cpu_efficiency < 50 THEN 1 END) as low_cpu_efficiency_jobs
                FROM job_resources
                WHERE mem_efficiency IS NOT NULL
            """
            )

            row = cursor.fetchone()
            if row:
                return {
                    "total_jobs": row[0],
                    "avg_mem_efficiency": round(row[1], 1) if row[1] else 0,
                    "avg_time_efficiency": round(row[2], 1) if row[2] else 0,
                    "avg_cpu_efficiency": round(row[3], 1) if row[3] else 0,
                    "avg_jobfs_efficiency": round(row[4], 1) if row[4] else 0,
                    "low_mem_efficiency_jobs": row[5],
                    "low_time_efficiency_jobs": row[6],
                    "low_cpu_efficiency_jobs": row[7],
                }
            return {}

    def get_inefficient_jobs(
        self, efficiency_threshold: float = 50.0, limit: int = 10
    ) -> List[Dict[str, Any]]:
        """Get jobs with low resource efficiency."""
        with self._connect() as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                """
                SELECT job_id, command, mem_efficiency, time_efficiency,
                       cpu_efficiency, mem_requested_mb, mem_used_mb,
                       time_requested_sec, time_used_sec, timestamp
                FROM job_resources
                WHERE (mem_efficiency < ? OR time_efficiency < ? OR cpu_efficiency < ?)
                  AND mem_efficiency IS NOT NULL
                ORDER BY timestamp DESC
                LIMIT ?
            """,
                (
                    efficiency_threshold,
                    efficiency_threshold,
                    efficiency_threshold,
                    limit,
                ),
            )

            return [dict(row) for row in cursor.fetchall()]

    def get_resource_trends(self, days: int = 30) -> List[Dict[str, Any]]:
        """Get resource usage trends over time."""
        with self._connect() as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                """
                SELECT
                    DATE(timestamp) as date,
                    COUNT(*) as job_count,
                    AVG(mem_efficiency) as avg_mem_eff,
                    AVG(time_efficiency) as avg_time_eff,
                    AVG(cpu_efficiency) as avg_cpu_eff,
                    SUM(mem_requested_mb) as total_mem_req_mb,
                    SUM(mem_used_mb) as total_mem_used_mb
                FROM job_resources
                WHERE timestamp > date('now', '-{} days')
                  AND mem_efficiency IS NOT NULL
                GROUP BY DATE(timestamp)
                ORDER BY date DESC
            """.format(
                    days
                )
            )

            return [dict(row) for row in cursor.fetchall()]

    # Job Status Tracking Methods (new for v3)

    def log_job_submitted(
        self,
        job_id: str,
        command: str,
        tags: Optional[Sequence[str]] = None,
        username: Optional[str] = None,
        joblog_path: Optional[str] = None,
        queue: Optional[str] = None,
        cpus_requested: Optional[int] = None,
        job_name: Optional[str] = None,
    ) -> bool:
        """Log initial job submission.

        Writes to the legacy ``job_resources`` table **and** updates the
        unified ``queue`` row (looked up by ``pbs_job_id``) so that both
        tables stay in sync during the transition period.
        """
        try:
            now = datetime.now().isoformat()
            clean_command = self._clean_command(command)
            tags_json = json.dumps(list(tags) if tags else [])
            if username is None:
                try:
                    import getpass

                    username = getpass.getuser()
                except Exception:
                    username = ""

            with self._connect() as conn:
                # Legacy table
                conn.execute(
                    """
                    INSERT OR REPLACE INTO job_resources
                    (job_id, timestamp, command, status, submitted_at, last_status_update,
                     tags, username, joblog_path, queue, cpus_requested, job_name)
                    VALUES (?, ?, ?, 'submitted', ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        job_id,
                        now,
                        clean_command,
                        now,
                        now,
                        tags_json,
                        username,
                        joblog_path,
                        queue,
                        cpus_requested,
                        job_name,
                    ),
                )

                # Unified queue table — update the row that was created
                # by create_queue_entry() before qsub.
                conn.execute(
                    """
                    UPDATE queue
                    SET status='submitted', last_status_update=?,
                        username=?, joblog_path=COALESCE(joblog_path, ?),
                        queue_name=COALESCE(queue_name, ?),
                        cpus_requested=COALESCE(cpus_requested, ?)
                    WHERE pbs_job_id=?
                    """,
                    (now, username, joblog_path, queue, cpus_requested, job_id),
                )

                conn.commit()

            logger.debug(
                "Logged job submission for %s (user: %s, tags: %s)",
                job_id,
                username,
                tags_json,
            )
            return True
        except Exception as e:
            logger.debug("Failed to log job submission for %s: %s", job_id, e)
            return False

    def update_job_status(self, job_id: str, status: str) -> bool:
        """Update job status (submitted -> running -> completed/failed/cancelled).

        Writes to both ``job_resources`` (legacy) and unified ``queue`` table.
        """
        try:
            now = datetime.now().isoformat()

            with self._connect() as conn:
                # Legacy table
                conn.execute(
                    "UPDATE job_resources SET status=?, last_status_update=? WHERE job_id=?",
                    (status, now, job_id),
                )

                # Update specific timestamp columns based on status
                if status == "running":
                    conn.execute(
                        "UPDATE job_resources SET started_at=? WHERE job_id=?",
                        (now, job_id),
                    )
                elif status in ("completed", "failed"):
                    conn.execute(
                        "UPDATE job_resources SET completed_at=? WHERE job_id=?",
                        (now, job_id),
                    )

                # Unified queue table
                conn.execute(
                    "UPDATE queue SET status=?, last_status_update=? WHERE pbs_job_id=?",
                    (status, now, job_id),
                )
                if status == "running":
                    conn.execute(
                        "UPDATE queue SET started_at=? WHERE pbs_job_id=?",
                        (now, job_id),
                    )
                elif status in ("completed", "failed"):
                    conn.execute(
                        "UPDATE queue SET completed_at=? WHERE pbs_job_id=?",
                        (now, job_id),
                    )

                conn.commit()

            logger.debug("Updated job %s status to %s", job_id, status)
            return True
        except Exception as e:
            logger.debug("Failed to update job status for %s: %s", job_id, e)
            return False

    def job_ids_starting_with(self, prefix: str) -> List[str]:
        """Job IDs in the DB that start with ``prefix`` (either table).

        Lets callers match ``12345`` to a stored ``12345.gadi-pbs``. Returns
        an empty list if the DB can't be read.
        """
        if not prefix:
            return []
        pattern = (
            prefix.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "%"
        )
        try:
            with self._connect(timeout=_STATUS_READ_TIMEOUT) as conn:
                rows = conn.execute(
                    "SELECT pbs_job_id FROM queue WHERE pbs_job_id LIKE ? ESCAPE '\\' "
                    "UNION SELECT job_id FROM job_resources "
                    "WHERE job_id LIKE ? ESCAPE '\\'",
                    (pattern, pattern),
                ).fetchall()
            return [r[0] for r in rows if r[0]]
        except Exception as e:  # pylint: disable=broad-except
            logger.debug("Failed to look up job IDs like %s: %s", prefix, e)
            return []

    def update_job_exit_code(self, job_id: str, exit_code: int) -> bool:
        """Update the exit code for a job."""
        try:
            with self._connect() as conn:
                conn.execute(
                    "UPDATE job_resources SET exit_code=? WHERE job_id=?",
                    (exit_code, job_id),
                )
                conn.commit()
            logger.debug("Updated job %s exit_code to %s", job_id, exit_code)
            return True
        except Exception as e:
            logger.debug("Failed to update exit code for %s: %s", job_id, e)
            return False

    def finalize_job(self, job_id: str, exit_code: int, joblog_path: str) -> bool:
        """Record job completion: update status and exit code.

        Job scripts no longer call this (see ``qxub.queue.runs``).
        Resource data (mem, walltime, cpu) is NOT collected here because PBS
        only appends the resource section to the joblog AFTER the script exits.
        The joblog_path is stored so that ``backfill_resources()`` can pick it
        up the next time ``qxub resources list`` is run.

        Writes to both ``job_resources`` (legacy) and unified ``queue`` table.
        """
        now = datetime.now().isoformat()
        status = "completed" if exit_code == 0 else "failed"
        try:
            with self._connect() as conn:
                # Legacy table
                conn.execute(
                    """UPDATE job_resources
                       SET status=?, exit_code=?, completed_at=?, last_status_update=?,
                           joblog_path=COALESCE(joblog_path, ?)
                       WHERE job_id=?""",
                    (status, exit_code, now, now, joblog_path, job_id),
                )

                # Unified queue table
                conn.execute(
                    """UPDATE queue
                       SET status=?, exit_code=?, completed_at=?, last_status_update=?,
                           joblog_path=COALESCE(joblog_path, ?)
                       WHERE pbs_job_id=?""",
                    (status, exit_code, now, now, joblog_path, job_id),
                )

                conn.commit()
            logger.debug(
                "Finalized job %s (exit=%s, status=%s)", job_id, exit_code, status
            )
            return True
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("Failed to finalize job %s: %s", job_id, exc)
            return False

    # ------------------------------------------------------------------
    # Joblog discovery helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_job_id_from_log(path: str) -> Optional[str]:
        """Quickly extract the PBS Job Id from the tail of a .log file.

        PBS appends the resource section at the very end of the file, so
        reading the last 2 KB is enough to find "Job Id: <id>".
        """
        import re  # pylint: disable=import-outside-toplevel

        try:
            with open(path, "rb") as fh:
                fh.seek(0, 2)
                size = fh.tell()
                fh.seek(max(0, size - 2048))
                tail = fh.read().decode("utf-8", errors="replace")
            m = re.search(r"Job Id:\s*(\S+)", tail)
            return m.group(1) if m else None
        except Exception:  # pylint: disable=broad-except
            return None

    def _discover_joblog_paths(self, log_dirs: List[str]) -> Dict[str, str]:
        """Scan *log_dirs* for ``*.log`` files and return ``{job_id: path}``."""
        import glob  # pylint: disable=import-outside-toplevel

        index: Dict[str, str] = {}
        for d in log_dirs:
            for path in glob.glob(os.path.join(d, "*.log")):
                job_id = self._extract_job_id_from_log(path)
                if job_id:
                    index[job_id] = path
        return index

    def _resolve_log_dirs(self) -> List[str]:
        """Return candidate log directories derived from already-stored paths."""
        try:
            with self._connect() as conn:
                rows = conn.execute(
                    "SELECT DISTINCT joblog_path FROM job_resources"
                    " WHERE joblog_path IS NOT NULL LIMIT 20"
                ).fetchall()
            dirs: List[str] = []
            seen: set = set()
            for (p,) in rows:
                d = os.path.dirname(p)
                if d and d not in seen and os.path.isdir(d):
                    dirs.append(d)
                    seen.add(d)
            return dirs
        except Exception:  # pylint: disable=broad-except
            return []

    # ------------------------------------------------------------------
    # Public backfill entry point
    # ------------------------------------------------------------------

    def backfill_resources(self, limit: int = 200) -> int:
        """Parse joblogs for recently completed jobs that have no resource data.

        Two-phase approach:
        1. Discovery — for completed jobs whose ``joblog_path`` is still NULL,
           scan known log directories to find the matching file.
        2. Parsing  — for all backfillable rows (now including newly discovered
           ones), parse the joblog and populate the resource columns.

        Designed to be called automatically on ``qxub resources list`` so data
        appears without any manual intervention.

        Returns the number of jobs successfully backfilled.
        """
        from .parser import (  # pylint: disable=import-outside-toplevel
            parse_joblog_resources,
        )

        # ---- Phase 1: discover joblog paths for jobs that lack them ----------
        try:
            with self._connect() as conn:
                missing = conn.execute(
                    """
                    SELECT job_id FROM job_resources
                    WHERE status IN ('completed', 'failed')
                      AND mem_used_mb IS NULL
                      AND joblog_path IS NULL
                      AND submitted_at > datetime('now', '-180 days')
                    ORDER BY completed_at DESC
                    LIMIT ?
                    """,
                    (limit,),
                ).fetchall()
            missing_ids = {row[0] for row in missing}
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("backfill phase-1 query failed: %s", exc)
            missing_ids = set()

        if missing_ids:
            log_dirs = self._resolve_log_dirs()
            if log_dirs:
                index = self._discover_joblog_paths(log_dirs)
                matched = {
                    jid: path for jid, path in index.items() if jid in missing_ids
                }
                if matched:
                    try:
                        with self._connect() as conn:
                            conn.executemany(
                                "UPDATE job_resources SET joblog_path=?"
                                " WHERE job_id=? AND joblog_path IS NULL",
                                [(path, jid) for jid, path in matched.items()],
                            )
                            conn.commit()
                        logger.debug(
                            "Discovered joblog paths for %d historical jobs",
                            len(matched),
                        )
                    except Exception as exc:  # pylint: disable=broad-except
                        logger.debug("backfill path-update failed: %s", exc)

        # ---- Phase 2: parse joblogs and fill resource columns ---------------
        try:
            with self._connect() as conn:
                conn.row_factory = sqlite3.Row
                rows = conn.execute(
                    """
                    SELECT job_id, joblog_path FROM job_resources
                    WHERE status IN ('completed', 'failed')
                      AND mem_used_mb IS NULL
                      AND joblog_path IS NOT NULL
                    ORDER BY completed_at DESC
                    LIMIT ?
                    """,
                    (limit,),
                ).fetchall()
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("backfill phase-2 query failed: %s", exc)
            return 0

        filled = 0
        for row in rows:
            try:
                data = parse_joblog_resources(row["joblog_path"])
                if data:
                    ok = self.update_job_resources(row["job_id"], data)
                    if ok:
                        filled += 1
                        logger.debug("Backfilled resources for %s", row["job_id"])
            except Exception as exc:  # pylint: disable=broad-except
                logger.debug("Backfill failed for %s: %s", row["job_id"], exc)
        return filled

    def sync_from_run_files(
        self, pbs_job_ids: Optional[Sequence[str]] = None, limit: int = 500
    ) -> int:
        """Apply run-file state to unfinished jobs in the database.

        Jobs no longer write the database themselves; they update their run
        file (``qxub.queue.runs``). Readers call this first so ``list``,
        ``summary``, ``show`` and ``cancel`` see current states. Only rows
        not yet finished are looked at (the ``limit`` most recent, or just
        ``pbs_job_ids``), and only changed rows are written, in one short
        transaction from this (login) node.

        Returns the number of jobs updated. Never raises.
        """
        from ..queue import runs  # pylint: disable=import-outside-toplevel

        try:
            with self._connect() as conn:
                query = (
                    "SELECT entry_id, pbs_job_id, status FROM queue"
                    " WHERE status IN ('initiated', 'dispatched', 'submitted', 'running')"
                )
                params: List[Any] = []
                if pbs_job_ids is not None:
                    if not pbs_job_ids:
                        return 0
                    query += " AND pbs_job_id IN (%s)" % ",".join(
                        "?" * len(pbs_job_ids)
                    )
                    params.extend(pbs_job_ids)
                query += " ORDER BY submitted_at DESC LIMIT ?"
                params.append(limit)
                rows = conn.execute(query, params).fetchall()
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("sync_from_run_files: query failed: %s", exc)
            return 0

        updates = []
        for entry_id, pbs_job_id, status in rows:
            record = runs.read_run(runs.run_file(entry_id))
            if not record:
                continue
            state = record.get("state")
            job_id = pbs_job_id or record.get("native_job_id") or None
            if state == "finished":
                code = runs.exit_code(record)
                new_status = "completed" if code == 0 else "failed"
            elif state == "running" and status != "running":
                new_status, code = "running", None
            else:
                continue
            updates.append(
                (
                    entry_id,
                    job_id,
                    new_status,
                    code,
                    record.get("started_at") or None,
                    record.get("finished_at") or None,
                )
            )

        if not updates:
            return 0
        now = datetime.now().isoformat()
        try:
            with self._connect() as conn:
                for entry_id, job_id, status, code, started, finished in updates:
                    conn.execute(
                        """UPDATE queue
                           SET status=?, exit_code=COALESCE(?, exit_code),
                               started_at=COALESCE(started_at, ?),
                               completed_at=COALESCE(?, completed_at),
                               pbs_job_id=COALESCE(pbs_job_id, ?),
                               last_status_update=?
                           WHERE entry_id=?""",
                        (status, code, started, finished, job_id, now, entry_id),
                    )
                    if not job_id:
                        continue
                    conn.execute(
                        """UPDATE job_resources
                           SET status=?, exit_code=COALESCE(?, exit_code),
                               started_at=COALESCE(started_at, ?),
                               completed_at=COALESCE(?, completed_at),
                               last_status_update=?
                           WHERE job_id=?""",
                        (status, code, started, finished, now, job_id),
                    )
                    conn.execute(
                        """UPDATE project_jobs
                           SET status=?, exit_code=COALESCE(?, exit_code),
                               started_at=COALESCE(started_at, ?),
                               completed_at=COALESCE(?, completed_at)
                           WHERE pbs_job_id=?""",
                        (status, code, started, finished, job_id),
                    )
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("sync_from_run_files: update failed: %s", exc)
            return 0
        return len(updates)

    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get current status of a specific job.

        Checks the unified ``queue`` table first (by ``pbs_job_id``), then
        falls back to the legacy ``job_resources`` table.

        Returns ``None`` when the job is genuinely not in the database.
        Raises :class:`DatabaseError` on transient sqlite failures so the
        caller can distinguish "not found" from "DB unreachable".

        Retries up to 3 times with exponential backoff to tolerate
        filesystem latency spikes and lock contention.
        """
        last_exc: Optional[Exception] = None
        for attempt in range(3):
            try:
                with self._connect(timeout=_STATUS_READ_TIMEOUT) as conn:
                    conn.row_factory = sqlite3.Row

                    # Try unified queue table first (new jobs)
                    row = conn.execute(
                        """
                        SELECT pbs_job_id AS job_id, status, command,
                               submitted_at, started_at, completed_at,
                               last_status_update, exit_code, joblog_path
                        FROM queue WHERE pbs_job_id=?
                        """,
                        (job_id,),
                    ).fetchone()
                    if row:
                        return dict(row)

                    # Fall back to legacy job_resources table
                    row = conn.execute(
                        """
                        SELECT job_id, status, command, submitted_at, started_at,
                               completed_at, last_status_update, exit_code, joblog_path
                        FROM job_resources WHERE job_id=?
                        """,
                        (job_id,),
                    ).fetchone()
                    return dict(row) if row else None
            except Exception as e:  # pylint: disable=broad-except
                last_exc = e
                logger.debug(
                    "get_job_status attempt %d failed for %s: %s",
                    attempt + 1,
                    job_id,
                    e,
                )
                if attempt < 2:
                    time.sleep(0.5 * (2**attempt))  # 0.5 s, 1 s

        # All retries exhausted — raise so the caller knows the DB failed.
        raise DatabaseError(
            f"Failed to query job status for {job_id} after 3 attempts: {last_exc}"
        ) from last_exc

    def get_jobs_by_status(
        self, status: str = None, limit: int = None
    ) -> List[Dict[str, Any]]:
        """Get jobs filtered by status."""
        try:
            with self._connect() as conn:
                conn.row_factory = sqlite3.Row

                if status:
                    query = """
                        SELECT job_id, status, command, submitted_at, started_at,
                               completed_at, last_status_update, exit_code,
                               job_name, tags
                        FROM job_resources
                        WHERE status=?
                        ORDER BY submitted_at DESC
                    """
                    params = (status,)
                else:
                    query = """
                        SELECT job_id, status, command, submitted_at, started_at,
                               completed_at, last_status_update, exit_code,
                               job_name, tags
                        FROM job_resources
                        ORDER BY submitted_at DESC
                    """
                    params = ()

                if limit:
                    query += f" LIMIT {limit}"

                cursor = conn.execute(query, params)
                return [dict(row) for row in cursor.fetchall()]
        except Exception as e:
            logger.debug("Failed to get jobs by status: %s", e)
            return []

    def settle_gone_jobs(self, job_ids: Sequence[str]) -> Dict[str, int]:
        """Finish unfinished rows for jobs the scheduler no longer lists.

        These are rows left ``submitted`` or ``running`` with no run file to
        say how they ended: jobs from before run files, jobs whose last update
        was lost, or jobs whose run file consolidation has already removed
        (``consolidate.GRACE_DAYS``). If the statistics snapshot knows how the
        job ended, the row takes its status, exit code and end time.
        Otherwise it becomes ``unknown``. Only unfinished rows change.

        Returns counts: ``{"from_snapshot": n, "unknown": m}``.
        """
        counts = {"from_snapshot": 0, "unknown": 0}
        if not job_ids:
            return counts
        ended = _snapshot_outcomes(job_ids)
        now = datetime.now().isoformat()
        unfinished = "('initiated', 'dispatched', 'submitted', 'running')"
        try:
            with self._connect() as conn:
                for job_id in job_ids:
                    status, code, finished = ended.get(job_id, ("unknown", None, None))
                    params = (status, code, finished, now, job_id)
                    cur = conn.execute(
                        "UPDATE job_resources SET status=?,"
                        " exit_code=COALESCE(?, exit_code),"
                        " completed_at=COALESCE(?, completed_at),"
                        " last_status_update=?"
                        f" WHERE job_id=? AND status IN {unfinished}",
                        params,
                    )
                    conn.execute(
                        "UPDATE queue SET status=?,"
                        " exit_code=COALESCE(?, exit_code),"
                        " completed_at=COALESCE(?, completed_at),"
                        " last_status_update=?"
                        f" WHERE pbs_job_id=? AND status IN {unfinished}",
                        params,
                    )
                    if cur.rowcount:
                        key = "unknown" if status == "unknown" else "from_snapshot"
                        counts[key] += 1
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("settle_gone_jobs failed: %s", exc)
            return {"from_snapshot": 0, "unknown": 0}
        return counts

    def get_status_summary(self) -> Dict[str, int]:
        """Get counts of jobs by status."""
        try:
            with self._connect() as conn:
                cursor = conn.execute(
                    """
                    SELECT status, COUNT(*) as count
                    FROM job_resources
                    WHERE submitted_at > datetime('now', '-7 days')
                    GROUP BY status
                    """
                )
                return {row[0]: row[1] for row in cursor.fetchall()}
        except Exception as e:
            logger.debug("Failed to get status summary: %s", e)
            return {}

    def cleanup_old_jobs(self, days_old: int = 30) -> int:
        """Remove job records older than specified days. Returns count of deleted jobs."""
        try:
            with self._connect() as conn:
                cursor = conn.execute(
                    """
                    DELETE FROM job_resources
                    WHERE completed_at < datetime('now', '-{} days')
                    OR (status = 'submitted' AND submitted_at < datetime('now', '-{} days'))
                    """.format(
                        days_old, days_old // 2
                    )  # Clean up stuck submitted jobs sooner
                )
                deleted_count = cursor.rowcount
                conn.commit()

            logger.debug("Cleaned up %d old job records", deleted_count)
            return deleted_count
        except Exception as e:
            logger.debug("Failed to cleanup old jobs: %s", e)
            return 0

    def export_csv(self, output_path: Path, limit: Optional[int] = None):
        """Export resource data to CSV."""
        import csv

        with self._connect() as conn:
            query = "SELECT * FROM job_resources ORDER BY timestamp DESC"
            if limit:
                query += f" LIMIT {limit}"

            cursor = conn.execute(query)

            with open(output_path, "w", newline="") as csvfile:
                # Get column names
                columns = [description[0] for description in cursor.description]
                writer = csv.DictWriter(csvfile, fieldnames=columns)

                writer.writeheader()
                for row in cursor:
                    writer.writerow(dict(zip(columns, row)))


class SnapshotReader(ResourceTracker):
    """The report queries of ``ResourceTracker``, over the published snapshot.

    Opens ``qxub-stats.db`` read-only with ``immutable=1`` (see
    ``qxub.queue.consolidate.open_snapshot``), so reports never write or lock
    anything. Only the read methods make sense here.
    """

    def __init__(self, db_path: Path):  # pylint: disable=super-init-not-called
        self.db_path = Path(db_path)

    def _connect(self, timeout: float = _SQLITE_TIMEOUT):
        from ..queue.consolidate import (  # pylint: disable=import-outside-toplevel
            open_snapshot,
        )

        return _Closing(open_snapshot(self.db_path))


class _Closing:
    """``with`` that closes a read-only connection on the way out."""

    def __init__(self, conn):
        self.conn = conn

    def __enter__(self):
        return self.conn

    def __exit__(self, *exc):
        self.conn.close()
        return False


# Lazy global instance — avoids running _init_database() at import time,
# which would fail fatally if the DB is corrupted or on a broken filesystem.
_resource_tracker = None


def _get_resource_tracker():
    global _resource_tracker  # noqa: PLW0603
    if _resource_tracker is None:
        _resource_tracker = ResourceTracker()
    return _resource_tracker


class _LazyResourceTracker:
    """Proxy that defers ResourceTracker construction until first use."""

    def __getattr__(self, name):
        return getattr(_get_resource_tracker(), name)


resource_tracker = _LazyResourceTracker()
