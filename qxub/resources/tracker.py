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


# Default timeout (seconds) for all sqlite connections.  The shared DB on
# /g/data can experience lock contention when many concurrent qxtat check
# processes race on the same file.  30 s is generous enough to ride out
# filesystem latency spikes on Gadi /scratch.
_SQLITE_TIMEOUT = 30


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
        """Return a new SQLite connection with mmap disabled.

        Disabling mmap (``PRAGMA mmap_size=0``) prevents SIGBUS crashes
        on shared/network filesystems (Lustre, GPFS) where concurrent
        writers can corrupt memory-mapped pages.
        """
        conn = sqlite3.connect(self.db_path, timeout=timeout)
        conn.execute("PRAGMA mmap_size=0")
        return conn

    def _init_database(self):
        """Initialize SQLite database with resource tracking table."""
        with self._connect() as conn:
            # DELETE mode avoids -shm sidecar which is always mmap'd.
            # WAL's -shm causes SIGBUS on shared filesystems (Lustre/GPFS).
            conn.execute("PRAGMA journal_mode=DELETE")
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
                    joblog_path TEXT
                )
            """
            )
            conn.commit()

            # Migrate existing database if needed
            self._migrate_database_schema(conn)

    def _migrate_database_schema(self, conn):
        """Add new status tracking columns to existing database."""
        try:
            # Check if status column exists
            cursor = conn.execute("PRAGMA table_info(job_resources)")
            columns = [row[1] for row in cursor.fetchall()]

            if "status" not in columns:
                logger.debug("Migrating database schema to add status tracking columns")

                # Add new columns
                conn.execute(
                    "ALTER TABLE job_resources ADD COLUMN status TEXT DEFAULT 'completed'"
                )
                conn.execute("ALTER TABLE job_resources ADD COLUMN submitted_at TEXT")
                conn.execute("ALTER TABLE job_resources ADD COLUMN started_at TEXT")
                conn.execute("ALTER TABLE job_resources ADD COLUMN completed_at TEXT")
                conn.execute(
                    "ALTER TABLE job_resources ADD COLUMN last_status_update TEXT"
                )

                # Update existing records to have completed status and timestamps
                conn.execute(
                    """
                    UPDATE job_resources
                    SET status='completed',
                        completed_at=timestamp,
                        last_status_update=timestamp
                    WHERE status IS NULL OR status = 'completed'
                """
                )

                conn.commit()
                logger.debug("Database schema migration completed")

            if "tags" not in columns:
                logger.debug("Migrating database schema to add tags column")
                conn.execute(
                    "ALTER TABLE job_resources ADD COLUMN tags TEXT DEFAULT '[]'"
                )
                conn.commit()
                logger.debug("Tags column migration completed")

            if "username" not in columns:
                logger.debug("Migrating database schema to add username column")
                conn.execute("ALTER TABLE job_resources ADD COLUMN username TEXT")
                conn.commit()
                logger.debug("Username column migration completed")

            if "joblog_path" not in columns:
                logger.debug("Migrating database schema to add joblog_path column")
                conn.execute("ALTER TABLE job_resources ADD COLUMN joblog_path TEXT")
                conn.commit()
                logger.debug("joblog_path column migration completed")
        except Exception as e:
            logger.debug("Database migration failed (may be normal): %s", e)

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
        # Examples: /g/data/a56/software/qsub_tools/venv/bin/qxub
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

    def get_resource_trends(self, days: int = 30) -> Dict[str, Any]:
        """Get resource usage trends over time."""
        with self._connect() as conn:
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
                     tags, username, joblog_path, queue, cpus_requested)
                    VALUES (?, ?, ?, 'submitted', ?, ?, ?, ?, ?, ?, ?)
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

        Called from within the PBS job script at end-of-job.
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
                with self._connect() as conn:
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
                               completed_at, last_status_update, exit_code
                        FROM job_resources
                        WHERE status=?
                        ORDER BY submitted_at DESC
                    """
                    params = (status,)
                else:
                    query = """
                        SELECT job_id, status, command, submitted_at, started_at,
                               completed_at, last_status_update, exit_code
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
