"""
Simple resource tracking system for qxub.

Tracks only key metrics: job_id, requested vs used resources, efficiency.
Uses SQLite for simple querying and analysis.
"""

import logging
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from .resource_parser import bytes_to_human, size_to_bytes, time_to_seconds


class ResourceTracker:
    """Simple resource tracking focused on efficiency metrics."""

    def __init__(self, db_path: Optional[Path] = None):
        """Initialize resource tracker with SQLite database."""
        if db_path is None:
            # Use XDG config directory
            config_dir = Path.home() / ".config" / "qxub"
            config_dir.mkdir(parents=True, exist_ok=True)
            db_path = config_dir / "resources.db"

        self.db_path = Path(db_path)
        self._init_database()

    def _init_database(self):
        """Initialize SQLite database with resource tracking table."""
        with sqlite3.connect(self.db_path) as conn:
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
                    execution_sec INTEGER
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
                logging.debug(
                    "Migrating database schema to add status tracking columns"
                )

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
                logging.debug("Database schema migration completed")
        except Exception as e:
            logging.debug("Database migration failed (may be normal): %s", e)

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
            with sqlite3.connect(self.db_path) as conn:
                placeholders = ", ".join(["?" for _ in record])
                columns = ", ".join(record.keys())

                conn.execute(
                    f"INSERT OR REPLACE INTO job_resources ({columns}) VALUES ({placeholders})",
                    list(record.values()),
                )
                conn.commit()

            logging.debug("Logged resource data for job %s", job_id)
            return True

        except Exception as e:
            logging.debug("Failed to log resource data for job %s: %s", job_id, e)
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
                logging.debug("No valid resource data to update for job %s", job_id)
                return False

            values.append(job_id)  # For WHERE clause

            with sqlite3.connect(self.db_path) as conn:
                sql = f"UPDATE job_resources SET {', '.join(set_clauses)} WHERE job_id = ?"
                result = conn.execute(sql, values)

                if result.rowcount == 0:
                    logging.debug(
                        "No existing record found to update for job %s", job_id
                    )
                    return False

                conn.commit()

            logging.debug("Updated resource data for job %s", job_id)
            return True

        except Exception as e:
            logging.debug("Failed to update resource data for job %s: %s", job_id, e)
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

    def get_recent_jobs(self, limit: int = 20) -> List[Dict[str, Any]]:
        """Get recent jobs with resource data."""
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                """
                SELECT * FROM job_resources
                ORDER BY timestamp DESC
                LIMIT ?
            """,
                (limit,),
            )

            return [dict(row) for row in cursor.fetchall()]

    def get_efficiency_stats(self) -> Dict[str, Any]:
        """Get overall efficiency statistics."""
        with sqlite3.connect(self.db_path) as conn:
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
        with sqlite3.connect(self.db_path) as conn:
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
        with sqlite3.connect(self.db_path) as conn:
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

    def log_job_submitted(self, job_id: str, command: str) -> bool:
        """Log initial job submission."""
        try:
            now = datetime.now().isoformat()
            clean_command = self._clean_command(command)

            with sqlite3.connect(self.db_path) as conn:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO job_resources
                    (job_id, timestamp, command, status, submitted_at, last_status_update)
                    VALUES (?, ?, ?, 'submitted', ?, ?)
                    """,
                    (job_id, now, clean_command, now, now),
                )
                conn.commit()

            logging.debug("Logged job submission for %s", job_id)
            return True
        except Exception as e:
            logging.debug("Failed to log job submission for %s: %s", job_id, e)
            return False

    def update_job_status(self, job_id: str, status: str) -> bool:
        """Update job status (submitted -> running -> completed/failed)."""
        try:
            now = datetime.now().isoformat()

            with sqlite3.connect(self.db_path) as conn:
                # Update status and timestamp
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

                conn.commit()

            logging.debug("Updated job %s status to %s", job_id, status)
            return True
        except Exception as e:
            logging.debug("Failed to update job status for %s: %s", job_id, e)
            return False

    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get current status of a specific job."""
        try:
            with sqlite3.connect(self.db_path) as conn:
                conn.row_factory = sqlite3.Row
                cursor = conn.execute(
                    """
                    SELECT job_id, status, command, submitted_at, started_at,
                           completed_at, last_status_update, exit_code
                    FROM job_resources WHERE job_id=?
                    """,
                    (job_id,),
                )
                row = cursor.fetchone()
                return dict(row) if row else None
        except Exception as e:
            logging.debug("Failed to get job status for %s: %s", job_id, e)
            return None

    def get_jobs_by_status(
        self, status: str = None, limit: int = None
    ) -> List[Dict[str, Any]]:
        """Get jobs filtered by status."""
        try:
            with sqlite3.connect(self.db_path) as conn:
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
            logging.debug("Failed to get jobs by status: %s", e)
            return []

    def get_status_summary(self) -> Dict[str, int]:
        """Get counts of jobs by status."""
        try:
            with sqlite3.connect(self.db_path) as conn:
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
            logging.debug("Failed to get status summary: %s", e)
            return {}

    def cleanup_old_jobs(self, days_old: int = 30) -> int:
        """Remove job records older than specified days. Returns count of deleted jobs."""
        try:
            with sqlite3.connect(self.db_path) as conn:
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

            logging.debug("Cleaned up %d old job records", deleted_count)
            return deleted_count
        except Exception as e:
            logging.debug("Failed to cleanup old jobs: %s", e)
            return 0

    def export_csv(self, output_path: Path, limit: Optional[int] = None):
        """Export resource data to CSV."""
        import csv

        with sqlite3.connect(self.db_path) as conn:
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


# Global instance
resource_tracker = ResourceTracker()
