"""
qxub.queue - Queue and Database Management Package

Central module for:
- Unified database path resolution (shared or per-user)
- Queue table for virtual job IDs and pending dispatch (Phase 3+)
- project_jobs table for PBS-level job tracking
- mark_running / mark_complete helpers used by job scripts

Public API:
    get_db_path()                  -> Path  # Resolve active DB path
    init_db(db_path=None)                  # Create/migrate all tables
    new_virtual_id()               -> str  # Generate qx-{uuid}
    is_virtual_id(job_id)          -> bool
    create_queue_entry(...)        -> str  # Register a dispatched job; returns virtual ID
    mark_running(pbs_job_id, ...)          # Update status on job start
    mark_complete(pbs_job_id, exit_code, ...) # Update status on job end
    resolve_virtual_id(virtual_id) -> dict | None
"""

from .db import (
    create_queue_entry,
    get_db_path,
    init_db,
    is_virtual_id,
    mark_complete,
    mark_running,
    new_virtual_id,
    resolve_virtual_id,
)

__all__ = [
    "get_db_path",
    "init_db",
    "new_virtual_id",
    "is_virtual_id",
    "create_queue_entry",
    "mark_running",
    "mark_complete",
    "resolve_virtual_id",
]
