"""
qxub.queue - Queue and Database Management Package

Central module for:
- Unified database path resolution (shared or per-user)
- Unified job table: virtual ID as primary key, PBS ID nullable
- project_jobs table for PBS-level job tracking
- mark_running / mark_complete helpers used by job scripts

Public API:
    get_db_path()                  -> Path  # Resolve active DB path
    init_db(db_path=None)                  # Create/migrate all tables
    new_virtual_id()               -> str  # Generate qx-{uuid}
    is_virtual_id(job_id)          -> bool
    create_queue_entry(...)        -> str  # Register a job; returns virtual ID
    update_queue_entry(...)        -> bool # Fill in PBS ID after qsub
    get_entries_bulk(ids)          -> dict # Bulk status fetch
    mark_running(pbs_job_id, ...)          # Update status on job start
    mark_complete(pbs_job_id, exit_code, ...) # Update status on job end
    resolve_virtual_id(virtual_id) -> dict | None
"""

from .db import (
    create_queue_entry,
    get_db_path,
    get_entries_bulk,
    init_db,
    is_virtual_id,
    mark_complete,
    mark_running,
    new_virtual_id,
    resolve_virtual_id,
    update_queue_entry,
)

__all__ = [
    "get_db_path",
    "init_db",
    "new_virtual_id",
    "is_virtual_id",
    "create_queue_entry",
    "update_queue_entry",
    "get_entries_bulk",
    "mark_running",
    "mark_complete",
    "resolve_virtual_id",
]
