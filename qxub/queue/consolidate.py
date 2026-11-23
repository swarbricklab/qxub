"""Occasional consolidation of run files into a published statistics DB.

Phase 3 of ``docs/dev/sqlite-robustness.md``. Jobs record their state in run
files (``qxub.queue.runs``); nothing on a compute node opens SQLite. Now and
then, one process per round folds finished runs into ``qxub-stats.db``:

1. It takes a lease: ``mkdir .qxub-stats.lease`` next to the snapshot.
   ``mkdir`` is atomic across Lustre clients even under ``localflock``. A
   lease older than ``LEASE_STALE_SECONDS`` is broken by renaming it away;
   only one rename wins.
2. It copies the published snapshot to a private file, adds the new runs and
   their resource usage (from the joblog, or one batched ``qstat -fx``), and
   checks ``PRAGMA integrity_check``.
3. It publishes the copy with ``rename``. The published file is never
   written in place, so readers (``open_snapshot``, ``immutable=1``) never
   see a journal and can't roll one back from another node.
4. It deletes run files that the published snapshot holds, once they are
   older than ``GRACE_DAYS``. Unfinished runs are never deleted.

Each user consolidates their own run files into the shared snapshot. A
round that doesn't get the lease skips; a round that dies leaves its run
files for the next one.

Run it with ``qxub db consolidate``, or let ``maybe_start_in_background``
start it after a submission when the last round was ``DUE_SECONDS`` ago.
"""

import base64
import json
import logging
import os
import shutil
import socket
import sqlite3
import stat
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional

from . import runs

logger = logging.getLogger(__name__)

STATS_DB_NAME = "qxub-stats.db"
# The snapshot gets a directory of its own: every user who publishes needs
# write access to it, which the config directory beside it shouldn't grant.
STATS_DIR_NAME = "stats"
SNAPSHOT_VERSION = 1

# A round that holds the lease longer than this is presumed dead.
LEASE_STALE_SECONDS = 3 * 3600
# Start a background round after a submission if the last one was this long ago.
DUE_SECONDS = 6 * 3600
# Keep consolidated run files this long, so `qxtat check` still sees them.
GRACE_DAYS = 7
# A finished run waits this long for its resource usage before it is
# recorded without it (PBS appends usage to the joblog after the job exits).
USAGE_WAIT = timedelta(days=2)
# A run with no sign of starting or ending after this long is given up on.
ABANDON_DAYS = 30
# Runs per round, so one round stays well inside the lease.
BATCH_SIZE = 2000

ATTEMPT_MARKER = ".consolidate.attempted"
LAST_LOG = ".consolidate.last.log"


# -- paths -----------------------------------------------------------------------


def stats_db_path() -> Path:
    """The published snapshot.

    ``QXUB_STATS_DB``, then ``stats_db.path`` from config, then
    ``stats/qxub-stats.db`` beside the old team DB (``shared_db.path``), or
    beside this user's job DB when there is none.
    """
    env = os.environ.get("QXUB_STATS_DB")
    if env:
        return Path(env)
    configured = _config_value("stats_db.path")
    if configured:
        return Path(str(configured))
    from .db import (  # pylint: disable=import-outside-toplevel
        get_db_path,
        legacy_shared_db_path,
    )

    beside = legacy_shared_db_path() or get_db_path()
    return beside.parent / STATS_DIR_NAME / STATS_DB_NAME


def prepare_dir(directory: Path) -> None:
    """Create the snapshot directory so the whole group can publish to it.

    Group-writable and setgid when this process creates it. An existing
    directory is left as it is; ``writable_problem`` reports on it.
    """
    if directory.is_dir():
        return
    directory.mkdir(parents=True, exist_ok=True)
    try:
        mode = directory.stat().st_mode
        os.chmod(directory, mode | stat.S_IRWXG | stat.S_ISGID)
    except OSError as exc:
        logger.debug("Could not make %s group-writable: %s", directory, exc)


def writable_problem(snapshot: Path) -> Optional[str]:
    """Why this user can't publish ``snapshot``, or None if they can."""
    directory = snapshot.parent
    while not directory.exists() and directory.parent != directory:
        directory = directory.parent
    if os.access(directory, os.W_OK | os.X_OK):
        return None
    return (
        f"{directory} is not writable by {os.environ.get('USER') or 'this user'}."
        f" Its owner can run `chmod g+ws {directory}`, or set stats_db.path"
        " to a directory you can write"
    )


def _config_value(key: str):
    try:
        from ..config import config_manager  # pylint: disable=import-outside-toplevel

        return config_manager.get_config_value(key)
    except Exception:  # pylint: disable=broad-except
        return None


def lease_path(snapshot: Path) -> Path:
    return snapshot.with_name(f".{snapshot.name}.lease")


def open_snapshot(path: Optional[Path] = None) -> sqlite3.Connection:
    """Open the published snapshot read-only, without touching journals."""
    path = Path(path or stats_db_path())
    if not path.exists():
        raise FileNotFoundError(path)
    conn = sqlite3.connect(f"{path.resolve().as_uri()}?immutable=1", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


# -- the lease ---------------------------------------------------------------------


class Lease:
    """One consolidator at a time, using ``mkdir`` instead of POSIX locks."""

    def __init__(self, path: Path, stale_seconds: float = LEASE_STALE_SECONDS):
        self.path = Path(path)
        self.stale_seconds = stale_seconds
        self.token = uuid.uuid4().hex
        self.held = False
        self.error: Optional[OSError] = None

    def _owner(self, directory: Path) -> Optional[str]:
        try:
            return (directory / "owner").read_text(encoding="utf-8")
        except OSError:
            return None

    def holder(self) -> Optional[Dict[str, str]]:
        """Who holds the lease, and since when; None if nobody does."""
        if not self.path.is_dir():
            return None
        info = dict(
            line.split("=", 1)
            for line in (self._owner(self.path) or "").splitlines()
            if "=" in line
        )
        info.setdefault("age_seconds", str(int(self._age())))
        return info

    def _age(self) -> float:
        try:
            return time.time() - self.path.stat().st_mtime
        except OSError:
            return 0.0

    def acquire(self) -> bool:
        for _ in range(2):
            try:
                self.path.mkdir()
            except FileExistsError:
                if not self._break_if_stale():
                    return False
                continue
            except OSError as exc:
                self.error = exc
                logger.warning("Cannot take the lease %s: %s", self.path, exc)
                return False
            self.held = True
            try:
                (self.path / "owner").write_text(
                    f"token={self.token}\nhost={socket.gethostname()}\n"
                    f"pid={os.getpid()}\nuser={os.environ.get('USER', '')}\n"
                    f"started={datetime.now().isoformat(timespec='seconds')}\n",
                    encoding="utf-8",
                )
            except OSError:
                pass
            return True
        return False

    def _break_if_stale(self) -> bool:
        if self._age() < self.stale_seconds:
            return False
        seen = self._owner(self.path)
        aside = self.path.with_name(f"{self.path.name}.stale.{uuid.uuid4().hex}")
        try:
            os.rename(self.path, aside)
        except OSError:
            return False  # someone else broke it first
        if self._owner(aside) != seen:
            # A new lease was taken between our look and our rename: put it back.
            try:
                os.rename(aside, self.path)
            except OSError:
                pass
            return False
        logger.warning("Broke a stale consolidation lease (%s)", (seen or "").strip())
        shutil.rmtree(aside, ignore_errors=True)
        return True

    def still_held(self) -> bool:
        """True unless the lease was broken as stale and taken by another."""
        owner = self._owner(self.path)
        return self.held and owner is not None and f"token={self.token}" in owner

    def release(self) -> None:
        if not self.held:
            return
        self.held = False
        owner = self._owner(self.path)
        if owner is not None and f"token={self.token}" not in owner:
            return  # broken and retaken by someone else
        shutil.rmtree(self.path, ignore_errors=True)


# -- resource usage ------------------------------------------------------------------


def _mb(value_bytes) -> Optional[float]:
    return round(value_bytes / (1024 * 1024), 2) if value_bytes else None


def _pct(used, requested) -> Optional[float]:
    if used is None or not requested:
        return None
    return round(used / requested * 100, 1)


def usage_from_joblog(data: Dict[str, Any]) -> Dict[str, Any]:
    """``job_resources`` columns from ``parse_joblog_resources`` output."""
    mem_req = data.get("memory_requested_bytes")
    mem_used = data.get("memory_used_bytes")
    wall_req = data.get("walltime_requested_seconds")
    wall_used = data.get("walltime_used_seconds")
    ncpus = data.get("ncpus_requested")
    cpu_time = data.get("cpu_time_used_seconds")
    jobfs_req = data.get("jobfs_requested_bytes")
    jobfs_used = data.get("jobfs_used_bytes")
    return {
        "mem_requested_mb": _mb(mem_req),
        "mem_used_mb": _mb(mem_used),
        "time_requested_sec": wall_req,
        "time_used_sec": wall_used,
        "cpus_requested": ncpus,
        "cpus_used": data.get("ncpus_used"),
        "jobfs_requested_mb": _mb(jobfs_req),
        "jobfs_used_mb": _mb(jobfs_used),
        "mem_efficiency": _pct(mem_used, mem_req),
        "time_efficiency": _pct(wall_used, wall_req),
        "cpu_efficiency": _pct(cpu_time, (wall_used or 0) * (ncpus or 0)),
        "jobfs_efficiency": _pct(jobfs_used, jobfs_req),
    }


def usage_from_qstat(data: Dict[str, Any]) -> Dict[str, Any]:
    """``job_resources`` columns from ``parse_qstat_fx_output`` output."""
    from ..resources.parser import (  # pylint: disable=import-outside-toplevel
        size_to_bytes,
        time_to_seconds,
    )

    def size(value):
        try:
            return size_to_bytes(value) if value else None
        except Exception:  # pylint: disable=broad-except
            return None

    def seconds(value):
        try:
            return time_to_seconds(value) if value else None
        except Exception:  # pylint: disable=broad-except
            return None

    def integer(value):
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    requested = data.get("resources_requested") or {}
    used = data.get("resources_used") or {}
    efficiency = data.get("efficiency") or {}
    timing = data.get("timing") or {}
    mem_req, mem_used = size(requested.get("mem")), size(used.get("mem"))
    wall_req, wall_used = seconds(requested.get("walltime")), seconds(
        used.get("walltime")
    )
    jobfs_req, jobfs_used = size(requested.get("jobfs")), size(used.get("jobfs"))
    try:
        cpu_percent = float(used["cpupercent"])
    except (KeyError, TypeError, ValueError):
        cpu_percent = None
    ncpus = integer(used.get("ncpus")) or integer(requested.get("ncpus"))
    cpu_time = seconds(used.get("cput"))
    return {
        "mem_requested_mb": _mb(mem_req),
        "mem_used_mb": _mb(mem_used),
        "time_requested_sec": wall_req,
        "time_used_sec": wall_used,
        "cpus_requested": integer(requested.get("ncpus")),
        "cpus_used": integer(used.get("ncpus")),
        "jobfs_requested_mb": _mb(jobfs_req),
        "jobfs_used_mb": _mb(jobfs_used),
        "cpu_percent": cpu_percent,
        "mem_efficiency": efficiency.get("memory_efficiency", _pct(mem_used, mem_req)),
        "time_efficiency": efficiency.get("time_efficiency", _pct(wall_used, wall_req)),
        # As for the joblog: CPU time / (walltime x CPUs). PBS's cpupercent,
        # which the tracker stores as cpu_efficiency, is not divided by CPUs.
        "cpu_efficiency": _pct(cpu_time, (wall_used or 0) * (ncpus or 0)),
        "jobfs_efficiency": efficiency.get(
            "jobfs_efficiency", _pct(jobfs_used, jobfs_req)
        ),
        "queue_wait_sec": timing.get("queue_wait_seconds"),
        "execution_sec": timing.get("execution_seconds"),
    }


def query_qstat(job_ids: Iterable[str]) -> Dict[str, Dict[str, Any]]:
    """What the scheduler remembers about many jobs. Never raises.

    Asks the platform's scheduler backend in batches (PBS: ``qstat -fx``,
    100 jobs per call). Returns parsed data by job ID, for the jobs the
    scheduler still remembers.
    """
    from .. import backends  # pylint: disable=import-outside-toplevel

    return backends.get_backend().job_details(job_ids)


def _qstat_finished(data: Optional[Dict[str, Any]]) -> bool:
    if not data:
        return False
    from .. import backends  # pylint: disable=import-outside-toplevel

    return backends.get_backend().is_finished(data)


# -- the snapshot ---------------------------------------------------------------------


def _create_schema(conn: sqlite3.Connection) -> None:
    from ..resources.tracker import (
        create_job_resources_schema,
    )  # pylint: disable=import-outside-toplevel

    create_job_resources_schema(conn)
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS consolidated_runs (
            run_id TEXT PRIMARY KEY,
            job_id TEXT,
            username TEXT,
            status TEXT,
            consolidated_at TEXT NOT NULL
        )
        """
    )
    conn.execute(
        "CREATE TABLE IF NOT EXISTS snapshot_meta (key TEXT PRIMARY KEY, value TEXT)"
    )
    conn.execute(f"PRAGMA user_version={SNAPSHOT_VERSION}")


@dataclass
class Result:
    """What one round did."""

    snapshot: Path
    skipped: Optional[str] = None
    published: bool = False
    added: int = 0
    settled: int = 0
    abandoned: int = 0
    waiting: int = 0
    in_progress: int = 0
    removed: int = 0
    errors: List[str] = field(default_factory=list)

    def summary(self) -> str:
        if self.skipped:
            return f"Skipped: {self.skipped}"
        parts = [
            f"{self.added} run(s) added",
            f"{self.settled} settled from the joblog or qstat",
            f"{self.abandoned} abandoned",
            f"{self.waiting} waiting for resource usage",
            f"{self.in_progress} still in progress",
            f"{self.removed} run file(s) removed",
        ]
        state = "published" if self.published else "unchanged"
        return f"Snapshot {state}: {self.snapshot}\n" + ", ".join(parts)


def _parse_time(value: Optional[str]) -> Optional[datetime]:
    if not value:
        return None
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return None
    # Run files use local time; `date -Iseconds` adds the offset.
    return parsed.astimezone().replace(tzinfo=None) if parsed.tzinfo else parsed


def _joblog_data(record: Dict[str, str], job_id: Optional[str]):
    """Parsed joblog usage, only if the joblog belongs to this job.

    Joblog names are reused (``{name}.log``), so one from another run of the
    same job name must not be credited to this one.
    """
    path = record.get("joblog")
    if not path or not job_id or not os.path.exists(path):
        return None
    from ..resources.parser import (
        parse_joblog_resources,
    )  # pylint: disable=import-outside-toplevel

    data = parse_joblog_resources(path)
    if data and data.get("job_id") == job_id:
        return data
    return None


def _command(record: Dict[str, str]) -> Optional[str]:
    try:
        return base64.b64decode(record.get("command_b64", "")).decode("utf-8")
    except Exception:  # pylint: disable=broad-except
        return None


def _run_files(runs_dir: Path) -> List[Path]:
    try:
        return sorted(
            Path(e.path)
            for e in os.scandir(runs_dir)
            if e.name.endswith(runs.SUFFIX) and not e.name.startswith(".")
        )
    except OSError:
        return []


def consolidate(
    runs_dir: Optional[Path] = None,
    snapshot: Optional[Path] = None,
    *,
    qstat: Callable[[Iterable[str]], Dict[str, Dict[str, Any]]] = query_qstat,
    now: Optional[datetime] = None,
    grace_days: float = GRACE_DAYS,
    batch_size: int = BATCH_SIZE,
    dry_run: bool = False,
) -> Result:
    """Run one consolidation round. See the module docstring."""
    runs_dir = Path(runs_dir or runs.runs_dir())
    snapshot = Path(snapshot or stats_db_path())
    now = now or datetime.now()
    result = Result(snapshot=snapshot)
    _touch(runs_dir / ATTEMPT_MARKER)
    consolidated: Dict[str, str] = {}

    def add(conn) -> bool:
        done = {
            r[0]: r[1]
            for r in conn.execute(
                "SELECT run_id, consolidated_at FROM consolidated_runs"
            )
        }
        changed = _add_runs(
            conn, runs_dir, done, now, qstat, batch_size, result, dry_run
        )
        consolidated.update(
            (r[0], r[1])
            for r in conn.execute(
                "SELECT run_id, consolidated_at FROM consolidated_runs"
            )
        )
        return changed

    if _round(snapshot, now, result, add, dry_run) and not dry_run:
        try:
            result.removed = _remove_consolidated(
                runs_dir, consolidated, now - timedelta(days=grace_days)
            )
        except OSError as exc:
            result.errors.append(str(exc))
    return result


def _round(snapshot: Path, now: datetime, result, add, dry_run: bool) -> bool:
    """Take the lease, apply ``add(conn)`` to a private copy, publish it.

    ``add`` returns True if it changed anything; only then is a copy
    published. Records a skip or an error on ``result``; returns True if the
    round finished (published or unchanged).
    """
    problem = writable_problem(snapshot)
    if problem:
        result.skipped = problem
        return False
    try:
        prepare_dir(snapshot.parent)
    except OSError as exc:
        result.skipped = f"cannot create {snapshot.parent}: {exc}"
        return False
    lease = Lease(lease_path(snapshot))
    if not lease.acquire():
        if lease.error is not None:
            result.skipped = f"cannot take the lease: {lease.error}"
            return False
        holder = lease.holder() or {}
        result.skipped = "another consolidation is running ({}@{})".format(
            holder.get("user", "?"), holder.get("host", "?")
        )
        return False

    private = snapshot.with_name(f".{snapshot.name}.{uuid.uuid4().hex}.tmp")
    try:
        _remove_old_private_copies(snapshot)
        if snapshot.exists():
            shutil.copyfile(snapshot, private)
        conn = sqlite3.connect(private, isolation_level=None)
        try:
            conn.execute("PRAGMA journal_mode=MEMORY")
            conn.execute("BEGIN")
            _create_schema(conn)
            changed = add(conn)
            check = conn.execute("PRAGMA integrity_check").fetchone()[0]
            if check != "ok":
                raise sqlite3.DatabaseError(f"integrity_check: {check}")
            if changed:
                conn.execute(
                    "INSERT OR REPLACE INTO snapshot_meta VALUES ('published_at', ?)",
                    (now.isoformat(timespec="seconds"),),
                )
            conn.execute("COMMIT")
        finally:
            conn.close()

        if changed and not dry_run:
            if not lease.still_held():
                # This round outlived LEASE_STALE_SECONDS and another took
                # over; publishing now would drop that round's work.
                raise OSError(f"lost the lease {lease.path}; not publishing")
            _publish(private, snapshot)
            result.published = True
        return True
    except (OSError, sqlite3.Error) as exc:
        result.errors.append(str(exc))
        logger.warning("Publishing %s failed: %s", snapshot, exc)
        return False
    finally:
        try:
            private.unlink()
        except OSError:
            pass
        lease.release()


def _add_runs(conn, runs_dir, done, now, qstat, batch_size, result, dry_run) -> bool:
    """Add finished runs to the private copy. True if anything was added."""
    pending = []
    for path in _run_files(runs_dir):
        run_id = path.name[: -len(runs.SUFFIX)]
        if run_id in done:
            continue
        record = runs.read_run(path)
        if record is None:
            continue
        pending.append((run_id, path, record))
        if len(pending) >= batch_size:
            break

    # Usage from joblogs first; one batched qstat for the rest.
    joblogs = {}
    need_qstat = []
    for run_id, _path, record in pending:
        job_id = record.get("native_job_id") or None
        joblogs[run_id] = _joblog_data(record, job_id)
        if job_id and not joblogs[run_id]:
            need_qstat.append(job_id)
    pbs = qstat(need_qstat) if need_qstat else {}

    changed = False
    stamp = now.isoformat(timespec="seconds")
    for run_id, path, record in pending:
        job_id = record.get("native_job_id") or None
        joblog = joblogs[run_id]
        from_pbs = pbs.get(job_id) if job_id else None
        state = record.get("state")
        status = None

        if state != "finished":
            if joblog and joblog.get("exit_status") is not None:
                code = joblog["exit_status"]
            elif _qstat_finished(from_pbs) and from_pbs.get("exit_status") is not None:
                code = from_pbs["exit_status"]
            else:
                code = None
            if code is not None:
                # Killed without running the EXIT trap (e.g. SIGKILL).
                if not dry_run:
                    _settle_run_file(path, code, stamp)
                record = {**record, "state": "finished", "exit_code": str(code)}
                result.settled += 1
            else:
                started = _parse_time(record.get("submitted_at")) or _mtime(path)
                if started and now - started > timedelta(days=ABANDON_DAYS):
                    status = "abandoned"
                else:
                    result.in_progress += 1
                    continue

        if status is None:
            usage = None
            if joblog:
                usage = usage_from_joblog(joblog)
            elif from_pbs and from_pbs.get("resources_used"):
                usage = usage_from_qstat(from_pbs)
            finished = _parse_time(record.get("finished_at")) or _mtime(path)
            if usage is None and finished and now - finished < USAGE_WAIT:
                result.waiting += 1
                continue
            code = runs.exit_code(record)
            status = "completed" if code == 0 else "failed"
            if job_id:
                _insert_job(conn, job_id, record, status, code, usage or {}, stamp)

        conn.execute(
            "INSERT OR REPLACE INTO consolidated_runs VALUES (?, ?, ?, ?, ?)",
            (run_id, job_id, record.get("user"), status, stamp),
        )
        if status == "abandoned":
            result.abandoned += 1
        else:
            result.added += 1
        changed = True
    return changed


def _insert_job(conn, job_id, record, status, code, usage, stamp) -> None:
    try:
        tags = json.dumps(json.loads(record.get("tags") or "[]"))
    except ValueError:
        tags = "[]"
    row = {
        "job_id": job_id,
        "timestamp": record.get("finished_at") or stamp,
        "command": _command(record),
        "exit_code": code,
        "status": status,
        "submitted_at": record.get("submitted_at") or None,
        "started_at": record.get("started_at") or None,
        "completed_at": record.get("finished_at") or None,
        "last_status_update": stamp,
        "queue": record.get("queue") or None,
        "hostname": record.get("host") or None,
        "tags": tags,
        "username": record.get("user") or None,
        "joblog_path": record.get("joblog") or None,
        **{k: v for k, v in usage.items() if v is not None},
    }
    columns = ", ".join(row)
    marks = ", ".join("?" * len(row))
    conn.execute(
        f"INSERT OR REPLACE INTO job_resources ({columns}) VALUES ({marks})",
        list(row.values()),
    )


def _settle_run_file(path: Path, code: int, stamp: str) -> None:
    """Record an exit the job couldn't record itself. Atomic, like the job's."""
    tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        text = path.read_text(encoding="utf-8")
        tmp.write_text(
            text
            + runs.format_fields(
                {
                    "state": "finished",
                    "exit_code": code,
                    "settled_by": "consolidate",
                    "settled_at": stamp,
                }
            ),
            encoding="utf-8",
        )
        os.replace(tmp, path)
    except OSError as exc:
        logger.debug("Could not settle %s: %s", path, exc)
        try:
            tmp.unlink()
        except OSError:
            pass


def _publish(private: Path, snapshot: Path) -> None:
    try:
        if snapshot.exists():
            shutil.copymode(snapshot, private)
        else:
            os.chmod(private, private.stat().st_mode | stat.S_IRGRP)
    except OSError:
        pass
    with open(private, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(private, snapshot)


def _remove_consolidated(runs_dir: Path, consolidated: Dict[str, str], cutoff) -> int:
    """Remove this user's run files that the snapshot has held since ``cutoff``.

    Walks the run files that exist, not the snapshot, so the cost stays
    proportional to what is on disk.
    """
    removed = 0
    for path in _run_files(runs_dir):
        at = _parse_time(consolidated.get(path.name[: -len(runs.SUFFIX)]))
        if not at or at > cutoff:
            continue
        try:
            path.unlink()
            removed += 1
        except FileNotFoundError:
            pass
        except OSError as exc:
            logger.debug("Could not remove run file %s: %s", path, exc)
    return removed


def _remove_old_private_copies(snapshot: Path) -> None:
    """Private copies left by rounds that died (only ours can be removed)."""
    cutoff = time.time() - LEASE_STALE_SECONDS
    for path in snapshot.parent.glob(f".{snapshot.name}.*.tmp"):
        try:
            if path.stat().st_mtime < cutoff:
                path.unlink()
        except OSError:
            pass


def _mtime(path: Path) -> Optional[datetime]:
    try:
        return datetime.fromtimestamp(path.stat().st_mtime)
    except OSError:
        return None


def _touch(path: Path) -> bool:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        return True
    except OSError:
        return False


# -- when it runs ---------------------------------------------------------------------


def auto_enabled() -> bool:
    if os.environ.get("QXUB_AUTO_CONSOLIDATE", "").lower() in ("0", "false", "no"):
        return False
    return _config_value("stats_db.auto_consolidate") is not False


def is_due(runs_dir: Path, due_seconds: float = DUE_SECONDS) -> bool:
    """True if the last attempt was ``due_seconds`` ago.

    With no attempt yet, the clock starts now instead, so the first
    submission doesn't start a round.
    """
    marker = runs_dir / ATTEMPT_MARKER
    try:
        return time.time() - marker.stat().st_mtime >= due_seconds
    except FileNotFoundError:
        _touch(marker)
        return False
    except OSError:
        return False


def maybe_start_in_background(runs_dir: Optional[Path] = None) -> bool:
    """Start a detached round if one is due. Cheap; never raises or waits."""
    try:
        if not auto_enabled():
            return False
        runs_dir = Path(runs_dir or runs.runs_dir())
        if not is_due(runs_dir) or not _touch(runs_dir / ATTEMPT_MARKER):
            return False
        with open(runs_dir / LAST_LOG, "w", encoding="utf-8") as log:
            subprocess.Popen(  # pylint: disable=consider-using-with
                [sys.executable, "-m", "qxub.queue.consolidate"],
                stdin=subprocess.DEVNULL,
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
                close_fds=True,
                # Not the user's cwd: `-m` would import a qxub checkout there.
                cwd=str(runs_dir),
                env={**os.environ, "QXUB_AUTO_CONSOLIDATE": "0"},
            )
        return True
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("Could not start background consolidation: %s", exc)
        return False


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    result = consolidate()
    print(f"{datetime.now().isoformat(timespec='seconds')} {result.summary()}")
    return 1 if result.errors else 0


if __name__ == "__main__":
    sys.exit(main())
