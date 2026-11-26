"""Occasional consolidation of run files into a published statistics DB.

Phase 3 of ``docs/dev/sqlite-robustness.md``. Jobs record their state in run
files (``qxub.queue.runs``); nothing on a compute node opens SQLite. Now and
then, one process per round folds finished runs into ``qxub-stats.db``:

1. It takes a lease: ``mkdir .qxub-stats.lease`` next to the snapshot.
   ``mkdir`` is atomic across Lustre clients even under ``localflock``. A
   lease older than ``LEASE_STALE_SECONDS`` is broken by renaming it away;
   only one rename wins.
2. It copies the published snapshot to a private file (in a new ``0700``
   temporary directory, never in the shared directory, and only if the
   snapshot is a regular file, not a link), adds the new runs and
   their resource usage (from the joblog, or one batched query of the scheduler), and
   checks ``PRAGMA integrity_check``.
3. It publishes the copy: written to a new file in the shared directory
   through the descriptor that created it, then ``rename``d over the
   snapshot. The published file is never written in place, so readers (``open_snapshot``, ``immutable=1``) never
   see a journal and can't roll one back from another node.
4. It marks the run files it published (``consolidated_at``), and deletes
   run files marked more than ``GRACE_DAYS`` ago that the snapshot holds.
   Unfinished runs are never deleted. The snapshot can be replaced by
   anyone in the group, so its rows alone never skip or delete a run file:
   only the mark this user wrote in their own run file does.

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
import tempfile
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
# The published snapshot's mode: the group reads it, only its publisher
# writes it (others replace it by rename). Never copied from the old one.
SNAPSHOT_MODE = 0o640

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
# Written into a run file once a round of this user's has published it.
CONSOLIDATED_AT = "consolidated_at"
CONSOLIDATED_AS = "consolidated_as"
LAST_LOG = ".consolidate.last.log"

_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)


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
        # Others can write the shared directory: never follow a link here.
        try:
            fd = os.open(directory / "owner", os.O_RDONLY | _NOFOLLOW | os.O_NONBLOCK)
        except OSError:
            return None
        try:
            if not stat.S_ISREG(os.fstat(fd).st_mode):
                return None
            return os.read(fd, 4096).decode("utf-8", "replace")
        except OSError:
            return None
        finally:
            os.close(fd)

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
                # A new file only: O_EXCL never follows a link put in its place.
                fd = os.open(
                    self.path / "owner",
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW,
                    0o644,
                )
                with os.fdopen(fd, "w", encoding="utf-8") as fh:
                    fh.write(
                        f"token={self.token}\nhost={socket.gethostname()}\n"
                        f"pid={os.getpid()}\nuser={os.environ.get('USER', '')}\n"
                        f"started={datetime.now().isoformat(timespec='seconds')}\n"
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


def usage_from_info(info) -> Dict[str, Any]:
    """``job_resources`` columns from a backend's :class:`~qxub.backends.JobInfo`."""
    usage = info.usage
    wall_used = usage.walltime_used
    ncpus = usage.cpus_used or usage.cpus_requested
    return {
        "mem_requested_mb": _mb(usage.memory_requested),
        "mem_used_mb": _mb(usage.memory_used),
        "time_requested_sec": usage.walltime_requested,
        "time_used_sec": wall_used,
        "cpus_requested": usage.cpus_requested,
        "cpus_used": usage.cpus_used,
        "jobfs_requested_mb": _mb(usage.disk_requested),
        "jobfs_used_mb": _mb(usage.disk_used),
        "cpu_percent": usage.cpu_percent,
        "mem_efficiency": _pct(usage.memory_used, usage.memory_requested),
        "time_efficiency": _pct(wall_used, usage.walltime_requested),
        # As for the joblog: CPU time / (walltime x CPUs). PBS's cpupercent,
        # which the tracker stores as cpu_efficiency, is not divided by CPUs.
        "cpu_efficiency": _pct(usage.cpu_time, (wall_used or 0) * (ncpus or 0)),
        "jobfs_efficiency": _pct(usage.disk_used, usage.disk_requested),
        "queue_wait_sec": info.queue_wait,
        "execution_sec": info.run_time,
    }


def query_jobs(job_ids: Iterable[str]) -> Dict[str, Any]:
    """What the schedulers remember about many jobs (:class:`JobInfo` by
    job ID). Never raises.

    Asks each job's backend in batches (PBS: ``qstat -fx``, 100 jobs per
    call). ``job_ids`` may be a dict of job ID to the backend its run file
    recorded; otherwise the backend is inferred from the ID. Returns parsed
    data by job ID, for the jobs the scheduler still remembers.
    """
    from .. import backends  # pylint: disable=import-outside-toplevel

    # Each job asks its own backend (PBS: qstat -fx; local: its records)
    by_backend: Dict[Any, List[str]] = {}
    owners: Dict[str, Any] = {}
    recorded = job_ids if isinstance(job_ids, dict) else {}
    for job_id in job_ids:
        backend = backends.backend_for_job(job_id, recorded.get(job_id))
        key = type(backend)
        owners.setdefault(key, backend)
        by_backend.setdefault(key, []).append(job_id)
    found: Dict[str, Dict[str, Any]] = {}
    for key, ids in by_backend.items():
        found.update(owners[key].job_info(ids))
    return found


# -- the snapshot ---------------------------------------------------------------------


def _create_schema(conn: sqlite3.Connection) -> None:
    from ..resources.tracker import (  # pylint: disable=import-outside-toplevel
        create_job_resources_schema,
    )

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
            f"{self.settled} settled from the joblog or scheduler",
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
    from ..resources.parser import (  # pylint: disable=import-outside-toplevel
        parse_joblog_resources,
    )

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
    lookup: Callable[[Iterable[str]], Dict[str, Any]] = query_jobs,
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
    # Run files added this round, to mark once they are published, and
    # run files the snapshot already held, to mark without re-adding.
    added: Dict[Path, str] = {}
    adopted: Dict[Path, str] = {}

    def add(conn) -> bool:
        done = {
            r[0]: r[1]
            for r in conn.execute("SELECT run_id, job_id FROM consolidated_runs")
        }
        changed = _add_runs(
            conn,
            runs_dir,
            done,
            now,
            lookup,
            batch_size,
            result,
            dry_run,
            added,
            adopted,
        )
        consolidated.update(
            (r[0], r[1])
            for r in conn.execute(
                "SELECT run_id, consolidated_at FROM consolidated_runs"
            )
        )
        return changed

    def mark() -> None:
        # Still under the lease, so an overlapping round can't add these
        # runs again between the publish and the mark.
        stamp = now.isoformat(timespec="seconds")
        marks = {**adopted, **(added if result.published else {})}
        for path, status in marks.items():
            runs.update(path, {CONSOLIDATED_AT: stamp, CONSOLIDATED_AS: status})

    if _round(snapshot, now, result, add, dry_run, mark) and not dry_run:
        try:
            result.removed = _remove_consolidated(
                runs_dir, consolidated, now - timedelta(days=grace_days)
            )
        except OSError as exc:
            result.errors.append(str(exc))
    return result


def _round(
    snapshot: Path, now: datetime, result, add, dry_run: bool, done=None
) -> bool:
    """Take the lease, apply ``add(conn)`` to a private copy, publish it.

    ``add`` returns True if it changed anything; only then is a copy
    published. ``done()``, if given, runs once the round has finished
    (published or unchanged), before the lease is released. Records a skip
    or an error on ``result``; returns True if the round finished.
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

    # Work where only this user can see or rename anything: the shared
    # directory is writable by the whole group and has no sticky bit.
    workdir = None
    try:
        _remove_old_private_copies(snapshot)
        workdir = Path(tempfile.mkdtemp(prefix="qxub-consolidate-"))
        private = workdir / STATS_DB_NAME
        _copy_snapshot(snapshot, private)
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
        if done is not None and not dry_run:
            done()
        return True
    except (OSError, sqlite3.Error) as exc:
        result.errors.append(str(exc))
        logger.warning("Publishing %s failed: %s", snapshot, exc)
        return False
    finally:
        if workdir is not None:
            shutil.rmtree(workdir, ignore_errors=True)
        lease.release()


def _copy_snapshot(snapshot: Path, private: Path) -> None:
    """Copy the published snapshot to ``private``, if there is one.

    Others can replace the snapshot, so it is opened without following a
    link and read only if that open file is a regular file: never the file
    a link points to, nor a FIFO or device.
    """
    try:
        fd = os.open(snapshot, os.O_RDONLY | _NOFOLLOW | os.O_NONBLOCK)
    except FileNotFoundError:
        return
    except OSError as exc:
        if os.path.islink(snapshot):
            raise OSError(f"{snapshot} is a symbolic link; not reading it") from exc
        raise
    try:
        if not stat.S_ISREG(os.fstat(fd).st_mode):
            raise OSError(f"{snapshot} is not a regular file; not reading it")
        out = os.open(private, os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW, 0o600)
        with os.fdopen(fd, "rb", closefd=False) as src, os.fdopen(out, "wb") as dst:
            shutil.copyfileobj(src, dst)
    finally:
        os.close(fd)


def _add_runs(
    conn, runs_dir, done, now, lookup, batch_size, result, dry_run, added, adopted
) -> bool:
    """Add finished runs to the private copy. True if anything was added.

    ``done`` maps the snapshot's consolidated run IDs to their job IDs.
    Records in ``added`` the run files it added (path -> status). A run is
    skipped as already published only if the snapshot holds it *and* its
    own file carries this user's mark: anyone in the group can put rows in
    the snapshot.

    A run published before run files were marked (by an older qxub) is not
    added again if its own file says it finished and the snapshot holds its
    job: re-adding it now, after its joblog may be gone and the scheduler
    has forgotten it, would lose its usage. It goes in ``adopted``, to be
    marked. Such a run is past running, so a planted row can at most keep
    its usage out of the statistics, as replacing the snapshot could anyway.
    """
    pending = []
    for path in _run_files(runs_dir):
        run_id = path.name[: -len(runs.SUFFIX)]
        record = runs.read_run(path)
        if record is None:
            continue
        if run_id in done and record.get(CONSOLIDATED_AT):
            continue
        if run_id in done and _held_finished(conn, done[run_id], record):
            adopted[path] = "finished"
            continue
        pending.append((run_id, path, record))
        if len(pending) >= batch_size:
            break

    # Usage from joblogs first; one batched scheduler query for the rest.
    joblogs = {}
    ask: Dict[str, Optional[str]] = {}  # job ID -> recorded backend
    for run_id, _path, record in pending:
        job_id = record.get("native_job_id") or None
        joblogs[run_id] = _joblog_data(record, job_id)
        if job_id and not joblogs[run_id]:
            ask[job_id] = record.get("backend") or None
    known = lookup(ask) if ask else {}

    changed = False
    stamp = now.isoformat(timespec="seconds")
    for run_id, path, record in pending:
        job_id = record.get("native_job_id") or None
        joblog = joblogs[run_id]
        info = known.get(job_id) if job_id else None
        state = record.get("state")
        status = None

        if state != "finished":
            if joblog and joblog.get("exit_status") is not None:
                code = joblog["exit_status"]
            elif info and info.is_finished and info.exit_status is not None:
                code = info.exit_status
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
            elif info and info.usage.has_used():
                usage = usage_from_info(info)
            finished = _parse_time(record.get("finished_at")) or _mtime(path)
            if usage is None and finished and now - finished < USAGE_WAIT:
                result.waiting += 1
                continue
            code = runs.exit_code(record)
            status = "completed" if code == 0 else "failed"
            if job_id:
                usage = {
                    **(usage or {}),
                    "ended_by": _ended_by(job_id, record, joblog, info),
                }
                _insert_job(conn, job_id, record, status, code, usage, stamp)

        conn.execute(
            "INSERT OR REPLACE INTO consolidated_runs VALUES (?, ?, ?, ?, ?)",
            (run_id, job_id, record.get("user"), status, stamp),
        )
        added[path] = "abandoned" if status == "abandoned" else "finished"
        if status == "abandoned":
            result.abandoned += 1
        else:
            result.added += 1
        changed = True
    return changed


def _held_finished(conn, held_job_id, record) -> bool:
    """True if the run file says it finished and the snapshot already has a
    ``job_resources`` row for the same job as the run file names."""
    if record.get("state") != "finished":
        return False
    job_id = record.get("native_job_id") or None
    if job_id is None or held_job_id != job_id:
        return False
    return (
        conn.execute(
            "SELECT 1 FROM job_resources WHERE job_id = ?", (job_id,)
        ).fetchone()
        is not None
    )


def _ended_by(job_id, record, joblog, info) -> Optional[str]:
    """Why the job ended: as its run file says (``qxub cancel``, or qx
    from the job log), the scheduler said, or its backend reads the job
    log's exit status (PBS: -29 for walltime; the run file then holds 143)."""
    if record.get("ended_by"):
        return record["ended_by"]
    if info is not None and info.ended_by:
        return info.ended_by
    from .. import backends  # pylint: disable=import-outside-toplevel

    backend = backends.backend_for_job(job_id, record.get("backend") or None)
    for code in ((joblog or {}).get("exit_status"), runs.exit_code(record)):
        if code:
            reason = backend.ended_by(code)
            if reason:
                return reason
    return None


def _backend_name(job_id: str) -> str:
    """For run files written before they recorded the backend."""
    from .. import backends  # pylint: disable=import-outside-toplevel

    return backends.backend_for_job(job_id).name


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
        "backend": record.get("backend") or _backend_name(job_id),
        **{k: v for k, v in usage.items() if v is not None},
    }
    columns = ", ".join(row)
    marks = ", ".join("?" * len(row))
    # A job added again (its usage now unknown) keeps the usage it had.
    keep = ", ".join(f"{c} = COALESCE(excluded.{c}, {c})" for c in row if c != "job_id")
    conn.execute(
        f"INSERT INTO job_resources ({columns}) VALUES ({marks})"
        f" ON CONFLICT(job_id) DO UPDATE SET {keep}",
        list(row.values()),
    )


def _settle_run_file(path: Path, code: int, stamp: str) -> None:
    """Record an exit the job couldn't record itself. Atomic, like the job's."""
    runs.update(
        path,
        {
            "state": "finished",
            "exit_code": code,
            "settled_by": "consolidate",
            "settled_at": stamp,
        },
    )


def _publish(private: Path, snapshot: Path) -> None:
    """Put the finished ``private`` copy in place of ``snapshot``.

    The copy is written to a new file in the shared directory, created
    exclusively (so never through a link), and given ``SNAPSHOT_MODE`` and
    synced through the descriptor that created it; then renamed over the
    snapshot. The mode is qxub's own, never the old snapshot's: anyone in
    the group may own that. If the staged name no longer refers to the
    file this process created, nothing is published.
    """
    staged = snapshot.with_name(f".{snapshot.name}.{uuid.uuid4().hex}.tmp")
    fd = os.open(staged, os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW, 0o600)
    try:
        mine = os.fstat(fd)
    except OSError:
        os.close(fd)
        raise
    renamed = False
    try:
        with open(private, "rb") as src, os.fdopen(fd, "wb", closefd=False) as dst:
            shutil.copyfileobj(src, dst)
        _give_to_dir_group(fd, mine, snapshot.parent)
        os.fchmod(fd, SNAPSHOT_MODE)
        os.fsync(fd)
        if not _same_file(staged, mine):
            raise OSError(f"{staged} was replaced by someone else; not publishing")
        os.replace(staged, snapshot)
        renamed = True
    finally:
        os.close(fd)
        if not renamed and _same_file(staged, mine):
            try:
                staged.unlink()
            except OSError:
                pass


# Whether _give_to_dir_group has warned that it can't set the group.
_GROUP_WARNED = False


def _give_to_dir_group(fd: int, mine: os.stat_result, directory: Path) -> None:
    """Give the staged snapshot ``fd`` the group of ``directory``.

    ``SNAPSHOT_MODE`` lets the group read it, which is only the team's group
    if the directory is setgid; without that the file has the user's own
    group. If the user isn't in the directory's group (``EPERM``), the file
    keeps its group and a warning is logged once.
    """
    global _GROUP_WARNED  # pylint: disable=global-statement
    gid = os.stat(directory).st_gid
    if mine.st_gid == gid:
        return
    try:
        os.fchown(fd, -1, gid)
    except PermissionError as exc:
        if not _GROUP_WARNED:
            _GROUP_WARNED = True
            logger.warning(
                "Cannot give the statistics snapshot the group of %s (gid %d): "
                "%s; others in that group may not be able to read it",
                directory,
                gid,
                exc,
            )


def _same_file(path: Path, known: os.stat_result) -> bool:
    """True if ``path`` itself (not a link's target) is the file ``known``."""
    try:
        there = os.lstat(path)
    except OSError:
        return False
    return (there.st_dev, there.st_ino) == (known.st_dev, known.st_ino)


def _remove_consolidated(runs_dir: Path, consolidated: Dict[str, str], cutoff) -> int:
    """Remove this user's run files consolidated before ``cutoff``.

    Only a run file that this user's own round marked as published before
    ``cutoff``, that says it finished (or was given up on), and that the
    snapshot still holds. The snapshot alone is never enough: anyone in
    the group can replace it.

    Walks the run files that exist, not the snapshot, so the cost stays
    proportional to what is on disk.
    """
    removed = 0
    for path in _run_files(runs_dir):
        if path.name[: -len(runs.SUFFIX)] not in consolidated:
            continue
        record = runs.read_run(path) or {}
        at = _parse_time(record.get(CONSOLIDATED_AT))
        if not at or at > cutoff:
            continue
        if record.get("state") != "finished" and (
            record.get(CONSOLIDATED_AS) != "abandoned"
        ):
            continue
        try:
            path.unlink()
            removed += 1
        except FileNotFoundError:
            pass
        except OSError as exc:
            logger.debug("Could not remove run file %s: %s", path, exc)
        runs.discard(path.with_suffix(runs.SETTINGS_SUFFIX))
    return removed


def _remove_old_private_copies(snapshot: Path) -> None:
    """Staged copies left by publishes that died (only ours can be removed)."""
    cutoff = time.time() - LEASE_STALE_SECONDS
    for path in snapshot.parent.glob(f".{snapshot.name}.*.tmp"):
        try:
            if path.lstat().st_mtime < cutoff:
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
