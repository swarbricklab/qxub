"""Import the old team-wide job DB into the statistics snapshot, once.

Phase 4 of ``docs/dev/sqlite-robustness.md``. Before run files, every login
node wrote ``shared_db.path`` directly, and on Lustre it ended up corrupt.
``recover`` reads the job rows it still holds and adds them to the published
snapshot with the same lease and publish steps as a consolidation round
(``consolidate._round``):

1. The old DB is copied to a private temporary directory and opened with
   ``immutable=1``, so nothing here can write it or roll back a journal.
2. If SQLite can't read it all ("database disk image is malformed"), the
   ``sqlite3`` command's ``.recover`` rebuilds what it can into a new file,
   and the rows come from that.
3. Rows the corruption garbled are left out (``malformed``): a job ID that
   isn't a PBS ID, a time that isn't ISO, or text in a numeric column.
4. A job still marked submitted or running but with an exit code lost its
   last status update; it is stored as completed or failed from that code.
5. Only finished jobs, or jobs with resource usage, are added: other jobs
   still marked submitted or running are stale or will arrive from run
   files. A job the snapshot already has is never overwritten.

Running it again adds nothing new, so it is safe to repeat.
"""

import logging
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .consolidate import _round, stats_db_path

logger = logging.getLogger(__name__)

FINISHED = ("completed", "failed", "cancelled")
# 12345678.gadi-pbs, 123[].gadi-pbs, 123[4].gadi-pbs, or a bare number.
PBS_JOB_ID = re.compile(r"\d+(\[\d*\])?(\.[A-Za-z0-9][A-Za-z0-9.-]*)?")
TIME_COLUMNS = ("timestamp", "submitted_at", "started_at", "completed_at")
RECOVERED_FROM = "recovered_from"


@dataclass
class RecoverResult:
    """What one recovery did."""

    snapshot: Path
    source: Optional[Path] = None
    method: str = ""
    read: int = 0
    added: int = 0
    present: int = 0
    unfinished: int = 0
    malformed: int = 0
    status_from_exit: int = 0
    skipped: Optional[str] = None
    published: bool = False
    errors: List[str] = field(default_factory=list)

    def summary(self) -> str:
        if self.skipped:
            return f"Skipped: {self.skipped}"
        state = "published" if self.published else "unchanged"
        return (
            f"Snapshot {state}: {self.snapshot}\n"
            f"Read {self.read} job(s) from {self.source} ({self.method}): "
            f"{self.added} added, {self.present} already there, "
            f"{self.malformed} malformed and {self.unfinished} unfinished left out"
            f"; {self.status_from_exit} marked completed/failed from their exit code"
        )


def recover(
    source: Path,
    snapshot: Optional[Path] = None,
    *,
    now: Optional[datetime] = None,
    dry_run: bool = False,
) -> RecoverResult:
    """Add the job rows ``source`` still holds to the snapshot."""
    source = Path(source)
    snapshot = Path(snapshot or stats_db_path())
    now = now or datetime.now()
    result = RecoverResult(snapshot=snapshot, source=source)
    if not source.is_file():
        result.skipped = f"{source} does not exist"
        return result
    if source.resolve() == snapshot.resolve():
        result.skipped = f"{source} is the snapshot itself"
        return result

    with tempfile.TemporaryDirectory(prefix="qxub-recover-") as tmp:
        try:
            columns, rows, result.method = read_jobs(source, Path(tmp))
        except (OSError, sqlite3.Error, RuntimeError) as exc:
            result.errors.append(f"cannot read {source}: {exc}")
            return result
    result.read = len(rows)

    def add(conn) -> bool:
        types = _snapshot_columns(conn)
        wanted = [c for c in columns if c in types]
        if "job_id" not in wanted:
            raise sqlite3.DatabaseError(f"{source} has no job_id column")
        numeric = {c for c in wanted if types[c] in ("INTEGER", "REAL")}
        names = ", ".join(wanted)
        marks = ", ".join("?" * len(wanted))
        for raw in rows:
            row = {c: v for c, v in zip(columns, raw) if c in types}
            if not _well_formed(row, numeric):
                result.malformed += 1
                continue
            if _status_from_exit_code(row):
                result.status_from_exit += 1
            if not _worth_keeping(row):
                result.unfinished += 1
                continue
            cur = conn.execute(
                f"INSERT OR IGNORE INTO job_resources ({names}) VALUES ({marks})",
                [row[c] for c in wanted],
            )
            if cur.rowcount:
                result.added += 1
            else:
                result.present += 1
        if not result.added:
            return False
        conn.execute(
            "INSERT OR REPLACE INTO snapshot_meta VALUES (?, ?)",
            (
                RECOVERED_FROM,
                f"{source} at {now.isoformat(timespec='seconds')}"
                f" ({result.added} jobs, {result.method})",
            ),
        )
        return True

    _round(snapshot, now, result, add, dry_run)
    return result


def read_jobs(source: Path, workdir: Path) -> Tuple[List[str], List[tuple], str]:
    """Columns and rows of ``job_resources`` in ``source``, and how they were read.

    Reads a private copy. Falls back to ``sqlite3 .recover`` if the copy is
    too damaged to read straight through.
    """
    copy = workdir / "source.db"
    shutil.copyfile(source, copy)
    try:
        columns, rows = _read_table(copy)
        return columns, rows, "read directly"
    except sqlite3.DatabaseError as exc:
        logger.info("Reading %s failed (%s); trying sqlite3 .recover", source, exc)
        first_error = exc
    rebuilt = workdir / "recovered.db"
    _sqlite3_recover(copy, rebuilt, first_error)
    columns, rows = _read_table(rebuilt)
    return columns, rows, "rebuilt with sqlite3 .recover"


def _read_table(path: Path) -> Tuple[List[str], List[tuple]]:
    conn = sqlite3.connect(f"{path.resolve().as_uri()}?immutable=1", uri=True)
    try:
        columns = [r[1] for r in conn.execute("PRAGMA table_info(job_resources)")]
        if not columns:
            raise sqlite3.DatabaseError("no job_resources table")
        rows = conn.execute("SELECT * FROM job_resources").fetchall()
        return columns, rows
    finally:
        conn.close()


def sqlite3_cli() -> Optional[str]:
    """The ``sqlite3`` command: the one beside this Python first, then PATH.

    qxub's conda environment pins ``sqlite>=3.29`` (``.recover`` arrived in
    3.29), and its ``bin`` may not be on PATH when qxub is run by full path.
    """
    beside = Path(sys.executable).parent / "sqlite3"
    if beside.is_file() and os.access(beside, os.X_OK):
        return str(beside)
    return shutil.which("sqlite3")


def _sqlite3_recover(copy: Path, rebuilt: Path, first_error: Exception) -> None:
    exe = sqlite3_cli()
    if not exe:
        raise RuntimeError(
            f"{first_error}, and there is no sqlite3 command to run .recover."
            " Install qxub's conda environment (it includes sqlite>=3.29)"
            " or put a sqlite3 of 3.29 or later on PATH"
        )
    dump = subprocess.run(
        [exe, str(copy), ".recover"],
        capture_output=True,
        text=True,
        check=False,
    )
    if dump.returncode != 0 or not dump.stdout:
        raise RuntimeError(
            f"{first_error}; sqlite3 .recover failed: "
            f"{(dump.stderr or '').strip() or 'no output'}"
        )
    load = subprocess.run(
        [exe, str(rebuilt)],
        input=dump.stdout,
        capture_output=True,
        text=True,
        check=False,
    )
    if load.returncode != 0 and not rebuilt.exists():
        raise RuntimeError(
            f"{first_error}; loading the .recover output failed: "
            f"{load.stderr.strip()}"
        )


def _snapshot_columns(conn) -> Dict[str, str]:
    """Column name -> declared type of the snapshot's ``job_resources``."""
    return {
        r[1]: (r[2] or "").upper()
        for r in conn.execute("PRAGMA table_info(job_resources)")
    }


def _well_formed(row: Dict[str, Any], numeric) -> bool:
    """False for rows the corruption mixed up with other tables or pages."""
    job_id = row.get("job_id")
    if not isinstance(job_id, str) or not PBS_JOB_ID.fullmatch(job_id):
        return False
    if not isinstance(row.get("timestamp"), str):
        return False
    for column in TIME_COLUMNS:
        value = row.get(column)
        if value is None:
            continue
        if not isinstance(value, str):
            return False
        try:
            datetime.fromisoformat(value)
        except ValueError:
            return False
    return all(
        row.get(c) is None
        or (isinstance(row[c], (int, float)) and not isinstance(row[c], bool))
        for c in numeric
    )


def _status_from_exit_code(row: Dict[str, Any]) -> bool:
    """Finish a job whose exit code was recorded but whose status wasn't."""
    code = row.get("exit_code")
    if row.get("status") in FINISHED or not isinstance(code, int):
        return False
    row["status"] = "completed" if code == 0 else "failed"
    return True


def _worth_keeping(row: Dict[str, Any]) -> bool:
    return row.get("status") in FINISHED or row.get("mem_used_mb") is not None


def recovered_from(snapshot: Optional[Path] = None) -> Optional[str]:
    """What the snapshot says was recovered into it, if anything."""
    from .consolidate import open_snapshot  # pylint: disable=import-outside-toplevel

    try:
        conn = open_snapshot(snapshot)
    except (OSError, sqlite3.Error):
        return None
    try:
        row = conn.execute(
            "SELECT value FROM snapshot_meta WHERE key=?", (RECOVERED_FROM,)
        ).fetchone()
        return row[0] if row else None
    except sqlite3.Error:
        return None
    finally:
        conn.close()
