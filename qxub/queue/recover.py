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
3. Only finished jobs, or jobs with resource usage, are added: jobs still
   marked submitted or running there are stale or will arrive from run files.
   A job the snapshot already has is never overwritten.

Running it again adds nothing new, so it is safe to repeat.
"""

import logging
import shutil
import sqlite3
import subprocess
import tempfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .consolidate import _round, stats_db_path

logger = logging.getLogger(__name__)

FINISHED = ("completed", "failed", "cancelled")
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
            f"{self.unfinished} unfinished and left out"
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
        wanted = [c for c in columns if c in _snapshot_columns(conn)]
        if "job_id" not in wanted:
            raise sqlite3.DatabaseError(f"{source} has no job_id column")
        index = {c: columns.index(c) for c in wanted}
        names = ", ".join(wanted)
        marks = ", ".join("?" * len(wanted))
        for row in rows:
            job_id = row[columns.index("job_id")]
            if not job_id:
                continue
            if not _worth_keeping(dict(zip(columns, row))):
                result.unfinished += 1
                continue
            cur = conn.execute(
                f"INSERT OR IGNORE INTO job_resources ({names}) VALUES ({marks})",
                [row[index[c]] for c in wanted],
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


def _sqlite3_recover(copy: Path, rebuilt: Path, first_error: Exception) -> None:
    exe = shutil.which("sqlite3")
    if not exe:
        raise RuntimeError(
            f"{first_error}, and the sqlite3 command isn't on PATH to run"
            " .recover (try `module load sqlite`)"
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


def _snapshot_columns(conn) -> set:
    return {r[1] for r in conn.execute("PRAGMA table_info(job_resources)")}


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
