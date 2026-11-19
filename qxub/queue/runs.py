"""One file per run: the job's record, written without SQLite or locks.

Phase 2 of ``docs/dev/sqlite-robustness.md``. On Gadi, ``/g/data`` and
``/scratch`` are mounted with ``localflock``, so SQLite processes on
different nodes don't see each other's locks, and jobs writing the shared
database from compute nodes corrupted it. Instead, each run has one small
text file:

    <status dir>/runs/<run_id>.run

- The submitting process writes it before ``qsub`` (``state=submitted``).
- The job adds to it at start (``state=running``) and at exit
  (``state=finished``, ``exit_code``) with plain shell; see the job
  templates. Its path is passed to the job as ``QXUB_RUN_FILE``.
- Every write replaces the whole file through a temporary file and
  ``rename``, which is atomic across nodes on Lustre. The writers take turns
  (submission happens before the job can start), so no locks are needed.

The format is one ``key=value`` per line so the job script can extend it
without a JSON parser. A key may appear more than once; the last one wins.
Values are single-line; lists and dicts are stored as compact JSON.
"""

import json
import logging
import os
import uuid
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

SUFFIX = ".run"
_warned = False


def runs_dir() -> Path:
    """Directory holding run files; visible from compute nodes."""
    from ..config.site import default_status_dir

    return default_status_dir() / "runs"


def run_file(run_id: str) -> Path:
    return runs_dir() / f"{run_id}{SUFFIX}"


def _value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        text = json.dumps(value, separators=(",", ":"), default=str)
    else:
        text = str(value)
    return text.replace("\r", " ").replace("\n", " ")


def format_fields(fields: Dict[str, Any]) -> str:
    return "".join(f"{key}={_value(value)}\n" for key, value in fields.items())


def write_submitted(run_id: str, fields: Dict[str, Any]) -> Optional[Path]:
    """Write the run file at submission. Never raises; warns once if it fails.

    Returns the file's path, or None if it could not be written. The job
    still creates the file itself at start if this failed.
    """
    path = run_file(run_id)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        try:
            with open(tmp, "w", encoding="utf-8") as f:
                f.write(format_fields({"run_id": run_id, **fields}))
            os.replace(tmp, path)
        except BaseException:
            try:
                tmp.unlink()
            except OSError:
                pass
            raise
        return path
    except OSError as exc:
        global _warned  # noqa: PLW0603
        level = logging.DEBUG if _warned else logging.WARNING
        _warned = True
        logger.log(level, "Could not write run file %s: %s", path, exc)
        return None


def read_run(path) -> Optional[Dict[str, str]]:
    """Parse a run file; None if missing or unreadable. Last value wins."""
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError:
        return None
    record: Dict[str, str] = {}
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if sep and key:
            record[key] = value
    return record


def exit_code(record: Dict[str, str]) -> Optional[int]:
    try:
        return int(record.get("exit_code", ""))
    except ValueError:
        return None


def submission_vars(run_id: str, path: Optional[Path]) -> str:
    """``qsub -v`` additions that tell the job where its run file is."""
    path = path or run_file(run_id)
    return f",QXUB_RUN_ID={run_id},QXUB_RUN_FILE={path}"
