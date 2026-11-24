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
import re
import shlex
import uuid
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

SUFFIX = ".run"
# The job's settings, next to its run file (see attach_settings).
SETTINGS_SUFFIX = ".env"
_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_warned = False


def runs_dir() -> Path:
    """Directory holding run files; visible from compute nodes."""
    from ..config.site import default_status_dir

    return default_status_dir() / "runs"


def run_file(run_id: str) -> Path:
    return run_file_in(runs_dir(), run_id)


def run_file_in(directory, run_id: str) -> Path:
    return Path(directory) / f"{run_id}{SUFFIX}"


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


def discard(path: Optional[Path]) -> None:
    """Remove a run or settings file whose job was never submitted."""
    if path is None:
        return
    try:
        Path(path).unlink()
    except OSError:
        pass


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


def settings_file(run_id: str) -> Path:
    return runs_dir() / f"{run_id}{SETTINGS_SUFFIX}"


def format_settings(env: Dict[str, Any]) -> str:
    """``env`` as a bash file the job script sources: ``export K='v'``."""
    lines = ["# qxub job settings, read by the job script (QXUB_JOB_FILE)"]
    for name, value in env.items():
        if not _NAME_RE.match(name):
            raise ValueError(f"not a variable name: {name!r}")
        lines.append(f"export {name}={shlex.quote(str(value))}")
    return "\n".join(lines) + "\n"


def reads_settings_file(script) -> bool:
    """Whether the job script loads ``$QXUB_JOB_FILE``.

    qxub's own job scripts do. A custom ``--template`` written for older
    qxub expects its variables from the scheduler (``qsub -v``) instead.
    """
    try:
        return "QXUB_JOB_FILE" in Path(script).read_text(errors="replace")
    except OSError:
        return False


def attach_settings(spec, run_id: str, write: bool = True):
    """``spec`` with its variables in a settings file instead of on the
    scheduler's command line.

    ``qsub -v`` splits on commas and limits the list's length, and each
    scheduler passes variables differently. The job's variables go to
    ``<runs dir>/<run_id>.env`` instead (readable by the user only), and
    the backend passes the job only its path, as ``QXUB_JOB_FILE``.

    With ``write=False`` (a dry run) the file is named but not written. If
    it can't be written, ``spec`` is returned unchanged and the variables
    go on the command line as before. So are they for a job script that
    doesn't read ``QXUB_JOB_FILE``, such as a ``--template`` written for
    older qxub.
    """
    from dataclasses import replace  # pylint: disable=import-outside-toplevel

    if not spec.env or not reads_settings_file(spec.script):
        return spec
    path = settings_file(run_id)
    if write:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
            fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as f:
                    f.write(format_settings(spec.env))
                os.replace(tmp, path)
            except BaseException:
                discard(tmp)
                raise
        except (OSError, ValueError) as exc:
            logger.warning(
                "Could not write job settings %s (%s); passing them on the "
                "command line",
                path,
                exc,
            )
            return spec
    return replace(spec, env_file=str(path))


def job_env(run_id: str, path: Optional[Path]) -> Dict[str, str]:
    """Job variables that tell the job its run ID and where to record status.

    The job writes its status files to ``QXUB_STATUS_DIR``, the same
    directory ``qxub monitor`` and ``qxub status`` read, named after
    ``QXUB_RUN_ID``.
    """
    from ..config.site import default_status_dir

    return {
        "QXUB_RUN_ID": run_id,
        "QXUB_RUN_FILE": str(path or run_file(run_id)),
        "QXUB_STATUS_DIR": str(default_status_dir()),
    }
