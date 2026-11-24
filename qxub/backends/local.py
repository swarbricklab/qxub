"""Local backend: each job is a background process on this machine.

No scheduler is involved. ``submit`` starts the job script detached from
the terminal, with its output going to the job log, and records the job in
``<status dir>/local/<job id>/``:

- ``pid``: the process group running the job;
- ``meta.json``: name, queue, job log, submission time;
- ``exit_status``: written when the job ends.

When the job ends, its exit status is also appended to the job log as
``Exit Status: N``, the line ``qxub exec`` looks for in a PBS job log.
Resource requests are not enforced. A platform with ``type: local`` selects
this backend (see ``docs/platforms/local/local.yaml``).
"""

import json
import logging
import os
import signal
import subprocess
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from .base import JobSpec, LiveJob, SchedulerBackend, SubmitError

logger = logging.getLogger(__name__)

JOB_ID_VARIABLE = "QXUB_LOCAL_JOBID"
JOB_ID_PREFIX = "local-"
QUEUE = "local"

# Runs the job script, then records its exit status. Arguments: the job
# script, the job's directory, the job log. `cancel` signals the whole
# process group: this shell catches the signal (a no-op handler, not an
# ignored signal, which the job script could not trap) and so lives to
# record the job's exit.
_RUNNER = """
trap ':' TERM INT HUP
bash "$1"
code=$?
printf '\\nExit Status: %s\\n' "$code" >> "$3"
printf '%s\\n' "$code" > "$2/.exit_status.tmp" && mv -f "$2/.exit_status.tmp" "$2/exit_status"
exit "$code"
"""


def jobs_dir() -> Path:
    """Where local jobs are recorded."""
    from ..config.site import (  # pylint: disable=import-outside-toplevel
        default_status_dir,
    )

    return default_status_dir() / "local"


def _job_env(spec: JobSpec, job_id: str) -> Dict[str, str]:
    """The job's variables: the settings file, or ``spec.env`` itself."""
    env = {JOB_ID_VARIABLE: job_id}
    if spec.env_file:
        env["QXUB_JOB_FILE"] = spec.env_file
    else:
        env.update({k: str(v) for k, v in spec.env.items()})
    return env


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class LocalBackend(SchedulerBackend):
    """Jobs as background processes on this machine."""

    name = "local"
    job_id_variable = JOB_ID_VARIABLE

    # -- submission -------------------------------------------------------

    def submit_command(self, spec: JobSpec) -> List[str]:
        env = _job_env(spec, f"{JOB_ID_PREFIX}<id>")
        return ["env", *(f"{k}={v}" for k, v in env.items()), "bash", spec.script]

    def submit(self, spec: JobSpec) -> str:
        job_id = f"{JOB_ID_PREFIX}{uuid.uuid4().hex[:12]}"
        job_dir = jobs_dir() / job_id
        job_log = Path(spec.job_log) if spec.job_log else job_dir / "job.log"
        try:
            job_dir.mkdir(parents=True)
            job_log.parent.mkdir(parents=True, exist_ok=True)
            (job_dir / "meta.json").write_text(
                json.dumps(
                    {
                        "name": spec.name,
                        "queue": QUEUE,
                        "job_log": str(job_log),
                        "submitted": datetime.now().isoformat(timespec="seconds"),
                    }
                )
            )
            with open(job_log, "ab") as log:
                process = subprocess.Popen(  # pylint: disable=consider-using-with
                    [
                        "bash",
                        "-c",
                        _RUNNER,
                        "qxub-local",
                        spec.script,
                        str(job_dir),
                        str(job_log),
                    ],
                    stdin=subprocess.DEVNULL,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    env={**os.environ, **_job_env(spec, job_id)},
                    start_new_session=True,  # its own process group, off the terminal
                )
            (job_dir / "pid").write_text(f"{process.pid}\n")
        except OSError as exc:
            raise SubmitError(f"cannot start the job: {exc}") from exc
        logger.info("Started local job %s (pid %s)", job_id, process.pid)
        return job_id

    def interactive_command(self, spec: JobSpec) -> List[str]:
        job_id = f"{JOB_ID_PREFIX}{uuid.uuid4().hex[:12]}"
        return ["env", f"{JOB_ID_VARIABLE}={job_id}", "bash", spec.script]

    # -- queries ----------------------------------------------------------

    def _read(self, job_id: str) -> Optional[Dict[str, Any]]:
        job_dir = jobs_dir() / job_id
        try:
            meta = json.loads((job_dir / "meta.json").read_text())
            pid = int((job_dir / "pid").read_text().strip())
        except (OSError, ValueError):
            return None
        try:
            exit_status: Optional[int] = int(
                (job_dir / "exit_status").read_text().strip()
            )
        except (OSError, ValueError):
            exit_status = None
        if exit_status is None and not _alive(pid):
            exit_status = -1  # ended without recording a status (killed)
        return {**meta, "pid": pid, "exit_status": exit_status}

    def list_jobs(
        self, user: Optional[str] = None, full_names: bool = False
    ) -> List[LiveJob]:
        jobs = []
        try:
            entries = sorted(jobs_dir().iterdir())
        except OSError:
            return []
        for entry in entries:
            job = self._read(entry.name)
            if job and job["exit_status"] is None:
                jobs.append(LiveJob(entry.name, job["queue"], "R", job["name"]))
        return jobs

    def job_details(self, job_ids: Iterable[str]) -> Dict[str, Dict[str, Any]]:
        found = {}
        for job_id in dict.fromkeys(job_ids):
            job = self._read(job_id)
            if not job:
                continue
            finished = job["exit_status"] is not None
            found[job_id] = {
                "exit_status": job["exit_status"],
                "metadata": {
                    "job_state": "F" if finished else "R",
                    "job_name": job["name"],
                    "queue": job["queue"],
                },
                "timing": {"submitted": job["submitted"]},
                "resources_requested": {},
                "resources_used": {},
            }
        return found

    def is_finished(self, details: Optional[Dict[str, Any]]) -> bool:
        return ((details or {}).get("metadata") or {}).get("job_state") == "F"

    def looks_like_job_id(self, value: str) -> bool:
        return (value or "").startswith(JOB_ID_PREFIX)

    # -- control ----------------------------------------------------------

    def cancel(self, job_id: str) -> Tuple[bool, str]:
        job = self._read(job_id)
        if not job:
            return False, f"no such local job: {job_id}"
        if job["exit_status"] is not None:
            return False, f"job {job_id} has already finished"
        try:
            # The job script's TERM trap records the exit (143)
            os.killpg(job["pid"], signal.SIGTERM)
        except ProcessLookupError:
            return False, f"job {job_id} has already finished"
        except OSError as exc:
            return False, f"cannot stop job {job_id}: {exc}"
        return True, ""


__all__ = ["LocalBackend", "jobs_dir", "JOB_ID_VARIABLE"]
