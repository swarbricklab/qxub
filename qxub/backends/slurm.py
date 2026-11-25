"""Slurm backend: ``sbatch``, ``squeue``, ``sacct``, ``scancel`` and ``srun``.

A platform with ``type: slurm`` selects it (see
``docs/platforms/slurm/slurm.yaml``). qxub's queues are Slurm partitions,
and the account is ``--account``.

- **Submitting:** ``sbatch --parsable``, with the job's variables in
  sbatch's own environment and ``--export=ALL``, Slurm's default, so values
  with commas reach the job unchanged. Requests become ``--cpus-per-task``,
  ``--mem``, ``--time``, ``--gres=gpu:N`` and ``--tmp``; any other ``-l
  name=value`` is passed as ``--name=value`` (``-l constraint=skylake``).
  The job log is ``--output``.
- **Interactive sessions:** ``srun --pty … bash script``.
- **Live jobs:** ``squeue``.
- **Finished jobs:** ``sacct``. Slurm says why a job ended (``TIMEOUT``,
  ``CANCELLED``, ``OUT_OF_MEMORY``), so ``ended_by`` comes from the state,
  not from the exit status.

The parsers are plain functions, tested on captured output.
"""

import getpass
import logging
import math
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import replace
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Tuple

from ..resources.spec import format_duration, parse_size
from .base import (
    ENDED_BY_CANCEL,
    ENDED_BY_MEMORY,
    ENDED_BY_WALLTIME,
    JobEnd,
    JobInfo,
    JobSpec,
    JobState,
    JobUsage,
    LiveJob,
    SchedulerBackend,
    SchedulerQueryError,
    SubmitError,
)

logger = logging.getLogger(__name__)

# A Slurm job ID: "12345", an array task "12345_7", a pending array
# "12345_[1-10]", or a heterogeneous component "12345+1".
JOB_ID_RE = re.compile(r"^\d+(_\d+|_\[[^\]]*\])?(\+\d+)?$")

# Job IDs per ``sacct`` call.
SACCT_CHUNK = 100

# What ``sacct`` reports, in this order (``--parsable2``: "|"-separated).
SACCT_FIELDS = [
    "JobID",
    "JobName",
    "Partition",
    "State",
    "ExitCode",
    "Submit",
    "Start",
    "End",
    "NodeList",
    "ReqTRES",
    "AllocTRES",
    "TimelimitRaw",
    "ElapsedRaw",
    "TotalCPU",
    "MaxRSS",
]

# ``squeue`` columns: job ID, partition, compact state, whole name.
SQUEUE_FORMAT = "%i|%P|%t|%j"

# squeue's compact states
_LIVE_STATES = {
    "PD": JobState.QUEUED,  # pending
    "RQ": JobState.QUEUED,  # requeued
    "RF": JobState.QUEUED,  # requeued after a node failure
    "RS": JobState.QUEUED,  # being resized
    "CF": JobState.RUNNING,  # configuring: nodes allocated, booting
    "R": JobState.RUNNING,
    "CG": JobState.RUNNING,  # completing
    "SO": JobState.RUNNING,  # staging out
    "S": JobState.HELD,  # suspended
    "ST": JobState.HELD,  # stopped
    "RH": JobState.HELD,  # requeued and held
    "RD": JobState.HELD,  # held for a reservation being deleted
}

# sacct's states (the first word: "CANCELLED by 1234")
_STATES = {
    "PENDING": JobState.QUEUED,
    "REQUEUED": JobState.QUEUED,
    "RESIZING": JobState.QUEUED,
    "RUNNING": JobState.RUNNING,
    "COMPLETING": JobState.RUNNING,
    "STAGE_OUT": JobState.RUNNING,
    "SUSPENDED": JobState.HELD,
    "STOPPED": JobState.HELD,
    "REQUEUE_HOLD": JobState.HELD,
}
_ENDED_BY = {
    "TIMEOUT": ENDED_BY_WALLTIME,
    "DEADLINE": ENDED_BY_WALLTIME,
    "CANCELLED": ENDED_BY_CANCEL,
    "OUT_OF_MEMORY": ENDED_BY_MEMORY,
}
# The exit status to report when Slurm stopped a job and sacct's ExitCode
# says nothing (0:0, e.g. cancelled while queued): what the job script's
# TERM trap records when it's stopped (143), or a SIGKILL for memory.
_STOPPED_EXIT = {
    ENDED_BY_WALLTIME: 143,
    ENDED_BY_CANCEL: 143,
    ENDED_BY_MEMORY: 137,
}

# PBS mail letters (qxub's --email-opts) as Slurm --mail-type values
_MAIL_TYPES = {"a": "FAIL", "b": "BEGIN", "e": "END", "n": "NONE"}


def _current_user() -> str:
    return os.environ.get("USER") or getpass.getuser()


def _megabytes(num_bytes: int) -> str:
    """A size for ``--mem``/``--tmp``, whose plain unit is MB."""
    return f"{math.ceil(num_bytes / 1024**2)}M"


def mail_types(letters: Optional[str]) -> str:
    """qxub's PBS-style mail letters (``ae``) as ``--mail-type`` (``FAIL,END``)."""
    types = [_MAIL_TYPES[c] for c in (letters or "ae") if c in _MAIL_TYPES]
    return ",".join(dict.fromkeys(types)) or "FAIL,END"


def slurm_options(spec: JobSpec, interactive: bool = False) -> List[str]:
    """The ``sbatch``/``srun`` options for ``spec``.

    An interactive session gets no job log, mail or site options: the
    session script carries what it needs.
    """
    opts = [f"--job-name={spec.name}"]
    if spec.queue:
        opts.append(f"--partition={spec.queue}")
    if spec.account:
        opts.append(f"--account={spec.account}")
    res = spec.resources
    if res.cpus is not None:
        opts.append(f"--cpus-per-task={res.cpus}")
    if res.memory is not None:
        opts.append(f"--mem={_megabytes(res.memory)}")
    if res.walltime is not None:
        opts.append(f"--time={format_duration(res.walltime)}")
    if res.gpus:
        opts.append(f"--gres=gpu:{res.gpus}")
    if res.disk is not None:
        opts.append(f"--tmp={_megabytes(res.disk)}")
    opts += [f"--{name}={value}" for name, value in res.native.items()]
    if interactive:
        return opts
    if spec.job_log:
        opts.append(f"--output={spec.job_log}")
    if spec.mail_to:
        opts += [
            f"--mail-user={spec.mail_to}",
            f"--mail-type={mail_types(spec.mail_events)}",
        ]
    for arg in spec.native_args:
        opts += shlex.split(arg)
    return opts


def job_variables(spec: JobSpec) -> Dict[str, str]:
    """The variables the job gets: the settings file, or ``spec.env``."""
    if spec.env_file:
        return {"QXUB_JOB_FILE": str(spec.env_file)}
    return {k: str(v) for k, v in spec.env.items()}


def sbatch_argv(spec: JobSpec) -> List[str]:
    """The ``sbatch`` command line for ``spec`` (its variables go in
    sbatch's environment, :func:`job_variables`)."""
    return [
        "sbatch",
        "--parsable",
        "--export=ALL",
        *slurm_options(spec),
        str(spec.script),
    ]


def srun_interactive_argv(spec: JobSpec) -> List[str]:
    """The ``srun --pty`` command that runs ``spec.script`` interactively."""
    spec = replace(
        spec, env={}, env_file=None, job_log=None, mail_to=None, native_args=[]
    )
    return [
        "srun",
        "--pty",
        *slurm_options(spec, interactive=True),
        "bash",
        str(spec.script),
    ]


class SlurmBackend(SchedulerBackend):
    """Slurm."""

    name = "slurm"
    job_id_variable = "SLURM_JOB_ID"

    def available(self) -> bool:
        return shutil.which("squeue") is not None

    # -- submission -------------------------------------------------------

    def submit_command(self, spec: JobSpec) -> List[str]:
        variables = [f"{k}={v}" for k, v in job_variables(spec).items()]
        return (
            ["env", *variables, *sbatch_argv(spec)] if variables else sbatch_argv(spec)
        )

    def submit(self, spec: JobSpec) -> str:
        argv = sbatch_argv(spec)
        logger.info("Submitting: %s", shlex.join(argv))
        try:
            result = subprocess.run(
                argv,
                capture_output=True,
                text=True,
                check=False,
                env={**os.environ, **job_variables(spec)},
            )
        except OSError as exc:
            raise SubmitError(f"cannot run sbatch: {exc}", exit_code=1) from exc
        if result.returncode != 0:
            raise SubmitError(
                f"sbatch failed (exit {result.returncode}): {result.stderr.strip()}",
                exit_code=result.returncode,
                stderr=result.stderr,
            )
        return parse_sbatch_output(result.stdout)

    def interactive_command(self, spec: JobSpec) -> List[str]:
        return srun_interactive_argv(spec)

    # -- queries ----------------------------------------------------------

    def list_jobs(
        self, user: Optional[str] = None, full_names: bool = False
    ) -> List[LiveJob]:
        # squeue never cuts a name printed with a bare %j
        argv = [
            "squeue",
            "--noheader",
            f"--user={user or _current_user()}",
            f"--format={SQUEUE_FORMAT}",
        ]
        try:
            result = subprocess.run(argv, capture_output=True, text=True, check=False)
        except OSError as exc:
            raise SchedulerQueryError(f"cannot run squeue: {exc}") from exc
        if result.returncode != 0:
            raise SchedulerQueryError(
                f"{' '.join(argv)} failed (exit {result.returncode}): "
                f"{result.stderr.strip()}"
            )
        return parse_squeue_output(result.stdout)

    def live_state(self, state: str) -> JobState:
        return _LIVE_STATES.get(state, JobState.UNKNOWN)

    def job_details(self, job_ids: Iterable[str]) -> Dict[str, Dict[str, Any]]:
        ids = list(dict.fromkeys(job_ids))
        found: Dict[str, Dict[str, Any]] = {}
        for start in range(0, len(ids), SACCT_CHUNK):
            chunk = ids[start : start + SACCT_CHUNK]
            argv = [
                "sacct",
                "--noheader",
                "--parsable2",
                f"--jobs={','.join(chunk)}",
                f"--format={','.join(SACCT_FIELDS)}",
            ]
            try:
                result = subprocess.run(
                    argv, capture_output=True, text=True, timeout=120, check=False
                )
            except (OSError, subprocess.SubprocessError) as exc:
                logger.debug("sacct failed: %s", exc)
                continue
            if result.returncode != 0:
                logger.debug("sacct failed: %s", result.stderr.strip())
                continue
            found.update(parse_sacct_output(result.stdout))
        return found

    def job_info(self, job_ids: Iterable[str]) -> Dict[str, JobInfo]:
        return {
            job_id: job_info_from_sacct(job_id, details)
            for job_id, details in self.job_details(job_ids).items()
        }

    def job_end(self, job_id: str, job_log=None) -> Optional[JobEnd]:
        info = self.find_info(job_id)
        if info is None or not info.is_finished or info.exit_status is None:
            return None
        return JobEnd(info.exit_status, info.ended_by)

    def looks_like_job_id(self, value: str) -> bool:
        return bool(JOB_ID_RE.match(value or ""))

    # -- control ----------------------------------------------------------

    def cancel(self, job_id: str) -> Tuple[bool, str]:
        try:
            result = subprocess.run(
                ["scancel", job_id], capture_output=True, text=True, check=False
            )
        except OSError as exc:
            return False, f"cannot run scancel: {exc}"
        if result.returncode == 0:
            return True, ""
        return False, result.stderr.strip() or f"scancel exited {result.returncode}"


# -- parsers -------------------------------------------------------------------


def parse_sbatch_output(output: str) -> str:
    """The job ID from ``sbatch --parsable`` (``12345`` or ``12345;cluster``)."""
    line = output.strip().splitlines()[-1] if output.strip() else ""
    job_id = line.split(";", 1)[0].strip()
    if not JOB_ID_RE.match(job_id):
        raise SubmitError(f"sbatch printed no job ID: {output.strip()!r}")
    return job_id


def parse_squeue_output(output: str) -> List[LiveJob]:
    """Parse ``squeue --noheader --format='%i|%P|%t|%j'`` output."""
    jobs = []
    for line in output.splitlines():
        parts = line.split("|", 3)
        if len(parts) == 4 and JOB_ID_RE.match(parts[0].strip()):
            job_id, partition, state, name = (p.strip() for p in parts)
            jobs.append(LiveJob(job_id, partition, state, name))
    return jobs


def parse_sacct_output(output: str) -> Dict[str, Dict[str, Any]]:
    """Parse ``sacct --parsable2`` output (:data:`SACCT_FIELDS`), by job ID.

    Each job is its own line's fields, with its steps (``12345.batch``,
    ``12345.extern``, ``12345.0``) under ``"steps"``: their ``MaxRSS``
    is the job's memory use.
    """
    jobs: Dict[str, Dict[str, Any]] = {}
    for line in output.splitlines():
        values = line.split("|")
        if len(values) != len(SACCT_FIELDS):
            continue
        row = dict(zip(SACCT_FIELDS, values))
        job_id, _, step = row["JobID"].partition(".")
        if step:
            if job_id in jobs:
                jobs[job_id]["steps"][step] = row
        elif JOB_ID_RE.match(job_id):
            jobs[job_id] = {**row, "steps": {}}
    return jobs


def _when(value: str) -> Optional[datetime]:
    try:
        return datetime.fromisoformat(value)
    except (TypeError, ValueError):
        return None  # "Unknown", "None"


def _tres(value: str) -> Dict[str, str]:
    """``cpu=4,mem=16G,gres/gpu=1`` as a dict."""
    found = {}
    for item in (value or "").split(","):
        name, sep, amount = item.partition("=")
        if sep:
            found[name.strip()] = amount.strip()
    return found


def _int(value) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _bytes(value) -> Optional[int]:
    """A Slurm size (``16G``, ``4000M``, ``123456K``; plain is bytes)."""
    try:
        return parse_size(value, bare_unit="B") if value else None
    except ValueError:
        return None


def parse_slurm_duration(value: str) -> Optional[int]:
    """Seconds from a Slurm duration: ``[D-][HH:]MM:SS[.mmm]``."""
    if not value:
        return None
    days, _, rest = value.rpartition("-")
    try:
        parts = [float(p) for p in rest.split(":")]
        day_count = int(days) if days else 0
    except ValueError:
        return None
    if not 1 <= len(parts) <= 3:
        return None
    seconds = 0.0
    for part in parts:
        seconds = seconds * 60 + part
    return int(day_count * 86400 + seconds)


def _exit_status(state: str, exit_code: str, ended_by: Optional[str]) -> Optional[int]:
    """The job's exit status from sacct's ``ExitCode`` (``code:signal``):
    128 + the signal if one killed it, as a shell reports it."""
    code, _, signal = (exit_code or "").partition(":")
    code_n, signal_n = _int(code), _int(signal) or 0
    if code_n is None:
        return None
    # Not a signal number: OUT_OF_MEMORY jobs show 0:125
    status = 128 + signal_n if 0 < signal_n < 65 else code_n
    if status == 0 and ended_by:
        status = _STOPPED_EXIT[ended_by]
    elif status == 0 and state != "COMPLETED":
        status = 1  # FAILED, NODE_FAIL, BOOT_FAIL, PREEMPTED, ... with no code
    return status


def job_info_from_sacct(job_id: str, details: Dict[str, Any]) -> JobInfo:
    """A :class:`JobInfo` from one job of :func:`parse_sacct_output`."""
    state_word = (details.get("State") or "").split(" ")[0].rstrip("+")
    state = _STATES.get(
        state_word, JobState.FINISHED if state_word else JobState.UNKNOWN
    )
    finished = state is JobState.FINISHED
    ended_by = _ENDED_BY.get(state_word) if finished else None
    requested = _tres(details.get("ReqTRES"))
    allocated = _tres(details.get("AllocTRES"))
    rss = [_bytes(s.get("MaxRSS")) for s in details.get("steps", {}).values()]
    rss = [r for r in rss if r is not None]
    limit = _int(details.get("TimelimitRaw"))  # minutes
    host = (details.get("NodeList") or "").split(",")[0]
    gpus = requested.get("gres/gpu")
    return JobInfo(
        job_id=job_id,
        state=state,
        exit_status=(
            _exit_status(state_word, details.get("ExitCode"), ended_by)
            if finished
            else None
        ),
        ended_by=ended_by,
        name=details.get("JobName") or None,
        queue=details.get("Partition") or None,
        host=host if host and host != "None assigned" else None,
        submitted=_when(details.get("Submit")),
        started=_when(details.get("Start")),
        finished=_when(details.get("End")) if finished else None,
        usage=JobUsage(
            cpus_requested=_int(requested.get("cpu")),
            cpus_used=_int(allocated.get("cpu")),
            cpu_time=parse_slurm_duration(details.get("TotalCPU")),
            memory_requested=_bytes(requested.get("mem")),
            memory_used=max(rss) if rss else None,
            walltime_requested=limit * 60 if limit is not None else None,
            walltime_used=_int(details.get("ElapsedRaw")),
            gpus_requested=_int(gpus) if gpus is not None else None,
        ),
        native=details,
    )


__all__ = ["SlurmBackend", "parse_sacct_output", "parse_squeue_output"]
