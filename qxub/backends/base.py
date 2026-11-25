"""The interface qxub uses to ask a scheduler about jobs and to cancel them.

Each scheduler (PBS Pro, Slurm, a local runner, and later others) gets
one :class:`SchedulerBackend`. The platform definition's ``type`` picks
which one (see :func:`qxub.backends.get_backend`).

A job to submit, or an interactive session to start, is described by a
:class:`JobSpec`; the backend turns it into its own command line, including the :class:`ResourceSpec` (PBS writes
``-l ncpus=…``).
"""

import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, fields
from datetime import datetime
from enum import Enum
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, Tuple

from ..resources.spec import ResourceSpec


class LiveJob(NamedTuple):
    """A job the scheduler still lists (queued, running, held, ...).

    ``state`` is the scheduler's own code (PBS: ``Q``, ``R``, ``H``, ...).
    ``name`` may be truncated by the scheduler's listing.
    """

    job_id: str
    queue: str
    state: str
    name: str


class JobState(str, Enum):
    """Where a job is, in terms every scheduler can map to."""

    QUEUED = "queued"  # waiting to start (PBS Q, W, T)
    RUNNING = "running"  # started, including ending (PBS R, B, E)
    HELD = "held"  # won't start until released (PBS H, S, U)
    FINISHED = "finished"  # ended, however it ended (PBS F, X)
    UNKNOWN = "unknown"


#: :attr:`JobInfo.ended_by` values, when the scheduler can tell.
ENDED_BY_WALLTIME = "walltime"  # killed for running past its walltime
ENDED_BY_CANCEL = "cancelled"  # cancelled (qxub cancel, qdel, ...)
ENDED_BY_MEMORY = "memory"  # killed for using more memory than it asked for

_ENDED_BY_TEXT = {
    ENDED_BY_WALLTIME: "stopped at its walltime limit",
    ENDED_BY_CANCEL: "cancelled",
    ENDED_BY_MEMORY: "stopped for exceeding its memory limit",
}


def describe_ended_by(ended_by: str) -> str:
    """How a job ended, to follow "Job was" / "Job": ``stopped at its
    walltime limit``, ``cancelled``."""
    return _ENDED_BY_TEXT.get(ended_by, f"ended by {ended_by}")


@dataclass
class JobUsage:
    """What a job asked for and used. Sizes in bytes, times in seconds.

    Anything the scheduler doesn't report is ``None``.
    """

    cpus_requested: Optional[int] = None
    cpus_used: Optional[int] = None
    #: CPU time summed over the job's CPUs.
    cpu_time: Optional[int] = None
    #: The scheduler's own CPU percentage (PBS cpupercent: 100 per busy CPU).
    cpu_percent: Optional[float] = None
    memory_requested: Optional[int] = None
    memory_used: Optional[int] = None
    walltime_requested: Optional[int] = None
    walltime_used: Optional[int] = None
    #: Node-local disk (PBS jobfs).
    disk_requested: Optional[int] = None
    disk_used: Optional[int] = None
    gpus_requested: Optional[int] = None

    def has_used(self) -> bool:
        """Whether any usage (not only requests) was reported."""
        return any(
            getattr(self, f.name) is not None
            for f in fields(self)
            if f.name.endswith("_used") or f.name in ("cpu_time", "cpu_percent")
        )


@dataclass
class JobInfo:
    """What a scheduler knows about one job, in scheduler-neutral terms.

    ``exit_status`` is the scheduler's own value (PBS: negative if PBS
    ended the job, 256 + N if it was killed by signal N). ``ended_by`` says
    why the scheduler ended it, when it can tell (``ENDED_BY_*``).
    ``native`` is the scheduler's own record (PBS: ``parse_qstat_fx_output``).
    """

    job_id: str
    state: JobState
    exit_status: Optional[int] = None
    ended_by: Optional[str] = None
    name: Optional[str] = None
    queue: Optional[str] = None
    host: Optional[str] = None
    submitted: Optional[datetime] = None
    started: Optional[datetime] = None
    finished: Optional[datetime] = None
    usage: JobUsage = field(default_factory=JobUsage)
    native: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_finished(self) -> bool:
        return self.state is JobState.FINISHED

    @property
    def queue_wait(self) -> Optional[int]:
        """Seconds from submission to start."""
        if self.submitted and self.started:
            return int((self.started - self.submitted).total_seconds())
        return None

    @property
    def run_time(self) -> Optional[int]:
        """Seconds from start to end."""
        if self.started and self.finished:
            return int((self.finished - self.started).total_seconds())
        return None


class JobEnd(NamedTuple):
    """How a job ended, as its scheduler reports it."""

    exit_status: int
    ended_by: Optional[str] = None


_JOB_LOG_EXIT_RE = re.compile(r"Exit Status:\s*(-?\d+)")


def read_job_log_exit(path) -> Optional[int]:
    """The exit status in a job log's ``Exit Status:`` line, if it has one.

    Gadi's PBS appends that line when a job ends, however it ended; the
    local backend writes it too. One read of one file.
    """
    try:
        with open(str(path), "r", encoding="utf-8") as f:
            found = _JOB_LOG_EXIT_RE.search(f.read())
    except (OSError, ValueError):
        return None
    return int(found.group(1)) if found else None


# The state letters `qxub cancel all --status` takes (PBS's)
_LIVE_STATES = {"Q": JobState.QUEUED, "R": JobState.RUNNING, "H": JobState.HELD}


class SchedulerQueryError(Exception):
    """The scheduler couldn't be asked which jobs it still has."""


class SubmitError(Exception):
    """The scheduler refused a job, or couldn't be asked to take it.

    Attributes:
        exit_code: What the CLI should exit with.
        stderr:    The scheduler's own message, if any.
    """

    def __init__(self, message: str, exit_code: int = 1, stderr: str = ""):
        super().__init__(message)
        self.exit_code = exit_code
        self.stderr = stderr


@dataclass
class JobSpec:
    """A job to submit, in terms that don't depend on the scheduler.

    ``env`` holds the variables the job script reads (``cmd_b64``, ``out``,
    ``QXUB_RUN_ID``, ...), in the order they should be passed.
    """

    script: str
    name: str
    queue: Optional[str] = None
    #: The account or project the job is charged to (PBS ``-P``).
    account: Optional[str] = None
    resources: ResourceSpec = field(default_factory=ResourceSpec)
    #: Where the scheduler writes its own log of the job (PBS ``-o``).
    job_log: Optional[str] = None
    mail_to: Optional[str] = None
    #: When to send mail, in the scheduler's own letters (PBS ``-m ae``).
    mail_events: Optional[str] = None
    env: Dict[str, str] = field(default_factory=dict)
    #: A file holding ``env`` (``qxub.queue.runs.attach_settings``). When
    #: set, the backend passes the job only its path, as ``QXUB_JOB_FILE``.
    env_file: Optional[str] = None
    #: Extra scheduler arguments the platform adds to every job.
    native_args: List[str] = field(default_factory=list)


class SchedulerBackend(ABC):
    """Queries and control for one kind of scheduler."""

    #: Short name, e.g. ``"pbs"``.
    name: str = ""
    #: The environment variable that holds a running job's own ID.
    job_id_variable: str = ""
    #: How the job's variables reach it when there is no settings file, if
    #: other users may read them there (PBS: ``"qsub -v"``, shown by
    #: ``qstat -f``). Empty when they go privately, e.g. in the submitting
    #: process's environment.
    exposed_variables_via: str = ""

    @abstractmethod
    def list_jobs(
        self, user: Optional[str] = None, full_names: bool = False
    ) -> List[LiveJob]:
        """The jobs the scheduler still lists for ``user`` (default: you).

        Finished jobs are not included. With ``full_names`` the names are
        never truncated, which may cost a slower query.

        Raises :class:`SchedulerQueryError` if the scheduler can't be asked.
        """

    @abstractmethod
    def job_info(self, job_ids: Iterable[str]) -> Dict[str, JobInfo]:
        """What the scheduler knows about each job, by job ID.

        Includes finished jobs while the scheduler remembers them. Jobs it
        doesn't know are left out. Never raises. Use this rather than
        :meth:`job_details` unless you need the scheduler's own record.
        """

    @abstractmethod
    def job_details(self, job_ids: Iterable[str]) -> Dict[str, Dict[str, Any]]:
        """The scheduler's own record of each job, by job ID.

        Like :meth:`job_info`, but in the scheduler's own layout (PBS:
        ``parse_qstat_fx_output``, which ``qxub history`` stores). Each
        :class:`JobInfo` carries it as ``native``.
        """

    @abstractmethod
    def cancel(self, job_id: str) -> Tuple[bool, str]:
        """Ask the scheduler to cancel ``job_id``.

        Returns ``(True, "")`` on success, else ``(False, reason)``.
        """

    @abstractmethod
    def looks_like_job_id(self, value: str) -> bool:
        """Whether ``value`` is shaped like one of this scheduler's job IDs.

        Used to tell scheduler jobs from rows qxub records for other things
        (``qxi_...`` interactive sessions, ``qx-...`` run IDs).
        """

    @abstractmethod
    def submit_command(self, spec: JobSpec) -> List[str]:
        """The command line that would submit ``spec`` (for dry runs)."""

    @abstractmethod
    def submit(self, spec: JobSpec) -> str:
        """Submit ``spec`` and return the scheduler's job ID.

        Raises :class:`SubmitError` if the job isn't accepted.
        """

    @abstractmethod
    def interactive_command(self, spec: JobSpec) -> List[str]:
        """The command that runs ``spec.script`` as an interactive session
        on this terminal (PBS: ``qsub -I … -- script``).

        The caller runs it in the foreground, or inside tmux.
        """

    def available(self) -> bool:
        """Whether this machine can ask this scheduler about jobs.

        Jobs of a backend other than the current platform's are listed
        (``backends.live_jobs``) only if it is available.
        """
        return True

    def same_job(self, job_id: str) -> str:
        """A key under which different spellings of one job ID are equal."""
        return job_id

    def run_file(self, job_id: str) -> Optional[str]:
        """Where the job's run file is, if this backend recorded it.

        Otherwise it's looked for in the current platform's status
        directory, which may not be the one the job was submitted with.
        """
        return None

    def live_state(self, state: str) -> JobState:
        """What a :class:`LiveJob`'s ``state`` (the scheduler's own code)
        means, so jobs can be filtered the same way on every scheduler."""
        return _LIVE_STATES.get(state, JobState.UNKNOWN)

    def job_end(self, job_id: str, job_log=None) -> Optional[JobEnd]:
        """How ``job_id`` ended, once the scheduler has ended it; ``None``
        while it is queued or running, or if that can't be told.

        ``qxub exec`` asks this every so often, for a job that ends without
        recording its exit in qxub's status files (killed at its walltime,
        cancelled outside qxub, or unable to start). It should be cheap. By
        default it reads the job log's ``Exit Status:`` line
        (:func:`read_job_log_exit`).
        """
        exit_status = read_job_log_exit(job_log) if job_log else None
        if exit_status is None:
            return None
        return JobEnd(exit_status, self.ended_by(exit_status))

    def ended_by(self, exit_status: Optional[int]) -> Optional[str]:
        """Why the scheduler ended a job, if its exit status says
        (``ENDED_BY_*``); ``None`` if it doesn't.

        For an exit status read outside :meth:`job_info`, such as from the
        job log (PBS: Exit_status -29 is a walltime kill, 271 a qdel).
        """
        return None

    def find_info(self, job_id: str) -> Optional[JobInfo]:
        """:meth:`job_info` for one job, also matching another spelling of
        its ID (PBS: ``123`` for ``123.gadi-pbs``)."""
        return self.find_details(self.job_info([job_id]), job_id)

    def find_details(self, details: Dict[str, Any], job_id: str) -> Optional[Any]:
        """``details[job_id]``, also matching another spelling of the ID."""
        if job_id in details:
            return details[job_id]
        key = self.same_job(job_id)
        for found_id, data in details.items():
            if self.same_job(found_id) == key:
                return data
        return None
