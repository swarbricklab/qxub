"""The interface qxub uses to ask a scheduler about jobs and to cancel them.

Each scheduler (PBS Pro, and later Slurm, a local runner and others) gets
one :class:`SchedulerBackend`. The platform definition's ``type`` picks
which one (see :func:`qxub.backends.get_backend`).

A job to submit, or an interactive session to start, is described by a
:class:`JobSpec`; the backend turns it into its own command line, including the :class:`ResourceSpec` (PBS writes
``-l ncpus=…``).
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
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
    def job_details(self, job_ids: Iterable[str]) -> Dict[str, Dict[str, Any]]:
        """Everything the scheduler knows about each job, by job ID.

        Includes finished jobs while the scheduler remembers them. Jobs it
        doesn't know are left out. Never raises.

        The value has the shape ``parse_qstat_fx_output`` returns
        (``exit_status``, ``resources_requested``, ``resources_used``,
        ``timing``, ``execution``, ``metadata``, ``efficiency``). It is still
        PBS-flavoured; a neutral ``JobUsage`` is later work.
        """

    @abstractmethod
    def is_finished(self, details: Optional[Dict[str, Any]]) -> bool:
        """Whether ``details`` (from :meth:`job_details`) is a finished job."""

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

    def same_job(self, job_id: str) -> str:
        """A key under which different spellings of one job ID are equal."""
        return job_id

    def find_details(
        self, details: Dict[str, Dict[str, Any]], job_id: str
    ) -> Optional[Dict[str, Any]]:
        """``details[job_id]``, also matching another spelling of the ID."""
        if job_id in details:
            return details[job_id]
        key = self.same_job(job_id)
        for found_id, data in details.items():
            if self.same_job(found_id) == key:
                return data
        return None
