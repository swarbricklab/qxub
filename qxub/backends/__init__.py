"""Scheduler backends: how qxub asks a scheduler about jobs and cancels them.

Use :func:`get_backend` rather than a backend class, so the platform
definition decides which scheduler qxub talks to::

    from qxub import backends

    for job in backends.get_backend().list_jobs():
        ...
"""

import logging
from typing import Dict, List, Optional, Tuple, Type

from ..resources.spec import ResourceSpec
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
    describe_ended_by,
    read_job_log_exit,
)
from .local import LocalBackend
from .pbs import PbsBackend
from .slurm import SlurmBackend

logger = logging.getLogger(__name__)

# Platform ``type`` values and the backend each selects.
BACKENDS: Dict[str, Type[SchedulerBackend]] = {
    "pbs": PbsBackend,
    "pbs_pro": PbsBackend,
    "pbspro": PbsBackend,
    "torque": PbsBackend,
    "slurm": SlurmBackend,
    # Jobs as background processes on this machine
    "local": LocalBackend,
}
DEFAULT_TYPE = "pbs"


def backend_for_type(scheduler_type: Optional[str]) -> SchedulerBackend:
    """The backend for a platform ``type`` (PBS if empty or unknown)."""
    kind = (scheduler_type or DEFAULT_TYPE).strip().lower()
    cls = BACKENDS.get(kind)
    if cls is None:
        logger.warning(
            "No scheduler backend for platform type %r (known: %s); using PBS",
            scheduler_type,
            ", ".join(sorted(BACKENDS)),
        )
        cls = BACKENDS[DEFAULT_TYPE]
    return cls()


def get_backend(platform=None) -> SchedulerBackend:
    """The backend for ``platform``, or for the current platform.

    The platform definition's ``type`` (``pbs_pro`` for NCI Gadi) picks the
    backend. Without a platform definition qxub assumes PBS, as it always
    has.
    """
    if platform is None:
        try:
            from ..platforms import (  # pylint: disable=import-outside-toplevel
                get_current_platform,
            )

            platform = get_current_platform()
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("No current platform (%s); assuming PBS", exc)
    return backend_for_type(getattr(platform, "type", None))


def backend_for_job(
    job_id: str, backend_name: Optional[str] = None
) -> SchedulerBackend:
    """The backend that owns ``job_id``.

    ``backend_name`` is the backend recorded when the job was submitted (the
    job database's and run file's ``backend``); use it when known. Otherwise
    the current platform's backend if it recognises the ID, else the first
    backend that does (``local-...`` is the local backend's), else the
    current one. IDs that look alike across schedulers (PBS and Slurm are
    both numbers) need the recorded name.
    """
    if backend_name:
        return backend_for_type(backend_name)
    current = get_backend()
    if current.looks_like_job_id(job_id):
        return current
    for cls in dict.fromkeys(BACKENDS.values()):
        if not isinstance(current, cls):
            candidate = cls()
            if candidate.looks_like_job_id(job_id):
                return candidate
    return current


def live_jobs(
    user: Optional[str] = None, full_names: bool = False
) -> List[Tuple[SchedulerBackend, LiveJob]]:
    """Queued and running jobs of every backend this machine can ask.

    The current platform's backend first: if it can't be asked,
    :class:`SchedulerQueryError` is raised. Then each other backend that is
    :meth:`~SchedulerBackend.available` here (local jobs under the Gadi
    platform, PBS jobs under the local one); one that fails is skipped.
    """
    current = get_backend()
    found = [(current, job) for job in current.list_jobs(user, full_names)]
    for cls in dict.fromkeys(BACKENDS.values()):
        if isinstance(current, cls):
            continue
        other = cls()
        try:
            if other.available():
                found += [(other, job) for job in other.list_jobs(user, full_names)]
        except SchedulerQueryError as exc:
            logger.debug("Skipping %s jobs: %s", other.name, exc)
    return found


__all__ = [
    "BACKENDS",
    "ENDED_BY_CANCEL",
    "ENDED_BY_MEMORY",
    "ENDED_BY_WALLTIME",
    "JobEnd",
    "describe_ended_by",
    "JobInfo",
    "JobSpec",
    "JobState",
    "JobUsage",
    "LiveJob",
    "LocalBackend",
    "PbsBackend",
    "ResourceSpec",
    "SchedulerBackend",
    "SchedulerQueryError",
    "SlurmBackend",
    "SubmitError",
    "backend_for_job",
    "backend_for_type",
    "get_backend",
    "live_jobs",
    "read_job_log_exit",
]
