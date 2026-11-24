"""Scheduler backends: how qxub asks a scheduler about jobs and cancels them.

Use :func:`get_backend` rather than a backend class, so the platform
definition decides which scheduler qxub talks to::

    from qxub import backends

    for job in backends.get_backend().list_jobs():
        ...
"""

import logging
from typing import Dict, Optional, Type

from ..resources.spec import ResourceSpec
from .base import JobSpec, LiveJob, SchedulerBackend, SchedulerQueryError, SubmitError
from .pbs import PbsBackend

logger = logging.getLogger(__name__)

# Platform ``type`` values and the backend each selects.
BACKENDS: Dict[str, Type[SchedulerBackend]] = {
    "pbs": PbsBackend,
    "pbs_pro": PbsBackend,
    "pbspro": PbsBackend,
    "torque": PbsBackend,
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


__all__ = [
    "BACKENDS",
    "JobSpec",
    "LiveJob",
    "PbsBackend",
    "ResourceSpec",
    "SchedulerBackend",
    "SchedulerQueryError",
    "SubmitError",
    "backend_for_type",
    "get_backend",
]
