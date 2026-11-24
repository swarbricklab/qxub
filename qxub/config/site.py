"""Site conventions taken from the current platform definition.

Values such as the user's default project, the shared scratch area and extra
scheduler directives differ between sites. They used to be hard-coded for
NCI Gadi (project ``a56``, ``/scratch/$PROJECT/$USER``, ``-W umask=0027``).
They now come from the platform file for the current host (see
``docs/platforms/nci_gadi.yaml``), and each helper has a neutral fallback
for when no platform is defined.
"""

import getpass
import logging
import os
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


def current_platform():
    """The platform for this host, or ``None`` if none is defined."""
    try:
        from ..platforms import get_current_platform

        return get_current_platform()
    except Exception as e:  # pylint: disable=broad-except
        logger.debug("Could not determine the current platform: %s", e)
        return None


def _user() -> str:
    return os.getenv("USER") or getpass.getuser()


def default_project() -> Optional[str]:
    """The project/account from the platform's ``account_env`` variable."""
    platform = current_platform()
    if platform and platform.account_env:
        return os.getenv(platform.account_env) or None
    return None


def scratch_dir() -> Optional[Path]:
    """The per-user shared scratch directory, if the platform defines one.

    Returns ``None`` when there is no platform, it has no ``paths.scratch``,
    or the template needs a project and none is known.
    """
    platform = current_platform()
    template = platform.paths.get("scratch") if platform else None
    if not template:
        return None
    project = default_project()
    if "{project}" in template and not project:
        return None
    return Path(template.format(project=project or "", user=_user()))


def status_dir_storage(job_project: Optional[str]) -> Optional[str]:
    """The storage request a job charged to ``job_project`` needs to reach
    :func:`default_status_dir`, or ``None`` if it needs none.

    On Gadi the status directory is under the *default* project's scratch,
    and a job only sees its own project's scratch plus what ``storage=``
    asks for. The platform's ``paths.scratch_storage`` names the request.
    """
    platform = current_platform()
    template = platform.paths.get("scratch_storage") if platform else None
    project = default_project()
    if not template or not project or not job_project or job_project == project:
        return None
    if scratch_dir() is None:
        return None
    return template.format(project=project, user=_user())


def default_status_dir() -> Path:
    """Where status files for interactive sessions and monitoring live.

    ``{scratch}/qxub/status`` when the platform has a scratch area (it must be
    visible from compute nodes), otherwise ``$XDG_STATE_HOME/qxub/status``.
    """
    scratch = scratch_dir()
    if scratch is not None:
        return scratch / "qxub" / "status"
    return state_dir() / "status"


def state_dir() -> Path:
    """``$XDG_STATE_HOME/qxub``: this machine's state, whatever the platform."""
    state_home = os.getenv("XDG_STATE_HOME") or str(Path.home() / ".local" / "state")
    return Path(state_home) / "qxub"


def default_directives() -> List[str]:
    """Scheduler arguments the platform adds to every submission."""
    platform = current_platform()
    return list(platform.default_directives) if platform else []
