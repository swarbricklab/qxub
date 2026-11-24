"""How a remote qxub tells the local one that it has submitted a job.

If the SSH connection drops, the local side retries the remote command
unless a job was already submitted, which a retry would submit again. The
remote qxub says so with a line on stderr::

    QXUB-SUBMITTED job_id=<ID> backend=<name>

whatever the scheduler's IDs look like. The local side asks for it by
setting :data:`ENV` in the remote command, and hides the line from the
user. Remote qxub versions that don't announce are recognised by their
PBS-style job IDs, as before.
"""

import os
import re
import sys
from typing import Optional, Tuple

#: Set (to ``1``) in the remote command's environment to ask for the line.
ENV = "QXUB_REMOTE_HANDSHAKE"
MARKER = "QXUB-SUBMITTED"

_LINE_RE = re.compile(rf"^{MARKER} job_id=(\S+)(?: backend=(\S+))?\s*$", re.M)


def requested() -> bool:
    """Whether to announce submissions, which only the first may.

    Clears the request, so jobs started from this process (whose
    environment the local backend passes on) don't announce theirs.
    """
    return os.environ.pop(ENV, "") == "1"


def announce(job_id: str, backend: str) -> None:
    """Tell the local qxub that ``job_id`` has been submitted."""
    print(f"{MARKER} job_id={job_id} backend={backend}", file=sys.stderr, flush=True)


def is_marker(line: str) -> bool:
    return line.startswith(MARKER + " ")


def find(text: str) -> Optional[Tuple[str, Optional[str]]]:
    """``(job_id, backend)`` from the first announcement in ``text``."""
    match = _LINE_RE.search(text or "")
    return (match.group(1), match.group(2)) if match else None
