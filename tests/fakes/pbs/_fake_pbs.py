"""Shared logic for the fake PBS commands in this directory.

The ``qsub``, ``qstat`` and ``qdel`` scripts next to this file stand in for
the real PBS Pro client tools so the test suite can exercise qxub's
submission, status and cancel paths without a scheduler.

Behaviour is driven by environment variables:

``FAKE_PBS_STATE``   Directory for state (required). Every call is appended
                     to ``calls.jsonl`` there; job IDs come from a counter.
``FAKE_QSUB_EXIT``   Exit code for ``qsub`` (default 0). On a non-zero exit,
                     ``FAKE_QSUB_STDERR`` is written to stderr.
``FAKE_QDEL_EXIT``   Exit code for ``qdel`` (default 0); stderr is
                     ``FAKE_QDEL_STDERR``.

``qstat -fx JOBID`` prints ``$FAKE_PBS_STATE/qstat/JOBID.txt`` if it exists,
otherwise it fails the way PBS does for an unknown job.
"""

import json
import os
import sys
from pathlib import Path

JOB_SUFFIX = "gadi-pbs"
FIRST_JOB_NUMBER = 1000


def _state_dir() -> Path:
    state = os.environ.get("FAKE_PBS_STATE")
    if not state:
        sys.stderr.write("fake PBS: FAKE_PBS_STATE is not set\n")
        sys.exit(2)
    path = Path(state)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _record(state: Path, command: str, argv: list) -> None:
    entry = {"cmd": command, "argv": argv, "cwd": os.getcwd()}
    with open(state / "calls.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry) + "\n")


def _next_job_id(state: Path) -> str:
    counter = state / "next_job_number"
    number = int(counter.read_text()) if counter.exists() else FIRST_JOB_NUMBER
    counter.write_text(str(number + 1))
    return f"{number}.{JOB_SUFFIX}"


def _fail(env_prefix: str) -> None:
    code = int(os.environ.get(f"{env_prefix}_EXIT", "0"))
    if code:
        sys.stderr.write(os.environ.get(f"{env_prefix}_STDERR", "") + "\n")
        sys.exit(code)


def qsub(argv: list) -> None:
    state = _state_dir()
    _record(state, "qsub", argv)
    _fail("FAKE_QSUB")
    print(_next_job_id(state))


def qstat(argv: list) -> None:
    state = _state_dir()
    _record(state, "qstat", argv)
    job_ids = [a for a in argv if not a.startswith("-")]
    if not job_ids:
        return
    canned = state / "qstat" / f"{job_ids[-1]}.txt"
    if canned.exists():
        sys.stdout.write(canned.read_text())
        return
    sys.stderr.write(f"qstat: Unknown Job Id {job_ids[-1]}\n")
    sys.exit(153)


def qdel(argv: list) -> None:
    state = _state_dir()
    _record(state, "qdel", argv)
    _fail("FAKE_QDEL")


COMMANDS = {"qsub": qsub, "qstat": qstat, "qdel": qdel}


def main() -> None:
    name = Path(sys.argv[0]).name
    COMMANDS[name](sys.argv[1:])
