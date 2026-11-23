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

``qstat -fx JOBID...`` prints ``$FAKE_PBS_STATE/qstat/JOBID.txt`` for each
job that has one, and fails the way PBS does if any job is unknown.

``qsub`` adds the job to ``live.json`` (queued, with its ``-N`` name and
``-q`` queue); ``qdel`` removes it. ``qstat -w -u USER`` lists those jobs
the way PBS Pro does, with names cut to 15 characters. ``FAKE_QSTAT_EXIT``
makes that listing fail.
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


def _option(argv: list, flag: str, default: str = "") -> str:
    return argv[argv.index(flag) + 1] if flag in argv[:-1] else default


def _live(state: Path) -> dict:
    path = state / "live.json"
    return json.loads(path.read_text()) if path.exists() else {}


def _save_live(state: Path, live: dict) -> None:
    (state / "live.json").write_text(json.dumps(live))


def qsub(argv: list) -> None:
    state = _state_dir()
    _record(state, "qsub", argv)
    _fail("FAKE_QSUB")
    job_id = _next_job_id(state)
    live = _live(state)
    live[job_id] = {
        "name": _option(argv, "-N", "STDIN"),
        "queue": _option(argv, "-q", "normal"),
        "state": "Q",
    }
    _save_live(state, live)
    print(job_id)


def _user_listing(user: str, live: dict) -> None:
    print("\ngadi-pbs:")
    print(" " * 99 + "Req'd  Req'd   Elap")
    print(
        "Job ID                         Username        Queue           Jobname"
        "         SessID   NDS  TSK   Memory Time  S Time"
    )
    print("-" * 30 + " " + "-" * 15 + " " + "-" * 15 + " " + "-" * 15 + " ...")
    for job_id, job in live.items():
        name = job["name"] if len(job["name"]) <= 15 else job["name"][:14] + "*"
        print(
            f"{job_id:<30} {user:<15} {job['queue'] + '-exec':<15} {name:<15}"
            f"       --    1     1    4gb 01:00 {job['state']}   --"
        )


def qstat(argv: list) -> None:
    state = _state_dir()
    _record(state, "qstat", argv)
    if "-u" in argv:
        _fail("FAKE_QSTAT")
        _user_listing(_option(argv, "-u"), _live(state))
        return
    job_ids = [a for a in argv if not a.startswith("-")]
    if not job_ids:
        if "-f" in argv:  # every live job, full format (names not cut)
            for job_id, job in _live(state).items():
                print(f"Job Id: {job_id}\n    Job_Name = {job['name']}\n")
        return
    unknown = False
    for job_id in job_ids:
        canned = state / "qstat" / f"{job_id}.txt"
        if canned.exists():
            sys.stdout.write(canned.read_text().rstrip("\n") + "\n\n")
        else:
            sys.stderr.write(f"qstat: Unknown Job Id {job_id}\n")
            unknown = True
    if unknown:
        sys.exit(153)


def qdel(argv: list) -> None:
    state = _state_dir()
    _record(state, "qdel", argv)
    _fail("FAKE_QDEL")
    live = _live(state)
    for job_id in argv:
        live.pop(job_id, None)
    _save_live(state, live)


COMMANDS = {"qsub": qsub, "qstat": qstat, "qdel": qdel}


def main() -> None:
    name = Path(sys.argv[0]).name
    COMMANDS[name](sys.argv[1:])
