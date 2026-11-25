"""Shared logic for the fake Slurm commands in this directory.

``sbatch``, ``squeue``, ``sacct`` and ``scancel`` stand in for Slurm's
client tools, so the tests can submit, list, look up and cancel Slurm jobs
without a cluster. Like ``tests/fakes/pbs``, behaviour is driven by
environment variables:

``FAKE_SLURM_STATE``  Directory for state (required). Every call is
                      appended to ``calls.jsonl`` there, with the
                      ``QXUB_JOB_FILE`` sbatch was given in its environment.
``FAKE_SBATCH_EXIT``  Exit code for ``sbatch`` (default 0); stderr is
                      ``FAKE_SBATCH_STDERR``. Likewise ``FAKE_SCANCEL_*``
                      and ``FAKE_SQUEUE_*``.

``sbatch`` adds the job to ``live.json`` (pending, with its
``--job-name`` and ``--partition``) and prints its ID; ``scancel`` removes
it; ``squeue`` lists live jobs as ``%i|%P|%t|%j``. ``sacct --jobs=ID,...``
prints ``$FAKE_SLURM_STATE/sacct/ID.txt`` for each job that has one.
"""

import json
import os
import sys
from pathlib import Path

FIRST_JOB_NUMBER = 5000


def _state_dir() -> Path:
    state = os.environ.get("FAKE_SLURM_STATE")
    if not state:
        sys.stderr.write("fake Slurm: FAKE_SLURM_STATE is not set\n")
        sys.exit(2)
    path = Path(state)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _record(state: Path, command: str, argv: list) -> None:
    entry = {"cmd": command, "argv": argv, "cwd": os.getcwd()}
    if "QXUB_JOB_FILE" in os.environ:
        entry["QXUB_JOB_FILE"] = os.environ["QXUB_JOB_FILE"]
    with open(state / "calls.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry) + "\n")


def _fail(env_prefix: str) -> None:
    code = int(os.environ.get(f"{env_prefix}_EXIT", "0"))
    if code:
        sys.stderr.write(os.environ.get(f"{env_prefix}_STDERR", "") + "\n")
        sys.exit(code)


def _option(argv: list, name: str, default: str = "") -> str:
    for arg in argv:
        if arg.startswith(f"--{name}="):
            return arg.split("=", 1)[1]
    return default


def _live(state: Path) -> dict:
    path = state / "live.json"
    return json.loads(path.read_text()) if path.exists() else {}


def _save_live(state: Path, live: dict) -> None:
    (state / "live.json").write_text(json.dumps(live))


def sbatch(argv: list) -> None:
    state = _state_dir()
    _record(state, "sbatch", argv)
    _fail("FAKE_SBATCH")
    counter = state / "next_job_number"
    number = int(counter.read_text()) if counter.exists() else FIRST_JOB_NUMBER
    counter.write_text(str(number + 1))
    live = _live(state)
    live[str(number)] = {
        "name": _option(argv, "job-name", "sbatch"),
        "queue": _option(argv, "partition", "cpu"),
        "state": "PD",
    }
    _save_live(state, live)
    print(number)


def squeue(argv: list) -> None:
    state = _state_dir()
    _record(state, "squeue", argv)
    _fail("FAKE_SQUEUE")
    for job_id, job in _live(state).items():
        print(f"{job_id}|{job['queue']}|{job['state']}|{job['name']}")


def sacct(argv: list) -> None:
    state = _state_dir()
    _record(state, "sacct", argv)
    for job_id in _option(argv, "jobs").split(","):
        canned = state / "sacct" / f"{job_id}.txt"
        if job_id and canned.exists():
            sys.stdout.write(canned.read_text().rstrip("\n") + "\n")


def scancel(argv: list) -> None:
    state = _state_dir()
    _record(state, "scancel", argv)
    _fail("FAKE_SCANCEL")
    live = _live(state)
    for job_id in argv:
        live.pop(job_id, None)
    _save_live(state, live)


COMMANDS = {"sbatch": sbatch, "squeue": squeue, "sacct": sacct, "scancel": scancel}


def main() -> None:
    COMMANDS[Path(sys.argv[0]).name](sys.argv[1:])
