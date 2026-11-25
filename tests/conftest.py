"""Shared pytest fixtures for the qxub test suite."""

import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import pytest

# Upper bound for a cold `import qxub`, in seconds.
IMPORT_TIME_LIMIT = 1.0


# ---------------------------------------------------------------------------
# The whole session runs away from the machine's own qxub setup
# ---------------------------------------------------------------------------


def pytest_configure(config):
    """Give the test session its own home, config and ``qxub`` command.

    Runs before any test module imports qxub, so qxub's global config sees
    only this. Without it, a run on a cluster login node (Gadi) picked up
    the user's config, ``$PROJECT``, the site's platform files and
    ``QXUB_SHARED_DB``, wrote history and job-database rows into the user's
    real ``~/.config/qxub``, and tests that run ``qxub`` by name failed when
    the CLI wasn't installed on ``PATH``. :class:`QxubSandbox` isolates
    each subprocess run the same way.
    """
    import atexit
    import shutil
    import site
    import tempfile

    root = Path(tempfile.mkdtemp(prefix="qxub-tests-"))
    atexit.register(shutil.rmtree, root, ignore_errors=True)
    home = root / "home"
    bin_dir = root / "bin"
    for path in (home / ".config", root / "xdg", bin_dir):
        path.mkdir(parents=True)

    # Subprocess Pythons still find packages installed with --user.
    os.environ.setdefault("PYTHONUSERBASE", site.getuserbase())
    # The site's settings go; test switches (QXUB_TEST_REMOTE=1) stay.
    for name in list(os.environ):
        if name.startswith("QXUB_TEST_"):
            continue
        if name.startswith(("QXUB_", "PBS_", "SLURM_")) or name == "PROJECT":
            del os.environ[name]
    os.environ.update(
        HOME=str(home),
        XDG_CONFIG_HOME=str(home / ".config"),
        XDG_CONFIG_DIRS=str(root / "xdg"),
        XDG_STATE_HOME=str(home / ".local" / "state"),
        XDG_DATA_HOME=str(home / ".local" / "share"),
        # Never match a platform's host pattern by accident (e.g. on Gadi)
        HOSTNAME="qxub-test-host",
        # No background consolidation rounds (qxub.queue.consolidate)
        QXUB_AUTO_CONSOLIDATE="0",
    )

    # Every Python the tests start imports this checkout's qxub first, not
    # an installed one the interpreter can also see (#163 review: a venv on
    # Gadi saw qxub 3.9.1 in another environment's site-packages).
    os.environ["PYTHONPATH"] = os.pathsep.join(
        p for p in [str(REPO_ROOT), os.environ.get("PYTHONPATH", "")] if p
    )

    # `qxub` on PATH is this checkout's CLI, run by the interpreter under
    # test, whether or not (and whichever) qxub is installed.
    shim = bin_dir / "qxub"
    shim.write_text(
        f"#!{sys.executable}\n"
        f"import sys; sys.path.insert(0, {str(REPO_ROOT)!r})\n"
        f"{_CLI_SHIM.replace('; ', chr(10))}\n"
    )
    shim.chmod(0o755)
    os.environ["PATH"] = os.pathsep.join([str(bin_dir), os.environ.get("PATH", "")])


def measure_cold_import_time(module: str = "qxub") -> float:
    """Time importing ``module`` in a fresh interpreter.

    Timing ``import qxub`` inside the test process only measures a
    ``sys.modules`` cache hit, so the import runs in a subprocess instead.
    """
    code = (
        "import time; t = time.perf_counter(); "
        f"import {module}; print(time.perf_counter() - t)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    )
    return float(result.stdout.strip().splitlines()[-1])


def fastest_of(measure, limit: float, attempts: int = 3) -> float:
    """The fastest of up to ``attempts`` timings from ``measure()``.

    Stops at the first one under ``limit``. A single wall-clock timing on a
    shared login node fails whenever the node is busy (#163 review: 3.02 s
    against a 3.0 s limit at load 11); a real slowdown is slow every time.
    """
    best = float("inf")
    for _ in range(attempts):
        best = min(best, measure())
        if best < limit:
            break
    return best


def timed_run(argv, limit: float, timeout=None, attempts: int = 3):
    """Run ``argv`` up to ``attempts`` times (see :func:`fastest_of`).

    Returns the fastest run's result and time. A run that hits ``timeout``
    (default: twice ``limit``) counts as a slow try, not an error; only if
    every try times out does the test fail.
    """
    timeout = timeout if timeout is not None else 2 * limit
    runs = []  # (seconds, result or None if it timed out)

    def measure():
        start = time.perf_counter()
        try:
            result = subprocess.run(
                argv, capture_output=True, text=True, timeout=timeout
            )
        except subprocess.TimeoutExpired:
            runs.append((float(timeout), None))
        else:
            runs.append((time.perf_counter() - start, result))
        return runs[-1][0]

    fastest_of(measure, limit, attempts)
    elapsed, result = min(runs, key=lambda run: run[0])
    if result is None:
        pytest.fail(f"{' '.join(argv)} timed out after {timeout}s on every try")
    return result, elapsed


@pytest.fixture
def cold_import_time():
    """Seconds taken by a cold ``import qxub`` (the fastest of a few tries)."""
    return fastest_of(measure_cold_import_time, IMPORT_TIME_LIMIT)


# ---------------------------------------------------------------------------
# Sandboxed qxub runs against a fake PBS
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
FAKE_PBS_BIN = Path(__file__).resolve().parent / "fakes" / "pbs"
FAKE_SLURM_BIN = Path(__file__).resolve().parent / "fakes" / "slurm"
GADI_PLATFORM_DIR = REPO_ROOT / "docs" / "platforms"
SLURM_PLATFORM = GADI_PLATFORM_DIR / "slurm" / "slurm.yaml"

# Run the CLI the way the console script does, but with a stable prog name.
_CLI_SHIM = "import sys; sys.argv[0] = 'qxub'; from qxub.cli import qxub; qxub()"


class QxubSandbox:
    """Runs ``qxub`` in a subprocess with an isolated home and a fake PBS.

    Each sandbox has its own ``HOME``, XDG config directories and working
    directory, so no user or system qxub config leaks in. The fake ``qsub``,
    ``qstat`` and ``qdel`` from ``tests/fakes/pbs``, and the fake Slurm
    commands from ``tests/fakes/slurm``, come first on ``PATH``, so nothing
    can reach a real scheduler even on a cluster login node.
    """

    def __init__(self, root: Path):
        self.root = root.resolve()
        self.home = self.root / "home"
        self.cwd = self.root / "work"
        self.pbs_state = self.root / "pbs"
        self.slurm_state = self.root / "slurm"
        for path in (self.home, self.cwd, self.pbs_state, self.root / "xdg"):
            path.mkdir(parents=True, exist_ok=True)

        env = {
            k: v
            for k, v in os.environ.items()
            if not k.startswith(("QXUB_", "PBS_", "SLURM_", "FAKE_"))
            and k not in ("PROJECT", "HOSTNAME", "XDG_STATE_HOME")
        }
        env.update(
            HOME=str(self.home),
            XDG_CONFIG_HOME=str(self.home / ".config"),
            XDG_CONFIG_DIRS=str(self.root / "xdg"),
            PATH=os.pathsep.join(
                [str(FAKE_PBS_BIN), str(FAKE_SLURM_BIN), env.get("PATH", "")]
            ),
            FAKE_PBS_STATE=str(self.pbs_state),
            FAKE_SLURM_STATE=str(self.slurm_state),
            NO_COLOR="1",
            XDG_STATE_HOME=str(self.home / ".local" / "state"),
            # Never match a platform's host pattern by accident (e.g. on Gadi)
            HOSTNAME="qxub-test-host",
            # No background consolidation rounds (qxub.queue.consolidate)
            QXUB_AUTO_CONSOLIDATE="0",
        )
        self.env = env

    @property
    def config_dir(self) -> Path:
        return self.home / ".config" / "qxub"

    def write_config(self, text: str) -> Path:
        """Write the user ``config.yaml``."""
        self.config_dir.mkdir(parents=True, exist_ok=True)
        path = self.config_dir / "config.yaml"
        path.write_text(text)
        return path

    def use_gadi_platform(self, extra_config: str = "") -> None:
        """Select the bundled ``nci_gadi`` platform as the default platform."""
        self.env["QXUB_PLATFORM_PATHS"] = str(GADI_PLATFORM_DIR)
        definition = GADI_PLATFORM_DIR / "nci_gadi.yaml"
        self.write_config(
            "defaults:\n"
            "  platform: nci_gadi\n"
            "platforms:\n"
            "  nci_gadi:\n"
            f"    definition: {definition}\n" + extra_config
        )

    def use_slurm_platform(self) -> None:
        """Select the example Slurm platform (``docs/platforms/slurm``)."""
        self.env["QXUB_PLATFORM_FILE"] = str(SLURM_PLATFORM)

    def run(self, *args: str, env: dict = None, input: str = None, timeout=60):
        """Run ``qxub *args`` and return the ``CompletedProcess``."""
        return subprocess.run(
            [sys.executable, "-c", _CLI_SHIM, *args],
            cwd=self.cwd,
            env={**self.env, **(env or {})},
            input=input,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

    # -- fake PBS controls -------------------------------------------------

    def pbs_calls(self, command: str = None) -> list:
        """Calls the fake PBS commands received, oldest first."""
        log = self.pbs_state / "calls.jsonl"
        if not log.exists():
            return []
        calls = [json.loads(line) for line in log.read_text().splitlines()]
        return [c for c in calls if command is None or c["cmd"] == command]

    def slurm_calls(self, command: str = None) -> list:
        """Calls the fake Slurm commands received, oldest first."""
        log = self.slurm_state / "calls.jsonl"
        if not log.exists():
            return []
        calls = [json.loads(line) for line in log.read_text().splitlines()]
        return [c for c in calls if command is None or c["cmd"] == command]

    def slurm_live_jobs(self) -> dict:
        """Jobs the fake Slurm still lists: ``{job_id: {name, queue, state}}``."""
        path = self.slurm_state / "live.json"
        return json.loads(path.read_text()) if path.exists() else {}

    def set_sacct(self, job_id: str, text: str) -> None:
        """Set what ``sacct --jobs=job_id`` prints."""
        sacct_dir = self.slurm_state / "sacct"
        sacct_dir.mkdir(parents=True, exist_ok=True)
        (sacct_dir / f"{job_id}.txt").write_text(text)

    def fail_qsub(self, exit_code: int, stderr: str) -> None:
        self.env["FAKE_QSUB_EXIT"] = str(exit_code)
        self.env["FAKE_QSUB_STDERR"] = stderr

    def fail_qdel(self, exit_code: int, stderr: str) -> None:
        self.env["FAKE_QDEL_EXIT"] = str(exit_code)
        self.env["FAKE_QDEL_STDERR"] = stderr

    def finish_job(self, job_id: str) -> None:
        """The scheduler forgets ``job_id``, as when it ends."""
        path = self.pbs_state / "live.json"
        live = json.loads(path.read_text()) if path.exists() else {}
        live.pop(job_id, None)
        path.write_text(json.dumps(live))

    def live_jobs(self) -> dict:
        """Jobs the fake scheduler still lists: ``{job_id: {name, queue, state}}``."""
        path = self.pbs_state / "live.json"
        return json.loads(path.read_text()) if path.exists() else {}

    def set_qstat(self, job_id: str, text: str) -> None:
        """Set what ``qstat -fx job_id`` prints."""
        qstat_dir = self.pbs_state / "qstat"
        qstat_dir.mkdir(exist_ok=True)
        (qstat_dir / f"{job_id}.txt").write_text(text)


@pytest.fixture
def sandbox(tmp_path):
    """A :class:`QxubSandbox` rooted in a fresh temporary directory."""
    return QxubSandbox(tmp_path)


# ---------------------------------------------------------------------------
# Golden files
# ---------------------------------------------------------------------------

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Rewrite golden files under tests/golden from current output.",
    )


@pytest.fixture
def golden(request):
    """Compare text with ``tests/golden/<name>``.

    Run ``pytest --update-golden`` to (re)write the files after an
    intentional change, then review the diff before committing.
    """
    update = request.config.getoption("--update-golden")

    def check(name: str, actual: str) -> None:
        path = GOLDEN_DIR / name
        if update:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(actual)
            return
        if not path.exists():
            pytest.fail(f"Missing golden file {path}; run pytest --update-golden")
        assert actual == path.read_text(), f"Output differs from golden file {path}"

    return check


# -- the variables a job is given ----------------------------------------------

SETTINGS_HEADER = "Job settings ("


def split_qsub_vars(value: str) -> list:
    """Split a ``qsub -v`` list on the commas PBS treats as separators.

    Commas inside single or double quotes belong to the value, as in
    ``bind='/a,/b'``.
    """
    items, current, quote = [], "", None
    for char in value:
        if quote:
            quote = None if char == quote else quote
        elif char in "'\"":
            quote = char
        elif char == ",":
            items.append(current)
            current = ""
            continue
        current += char
    items.append(current)
    return items


def parse_settings(text: str) -> dict:
    """The variables in a job settings file (``export NAME=value`` lines)."""
    found = {}
    tokens = iter(shlex.split(text, comments=True))  # values may span lines
    for token in tokens:
        if token == "export":
            name, _, value = next(tokens).partition("=")
            found[name] = value
    return found


def dry_run_settings(stdout: str) -> str:
    """The job settings a dry run prints after its command line."""
    lines = stdout.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(SETTINGS_HEADER):
            return "\n".join(lines[i + 1 :]) + "\n"
    return ""


def job_variables(argv, dry_run_stdout=None) -> dict:
    """The variables a ``qsub`` argv gives the job.

    Either the ``-v`` list itself or, when it names a settings file
    (``QXUB_JOB_FILE``), that file's variables. A dry run doesn't write the
    file, so pass its stdout instead.
    """
    if "-v" not in argv:
        return {}
    found = {}
    for item in split_qsub_vars(argv[argv.index("-v") + 1]):
        name, _, value = item.partition("=")
        if len(value) > 1 and value[0] == value[-1] and value[0] in "'\"":
            value = value[1:-1]
        found[name] = value
    if "QXUB_JOB_FILE" in found:
        if dry_run_stdout is not None:
            text = dry_run_settings(dry_run_stdout)
        else:
            text = Path(found["QXUB_JOB_FILE"]).read_text()
        found.update(parse_settings(text))
    return found
