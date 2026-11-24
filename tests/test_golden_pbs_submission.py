"""Golden tests for the ``qsub`` command line that ``qxub exec`` builds.

These pin today's PBS output so the scheduler-neutral refactor (see
``docs/dev/pbs-coupling-audit``) can show it stays byte-identical. Each case
runs ``qxub exec --dry-run`` in a sandbox, normalises the paths, timestamps
and host details that vary between runs, and compares the result with
``tests/golden/exec/<case>.txt``.

The golden files record current behaviour, including known bugs. Fixing
one means updating its golden file in the same change: run
``pytest --update-golden`` and review the diff. Known bugs pinned here:

* ``gadi_hugemem``: queue auto-selection only runs for ``--gpus``,
  ``--internet`` or ``-q auto``, so a 600GB job goes to ``normal``.
"""

import re
import shlex

import pytest
from conftest import REPO_ROOT, dry_run_settings, split_qsub_vars

DRY_RUN_PREFIX = "Dry run - would execute: "
JOBSCRIPTS_DIR = REPO_ROOT / "qxub" / "jobscripts"

# One entry per case: (golden name, extra setup, qxub exec arguments).
CASES = [
    (
        "default_resources",
        None,
        ["-P", "ab12", "-q", "normal", "--mem", "4GB", "--cpus", "2"]
        + ["--runtime", "1h", "--", "echo", "hello"],
    ),
    ("no_project", None, ["--", "echo", "hello"]),
    ("conda_env", None, ["-P", "ab12", "--env", "myenv", "--", "python", "x.py"]),
    (
        "modules",
        None,
        ["-P", "ab12", "--mod", "samtools", "--mod", "bcftools/1.19"]
        + ["--", "samtools", "--version"],
    ),
    ("modules_list", None, ["-P", "ab12", "--mods", "samtools,bcftools", "--", "ls"]),
    (
        "conda_and_modules",
        None,
        ["-P", "ab12", "--env", "myenv", "--mod", "cuda", "--", "nvcc", "-V"],
    ),
    (
        "singularity",
        None,
        ["-P", "ab12", "--sif", "/containers/tool.sif"]
        + ["--bind", "/g/data/ab12,/scratch/ab12", "--", "tool", "run"],
    ),
    (
        "raw_resources",
        None,
        ["-P", "ab12", "-l", "ncpus=8", "-l", "mem=32GB"]
        + ["--disk", "100GB", "--volumes", "gdata/ab12+scratch/ab12", "--", "ls"],
    ),
    (
        "pre_post_vars",
        None,
        ["-P", "ab12", "--pre", "echo start", "--post", "echo done"]
        + ["--var", "FOO=bar", "--vars", "A=1,B=2", "--", "env"],
    ),
    (
        "names_and_paths",
        None,
        ["-P", "ab12", "-N", "myjob", "--out", "logs/out.txt"]
        + ["--err", "logs/err.txt", "--execdir", "/tmp/elsewhere", "--", "pwd"],
    ),
    (
        "email",
        None,
        ["-P", "ab12", "--email", "me@example.org", "--notify", "--", "ls"],
    ),
    (
        "quoted_command",
        None,
        ["-P", "ab12", "--", "bash", "-c", "echo 'a b' | wc -w && echo \"$HOME\""],
    ),
    (
        "gadi_gpu",
        "gadi",
        ["-P", "ab12", "--gpus", "1", "--cpus", "12", "--", "nvidia-smi"],
    ),
    (
        "gadi_hugemem",
        "gadi",
        ["-P", "ab12", "--mem", "600GB", "--cpus", "4", "--", "big_job"],
    ),
    ("gadi_internet", "gadi", ["-P", "ab12", "--internet", "--", "wget", "x"]),
    (
        "config_defaults",
        "config",
        ["--", "echo", "configured"],
    ),
]

CONFIG_YAML = """\
defaults:
  project: cd34
  queue: express
  resources:
    - mem=8GB
    - ncpus=1
    - walltime=00:30:00
"""


def normalise(text: str, sandbox) -> str:
    """Replace the run-specific parts of qxub output with placeholders."""
    replacements = [
        (str(JOBSCRIPTS_DIR), "<JOBSCRIPTS>"),
        (str(sandbox.cwd), "<CWD>"),
        (str(sandbox.home), "<HOME>"),
        (str(sandbox.root), "<SANDBOX>"),
        (str(REPO_ROOT), "<REPO>"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    text = re.sub(r"\d{8}[_-]\d{6}", "<TIMESTAMP>", text)
    text = re.sub(r"qxub_version=[^,\s]+", "qxub_version=<VERSION>", text)
    text = re.sub(r"qxub_path=[^,\s]+", "qxub_path=<QXUB_PATH>", text)
    text = re.sub(r"qxub_python=[^,\s]+", "qxub_python=<PYTHON>", text)
    text = re.sub(r"qxub_submit_host=[^,\s]+", "qxub_submit_host=<HOST>", text)
    text = re.sub(r"qx-[0-9a-f]{8}-[0-9a-f-]{27}", "<RUN_ID>", text)
    return text


def render_qsub_argv(argv: list) -> str:
    """One argument per line, with the ``-v`` list split onto its own lines.

    This keeps golden-file diffs readable. The rendering is only for
    display; any change to the real argument list still changes it.
    """
    lines = []
    args = iter(argv)
    for arg in args:
        if arg == "-v":
            lines.append("-v")
            lines.extend(f"    {item}" for item in split_qsub_vars(next(args)))
        else:
            lines.append(arg)
    return "\n".join(lines) + "\n"


def dry_run_qsub_argv(result) -> list:
    lines = [l for l in result.stdout.splitlines() if l.startswith(DRY_RUN_PREFIX)]
    assert len(lines) == 1, f"no dry-run line in output:\n{result.stdout}"
    return shlex.split(lines[0][len(DRY_RUN_PREFIX) :])


def _prepare(sandbox, setup):
    if setup == "gadi":
        sandbox.use_gadi_platform()
    elif setup == "config":
        sandbox.write_config(CONFIG_YAML)


@pytest.mark.parametrize("name, setup, args", CASES, ids=[c[0] for c in CASES])
def test_exec_dry_run_matches_golden(sandbox, golden, name, setup, args):
    _prepare(sandbox, setup)
    result = sandbox.run("exec", "--dry-run", *args)
    assert result.returncode == 0, result.stdout + result.stderr

    argv = dry_run_qsub_argv(result)
    assert argv[0] == "qsub"
    rendered = render_qsub_argv(argv)
    settings = dry_run_settings(result.stdout)
    if settings:
        rendered += "--- job settings (QXUB_JOB_FILE) ---\n" + settings
    golden(f"exec/{name}.txt", normalise(rendered, sandbox))
    assert sandbox.pbs_calls() == [], "a dry run must not call PBS"


@pytest.mark.parametrize("_name, setup, args", CASES, ids=[c[0] for c in CASES])
def test_terse_submission_sends_the_dry_run_command(sandbox, _name, setup, args):
    """``--terse`` goes through a separate submission path; both must agree."""
    _prepare(sandbox, setup)
    dry = dry_run_qsub_argv(sandbox.run("exec", "--dry-run", *args))

    result = sandbox.run("exec", "--terse", *args)
    assert result.returncode == 0, result.stdout + result.stderr
    (call,) = sandbox.pbs_calls("qsub")

    def strip_timestamps(argv):
        return [normalise(a, sandbox) for a in argv]

    assert strip_timestamps(["qsub", *call["argv"]]) == strip_timestamps(dry)
