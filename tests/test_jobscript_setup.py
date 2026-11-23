"""Job scripts fail loudly when a setup step fails (#109).

A missing execution directory used to end the job with a bare ``cd`` error
in the joblog and an opaque exit code 1. Every setup failure now says what
failed and how to fix it, on stderr and in the job's ``$err`` file, which
``qxub exec`` and CI steps stream.
"""

import base64
import subprocess
from pathlib import Path

import pytest
from conftest import REPO_ROOT

JOBSCRIPTS = REPO_ROOT / "qxub" / "jobscripts"
TEMPLATES = ["qdefault", "qmod", "qconda", "qconda_mod", "qsing"]

FAKE_MODULE = """#!/bin/bash
# module load NAME: fails for names starting with "bad"
[ "$1" = load ] && case "$2" in bad*) echo "ERROR: Unable to locate a modulefile for '$2'" >&2; exit 1;; esac
exit 0
"""

FAKE_CONDA = """#!/bin/bash
if [ "$1 $2" = "shell.bash hook" ]; then
    cat <<'HOOK'
conda() {
    if [ "$1" = activate ]; then
        [ "$2" = missing ] && { echo "EnvironmentNameNotFound: $2" >&2; return 1; }
        export CONDA_DEFAULT_ENV="$2"
    fi
}
HOOK
fi
"""

FAKE_SINGULARITY = """#!/bin/bash
# singularity exec [--bind X] IMAGE CMD...: the image must be one argument
shift
[ "$1" = --bind ] && shift 2
case "$1" in
    *://*) ;;
    *) [ -e "$1" ] || { echo "FATAL: could not open image $1" >&2; exit 255; } ;;
esac
shift
"$@"
"""


@pytest.fixture
def job(tmp_path):
    """Run a job script with fake module/conda/singularity commands."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, text in (
        ("module", FAKE_MODULE),
        ("conda", FAKE_CONDA),
        ("singularity", FAKE_SINGULARITY),
    ):
        path = bin_dir / name
        path.write_text(text)
        path.chmod(0o755)
    image = tmp_path / "image.sif"
    image.write_text("")
    logs = tmp_path / "logs"

    def run(template, command="touch ran", **overrides):
        env = {
            "PATH": f"{bin_dir}:/usr/bin:/bin",
            "HOME": str(tmp_path),
            "cmd_b64": base64.b64encode(command.encode()).decode(),
            "cwd": str(tmp_path),
            "out": str(logs / "job.out"),
            "err": str(logs / "job.err"),
            "quiet": "false",
            "PBS_JOBID": "123.gadi-pbs",
            "mods": "python3",
            "env": "myenv",
            "sif": str(image),
            "QXUB_RUN_FILE": str(tmp_path / "runs" / "qx-1.run"),
        }
        env.update(overrides)
        result = subprocess.run(
            ["bash", str(JOBSCRIPTS / f"{template}.pbs")],
            env=env,
            capture_output=True,
            text=True,
            timeout=30,
            cwd=tmp_path,
        )
        err = logs / "job.err"
        result.err_file = err.read_text() if err.exists() else ""
        result.run_file = (tmp_path / "runs" / "qx-1.run").read_text()
        return result

    return run


@pytest.mark.parametrize("template", TEMPLATES)
def test_missing_execdir_says_so_where_callers_look(job, tmp_path, template):
    missing = tmp_path / "no such dir"
    result = job(template, cwd=str(missing), create_execdir="false")

    assert result.returncode == 1
    message = f"execution directory does not exist or can't be entered: {missing}"
    assert message in result.stderr
    assert message in result.err_file
    assert "--create-execdir" in result.err_file
    assert not (tmp_path / "ran").exists()
    assert not (missing / "ran").exists()
    assert "state=finished" in result.run_file
    assert "exit_code=1" in result.run_file


@pytest.mark.parametrize("template", TEMPLATES)
def test_execdir_with_spaces_works(job, tmp_path, template):
    workdir = tmp_path / "a dir with spaces"
    workdir.mkdir()
    result = job(template, cwd=str(workdir))

    assert result.returncode == 0, result.stdout + result.stderr
    assert (workdir / "ran").exists()


@pytest.mark.parametrize("template", TEMPLATES)
def test_create_execdir_with_spaces(job, tmp_path, template):
    workdir = tmp_path / "new dir"
    result = job(template, cwd=str(workdir), create_execdir="true")

    assert result.returncode == 0, result.stdout + result.stderr
    assert (workdir / "ran").exists()


@pytest.mark.parametrize("template", ["qmod", "qconda_mod"])
def test_module_that_wont_load_is_named(job, tmp_path, template):
    result = job(template, mods="python3 bad-module")

    assert result.returncode == 1
    assert "could not load module: bad-module" in result.err_file
    assert "module avail bad-module" in result.err_file
    assert not (tmp_path / "ran").exists()


@pytest.mark.parametrize("template", ["qconda", "qconda_mod"])
def test_conda_env_that_wont_activate_is_named(job, tmp_path, template):
    result = job(template, env="missing")

    assert result.returncode == 1
    assert "Failed to activate conda environment: missing" in result.err_file
    assert "Failed to activate conda environment: missing" in result.stderr
    assert not (tmp_path / "ran").exists()


def test_missing_singularity_image_is_named(job, tmp_path):
    result = job("qsing", sif=str(tmp_path / "nope.sif"))

    assert result.returncode == 1
    assert f"Singularity image not found: {tmp_path / 'nope.sif'}" in result.err_file
    assert not (tmp_path / "ran").exists()


def test_singularity_image_path_with_spaces_reaches_singularity(job, tmp_path):
    """#140 retest: the inner ``bash -c`` split the image path."""
    image = tmp_path / "img dir" / "my image.sif"
    image.parent.mkdir()
    image.write_text("")
    for quiet in ("false", "true"):
        result = job("qsing", sif=str(image), quiet=quiet)
        assert result.returncode == 0, result.stdout + result.stderr + result.err_file
        assert (tmp_path / "ran").exists()
        (tmp_path / "ran").unlink()


def test_singularity_uri_is_not_checked_as_a_file(job, tmp_path):
    result = job("qsing", sif="docker://ubuntu:22.04")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "ran").exists()


def test_quiet_mode_reports_on_stderr_only(job, tmp_path):
    missing = tmp_path / "gone"
    result = job("qdefault", cwd=str(missing), quiet="true")

    assert result.returncode == 1
    assert "execution directory does not exist" in result.stderr
    assert result.err_file == ""
