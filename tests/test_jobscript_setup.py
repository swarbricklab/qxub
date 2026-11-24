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
JOB_SCRIPT = JOBSCRIPTS / "job.sh"
# Every execution context runs the same job script with different settings.
TEMPLATES = ["default", "module", "conda", "conda_module", "singularity"]

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
# singularity exec [--bind X] IMAGE CMD...: the image must be one argument.
# Records the bind list, and marks what runs "inside" the container.
shift
if [ "$1" = --bind ]; then
    printf '%s' "$2" > "${FAKE_SING_BIND:-/dev/null}"
    shift 2
fi
case "$1" in
    *://*) ;;
    *) [ -e "$1" ] || { echo "FATAL: could not open image $1" >&2; exit 255; } ;;
esac
shift
if [ -n "${FAKE_IMAGE_WITHOUT_BASH:-}" ]; then
    # An image with sh but no bash: only sh (and what it needs) on PATH
    mkdir -p "$HOME/nobash"
    for tool in sh cat; do ln -sf "$(command -v $tool)" "$HOME/nobash/$tool"; done
    INSIDE_CONTAINER=yes PATH="$HOME/nobash" "$@"
    exit
fi
INSIDE_CONTAINER=yes "$@"
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
            "QXUB_RUN_FILE": str(tmp_path / "runs" / "qx-1.run"),
        }
        # The settings qxub passes for each context
        if template in ("module", "conda_module"):
            env["mods"] = "python3"
        if template in ("conda", "conda_module"):
            env["env"] = "myenv"
        if template == "singularity":
            env["sif"] = str(image)
        env.update(overrides)
        result = subprocess.run(
            ["bash", str(JOB_SCRIPT)],
            env=env,
            capture_output=True,
            text=True,
            timeout=30,
            cwd=tmp_path,
        )
        err = logs / "job.err"
        result.err_file = err.read_text() if err.exists() else ""
        run_file = Path(env["QXUB_RUN_FILE"])
        result.run_file = run_file.read_text() if run_file.exists() else ""
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


@pytest.mark.parametrize("template", ["module", "conda_module"])
def test_module_that_wont_load_is_named(job, tmp_path, template):
    result = job(template, mods="python3 bad-module")

    assert result.returncode == 1
    assert "could not load module: bad-module" in result.err_file
    assert "module avail bad-module" in result.err_file
    assert not (tmp_path / "ran").exists()


@pytest.mark.parametrize("template", ["conda", "conda_module"])
def test_conda_env_that_wont_activate_is_named(job, tmp_path, template):
    result = job(template, env="missing")

    assert result.returncode == 1
    assert "Failed to activate conda environment: missing" in result.err_file
    assert "Failed to activate conda environment: missing" in result.stderr
    assert not (tmp_path / "ran").exists()


def test_missing_singularity_image_is_named(job, tmp_path):
    result = job("singularity", sif=str(tmp_path / "nope.sif"))

    assert result.returncode == 1
    assert f"Singularity image not found: {tmp_path / 'nope.sif'}" in result.err_file
    assert not (tmp_path / "ran").exists()


def test_singularity_image_path_with_spaces_reaches_singularity(job, tmp_path):
    """#140 retest: the inner ``bash -c`` split the image path."""
    image = tmp_path / "img dir" / "my image.sif"
    image.parent.mkdir()
    image.write_text("")
    for quiet in ("false", "true"):
        result = job("singularity", sif=str(image), quiet=quiet)
        assert result.returncode == 0, result.stdout + result.stderr + result.err_file
        assert (tmp_path / "ran").exists()
        (tmp_path / "ran").unlink()


def test_singularity_uri_is_not_checked_as_a_file(job, tmp_path):
    result = job("singularity", sif="docker://ubuntu:22.04")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "ran").exists()


def test_quiet_mode_reports_on_stderr_only(job, tmp_path):
    missing = tmp_path / "gone"
    result = job("default", cwd=str(missing), quiet="true")

    assert result.returncode == 1
    assert "execution directory does not exist" in result.stderr
    assert result.err_file == ""


# -- where each step's output goes (#150) --------------------------------------


def _b64(text):
    return base64.b64encode(text.encode()).decode()


@pytest.mark.parametrize("template", TEMPLATES)
def test_pre_and_post_keep_the_main_commands_output(job, tmp_path, template):
    """``$out``/``$err`` hold the main command's output; the rest goes to the joblog."""
    result = job(
        template,
        command="echo MAIN-OUT; echo MAIN-ERR >&2",
        pre_cmd_b64=_b64("echo PRE-OUT; echo PRE-ERR >&2"),
        post_cmd_b64=_b64("echo POST-OUT; echo POST-ERR >&2"),
    )

    out = (tmp_path / "logs" / "job.out").read_text()
    assert result.returncode == 0, result.stdout + result.stderr
    assert out == "MAIN-OUT\n"
    assert result.err_file == "MAIN-ERR\n"
    joblog = result.stdout + result.stderr
    for text in ("PRE-OUT", "PRE-ERR", "POST-OUT", "POST-ERR"):
        assert text in joblog


def test_user_vars_reach_the_main_command_with_post(job, tmp_path):
    job(
        "default",
        command='echo "FOO=$FOO"',
        user_vars_b64=_b64("FOO=1"),
        post_cmd_b64=_b64("echo POST"),
    )

    assert (tmp_path / "logs" / "job.out").read_text() == "FOO=1\n"


@pytest.mark.parametrize("template", TEMPLATES)
def test_out_and_err_start_empty_on_a_rerun(job, tmp_path, template):
    logs = tmp_path / "logs"
    logs.mkdir()
    (logs / "job.out").write_text("old output\n")
    (logs / "job.err").write_text("old error\n")

    result = job(template, command="echo new")

    assert (logs / "job.out").read_text() == "new\n"
    assert result.err_file == ""


# -- one job script: how the steps run (step 6b) -------------------------------


@pytest.mark.parametrize("template", TEMPLATES)
def test_what_pre_sets_up_reaches_the_command(job, tmp_path, template):
    (tmp_path / "sub").mkdir()
    result = job(
        template,
        'echo "$GREETING from ${PWD##*/}" > "$HOME/seen"',
        pre_cmd_b64=_b64("export GREETING=hello; cd sub"),
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "seen").read_text() == "hello from sub\n"


@pytest.mark.parametrize("template", TEMPLATES)
def test_steps_run_like_a_plain_shell(job, tmp_path, template):
    """A failure part-way through a step doesn't stop it."""
    result = job(
        template,
        "false; echo main-went-on",
        pre_cmd_b64=_b64("false; touch pre-went-on"),
        post_cmd_b64=_b64("false; touch post-went-on"),
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "logs" / "job.out").read_text() == "main-went-on\n"
    assert (tmp_path / "pre-went-on").exists()
    assert (tmp_path / "post-went-on").exists()


def test_the_whole_command_runs_inside_the_container(job, tmp_path):
    """Pipes and && used to be run by the host shell around singularity."""
    result = job(
        "singularity",
        'echo x | (read; echo "${INSIDE_CONTAINER:-host}") && echo "${INSIDE_CONTAINER:-host}"',
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "logs" / "job.out").read_text() == "yes\nyes\n"


def test_bind_list_reaches_singularity_as_one_argument(job, tmp_path):
    bound = tmp_path / "bound"
    result = job(
        "singularity",
        "true",
        bind="/g/data/ab12,/scratch/a dir",
        FAKE_SING_BIND=str(bound),
    )

    assert result.returncode == 0, result.stderr
    assert bound.read_text() == "/g/data/ab12,/scratch/a dir"


def test_image_without_bash_runs_the_command_with_sh(job, tmp_path):
    """#161 review: minimal images have /bin/sh but no bash."""
    result = job(
        "singularity",
        'echo "in ${INSIDE_CONTAINER:-host}" && echo a,b | cat',
        FAKE_IMAGE_WITHOUT_BASH="1",
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "logs" / "job.out").read_text() == "in yes\na,b\n"
