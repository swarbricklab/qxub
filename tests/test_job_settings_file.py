"""Job settings go in a file the job sources, not on qsub's command line.

``qsub -v`` splits its list on commas and limits its length, and each
scheduler passes variables differently. qxub writes the job's variables to
``<runs dir>/<run ID>.env`` and passes the job only ``QXUB_JOB_FILE``.
"""

import stat
import subprocess
from pathlib import Path

from conftest import REPO_ROOT, job_variables, parse_settings

from qxub.backends import JobSpec
from qxub.backends.pbs import qsub_argv
from qxub.queue import runs

JOBSCRIPTS = REPO_ROOT / "qxub" / "jobscripts"
AWKWARD = "a,b 'quoted' \"dq\" $HOME `x` ; c\nsecond line"


def test_settings_round_trip_through_bash(tmp_path):
    env = {"plain": "x", "awkward": AWKWARD, "empty": ""}
    path = tmp_path / "s.env"
    path.write_text(runs.format_settings(env))

    script = '. "$1"; printf "%s\\0" "$plain" "$awkward" "$empty"'
    out = subprocess.run(
        ["bash", "-c", script, "bash", str(path)], capture_output=True, text=True
    )
    assert out.stdout.split("\0")[:3] == ["x", AWKWARD, ""]
    assert parse_settings(path.read_text()) == env


def test_qsub_gets_only_the_settings_path():
    spec = JobSpec(script="s", name="n", env={"a": "1,2"}, env_file="/r/qx-1.env")
    assert qsub_argv(spec) == [
        "qsub",
        "-v",
        "QXUB_JOB_FILE=/r/qx-1.env",
        "-N",
        "n",
        "s",
    ]
    # Without a file, the variables go on the command line as before.
    spec = JobSpec(script="s", name="n", env={"a": "1,2"})
    assert qsub_argv(spec)[:3] == ["qsub", "-v", "a='1,2'"]


def test_unwritable_runs_dir_falls_back_to_qsub_v(monkeypatch, tmp_path):
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("")
    monkeypatch.setattr(runs, "runs_dir", lambda: blocker / "runs")
    script = JOBSCRIPTS / "qdefault.pbs"  # reads QXUB_JOB_FILE
    spec = JobSpec(script=str(script), name="n", env={"a": "1"})

    assert runs.attach_settings(spec, "qx-1") is spec
    monkeypatch.setattr(runs, "runs_dir", lambda: tmp_path / "runs")
    assert runs.attach_settings(spec, "qx-1").env_file  # writable: used


def test_exec_writes_the_settings_file(sandbox):
    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--sif", "/c/t.sif",
        "--bind", "/a,/b", "--", "echo", "hi, there",
    )  # fmt: skip
    assert result.returncode == 0, result.stdout + result.stderr

    (call,) = sandbox.pbs_calls("qsub")
    argv = call["argv"]
    (setting,) = [argv[i + 1] for i, a in enumerate(argv) if a == "-v"]
    name, _, path = setting.partition("=")
    assert name == "QXUB_JOB_FILE"
    path = Path(path)
    assert path.name == job_variables(argv)["QXUB_RUN_ID"] + ".env"
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert job_variables(argv)["bind"] == "/a,/b"


def _run_template(tmp_path, template, env_file, **env):
    logs = tmp_path / "logs"
    base = {"PATH": "/usr/bin:/bin", "HOME": str(tmp_path), "PBS_JOBID": "1.x"}
    return (
        subprocess.run(
            ["bash", str(JOBSCRIPTS / "job.sh")],
            env={**base, "QXUB_JOB_FILE": str(env_file), **env},
            capture_output=True,
            text=True,
            timeout=30,
            cwd=tmp_path,
        ),
        logs,
    )


def test_job_reads_everything_from_the_settings_file(tmp_path):
    import base64

    logs = tmp_path / "logs"
    settings = {
        "cmd_b64": base64.b64encode(b'printf "%s" "$awkward" > seen').decode(),
        "cwd": str(tmp_path),
        "out": str(logs / "job.out"),
        "err": str(logs / "job.err"),
        "quiet": "false",
        "awkward": AWKWARD,
        "QXUB_RUN_ID": "qx-1",
        "QXUB_RUN_FILE": str(tmp_path / "runs" / "qx-1.run"),
        "QXUB_STATUS_DIR": str(tmp_path / "status"),
    }
    env_file = tmp_path / "qx-1.env"
    env_file.write_text(runs.format_settings(settings))

    result, _ = _run_template(tmp_path, "default", env_file)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "seen").read_text() == AWKWARD
    assert (tmp_path / "status" / "final_exit_code_qx-1").exists()


def test_missing_settings_file_fails_the_job_clearly(tmp_path):
    result, _ = _run_template(tmp_path, "default", tmp_path / "gone.env")

    assert result.returncode == 1
    assert "Cannot read the job settings file" in result.stderr


def test_custom_template_without_the_file_still_gets_qsub_v(sandbox):
    """A --template written for older qxub never loads QXUB_JOB_FILE."""
    template = sandbox.cwd / "old.pbs"
    template.write_text('#!/bin/bash\necho "$cmd_b64" | base64 -d | bash\n')
    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--template", str(template),
        "--", "echo", "a,b",
    )  # fmt: skip
    assert result.returncode == 0, result.stdout + result.stderr

    (call,) = sandbox.pbs_calls("qsub")
    variables = job_variables(call["argv"])
    assert "QXUB_JOB_FILE" not in variables
    assert variables["cmd_b64"] and variables["QXUB_RUN_ID"]


def test_custom_template_that_reads_the_file_gets_it(sandbox):
    template = sandbox.cwd / "new.pbs"
    template.write_text('#!/bin/bash\n. "$QXUB_JOB_FILE"\neval "$cmd"\n')
    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--template", str(template), "--", "true"
    )
    assert result.returncode == 0, result.stdout + result.stderr

    (call,) = sandbox.pbs_calls("qsub")
    assert "QXUB_JOB_FILE" in job_variables(call["argv"])
