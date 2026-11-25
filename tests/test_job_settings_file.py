"""Job settings go in a file the job sources, not on qsub's command line.

``qsub -v`` splits its list on commas and limits its length, and each
scheduler passes variables differently. qxub writes the job's variables to
``<runs dir>/<run ID>.env`` and passes the job only ``QXUB_JOB_FILE``.
"""

import stat
import subprocess
from pathlib import Path

import pytest
from conftest import REPO_ROOT, job_variables, parse_settings

from qxub.backends import JobSpec, SubmitError
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


def parse_pbs_v(value: str) -> dict:
    """A ``qsub -v`` list as PBS reads it: a quote anywhere opens a region
    that runs to the matching quote (the quotes are dropped), commas outside
    one separate items, and a backslash escapes the next character."""
    items, current, quote, chars = [], "", None, iter(value)
    for char in chars:
        if quote:
            if char == quote:
                quote = None
            else:
                current += char
        elif char in "'\"":
            quote = char
        elif char == "\\":
            current += next(chars, "")
        elif char == ",":
            items.append(current)
            current = ""
        else:
            current += char
    items.append(current)
    return dict(item.partition("=")[::2] for item in items)


HOSTILE = [
    "plain",
    "/a,/b",
    "/o'brien",
    'say "hi"',
    "x,QXUB_JOB_FILE=/tmp/evil",
    "x',QXUB_JOB_FILE=/tmp/evil,'",
    'x",QXUB_JOB_FILE=/tmp/evil,"',
    "/w/a,b'c\"d,cmd_b64=ZXZpbA==",
    "/w/a'b\"c",
    "ends with \\",
    "back\\,slash",
    "new\nline,x=1",
    "",
]


def test_qsub_v_values_round_trip_or_are_refused():
    for value in HOSTILE:
        env = {"cwd": value, "out": "/o", "QXUB_RUN_ID": "qx-1"}
        try:
            argv = qsub_argv(JobSpec(script="s", name="n", env=env))
        except SubmitError:
            continue
        assert parse_pbs_v(argv[argv.index("-v") + 1]) == env, value
    # The settings file's path gets the same care
    for value in filter(None, HOSTILE):
        try:
            argv = qsub_argv(JobSpec(script="s", name="n", env={}, env_file=value))
        except SubmitError:
            continue
        assert parse_pbs_v(argv[argv.index("-v") + 1]) == {"QXUB_JOB_FILE": value}


def test_qsub_v_refuses_what_it_cannot_carry():
    for value in ("a'b\"c", "a\\b", "a\nb"):
        with pytest.raises(SubmitError, match="job variable cwd"):
            qsub_argv(JobSpec(script="s", name="n", env={"cwd": value}))
    with pytest.raises(SubmitError, match="not a variable name"):
        qsub_argv(JobSpec(script="s", name="n", env={"a,b": "1"}))


def test_unwritable_runs_dir_refuses_the_job(monkeypatch, tmp_path):
    """The settings never fall back to qsub -v: other users may read them
    there (qstat -f), and its quoting is weaker than the file's."""
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("")
    monkeypatch.setattr(runs, "runs_dir", lambda: blocker / "runs")
    script = JOBSCRIPTS / "job.sh"  # reads QXUB_JOB_FILE
    spec = JobSpec(script=str(script), name="n", env={"a": "1"})

    with pytest.raises(SubmitError, match="Cannot write the job settings file"):
        runs.attach_settings(spec, "qx-1")
    monkeypatch.setattr(runs, "runs_dir", lambda: tmp_path / "runs")
    assert runs.attach_settings(spec, "qx-1").env_file  # writable: used


def test_exec_with_unwritable_runs_dir_submits_nothing(sandbox):
    runs_path = sandbox.home / ".local" / "state" / "qxub" / "status" / "runs"
    runs_path.parent.mkdir(parents=True)
    runs_path.write_text("")  # a file, so nothing can be written under it
    execdir = sandbox.cwd / "x,QXUB_JOB_FILE=evil,'\"y"
    execdir.mkdir()
    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--execdir", str(execdir),
        "--var", "TOKEN=s3cret", "--", "echo", "hi",
    )  # fmt: skip

    assert result.returncode != 0
    assert "Cannot write the job settings file" in result.stderr
    assert sandbox.pbs_calls("qsub") == []


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

    assert "doesn't read $QXUB_JOB_FILE" in result.stderr  # said, visibly
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


def test_only_scripts_that_source_the_file_get_it(tmp_path):
    cases = {
        '. "$QXUB_JOB_FILE"\n': True,
        'if ! . "$QXUB_JOB_FILE"; then exit 1; fi\n': True,
        "source ${QXUB_JOB_FILE}\n": True,
        "# reads QXUB_JOB_FILE? not yet\necho $cmd_b64\n": False,
        'echo "$QXUB_JOB_FILE"\n': False,
    }
    for i, (text, expected) in enumerate(cases.items()):
        script = tmp_path / f"t{i}.sh"
        script.write_text("#!/bin/bash\n" + text)
        assert runs.reads_settings_file(script) is expected, text
    assert runs.reads_settings_file(JOBSCRIPTS / "job.sh")
