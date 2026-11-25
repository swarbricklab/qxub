"""Values pasted into the ``qxi`` session script never run as code.

The working directory, conda env, modules, container, bind mounts, record
directory and container prompt used to be pasted into the script bare or in
double quotes, so a directory name or config value holding ``$(...)``,
backticks or quotes ran commands as the user when the session started.
"""

import os
import stat
import subprocess
import textwrap

import pytest

from qxub import interactive_cli
from qxub.config import site

EVIL = (
    "x$(touch PWNED1)`touch PWNED2`\";touch PWNED3;\"';touch PWNED4;'"
    "\nRCEOF\ntouch PWNED5\n"
)


@pytest.fixture
def env(tmp_path, monkeypatch):
    """Stub module/conda/singularity/script that log their arguments."""
    monkeypatch.setattr(site, "current_platform", lambda: None)
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "state"))
    monkeypatch.setattr(interactive_cli, "_get_config_manager", lambda: None)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "calls.log"
    for name in ("module", "conda", "singularity", "script"):
        stub = bin_dir / name
        stub.write_text(
            textwrap.dedent(
                f"""\
                #!/bin/bash
                printf '%s\\0' {name} "$@" >> {log}
                printf '\\n' >> {log}
                """
            )
        )
        stub.chmod(0o755)
    return {
        "PATH": f"{bin_dir}:/usr/bin:/bin",
        "HOME": str(tmp_path),
        "PBS_JOBID": "9.pbs",
        "PROJECT": "ab12",
        "log": log,
    }


def run(tmp_path, env, **options):
    script = interactive_cli._build_interactive_script(
        shell="/bin/bash",
        pre_cmd=None,
        post_cmd=None,
        session_id="abc123",
        **{"working_dir": str(tmp_path), **options},
    )
    rc_tmp = tmp_path / "rctmp"
    rc_tmp.mkdir(exist_ok=True)
    script = script.replace("/tmp/qxub_rc.", f"{rc_tmp}/qxub_rc.")
    path = tmp_path / "init.sh"
    path.write_text(script)
    run_env = {k: v for k, v in env.items() if k != "log"}
    result = subprocess.run(
        ["bash", str(path)],
        input="exit 0\n",
        capture_output=True,
        text=True,
        timeout=30,
        cwd=tmp_path,
        env=run_env,
    )
    pwned = [p for p in tmp_path.rglob("PWNED*")]
    assert not pwned, (pwned, script)
    calls = []
    if env["log"].exists():
        for line in env["log"].read_text().split("\0\n"):
            if line:
                calls.append(line.split("\0"))
    return result, calls


def test_hostile_working_directory(tmp_path, env):
    evil_dir = tmp_path / EVIL.replace("\n", "_")
    evil_dir.mkdir()
    result, _ = run(tmp_path, env, working_dir=str(evil_dir))
    assert result.returncode == 0, result.stdout + result.stderr
    assert str(evil_dir) in result.stdout


def test_hostile_conda_env_and_modules(tmp_path, env):
    result, calls = run(tmp_path, env, conda_env=EVIL, modules=[EVIL, "ok/1.0"])
    assert ["module", "load", EVIL] in calls
    assert ["module", "load", "ok/1.0"] in calls
    assert ["conda", "activate", EVIL] in calls
    assert result.returncode == 0, result.stdout + result.stderr


def test_hostile_container_and_bind(tmp_path, env):
    container = f"/images/{EVIL}.sif"
    result, calls = run(tmp_path, env, container=container, bind=EVIL)
    assert ["singularity", "shell", "--bind", EVIL, container] in calls
    assert result.returncode == 0, result.stdout + result.stderr


def test_hostile_record_dir_and_command(tmp_path, env):
    record_dir = tmp_path / "rec" / EVIL.replace("\n", "_")
    result, calls = run(
        tmp_path, env, run_cmd="echo hi", record=True, record_dir=str(record_dir)
    )
    (call,) = [c for c in calls if c[0] == "script"]
    assert call[1:3] == ["-eqc", "echo hi"]
    assert call[3].startswith(f"-t{record_dir}/qxub-timing-")


def test_container_prompt_keeps_container_name_inert():
    prompt = interactive_cli._get_container_prompt("a$(id)`id`\\w")
    assert "$(" not in prompt and "`" not in prompt and "\\w" not in prompt


def test_record_dir_variables_are_expanded(tmp_path, monkeypatch):
    monkeypatch.setenv("QXUB_TEST_REC", "/scratch/x")
    script = interactive_cli._build_interactive_script(
        shell="/bin/bash",
        working_dir="/tmp",
        pre_cmd=None,
        post_cmd=None,
        run_cmd="R",
        record=True,
        record_dir="$QXUB_TEST_REC/rec",
    )
    assert "RECORD_DIR=/scratch/x/rec" in script
