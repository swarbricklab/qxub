"""Site conventions come from the platform definition, not from code.

Covers ``qxub/config/site.py`` and how ``qxub exec`` uses it: the default
project (``account_env``), the shared scratch area (``paths.scratch``),
extra scheduler directives (``default_directives``) and matching the
current host against each platform's ``host``.
"""

import re
import shlex
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import REPO_ROOT, job_variables

from qxub.config import site
from qxub.platforms import core as platform_core

JOBSCRIPTS = sorted((REPO_ROOT / "qxub" / "jobscripts").glob("*.pbs"))

PLATFORM_YAML = """\
platform:
  name: {name}
  type: pbs_pro
  host: {host}
  account_env: TEST_ACCOUNT
  paths:
    scratch: "{scratch}/{{project}}/{{user}}"
  default_directives:
    - "-W umask=0027"
  queues:
    - name: normal
      type: standard
"""


def write_platform(directory: Path, name: str, host: str, scratch: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{name}.yaml").write_text(
        PLATFORM_YAML.format(name=name, host=host, scratch=scratch)
    )


@pytest.fixture
def platforms(tmp_path, monkeypatch):
    """Two platform files in a search path, and a fresh platform loader."""
    platform_dir = tmp_path / "platforms"
    scratch = tmp_path / "scratch"
    write_platform(platform_dir, "alpha", "alpha-login", scratch)
    write_platform(platform_dir, "beta", "beta", scratch)

    monkeypatch.setenv("QXUB_PLATFORM_PATHS", str(platform_dir))
    monkeypatch.setenv("USER", "tester")
    monkeypatch.delenv("QXUB_PLATFORM_OVERRIDE", raising=False)
    monkeypatch.delenv("TEST_ACCOUNT", raising=False)
    monkeypatch.setattr(platform_core, "_platform_loader", None)
    return scratch


# -- platform detection --------------------------------------------------------


@pytest.mark.parametrize(
    "hostname, expected",
    [("alpha-login-02.example.org", "alpha"), ("BETA-07", "beta"), ("other", None)],
)
def test_platform_is_detected_from_its_host_field(
    platforms, monkeypatch, hostname, expected
):
    monkeypatch.setenv("HOSTNAME", hostname)
    assert platform_core.detect_platform() == expected


# -- site helpers ------------------------------------------------------------


def test_default_project_comes_from_account_env(platforms, monkeypatch):
    monkeypatch.setenv("HOSTNAME", "beta-1")
    assert site.default_project() is None

    monkeypatch.setenv("TEST_ACCOUNT", "ef56")
    assert site.default_project() == "ef56"


def test_scratch_and_status_dirs_follow_the_platform(platforms, monkeypatch):
    monkeypatch.setenv("HOSTNAME", "beta-1")
    monkeypatch.setenv("TEST_ACCOUNT", "ef56")

    assert site.scratch_dir() == platforms / "ef56" / "tester"
    assert (
        site.default_status_dir() == platforms / "ef56" / "tester" / "qxub" / "status"
    )


def test_scratch_needs_a_project_when_the_template_uses_one(platforms, monkeypatch):
    monkeypatch.setenv("HOSTNAME", "beta-1")
    assert site.scratch_dir() is None


def test_without_a_platform_status_files_go_under_xdg_state(tmp_path, monkeypatch):
    monkeypatch.setattr(site, "current_platform", lambda: None)
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "state"))

    assert site.default_project() is None
    assert site.scratch_dir() is None
    assert site.default_directives() == []
    assert site.default_status_dir() == tmp_path / "state" / "qxub" / "status"


# -- qxub exec ---------------------------------------------------------------


@pytest.fixture
def sandbox_platform(sandbox):
    """The sandbox with platform ``site`` matching its host name."""
    scratch = sandbox.root / "scratch"
    write_platform(sandbox.root / "platforms", "site", "qxub-test-host", scratch)
    sandbox.env["QXUB_PLATFORM_PATHS"] = str(sandbox.root / "platforms")
    return sandbox, scratch


def dry_run_argv(sandbox, *args, env=None):
    result = sandbox.run("exec", "--dry-run", *args, "--", "true", env=env)
    assert result.returncode == 0, result.stdout + result.stderr
    (line,) = [l for l in result.stdout.splitlines() if l.startswith("Dry run")]
    return shlex.split(line.split(": ", 1)[1])


def dry_run_variables(sandbox, *args, env=None):
    result = sandbox.run("exec", "--dry-run", *args, "--", "true", env=env)
    assert result.returncode == 0, result.stdout + result.stderr
    (line,) = [l for l in result.stdout.splitlines() if l.startswith("Dry run")]
    return job_variables(shlex.split(line.split(": ", 1)[1]), result.stdout)


def test_exec_uses_the_platform_project_and_directives(sandbox_platform):
    sandbox, _ = sandbox_platform

    argv = dry_run_argv(sandbox, env={"TEST_ACCOUNT": "ef56"})

    assert argv[argv.index("-P") + 1] == "ef56"
    assert "umask=0027" in argv


def test_explicit_project_wins_over_account_env(sandbox_platform):
    sandbox, _ = sandbox_platform

    argv = dry_run_argv(sandbox, "-P", "ab12", env={"TEST_ACCOUNT": "ef56"})

    assert argv[argv.index("-P") + 1] == "ab12"


def test_exec_without_any_project_omits_minus_p(sandbox):
    assert "-P" not in dry_run_argv(sandbox)


def test_default_output_goes_to_platform_scratch_when_it_exists(sandbox_platform):
    sandbox, scratch = sandbox_platform
    user_scratch = scratch / "ef56" / sandbox.env.get("USER", "")
    user_scratch.mkdir(parents=True)

    variables = dry_run_variables(sandbox, env={"TEST_ACCOUNT": "ef56"})
    assert variables["out"].startswith(f"{user_scratch}/qt/")


def test_default_output_falls_back_to_cwd(sandbox_platform):
    sandbox, _ = sandbox_platform

    variables = dry_run_variables(sandbox, env={"TEST_ACCOUNT": "ef56"})
    assert variables["out"].startswith(f"{sandbox.cwd}/qt/")


# -- job scripts ---------------------------------------------------------------


@pytest.mark.parametrize("script", JOBSCRIPTS, ids=[p.name for p in JOBSCRIPTS])
def test_job_scripts_have_no_site_specific_values(script):
    text = script.read_text()
    for value in ("a56", "/g/data", "#PBS -P", "#PBS -l storage", "umask"):
        assert value not in text, f"{script.name} contains {value!r}"


@pytest.mark.parametrize("script", JOBSCRIPTS, ids=[p.name for p in JOBSCRIPTS])
def test_job_scripts_never_open_the_job_database(script):
    """Jobs update their run file instead (docs/dev/sqlite-robustness.md)."""
    text = script.read_text()
    for value in ("python", "qxub.queue", "qxub.resources", "QXUB_SHARED_DB"):
        assert value not in text, f"{script.name} contains {value!r}"
    assert "qxub_run_update state=running" in text
    # The finished state is recorded from the EXIT trap (directly, or from
    # the function it calls).
    assert "qxub_run_update state=finished" in text
    assert re.search(r"^trap .+ EXIT$", text, re.MULTILINE)


# -- every submission path passes the submitting Python ------------------------


def test_job_env_passes_qxub_python():
    from qxub.execution.context import ExecutionContext

    env = ExecutionContext("default", None, "default").job_env(
        ("true",), {"out": "o", "err": "e", "execdir": "/w", "quiet": True}
    )

    assert env["qxub_python"] == sys.executable


def test_submit_job_api_passes_qxub_python(sandbox):
    script = (
        "from qxub.execution.submit import submit_job; "
        "submit_job('true', env='base', project='ab12')"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=sandbox.cwd,
        env=sandbox.env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    (call,) = sandbox.pbs_calls("qsub")
    assert job_variables(call["argv"])["qxub_python"] == sys.executable
