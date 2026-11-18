"""``qxi`` records how a session ended (#125).

The generated init script used to ``exec`` the interactive shell, which
replaced the process holding the EXIT trap, so ``exit_code`` and
``end_time`` were never written. ``--post`` also set its own EXIT trap,
replacing the status one. After the session, qxub now also asks
``qstat -fx`` for the final exit status, which covers sessions killed at
walltime, where no trap runs.
"""

import subprocess
from pathlib import Path

import pytest

from qxub import interactive_cli
from qxub.config import site

FIXTURES = Path(__file__).parent / "fixtures"
SESSION = "abc123def456"


@pytest.fixture
def status_root(tmp_path, monkeypatch):
    """No platform, so status files go under ``$XDG_STATE_HOME``."""
    monkeypatch.setattr(site, "current_platform", lambda: None)
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "state"))
    return tmp_path / "state" / "qxub" / "status" / f"session_{SESSION}"


def run_session(tmp_path, typed: str, **options) -> subprocess.CompletedProcess:
    """Run the init script; ``typed`` is what the user types in the session.

    Piped stdin makes the inner bash read commands from it rather than a
    terminal, which is enough to drive the session from a test.
    """
    script = interactive_cli._build_interactive_script(
        shell="/bin/bash",
        working_dir=str(tmp_path),
        session_id=SESSION,
        **{"pre_cmd": None, "post_cmd": None, **options},
    )
    path = tmp_path / "qxi_init.sh"
    path.write_text(script)
    return subprocess.run(
        ["bash", str(path)],
        input=typed,
        capture_output=True,
        text=True,
        timeout=30,
        env={"PATH": "/usr/bin:/bin", "HOME": str(tmp_path), "PBS_JOBID": "9.pbs"},
    )


def test_exit_code_and_end_time_are_written(tmp_path, status_root):
    result = run_session(tmp_path, "exit 7\n")

    assert result.returncode == 7, result.stderr
    assert (status_root / "job_id").read_text().strip() == "9.pbs"
    assert (status_root / "exit_code").read_text().strip() == "7"
    assert (status_root / "end_time").exists()


def test_post_command_and_status_both_run(tmp_path, status_root):
    marker = tmp_path / "post_ran"

    result = run_session(tmp_path, "exit 3\n", post_cmd=f"touch {marker}")

    assert result.returncode == 3, result.stderr
    assert marker.exists()
    assert (status_root / "exit_code").read_text().strip() == "3"


def test_failed_pre_command_records_status_but_skips_post(tmp_path, status_root):
    marker = tmp_path / "post_ran"

    result = run_session(tmp_path, "", pre_cmd="false", post_cmd=f"touch {marker}")

    assert result.returncode == 1
    assert not marker.exists()
    assert (status_root / "exit_code").read_text().strip() == "1"


@pytest.mark.parametrize(
    "options",
    [{}, {"run_cmd": "R"}, {"container": "/c/tool.sif"}, {"conda_env": "dt3"}],
    ids=["shell", "run_cmd", "container", "conda"],
)
def test_the_session_never_replaces_the_script(tmp_path, options, status_root):
    script = interactive_cli._build_interactive_script(
        shell="/bin/bash",
        working_dir=str(tmp_path),
        pre_cmd=None,
        post_cmd="echo bye",
        session_id=SESSION,
        **options,
    )
    lines = [l.strip() for l in script.splitlines()]

    assert not [l for l in lines if l.startswith("exec ")]
    assert [l for l in lines if l.startswith("trap ")] == ["trap _qxub_on_exit EXIT"]


# -- after the session: qstat -fx ----------------------------------------------


@pytest.fixture
def fake_pbs(sandbox, monkeypatch):
    monkeypatch.setenv("PATH", sandbox.env["PATH"])
    monkeypatch.setenv("FAKE_PBS_STATE", str(sandbox.pbs_state))
    return sandbox


def finished_job(exit_status: int) -> str:
    text = (FIXTURES / "qstat_fx_finished.txt").read_text()
    return text.replace("Exit_status = 0", f"Exit_status = {exit_status}")


def test_session_killed_at_walltime_gets_exit_status_from_qstat(fake_pbs, status_root):
    # Killed at walltime: only the files written at start exist.
    status_root.mkdir(parents=True)
    (status_root / "job_id").write_text("1000.gadi-pbs\n")
    fake_pbs.set_qstat("1000.gadi-pbs", finished_job(-29))

    info = interactive_cli._process_session_completion(SESSION)

    assert info["job_id"] == "1000.gadi-pbs"
    assert info["exit_code"] == -29
    assert info["resources"]["exit_status"] == -29
    assert fake_pbs.pbs_calls("qstat")[0]["argv"] == ["-fx", "1000.gadi-pbs"]


def test_status_files_are_used_when_qstat_has_no_record(fake_pbs, status_root):
    status_root.mkdir(parents=True)
    (status_root / "job_id").write_text("424242.gadi-pbs\n")
    (status_root / "exit_code").write_text("5\n")

    info = interactive_cli._process_session_completion(SESSION)

    assert info["exit_code"] == 5
    assert "resources" not in info
    assert len(fake_pbs.pbs_calls("qstat")) == 1, "no retries for an unknown job"
