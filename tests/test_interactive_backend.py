"""``qxub interactive`` starts its session through the scheduler backend."""

import os
import stat
from pathlib import Path

from qxub import interactive_cli


def test_qxi_runs_qsub_i_with_its_session_script(sandbox):
    result = sandbox.run(
        "interactive", "-P", "ab12", "--mem", "4GB", "--runtime", "30m",
        "--volumes", "gdata/ab12",
    )  # fmt: skip
    assert result.returncode == 0, result.stdout + result.stderr

    (call,) = sandbox.pbs_calls("qsub")
    argv = call["argv"]
    assert argv[0] == "-I"
    assert argv[argv.index("-P") + 1] == "ab12"
    assert [argv[i + 1] for i, a in enumerate(argv) if a == "-l"] == [
        "mem=4GB",
        "walltime=00:30:00",
        "storage=gdata/ab12",
    ]
    assert argv[-2] == "--" and argv[-1].endswith(".sh")
    assert "-v" not in argv  # sessions get no job variables


def test_session_script_uses_the_backends_job_id_and_no_qstat():
    script = interactive_cli._build_interactive_script(
        shell="/bin/bash",
        working_dir="/tmp",
        pre_cmd=None,
        post_cmd=None,
        session_id="abc",
        job_id_variable="SLURM_JOB_ID",
    )
    assert 'echo "${SLURM_JOB_ID}" > "$QXUB_SESSION_DIR/job_id"' in script
    assert "qstat" not in script and "PBS_JOBID" not in script


def test_session_script_is_private_and_removed_afterwards(sandbox):
    result = sandbox.run("interactive", "-P", "ab12")
    assert result.returncode == 0, result.stdout + result.stderr

    (call,) = sandbox.pbs_calls("qsub")
    script = Path(call["argv"][-1])
    # A fresh directory of its own, not a fixed name shared by every user
    assert script.parent.name.startswith("qxub_interactive_")
    assert not script.parent.exists()


def test_write_session_script_is_owner_only(tmp_path, monkeypatch):
    import tempfile

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    # Another user's fixed-name directory is never used
    (tmp_path / "qxub_interactive").mkdir(mode=0o777)
    old = os.umask(0o022)
    try:
        directory, script = interactive_cli._write_session_script("echo hi\n")
        other, _ = interactive_cli._write_session_script("echo hi\n")
    finally:
        os.umask(old)

    assert directory != other and directory.parent == tmp_path
    assert directory.name != "qxub_interactive"
    assert stat.S_IMODE(directory.stat().st_mode) & 0o777 == 0o700
    assert stat.S_IMODE(script.stat().st_mode) & 0o777 == 0o700
    assert script.read_text() == "echo hi\n"

    interactive_cli._remove_session_script(directory)
    assert not directory.exists()
