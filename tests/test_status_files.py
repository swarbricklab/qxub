"""Status files live in QXUB_STATUS_DIR and are named after the run ID.

Job scripts used to derive their status directory from ``$out``
(``${out%/*}/status``) while ``qxub monitor`` and ``qxub status`` read
``site.default_status_dir()``; the two only agreed when the paths happened to
line up. Per-step files (``pre_exit_code`` …) also had no job in their name,
so jobs sharing an output directory overwrote each other's.
"""

import pytest
from test_golden_pbs_submission import dry_run_qsub_argv
from test_jobscript_setup import TEMPLATES, _b64, job  # noqa: F401 (fixture)

from qxub.core.scheduler import (
    _format_completion_message,
    job_status_from_files,
    monitor_job_single_thread,
    status_file_key,
)
from qxub.queue import runs

PBS_ID = "123.gadi-pbs"


def _first_line(path):
    return path.read_text().splitlines()[0]


@pytest.mark.parametrize("template", TEMPLATES)
def test_status_files_use_the_run_id_and_the_given_dir(job, tmp_path, template):
    status = tmp_path / "status"
    result = job(
        template,
        QXUB_RUN_ID="qx-1",
        QXUB_STATUS_DIR=str(status),
        pre_cmd_b64=_b64("true"),
        post_cmd_b64=_b64("true"),
    )

    assert result.returncode == 0, result.stderr
    for name in ("job_started", "pre_exit_code", "main_exit_code", "post_exit_code"):
        assert (status / f"{name}_qx-1").exists(), name
    assert _first_line(status / "final_exit_code_qx-1") == "0"
    assert (status / f"native_{PBS_ID}").read_text().strip() == "qx-1"
    assert not (tmp_path / "logs" / "status").exists()


@pytest.mark.parametrize("template", TEMPLATES)
def test_jobs_from_older_qxub_keep_the_old_names(job, tmp_path, template):
    """Without QXUB_RUN_ID/QXUB_STATUS_DIR: next to $out, by job ID."""
    result = job(template)

    assert result.returncode == 0, result.stderr
    status = tmp_path / "logs" / "status"
    assert (status / f"job_started_{PBS_ID}").exists()
    assert _first_line(status / f"final_exit_code_{PBS_ID}") == "0"
    assert (status / f"main_exit_code_{PBS_ID}").exists()


def test_readers_find_a_job_by_its_pbs_id_or_run_id(job, tmp_path):
    status = tmp_path / "status"
    job("qdefault", "sh -c 'exit 3'", QXUB_RUN_ID="qx-1", QXUB_STATUS_DIR=str(status))

    assert status_file_key(PBS_ID, status) == "qx-1"
    assert status_file_key("qx-1", status) == "qx-1"
    assert status_file_key("999.gadi-pbs", status) == "999.gadi-pbs"
    assert job_status_from_files(PBS_ID, status) == ("F", 3)
    assert job_status_from_files("qx-1", status) == ("F", 3)
    assert job_status_from_files("999.gadi-pbs", status) == ("Q", None)


def test_jobs_sharing_a_dir_report_their_own_steps(job, tmp_path):
    status = tmp_path / "status"
    job("qdefault", "true", QXUB_RUN_ID="qx-1", QXUB_STATUS_DIR=str(status))
    job(
        "qdefault",
        "true",
        QXUB_RUN_ID="qx-2",
        QXUB_STATUS_DIR=str(status),
        PBS_JOBID="124.gadi-pbs",
        pre_cmd_b64=_b64("true"),
        post_cmd_b64=_b64("true"),
    )

    first = _format_completion_message(PBS_ID, 0, status)
    second = _format_completion_message("124.gadi-pbs", 0, status)
    assert first.endswith("(main ✅)")
    assert second.endswith("(pre ✅, main ✅, post ✅)")


def test_exec_monitor_reads_the_run_id_files(tmp_path):
    status = tmp_path / "status"
    status.mkdir()
    (status / "final_exit_code_qx-9").write_text("4\n2026-01-01T00:00:00\n")

    exit_status = monitor_job_single_thread(
        PBS_ID,
        tmp_path / "job.out",
        tmp_path / "job.err",
        quiet=True,
        run_id="qx-9",
        status_dir=status,
    )
    assert exit_status == 4


def test_job_env_names_the_status_dir(monkeypatch, tmp_path):
    from qxub.config import site

    monkeypatch.setattr(site, "default_status_dir", lambda: tmp_path / "st")
    env = runs.job_env("qx-1", tmp_path / "st" / "runs" / "qx-1.run")
    assert env["QXUB_STATUS_DIR"] == str(tmp_path / "st")
    assert env["QXUB_RUN_ID"] == "qx-1"


# -- jobs charged to another project (review of #154) -------------------------


def _storage(argv):
    values = [argv[i + 1] for i, a in enumerate(argv) if a == "-l"]
    return [v for v in values if v.startswith("storage=")]


@pytest.mark.parametrize(
    "args, expected",
    [
        (["-P", "px14"], ["storage=scratch/a56"]),
        (["-P", "px14", "--volumes", "gdata/px14"], ["storage=gdata/px14+scratch/a56"]),
        (["-P", "px14", "--volumes", "scratch/a56"], ["storage=scratch/a56"]),
        (["-P", "a56"], []),  # the default project: its scratch is mounted
    ],
)
def test_another_project_gets_the_status_dirs_scratch(sandbox, args, expected):
    """Gadi mounts only the job's own scratch; qxub's status dir is under $PROJECT's."""
    sandbox.use_gadi_platform()
    result = sandbox.run(
        "exec", "--dry-run", *args, "--", "true", env={"PROJECT": "a56"}
    )
    assert _storage(dry_run_qsub_argv(result)) == expected


def test_unreachable_status_dir_does_not_fail_the_job(job, tmp_path):
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("")
    result = job(
        "qdefault",
        QXUB_RUN_ID="qx-1",
        QXUB_STATUS_DIR=str(blocker / "status"),
        QXUB_RUN_FILE=str(blocker / "runs" / "qx-1.run"),
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "ran").exists()
    assert "Cannot create the status directory" in result.stderr
    assert _first_line(tmp_path / "logs" / "status" / "final_exit_code_qx-1") == "0"


@pytest.mark.parametrize("pbs_exit", [1, 0])
def test_exec_monitor_stops_when_the_joblog_ends_first(tmp_path, pbs_exit):
    """No status file ever appears: the joblog's exit status ends the wait."""
    joblog = tmp_path / "job.log"
    joblog.write_text(f"=====\nExit Status:        {pbs_exit}\n")

    exit_status = monitor_job_single_thread(
        PBS_ID,
        tmp_path / "job.out",
        tmp_path / "job.err",
        quiet=True,
        joblog_file=joblog,
        joblog_check_interval=0,
        run_id="qx-9",
        status_dir=tmp_path / "status",
    )
    assert exit_status == pbs_exit


# -- a command that calls exit ------------------------------------------------


@pytest.mark.parametrize("template", TEMPLATES)
def test_exit_in_the_command_is_recorded(job, tmp_path, template):
    """``exit N`` ends the command, not the job script before it records N."""
    status = tmp_path / "status"
    result = job(template, "exit 3", QXUB_RUN_ID="qx-1", QXUB_STATUS_DIR=str(status))

    assert result.returncode == 3
    assert _first_line(status / "final_exit_code_qx-1") == "3"
    assert "exit_code=3" in result.run_file


@pytest.mark.parametrize("template", TEMPLATES)
def test_exit_0_in_the_command_still_runs_post(job, tmp_path, template):
    status = tmp_path / "status"
    result = job(
        template,
        "exit 0",
        QXUB_RUN_ID="qx-1",
        QXUB_STATUS_DIR=str(status),
        post_cmd_b64=_b64("touch post_ran"),
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "post_ran").exists()
    assert _first_line(status / "final_exit_code_qx-1") == "0"


@pytest.mark.parametrize("step", ["pre", "post"])
def test_exit_in_pre_or_post_is_recorded(job, tmp_path, step):
    """qdefault evals --pre/--post in the job's shell, so exit ends the job."""
    status = tmp_path / "status"
    result = job(
        "qdefault",
        "true",
        QXUB_RUN_ID="qx-1",
        QXUB_STATUS_DIR=str(status),
        **{f"{step}_cmd_b64": _b64("exit 4")},
    )

    assert result.returncode == 4
    assert _first_line(status / "final_exit_code_qx-1") == "4"
