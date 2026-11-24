"""Status files live in QXUB_STATUS_DIR and are named after the run ID.

Job scripts used to derive their status directory from ``$out``
(``${out%/*}/status``) while ``qxub monitor`` and ``qxub status`` read
``site.default_status_dir()``; the two only agreed when the paths happened to
line up. Per-step files (``pre_exit_code`` …) also had no job in their name,
so jobs sharing an output directory overwrote each other's.
"""

import pytest
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
