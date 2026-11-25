"""Backends describe jobs as JobInfo/JobUsage, not in qstat -fx's layout (#65).

Consolidation, ``resources show`` and the end of an interactive session
read these, so a new backend only has to fill them in.
"""

from datetime import datetime

import pytest
from conftest import REPO_ROOT

from qxub.backends import (
    ENDED_BY_CANCEL,
    ENDED_BY_WALLTIME,
    JobInfo,
    JobState,
    JobUsage,
    LocalBackend,
)
from qxub.backends.pbs import job_info_from_details, parse_qstat_fx_output

LOCAL_PLATFORM = REPO_ROOT / "docs" / "platforms" / "local" / "local.yaml"

FINISHED = """\
Job Id: 185107076.gadi-pbs
    Job_Name = align
    Job_Owner = abc123@gadi-login-01.gadi.nci.org.au
    resources_used.cpupercent = 380
    resources_used.cput = 03:10:00
    resources_used.jobfs = 0b
    resources_used.mem = 6291456kb
    resources_used.ncpus = 4
    resources_used.walltime = 00:50:00
    job_state = F
    queue = normal-exec
    exec_host = gadi-cpu-clx-0123/4*4
    Exit_status = {code}
    Resource_List.jobfs = 100mb
    Resource_List.mem = 16gb
    Resource_List.ncpus = 4
    Resource_List.ngpus = 0
    Resource_List.walltime = 01:00:00
    qtime = Tue Nov 24 09:00:00 2026
    stime = Tue Nov 24 09:10:00 2026
    mtime = Tue Nov 24 10:00:05 2026
    obittime = Tue Nov 24 10:00:00 2026
"""


def pbs(text):
    return job_info_from_details(
        "185107076.gadi-pbs", parse_qstat_fx_output(text.format(code=0))
    )


def test_pbs_job_in_neutral_terms():
    info = pbs(FINISHED)

    assert info.state is JobState.FINISHED and info.is_finished
    assert info.exit_status == 0 and info.ended_by is None
    assert (info.name, info.queue, info.host) == (
        "align",
        "normal-exec",
        "gadi-cpu-clx-0123",
    )
    assert info.submitted == datetime(2026, 11, 24, 9, 0)
    assert info.queue_wait == 600
    assert info.run_time == 50 * 60, "from obittime, not mtime"
    assert info.usage == JobUsage(
        cpus_requested=4,
        cpus_used=4,
        cpu_time=3 * 3600 + 600,
        cpu_percent=380.0,
        memory_requested=16 * 1024**3,
        memory_used=6 * 1024**3,
        walltime_requested=3600,
        walltime_used=3000,
        disk_requested=100 * 1024**2,
        disk_used=0,
        gpus_requested=0,
    )
    assert info.native["resources_used"]["mem"] == "6291456kb"


@pytest.mark.parametrize(
    "letter, state",
    [
        ("Q", JobState.QUEUED),
        ("W", JobState.QUEUED),
        ("R", JobState.RUNNING),
        ("E", JobState.RUNNING),
        ("H", JobState.HELD),
        ("S", JobState.HELD),
        ("F", JobState.FINISHED),
        ("X", JobState.FINISHED),
        ("?", JobState.UNKNOWN),
    ],
)
def test_pbs_states(letter, state):
    info = pbs(FINISHED.replace("job_state = F", f"job_state = {letter}"))
    assert info.state is state
    assert (info.finished is None) == (state is not JobState.FINISHED)


def test_pbs_walltime_kill_is_told_apart():
    """Killed at walltime: PBS's Exit_status -29; qdel'd: 271 (#170 review)."""
    walltime = job_info_from_details(
        "1", parse_qstat_fx_output(FINISHED.format(code=-29))
    )
    qdel = job_info_from_details("1", parse_qstat_fx_output(FINISHED.format(code=271)))

    assert walltime.ended_by == ENDED_BY_WALLTIME
    assert qdel.ended_by == ENDED_BY_CANCEL
    other = job_info_from_details("1", parse_qstat_fx_output(FINISHED.format(code=1)))
    assert other.ended_by is None


def test_usage_knows_when_nothing_was_used():
    assert not JobUsage(cpus_requested=4, memory_requested=1).has_used()
    assert JobUsage(cpu_time=0).has_used()
    assert JobInfo("1", JobState.QUEUED).queue_wait is None


# -- local jobs ---------------------------------------------------------------


def test_local_job_info(sandbox, monkeypatch):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    job_id = sandbox.run("exec", "--terse", "--name", "li", "--", "true").stdout.strip()
    monkeypatch.setenv("XDG_STATE_HOME", sandbox.env["XDG_STATE_HOME"])
    backend = LocalBackend()
    sandbox.run("monitor", "--quiet", job_id)

    info = backend.find_info(job_id)
    assert info.is_finished and info.exit_status == 0 and info.ended_by is None
    assert info.name == "li" and info.queue == "local"
    assert info.finished >= info.started
    assert backend.job_details([job_id])[job_id] == info.native


def test_cancelled_local_job_says_so(sandbox, monkeypatch):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    job_id = sandbox.run("exec", "--terse", "--", "sleep 30").stdout.strip()
    monkeypatch.setenv("XDG_STATE_HOME", sandbox.env["XDG_STATE_HOME"])

    assert sandbox.run("cancel", "id", "-f", job_id).returncode == 0
    sandbox.run("monitor", "--quiet", job_id)

    info = LocalBackend().find_info(job_id)
    assert info.is_finished and info.ended_by == ENDED_BY_CANCEL


def test_failed_cancel_leaves_no_marker(monkeypatch, tmp_path):
    import qxub.backends.local as local_mod

    monkeypatch.setattr(local_mod, "jobs_dir", lambda: tmp_path)
    job = tmp_path / "local-gone"
    job.mkdir()
    (job / "meta.json").write_text(
        '{"name": "j", "queue": "local", "job_log": "/x", "submitted": "2026-11-24T00:00:00"}'
    )
    (job / "pid").write_text("1\n")  # not ours to signal
    monkeypatch.setattr(
        local_mod.os, "killpg", lambda *a: (_ for _ in ()).throw(ProcessLookupError())
    )
    monkeypatch.setattr(local_mod, "_alive", lambda pid: True)

    ok, _ = LocalBackend().cancel("local-gone")
    assert not ok and not (job / "cancelled").exists()


# -- resources show -----------------------------------------------------------


def test_resources_show_reads_job_info(sandbox):
    sandbox.set_qstat("185107076.gadi-pbs", FINISHED.format(code=-29))

    result = sandbox.run("resources", "show", "185107076.gadi-pbs")

    assert result.returncode == 0, result.stderr
    out = result.stdout
    assert "State: finished (walltime)" in out
    assert "Exit Status: -29" in out
    assert "Host: gadi-cpu-clx-0123" in out
    assert "16.0GB" in out and "6.0GB" in out and "37.5%" in out
    assert "01:00:00" in out and "00:50:00" in out and "83.3%" in out
    assert "95.0%" in out  # 3h10m CPU time / (50 min x 4 CPUs)
    assert "Queue wait: 600 seconds" in out


def test_resources_show_a_local_job(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    job_id = sandbox.run("exec", "--terse", "--", "true").stdout.strip()
    sandbox.run("monitor", "--quiet", job_id)

    result = sandbox.run("resources", "show", job_id)

    assert result.returncode == 0, result.stderr
    assert "State: finished" in result.stdout and "Queue: local" in result.stdout


# -- how qx learns a PBS job ended ---------------------------------------------


def test_pbs_job_end_reads_the_job_log_first(tmp_path, monkeypatch):
    from qxub.backends import JobEnd, PbsBackend

    asked = []
    backend = PbsBackend()
    monkeypatch.setattr(backend, "job_details", lambda ids: asked.append(ids) or {})
    log = tmp_path / "j.log"

    assert backend.job_end("1.pbs", log) is None, "no log yet: running"
    log.write_text("output\nExit Status:        -29 (walltime)\n")
    assert backend.job_end("1.pbs", log) == JobEnd(-29, ENDED_BY_WALLTIME)
    assert asked == [], "Gadi: the log is enough"


def test_pbs_job_end_asks_pbs_when_the_log_has_no_exit_line(tmp_path, monkeypatch):
    """Plain PBS (OpenPBS, UNSW Katana) writes no Exit Status line."""
    from qxub.backends import JobEnd, PbsBackend

    log = tmp_path / "j.log"
    log.write_text("=>> PBS: job killed: walltime 68 exceeded limit 60\n")
    backend = PbsBackend()
    details = {"1.pbs": parse_qstat_fx_output(FINISHED.format(code=-29))}
    monkeypatch.setattr(backend, "job_details", lambda ids: details)

    assert backend.job_end("1.pbs", log) == JobEnd(-29, ENDED_BY_WALLTIME)

    exiting = FINISHED.format(code=0).replace("job_state = F", "job_state = E")
    details["1.pbs"] = parse_qstat_fx_output(exiting)
    assert backend.job_end("1.pbs", log) is None, "not finished yet"
