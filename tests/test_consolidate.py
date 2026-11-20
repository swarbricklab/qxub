"""Consolidating run files into the published stats snapshot (Phase 3).

See ``docs/dev/sqlite-robustness.md`` and ``qxub/queue/consolidate.py``.
"""

import os
import sqlite3
import subprocess
import sys
import textwrap
import time
from datetime import datetime, timedelta
from pathlib import Path

import pytest

from qxub.queue import consolidate as cons
from qxub.queue import runs

NOW = datetime(2026, 11, 20, 12, 0, 0)

JOBLOG = """\
some job output
======================================================================================
                  Resource Usage on {when}:
   Job Id:             {job_id}
   Project:            ab12
   Exit Status:        {code}
   Service Units:      0.50
   NCPUs Requested:    4                      NCPUs Used: 4
                                           CPU Time Used: 00:20:00
   Memory Requested:   8.0GB                 Memory Used: 2.0GB
   Walltime requested: 02:00:00            Walltime Used: 00:10:00
   JobFS requested:    100.0MB                JobFS used: 0B
======================================================================================
"""

QSTAT = """\
Job Id: {job_id}
    Job_Name = thing
    job_state = {state}
    queue = normal
    Resource_List.mem = 16gb
    Resource_List.ncpus = 2
    Resource_List.walltime = 01:00:00
    resources_used.mem = 4gb
    resources_used.ncpus = 2
    resources_used.walltime = 00:30:00
    resources_used.cpupercent = 150
    resources_used.cput = 00:45:00
    Exit_status = {code}
"""


@pytest.fixture
def area(tmp_path):
    """A runs directory and where the snapshot goes."""
    runs_dir = tmp_path / "runs"
    runs_dir.mkdir()
    return runs_dir, tmp_path / "shared" / "qxub-stats.db"


def write_run(runs_dir, run_id, **fields):
    base = {
        "run_id": run_id,
        "state": "submitted",
        "submitted_at": (NOW - timedelta(hours=3)).isoformat(timespec="seconds"),
        "user": "alice",
        "command_b64": "ZWNobyBoaQ==",  # echo hi
        "queue": "normal",
        "tags": '["workflow=w","rule=r"]',
    }
    base.update(fields)
    path = runs_dir / f"{run_id}.run"
    path.write_text(runs.format_fields(base))
    return path


def finished(runs_dir, run_id, job_id, code=0, ago=timedelta(hours=1), **fields):
    when = (NOW - ago).isoformat(timespec="seconds")
    return write_run(
        runs_dir,
        run_id,
        native_job_id=job_id,
        started_at=when,
        finished_at=when,
        host="gadi-cpu-1",
        exit_code=code,
        state="finished",
        **fields,
    )


def no_qstat(job_ids):
    raise AssertionError(f"qstat should not be called: {job_ids}")


def run(area, **kwargs):
    runs_dir, snapshot = area
    kwargs.setdefault("qstat", no_qstat)
    kwargs.setdefault("now", NOW)
    return cons.consolidate(runs_dir, snapshot, **kwargs)


def rows(snapshot, sql="SELECT * FROM job_resources ORDER BY job_id"):
    conn = cons.open_snapshot(snapshot)
    try:
        return [dict(r) for r in conn.execute(sql)]
    finally:
        conn.close()


# -- what goes in ---------------------------------------------------------------------


def test_finished_run_with_joblog_is_published_with_usage(area, tmp_path):
    runs_dir, snapshot = area
    joblog = tmp_path / "job.log"
    joblog.write_text(JOBLOG.format(when=NOW, job_id="1.gadi-pbs", code=3))
    finished(runs_dir, "qx-a", "1.gadi-pbs", code=3, joblog=joblog)

    result = run(area)

    assert result.published and result.added == 1, result.summary()
    (job,) = rows(snapshot)
    assert job["job_id"] == "1.gadi-pbs"
    assert job["status"] == "failed" and job["exit_code"] == 3
    assert job["command"] == "echo hi"
    assert job["username"] == "alice"
    assert job["tags"] == '["workflow=w", "rule=r"]'
    assert job["mem_used_mb"] == 2048 and job["mem_efficiency"] == 25.0
    assert job["time_efficiency"] == pytest.approx(8.3)
    assert job["cpu_efficiency"] == 50.0  # 20 min CPU / (10 min x 4 cpus)
    assert (runs_dir / "qx-a.run").exists(), "kept for the grace period"


def test_usage_comes_from_one_batched_qstat_when_there_is_no_joblog(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs")
    finished(runs_dir, "qx-b", "2.gadi-pbs")
    calls = []

    def qstat(job_ids):
        calls.append(list(job_ids))
        return {
            j: cons._qstat_parse(QSTAT.format(job_id=j, state="F", code=0))
            for j in job_ids
        }

    result = run(area, qstat=qstat)

    assert calls == [["1.gadi-pbs", "2.gadi-pbs"]]
    assert result.added == 2
    assert [r["mem_used_mb"] for r in rows(snapshot)] == [4096, 4096]
    assert rows(snapshot)[0]["cpu_percent"] == 150.0
    # Review of #135: the same measure as from a joblog, not cpupercent.
    assert rows(snapshot)[0]["cpu_efficiency"] == 75.0  # 45 min / (30 min x 2)


def test_query_qstat_parses_several_jobs_from_one_call(sandbox, monkeypatch):
    for job_id in ("1.gadi-pbs", "2.gadi-pbs"):
        sandbox.set_qstat(job_id, QSTAT.format(job_id=job_id, state="F", code=0))
    monkeypatch.setenv("PATH", sandbox.env["PATH"])
    monkeypatch.setenv("FAKE_PBS_STATE", sandbox.env["FAKE_PBS_STATE"])

    found = cons.query_qstat(["1.gadi-pbs", "2.gadi-pbs", "3.gadi-pbs"])

    assert sorted(found) == ["1.gadi-pbs", "2.gadi-pbs"]
    assert len(sandbox.pbs_calls("qstat")) == 1


def test_recent_run_without_usage_waits_and_old_one_goes_in_without(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-new", "1.gadi-pbs", ago=timedelta(hours=1))
    finished(runs_dir, "qx-old", "2.gadi-pbs", ago=timedelta(days=3))

    result = run(area, qstat=lambda ids: {})

    assert (result.added, result.waiting) == (1, 1)
    (job,) = rows(snapshot)
    assert job["job_id"] == "2.gadi-pbs" and job["mem_used_mb"] is None


def test_joblog_of_another_run_with_the_same_name_is_ignored(area, tmp_path):
    runs_dir, snapshot = area
    joblog = tmp_path / "job.log"
    joblog.write_text(JOBLOG.format(when=NOW, job_id="99.gadi-pbs", code=0))
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3), joblog=joblog)

    run(area, qstat=lambda ids: {})

    assert rows(snapshot)[0]["mem_used_mb"] is None


# -- runs that never recorded their end -----------------------------------------------------


def test_killed_job_is_settled_from_its_joblog(area, tmp_path):
    """SIGKILL skips the EXIT trap; the joblog still has the exit status."""
    runs_dir, snapshot = area
    joblog = tmp_path / "job.log"
    joblog.write_text(JOBLOG.format(when=NOW, job_id="1.gadi-pbs", code=137))
    path = write_run(
        runs_dir, "qx-a", state="running", native_job_id="1.gadi-pbs", joblog=joblog
    )

    result = run(area)

    assert result.settled == 1 and result.added == 1
    record = runs.read_run(path)
    assert record["state"] == "finished" and record["exit_code"] == "137"
    assert record["settled_by"] == "consolidate"
    assert rows(snapshot)[0]["status"] == "failed"


def test_job_pbs_reports_finished_is_settled_from_qstat(area):
    runs_dir, snapshot = area
    path = write_run(runs_dir, "qx-a", state="running", native_job_id="1.gadi-pbs")

    run(
        area,
        qstat=lambda ids: {
            "1.gadi-pbs": cons._qstat_parse(
                QSTAT.format(job_id="1.gadi-pbs", state="F", code=271)
            )
        },
    )

    assert runs.exit_code(runs.read_run(path)) == 271
    assert rows(snapshot)[0]["exit_code"] == 271


def test_runs_in_progress_are_left_alone(area):
    runs_dir, snapshot = area
    queued = write_run(runs_dir, "qx-q")
    running = write_run(runs_dir, "qx-r", state="running", native_job_id="1.gadi-pbs")
    before = {p: p.read_text() for p in (queued, running)}

    result = run(
        area,
        qstat=lambda ids: {
            "1.gadi-pbs": cons._qstat_parse(
                QSTAT.format(job_id="1.gadi-pbs", state="R", code="")
            )
        },
        grace_days=0,
    )

    assert result.in_progress == 2 and not result.published
    assert {p: p.read_text() for p in before} == before
    assert not snapshot.exists()


def test_run_that_never_started_is_abandoned_and_later_removed(area):
    runs_dir, snapshot = area
    path = write_run(
        runs_dir,
        "qx-a",
        submitted_at=(NOW - timedelta(days=40)).isoformat(timespec="seconds"),
    )

    result = run(area)
    assert result.abandoned == 1 and result.published
    assert rows(snapshot) == []
    (entry,) = rows(snapshot, "SELECT * FROM consolidated_runs")
    assert entry["status"] == "abandoned"

    run(area, now=NOW + timedelta(days=cons.GRACE_DAYS + 1))
    assert not path.exists()


# -- publishing and cleanup ---------------------------------------------------------------


def test_rounds_are_idempotent_and_unchanged_snapshots_are_not_republished(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    run(area, qstat=lambda ids: {})
    inode = snapshot.stat().st_ino

    again = run(area, qstat=lambda ids: {})

    assert not again.published and again.added == 0
    assert snapshot.stat().st_ino == inode
    assert len(rows(snapshot)) == 1


def test_new_runs_are_added_to_what_is_already_published(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    run(area, qstat=lambda ids: {})
    finished(runs_dir, "qx-b", "2.gadi-pbs", ago=timedelta(days=3))

    run(area, qstat=lambda ids: {})

    assert [r["job_id"] for r in rows(snapshot)] == ["1.gadi-pbs", "2.gadi-pbs"]


def test_run_files_are_removed_only_after_the_grace_period(area):
    runs_dir, snapshot = area
    path = finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    run(area, qstat=lambda ids: {})

    run(area, now=NOW + timedelta(days=cons.GRACE_DAYS - 1))
    assert path.exists()

    result = run(area, now=NOW + timedelta(days=cons.GRACE_DAYS + 1))
    assert result.removed == 1 and not path.exists()
    assert len(rows(snapshot)) == 1, "the snapshot keeps the job"


def test_dry_run_changes_nothing(area, tmp_path):
    runs_dir, snapshot = area
    joblog = tmp_path / "job.log"
    joblog.write_text(JOBLOG.format(when=NOW, job_id="1.gadi-pbs", code=9))
    path = write_run(
        runs_dir, "qx-a", state="running", native_job_id="1.gadi-pbs", joblog=joblog
    )
    before = path.read_text()

    result = run(area, dry_run=True, grace_days=0)

    assert result.added == 1 and not result.published
    assert not snapshot.exists() and path.read_text() == before


def test_published_file_is_private_copy_renamed(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    run(area, qstat=lambda ids: {})

    names = sorted(p.name for p in snapshot.parent.iterdir())
    assert names == ["qxub-stats.db"], "no journal, temp copy or lease left behind"
    conn = sqlite3.connect(snapshot)
    assert conn.execute("PRAGMA integrity_check").fetchone() == ("ok",)
    assert conn.execute("PRAGMA journal_mode").fetchone() == ("delete",)


def test_open_snapshot_is_read_only(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    run(area, qstat=lambda ids: {})

    conn = cons.open_snapshot(snapshot)
    with pytest.raises(sqlite3.OperationalError, match="readonly"):
        conn.execute("DELETE FROM job_resources")


# -- the lease -------------------------------------------------------------------------------


def test_round_skips_while_another_holds_the_lease(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    snapshot.parent.mkdir()
    other = cons.Lease(cons.lease_path(snapshot))
    assert other.acquire()

    result = run(area, qstat=lambda ids: {})

    assert result.skipped and "another consolidation" in result.skipped
    assert not snapshot.exists()
    other.release()
    assert not cons.lease_path(snapshot).exists()


def test_stale_lease_is_broken(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    snapshot.parent.mkdir()
    dead = cons.Lease(cons.lease_path(snapshot))
    assert dead.acquire()
    old = time.time() - cons.LEASE_STALE_SECONDS - 60
    os.utime(cons.lease_path(snapshot), (old, old))

    result = run(area, qstat=lambda ids: {})

    assert result.published
    assert [p.name for p in snapshot.parent.iterdir()] == ["qxub-stats.db"]


def test_round_that_lost_its_lease_does_not_publish(area):
    """Review of #135: a round past LEASE_STALE_SECONDS whose lease was taken."""
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    lease = cons.lease_path(snapshot)

    def qstat_while_lease_is_taken_over(job_ids):
        os.rename(lease, lease.with_name("gone"))
        assert cons.Lease(lease).acquire()
        return {}

    result = run(area, qstat=qstat_while_lease_is_taken_over)

    assert not result.published and "lost the lease" in result.errors[0]
    assert not snapshot.exists()
    assert lease.is_dir(), "the other round's lease is left alone"


def test_unwritable_snapshot_directory_is_reported(area):
    """Review of #135: other users couldn't take the lease, silently."""
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))
    snapshot.parent.mkdir()
    snapshot.parent.chmod(0o555)
    try:
        result = run(area, qstat=lambda ids: {})
    finally:
        snapshot.parent.chmod(0o755)

    assert "chmod g+ws" in result.skipped
    assert cons.writable_problem(snapshot) is None


def test_snapshot_directory_is_created_group_writable(area):
    runs_dir, snapshot = area
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=3))

    run(area, qstat=lambda ids: {})

    assert snapshot.parent.stat().st_mode & 0o070 == 0o070
    assert snapshot.stat().st_mode & 0o040


def test_default_snapshot_has_its_own_directory(tmp_path, monkeypatch):
    monkeypatch.delenv("QXUB_STATS_DB", raising=False)
    monkeypatch.setenv("QXUB_DB_PATH", str(tmp_path / "config" / "qxub.db"))
    # A site config (on Gadi, shared_db.path) must not leak in.
    monkeypatch.setattr(cons, "_config_value", lambda key: None)
    monkeypatch.setattr("qxub.queue.db.legacy_shared_db_path", lambda: None)

    assert cons.stats_db_path() == tmp_path / "config" / "stats" / "qxub-stats.db"


WORKER = textwrap.dedent(
    """
    import sys
    from qxub.queue import consolidate as cons
    result = cons.consolidate(sys.argv[1], sys.argv[2], qstat=lambda ids: {})
    print("published" if result.published else (result.skipped or "unchanged"))
    """
)


def test_concurrent_rounds_from_many_users_lose_nothing(tmp_path):
    snapshot = tmp_path / "shared" / "qxub-stats.db"
    users = []
    for u in range(6):
        runs_dir = tmp_path / f"user{u}" / "runs"
        runs_dir.mkdir(parents=True)
        for i in range(20):
            finished(
                runs_dir, f"qx-{u}-{i}", f"{u}{i:03d}.gadi-pbs", ago=timedelta(days=3)
            )
        users.append(runs_dir)

    def start(runs_dir):
        return subprocess.Popen(
            [sys.executable, "-c", WORKER, str(runs_dir), str(snapshot)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    procs = [start(d) for d in users]
    outputs = [p.communicate(timeout=120) for p in procs]
    assert [p.returncode for p in procs] == [0] * len(procs), outputs
    # Rounds that found the lease taken skipped; a second pass catches up.
    for runs_dir in users:
        start(runs_dir).communicate(timeout=120)

    conn = sqlite3.connect(snapshot)
    assert conn.execute("PRAGMA integrity_check").fetchone() == ("ok",)
    assert conn.execute("SELECT COUNT(*) FROM job_resources").fetchone() == (120,)


# -- when it runs -------------------------------------------------------------------------------


def test_first_submission_starts_the_clock_instead_of_a_round(tmp_path, monkeypatch):
    monkeypatch.delenv("QXUB_AUTO_CONSOLIDATE", raising=False)
    started = []
    monkeypatch.setattr(
        cons.subprocess, "Popen", lambda *a, **k: started.append((a, k))
    )

    assert not cons.maybe_start_in_background(tmp_path)
    assert (tmp_path / cons.ATTEMPT_MARKER).exists()
    assert not cons.maybe_start_in_background(tmp_path)

    old = time.time() - cons.DUE_SECONDS - 60
    os.utime(tmp_path / cons.ATTEMPT_MARKER, (old, old))
    assert cons.maybe_start_in_background(tmp_path)
    ((args, kwargs),) = started
    assert args[0][1:] == ["-m", "qxub.queue.consolidate"]
    assert kwargs["cwd"] == str(tmp_path), "not the user's cwd (review of #135)"
    assert not cons.maybe_start_in_background(tmp_path), "the clock restarted"


def test_auto_consolidation_can_be_turned_off(tmp_path, monkeypatch):
    monkeypatch.setenv("QXUB_AUTO_CONSOLIDATE", "0")
    (tmp_path / cons.ATTEMPT_MARKER).touch()
    old = time.time() - cons.DUE_SECONDS - 60
    os.utime(tmp_path / cons.ATTEMPT_MARKER, (old, old))
    monkeypatch.setattr(
        cons.subprocess, "Popen", lambda *a, **k: pytest.fail("started")
    )

    assert not cons.maybe_start_in_background(tmp_path)


def test_failed_qsub_leaves_no_run_file(sandbox):
    sandbox.fail_qsub(38, "qsub: Job violates queue limits")

    result = sandbox.run("exec", "--terse", "--", "true")

    assert result.returncode == 38
    runs_dir = sandbox.home / ".local" / "state" / "qxub" / "status" / "runs"
    assert list(runs_dir.glob("*.run")) == []


# -- qxub db -----------------------------------------------------------------------------


def test_qxub_db_consolidate_and_status(sandbox, tmp_path):
    runs_dir = sandbox.home / ".local" / "state" / "qxub" / "status" / "runs"
    runs_dir.mkdir(parents=True)
    finished(runs_dir, "qx-a", "1.gadi-pbs", ago=timedelta(days=30))
    snapshot = tmp_path / "stats" / "qxub-stats.db"
    env = {"QXUB_STATS_DB": str(snapshot)}

    result = sandbox.run("db", "consolidate", env=env)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Snapshot published" in result.stdout and "1 run(s) added" in result.stdout

    result = sandbox.run("db", "status", env=env)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Jobs:       1 (0 with resource usage)" in result.stdout
    assert "Run files:  1" in result.stdout
