"""``qxub cancel`` searches the job DB plus one scheduler listing (#143).

Before, ``cancel filter --job-name`` ran ``qstat -fx`` once per row of the
job DB, and ``cancel status``/``filter`` targeted rows left ``submitted``
long ago. Now candidates are the rows whose job the scheduler still lists,
from one ``qstat -w -u`` call, and job names come from the DB.
"""

import os
import sqlite3
from datetime import datetime, timedelta

from qxub.resources.tracker import create_job_resources_schema
from qxub.status_cli import _format_status


def submit(sandbox, name, *extra):
    result = sandbox.run(
        "exec", "--terse", "-P", "ab12", "--name", name, *extra, "--", "sleep", "1"
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout.strip()


def job_db(sandbox):
    (path,) = sandbox.home.rglob("qxub.db")
    return sqlite3.connect(path)


def status_of(sandbox, job_id):
    with job_db(sandbox) as conn:
        return conn.execute(
            "SELECT status FROM job_resources WHERE job_id=?", (job_id,)
        ).fetchone()[0]


def scheduler_queries(sandbox):
    return [c["argv"] for c in sandbox.pbs_calls("qstat")]


def test_job_name_filter_asks_the_scheduler_once(sandbox):
    # Longer than the 15 characters qstat -w -u shows: matched from the DB.
    wanted = [submit(sandbox, f"merge_captures_{i}") for i in range(3)]
    other = submit(sandbox, "align_sample")
    sandbox.finish_job(submit(sandbox, "merge_captures_old"))

    result = sandbox.run(
        "cancel", "filter", "--job-name", "^merge_captures", "-n", "-q"
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.split() == list(reversed(wanted))
    assert other not in result.stdout
    assert scheduler_queries(sandbox) == [["-w", "-u", sandbox.env["USER"]]]
    assert sandbox.pbs_calls("qdel") == []


def test_tag_filter_needs_every_tag(sandbox):
    merge = submit(sandbox, "a", "--tags", "workflow=w1,rule=merge")
    submit(sandbox, "b", "--tags", "workflow=w1,rule=align")
    submit(sandbox, "c", "--tags", "workflow=w2,rule=merge")

    result = sandbox.run(
        "cancel", "filter", "--tag", "workflow=w1", "--tag", "rule=merge", "-n", "-q"
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.split() == [merge]

    bad = sandbox.run("cancel", "filter", "--tag", "merge", "-n")
    assert bad.returncode == 1
    assert "KEY=VALUE" in bad.stderr


def test_filter_cancels_only_what_matched(sandbox):
    wanted = submit(sandbox, "teardown", "--tags", "workflow=w1")
    keep = submit(sandbox, "keep", "--tags", "workflow=w2")

    result = sandbox.run("cancel", "filter", "--tag", "workflow=w1", "-f", "-q")
    assert result.returncode == 0, result.stdout + result.stderr
    assert [c["argv"] for c in sandbox.pbs_calls("qdel")] == [[wanted]]
    assert set(sandbox.live_jobs()) == {keep}


def backdate(sandbox, *job_ids):
    long_ago = (datetime.now() - timedelta(days=30)).isoformat()
    with job_db(sandbox) as conn:
        for job_id in job_ids:
            conn.execute(
                "UPDATE job_resources SET submitted_at=? WHERE job_id=?",
                (long_ago, job_id),
            )


def test_old_rows_the_scheduler_forgot_are_not_targeted(sandbox):
    live = submit(sandbox, "still_queued")
    old = submit(sandbox, "long_gone")
    young = submit(sandbox, "just_finished")
    sandbox.finish_job(old)
    sandbox.finish_job(young)
    backdate(sandbox, old)

    result = sandbox.run("cancel", "status", "-s", "submitted", "-n", "-q")
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.split() == [live]
    # A dry run changes nothing, but says what it saw.
    assert status_of(sandbox, old) == "submitted"
    loud = sandbox.run("cancel", "status", "-s", "submitted", "-n")
    assert "1 old job(s)" in loud.stdout and "dry run" in loud.stdout
    assert status_of(sandbox, old) == "submitted"

    result = sandbox.run("cancel", "status", "-s", "submitted", "-f")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "0 from their joblog, 1 as unknown" in result.stdout
    assert status_of(sandbox, old) == "unknown"
    # Too recent to judge: qstat may simply not list it yet.
    assert status_of(sandbox, young) == "submitted"
    assert [c["argv"] for c in sandbox.pbs_calls("qdel")] == [[live]]


JOBLOG = """\
output
======================================================================================
                  Resource Usage on 2026-10-01 00:00:00:
   Job Id:             {job_id}
   Project:            ab12
   Exit Status:        {code}
   Service Units:      0.50
   NCPUs Requested:    1                      NCPUs Used: 1
                                           CPU Time Used: 00:00:01
   Memory Requested:   1.0GB                 Memory Used: 10.0MB
   Walltime requested: 00:10:00            Walltime Used: 00:00:05
   JobFS requested:    100.0MB                JobFS used: 0B
======================================================================================
"""


def test_tidying_covers_every_status_and_uses_the_joblog(sandbox):
    finished = submit(sandbox, "ended_long_ago")
    lost = submit(sandbox, "lost_running")
    for job_id in (finished, lost):
        sandbox.finish_job(job_id)
    backdate(sandbox, finished, lost)
    # The run file is long gone, but the job's own joblog says how it ended.
    joblog = sandbox.root / "ended_long_ago.log"
    joblog.write_text(JOBLOG.format(job_id=finished, code=3))
    with job_db(sandbox) as conn:
        conn.execute(
            "UPDATE job_resources SET status='running' WHERE job_id=?", (lost,)
        )
        conn.execute(
            "UPDATE job_resources SET joblog_path=? WHERE job_id=?",
            (str(joblog), finished),
        )

    # Asking about 'submitted' jobs tidies the 'running' row too.
    result = sandbox.run("cancel", "status", "-s", "submitted", "-f")
    assert "1 from their joblog, 1 as unknown" in result.stdout
    with job_db(sandbox) as conn:
        row = conn.execute(
            "SELECT status, exit_code FROM job_resources WHERE job_id=?",
            (finished,),
        ).fetchone()
    assert row == ("failed", 3)
    assert status_of(sandbox, lost) == "unknown"


def test_planted_snapshot_rows_do_not_decide_how_jobs_ended(sandbox):
    """Security: anyone in the group can replace the stats snapshot."""
    gone = submit(sandbox, "killed_without_trap")
    sandbox.finish_job(gone)
    backdate(sandbox, gone)
    snapshot = sandbox.root / "stats" / "qxub-stats.db"
    snapshot.parent.mkdir()
    with sqlite3.connect(snapshot) as conn:
        create_job_resources_schema(conn)
        conn.execute(
            "INSERT INTO job_resources (job_id, timestamp, status, exit_code,"
            " completed_at) VALUES (?, 'x', 'completed', 0, '2026-10-01T00:00:00')",
            (gone,),
        )
    sandbox.env["QXUB_STATS_DB"] = str(snapshot)

    result = sandbox.run("cancel", "filter", "--job-name", "nothing", "-f")
    assert "0 from their joblog, 1 as unknown" in result.stdout
    with job_db(sandbox) as conn:
        row = conn.execute(
            "SELECT status, exit_code FROM job_resources WHERE job_id=?", (gone,)
        ).fetchone()
        queue = conn.execute(
            "SELECT status FROM queue WHERE pbs_job_id=?", (gone,)
        ).fetchall()
    assert row == ("unknown", None)
    assert all(status != "completed" for (status,) in queue)
    check = sandbox.run("status", "check", gone)
    assert "completed" not in check.stdout.lower()


def test_nothing_is_cancelled_if_the_scheduler_cant_be_asked(sandbox):
    submit(sandbox, "x")
    sandbox.env["FAKE_QSTAT_EXIT"] = "1"
    sandbox.env["FAKE_QSTAT_STDERR"] = "qstat: cannot connect to server"

    result = sandbox.run("cancel", "filter", "--job-name", "x", "-f")
    assert result.returncode == 1
    assert "cannot connect to server" in result.stderr
    assert sandbox.pbs_calls("qdel") == []


def test_cancel_name_dry_run(sandbox):
    job = submit(sandbox, "merge_captures_1")

    result = sandbox.run("cancel", "name", "-p", "^merge_", "--dry-run", "-q")
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.split() == [job]
    assert sandbox.pbs_calls("qdel") == []


def test_quoted_status_is_normalised():
    conn = sqlite3.connect(":memory:")
    create_job_resources_schema(conn)
    conn.execute(
        "INSERT INTO job_resources (job_id, timestamp, status)"
        " VALUES ('1.gadi-pbs', 'now', '\"running\"')"
    )
    create_job_resources_schema(conn)
    assert conn.execute("SELECT status FROM job_resources").fetchone() == ("running",)


def test_cancelled_has_a_label():
    assert _format_status("cancelled") == "🚫 Cancelled"


def test_joblog_outcome_reads_only_this_users_regular_file(tmp_path):
    """Security: the name is opened once, without following a link, and
    that open file is what is checked and read."""
    from qxub.resources.tracker import _joblog_outcome

    joblog = tmp_path / "job.log"
    joblog.write_text(JOBLOG.format(job_id="1.gadi-pbs", code=3))
    assert _joblog_outcome("1.gadi-pbs", str(joblog)) == ("failed", 3, None)
    assert _joblog_outcome("2.gadi-pbs", str(joblog)) is None

    link = tmp_path / "link.log"
    link.symlink_to(joblog)
    assert _joblog_outcome("1.gadi-pbs", str(link)) is None

    fifo = tmp_path / "fifo.log"
    os.mkfifo(fifo)
    assert _joblog_outcome("1.gadi-pbs", str(fifo)) is None  # doesn't block


def test_joblog_outcome_uses_the_last_usage_block(tmp_path):
    """The job's own output can imitate the scheduler's footer."""
    from qxub.resources.tracker import _joblog_outcome

    joblog = tmp_path / "job.log"
    joblog.write_text(
        JOBLOG.format(job_id="1.gadi-pbs", code=0)
        + JOBLOG.format(job_id="1.gadi-pbs", code=137)
    )
    assert _joblog_outcome("1.gadi-pbs", str(joblog)) == ("failed", 137, None)


def test_only_the_end_of_a_joblog_is_read(tmp_path, monkeypatch):
    from qxub.resources import parser

    monkeypatch.setattr(parser, "JOBLOG_TAIL_BYTES", 2000)
    joblog = tmp_path / "job.log"
    early = JOBLOG.format(job_id="1.gadi-pbs", code=0)
    joblog.write_text(
        early + "x" * 5000 + "\n" + JOBLOG.format(job_id="1.gadi-pbs", code=2)
    )
    assert parser.parse_joblog_resources(str(joblog))["exit_status"] == 2
    joblog.write_text(early + "x" * 5000 + "\n")
    assert parser.parse_joblog_resources(str(joblog)) is None


def _local_runner_log(tmp_path, output: str, code: int):
    """A job log written as the local backend's runner writes it."""
    import subprocess

    from qxub.backends.local import _RUNNER

    script = tmp_path / "job.sh"
    script.write_text(f"printf '%s' {output!r}\nexit {code}\n")
    joblog = tmp_path / "job.log"
    with open(joblog, "ab") as log:
        subprocess.run(
            ["bash", "-c", _RUNNER, "runner", str(script), str(tmp_path), str(joblog)],
            stdout=log,
            stderr=log,
            check=False,
        )
    return joblog


def test_joblog_outcome_reads_a_local_backend_log(tmp_path):
    """Local job logs have no PBS footer, just the runner's last line."""
    from qxub.resources.tracker import _joblog_outcome

    joblog = _local_runner_log(tmp_path, "hello\n", 2)
    assert joblog.read_text().rstrip().endswith("Exit Status: 2")
    assert _joblog_outcome("local-abc", str(joblog)) == ("failed", 2, None)
    # Without a job ID to check, only a local job's own log is read this way
    assert _joblog_outcome("1.gadi-pbs", str(joblog)) is None

    link = tmp_path / "link.log"
    link.symlink_to(joblog)
    assert _joblog_outcome("local-abc", str(link)) is None


def test_joblog_outcome_local_log_needs_the_runners_final_line(tmp_path):
    """The job's output may print an Exit Status line of its own."""
    from qxub.resources.tracker import _joblog_outcome

    joblog = tmp_path / "job.log"
    joblog.write_text("Exit Status: 0\nstill running...\n")
    assert _joblog_outcome("local-abc", str(joblog)) is None

    joblog = _local_runner_log(tmp_path, "Exit Status: 0\n", 1)
    assert _joblog_outcome("local-abc", str(joblog)) == ("failed", 1, None)


def test_settle_gone_local_job_from_its_joblog(tmp_path):
    from qxub.resources.tracker import ResourceTracker

    joblog = _local_runner_log(tmp_path, "", 0)
    tracker = ResourceTracker(db_path=tmp_path / "qxub.db")
    with sqlite3.connect(tmp_path / "qxub.db") as conn:
        conn.execute(
            "INSERT INTO job_resources (job_id, timestamp, status, joblog_path)"
            " VALUES (?,?,?,?)",
            ("local-abc", datetime.now().isoformat(), "running", str(joblog)),
        )
    assert tracker.settle_gone_jobs(["local-abc"])["from_joblog"] == 1
    with sqlite3.connect(tmp_path / "qxub.db") as conn:
        row = conn.execute(
            "SELECT status, exit_code FROM job_resources WHERE job_id='local-abc'"
        ).fetchone()
    assert row == ("completed", 0)
