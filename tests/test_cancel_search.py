"""``qxub cancel`` searches the job DB plus one scheduler listing (#143).

Before, ``cancel filter --job-name`` ran ``qstat -fx`` once per row of the
job DB, and ``cancel status``/``filter`` targeted rows left ``submitted``
long ago. Now candidates are the rows whose job the scheduler still lists,
from one ``qstat -w -u`` call, and job names come from the DB.
"""

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


def test_old_rows_the_scheduler_forgot_are_not_targeted(sandbox):
    live = submit(sandbox, "still_queued")
    old = submit(sandbox, "long_gone")
    young = submit(sandbox, "just_finished")
    sandbox.finish_job(old)
    sandbox.finish_job(young)
    long_ago = (datetime.now() - timedelta(days=30)).isoformat()
    with job_db(sandbox) as conn:
        conn.execute(
            "UPDATE job_resources SET submitted_at=? WHERE job_id=?", (long_ago, old)
        )

    result = sandbox.run("cancel", "status", "-s", "submitted", "-n", "-q")
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.split() == [live]
    assert status_of(sandbox, old) == "unknown"
    # Too recent to judge: qstat may simply not list it yet.
    assert status_of(sandbox, young) == "submitted"
    assert status_of(sandbox, live) == "submitted"


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
