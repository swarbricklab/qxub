"""qxub against a real OpenPBS (the ci/pbs container).

Unlike Gadi's PBS, nothing here appends an ``Exit Status:`` line to the
job log, as at other PBS sites (UNSW Katana): these show how qx learns how
a job ended without it.
"""

import subprocess

import pytest
from cluster_helpers import cluster, needs, start_qx, wait_for

pytestmark = needs("pbs")
__all__ = ["cluster"]  # the fixture, used by name


def live_jobs():
    """IDs of the jobs PBS hasn't finished (queued, running, exiting)."""
    out = subprocess.run(["qselect"], capture_output=True, text=True, check=False)
    return out.stdout.split()


@pytest.fixture(autouse=True)
def empty_queue():
    """Each test starts with no jobs left over: PBS lists a job for a few
    seconds after qx has seen it end, and a cancelled one for its kill
    delay, so a test could otherwise pick up the last test's job."""
    wait_for(lambda: not live_jobs(), timeout=120, what="the queue to empty")


def running_job():
    """The one running job's ID, once there is one."""
    out = subprocess.run(
        ["qstat", "-r", "-w"], capture_output=True, text=True, check=False
    ).stdout
    jobs = [line.split()[0] for line in out.splitlines() if line[:1].isdigit()]
    return jobs[0] if len(jobs) == 1 else None


def test_qx_runs_a_job(cluster):
    result = cluster.run("exec", "--", 'echo "hello from job $PBS_JOBID"')

    assert result.returncode == 0, result.stdout + result.stderr
    assert "hello from job " in result.stdout
    assert "✅ Job completed successfully" in result.stderr


def test_exit_code_reaches_qx(cluster):
    result = cluster.run("exec", "--", "echo failing; exit 3")

    assert result.returncode == 3, result.stdout + result.stderr
    assert "Job failed with exit code 3" in result.stderr


def test_queue_and_resources_reach_the_job(cluster):
    result = cluster.run(
        "exec",
        "--queue",
        "express",
        "--cpus",
        "2",
        "--mem",
        "100MB",
        "--",
        "echo $PBS_QUEUE $NCPUS",
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "express 2"


def test_qxub_cancel_stops_qx(cluster):
    qx = start_qx(cluster, "sleep 300")
    try:
        job_id = wait_for(running_job, what="the job to start")
        cancel = cluster.run("cancel", "id", "-f", job_id)
        assert cancel.returncode == 0, cancel.stdout + cancel.stderr
        _, err = qx.communicate(timeout=60)
    finally:
        qx.kill()

    assert qx.returncode == 143, err
    assert "🛑 Job was cancelled (exit code 143)" in err


def test_qdel_outside_qxub_is_noticed(cluster):
    """No Exit Status line in the job log: qx asks PBS (Exit_status 271)."""
    qx = start_qx(cluster, "sleep 300")
    try:
        job_id = wait_for(running_job, what="the job to start")
        subprocess.run(["qdel", job_id], check=True)
        _, err = qx.communicate(timeout=90)
    finally:
        qx.kill()

    assert "🛑 Job was cancelled (exit code 271)" in err, err


def test_walltime_kill_says_so(cluster):
    result = cluster.run("exec", "--runtime", "1m", "--", "sleep 600", timeout=300)

    assert result.returncode != 0, result.stdout + result.stderr
    assert "Job was stopped at its walltime limit (exit code -29)" in result.stderr


def test_cancel_all_by_state(cluster):
    qx = start_qx(cluster, "sleep 300")
    try:
        job_id = wait_for(running_job, what="the job to start")
        assert (
            cluster.run("cancel", "all", "-s", "Q", "-f").returncode == 1
        ), "none queued"
        assert cluster.run("cancel", "all", "-s", "R", "-f").returncode == 0
        qx.communicate(timeout=90)
    finally:
        qx.kill()

    assert qx.returncode == 143
    wait_for(lambda: job_id not in live_jobs(), what="PBS to end the job")


def test_resources_show_reads_qstat(cluster):
    job_id = cluster.run("exec", "--terse", "--", "sleep 2").stdout.strip()
    assert job_id
    cluster.run("monitor", "--quiet", job_id, timeout=120)

    def finished():
        out = cluster.run("resources", "show", job_id).stdout
        return out if "State: finished" in out else None

    out = wait_for(finished, timeout=30, what="PBS to record the job")
    assert "Queue: normal" in out and "Exit Status: 0" in out, out
