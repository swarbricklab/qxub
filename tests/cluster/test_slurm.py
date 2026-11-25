"""qxub against a real Slurm (the ci/slurm container).

Everything here runs real jobs: submission with sbatch, streaming, the end
of a job as qx sees it (normal exit, qxub cancel, scancel, cancelled while
queued, walltime), and sacct's record of it.
"""

import subprocess

import pytest
from cluster_helpers import cluster, needs, start_qx, wait_for

pytestmark = needs("slurm")
__all__ = ["cluster"]  # the fixture, used by name


def squeue():
    """``{job_id: state}`` of the queue's jobs (``PENDING``, ``RUNNING``)."""
    out = subprocess.run(
        ["squeue", "-h", "-o", "%i %T"], capture_output=True, text=True, check=True
    ).stdout
    return dict(line.split() for line in out.splitlines())


@pytest.fixture(autouse=True)
def empty_queue():
    """Each test starts with no jobs left over from the last."""
    wait_for(lambda: not squeue(), timeout=120, what="the queue to empty")


def only_job(state):
    """The one job in ``state``, once there is one."""

    def find():
        jobs = [job_id for job_id, s in squeue().items() if s == state]
        return jobs[0] if len(jobs) == 1 else None

    return find


def test_qx_runs_a_job(cluster):
    result = cluster.run("exec", "--", 'echo "hello from job $SLURM_JOB_ID"')

    assert result.returncode == 0, result.stdout + result.stderr
    assert "hello from job " in result.stdout
    assert "✅ Job completed successfully" in result.stderr


def test_exit_code_reaches_qx(cluster):
    result = cluster.run("exec", "--", "echo failing; exit 3")

    assert result.returncode == 3, result.stdout + result.stderr
    assert "Job failed with exit code 3" in result.stderr


def test_partition_and_resources_reach_the_job(cluster):
    result = cluster.run(
        "exec",
        "--queue",
        "debug",
        "--cpus",
        "2",
        "--mem",
        "100MB",
        "--",
        "echo $SLURM_JOB_PARTITION $SLURM_CPUS_PER_TASK $SLURM_MEM_PER_NODE",
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "debug 2 100"


def test_qxub_cancel_stops_qx(cluster):
    qx = start_qx(cluster, "sleep 300")
    try:
        job_id = wait_for(only_job("RUNNING"), what="the job to start")
        assert cluster.run("cancel", "id", "-f", job_id).returncode == 0
        _, err = qx.communicate(timeout=60)
    finally:
        qx.kill()

    assert qx.returncode == 143, err
    assert "🛑 Job was cancelled (exit code 143)" in err


def test_scancel_outside_qxub_is_noticed(cluster):
    """No run file mark: qx learns it from sacct (CANCELLED)."""
    qx = start_qx(cluster, "sleep 300")
    try:
        job_id = wait_for(only_job("RUNNING"), what="the job to start")
        subprocess.run(["scancel", job_id], check=True)
        _, err = qx.communicate(timeout=60)
    finally:
        qx.kill()

    assert qx.returncode == 143, err
    assert "🛑 Job was cancelled (exit code 143)" in err


def test_job_cancelled_while_queued_is_noticed(cluster):
    """It never starts, so it writes no status files at all."""
    qx = start_qx(cluster, "true", args=["-l", "begin=now+3600"])
    try:
        job_id = wait_for(only_job("PENDING"), what="the job to be queued")
        subprocess.run(["scancel", job_id], check=True)
        _, err = qx.communicate(timeout=60)
    finally:
        qx.kill()

    assert qx.returncode == 143, err
    assert "🛑 Job was cancelled (exit code 143) (from slurm;" in err


def test_walltime_kill_says_so(cluster):
    """Slurm's shortest limit is a minute; it ends the job soon after."""
    result = cluster.run("exec", "--runtime", "1m", "--", "sleep 600", timeout=300)

    assert result.returncode == 143, result.stdout + result.stderr
    assert "Job was stopped at its walltime limit (exit code 143)" in result.stderr


def test_cancel_all_by_state(cluster):
    queued = cluster.run("exec", "--terse", "-l", "begin=now+3600", "--", "true")
    job_id = queued.stdout.strip()
    assert queued.returncode == 0 and job_id, queued.stderr
    wait_for(only_job("PENDING"), what="the job to be queued")

    assert cluster.run("cancel", "all", "-s", "R", "-f").returncode == 1, "none running"
    assert cluster.run("cancel", "all", "-s", "Q", "-f").returncode == 0

    wait_for(lambda: job_id not in squeue(), what="the job to go")


def test_resources_show_reads_sacct(cluster):
    job_id = cluster.run("exec", "--terse", "--", "sleep 2").stdout.strip()
    assert job_id
    cluster.run("monitor", "--quiet", job_id, timeout=120)

    def finished():
        out = cluster.run("resources", "show", job_id).stdout
        return out if "State: finished" in out else None

    out = wait_for(finished, timeout=30, what="sacct to record the job")
    assert "Queue: cpu" in out and "Exit Status: 0" in out, out
