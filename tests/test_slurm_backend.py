"""The Slurm backend: sbatch, squeue, sacct, scancel and srun (#65, step 8).

The parsers run on captured output; the rest runs ``qxub`` against the fake
Slurm commands in ``tests/fakes/slurm``.
"""

from datetime import datetime

import pytest

from qxub import backends
from qxub.backends import (
    ENDED_BY_CANCEL,
    ENDED_BY_MEMORY,
    ENDED_BY_WALLTIME,
    JobEnd,
    JobSpec,
    JobState,
    SlurmBackend,
    SubmitError,
)
from qxub.backends.slurm import (
    job_info_from_sacct,
    mail_types,
    parse_sacct_output,
    parse_sbatch_output,
    parse_slurm_duration,
    parse_squeue_output,
)
from qxub.config import site
from qxub.core import scheduler
from qxub.queue import runs
from qxub.resources.spec import ResourceSpec

GB = 1024**3


def _sacct(*lines):
    return "\n".join(lines) + "\n"


# JobID|JobName|Partition|State|ExitCode|Submit|Start|End|NodeList|ReqTRES|
# AllocTRES|TimelimitRaw|ElapsedRaw|TotalCPU|MaxRSS
COMPLETED = _sacct(
    "5000|align|cpu|COMPLETED|0:0|2026-11-24T09:00:00|2026-11-24T09:10:00|"
    "2026-11-24T10:00:00|node01|billing=4,cpu=4,gres/gpu=1,mem=16G,node=1|"
    "billing=4,cpu=4,mem=16G,node=1|60|3000|03:10:00|",
    "5000.batch|batch||COMPLETED|0:0|2026-11-24T09:10:00|2026-11-24T09:10:00|"
    "2026-11-24T10:00:00|node01||cpu=4,mem=16G,node=1||3000|03:09:59|6291456K",
    "5000.extern|extern||COMPLETED|0:0|2026-11-24T09:10:00|2026-11-24T09:10:00|"
    "2026-11-24T10:00:00|node01||cpu=4,mem=16G,node=1||3000|00:00:01|1024K",
)


def _job(state, code, job_id="5001"):
    return _sacct(
        f"{job_id}|j|cpu|{state}|{code}|2026-11-24T09:00:00|Unknown|"
        "2026-11-24T09:05:00|None assigned|cpu=1,mem=1G,node=1||10|0|00:00:00|"
    )


def _info(text, job_id="5001"):
    return job_info_from_sacct(job_id, parse_sacct_output(text)[job_id])


# -- commands -----------------------------------------------------------------


def test_sbatch_command_from_a_job_spec():
    spec = JobSpec(
        script="/q/job.sh",
        name="align",
        queue="gpu",
        account="abc",
        resources=ResourceSpec(
            cpus=4,
            memory=16 * GB,
            walltime=2 * 3600,
            gpus=1,
            disk=100 * 1024**2,
            native={"constraint": "skylake"},
        ),
        job_log="/q/align.log",
        mail_to="me@example.org",
        mail_events="ae",
        env_file="/s/qx-1.env",
        native_args=["--qos=high"],
    )

    assert SlurmBackend().submit_command(spec) == [
        "env",
        "QXUB_JOB_FILE=/s/qx-1.env",
        "sbatch",
        "--parsable",
        "--export=ALL",
        "--job-name=align",
        "--partition=gpu",
        "--account=abc",
        "--cpus-per-task=4",
        "--mem=16384M",
        "--time=02:00:00",
        "--gres=gpu:1",
        "--tmp=100M",
        "--constraint=skylake",
        "--output=/q/align.log",
        "--mail-user=me@example.org",
        "--mail-type=FAIL,END",
        "--qos=high",
        "/q/job.sh",
    ]


def test_interactive_session_is_srun_pty():
    spec = JobSpec(
        script="/tmp/session.sh",
        name="qxub-interactive",
        queue="cpu",
        resources=ResourceSpec(cpus=2),
        job_log="/q/x.log",
        native_args=["--qos=high"],
    )
    assert SlurmBackend().interactive_command(spec) == [
        "srun",
        "--pty",
        "--job-name=qxub-interactive",
        "--partition=cpu",
        "--cpus-per-task=2",
        "bash",
        "/tmp/session.sh",
    ]


def test_mail_letters_become_mail_types():
    assert mail_types("abe") == "FAIL,BEGIN,END"
    assert mail_types("n") == "NONE"
    assert mail_types(None) == "FAIL,END"


def test_sbatch_output():
    assert parse_sbatch_output("12345\n") == "12345"
    assert parse_sbatch_output("12345;cluster2\n") == "12345"
    with pytest.raises(SubmitError):
        parse_sbatch_output("sbatch: error: something\n")


def test_job_ids():
    backend = SlurmBackend()
    for job_id in ("12345", "12345_7", "12345_[1-10]", "12345+1"):
        assert backend.looks_like_job_id(job_id), job_id
    for value in ("qx-abc", "local-abc", "12345.gadi-pbs", ""):
        assert not backend.looks_like_job_id(value), value


# -- listing ------------------------------------------------------------------


def test_squeue_listing_and_states():
    jobs = parse_squeue_output(
        "5000|cpu|R|align with spaces|and a bar\n5001|gpu|PD|train\nnot a job\n"
    )
    assert [(j.job_id, j.queue, j.state, j.name) for j in jobs] == [
        ("5000", "cpu", "R", "align with spaces|and a bar"),
        ("5001", "gpu", "PD", "train"),
    ]
    backend = SlurmBackend()
    assert backend.live_state("PD") is JobState.QUEUED
    assert backend.live_state("CG") is JobState.RUNNING
    assert backend.live_state("S") is JobState.HELD
    assert backend.live_state("??") is JobState.UNKNOWN


# -- sacct --------------------------------------------------------------------


def test_finished_job_in_neutral_terms():
    info = _info(COMPLETED, "5000")

    assert info.state is JobState.FINISHED
    assert info.exit_status == 0 and info.ended_by is None
    assert (info.name, info.queue, info.host) == ("align", "cpu", "node01")
    assert info.submitted == datetime(2026, 11, 24, 9, 0)
    assert info.queue_wait == 600 and info.run_time == 50 * 60
    usage = info.usage
    assert (usage.cpus_requested, usage.cpus_used) == (4, 4)
    assert usage.cpu_time == 3 * 3600 + 600
    assert usage.memory_requested == 16 * GB
    assert usage.memory_used == 6 * GB, "the largest step's MaxRSS"
    assert (usage.walltime_requested, usage.walltime_used) == (3600, 3000)
    assert usage.gpus_requested == 1
    assert info.native["steps"]["batch"]["MaxRSS"] == "6291456K"


@pytest.mark.parametrize(
    "state, code, exit_status, ended_by",
    [
        ("TIMEOUT", "0:15", 143, ENDED_BY_WALLTIME),
        ("TIMEOUT", "0:0", 143, ENDED_BY_WALLTIME),
        ("CANCELLED by 1234", "0:0", 143, ENDED_BY_CANCEL),  # while queued
        ("CANCELLED by 1234", "0:15", 143, ENDED_BY_CANCEL),
        ("OUT_OF_MEMORY", "0:125", 137, ENDED_BY_MEMORY),
        ("FAILED", "3:0", 3, None),
        ("NODE_FAIL", "0:0", 1, None),
        ("COMPLETED", "0:0", 0, None),
    ],
)
def test_why_a_job_ended_comes_from_its_state(state, code, exit_status, ended_by):
    info = _info(_job(state, code))
    assert info.is_finished
    assert (info.exit_status, info.ended_by) == (exit_status, ended_by)


def test_unfinished_jobs():
    for state in ("PENDING", "RUNNING", "COMPLETING"):
        info = _info(_job(state, "0:0"))
        assert not info.is_finished and info.exit_status is None, state
    assert _info(_job("SUSPENDED", "0:0")).state is JobState.HELD


def test_slurm_durations():
    assert parse_slurm_duration("01:02:03") == 3723
    assert parse_slurm_duration("2-00:00:01") == 2 * 86400 + 1
    assert parse_slurm_duration("05:30.123") == 330
    assert parse_slurm_duration("") is None
    assert parse_slurm_duration("Unknown") is None


# -- qxub against the fake Slurm ----------------------------------------------


def test_exec_submits_with_sbatch(sandbox):
    sandbox.use_slurm_platform()

    result = sandbox.run("exec", "--terse", "--cpus", "2", "--", "echo hi")

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "5000"
    assert sandbox.pbs_calls() == []
    (call,) = sandbox.slurm_calls("sbatch")
    assert "--cpus-per-task=2" in call["argv"]
    assert not any(
        a.startswith("--partition") for a in call["argv"]
    ), "no queue asked for and none named normal: Slurm's default partition"
    settings = open(call["QXUB_JOB_FILE"], encoding="utf-8").read()
    assert "QXUB_JOB_ID_VAR=SLURM_JOB_ID" in settings


def test_exec_names_a_partition_and_dry_runs(sandbox):
    sandbox.use_slurm_platform()

    dry = sandbox.run("exec", "--dry", "--queue", "gpu", "--", "true")

    assert dry.returncode == 0, dry.stderr
    assert "env QXUB_JOB_FILE=" in dry.stdout and "sbatch --parsable" in dry.stdout
    assert "--partition=gpu" in dry.stdout
    assert sandbox.slurm_calls("sbatch") == []


def test_cancel_and_cancel_all(sandbox):
    sandbox.use_slurm_platform()
    first = sandbox.run("exec", "--terse", "--", "true").stdout.strip()
    second = sandbox.run("exec", "--terse", "--", "true").stdout.strip()

    assert sandbox.run("cancel", "id", "-f", first).returncode == 0
    assert [c["argv"] for c in sandbox.slurm_calls("scancel")] == [[first]]

    queued = sandbox.run("cancel", "all", "-s", "Q", input="n\n")
    assert second in queued.stdout, queued.stdout + queued.stderr
    assert sandbox.run("cancel", "all", "-s", "R", "-f").returncode == 1, "none running"

    assert sandbox.run("cancel", "all", "-f").returncode == 0
    assert sandbox.slurm_live_jobs() == {}


def test_qx_notices_a_job_cancelled_while_queued(sandbox):
    """It never runs, so it writes no status files; sacct says it ended."""
    sandbox.use_slurm_platform()
    sandbox.write_config("monitoring:\n  joblog_check_interval_sec: 1\n")
    sandbox.set_sacct("5000", _job("CANCELLED by 1234", "0:0", job_id="5000"))

    result = sandbox.run("exec", "--", "echo hi", timeout=30)

    assert result.returncode == 143, result.stdout + result.stderr
    assert "🛑 Job was cancelled (exit code 143) (from slurm;" in result.stderr


def test_resources_show_a_slurm_job(sandbox):
    sandbox.use_slurm_platform()
    sandbox.set_sacct("5000", COMPLETED)

    result = sandbox.run("resources", "show", "5000")

    assert result.returncode == 0, result.stderr
    assert "State: finished" in result.stdout and "Host: node01" in result.stdout
    assert "16.0GB" in result.stdout and "6.0GB" in result.stdout


# -- pieces qx and status use -------------------------------------------------


def test_job_end_from_sacct(monkeypatch):
    backend = SlurmBackend()
    monkeypatch.setattr(
        backend, "job_details", lambda ids: parse_sacct_output(_job("TIMEOUT", "0:15"))
    )
    assert backend.job_end("5001") == JobEnd(143, ENDED_BY_WALLTIME)

    monkeypatch.setattr(
        backend, "job_details", lambda ids: parse_sacct_output(_job("RUNNING", "0:0"))
    )
    assert backend.job_end("5001") is None


def test_reason_the_scheduler_gave_is_kept(tmp_path):
    path = tmp_path / "qx-s.run"
    path.write_text(runs.format_fields({"state": "finished", "exit_code": 137}))

    reason = scheduler.job_end_reason("5001", 137, path, ENDED_BY_MEMORY)

    assert reason == ENDED_BY_MEMORY
    assert runs.read_run(path)["ended_by"] == ENDED_BY_MEMORY
    assert (
        backends.describe_ended_by(reason) == "stopped for exceeding its memory limit"
    )


class _Platform:
    def __init__(self, queues, default_queue=None):
        self.queues = dict.fromkeys(queues)
        self.default_queue = default_queue


@pytest.mark.parametrize(
    "platform, queue",
    [
        (None, "normal"),
        (_Platform(["normal", "express"]), "normal"),
        (_Platform(["cpu", "gpu"]), None),
        (_Platform(["cpu", "gpu"], default_queue="cpu"), "cpu"),
    ],
)
def test_default_queue_comes_from_the_platform(monkeypatch, platform, queue):
    monkeypatch.setattr(site, "current_platform", lambda: platform)
    assert site.default_queue() == queue


def test_status_adds_gadis_suffix_only_under_pbs(monkeypatch):
    from qxub import status_cli

    monkeypatch.setattr(backends, "get_backend", lambda platform=None: SlurmBackend())
    assert status_cli._full_job_id("5000") == "5000"
    monkeypatch.setattr(
        backends, "get_backend", lambda platform=None: backends.PbsBackend()
    )
    assert status_cli._full_job_id("5000") == "5000.gadi-pbs"
