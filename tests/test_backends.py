"""Scheduler backends (#65): the one place qxub asks a scheduler about jobs.

The platform definition's ``type`` picks the backend. ``cancel``, ``status``
lookups, ``resources show``, interactive completion and consolidation all go
through it, so a second scheduler only needs a new backend class. The tests
that swap in :class:`FakeBackend` show the CLI no longer runs PBS commands
itself.
"""

import logging
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from qxub import backends, cancel_cli
from qxub.backends import (
    JobInfo,
    JobSpec,
    JobState,
    LiveJob,
    ResourceSpec,
    SchedulerBackend,
    SchedulerQueryError,
    SubmitError,
)
from qxub.backends.pbs import PbsBackend, parse_qstat_full_listing, qsub_argv
from qxub.core import scheduler
from qxub.queue import consolidate as cons


class FakeBackend(SchedulerBackend):
    """A scheduler in memory: no commands are run."""

    name = "fake"

    def __init__(self, jobs=(), details=None, fail=None):
        self.jobs = list(jobs)
        self.details = details or {}
        self.fail = fail
        self.calls = []

    def list_jobs(self, user=None, full_names=False):
        self.calls.append(("list_jobs", user, full_names))
        if self.fail:
            raise SchedulerQueryError(self.fail)
        return list(self.jobs)

    def job_details(self, job_ids):
        job_ids = list(job_ids)
        self.calls.append(("job_details", job_ids))
        return {j: self.details[j] for j in job_ids if j in self.details}

    def job_info(self, job_ids):
        return {
            job_id: JobInfo(
                job_id,
                JobState.FINISHED if details.get("done") else JobState.RUNNING,
                exit_status=details.get("exit_status"),
                native=details,
            )
            for job_id, details in self.job_details(job_ids).items()
        }

    def cancel(self, job_id):
        self.calls.append(("cancel", job_id))
        before = len(self.jobs)
        self.jobs = [j for j in self.jobs if j.job_id != job_id]
        return (True, "") if len(self.jobs) < before else (False, "no such job")

    def looks_like_job_id(self, value):
        return value.startswith("fake-")

    def submit_command(self, spec):
        return ["fake-submit", spec.name, spec.script]

    def interactive_command(self, spec):
        return ["fake-interactive", spec.name, spec.script]

    def submit(self, spec):
        self.calls.append(("submit", spec))
        if self.fail:
            raise SubmitError(self.fail, exit_code=3, stderr=self.fail)
        job_id = f"fake-{len(self.jobs) + 1}"
        self.jobs.append(LiveJob(job_id, spec.queue or "", "Q", spec.name))
        return job_id


@pytest.fixture
def fake_backend(monkeypatch, tmp_path):
    """Use a :class:`FakeBackend`, and a throwaway job DB, in this process."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("QXUB_DB_PATH", str(tmp_path / "qxub.db"))
    backend = FakeBackend()
    monkeypatch.setattr(backends, "get_backend", lambda platform=None: backend)
    return backend


# -- choosing a backend --------------------------------------------------------


@pytest.mark.parametrize("kind", ["pbs_pro", "pbs", "PBSPro", "torque", None, ""])
def test_pbs_platform_types_select_the_pbs_backend(kind):
    assert isinstance(backends.backend_for_type(kind), PbsBackend)


def test_unknown_platform_type_warns_and_uses_pbs(caplog):
    with caplog.at_level(logging.WARNING, logger="qxub.backends"):
        backend = backends.backend_for_type("lsf")
    assert isinstance(backend, PbsBackend)
    assert "No scheduler backend for platform type 'lsf'" in caplog.text


def test_platform_type_picks_the_backend(monkeypatch):
    monkeypatch.setitem(backends.BACKENDS, "fake", FakeBackend)
    platform = SimpleNamespace(type="fake")
    assert isinstance(backends.get_backend(platform), FakeBackend)


def test_without_a_platform_qxub_assumes_pbs(monkeypatch):
    import qxub.platforms

    def no_platform():
        raise RuntimeError("no platform files")

    monkeypatch.setattr(qxub.platforms, "get_current_platform", no_platform)
    assert isinstance(backends.get_backend(), PbsBackend)


def test_gadi_platform_definition_selects_pbs():
    from conftest import GADI_PLATFORM_DIR

    from qxub.platforms.core import PlatformLoader

    gadi = PlatformLoader([GADI_PLATFORM_DIR]).get_platform("nci_gadi")
    assert gadi.type == "pbs_pro"
    assert isinstance(backends.get_backend(gadi), PbsBackend)


# -- the PBS backend -----------------------------------------------------------

QSTAT_F = """\
Job Id: 101.gadi-pbs
    Job_Name = merge_captures_long_name
    Job_Owner = alice@gadi-login-01.gadi.nci.org.au
    job_state = R
    queue = normal-exec

Job Id: 102.gadi-pbs
    Job_Name = someone_elses
    Job_Owner = bob@gadi-login-02.gadi.nci.org.au
    job_state = Q
    queue = express-exec
"""


def test_full_listing_keeps_whole_names_and_only_your_jobs():
    assert parse_qstat_full_listing(QSTAT_F, "alice") == [
        LiveJob("101.gadi-pbs", "normal-exec", "R", "merge_captures_long_name")
    ]
    assert len(parse_qstat_full_listing(QSTAT_F)) == 2


def test_pbs_job_ids():
    pbs = PbsBackend()
    for good in ("123", "123.gadi-pbs", "123[].gadi-pbs", "123[4].gadi-pbs"):
        assert pbs.looks_like_job_id(good)
    for bad in ("qxi_abc", "qx-1234", "", "gadi-pbs:"):
        assert not pbs.looks_like_job_id(bad)
    assert pbs.same_job("123.gadi-pbs") == pbs.same_job("123") == "123"


@pytest.fixture
def fake_pbs(sandbox, monkeypatch):
    """The fake PBS commands on PATH for code called from the test."""
    monkeypatch.setenv("PATH", sandbox.env["PATH"])
    monkeypatch.setenv("FAKE_PBS_STATE", str(sandbox.pbs_state))
    return sandbox


FX = "Job Id: {id}\n    job_state = F\n    Exit_status = 0\n"


def test_job_details_asks_in_batches_and_skips_unknown_jobs(fake_pbs):
    ids = [f"{n}.gadi-pbs" for n in range(150)]
    for job_id in ids[:-1]:
        fake_pbs.set_qstat(job_id, FX.format(id=job_id))

    found = PbsBackend().job_info(ids)

    assert sorted(found) == sorted(ids[:-1])
    assert [len(c["argv"]) - 1 for c in fake_pbs.pbs_calls("qstat")] == [100, 50]
    assert found["0.gadi-pbs"].is_finished and found["0.gadi-pbs"].exit_status == 0


def test_job_resource_data_matches_a_short_job_id(fake_pbs):
    fake_pbs.set_qstat("123", FX.format(id="123.gadi-pbs"))
    assert scheduler.get_job_resource_data("123")["exit_status"] == 0
    assert scheduler.get_job_info("123").job_id == "123.gadi-pbs"


def test_cancel_reports_why_it_failed(fake_pbs, monkeypatch):
    monkeypatch.setenv("FAKE_QDEL_EXIT", "35")
    monkeypatch.setenv("FAKE_QDEL_STDERR", "qdel: Unknown Job Id 9.gadi-pbs")
    assert PbsBackend().cancel("9.gadi-pbs") == (
        False,
        "qdel: Unknown Job Id 9.gadi-pbs",
    )


def test_listing_failure_is_a_query_error(fake_pbs, monkeypatch):
    monkeypatch.setenv("FAKE_QSTAT_EXIT", "2")
    monkeypatch.setenv("FAKE_QSTAT_STDERR", "qstat: cannot connect to server")
    with pytest.raises(SchedulerQueryError, match="cannot connect to server"):
        PbsBackend().list_jobs("alice")


# -- callers go through the backend --------------------------------------------


def test_cancel_all_lists_and_cancels_through_the_backend(fake_backend):
    fake_backend.jobs = [
        LiveJob("fake-1", "normal", "R", "a"),
        LiveJob("fake-2", "normal", "Q", "b"),
        LiveJob("fake-3", "express", "R", "c"),
    ]

    result = CliRunner().invoke(
        cancel_cli.cancel_cli, ["all", "-u", "alice", "-s", "R", "-q", "normal", "-f"]
    )

    assert result.exit_code == 0, result.output
    assert ("list_jobs", "alice", False) in fake_backend.calls
    assert [c for c in fake_backend.calls if c[0] == "cancel"] == [("cancel", "fake-1")]


def test_cancel_name_asks_for_whole_names(fake_backend):
    fake_backend.jobs = [
        LiveJob("fake-1", "normal", "R", "merge_captures_long_name"),
        LiveJob("fake-2", "normal", "R", "align"),
    ]

    result = CliRunner().invoke(
        cancel_cli.cancel_cli, ["name", "^merge", "--pattern", "-n", "-q"]
    )

    assert result.exit_code == 0, result.output
    assert result.output.split() == ["fake-1"]
    assert fake_backend.calls == [("list_jobs", None, True)]


def test_cancel_wait_polls_the_backend_listing(fake_backend):
    fake_backend.jobs = [LiveJob("fake-1", "normal", "R", "a")]
    assert not cancel_cli._wait_for_deletion("fake-1", timeout=0)
    fake_backend.jobs = []
    assert cancel_cli._wait_for_deletion("fake-1", timeout=5)


def test_cancel_wait_gives_up_if_the_scheduler_cannot_be_asked(fake_backend):
    fake_backend.fail = "down"
    assert not cancel_cli._wait_for_deletion("fake-1", timeout=5)


def test_consolidation_and_resource_lookups_use_the_backend(fake_backend):
    fake_backend.details = {"fake-1": {"exit_status": 3, "done": True}}

    (info,) = cons.query_jobs(["fake-1", "fake-2"]).values()
    assert info.is_finished and info.exit_status == 3
    assert scheduler.get_job_info("fake-1") == info
    assert scheduler.get_job_info("fake-2") is None
    assert scheduler.get_job_resource_data("fake-1") == {"exit_status": 3, "done": True}
    assert scheduler.get_job_resource_data("fake-2") is None
    assert scheduler.qdel("fake-9", quiet=True) is False


class ShortIdBackend(FakeBackend):
    """Like PBS: ``fake-1`` and ``fake-1.server`` are the same job."""

    def same_job(self, job_id):
        return job_id.split(".", 1)[0]


def test_cancel_id_marks_the_db_row_for_a_short_id(fake_backend, monkeypatch):
    from qxub.resources import ResourceTracker

    backend = ShortIdBackend(jobs=[LiveJob("fake-1", "normal", "R", "a")])
    monkeypatch.setattr(backends, "get_backend", lambda platform=None: backend)
    tracker = ResourceTracker()
    for job_id in ("fake-1.server", "fake-10.server"):
        tracker.log_job_submitted(job_id, command="sleep 60")

    result = CliRunner().invoke(cancel_cli.cancel_cli, ["id", "fake-1", "-f"])

    assert result.exit_code == 0, result.output
    assert tracker.get_job_status("fake-1.server")["status"] == "cancelled"
    assert tracker.get_job_status("fake-10.server")["status"] == "submitted"


def test_cancel_id_shows_why_it_failed(fake_backend):
    result = CliRunner().invoke(cancel_cli.cancel_cli, ["id", "fake-9", "-f"])

    assert result.exit_code == 1
    assert "Could not cancel fake-9: no such job" in result.output


# -- submission ----------------------------------------------------------------


def test_qsub_argv_renders_a_job_spec():
    spec = JobSpec(
        script="/t/job.pbs",
        name="j",
        queue="normal",
        account="ab12",
        resources=ResourceSpec(cpus=2, memory=4 * 1024**3),
        job_log="/logs/my job.log",
        mail_to="me@example.org",
        env={"cmd_b64": "ZWNobw==", "bind": "/a,/b", "odd": "it's,here", "p": "/x y"},
        native_args=["-W umask=0027"],
    )

    assert qsub_argv(spec) == [
        "qsub",
        "-v",
        "cmd_b64=ZWNobw==,bind='/a,/b',odd=\"it's,here\",p=/x y",
        "-N", "j", "-q", "normal", "-P", "ab12",
        "-l", "ncpus=2", "-l", "mem=4GB",
        "-o", "/logs/my job.log",
        "-M", "me@example.org", "-m", "ae",
        "-W", "umask=0027",
        "/t/job.pbs",
    ]  # fmt: skip
    assert qsub_argv(JobSpec(script="s", name="n")) == ["qsub", "-N", "n", "s"]


def test_rejected_job_exits_2(fake_pbs, monkeypatch):
    monkeypatch.setenv("FAKE_QSUB_EXIT", "166")
    monkeypatch.setenv("FAKE_QSUB_STDERR", "qsub: Illegal attribute or resource value")

    with pytest.raises(SubmitError) as excinfo:
        PbsBackend().submit(JobSpec(script="s", name="n"))

    assert excinfo.value.exit_code == 2
    assert "Illegal attribute" in excinfo.value.stderr
    assert issubclass(scheduler.QsubError, SubmitError)


@pytest.fixture
def submit_env(fake_backend, monkeypatch, tmp_path):
    """Run submissions from a scratch directory with no platform."""
    from qxub.config import site

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(site, "current_platform", lambda: None)
    return fake_backend


def test_submit_job_goes_through_the_backend(submit_env):
    from qxub.execution import submit_job
    from qxub.resources import ResourceTracker

    result = submit_job(
        ["echo", "hi"],
        env="myenv",
        queue="normal",
        project="ab12",
        name="greet",
        resources={"ncpus": 2, "mem": "4GB"},
    )

    (spec,) = [c[1] for c in submit_env.calls if c[0] == "submit"]
    assert result.job_id == "fake-1"
    assert (spec.name, spec.queue, spec.account) == ("greet", "normal", "ab12")
    assert spec.resources == ResourceSpec(cpus=2, memory=4 * 1024**3)
    assert spec.env["env"] == "myenv"
    assert spec.env["QXUB_RUN_ID"] == result.virtual_id
    assert result.submission_command == "fake-submit greet " + spec.script
    assert ResourceTracker().get_job_status("fake-1")["status"] == "submitted"


def test_refused_submit_job_raises_and_removes_the_run_file(submit_env):
    from qxub.execution import QsubError, submit_job

    submit_env.fail = "queue is closed"
    with pytest.raises(QsubError, match="queue is closed"):
        submit_job(["echo", "hi"], queue="normal", project="ab12")

    (spec,) = [c[1] for c in submit_env.calls if c[0] == "submit"]
    from pathlib import Path

    assert not Path(spec.env["QXUB_RUN_FILE"]).exists()
    assert spec.env_file and not Path(spec.env_file).exists()


def test_exec_submits_through_the_backend(submit_env):
    from qxub.exec_cli import exec_cli

    result = CliRunner().invoke(
        exec_cli, ["--terse", "-P", "ab12", "-q", "normal", "--", "echo", "hi"]
    )

    assert result.exit_code == 0, result.output
    assert result.output.strip().splitlines()[-1] == "fake-1"
    (spec,) = [c[1] for c in submit_env.calls if c[0] == "submit"]
    assert spec.account == "ab12"


def test_exec_shows_why_the_scheduler_refused(submit_env):
    from qxub.exec_cli import exec_cli

    submit_env.fail = "queue is closed"
    result = CliRunner().invoke(
        exec_cli, ["--terse", "-P", "ab12", "-q", "normal", "--", "echo", "hi"]
    )

    assert result.exit_code == 3
    assert "queue is closed" in result.output


# -- interactive sessions ----------------------------------------------------


def test_qsub_interactive_argv():
    from qxub.backends.pbs import qsub_interactive_argv

    spec = JobSpec(
        script="/t/session.sh",
        name="qxub-interactive",
        queue="normal",
        account="ab12",
        resources=ResourceSpec(cpus=2, native={"storage": "gdata/ab12"}),
        # Batch-only settings a session doesn't get:
        job_log="/l/x.log",
        mail_to="me@example.org",
        env={"cmd_b64": "x"},
        env_file="/r/x.env",
        native_args=["-W umask=0027"],
    )
    assert qsub_interactive_argv(spec) == [
        "qsub", "-I", "-N", "qxub-interactive", "-q", "normal", "-P", "ab12",
        "-l", "ncpus=2", "-l", "storage=gdata/ab12", "--", "/t/session.sh",
    ]  # fmt: skip
    assert PbsBackend().interactive_command(spec) == qsub_interactive_argv(spec)
    assert PbsBackend.job_id_variable == "PBS_JOBID"
