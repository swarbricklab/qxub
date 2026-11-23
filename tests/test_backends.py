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
from qxub.backends import LiveJob, SchedulerBackend, SchedulerQueryError
from qxub.backends.pbs import PbsBackend, parse_qstat_full_listing
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

    def is_finished(self, details):
        return bool(details) and details.get("done", False)

    def cancel(self, job_id):
        self.calls.append(("cancel", job_id))
        before = len(self.jobs)
        self.jobs = [j for j in self.jobs if j.job_id != job_id]
        return (True, "") if len(self.jobs) < before else (False, "no such job")

    def looks_like_job_id(self, value):
        return value.startswith("fake-")


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
        backend = backends.backend_for_type("slurm")
    assert isinstance(backend, PbsBackend)
    assert "No scheduler backend for platform type 'slurm'" in caplog.text


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

    found = PbsBackend().job_details(ids)

    assert sorted(found) == sorted(ids[:-1])
    assert [len(c["argv"]) - 1 for c in fake_pbs.pbs_calls("qstat")] == [100, 50]
    assert PbsBackend().is_finished(found["0.gadi-pbs"])


def test_job_resource_data_matches_a_short_job_id(fake_pbs):
    fake_pbs.set_qstat("123", FX.format(id="123.gadi-pbs"))
    assert scheduler.get_job_resource_data("123")["exit_status"] == 0


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

    assert cons.query_qstat(["fake-1", "fake-2"]) == fake_backend.details
    assert cons._qstat_finished(fake_backend.details["fake-1"])
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
