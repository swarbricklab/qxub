"""``qxub monitor`` config defaults and ``qxub status list --all``."""

import pytest
from click.testing import CliRunner

from qxub import monitor_cli, status_cli


@pytest.fixture
def config(monkeypatch):
    """Set config values that ``monitor`` reads."""
    values = {}
    monkeypatch.setattr(
        monitor_cli.config_manager, "get_config_value", lambda key: values.get(key)
    )
    return values


def test_monitor_uses_configured_defaults(config):
    config.update(
        {"monitor.default_suffix": ".gadi-pbs", "monitor.default_interval": 15}
    )
    assert monitor_cli._monitor_defaults(None, None) == (".gadi-pbs", 15)


def test_command_line_beats_config(config):
    config.update(
        {"monitor.default_suffix": ".gadi-pbs", "monitor.default_interval": 15}
    )
    # -i 5 is the old default value, and was overridden by config before
    assert monitor_cli._monitor_defaults(".x", 5) == (".x", 5)


def test_monitor_without_config_uses_5_seconds(config):
    assert monitor_cli._monitor_defaults(None, None) == (None, 5)


def test_bad_configured_interval_warns_and_uses_5(config, capsys):
    config["monitor.default_interval"] = "soon"
    assert monitor_cli._monitor_defaults(None, None) == (None, 5)
    assert "monitor.default_interval" in capsys.readouterr().err


def test_status_list_all_ignores_the_limit(monkeypatch):
    jobs = [{"job_id": f"{n}.gadi-pbs", "status": "completed"} for n in range(30)]
    asked = []

    def get_jobs_by_status(status=None, limit=None):
        asked.append(limit)
        return jobs[:limit] if limit else jobs

    tracker = status_cli.resource_tracker
    monkeypatch.setattr(tracker, "sync_from_run_files", lambda *a, **k: 0)
    monkeypatch.setattr(tracker, "get_jobs_by_status", get_jobs_by_status)
    runner = CliRunner()

    some = runner.invoke(status_cli.status_cli, ["list", "--id-only"])
    every = runner.invoke(status_cli.status_cli, ["list", "--all", "--id-only"])

    assert asked == [20, None]
    assert len(some.output.split()) == 20
    assert len(every.output.split()) == 30
