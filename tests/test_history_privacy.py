"""History and the job database are private to their user.

Records hold full command lines and the values of ``--var``/``--vars``, which
often carry credentials. Files must be owner-only whatever the umask, and
variable values must not be stored at all.
"""

import os
import stat
import sys
from types import SimpleNamespace

import pytest

from qxub.history.manager import HistoryManager, redact_argv, redact_params

SECRET = "s3cr3t-T0ken"
ARGV = [
    "qxub",
    "exec",
    "--env",
    "dt3",
    "--var",
    f"API_TOKEN={SECRET}",
    f"--var=OTHER={SECRET}",
    "--vars",
    f"A={SECRET},B=x,{SECRET}",
    "--",
    "echo",
    "hi",
]


@pytest.fixture
def umask_022():
    old = os.umask(0o022)
    yield
    os.umask(old)


def _modes(root):
    for dirpath, dirnames, filenames in os.walk(root):
        for name in dirnames + filenames:
            path = os.path.join(dirpath, name)
            yield path, stat.S_IMODE(os.lstat(path).st_mode)


def _all_text(root):
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            with open(os.path.join(dirpath, name), "rb") as f:
                yield f.read()


def test_redact_argv_keeps_names_and_command():
    assert redact_argv(ARGV) == [
        "qxub",
        "exec",
        "--env",
        "dt3",
        "--var",
        "API_TOKEN=***",
        "--var=OTHER=***",
        "--vars",
        "A=***,B=***",
        "--",
        "echo",
        "hi",
    ]
    # After "--" it is the job's own command
    assert redact_argv(["qxub", "--", "x", "--var", "K=v"])[-1] == "K=v"
    assert redact_argv(["qxub", "--var", "BARE"]) == ["qxub", "--var", "BARE"]


def test_redact_params():
    params = {"var": (f"K={SECRET}",), "vars": f"A={SECRET}", "env": "dt3"}
    assert redact_params(params) == {"var": ["K=***"], "vars": "A=***", "env": "dt3"}


def test_records_are_owner_only_and_hold_no_var_values(
    tmp_path, monkeypatch, umask_022
):
    monkeypatch.setattr(sys, "argv", ARGV)
    config_dir = tmp_path / "qxub"
    h = HistoryManager(config_dir=config_dir)
    ctx = SimpleNamespace(
        params={"var": (f"API_TOKEN={SECRET}",), "vars": f"A={SECRET}"}
    )
    ts = h.log_execution(ctx, job_id=None)
    assert h.update_execution_with_resources(ts, "1.pbs", {"exit_status": 0})

    history = config_dir / "history"
    modes = dict(_modes(history))
    assert modes  # records were written
    assert stat.S_IMODE(history.stat().st_mode) & 0o777 == 0o700
    for path, mode in modes.items():
        assert not mode & 0o077, (path, oct(mode))
    for content in _all_text(config_dir):
        assert SECRET.encode() not in content

    record = h.get_executions(limit=1)[0]
    assert "API_TOKEN=***" in record["context"]["command_line"]
    assert "echo hi" in record["context"]["command_line"]


def test_existing_loose_history_is_tightened(tmp_path, monkeypatch, umask_022):
    monkeypatch.setattr(sys, "argv", ["qxub", "exec", "--env", "dt3", "--", "true"])
    config_dir = tmp_path / "qxub"
    executions = config_dir / "history" / "executions"
    executions.mkdir(parents=True, mode=0o755)
    os.chmod(config_dir / "history", 0o755)
    old = executions / "1.json"
    old.write_text('{"recipe_hash": "x"}')
    os.chmod(old, 0o644)
    legacy = config_dir / "executions.yaml.imported"
    legacy.write_text("executions: {}\n")
    os.chmod(legacy, 0o644)

    h = HistoryManager(config_dir=config_dir)
    h.get_executions()  # reading is first use too

    for path, mode in _modes(config_dir / "history"):
        assert not mode & 0o077, (path, oct(mode))
    assert stat.S_IMODE(old.stat().st_mode) & 0o777 == 0o600
    assert stat.S_IMODE(legacy.stat().st_mode) & 0o777 == 0o600


def test_legacy_import_output_is_owner_only(tmp_path, monkeypatch, umask_022):
    config_dir = tmp_path / "qxub"
    config_dir.mkdir()
    (config_dir / "executions.yaml").write_text(
        "executions:\n  '123':\n    recipe_hash: abc\n"
    )
    h = HistoryManager(config_dir=config_dir)
    assert h.get_executions()
    for path, mode in _modes(config_dir / "history"):
        assert not mode & 0o077, (path, oct(mode))
    imported = config_dir / "executions.yaml.imported"
    assert stat.S_IMODE(imported.stat().st_mode) & 0o777 == 0o600


def test_job_db_is_owner_only(tmp_path, monkeypatch, umask_022):
    from qxub.queue import db

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("QXUB_DB_PATH", raising=False)
    monkeypatch.setattr(db, "_PRIVATE_CHECKED", set())

    path = db.get_db_path()
    assert path == home / ".config" / "qxub" / "qxub.db"
    assert stat.S_IMODE(path.parent.stat().st_mode) & 0o777 == 0o700
    db.init_db(path)
    assert stat.S_IMODE(path.stat().st_mode) & 0o777 == 0o600

    # A database left readable by older qxub is tightened on first use
    os.chmod(path, 0o644)
    monkeypatch.setattr(db, "_PRIVATE_CHECKED", set())
    db.get_db_path()
    assert stat.S_IMODE(path.stat().st_mode) & 0o777 == 0o600
