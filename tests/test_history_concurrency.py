"""History survives concurrent qxub processes (#129).

``executions.yaml`` used to be rewritten in full by every process, so
parallel ``qx``/``qxi`` runs lost each other's records or wiped the file.
Each record is now its own file.
"""

import subprocess
import sys
import textwrap
from types import SimpleNamespace

import pytest
import yaml

from qxub.history import manager as manager_mod
from qxub.history.manager import HistoryManager

ARGV = ["qxub", "exec", "--env", "dt3", "--", "echo", "hi"]


@pytest.fixture
def history(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ARGV)
    return HistoryManager(config_dir=tmp_path / "qxub")


def ctx(**params):
    return SimpleNamespace(params=params)


WORKER = textwrap.dedent(
    """
    import sys
    from types import SimpleNamespace
    from qxub.history.manager import HistoryManager

    config_dir, worker, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    sys.argv = ["qxub", "exec", "--env", "dt3", "--", "echo", worker]
    h = HistoryManager(config_dir=config_dir)
    for i in range(n):
        ts = h.log_execution(SimpleNamespace(params={"worker": worker, "i": i}))
        assert h.update_execution_with_resources(
            ts, f"{worker}{i}.pbs", {"exit_status": i}
        ), ts
    """
)


def test_parallel_processes_keep_every_record(tmp_path):
    config_dir = tmp_path / "qxub"
    workers, per_worker = 6, 8

    procs = [
        subprocess.Popen(
            [sys.executable, "-c", WORKER, str(config_dir), f"w{w}-", str(per_worker)]
        )
        for w in range(workers)
    ]
    assert [p.wait(timeout=120) for p in procs] == [0] * workers

    executions = HistoryManager(config_dir=config_dir).get_executions(limit=1000)
    job_ids = sorted(e["execution"]["job_id"] for e in executions)
    expected = sorted(
        f"w{w}-{i}.pbs" for w in range(workers) for i in range(per_worker)
    )
    assert job_ids == expected
    assert all(
        e["execution"]["exit_code"] == int(e["execution"]["job_id"].split("-")[1][0])
        for e in executions
    )


def test_same_microsecond_gets_distinct_records(history, monkeypatch):
    monkeypatch.setattr(manager_mod.time, "time", lambda: 1_700_000_000.0)

    first = history.log_execution(ctx(i=1))
    second = history.log_execution(ctx(i=2))

    assert first != second
    assert len(history.get_executions(limit=10)) == 2


def test_update_and_lookup_by_job_id(history):
    ts = history.log_execution(ctx())

    assert history.update_execution_with_resources(
        ts, "123.gadi-pbs", {"exit_status": 7, "resources_used": {"cput": "00:01:00"}}
    )
    record = history.get_execution_by_job_id("123")
    assert record["timestamp"] == ts
    assert record["execution"]["exit_code"] == 7
    assert record["resources"]["used"] == {"cput": "00:01:00"}


def test_update_of_missing_record_is_reported(history, caplog):
    assert not history.update_execution_with_resources("42", "1.pbs", {})
    assert "not found" in caplog.text


def test_half_written_files_are_ignored(history):
    history.log_execution(ctx())
    executions_dir = history.config_dir / "history" / "executions"
    (executions_dir / "99999999999999999.json").write_text("")  # claimed, not written
    (executions_dir / ".1.abc.tmp").write_text('{"partial":')

    assert len(history.get_executions(limit=10)) == 1


def test_recipe_counts_come_from_executions(history):
    for _ in range(3):
        history.log_execution(ctx())

    (recipe,) = history.get_recipes().values()
    assert recipe["executor"] == {"type": "conda", "env": "dt3"}
    assert recipe["metadata"]["run_count"] == 3
    assert recipe["metadata"]["first_seen"]
    assert recipe["metadata"]["last_seen"]


def test_old_records_are_pruned(history, monkeypatch):
    monkeypatch.setattr(manager_mod, "MAX_EXECUTIONS", 5)
    monkeypatch.setattr(manager_mod, "PRUNE_SLACK", 2)
    keys = [history.log_execution(ctx(i=i)) for i in range(8)]

    kept = [e["timestamp"] for e in history.get_executions(limit=100)]
    assert kept == list(reversed(keys[3:]))


def test_yaml_history_is_imported_once(history):
    config_dir = history.config_dir
    (config_dir / "executions.yaml").write_text(
        yaml.safe_dump(
            {
                "executions": {
                    "1700000000000001": {
                        "recipe_hash": "abc123abc123",
                        "execution": {"job_id": "1.gadi-pbs", "status": "completed"},
                    }
                }
            }
        )
    )
    (config_dir / "recipes.yaml").write_text(
        yaml.safe_dump(
            {
                "recipes": {
                    "abc123abc123": {
                        "executor": {"type": "default"},
                        "target": {"cmd": "ls"},
                        "metadata": {
                            "first_seen": "2026-01-01T00:00:00",
                            "last_seen": "2026-01-02T00:00:00",
                            "run_count": 1,
                        },
                    }
                }
            }
        )
    )

    record = history.get_execution_by_job_id("1.gadi-pbs")

    assert record["timestamp"] == "1700000000000001"
    assert history.get_recipe_by_hash("abc123abc123")["metadata"]["run_count"] == 1
    assert not (config_dir / "executions.yaml").exists()
    assert (config_dir / "executions.yaml.imported").exists()

    # An older qxub writes the YAML again: import is idempotent.
    (config_dir / "executions.yaml.imported").rename(config_dir / "executions.yaml")
    history.log_execution(ctx())
    assert len(history.get_executions(limit=10)) == 2


def test_clear_all_history(history):
    history.log_execution(ctx())
    history.clear_all_history()
    assert history.get_executions() == []
    assert history.get_recipes() == {}


def write_legacy_recipes(config_dir, recipes):
    (config_dir / "recipes.yaml").write_text(yaml.safe_dump({"recipes": recipes}))


def test_imported_recipes_keep_run_counts(history):
    """Review of #130: run counts were reset to the retained executions."""
    ours = history.compute_recipe_hash(
        {"executor": {"type": "conda", "env": "dt3"}, "target": {"cmd": "echo hi"}}
    )
    write_legacy_recipes(
        history.config_dir,
        {
            ours: {
                "executor": {"type": "conda", "env": "dt3"},
                "target": {"cmd": "echo hi"},
                "metadata": {
                    "first_seen": "2026-11-12T09:00:00",
                    "last_seen": "2026-11-16T10:00:00",
                    "run_count": 850,
                },
            },
            "ef808e9ff02a": {
                "executor": {"type": "default"},
                "target": {"cmd": "dt-hook-pre-commit"},
                "metadata": {
                    "first_seen": "2026-11-07T09:00:00",
                    "last_seen": "2026-11-17T10:00:00",
                    "run_count": 40,
                },
            },
        },
    )
    # One of the 850 runs is still in executions.yaml: it must not count twice.
    (history.config_dir / "executions.yaml").write_text(
        yaml.safe_dump(
            {
                "executions": {
                    "1763283600000000": {
                        "recipe_hash": ours,
                        "execution": {"timestamp": "2026-11-16T10:00:00"},
                    }
                }
            }
        )
    )

    recipes = history.get_recipes()
    assert recipes[ours]["metadata"]["run_count"] == 850
    assert recipes[ours]["metadata"]["last_seen"] == "2026-11-16T10:00:00"
    assert recipes["ef808e9ff02a"]["metadata"]["run_count"] == 40
    assert recipes["ef808e9ff02a"]["metadata"]["last_seen"] == "2026-11-17T10:00:00"

    history.log_execution(ctx())

    recipes = history.get_recipes()
    assert recipes[ours]["metadata"]["run_count"] == 851
    assert recipes[ours]["metadata"]["last_seen"] > "2026-11-16T10:00:00"
    assert recipes["ef808e9ff02a"]["metadata"]["run_count"] == 40
