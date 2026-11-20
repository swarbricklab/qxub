"""Retiring the shared live DB (Phase 4 of docs/dev/sqlite-robustness.md).

The job DB is per user again; ``qxub db recover`` imports the old team DB
into the snapshot; ``qxub resources`` reports read the snapshot.
"""

import shutil
import sqlite3
from datetime import datetime
from pathlib import Path

import pytest
from click.testing import CliRunner

from qxub.queue import consolidate as cons
from qxub.queue import db as qdb
from qxub.queue import recover as rec
from qxub.resources.tracker import create_job_resources_schema

NOW = datetime(2026, 11, 20, 12, 0, 0)


def old_db(path: Path, rows) -> Path:
    """An old-style shared DB with ``job_resources`` rows."""
    conn = sqlite3.connect(path)
    create_job_resources_schema(conn)
    for row in rows:
        names = ", ".join(row)
        conn.execute(
            f"INSERT INTO job_resources ({names}) VALUES ({', '.join('?' * len(row))})",
            list(row.values()),
        )
    conn.commit()
    conn.close()
    return path


def job(job_id, status="completed", **extra):
    return {
        "job_id": job_id,
        "timestamp": "2026-06-01T10:00:00",
        "status": status,
        "username": "abc123",
        **extra,
    }


def snapshot_jobs(snapshot):
    conn = cons.open_snapshot(snapshot)
    try:
        return {
            r["job_id"]: dict(r) for r in conn.execute("SELECT * FROM job_resources")
        }
    finally:
        conn.close()


@pytest.fixture
def config_value(monkeypatch):
    """Set config values seen by ``config_manager.get_config_value``."""
    from qxub.config import config_manager

    values = {}
    monkeypatch.setattr(
        config_manager, "get_config_value", lambda key, *a, **k: values.get(key)
    )
    return values


# -- the job DB is per user ----------------------------------------------------


def test_job_db_ignores_shared_db_path(tmp_path, monkeypatch, config_value):
    monkeypatch.delenv("QXUB_DB_PATH", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("QXUB_SHARED_DB", str(tmp_path / "team" / "env.db"))
    config_value["shared_db.path"] = str(tmp_path / "team" / "config.db")

    assert qdb.get_db_path() == tmp_path / ".config" / "qxub" / "qxub.db"
    assert qdb.legacy_shared_db_path() == tmp_path / "team" / "env.db"
    monkeypatch.delenv("QXUB_SHARED_DB")
    assert qdb.legacy_shared_db_path() == tmp_path / "team" / "config.db"


def test_snapshot_still_defaults_beside_the_old_team_db(
    tmp_path, monkeypatch, config_value
):
    monkeypatch.delenv("QXUB_STATS_DB", raising=False)
    monkeypatch.delenv("QXUB_SHARED_DB", raising=False)
    config_value["shared_db.path"] = str(tmp_path / "team" / "qxub.db")
    assert cons.stats_db_path() == tmp_path / "team" / "stats" / "qxub-stats.db"


# -- recover ---------------------------------------------------------------------


def test_recover_adds_finished_jobs_once(tmp_path):
    source = old_db(
        tmp_path / "old.db",
        [
            job("1.gadi-pbs", mem_used_mb=100, mem_efficiency=50.0),
            job("2.gadi-pbs", status="failed"),
            job("3.gadi-pbs", status="submitted"),
            job("4.gadi-pbs", status="running", mem_used_mb=10),
        ],
    )
    snapshot = tmp_path / "stats" / "qxub-stats.db"

    result = rec.recover(source, snapshot, now=NOW)
    assert not result.errors and result.published
    assert result.method == "read directly"
    assert (result.read, result.added, result.unfinished) == (4, 3, 1)
    jobs = snapshot_jobs(snapshot)
    assert set(jobs) == {"1.gadi-pbs", "2.gadi-pbs", "4.gadi-pbs"}
    assert jobs["1.gadi-pbs"]["mem_efficiency"] == 50.0
    assert "old.db" in rec.recovered_from(snapshot)

    again = rec.recover(source, snapshot, now=NOW)
    assert (again.added, again.present, again.published) == (0, 3, False)


def test_recover_keeps_what_the_snapshot_already_has(tmp_path):
    snapshot = tmp_path / "stats" / "qxub-stats.db"
    newer = old_db(tmp_path / "newer.db", [job("1.gadi-pbs", mem_used_mb=999)])
    rec.recover(newer, snapshot, now=NOW)
    source = old_db(tmp_path / "old.db", [job("1.gadi-pbs", mem_used_mb=1)])

    result = rec.recover(source, snapshot, now=NOW)
    assert (result.added, result.present) == (0, 1)
    assert snapshot_jobs(snapshot)["1.gadi-pbs"]["mem_used_mb"] == 999


def test_recover_never_writes_the_source(tmp_path):
    source = old_db(tmp_path / "old.db", [job("1.gadi-pbs")])
    journal = source.with_name("old.db-journal")
    journal.write_bytes(b"not really a journal")
    before = (source.read_bytes(), journal.read_bytes())

    rec.recover(source, tmp_path / "stats" / "qxub-stats.db", now=NOW)
    assert (source.read_bytes(), journal.read_bytes()) == before


def test_recover_dry_run_publishes_nothing(tmp_path):
    source = old_db(tmp_path / "old.db", [job("1.gadi-pbs")])
    snapshot = tmp_path / "stats" / "qxub-stats.db"
    result = rec.recover(source, snapshot, now=NOW, dry_run=True)
    assert result.added == 1 and not result.published
    assert not snapshot.exists()


def corrupt(path: Path) -> None:
    """Overwrite a page in the middle of the table with garbage."""
    data = bytearray(path.read_bytes())
    page = 4096
    pages = len(data) // page
    assert pages > 6
    middle = (pages // 2) * page
    data[middle : middle + page] = b"\xff" * page
    path.write_bytes(bytes(data))


def big_old_db(path: Path) -> Path:
    rows = [job(f"{i}.gadi-pbs", command="x" * 200, mem_used_mb=i) for i in range(400)]
    return old_db(path, rows)


@pytest.mark.skipif(shutil.which("sqlite3") is None, reason="needs the sqlite3 CLI")
def test_recover_rebuilds_a_malformed_db(tmp_path):
    source = big_old_db(tmp_path / "old.db")
    corrupt(source)
    with pytest.raises(sqlite3.DatabaseError):
        rec._read_table(source)

    result = rec.recover(source, tmp_path / "stats" / "qxub-stats.db", now=NOW)
    assert not result.errors, result.errors
    assert result.method == "rebuilt with sqlite3 .recover"
    assert 0 < result.added < 400


def test_recover_without_sqlite3_cli_says_how_to_get_it(tmp_path, monkeypatch):
    source = big_old_db(tmp_path / "old.db")
    corrupt(source)
    monkeypatch.setattr(rec.shutil, "which", lambda name: None)

    result = rec.recover(source, tmp_path / "stats" / "qxub-stats.db", now=NOW)
    assert result.added == 0
    assert "module load sqlite" in result.errors[0]


def test_recover_missing_source_is_skipped(tmp_path):
    result = rec.recover(tmp_path / "nope.db", tmp_path / "s.db", now=NOW)
    assert "does not exist" in result.skipped


# -- the CLI -------------------------------------------------------------------


def test_db_recover_defaults_to_the_old_team_db(tmp_path, monkeypatch):
    from qxub.db_cli import db_cli

    source = old_db(tmp_path / "old.db", [job("1.gadi-pbs")])
    snapshot = tmp_path / "stats" / "qxub-stats.db"
    monkeypatch.setenv("QXUB_SHARED_DB", str(source))
    monkeypatch.setenv("QXUB_STATS_DB", str(snapshot))

    out = CliRunner().invoke(db_cli, ["recover"])
    assert out.exit_code == 0, out.output
    assert "1 added" in out.output
    assert "1.gadi-pbs" in snapshot_jobs(snapshot)

    status = CliRunner().invoke(db_cli, ["status"])
    assert f"Old job DB: {source}" in status.output
    assert "recovered:" in status.output


def test_resources_reports_read_the_snapshot(tmp_path, monkeypatch):
    from qxub.resources_cli import resources

    snapshot = tmp_path / "stats" / "qxub-stats.db"
    monkeypatch.setenv("QXUB_STATS_DB", str(snapshot))

    empty = CliRunner().invoke(resources, ["stats"])
    assert empty.exit_code == 0
    assert "qxub db consolidate" in empty.output

    source = old_db(
        tmp_path / "old.db",
        [
            job("1.gadi-pbs", mem_used_mb=10, mem_efficiency=20.0),
            job("2.gadi-pbs", mem_used_mb=10, mem_efficiency=80.0),
        ],
    )
    rec.recover(source, snapshot, now=NOW)

    stats = CliRunner().invoke(resources, ["stats"])
    assert stats.exit_code == 0, stats.output
    assert "2 jobs" in stats.output
    assert "50.0%" in stats.output

    listed = CliRunner().invoke(resources, ["list", "--csv"])
    assert "1.gadi-pbs" in listed.output and "2.gadi-pbs" in listed.output

    export = tmp_path / "out.csv"
    CliRunner().invoke(resources, ["export", str(export)])
    assert export.read_text().count("gadi-pbs") == 2
