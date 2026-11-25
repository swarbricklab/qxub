"""Config loading must not trust files that other people can plant or change."""

from pathlib import Path

import pytest

from qxub.config.files import xdg_config_dirs
from qxub.config.manager import ConfigManager
from qxub.config.shortcuts import ShortcutManager


@pytest.fixture
def home(tmp_path, monkeypatch):
    """A private home and XDG layout, with the working directory in ``work``."""
    home = tmp_path / "home"
    (home / ".config").mkdir(parents=True)
    xdg = tmp_path / "xdg"
    xdg.mkdir()
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(home / ".config"))
    monkeypatch.setenv("XDG_CONFIG_DIRS", str(xdg))
    monkeypatch.chdir(work)
    return tmp_path


# -- XDG_CONFIG_DIRS -----------------------------------------------------------


@pytest.mark.parametrize(
    "value, expected",
    [
        ("", ["/etc/xdg"]),
        (":", ["/etc/xdg"]),
        ("relative", ["/etc/xdg"]),
        (":/opt/a::/opt/b:", ["/opt/a", "/opt/b"]),
        ("rel:/opt/a", ["/opt/a"]),
    ],
)
def test_empty_and_relative_xdg_config_dirs_are_dropped(monkeypatch, value, expected):
    monkeypatch.setenv("XDG_CONFIG_DIRS", value)
    assert xdg_config_dirs() == [Path(p) for p in expected]


def test_unset_xdg_config_dirs_is_etc_xdg(monkeypatch):
    monkeypatch.delenv("XDG_CONFIG_DIRS", raising=False)
    assert xdg_config_dirs() == [Path("/etc/xdg")]


@pytest.mark.parametrize("value", ["", ":", "::", "qxub-rel:"])
def test_nothing_is_read_from_the_working_directory(home, monkeypatch, value):
    planted = Path.cwd() / "qxub"
    planted.mkdir()
    (planted / "config.yaml").write_text("defaults:\n  queue: planted\n")
    (planted / "shortcuts.json").write_text('{"ls": {"pre": "evil"}}')
    monkeypatch.setenv("XDG_CONFIG_DIRS", value)

    cm = ConfigManager()
    assert cm.get_config_value("defaults.queue") != "planted"
    assert all(d.is_absolute() for d in cm._get_xdg_config_dirs())

    sm = ShortcutManager()
    assert sm._system_shortcuts_file.is_absolute()
    assert sm.find_shortcut(["ls"]) is None


# -- shared system files -------------------------------------------------------


def test_group_writable_system_config_is_ignored(home, capsys):
    qxub_dir = home / "xdg" / "qxub"
    qxub_dir.mkdir()
    config = qxub_dir / "config.yaml"
    config.write_text("defaults:\n  queue: shared\n")
    assert ConfigManager().get_config_value("defaults.queue") == "shared"

    config.chmod(0o664)
    assert ConfigManager().get_config_value("defaults.queue") is None
    assert "writable by group or others" in capsys.readouterr().err


def test_system_config_in_a_writable_directory_is_ignored(home):
    qxub_dir = home / "xdg" / "qxub"
    qxub_dir.mkdir()
    (qxub_dir / "config.yaml").write_text("defaults:\n  queue: shared\n")
    qxub_dir.chmod(0o777)
    try:
        assert ConfigManager().get_config_value("defaults.queue") is None
    finally:
        qxub_dir.chmod(0o755)


def test_group_writable_system_shortcuts_are_ignored(home):
    qxub_dir = home / "xdg" / "qxub"
    qxub_dir.mkdir()
    shortcuts = qxub_dir / "shortcuts.json"
    shortcuts.write_text('{"ls": {"env": "base"}}')
    assert ShortcutManager().find_shortcut(["ls"]) is not None

    shortcuts.chmod(0o666)
    assert ShortcutManager().find_shortcut(["ls"]) is None
    assert ShortcutManager().list_shortcuts_with_origin() == {}
