"""Config loading must not trust files that other people can plant or change."""

import os
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


# -- user files hold secrets ---------------------------------------------------


def _mode(path):
    return path.stat().st_mode & 0o777


@pytest.fixture
def loose_umask():
    old = os.umask(0o022)
    yield
    os.umask(old)


def test_user_config_is_written_owner_only(home, loose_umask):
    cm = ConfigManager()
    cm.set_user_config_value("notifications.slack.bot_token", "xoxb-secret")
    config = home / "home" / ".config" / "qxub" / "config.yaml"
    assert _mode(config) == 0o600
    assert _mode(config.parent) == 0o700


def test_existing_loose_user_config_is_tightened(home, loose_umask):
    config_dir = home / "home" / ".config" / "qxub"
    config_dir.mkdir()
    config = config_dir / "config.yaml"
    config.write_text("defaults:\n  queue: normal\n")
    config.chmod(0o644)

    ConfigManager().set_user_config_value("defaults.name", "x")
    assert _mode(config) == 0o600
    assert not list(config_dir.glob(".config.yaml.*")), "no temp file left"


def test_user_config_template_is_owner_only(home, loose_umask):
    path = ConfigManager().create_user_config_template()
    assert _mode(path) == 0o600


def test_saving_an_alias_keeps_interpolations(home, loose_umask, monkeypatch):
    monkeypatch.setenv("QXUB_TEST_TOKEN", "s3cret")
    config_dir = home / "home" / ".config" / "qxub"
    config_dir.mkdir()
    config = config_dir / "config.yaml"
    config.write_text(
        "notifications:\n  slack:\n    bot_token: ${oc.env:QXUB_TEST_TOKEN}\n"
    )

    cm = ConfigManager()
    cm.save_alias("quick", {"queue": "express"})

    text = config.read_text()
    assert "${oc.env:QXUB_TEST_TOKEN}" in text
    assert "s3cret" not in text
    assert _mode(config) == 0o600
    assert cm.get_alias("quick") == {"queue": "express"}


def test_saving_an_alias_does_not_clobber_an_unreadable_config(home):
    config_dir = home / "home" / ".config" / "qxub"
    config_dir.mkdir()
    config = config_dir / "config.yaml"
    config.write_text("defaults: [unclosed\n")

    cm = ConfigManager()
    with pytest.raises(Exception):
        cm.save_alias("quick", {"queue": "express"})
    assert config.read_text() == "defaults: [unclosed\n"


def test_user_shortcuts_are_written_owner_only(home, loose_umask):
    sm = ShortcutManager()
    sm.add_shortcut("deploy", {"vars": ["TOKEN=s3cret"]})
    path = home / "home" / ".config" / "qxub" / "shortcuts.json"
    assert _mode(path) == 0o600
    path.chmod(0o644)
    sm.add_shortcut("other", {"env": "base"})
    assert _mode(path) == 0o600
    path.chmod(0o644)
    assert sm.remove_shortcut("other")
    assert _mode(path) == 0o600


def test_a_symlinked_user_config_keeps_its_link(home, loose_umask):
    config_dir = home / "home" / ".config" / "qxub"
    config_dir.mkdir()
    target = home / "dotfiles.yaml"
    target.write_text("defaults:\n  queue: normal\n")
    (config_dir / "config.yaml").symlink_to(target)

    ConfigManager().set_user_config_value("defaults.name", "x")
    assert (config_dir / "config.yaml").is_symlink()
    assert "name: x" in target.read_text()
    assert _mode(target) == 0o600
