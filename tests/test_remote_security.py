"""The ssh hop must not turn config or job values into ssh options or leaks."""

import pytest

from qxub.remote.platform_executor import PlatformRemoteExecutor, RemoteExecutionError


def _executor(host):
    return PlatformRemoteExecutor("p", {"host": host, "working_dir": "/w"})


@pytest.mark.parametrize(
    "host",
    [
        "gadi",
        "gadi.nci.org.au",
        "jr9959@gadi.nci.org.au",
        "first.last@login-1",
        "10.0.0.1",
        "user@[::1]",
        "my_host",
    ],
)
def test_ordinary_destinations_are_accepted(host):
    ssh = _executor(host)._build_ssh_command("qxub exec -- true")
    # The destination comes straight after "--", so ssh can't take it as
    # an option
    assert ssh[ssh.index("--") + 1] == host
    assert ssh.index("--") == len(ssh) - 3


@pytest.mark.parametrize(
    "host",
    [
        "-oProxyCommand=touch /tmp/pwned",
        "-F/tmp/evil",
        "--",
        "gadi -oProxyCommand=x",
        "gadi\n-oProxyCommand=x",
        "gadi\t",
        "user@-oProxyCommand=x",
        "-user@gadi",
        ["gadi"],
        123,
    ],
)
def test_option_like_or_malformed_destinations_are_rejected(host):
    with pytest.raises(RemoteExecutionError, match="not a valid ssh destination"):
        _executor(host)
