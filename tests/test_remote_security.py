"""The ssh hop must not turn config or job values into ssh options or leaks."""

import base64
import json
import re

import pytest
from conftest import REPO_ROOT

from qxub.execution.context import ExecutionContext
from qxub.remote.command_builder import build_remote_command, job_vars_payload
from qxub.remote.platform_executor import PlatformRemoteExecutor, RemoteExecutionError

LOCAL_PLATFORM = REPO_ROOT / "docs" / "platforms" / "local" / "local.yaml"


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


# -- job variables stay off command lines ----------------------------------------

SECRET = "s3cret-marker"


def _remote_command(options):
    ctx = ExecutionContext("default", None, "default")
    return build_remote_command(ctx, options, ["env"])


def test_variable_values_are_not_in_the_remote_command():
    options = {"var": (f"TOKEN={SECRET}",), "vars": f"A=1,B={SECRET}2", "verbose": 3}
    command = _remote_command(options)
    assert SECRET not in command
    assert "--vars-stdin" in command
    ssh = _executor("gadi")._build_ssh_command(command)
    assert not any(SECRET in arg for arg in ssh)
    assert json.loads(job_vars_payload(options)) == [
        f"TOKEN={SECRET}",
        "A=1",
        f"B={SECRET}2",
    ]


def test_no_variables_means_no_stdin_payload():
    assert "--vars-stdin" not in _remote_command({})
    assert job_vars_payload({}) is None


@pytest.mark.parametrize("stream", [True, False])
def test_the_payload_reaches_the_remote_stdin(monkeypatch, tmp_path, capfd, stream):
    received = tmp_path / "stdin"
    executor = _executor("gadi")
    monkeypatch.setattr(
        executor,
        "_build_ssh_command",
        lambda command: ["bash", "-c", f"cat > {received}"],
    )
    payload = json.dumps([f"TOKEN={SECRET}"]) + "\n"
    command = _remote_command({"var": (f"TOKEN={SECRET}",)})

    assert (
        executor.execute(command, stream_output=stream, verbose=3, stdin_data=payload)
        == 0
    )
    assert received.read_text() == payload
    out, err = capfd.readouterr()
    assert SECRET not in out + err


def test_remote_qxub_reads_the_variables_from_stdin(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    payload = json.dumps([f"TOKEN={SECRET}", "OTHER=2"]) + "\n"
    result = sandbox.run(
        "exec", "--dry-run", "--vars-stdin", "--var", "CLI=1", "--", "env",
        input=payload,
    )  # fmt: skip
    assert result.returncode == 0, result.stderr
    encoded = re.search(r"user_vars_b64='?([A-Za-z0-9+/=]+)", result.stdout).group(1)
    assert base64.b64decode(encoded).decode() == f"CLI=1;TOKEN={SECRET};OTHER=2"


def test_remote_qxub_rejects_a_bad_stdin_payload(sandbox):
    sandbox.env["QXUB_PLATFORM_FILE"] = str(LOCAL_PLATFORM)
    result = sandbox.run(
        "exec", "--dry-run", "--vars-stdin", "--", "env", input='{"a": 1}\n'
    )
    assert result.returncode != 0
    assert "JSON list" in result.stderr
