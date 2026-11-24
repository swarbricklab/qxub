"""``submit_job`` applies the configured resource defaults, as ``qxub exec`` does.

#157 review: a ``submit_job(env=...)`` job on Gadi ran with only
``storage=scratch/<project>``. ``qxub exec`` added the configured
``defaults.volumes`` (where the conda environments live), but ``submit_job``
skipped ``defaults.mem``, ``defaults.volumes`` and the rest, so conda wasn't
found.
"""

import subprocess
import sys

import pytest

CONFIG = """\
defaults:
  project: ab12
  mem: 8GB
  runtime: 2h
  volumes: gdata/ab12+gdata/cd34
"""


def _submit(sandbox, call):
    script = f"from qxub.execution.submit import submit_job; {call}"
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=sandbox.cwd,
        env=sandbox.env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    ((qsub,),) = [sandbox.pbs_calls("qsub")]
    argv = qsub["argv"]
    return [argv[i + 1] for i, a in enumerate(argv) if a == "-l"]


def _exec(sandbox, *args):
    result = sandbox.run("exec", "--terse", *args, "--", "true")
    assert result.returncode == 0, result.stdout + result.stderr
    argv = sandbox.pbs_calls("qsub")[-1]["argv"]
    return [argv[i + 1] for i, a in enumerate(argv) if a == "-l"]


@pytest.mark.parametrize(
    "call, expected",
    [
        (
            "submit_job('true', env='base')",
            ["mem=8GB", "walltime=02:00:00", "storage=gdata/ab12+gdata/cd34"],
        ),
        # What the caller asks for wins; defaults fill in the rest.
        (
            "submit_job('true', resources={'mem': '4GB', 'storage': 'gdata/xy99'})",
            ["mem=4GB", "walltime=02:00:00", "storage=gdata/xy99"],
        ),
    ],
)
def test_submit_job_applies_configured_defaults(sandbox, call, expected):
    sandbox.write_config(CONFIG)
    assert _submit(sandbox, call) == expected


def test_submit_job_and_exec_request_the_same(sandbox):
    sandbox.write_config(CONFIG)
    from_python = _submit(sandbox, "submit_job('true')")
    assert _exec(sandbox) == from_python
