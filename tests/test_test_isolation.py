"""The test session runs this checkout's qxub, whatever else is installed.

#163 review: on Gadi the session's interpreter could also see another
environment's site-packages with qxub 3.9.1, and the ``qxub`` shim and
sandbox runs imported that instead of the checkout.
"""

import os
import subprocess
import sys

from conftest import REPO_ROOT

import qxub

DECOY_VERSION = "0.0.0-installed-elsewhere"


def _env_with_a_decoy_qxub(tmp_path, base=None):
    """``base`` (default: this session's environment) where the interpreter
    can also see a different qxub, as with an installed copy."""
    decoy = tmp_path / "elsewhere" / "qxub"
    decoy.mkdir(parents=True)
    (decoy / "__init__.py").write_text(f"__version__ = {DECOY_VERSION!r}\n")
    (decoy / "cli.py").write_text(
        f"def qxub():\n    print('qxub', {DECOY_VERSION!r})\n"
    )
    env = dict(base if base is not None else os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        p for p in [env.get("PYTHONPATH", ""), str(decoy.parent)] if p
    )
    return env


def test_qxub_on_path_is_the_checkout(tmp_path):
    env = _env_with_a_decoy_qxub(tmp_path)
    out = subprocess.run(
        ["qxub", "--version"], env=env, capture_output=True, text=True, cwd=tmp_path
    )
    assert out.stdout.split() == ["qxub", qxub.__version__], out.stdout + out.stderr


def test_subprocess_pythons_import_the_checkout(tmp_path):
    env = _env_with_a_decoy_qxub(tmp_path)
    out = subprocess.run(
        [sys.executable, "-c", "import qxub; print(qxub.__file__)"],
        env=env,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert out.stdout.strip().startswith(str(REPO_ROOT / "qxub")), out.stdout


def test_sandbox_runs_import_the_checkout(sandbox):
    sandbox.env = _env_with_a_decoy_qxub(sandbox.root, sandbox.env)
    result = sandbox.run("--version")
    assert result.stdout.split() == ["qxub", qxub.__version__], result.stdout


# -- timing limits on busy machines ------------------------------------------


def test_fastest_of_retries_only_until_under_the_limit():
    from conftest import fastest_of

    timings = iter([3.02, 2.1, 9.9])  # a load spike, then a normal run
    calls = []

    def measure():
        calls.append(1)
        return next(timings)

    assert fastest_of(measure, limit=3.0) == 2.1
    assert len(calls) == 2  # stopped at the first run under the limit

    always_slow = iter([4.0, 3.5, 3.9])
    assert fastest_of(lambda: next(always_slow), limit=3.0) == 3.5  # still fails
