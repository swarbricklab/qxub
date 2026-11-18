"""``qxi`` resolves queue ``auto`` instead of passing it to ``qsub`` (#124)."""

import pytest
from conftest import GADI_PLATFORM_DIR


def gadi_config(sandbox, queue_defaults: str = "") -> None:
    """The bundled Gadi platform, plus extra lines under ``defaults``."""
    sandbox.env["QXUB_PLATFORM_PATHS"] = str(GADI_PLATFORM_DIR)
    sandbox.write_config(
        "defaults:\n"
        "  platform: nci_gadi\n"
        f"{queue_defaults}"
        "platforms:\n"
        "  nci_gadi:\n"
        f"    definition: {GADI_PLATFORM_DIR / 'nci_gadi.yaml'}\n"
    )


def qsub_queue(sandbox, *args) -> str:
    """The ``-q`` value of the ``qsub -I`` command ``qxi --dry-run`` prints."""
    result = sandbox.run("interactive", "--dry-run", "--runtime", "15m", *args)
    assert result.returncode == 0, result.stdout + result.stderr
    (line,) = [
        l for l in result.stdout.splitlines() if l.strip().startswith("qsub -I ")
    ]
    argv = line.split()
    assert "auto" not in argv, line
    return argv[argv.index("-q") + 1]


def test_auto_from_config_is_resolved(sandbox):
    gadi_config(sandbox, "  queue: auto\n")

    assert qsub_queue(sandbox, "--mem", "4GB", "--cpus", "1") == "normal"


def test_auto_picks_hugemem_for_large_memory(sandbox):
    gadi_config(sandbox)

    queue = qsub_queue(sandbox, "--queue", "auto", "--mem", "600GB", "--cpus", "4")

    assert queue == "hugemem"


def test_auto_with_internet_picks_an_internet_queue(sandbox):
    gadi_config(sandbox)

    queue = qsub_queue(sandbox, "--queue", "auto", "--internet", "--mem", "4GB")

    assert queue == "copyq"


@pytest.mark.parametrize(
    "queue_defaults, args, expected",
    [
        ("  queue: auto\n  interactive:\n    queue: express\n", [], "express"),
        ("  queue: auto\n", ["--queue", "normalbw"], "normalbw"),
    ],
    ids=["interactive_default", "explicit_queue"],
)
def test_named_queues_are_passed_through(sandbox, queue_defaults, args, expected):
    gadi_config(sandbox, queue_defaults)

    assert qsub_queue(sandbox, "--mem", "4GB", *args) == expected


def test_auto_without_a_platform_falls_back_to_normal(sandbox):
    sandbox.env["QXUB_PLATFORM_PATHS"] = str(sandbox.root / "no-platforms")
    sandbox.write_config("defaults:\n  queue: auto\n")

    result = sandbox.run("interactive", "--dry-run", "--runtime", "15m", "--mem", "4GB")

    assert "-q normal " in result.stdout
    assert "No platforms available" in result.stderr
    assert "using queue 'normal'" in result.stderr
