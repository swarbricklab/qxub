"""Tests for Platform.select_queue auto-selection rule ordering."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from qxub.platforms.core import (  # noqa: E402
    AutoSelectionRule,
    Platform,
    PlatformLoader,
    Queue,
    QueueLimits,
)

GADI_PLATFORM_DIR = Path(__file__).parent.parent / "docs" / "platforms"


@pytest.fixture(scope="module")
def gadi():
    loader = PlatformLoader(search_paths=[GADI_PLATFORM_DIR])
    platform = loader.get_platform("nci_gadi")
    assert platform is not None
    return platform


@pytest.mark.parametrize(
    "resources, expected",
    [
        # The default rule ("cpus <= 48 -> normal") sits before
        # "cpus <= 3200 -> normalsl" and must win for small jobs.
        ({"cpus": 4, "memory": "8GB", "walltime": "1:00:00"}, "normal"),
        ({"cpus": 48, "memory": "190GB"}, "normal"),
        ({"cpus": 1000, "memory": "100GB"}, "normalsl"),
        ({"cpus": 2000, "memory": "100GB"}, "normalbw"),
        ({"cpus": 4000, "memory": "100GB"}, "normalsr"),
        # Rules listed before the default still take precedence over it.
        ({"cpus": 4, "memory": "600GB"}, "hugemem"),
        ({"cpus": 12, "gpus": 1}, "gpuvolta"),
    ],
)
def test_gadi_rules_are_evaluated_in_order(gadi, resources, expected):
    assert gadi.select_queue(resources) == expected


def _queue(name, internet=False):
    return Queue(
        name=name,
        type="standard",
        limits=QueueLimits(),
        internet_connectivity=internet,
    )


def _platform(rules, queues=("small", "big", "fallback")):
    return Platform(
        name="test",
        type="pbs_pro",
        host="test",
        queues={q: _queue(q) for q in queues},
        auto_selection_rules=rules,
    )


def test_default_rule_condition_is_honoured_in_position():
    platform = _platform(
        [
            AutoSelectionRule("cpus <= 8", "small", is_default=True),
            AutoSelectionRule("cpus <= 100", "big"),
        ]
    )
    assert platform.select_queue({"cpus": 4}) == "small"
    assert platform.select_queue({"cpus": 50}) == "big"


def test_default_rule_is_fallback_when_nothing_matches():
    platform = _platform(
        [
            AutoSelectionRule("cpus <= 8", "small", is_default=True),
            AutoSelectionRule("cpus <= 100", "big"),
        ]
    )
    # Matches neither condition: fall back to the default rule's queue.
    assert platform.select_queue({"cpus": 500}) == "small"


def test_default_rule_without_condition_is_only_a_fallback():
    platform = _platform(
        [
            AutoSelectionRule("", "fallback", is_default=True),
            AutoSelectionRule("cpus > 8", "big"),
        ]
    )
    assert platform.select_queue({"cpus": 50}) == "big"
    assert platform.select_queue({"cpus": 4}) == "fallback"


def test_internet_selection_honours_default_rule_condition():
    platform = Platform(
        name="test",
        type="pbs_pro",
        host="test",
        queues={
            "copy": _queue("copy", internet=True),
            "copybig": _queue("copybig", internet=True),
            "compute": _queue("compute"),
        },
        auto_selection_rules=[
            AutoSelectionRule("cpus <= 1", "copy", is_default=True),
            AutoSelectionRule("cpus <= 8", "copybig"),
            AutoSelectionRule("cpus <= 100", "compute"),
        ],
    )
    assert platform.select_queue({"cpus": 1, "internet": True}) == "copy"
    assert platform.select_queue({"cpus": 4, "internet": True}) == "copybig"
    # "compute" matches but lacks internet; fall back to the default rule.
    assert platform.select_queue({"cpus": 50, "internet": True}) == "copy"
