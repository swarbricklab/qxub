"""
Tests for the platform model: queue limits, walltime rules, validation,
platform loading from YAML, queue selection and cost estimation.
"""

from pathlib import Path

import pytest
import yaml

import qxub.platforms.core as platform_core
from qxub.platforms.core import (
    PlatformLoader,
    Queue,
    QueueLimits,
    QueueSelector,
    WalltimeRule,
    select_best_queue,
)


def create_test_platform_file(temp_dir: Path) -> Path:
    """Create a test platform definition file."""
    platform_data = {
        "name": "test_platform",
        "type": "pbs",
        "host": "test.example.com",
        "description": "Test platform for unit testing",
        "queues": [
            {
                "name": "normal",
                "type": "compute",
                "priority": "normal",
                "su_billing_rate": 1.0,
                "limits": {
                    "max_cpus": 48,
                    "min_cpus": 1,
                    "max_memory": "192GB",
                    "min_memory": "1GB",
                },
                "walltime_rules": [{"cores": "1-48", "max_walltime": "48:00:00"}],
                "default_walltime": "1:00:00",
            },
            {
                "name": "express",
                "type": "compute",
                "priority": "high",
                "su_billing_rate": 2.0,
                "limits": {"max_cpus": 12, "min_cpus": 1, "max_memory": "48GB"},
                "walltime_rules": [{"cores": "1-12", "max_walltime": "2:00:00"}],
            },
            {
                "name": "gpu",
                "type": "gpu",
                "priority": "normal",
                "su_billing_rate": 3.0,
                "limits": {
                    "max_cpus": 12,
                    "min_cpus": 1,
                    "max_gpus": 4,
                    "min_gpus": 1,
                    "max_memory": "96GB",
                },
                "walltime_rules": [{"cores": "1-12", "max_walltime": "24:00:00"}],
            },
        ],
        "auto_selection_rules": [
            {"condition": "gpu_requested > 0", "queue": "gpu"},
            {"condition": "walltime <= 2:00:00", "queue": "express"},
            {"condition": "", "queue": "normal", "is_default": True},
        ],
    }

    platform_file = temp_dir / "test_platform.yaml"
    with open(platform_file, "w") as f:
        yaml.dump(platform_data, f, default_flow_style=False)

    return platform_file


@pytest.fixture
def platform_dir(tmp_path):
    create_test_platform_file(tmp_path)
    return tmp_path


@pytest.fixture
def test_platform(platform_dir):
    platform = PlatformLoader(search_paths=[platform_dir]).get_platform("test_platform")
    assert platform is not None
    return platform


def test_queue_limits():
    limits = QueueLimits(max_cpus=48, min_cpus=1, max_memory="192GB", min_memory="1GB")
    assert limits.max_cpus == 48
    assert limits.max_memory == "192GB"


@pytest.mark.parametrize(
    "cores, expected",
    [(24, True), (48, True), (1, True), (96, False)],
)
def test_walltime_rule(cores, expected):
    rule = WalltimeRule(cores="1-48", max_walltime="48:00:00")
    assert rule.matches_core_count(cores) is expected


@pytest.mark.parametrize(
    "resources, should_be_valid",
    [
        ({"cpus": 4, "memory": "8GB", "walltime": "2:00:00"}, True),
        ({"cpus": 1, "memory": "1GB", "walltime": "1:00:00"}, True),
        ({"cpus": 96, "memory": "8GB", "walltime": "2:00:00"}, False),
        ({"cpus": 4, "memory": "256GB", "walltime": "2:00:00"}, False),
        ({"cpus": 4, "memory": "8GB", "walltime": "72:00:00"}, False),
    ],
)
def test_queue_validation(resources, should_be_valid):
    queue = Queue(
        name="test_queue",
        type="compute",
        limits=QueueLimits(max_cpus=48, min_cpus=1, max_memory="192GB"),
        walltime_rules=[WalltimeRule(cores="1-48", max_walltime="48:00:00")],
    )
    result = queue.validate_resources(resources)
    assert result.is_valid is should_be_valid, result.errors


def test_platform_loading(platform_dir):
    loader = PlatformLoader(search_paths=[platform_dir])
    assert "test_platform" in loader.list_platforms()
    assert len(loader.get_platform("test_platform").queues) == 3


@pytest.mark.parametrize(
    "resources, expected_queue",
    [
        ({"cpus": 2, "memory": "8GB", "gpus": 1}, "gpu"),
        ({"cpus": 4, "memory": "8GB", "walltime": "1:00:00"}, "express"),
        ({"cpus": 8, "memory": "16GB", "walltime": "12:00:00"}, "normal"),
    ],
)
def test_queue_selection(test_platform, resources, expected_queue):
    result = QueueSelector(test_platform).select_queue(resources)
    assert result.best_queue == expected_queue, result.warnings


def test_cost_estimation(test_platform):
    # 4 CPUs x 4 hours on the normal queue (1.0 SU per CPU-hour)
    resources = {"cpus": 4, "memory": "8GB", "walltime": "4:00:00"}
    result = QueueSelector(test_platform).select_queue(resources)
    assert result.best_queue == "normal"
    assert result.estimated_cost == pytest.approx(4 * 4.0 * 1.0)


def test_high_level_function(platform_dir, monkeypatch):
    monkeypatch.setattr(
        platform_core,
        "_platform_loader",
        PlatformLoader(search_paths=[platform_dir]),
    )
    resources = {"cpus": 4, "memory": "8GB", "walltime": "4:00:00"}
    result = select_best_queue(resources, platform_name="test_platform")
    assert result.best_queue, result.warnings
