"""
qxub.resources - Resource Management Package

This package provides comprehensive resource management utilities for PBS job scheduling,
including parsing, validation, tracking, and conversion between different resource formats.

Package Structure:
    parser.py     - Parse resource specifications (memory, walltime, etc.)
    utils.py      - Resource formatting and manipulation utilities
    tracker.py    - Resource usage tracking and efficiency analysis
    spec.py       - ResourceSpec, and the one size and duration parser

Public API:
    # Resource Parsing
    parse_memory_size(spec) -> int        # Parse memory specs like "8GB" to bytes
    parse_walltime(spec) -> float         # Parse walltime specs to hours

    # Resource Formatting
    format_walltime(hours) -> str         # Format hours as HH:MM:SS
    bytes_to_human(bytes) -> str          # Format bytes as human-readable

    # Resource Tracking
    ResourceTracker                       # Track job resource usage

    # Resource requests
    ResourceSpec                          # What a job asks for, scheduler-neutral
    parse_size(spec) -> int               # "8GB" -> bytes
    parse_duration(spec) -> int           # "2h30m" -> seconds

Example Usage:
    from qxub.resources import ResourceSpec, parse_memory_size

    # Parse memory specification
    memory_bytes = parse_memory_size("8GB")  # Returns 8589934592

    # Describe a job's resources; the scheduler backend writes its own syntax
    spec = ResourceSpec.from_options(memory="8GB", walltime="2h30m", cpus=4)
"""

from .parser import (
    bytes_to_human,
    calculate_efficiency,
    parse_exec_host,
    parse_exec_vnode,
    parse_joblog_resources,
    parse_timestamp,
    size_to_bytes,
    time_to_seconds,
)
from .spec import ResourceSpec, format_duration, format_size, parse_duration, parse_size
from .tracker import ResourceTracker, resource_tracker

# Import all public APIs for convenient access
from .utils import (
    evaluate_condition,
    format_memory_size,
    format_walltime,
    parse_memory_size,
    parse_walltime,
    suggest_resource_adjustment,
)

# Define what gets imported with "from qxub.resources import *"
__all__ = [
    # Utility functions (most commonly used)
    "parse_memory_size",
    "parse_walltime",
    "format_memory_size",
    "format_walltime",
    "evaluate_condition",
    "suggest_resource_adjustment",
    # Parser functions (for advanced usage)
    "size_to_bytes",
    "bytes_to_human",
    "time_to_seconds",
    "parse_exec_host",
    "parse_exec_vnode",
    "parse_timestamp",
    "calculate_efficiency",
    "parse_joblog_resources",
    # Classes
    "ResourceTracker",
    "ResourceSpec",
    # Resource requests
    "parse_size",
    "parse_duration",
    "format_size",
    "format_duration",
    # Global instances
    "resource_tracker",
]
