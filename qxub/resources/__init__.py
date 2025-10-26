"""
qxub.resources - Resource Management Package

This package provides comprehensive resource management utilities for PBS job scheduling,
including parsing, validation, tracking, and conversion between different resource formats.

Package Structure:
    parser.py     - Parse resource specifications (memory, walltime, etc.)
    utils.py      - Resource formatting and manipulation utilities
    tracker.py    - Resource usage tracking and efficiency analysis
    mappers.py    - Convert between workflow engine and PBS resource formats

Public API:
    # Resource Parsing
    parse_memory_size(spec) -> int        # Parse memory specs like "8GB" to bytes
    parse_walltime(spec) -> float         # Parse walltime specs to hours

    # Resource Formatting
    format_walltime(hours) -> str         # Format hours as HH:MM:SS
    bytes_to_human(bytes) -> str          # Format bytes as human-readable

    # Resource Tracking
    ResourceTracker                       # Track job resource usage

    # Resource Mapping (Workflow Integration)
    ResourceMapper                        # Convert workflow resources to PBS

Example Usage:
    from qxub.resources import parse_memory_size, ResourceMapper

    # Parse memory specification
    memory_bytes = parse_memory_size("8GB")  # Returns 8589934592

    # Convert workflow resources to PBS
    mapper = ResourceMapper()
    mapper.add_memory("8GB")
    mapper.add_runtime("2h30m")
    pbs_resources = mapper.get_pbs_resources()  # ["mem=8GB", "walltime=2:30:00"]
"""

from .mappers import ResourceMapper
from .parser import (
    bytes_to_human,
    calculate_efficiency,
    parse_exec_host,
    parse_exec_vnode,
    parse_joblog_resources,
    parse_timestamp,
    seconds_to_time,
    size_to_bytes,
    time_to_seconds,
)
from .tracker import ResourceTracker, resource_tracker

# Import all public APIs for convenient access
from .utils import (
    compare_memory,
    compare_walltime,
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
    "compare_memory",
    "compare_walltime",
    "evaluate_condition",
    "suggest_resource_adjustment",
    # Parser functions (for advanced usage)
    "size_to_bytes",
    "bytes_to_human",
    "time_to_seconds",
    "seconds_to_time",
    "parse_exec_host",
    "parse_exec_vnode",
    "parse_timestamp",
    "calculate_efficiency",
    "parse_joblog_resources",
    # Classes
    "ResourceTracker",
    "ResourceMapper",
    # Global instances
    "resource_tracker",
]
