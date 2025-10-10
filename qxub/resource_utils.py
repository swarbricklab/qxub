"""
Resource parsing and validation utilities for qxub platform system.

Handles parsing of resource specifications like memory sizes, walltime durations,
and resource requirement expressions for queue selection and validation.
"""

import logging
import re
from typing import Optional, Union

logger = logging.getLogger(__name__)


def parse_memory_size(memory_str: str) -> Optional[int]:
    """
    Parse memory size string to bytes.

    Supports formats like: "4GB", "512MB", "2TB", "1024" (assumed MB)
    Returns memory size in bytes, or None if parsing fails.
    """
    if not memory_str:
        return None

    memory_str = memory_str.strip().upper()

    # Extract number and unit
    match = re.match(r"^(\d+(?:\.\d+)?)\s*([KMGT]?B?)$", memory_str)
    if not match:
        # Try just number (assume MB)
        try:
            return int(float(memory_str)) * 1024 * 1024
        except ValueError:
            logger.warning(f"Could not parse memory size: {memory_str}")
            return None

    value, unit = match.groups()
    value = float(value)

    # Convert to bytes
    multipliers = {
        "": 1024 * 1024,  # Default to MB
        "B": 1,
        "KB": 1024,
        "MB": 1024 * 1024,
        "GB": 1024 * 1024 * 1024,
        "TB": 1024 * 1024 * 1024 * 1024,
        "K": 1024,
        "M": 1024 * 1024,
        "G": 1024 * 1024 * 1024,
        "T": 1024 * 1024 * 1024 * 1024,
    }

    multiplier = multipliers.get(unit, 1024 * 1024)  # Default to MB
    return int(value * multiplier)


def format_memory_size(bytes_size: int) -> str:
    """Format bytes to human-readable memory size."""
    if bytes_size < 1024:
        return f"{bytes_size}B"
    elif bytes_size < 1024 * 1024:
        return f"{bytes_size // 1024}KB"
    elif bytes_size < 1024 * 1024 * 1024:
        return f"{bytes_size // (1024 * 1024)}MB"
    elif bytes_size < 1024 * 1024 * 1024 * 1024:
        return f"{bytes_size // (1024 * 1024 * 1024)}GB"
    else:
        return f"{bytes_size // (1024 * 1024 * 1024 * 1024)}TB"


def parse_walltime(walltime_str: str) -> Optional[float]:
    """
    Parse walltime string to hours.

    Supports formats like: "48:00:00", "2:30:00", "1h", "30m", "45s"
    Returns walltime in hours, or None if parsing fails.
    """
    if not walltime_str:
        return None

    walltime_str = walltime_str.strip()

    # Try HH:MM:SS format
    if ":" in walltime_str:
        parts = walltime_str.split(":")
        if len(parts) == 3:
            try:
                hours = int(parts[0])
                minutes = int(parts[1])
                seconds = int(parts[2])
                return hours + minutes / 60 + seconds / 3600
            except ValueError:
                pass
        elif len(parts) == 2:
            try:
                hours = int(parts[0])
                minutes = int(parts[1])
                return hours + minutes / 60
            except ValueError:
                pass

    # Try human-readable format (1h, 30m, 45s)
    match = re.match(r"^(\d+(?:\.\d+)?)\s*([hms]?)$", walltime_str.lower())
    if match:
        value, unit = match.groups()
        value = float(value)

        if unit == "h" or unit == "":
            return value
        elif unit == "m":
            return value / 60
        elif unit == "s":
            return value / 3600

    logger.warning(f"Could not parse walltime: {walltime_str}")
    return None


def format_walltime(hours: float) -> str:
    """Format hours to HH:MM:SS walltime string."""
    total_seconds = int(hours * 3600)
    hours_int = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return f"{hours_int:02d}:{minutes:02d}:{seconds:02d}"


def compare_memory(mem1_str: str, mem2_str: str) -> Optional[int]:
    """
    Compare two memory size strings.

    Returns:
        -1 if mem1 < mem2
         0 if mem1 == mem2
         1 if mem1 > mem2
         None if either cannot be parsed
    """
    mem1_bytes = parse_memory_size(mem1_str)
    mem2_bytes = parse_memory_size(mem2_str)

    if mem1_bytes is None or mem2_bytes is None:
        return None

    if mem1_bytes < mem2_bytes:
        return -1
    elif mem1_bytes > mem2_bytes:
        return 1
    else:
        return 0


def compare_walltime(time1_str: str, time2_str: str) -> Optional[int]:
    """
    Compare two walltime strings.

    Returns:
        -1 if time1 < time2
         0 if time1 == time2
         1 if time1 > time2
         None if either cannot be parsed
    """
    time1_hours = parse_walltime(time1_str)
    time2_hours = parse_walltime(time2_str)

    if time1_hours is None or time2_hours is None:
        return None

    if time1_hours < time2_hours:
        return -1
    elif time1_hours > time2_hours:
        return 1
    else:
        return 0


def evaluate_condition(condition: str, resources: dict) -> bool:
    """
    Evaluate a resource condition for queue selection.

    Supports conditions like:
    - "gpu_requested > 0"
    - "memory > 192GB"
    - "cpus <= 48"
    - "walltime > 10:00:00"

    Args:
        condition: Boolean condition string
        resources: Resource requirements dict

    Returns:
        True if condition is met, False otherwise
    """
    if not condition:
        return False

    try:
        # Extract variable, operator, and value
        match = re.match(r"(\w+)\s*([><=!]+)\s*(.+)", condition.strip())
        if not match:
            logger.warning(f"Could not parse condition: {condition}")
            return False

        variable, operator, value_str = match.groups()
        value_str = value_str.strip()

        # Get resource value
        if variable == "gpu_requested":
            resource_value = resources.get("gpus", 0)
        elif variable == "cpus":
            resource_value = resources.get("cpus", 1)
        elif variable == "memory":
            memory_str = resources.get("memory", "0")
            resource_value = parse_memory_size(memory_str)
            comparison_value = parse_memory_size(value_str)
            if resource_value is None or comparison_value is None:
                return False
        elif variable == "walltime":
            walltime_str = resources.get("walltime", "1:00:00")
            resource_value = parse_walltime(walltime_str)
            comparison_value = parse_walltime(value_str)
            if resource_value is None or comparison_value is None:
                return False
        else:
            # Unknown variable
            return False

        # For memory and walltime, we already parsed both values
        if variable in ["memory", "walltime"]:
            comparison_value = comparison_value
        else:
            try:
                comparison_value = float(value_str)
            except ValueError:
                logger.warning(f"Could not parse comparison value: {value_str}")
                return False

        # Apply operator
        if operator == ">":
            return resource_value > comparison_value
        elif operator == ">=":
            return resource_value >= comparison_value
        elif operator == "<":
            return resource_value < comparison_value
        elif operator == "<=":
            return resource_value <= comparison_value
        elif operator == "==" or operator == "=":
            return resource_value == comparison_value
        elif operator == "!=" or operator == "!=":
            return resource_value != comparison_value
        else:
            logger.warning(f"Unknown operator: {operator}")
            return False

    except Exception as e:
        logger.warning(f"Error evaluating condition '{condition}': {e}")
        return False

    return False


def suggest_resource_adjustment(
    resources: dict, queue_limits: dict, policy: str
) -> Optional[dict]:
    """
    Suggest resource adjustments based on queue limits and policy.

    Args:
        resources: Current resource requirements
        queue_limits: Queue limit constraints
        policy: Adjustment policy ("auto", "suggest", "user", "error")

    Returns:
        Dict with suggested adjustments, or None if no changes needed
    """
    suggestions = {}

    if policy == "user":
        return None  # Never suggest changes

    # CPU adjustments
    cpus = resources.get("cpus", 1)
    min_cpus = queue_limits.get("min_cpus")
    max_cpus = queue_limits.get("max_cpus")

    if min_cpus and cpus < min_cpus:
        suggestions["cpus"] = min_cpus
        suggestions["cpu_reason"] = f"Increased from {cpus} to meet queue minimum"
    elif max_cpus and cpus > max_cpus:
        suggestions["cpus"] = max_cpus
        suggestions["cpu_reason"] = f"Reduced from {cpus} to meet queue maximum"

    # Memory adjustments
    memory_str = resources.get("memory")
    if memory_str:
        memory_bytes = parse_memory_size(memory_str)
        max_memory_str = queue_limits.get("max_memory")
        min_memory_str = queue_limits.get("min_memory")

        if max_memory_str:
            max_memory_bytes = parse_memory_size(max_memory_str)
            if memory_bytes and max_memory_bytes and memory_bytes > max_memory_bytes:
                suggestions["memory"] = max_memory_str
                suggestions["memory_reason"] = (
                    f"Reduced from {memory_str} to meet queue maximum"
                )

        if min_memory_str:
            min_memory_bytes = parse_memory_size(min_memory_str)
            if memory_bytes and min_memory_bytes and memory_bytes < min_memory_bytes:
                suggestions["memory"] = min_memory_str
                suggestions["memory_reason"] = (
                    f"Increased from {memory_str} to meet queue minimum"
                )

    return suggestions if suggestions else None
