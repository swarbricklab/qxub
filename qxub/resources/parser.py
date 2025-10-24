"""
Resource parsing utilities for PBS job monitoring.

This module provides functions to parse PBS resource formats from qstat -fx output:
- Size formats: '4294967296b', '135580kb', '104857600b'
- Time formats: '00:30:00', '00:00:04'
- Host formats: 'gadi-cpu-clx-2841/3', '(hostname:ncpus=1:mem=4194304kb:jobfs=102400kb)'
"""

import re
from datetime import datetime
from typing import Any, Dict, Optional


def size_to_bytes(size_str: str) -> int:
    """
    Convert PBS size string to bytes.

    Args:
        size_str: Size string like '4294967296b', '135580kb', '100mb', '4gb'

    Returns:
        Size in bytes as integer

    Examples:
        >>> size_to_bytes('4294967296b')
        4294967296
        >>> size_to_bytes('135580kb')
        138833920
        >>> size_to_bytes('100mb')
        104857600
        >>> size_to_bytes('4gb')
        4294967296
    """
    if not size_str or size_str == "0b":
        return 0

    # Extract number and unit
    match = re.match(r"^(\d+(?:\.\d+)?)([a-zA-Z]+)$", size_str.strip())
    if not match:
        raise ValueError(f"Invalid size format: {size_str}")

    value, unit = match.groups()
    value = float(value)
    unit = unit.lower()

    # Convert to bytes
    multipliers = {
        "b": 1,
        "kb": 1024,
        "mb": 1024**2,
        "gb": 1024**3,
        "tb": 1024**4,
        "k": 1024,  # Sometimes PBS uses 'k' instead of 'kb'
        "m": 1024**2,  # Sometimes PBS uses 'm' instead of 'mb'
        "g": 1024**3,  # Sometimes PBS uses 'g' instead of 'gb'
    }

    if unit not in multipliers:
        raise ValueError(f"Unknown size unit: {unit}")

    return int(value * multipliers[unit])


def bytes_to_human(bytes_value: int) -> str:
    """
    Convert bytes to human-readable format.

    Args:
        bytes_value: Size in bytes

    Returns:
        Human-readable size string

    Examples:
        >>> bytes_to_human(4294967296)
        '4.0GB'
        >>> bytes_to_human(138833920)
        '132.4MB'
    """
    if bytes_value == 0:
        return "0B"

    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(bytes_value)

    for unit in units:
        if size < 1024.0:
            if unit == "B":
                return f"{int(size)}{unit}"
            else:
                return f"{size:.1f}{unit}"
        size /= 1024.0

    return f"{size:.1f}PB"


def time_to_seconds(time_str: str) -> int:
    """
    Convert PBS time string to seconds.

    Args:
        time_str: Time string like '00:30:00', '00:00:04', '01:45:30'

    Returns:
        Time in seconds as integer

    Examples:
        >>> time_to_seconds('00:30:00')
        1800
        >>> time_to_seconds('00:00:04')
        4
        >>> time_to_seconds('01:45:30')
        6330
    """
    if not time_str:
        return 0

    # Handle format HH:MM:SS
    match = re.match(r"^(\d+):(\d+):(\d+)$", time_str.strip())
    if not match:
        raise ValueError(f"Invalid time format: {time_str}")

    hours, minutes, seconds = map(int, match.groups())
    return hours * 3600 + minutes * 60 + seconds


def seconds_to_time(seconds: int) -> str:
    """
    Convert seconds to HH:MM:SS format.

    Args:
        seconds: Time in seconds

    Returns:
        Time string in HH:MM:SS format

    Examples:
        >>> seconds_to_time(1800)
        '00:30:00'
        >>> seconds_to_time(6330)
        '01:45:30'
    """
    if seconds < 0:
        return "00:00:00"

    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60

    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def parse_exec_host(exec_host: str) -> Dict[str, Any]:
    """
    Parse PBS exec_host string.

    Args:
        exec_host: Host string like 'gadi-cpu-clx-2841/3'

    Returns:
        Dictionary with hostname and core info

    Examples:
        >>> parse_exec_host('gadi-cpu-clx-2841/3')
        {'hostname': 'gadi-cpu-clx-2841', 'core': 3}
        >>> parse_exec_host('gadi-gpu-v100-0001/0*2')
        {'hostname': 'gadi-gpu-v100-0001', 'core': '0*2'}
    """
    if not exec_host:
        return {"hostname": None, "core": None}

    # Split on last '/' to separate hostname and core
    if "/" in exec_host:
        hostname, core_info = exec_host.rsplit("/", 1)

        # Try to convert core to int, but keep as string if it has special chars like '*'
        try:
            core = int(core_info)
        except ValueError:
            core = core_info

        return {"hostname": hostname, "core": core}
    else:
        return {"hostname": exec_host, "core": None}


def parse_exec_vnode(exec_vnode: str) -> Dict[str, Any]:
    """
    Parse PBS exec_vnode string.

    Args:
        exec_vnode: VNode string like '(gadi-cpu-clx-2841:ncpus=1:mem=4194304kb:jobfs=102400kb)'

    Returns:
        Dictionary with parsed vnode information

    Examples:
        >>> parse_exec_vnode('(gadi-cpu-clx-2841:ncpus=1:mem=4194304kb:jobfs=102400kb)')
        {'hostname': 'gadi-cpu-clx-2841', 'ncpus': 1, 'mem_kb': 4194304, 'jobfs_kb': 102400}
    """
    if not exec_vnode:
        return {}

    # Remove parentheses
    vnode_str = exec_vnode.strip("()")

    # Split by ':' to get components
    components = vnode_str.split(":")
    if not components:
        return {}

    result = {"hostname": components[0]}

    # Parse resource specifications
    for component in components[1:]:
        if "=" in component:
            key, value = component.split("=", 1)

            # Convert specific keys
            if key == "ncpus":
                result["ncpus"] = int(value)
            elif key == "mem":
                # Convert memory to KB for consistency
                result["mem_kb"] = size_to_bytes(value) // 1024
            elif key == "jobfs":
                # Convert jobfs to KB for consistency
                result["jobfs_kb"] = size_to_bytes(value) // 1024
            else:
                result[key] = value

    return result


def parse_timestamp(timestamp_str: str) -> Optional[datetime]:
    """
    Parse PBS timestamp string to datetime object.

    Args:
        timestamp_str: Timestamp like 'Sun Oct  5 10:25:10 2025'

    Returns:
        datetime object or None if parsing fails

    Examples:
        >>> parse_timestamp('Sun Oct  5 10:25:10 2025')
        datetime.datetime(2025, 10, 5, 10, 25, 10)
    """
    if not timestamp_str:
        return None

    try:
        # Handle PBS format: 'Sun Oct  5 10:25:10 2025'
        # Note the double space before single-digit day
        return datetime.strptime(timestamp_str, "%a %b %d %H:%M:%S %Y")
    except ValueError:
        try:
            # Try alternative format without day name
            return datetime.strptime(timestamp_str, "%b %d %H:%M:%S %Y")
        except ValueError:
            return None


def calculate_efficiency(used: int, requested: int) -> float:
    """
    Calculate efficiency percentage.

    Args:
        used: Amount used
        requested: Amount requested

    Returns:
        Efficiency percentage (0-100)

    Examples:
        >>> calculate_efficiency(135580, 4194304)  # KB
        3.2
        >>> calculate_efficiency(4, 1800)  # seconds
        0.2
    """
    if requested == 0:
        return 0.0

    efficiency = (used / requested) * 100
    return round(efficiency, 1)


def parse_joblog_resources(joblog_path: str) -> Optional[Dict[str, Any]]:
    """
    Parse resource usage information from a PBS joblog file.

    The joblog contains a resource usage section at the end that looks like:

    ========================================================================
    Resource Usage on 2025-10-21 11:14:48:
       Job Id:             152942530.gadi-pbs
       Project:            a56
       Exit Status:        0
       Service Units:      0.00
       NCPUs Requested:    1                      NCPUs Used: 1
       CPU Time Used: 00:00:02
       Memory Requested:   4.0GB                 Memory Used: 131.66MB
       Walltime requested: 02:00:00            Walltime Used: 00:00:05
       JobFS requested:    100.0MB                JobFS used: 0B
    ========================================================================

    Args:
        joblog_path: Path to the PBS joblog file

    Returns:
        Dictionary with parsed resource information, or None if parsing fails

    Example:
        >>> parse_joblog_resources("/path/to/joblog")
        {
            'job_id': '152942530.gadi-pbs',
            'project': 'a56',
            'exit_status': 0,
            'service_units': 0.0,
            'ncpus_requested': 1,
            'ncpus_used': 1,
            'cpu_time_used_seconds': 2,
            'memory_requested_bytes': 4294967296,
            'memory_used_bytes': 138131251,
            'walltime_requested_seconds': 7200,
            'walltime_used_seconds': 5,
            'jobfs_requested_bytes': 104857600,
            'jobfs_used_bytes': 0
        }
    """
    try:
        with open(joblog_path, "r") as f:
            content = f.read()

        # Look for the Resource Usage section - it's at the end of the file
        # Pattern to match the resource usage block - more flexible regex
        resource_section_pattern = r"Resource Usage on.*?:\s*\n(.*?)={10,}"
        match = re.search(resource_section_pattern, content, re.DOTALL)

        if not match:
            return None

        resource_text = match.group(1)

        # Parse individual fields using regex patterns
        result = {}

        # Job ID
        job_id_match = re.search(r"Job Id:\s*(\S+)", resource_text)
        if job_id_match:
            result["job_id"] = job_id_match.group(1)

        # Project
        project_match = re.search(r"Project:\s*(\S+)", resource_text)
        if project_match:
            result["project"] = project_match.group(1)

        # Exit Status
        exit_status_match = re.search(r"Exit Status:\s*(\d+)", resource_text)
        if exit_status_match:
            result["exit_status"] = int(exit_status_match.group(1))

        # Service Units
        service_units_match = re.search(r"Service Units:\s*([\d.]+)", resource_text)
        if service_units_match:
            result["service_units"] = float(service_units_match.group(1))

        # NCPUs Requested and Used
        ncpus_match = re.search(
            r"NCPUs Requested:\s*(\d+)\s+NCPUs Used:\s*(\d+)", resource_text
        )
        if ncpus_match:
            result["ncpus_requested"] = int(ncpus_match.group(1))
            result["ncpus_used"] = int(ncpus_match.group(2))

        # CPU Time Used
        cpu_time_match = re.search(r"CPU Time Used:\s*(\d+:\d+:\d+)", resource_text)
        if cpu_time_match:
            result["cpu_time_used_seconds"] = time_to_seconds(cpu_time_match.group(1))

        # Memory Requested and Used
        memory_match = re.search(
            r"Memory Requested:\s*([\d.]+\w+)\s+Memory Used:\s*([\d.]+\w+)",
            resource_text,
        )
        if memory_match:
            result["memory_requested_bytes"] = size_to_bytes(memory_match.group(1))
            result["memory_used_bytes"] = size_to_bytes(memory_match.group(2))

        # Walltime Requested and Used
        walltime_match = re.search(
            r"Walltime requested:\s*(\d+:\d+:\d+)\s+Walltime Used:\s*(\d+:\d+:\d+)",
            resource_text,
        )
        if walltime_match:
            result["walltime_requested_seconds"] = time_to_seconds(
                walltime_match.group(1)
            )
            result["walltime_used_seconds"] = time_to_seconds(walltime_match.group(2))

        # JobFS Requested and Used
        jobfs_match = re.search(
            r"JobFS requested:\s*([\d.]+\w+)\s+JobFS used:\s*([\d.]+\w+)", resource_text
        )
        if jobfs_match:
            result["jobfs_requested_bytes"] = size_to_bytes(jobfs_match.group(1))
            result["jobfs_used_bytes"] = size_to_bytes(jobfs_match.group(2))

        return result if result else None

    except Exception:
        return None


if __name__ == "__main__":
    # Test the functions
    print("Testing size_to_bytes:")
    print(f"4294967296b -> {size_to_bytes('4294967296b')} bytes")
    print(f"135580kb -> {size_to_bytes('135580kb')} bytes")
    print(f"4gb -> {size_to_bytes('4gb')} bytes")

    print("\nTesting bytes_to_human:")
    print(f"4294967296 -> {bytes_to_human(4294967296)}")
    print(f"138833920 -> {bytes_to_human(138833920)}")

    print("\nTesting time_to_seconds:")
    print(f"00:30:00 -> {time_to_seconds('00:30:00')} seconds")
    print(f"00:00:04 -> {time_to_seconds('00:00:04')} seconds")

    print("\nTesting parse_exec_host:")
    print(f"gadi-cpu-clx-2841/3 -> {parse_exec_host('gadi-cpu-clx-2841/3')}")

    print("\nTesting parse_exec_vnode:")
    vnode = "(gadi-cpu-clx-2841:ncpus=1:mem=4194304kb:jobfs=102400kb)"
    print(f"{vnode} -> {parse_exec_vnode(vnode)}")

    print("\nTesting calculate_efficiency:")
    print(f"Memory: {calculate_efficiency(135580, 4194304)}%")
    print(f"Time: {calculate_efficiency(4, 1800)}%")
