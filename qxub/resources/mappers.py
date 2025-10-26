"""
Resource mapping utilities for workflow engine compatibility.

This module provides functions to convert workflow-specific resource
specifications into PBS Pro format, enabling better integration with
Snakemake, NextFlow, CWL, and other workflow engines.
"""

import re
from typing import Dict, List, Tuple, Union


class ResourceMapper:
    """Convert workflow-friendly resource specs to PBS Pro format."""

    def __init__(self):
        self.pbs_resources = []

    def add_memory(self, memory: Union[str, int], unit: str = "MB") -> None:
        """Add memory specification.

        Args:
            memory: Memory amount (e.g., "4GB", "1000MB", 1000)
            unit: Default unit if memory is integer ("MB", "GB")
        """
        if isinstance(memory, int):
            self.pbs_resources.append(f"mem={memory}{unit}")
        elif isinstance(memory, str):
            # Parse memory strings like "4GB", "1000MB", "4 GB"
            memory_clean = memory.replace(" ", "").upper()
            if re.match(r"^\d+[KMGT]?B?$", memory_clean):
                self.pbs_resources.append(f"mem={memory_clean}")
            else:
                raise ValueError(f"Invalid memory format: {memory}")

    def add_runtime(self, runtime: Union[str, int], unit: str = "minutes") -> None:
        """Add runtime/walltime specification.

        Args:
            runtime: Runtime (e.g., "2h", "120m", "1:30:00", 120)
            unit: Default unit if runtime is integer ("minutes", "hours", "seconds")
        """
        if isinstance(runtime, int):
            if unit == "minutes":
                hours = runtime // 60
                mins = runtime % 60
                self.pbs_resources.append(f"walltime={hours}:{mins:02d}:00")
            elif unit == "hours":
                self.pbs_resources.append(f"walltime={runtime}:00:00")
            elif unit == "seconds":
                hours = runtime // 3600
                mins = (runtime % 3600) // 60
                secs = runtime % 60
                self.pbs_resources.append(f"walltime={hours}:{mins:02d}:{secs:02d}")
        elif isinstance(runtime, str):
            runtime_clean = runtime.replace(" ", "").lower()

            # Handle HH:MM:SS format
            if re.match(r"^\d{1,2}:\d{2}:\d{2}$", runtime_clean):
                self.pbs_resources.append(f"walltime={runtime_clean}")
            # Handle HH:MM format
            elif re.match(r"^\d{1,2}:\d{2}$", runtime_clean):
                self.pbs_resources.append(f"walltime={runtime_clean}:00")
            else:
                # Parse complex formats like "1h30m", "2h", "90m", "3600s"
                total_seconds = 0

                # Extract hours
                hours_match = re.search(r"(\d+)h", runtime_clean)
                if hours_match:
                    total_seconds += int(hours_match.group(1)) * 3600

                # Extract minutes
                minutes_match = re.search(r"(\d+)m", runtime_clean)
                if minutes_match:
                    total_seconds += int(minutes_match.group(1)) * 60

                # Extract seconds
                seconds_match = re.search(r"(\d+)s", runtime_clean)
                if seconds_match:
                    total_seconds += int(seconds_match.group(1))

                # Handle simple numeric strings (assume minutes if no unit)
                if not any([hours_match, minutes_match, seconds_match]):
                    if runtime_clean.isdigit():
                        total_seconds = int(runtime_clean) * 60  # Assume minutes
                    else:
                        raise ValueError(f"Invalid runtime format: {runtime}")

                if total_seconds == 0:
                    raise ValueError(f"Invalid runtime format: {runtime}")

                # Convert to HH:MM:SS
                hours = total_seconds // 3600
                minutes = (total_seconds % 3600) // 60
                seconds = total_seconds % 60
                self.pbs_resources.append(
                    f"walltime={hours}:{minutes:02d}:{seconds:02d}"
                )

    def add_cpus(self, cpus: int) -> None:
        """Add CPU/cores specification."""
        self.pbs_resources.append(f"ncpus={cpus}")

    def add_disk(self, disk: Union[str, int], unit: str = "MB") -> None:
        """Add disk/jobfs specification."""
        if isinstance(disk, int):
            self.pbs_resources.append(f"jobfs={disk}{unit}")
        elif isinstance(disk, str):
            disk_clean = disk.replace(" ", "").upper()
            if re.match(r"^\d+[KMGT]?B?$", disk_clean):
                self.pbs_resources.append(f"jobfs={disk_clean}")
            else:
                raise ValueError(f"Invalid disk format: {disk}")

    def add_custom(self, key: str, value: str) -> None:
        """Add custom PBS resource specification."""
        self.pbs_resources.append(f"{key}={value}")

    def get_pbs_resources(self) -> List[str]:
        """Get the PBS resource specifications."""
        return self.pbs_resources

    def get_pbs_string(self) -> str:
        """Get PBS resources as a single comma-separated string."""
        return ",".join(self.pbs_resources)


def map_snakemake_resources(resources: Dict) -> List[str]:
    """Convert Snakemake resources dict to PBS format.

    Args:
        resources: Snakemake resources dict with keys like:
                  mem_mb, runtime, threads, disk_mb, tmpdir, etc.

    Returns:
        List of PBS resource specifications

    Example:
        >>> resources = {'mem_mb': 4000, 'runtime': 120, 'threads': 4}
        >>> map_snakemake_resources(resources)
        ['mem=4000MB', 'walltime=2:00:00', 'ncpus=4']
    """
    mapper = ResourceMapper()

    # Memory
    if "mem_mb" in resources:
        mapper.add_memory(resources["mem_mb"], "MB")
    elif "mem_gb" in resources:
        mapper.add_memory(resources["mem_gb"], "GB")
    elif "mem" in resources:
        mapper.add_memory(resources["mem"])

    # Runtime
    if "runtime" in resources:
        mapper.add_runtime(resources["runtime"], "minutes")
    elif "walltime" in resources:
        mapper.add_runtime(resources["walltime"])

    # CPUs/Threads
    if "threads" in resources:
        mapper.add_cpus(resources["threads"])
    elif "cpus" in resources:
        mapper.add_cpus(resources["cpus"])

    # Disk
    if "disk_mb" in resources:
        mapper.add_disk(resources["disk_mb"], "MB")
    elif "disk_gb" in resources:
        mapper.add_disk(resources["disk_gb"], "GB")
    elif "jobfs" in resources:
        mapper.add_disk(resources["jobfs"])

    # Custom resources (passed through if they look like PBS resources)
    for key, value in resources.items():
        if key not in [
            "mem_mb",
            "mem_gb",
            "mem",
            "runtime",
            "walltime",
            "threads",
            "cpus",
            "disk_mb",
            "disk_gb",
            "jobfs",
        ]:
            # Check if it's a valid PBS resource name
            if isinstance(value, (str, int)):
                mapper.add_custom(key, str(value))

    return mapper.get_pbs_resources()


def map_nextflow_resources(
    memory: str = None, time: str = None, cpus: int = None, disk: str = None
) -> List[str]:
    """Convert NextFlow-style resources to PBS format.

    Args:
        memory: Memory specification (e.g., "8 GB", "4000 MB")
        time: Time specification (e.g., "2h", "120m", "1:30:00")
        cpus: Number of CPUs
        disk: Disk specification (e.g., "100 GB")

    Returns:
        List of PBS resource specifications
    """
    mapper = ResourceMapper()

    if memory:
        mapper.add_memory(memory)
    if time:
        mapper.add_runtime(time)
    if cpus:
        mapper.add_cpus(cpus)
    if disk:
        mapper.add_disk(disk)

    return mapper.get_pbs_resources()


def map_cwl_resources(
    cores_min: int = None,
    ram_min: int = None,
    tmpdir_min: int = None,
    time_limit: int = None,
) -> List[str]:
    """Convert CWL ResourceRequirement to PBS format.

    Args:
        cores_min: Minimum cores required
        ram_min: Minimum RAM in MB
        tmpdir_min: Minimum temp directory space in MB
        time_limit: Time limit in seconds

    Returns:
        List of PBS resource specifications
    """
    mapper = ResourceMapper()

    if cores_min:
        mapper.add_cpus(cores_min)
    if ram_min:
        mapper.add_memory(ram_min, "MB")
    if tmpdir_min:
        mapper.add_disk(tmpdir_min, "MB")
    if time_limit:
        mapper.add_runtime(time_limit, "seconds")

    return mapper.get_pbs_resources()


# Convenience functions for common patterns
def parse_resource_string(resource_str: str) -> Tuple[str, str]:
    """Parse resource strings like 'mem=4GB' or 'walltime=2:00:00'."""
    if "=" in resource_str:
        key, value = resource_str.split("=", 1)
        return key.strip(), value.strip()
    else:
        raise ValueError(f"Invalid resource format: {resource_str}")


def standardize_memory_unit(memory_str: str) -> str:
    """Standardize memory units to PBS format."""
    # Convert common variations to PBS format
    memory_str = memory_str.upper().replace(" ", "")

    # Handle variations like 'G', 'GB', 'GiB'
    memory_str = re.sub(r"GIB?$", "GB", memory_str)
    memory_str = re.sub(r"MIB?$", "MB", memory_str)
    memory_str = re.sub(r"KIB?$", "KB", memory_str)

    return memory_str
