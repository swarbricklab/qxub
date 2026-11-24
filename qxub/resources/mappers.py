"""
Resource mapping utilities for workflow engine compatibility.

This module provides functions to convert workflow-specific resource
specifications into PBS Pro format, enabling better integration with
Snakemake, NextFlow, CWL, and other workflow engines.
"""

import re
from typing import List, Union


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

    def add_gpus(self, gpus: int) -> None:
        """Add GPU count specification."""
        self.pbs_resources.append(f"ngpus={gpus}")

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

    def add_storage(self, storage: str) -> None:
        """Add storage volume specification for NCI storage mounts.

        Args:
            storage: Storage volumes to mount (e.g., 'gdata/ab12', 'gdata/ab12+gdata/cd34')
        """
        if isinstance(storage, str):
            # For NCI storage volumes, pass through directly to PBS
            # Examples: 'gdata/ab12', 'gdata/ab12+gdata/cd34', 'scratch/ab12+gdata/ab12'
            self.pbs_resources.append(f"storage={storage}")
        else:
            raise ValueError(f"Storage must be a string specifying volumes: {storage}")

    def get_pbs_resources(self) -> List[str]:
        """Get the PBS resource specifications."""
        return self.pbs_resources
