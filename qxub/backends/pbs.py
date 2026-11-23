"""PBS Pro (and Torque) backend: ``qstat`` and ``qdel``.

All of qxub's ``qstat``/``qdel`` calls live here. The parsers are plain
functions so they can be tested on captured output.
"""

import getpass
import logging
import os
import re
import subprocess
from datetime import datetime
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple

from ..resources import (
    bytes_to_human,
    calculate_efficiency,
    parse_exec_host,
    parse_exec_vnode,
    parse_timestamp,
    size_to_bytes,
    time_to_seconds,
)
from .base import LiveJob, SchedulerBackend, SchedulerQueryError

logger = logging.getLogger(__name__)

# Job IDs per ``qstat -fx`` call. PBS Pro takes many on one command line.
QSTAT_CHUNK = 100

# A PBS job ID: "12345", "12345.server", "12345[].server" (array parent) or
# "12345[7].server" (subjob). Header, separator and server-name lines of a
# ``qstat -u`` listing never start with a digit.
JOB_ID_RE = re.compile(r"^\d+(\[\d*\])?(\.\S+)?$")


def _current_user() -> str:
    return os.environ.get("USER") or getpass.getuser()


class PbsBackend(SchedulerBackend):
    """PBS Pro. Torque's ``qstat``/``qdel`` output has the same shape."""

    name = "pbs"

    def list_jobs(
        self, user: Optional[str] = None, full_names: bool = False
    ) -> List[LiveJob]:
        user = user or _current_user()
        # -w: don't truncate job IDs ("153392916.gadi-*"), which qdel would
        # reject. -f prints whole names, but for every job you can see.
        argv = ["qstat", "-f"] if full_names else ["qstat", "-w", "-u", user]
        try:
            result = subprocess.run(argv, capture_output=True, text=True, check=False)
        except OSError as exc:
            raise SchedulerQueryError(f"cannot run qstat: {exc}") from exc
        if result.returncode != 0:
            raise SchedulerQueryError(
                f"{' '.join(argv)} failed (exit {result.returncode}): "
                f"{result.stderr.strip()}"
            )
        if full_names:
            return parse_qstat_full_listing(result.stdout, user)
        return parse_qstat_user_listing(result.stdout)

    def job_details(self, job_ids: Iterable[str]) -> Dict[str, Dict[str, Any]]:
        ids = list(dict.fromkeys(job_ids))
        found: Dict[str, Dict[str, Any]] = {}
        for start in range(0, len(ids), QSTAT_CHUNK):
            chunk = ids[start : start + QSTAT_CHUNK]
            try:
                # Unknown IDs make qstat exit non-zero, but it still prints
                # the jobs it knows.
                out = subprocess.run(
                    ["qstat", "-fx", *chunk],
                    capture_output=True,
                    text=True,
                    timeout=120,
                    check=False,
                ).stdout
            except (OSError, subprocess.SubprocessError) as exc:
                logger.debug("qstat -fx failed: %s", exc)
                continue
            for job_id, text in qstat_blocks(out):
                data = parse_qstat_fx_output(text)
                if data:
                    found[job_id] = data
        return found

    def is_finished(self, details: Optional[Dict[str, Any]]) -> bool:
        return ((details or {}).get("metadata") or {}).get("job_state") == "F"

    def cancel(self, job_id: str) -> Tuple[bool, str]:
        try:
            result = subprocess.run(
                ["qdel", job_id], capture_output=True, text=True, check=False
            )
        except OSError as exc:
            return False, f"cannot run qdel: {exc}"
        if result.returncode == 0:
            return True, ""
        return False, result.stderr.strip() or f"qdel exited {result.returncode}"

    def looks_like_job_id(self, value: str) -> bool:
        return bool(JOB_ID_RE.match(value or ""))

    def same_job(self, job_id: str) -> str:
        """``12345`` and ``12345.gadi-pbs`` are the same job."""
        return job_id.split(".", 1)[0]


# -- parsers -------------------------------------------------------------------


def parse_qstat_user_listing(output: str) -> List[LiveJob]:
    """Parse ``qstat -w -u USER`` output.

    The listing looks like (PBS Pro; Torque is the same shape)::

        server:
                                                     Req'd  Req'd   Elap
        Job ID          Username Queue Jobname SessID NDS TSK Memory Time  S Time
        --------------- -------- ----- ------- ------ --- --- ------ ----- - -----
        12345.server    alice    normal myjob   4321   1   4  16gb 01:00 R 00:10

    Rows are recognised by a job-ID-shaped first token rather than by
    skipping a fixed number of header lines. The state is the second-to-last
    column (just before elapsed time), the queue the third and the name the
    fourth. PBS truncates long names in this listing.
    """
    jobs = []
    for line in output.splitlines():
        parts = line.split()
        if len(parts) < 4 or not JOB_ID_RE.match(parts[0]):
            continue
        jobs.append(LiveJob(parts[0], parts[2], parts[-2], parts[3]))
    return jobs


def qstat_blocks(output: str) -> Iterator[Tuple[str, str]]:
    """Split ``qstat -f``/``-fx`` output into ``(job_id, text)`` per job."""
    job_id, lines = None, []
    for line in output.splitlines():
        if line.startswith("Job Id:"):
            if job_id:
                yield job_id, "\n".join(lines)
            job_id, lines = line.split(":", 1)[1].strip(), []
        elif job_id:
            lines.append(line)
    if job_id:
        yield job_id, "\n".join(lines)


def _attributes(text: str) -> Dict[str, str]:
    attrs = {}
    for line in text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def parse_qstat_full_listing(output: str, user: Optional[str] = None) -> List[LiveJob]:
    """Parse ``qstat -f`` output (no job IDs) into :class:`LiveJob` rows.

    Names are whole. If ``user`` is given, jobs whose ``Job_Owner`` is
    someone else are left out; jobs without a ``Job_Owner`` are kept.
    """
    jobs = []
    for job_id, text in qstat_blocks(output):
        attrs = _attributes(text)
        owner = attrs.get("Job_Owner", "").split("@", 1)[0]
        if user and owner and owner != user:
            continue
        jobs.append(
            LiveJob(
                job_id,
                attrs.get("queue", ""),
                attrs.get("job_state", ""),
                attrs.get("Job_Name", ""),
            )
        )
    return jobs


def parse_qstat_fx_output(qstat_output):
    """
    Parse qstat -fx output to extract comprehensive job resource information.

    Args:
        qstat_output: Raw output from qstat -fx command

    Returns:
        dict: Parsed job resource data
    """
    if not qstat_output:
        return None

    data = {
        "exit_status": None,
        "resources_requested": {},
        "resources_used": {},
        "execution": {},
        "timing": {},
        "metadata": {},
    }

    # Parse line by line
    for line in qstat_output.split("\n"):
        line = line.strip()
        if not line or "=" not in line:
            continue

        # Split on first '=' to handle values with '=' in them
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()

        # Parse different categories of information
        if key == "Exit_status":
            try:
                data["exit_status"] = int(value)
            except ValueError:
                logger.warning("Could not parse exit status: %s", value)

        elif key.startswith("Resource_List."):
            # Requested resources
            resource_key = key[14:]  # Remove 'Resource_List.' prefix
            data["resources_requested"][resource_key] = value

        elif key.startswith("resources_used."):
            # Used resources
            resource_key = key[15:]  # Remove 'resources_used.' prefix
            data["resources_used"][resource_key] = value

        elif key in [
            "exec_host",
            "exec_vnode",
            "queue",
            "session_id",
            "project",
            "group_list",
        ]:
            # Execution environment
            data["execution"][key] = value

        elif key in ["qtime", "stime", "etime", "mtime", "ctime", "obittime"]:
            # Timing information
            data["timing"][key] = value

        elif key in ["Job_Id", "Job_Name", "Job_Owner", "job_state"]:
            # Job metadata
            data["metadata"][key] = value

    # Post-process resource data for easier consumption
    data = _enhance_resource_data(data)

    return data


def _enhance_resource_data(data):
    """
    Enhance parsed resource data with calculated fields and conversions.

    Args:
        data: Raw parsed data from parse_qstat_fx_output

    Returns:
        dict: Enhanced data with calculated efficiency metrics
    """
    if not data:
        return data

    # Convert resource values to standardized formats
    requested = data["resources_requested"]
    used = data["resources_used"]

    # Parse timing information
    timing = data["timing"]
    for time_key in list(timing.keys()):  # Convert to list to avoid iteration issues
        if timing[time_key]:
            parsed_time = parse_timestamp(timing[time_key])
            if parsed_time:
                timing[f"{time_key}_parsed"] = parsed_time.isoformat()

    # Calculate derived timing metrics
    if "qtime_parsed" in timing and "stime_parsed" in timing:
        try:
            qtime = datetime.fromisoformat(timing["qtime_parsed"])
            stime = datetime.fromisoformat(timing["stime_parsed"])
            timing["queue_wait_seconds"] = int((stime - qtime).total_seconds())
        except Exception:
            pass

    if "stime_parsed" in timing and "mtime_parsed" in timing:
        try:
            stime = datetime.fromisoformat(timing["stime_parsed"])
            mtime = datetime.fromisoformat(timing["mtime_parsed"])
            timing["execution_seconds"] = int((mtime - stime).total_seconds())
        except Exception:
            pass

    # Parse execution environment
    execution = data["execution"]
    if "exec_host" in execution:
        execution["exec_host_parsed"] = parse_exec_host(execution["exec_host"])
    if "exec_vnode" in execution:
        execution["exec_vnode_parsed"] = parse_exec_vnode(execution["exec_vnode"])

    # Calculate resource efficiency metrics
    efficiency = {}

    # Memory efficiency
    if "mem" in requested and "mem" in used:
        try:
            mem_requested_bytes = size_to_bytes(requested["mem"])
            mem_used_bytes = size_to_bytes(used["mem"])
            efficiency["memory_efficiency"] = calculate_efficiency(
                mem_used_bytes, mem_requested_bytes
            )
            efficiency["memory_requested_human"] = bytes_to_human(mem_requested_bytes)
            efficiency["memory_used_human"] = bytes_to_human(mem_used_bytes)
        except Exception as e:
            logger.warning("Could not calculate memory efficiency: %s", e)

    # Time efficiency
    if "walltime" in requested and "walltime" in used:
        try:
            time_requested_seconds = time_to_seconds(requested["walltime"])
            time_used_seconds = time_to_seconds(used["walltime"])
            efficiency["time_efficiency"] = calculate_efficiency(
                time_used_seconds, time_requested_seconds
            )
        except Exception as e:
            logger.warning("Could not calculate time efficiency: %s", e)

    # CPU efficiency (from cpupercent if available)
    if "cpupercent" in used:
        try:
            efficiency["cpu_efficiency"] = float(used["cpupercent"])
        except Exception:
            pass

    # Jobfs efficiency
    if "jobfs" in requested and "jobfs" in used:
        try:
            jobfs_requested_bytes = size_to_bytes(requested["jobfs"])
            jobfs_used_bytes = size_to_bytes(used["jobfs"])
            efficiency["jobfs_efficiency"] = calculate_efficiency(
                jobfs_used_bytes, jobfs_requested_bytes
            )
            efficiency["jobfs_requested_human"] = bytes_to_human(jobfs_requested_bytes)
            efficiency["jobfs_used_human"] = bytes_to_human(jobfs_used_bytes)
        except Exception as e:
            logger.warning("Could not calculate jobfs efficiency: %s", e)

    data["efficiency"] = efficiency

    return data
