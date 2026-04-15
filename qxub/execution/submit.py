"""
Programmatic job submission API for qxub.

Provides :func:`submit_job` — a pure-Python entrypoint that submits a PBS
job through qxub's full infrastructure (templates, status files, DB entries,
resource tracking) **without** depending on Click or the CLI layer.

Usage::

    from qxub.execution import submit_job

    result = submit_job(
        command=("python", "train.py"),
        resources={"walltime": "02:00:00", "mem": "8GB", "ncpus": 4},
        name="my_job",
        queue="normal",
        project="a56",
        env="myenv",          # conda environment
    )
    print(result.job_id)       # "12345.gadi-pbs"
    print(result.virtual_id)   # "qx-550e8400-..."
"""

import base64
import logging

logger = logging.getLogger(__name__)
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union


@dataclass
class SubmitResult:
    """Structured result from a programmatic job submission.

    Attributes:
        job_id:       Real PBS job ID (e.g. ``"12345.gadi-pbs"``).
        virtual_id:   qxub virtual ID (e.g. ``"qx-550e8400-..."``).
                      ``None`` if the queue entry could not be created.
        joblog_path:  Absolute path to the PBS joblog file.
        stdout_path:  Absolute path to the PBS stdout file.
        stderr_path:  Absolute path to the PBS stderr file.
        submission_command:
                      The raw ``qsub -v ...`` command that was executed
                      (useful for debugging).
    """

    job_id: str
    virtual_id: Optional[str] = None
    joblog_path: Optional[str] = None
    stdout_path: Optional[str] = None
    stderr_path: Optional[str] = None
    submission_command: str = ""


def submit_job(
    command: Union[Tuple[str, ...], List[str], str],
    *,
    # --- Execution context (mutually exclusive) ---
    env: Optional[str] = None,
    modules: Optional[List[str]] = None,
    sif: Optional[str] = None,
    bind: Optional[str] = None,
    # --- PBS resources ---
    resources: Optional[Dict[str, Any]] = None,
    queue: Optional[str] = None,
    name: Optional[str] = None,
    project: Optional[str] = None,
    # --- Paths ---
    log_dir: Optional[str] = None,
    execdir: Optional[str] = None,
    # --- Hooks ---
    pre: Optional[str] = None,
    post: Optional[str] = None,
    # --- Metadata ---
    tags: Optional[Sequence[str]] = None,
    # --- Options ---
    template: Optional[str] = None,
    config_path: Optional[str] = None,
) -> SubmitResult:
    """Submit a PBS job through qxub's full infrastructure.

    This is the programmatic equivalent of ``qxub exec --terse``.  It
    resolves configuration, selects a template, creates a queue entry,
    calls ``qsub``, updates the queue entry with the real PBS job ID,
    and logs to the resource tracker and execution history.

    Unlike the CLI path, this function:

    * **Never calls** ``sys.exit()`` — failures raise exceptions.
    * **Never monitors** the job — the caller handles status checking.
    * **Returns structured data** via :class:`SubmitResult`.

    Args:
        command:      The command to execute.  Can be a tuple/list of strings
                      or a single string.
        env:          Conda environment name (mutually exclusive with
                      *modules* and *sif*).
        modules:      List of environment module names.
        sif:          Path to a Singularity container (``.sif`` file).
        bind:         Singularity bind mounts (only used with *sif*).
        resources:    Resource requirements as a dict, e.g.
                      ``{"walltime": "02:00:00", "mem": "8GB", "ncpus": 4}``.
                      Keys are the same as PBS ``-l`` keys.
        queue:        PBS queue name (default: from config or ``"normal"``).
        name:         PBS job name (default: from config or ``"qt"``).
        project:      PBS project code (default: from config or ``$PROJECT``).
        log_dir:      Override log directory for .out / .err / .log files.
        execdir:      Working directory on the compute node (default: CWD).
        pre:          Shell command to run before the main command.
        post:         Shell command to run after the main command.
        tags:         List of tag strings for job metadata.
        template:     Custom job script template path.
        config_path:  Path to a qxub config file (overrides the default
                      config hierarchy).

    Returns:
        A :class:`SubmitResult` containing the PBS job ID, virtual ID,
        and file paths.

    Raises:
        QsubError:    If the ``qsub`` command fails.
        ValueError:   If multiple execution contexts are specified.
    """
    from ..config.handler import process_job_options
    from ..core.scheduler import QsubError, qsub
    from ..core.templates import get_template
    from ..queue import create_queue_entry, update_queue_entry
    from ..queue.db import get_db_path
    from ..resources.tracker import resource_tracker

    # ------------------------------------------------------------------
    # Validate execution context (mutually exclusive)
    # ------------------------------------------------------------------
    contexts = [env, modules, sif]
    if sum(bool(c) for c in contexts) > 1:
        raise ValueError(
            "Only one execution context may be specified: env, modules, or sif"
        )

    # ------------------------------------------------------------------
    # Normalise command
    # ------------------------------------------------------------------
    if isinstance(command, str):
        cmd_parts: Tuple[str, ...] = tuple(command.split())
    else:
        cmd_parts = tuple(command)

    cmd_str = " ".join(cmd_parts)

    # ------------------------------------------------------------------
    # Determine execution context + template
    # ------------------------------------------------------------------
    if env:
        context_type = "conda"
        context_value: Any = env
        template_type = "conda"
    elif modules:
        context_type = "module"
        context_value = modules
        template_type = "module"
    elif sif:
        context_type = "singularity"
        context_value = sif
        template_type = "singularity"
    else:
        context_type = "default"
        context_value = None
        template_type = "default"

    if not template:
        template = get_template(template_type)

    # ------------------------------------------------------------------
    # Build PBS resource list from dict
    # ------------------------------------------------------------------
    pbs_resources: List[str] = []
    if resources:
        for key, value in resources.items():
            pbs_resources.append(f"{key}={value}")

    # ------------------------------------------------------------------
    # Config override
    # ------------------------------------------------------------------
    if config_path:
        from ..config.manager import set_config_override

        set_config_override(config_path)

    from ..config import manager as config_mod

    config_manager = config_mod.config_manager

    # ------------------------------------------------------------------
    # Build params dict (mimics what exec_cli constructs)
    # ------------------------------------------------------------------
    params: Dict[str, Any] = {
        "command": cmd_parts,
        "resources": tuple(pbs_resources),
        "queue": queue,
        "name": name,
        "project": project,
        "log_dir": log_dir,
        "out": None,
        "err": None,
        "joblog": None,
        "execdir": execdir,
        "create_execdir": False,
        "email": None,
        "email_opts": None,
        "array": None,
        "dry": False,
        "quiet": True,
        "terse": True,
        "verbose": 0,
        "cpus_explicit": "ncpus" in (resources or {}),
        "internet": False,
        "tags": list(tags or []),
    }

    # ------------------------------------------------------------------
    # Process config defaults + template vars + qsub options string
    # ------------------------------------------------------------------
    processed_params, qsub_options = process_job_options(params, config_manager)

    # ------------------------------------------------------------------
    # Build submission variables (same logic as ExecutionContext)
    # ------------------------------------------------------------------
    cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")
    exec_dir = processed_params.get("execdir") or os.getcwd()

    base_vars = (
        f'cmd_b64="{cmd_b64}"'
        f",cwd={exec_dir}"
        f",out={processed_params['out']}"
        f",err={processed_params['err']}"
        f",quiet=true"
        f",create_execdir=false"
    )

    # Context-specific variables
    if context_type == "conda":
        context_vars = f"env={context_value}"
    elif context_type == "module":
        mods_str = " ".join(context_value)
        context_vars = f'mods="{mods_str}"'
    elif context_type == "singularity":
        context_vars = f"sif={context_value}"
        if bind:
            context_vars += f',bind="{bind}"'
    else:
        context_vars = ""

    if context_vars:
        submission_vars = f"{context_vars},{base_vars}"
    else:
        submission_vars = base_vars

    # Pre/post hooks
    if pre:
        pre_b64 = base64.b64encode(pre.encode("utf-8")).decode("ascii")
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_b64 = base64.b64encode(post.encode("utf-8")).decode("ascii")
        submission_vars += f',post_cmd_b64="{post_b64}"'

    # Shared DB path for job scripts
    db_path = str(get_db_path())
    submission_vars += f",QXUB_SHARED_DB={db_path}"

    # ------------------------------------------------------------------
    # Build qsub command
    # ------------------------------------------------------------------
    submission_command = f"qsub -v {submission_vars} {qsub_options} {template}"
    logger.info("Programmatic submission command: %s", submission_command)

    # ------------------------------------------------------------------
    # Pre-qsub: create queue entry (status = initiated)
    # ------------------------------------------------------------------
    exec_context_info = {
        "type": context_type,
        "value": (
            context_value
            if not isinstance(context_value, list)
            else " ".join(context_value)
        ),
    }

    tag_list = list(tags or [])

    virtual_id = None
    try:
        ncpus = int(resources.get("ncpus", 0)) if resources else None
        virtual_id = create_queue_entry(
            command=cmd_str,
            pbs_job_id=None,
            tags=tag_list,
            working_dir=exec_dir,
            exec_context=exec_context_info,
            joblog_path=processed_params.get("joblog"),
            queue_name=processed_params.get("queue"),
            cpus_requested=ncpus if ncpus else None,
        )
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("Failed to create queue entry: %s", exc)

    # ------------------------------------------------------------------
    # Submit to PBS (raises QsubError on failure — no sys.exit)
    # ------------------------------------------------------------------
    try:
        job_id = qsub(submission_command, quiet=True)
    except QsubError:
        if virtual_id:
            try:
                update_queue_entry(virtual_id, status="failed")
            except Exception:  # pylint: disable=broad-except
                pass
        raise  # Let QsubError propagate to the caller

    # ------------------------------------------------------------------
    # Post-qsub: update queue entry (status → dispatched)
    # ------------------------------------------------------------------
    if virtual_id:
        try:
            update_queue_entry(
                virtual_id,
                pbs_job_id=job_id,
                status="dispatched",
                joblog_path=processed_params.get("joblog"),
            )
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("Failed to update queue entry: %s", exc)

    # ------------------------------------------------------------------
    # Resource tracker + history logging
    # ------------------------------------------------------------------
    try:
        ncpus_req = int(resources.get("ncpus", 0)) if resources else None
        resource_tracker.log_job_submitted(
            job_id=job_id,
            command=cmd_str,
            tags=tag_list,
            joblog_path=processed_params.get("joblog"),
            queue=processed_params.get("queue"),
            cpus_requested=ncpus_req if ncpus_req else None,
        )
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("Failed to log job submission: %s", exc)

    # ------------------------------------------------------------------
    # Return structured result
    # ------------------------------------------------------------------
    return SubmitResult(
        job_id=job_id,
        virtual_id=virtual_id,
        joblog_path=processed_params.get("joblog"),
        stdout_path=str(processed_params["out"]),
        stderr_path=str(processed_params["err"]),
        submission_command=submission_command,
    )
