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
        project="ab12",
        env="myenv",          # conda environment
    )
    print(result.job_id)       # "12345.gadi-pbs"
    print(result.virtual_id)   # "qx-550e8400-..."
"""

import logging
import os
import shlex
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

logger = logging.getLogger(__name__)


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
                      The scheduler command that was run, shell-quoted
                      (useful for debugging).
    """

    job_id: str
    virtual_id: Optional[str] = None
    joblog_path: Optional[str] = None
    stdout_path: Optional[str] = None
    stderr_path: Optional[str] = None
    submission_command: str = ""


def exec_context_summary(context_type: str, context_value: Any) -> Dict[str, Any]:
    """How the queue table records an execution context."""
    if isinstance(context_value, dict):
        value = (
            f"{context_value.get('env', '')} "
            f"+ {' '.join(context_value.get('modules', []))}"
        )
    elif isinstance(context_value, list):
        value = " ".join(context_value)
    else:
        value = context_value
    return {"type": context_type, "value": value}


def dispatch(
    spec,
    *,
    run_id: str,
    run_path,
    command: str,
    tags: Sequence[str],
    working_dir: str,
    exec_context: Dict[str, Any],
    backend=None,
) -> Tuple[str, Optional[str]]:
    """Submit ``spec`` and record it: the one submission path.

    Creates the queue entry, submits through the platform's scheduler
    backend, then records the job ID in the queue table and the resource
    tracker. Returns ``(job_id, queue_entry_id)``; the entry ID is ``None``
    if the queue entry couldn't be written.

    If the scheduler refuses the job, the run file is removed, the queue
    entry is marked failed and :class:`~qxub.backends.SubmitError` is raised.
    """
    from .. import backends
    from ..queue import create_queue_entry, runs, update_queue_entry
    from ..queue.consolidate import maybe_start_in_background
    from ..remote import handshake
    from ..resources.tracker import resource_tracker

    backend = backend or backends.get_backend()
    # Run by a local qxub over SSH: tell it once the job is submitted
    announce = handshake.requested()
    cpus = spec.resources.cpus

    virtual_id = None
    try:
        virtual_id = create_queue_entry(
            command=command,
            pbs_job_id=None,
            tags=list(tags),
            working_dir=working_dir,
            exec_context=exec_context,
            joblog_path=spec.job_log,
            queue_name=spec.queue,
            cpus_requested=cpus,
            entry_id=run_id,
            backend=backend.name,
        )
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("Failed to create queue entry: %s", exc)

    try:
        job_id = backend.submit(spec)
    except BaseException:
        runs.discard(run_path)
        runs.discard(spec.env_file)
        if virtual_id:
            try:
                update_queue_entry(virtual_id, status="failed")
            except Exception:  # pylint: disable=broad-except
                pass
        raise
    if announce:
        handshake.announce(job_id, backend.name)

    if virtual_id:
        try:
            update_queue_entry(
                virtual_id,
                pbs_job_id=job_id,
                status="dispatched",
                joblog_path=spec.job_log,
            )
        except Exception as exc:  # pylint: disable=broad-except
            logger.debug("Failed to update queue entry: %s", exc)

    try:
        resource_tracker.log_job_submitted(
            job_id=job_id,
            command=command,
            tags=list(tags),
            joblog_path=spec.job_log,
            queue=spec.queue,
            cpus_requested=cpus,
            job_name=spec.name,
            backend=backend.name,
        )
    except Exception as exc:  # pylint: disable=broad-except
        logger.debug("Failed to log job submission: %s", exc)

    # Fold finished runs into the stats snapshot now and then (Phase 3).
    maybe_start_in_background(run_path.parent if run_path else None)

    return job_id, virtual_id


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
    resolves configuration, selects a template, and hands the job to
    :func:`dispatch` (the same pipeline ``qxub exec`` uses), which submits it
    through the platform's scheduler backend and records it.

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
        SubmitError:  If the scheduler refuses the job. ``QsubError`` is
                      the same class, so existing handlers still work.
        ValueError:   If multiple execution contexts are specified.
    """
    from .. import backends
    from ..config.handler import job_spec, process_job_options
    from ..core.templates import get_template
    from ..queue import runs
    from ..queue.db import new_virtual_id
    from .context import ExecutionContext

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
    # Resources: PBS names and values, e.g. {"ncpus": 4, "mem": "8GB"}
    # ------------------------------------------------------------------
    from ..backends.pbs import parse_pbs_resources

    resource_spec = parse_pbs_resources(
        f"{key}={value}" for key, value in (resources or {}).items()
    )

    # ------------------------------------------------------------------
    # Config override
    # ------------------------------------------------------------------
    if config_path:
        from ..config.manager import set_config_override

        set_config_override(config_path)

    from ..config import manager as config_mod
    from ..config.handler import resolve_resources

    config_manager = config_mod.config_manager

    # The configured defaults (defaults.mem, defaults.volumes, ...) apply
    # as they do for ``qxub exec``.
    resource_spec = resolve_resources(resource_spec, config_manager)

    # ------------------------------------------------------------------
    # Build params dict (mimics what exec_cli constructs)
    # ------------------------------------------------------------------
    params: Dict[str, Any] = {
        "command": cmd_parts,
        "resources": resource_spec,
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
    processed_params, _options = process_job_options(params, config_manager)

    # ------------------------------------------------------------------
    # Job variables (the same as ``qxub exec`` passes)
    # ------------------------------------------------------------------
    exec_dir = processed_params.get("execdir") or os.getcwd()
    job_vars = ExecutionContext(context_type, context_value, template_type).job_env(
        cmd_parts, processed_params, pre=pre, post=post, bind=bind
    )

    # Run file: the job updates it instead of writing the job database
    # (docs/dev/sqlite-robustness.md, Phase 2).
    run_id = new_virtual_id()
    backend = backends.get_backend()
    run_path = runs.write_submitted(
        run_id,
        {
            "backend": backend.name,
            "state": "submitted",
            "submitted_at": datetime.now().isoformat(timespec="seconds"),
            "user": os.environ.get("USER", ""),
            "submit_host": os.uname().nodename,
            "command_b64": job_vars["cmd_b64"],
            "working_dir": exec_dir,
            "exec_context": {"type": context_type, "value": context_value},
            "queue": processed_params.get("queue"),
            "job_name": processed_params.get("name"),
            "project": processed_params.get("project"),
            "resources": resource_spec.as_dict(),
            "tags": list(tags or []),
            "joblog": processed_params.get("joblog"),
            "out": processed_params.get("out"),
            "err": processed_params.get("err"),
        },
    )
    job_vars.update(
        runs.job_env(
            run_id,
            run_path,
            backend.job_id_variable,
            processed_params.get("name"),
        )
    )
    spec = runs.attach_settings(job_spec(processed_params, template, job_vars), run_id)
    submission_command = shlex.join(backend.submit_command(spec))
    logger.info("Programmatic submission command: %s", submission_command)

    # Raises SubmitError (also known as QsubError) if the job is refused
    job_id, virtual_id = dispatch(
        spec,
        run_id=run_id,
        run_path=run_path,
        command=cmd_str,
        tags=list(tags or []),
        working_dir=exec_dir,
        exec_context=exec_context_summary(context_type, context_value),
        backend=backend,
    )

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
