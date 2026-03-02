"""
Unified execution context handling for qxub.

This module consolidates the execution logic that was previously
duplicated across execute_conda, execute_module, execute_singularity,
and execute_default functions.
"""

import base64
import logging
import os
from pathlib import Path
from typing import List, Optional, Union

import click

from ..core.scheduler import QsubError, qsub
from ..history import history_manager
from ..queue import create_queue_entry, update_queue_entry
from ..queue.db import get_db_path
from ..resources import resource_tracker


def _parse_cpus_from_resources(resources: list) -> Optional[int]:
    """Extract ncpus value from a PBS resources list like ['ncpus=2', 'mem=4gb']."""
    import re  # pylint: disable=import-outside-toplevel

    for item in resources or []:
        for part in str(item).split(","):
            m = re.match(r"ncpus=(\d+)", part.strip())
            if m:
                return int(m.group(1))
    return None


class ExecutionContext:
    """Represents an execution context with its specific configuration."""

    def __init__(
        self,
        context_type: str,
        context_value: Union[str, List[str]],
        template_type: str,
    ):
        self.context_type = context_type  # 'conda', 'module', 'singularity', 'default'
        self.context_value = (
            context_value  # env name, module list, container path, or None
        )
        self.template_type = template_type

    def validate(self) -> None:
        """Validate the execution context."""
        if self.context_type == "conda" and not self.context_value:
            raise click.ClickException(
                "Error: No conda environment available. Either activate a conda environment or use --env to specify one."
            )
        elif self.context_type == "module" and not self.context_value:
            raise click.ClickException("Error: No environment modules specified.")
        elif self.context_type == "singularity" and not self.context_value:
            raise click.ClickException("Error: Singularity container path is required.")

    def get_default_template(self) -> str:
        """Get the default template for this execution context."""
        from ..core.templates import get_template

        return get_template(self.template_type)

    def build_submission_vars(
        self,
        command: tuple,
        ctx_obj: dict,
        pre: str = None,
        post: str = None,
        bind: str = None,
    ) -> str:
        """Build the submission variables string for qsub."""
        cmd_str = " ".join(command)
        cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")

        # Base submission variables
        out = Path(ctx_obj["out"])
        err = Path(ctx_obj["err"])
        # Default to current working directory if execdir is None
        exec_dir = ctx_obj["execdir"] or os.getcwd()
        create_execdir = str(ctx_obj.get("create_execdir", False)).lower()
        base_vars = f'cmd_b64="{cmd_b64}",cwd={exec_dir},out={out},err={err},quiet={str(ctx_obj["quiet"]).lower()},create_execdir={create_execdir}'

        # Context-specific variables
        if self.context_type == "conda":
            context_vars = f"env={self.context_value}"
        elif self.context_type == "module":
            mods_str = " ".join(self.context_value)
            context_vars = f'mods="{mods_str}"'
        elif self.context_type == "singularity":
            context_vars = f"sif={self.context_value}"
            if bind:
                context_vars += f',bind="{bind}"'
        else:  # default
            context_vars = ""

        # Combine base and context variables
        if context_vars:
            submission_vars = f"{context_vars},{base_vars}"
        else:
            submission_vars = base_vars

        # Add pre/post commands if provided
        if pre:
            pre_b64 = base64.b64encode(pre.encode("utf-8")).decode("ascii")
            submission_vars += f',pre_cmd_b64="{pre_b64}"'
        if post:
            post_b64 = base64.b64encode(post.encode("utf-8")).decode("ascii")
            submission_vars += f',post_cmd_b64="{post_b64}"'

        # Pass the resolved DB path so job scripts can write status to the right DB
        db_path = str(get_db_path())
        submission_vars += f",QXUB_SHARED_DB={db_path}"

        return submission_vars

    def log_debug_info(
        self,
        ctx_obj: dict,
        command: tuple,
        template: str,
        pre: str = None,
        post: str = None,
    ) -> None:
        """Log debug information for this execution context."""
        # Log context parameters
        for key, value in ctx_obj.items():
            logging.debug("Context: %s = %s", key, value)

        # Log context-specific information
        if self.context_type == "conda":
            logging.debug("Conda environment: %s", self.context_value)
        elif self.context_type == "module":
            logging.debug("Environment modules: %s", self.context_value)
        elif self.context_type == "singularity":
            logging.debug("Singularity container: %s", self.context_value)
        else:
            logging.debug("Default execution context")

        logging.debug("Jobscript template: %s", template)
        logging.debug("Command: %s", command)
        logging.debug("Pre-commands: %s", pre)
        logging.debug("Post-commands: %s", post)


def execute_unified(
    ctx,
    command: tuple,
    execution_context: ExecutionContext,
    template: str = None,
    pre: str = None,
    post: str = None,
    bind: str = None,
) -> None:
    """Unified execution function for all execution contexts.

    Args:
        ctx: Click context
        command: Command tuple to execute
        execution_context: ExecutionContext instance
        template: Optional template override
        pre: Pre-execution command
        post: Post-execution command
        bind: Singularity bind mounts (only used for singularity context)
    """
    # Validate execution context
    execution_context.validate()

    # Use default template if not provided
    if not template:
        template = execution_context.get_default_template()

    ctx_obj = ctx.obj

    # Log debug information
    execution_context.log_debug_info(ctx_obj, command, template, pre, post)

    # Build submission variables
    submission_vars = execution_context.build_submission_vars(
        command, ctx_obj, pre, post, bind
    )

    # Build final submission command
    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'
    logging.info("Submission command: %s", submission_command)

    # Progress message: Job command constructed (skip for terse mode)
    if not ctx_obj["quiet"] and not ctx_obj.get("terse", False):
        from ..core.scheduler import print_status

        print_status("🔧 Job command constructed", final=True)

    if ctx_obj["dry"]:
        print(f"Dry run - would execute: {submission_command}")
        # Log history even for dry runs
        try:
            # Prepare file paths for history (even for dry runs)
            file_paths = {
                "out": ctx_obj["out"],
                "err": ctx_obj["err"],
                "joblog": ctx_obj.get("joblog"),  # May be None initially
            }
            history_manager.log_execution(
                ctx,
                success=True,
                file_paths=file_paths,
                tags=ctx_obj.get("tags") or [],
            )
        except Exception as e:
            logging.debug("Failed to log execution history: %s", e)
        return

    # ------------------------------------------------------------------
    # Pre-qsub: create queue entry with virtual ID (status = initiated)
    # ------------------------------------------------------------------
    cmd_str = " ".join(command)
    tags = ctx_obj.get("tags") or []
    exec_context_info = {
        "type": execution_context.context_type,
        "value": (
            execution_context.context_value
            if not isinstance(execution_context.context_value, list)
            else " ".join(execution_context.context_value)
        ),
    }

    virtual_id = None
    try:
        virtual_id = create_queue_entry(
            command=cmd_str,
            pbs_job_id=None,
            tags=tags,
            working_dir=ctx_obj.get("execdir") or os.getcwd(),
            exec_context=exec_context_info,
            joblog_path=ctx_obj.get("joblog"),
            queue_name=ctx_obj.get("queue"),
            cpus_requested=_parse_cpus_from_resources(ctx_obj.get("resources", [])),
        )
    except Exception as e:
        logging.debug("Failed to create queue entry: %s", e)

    # ------------------------------------------------------------------
    # Submit job to PBS
    # ------------------------------------------------------------------
    try:
        job_id = qsub(submission_command, quiet=ctx_obj["quiet"])
    except QsubError as qsub_exc:
        # qsub failed — mark entry as failed and re-raise
        if virtual_id:
            try:
                update_queue_entry(virtual_id, status="failed")
            except Exception:
                pass
        # In CLI mode, display the error and exit with the PBS exit code
        click.echo(qsub_exc.stderr if qsub_exc.stderr else str(qsub_exc))
        import sys as _sys

        _sys.exit(qsub_exc.exit_code)
    except Exception as exc:
        # Unexpected error — mark entry as failed and re-raise
        if virtual_id:
            try:
                update_queue_entry(virtual_id, status="failed")
            except Exception:
                pass
        raise exc

    # ------------------------------------------------------------------
    # Post-qsub: fill in real PBS job ID (status → dispatched)
    # ------------------------------------------------------------------
    if virtual_id:
        try:
            update_queue_entry(
                virtual_id,
                pbs_job_id=job_id,
                status="dispatched",
                joblog_path=ctx_obj.get("joblog"),
            )
        except Exception as e:
            logging.debug("Failed to update queue entry: %s", e)

    # Log job submission for status and resource tracking (do this before terse return)
    try:
        resource_tracker.log_job_submitted(
            job_id=job_id,
            command=cmd_str,
            tags=tags,
            joblog_path=ctx_obj.get("joblog"),
            queue=ctx_obj.get("queue"),
            cpus_requested=_parse_cpus_from_resources(ctx_obj.get("resources", [])),
        )
    except Exception as e:
        logging.debug("Failed to log job submission: %s", e)

    # Log execution to history system (must happen before terse mode early return)
    try:
        # Prepare file paths for history
        file_paths = {
            "out": ctx_obj["out"],
            "err": ctx_obj["err"],
            "joblog": ctx_obj.get("joblog"),  # May be None initially
        }
        history_manager.log_execution(
            ctx,
            success=True,
            job_id=job_id,
            file_paths=file_paths,
            tags=ctx_obj.get("tags") or [],
        )
    except Exception as e:
        logging.debug("Failed to log execution history: %s", e)

    # Handle terse mode - emit real PBS job ID and return immediately (for pipeline use)
    # Virtual IDs will be emitted here once Phase 3 introduces pending-dispatch.
    if ctx_obj.get("terse", False):
        click.echo(job_id)
        logging.info(
            "Terse mode: emitted job ID %s and returning immediately",
            job_id,
        )
        return  # Return immediately for terse mode
    elif not ctx_obj["quiet"]:
        # Display job ID to user (only in normal mode, not quiet or terse)
        from ..core.scheduler import print_status

        success_msg = f"✅ Job submitted successfully! Job ID: {job_id}"
        print_status(success_msg, final=True)

    # All modes now continue to monitoring (removed the quiet mode early return)

    # Start concurrent monitoring of job and log files
    # Stream log files to STDOUT/STDERR as appropriate
    import sys
    from pathlib import Path

    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    out.parent.mkdir(parents=True, exist_ok=True)
    err.parent.mkdir(parents=True, exist_ok=True)
    out.touch()
    err.touch()

    # Use single-thread monitoring for simpler, more reliable operation
    from ..core.scheduler import monitor_job_single_thread

    try:
        # Terse mode should run monitoring silently like quiet mode
        monitoring_quiet = ctx_obj["quiet"] or ctx_obj.get("terse", False)
        verbose_level = ctx_obj.get("verbose", 0)
        exit_status = monitor_job_single_thread(
            job_id,
            out,
            err,
            quiet=monitoring_quiet,
            joblog_file=ctx_obj.get("joblog"),
            verbose=verbose_level,
            walltime_str=ctx_obj.get("runtime"),
            joblog_check_interval=ctx_obj.get("joblog_check_interval", 60),
            walltime_offset_sec=ctx_obj.get("walltime_offset_sec", 0),
        )
        # Exit with the job's exit status
        sys.exit(exit_status)
    except KeyboardInterrupt:
        # Handle Ctrl-C gracefully
        print("\n🛑 Interrupted! Cleaning up job...")
        from ..core.scheduler import qdel

        success = qdel(job_id, quiet=False)
        if success:
            print("✅ Job cleanup completed")
            sys.exit(130)  # Standard exit code for SIGINT
        else:
            print(
                f"⚠️  Job cleanup failed - you may need to manually run: qdel {job_id}"
            )
            sys.exit(1)


# Convenience functions for creating execution contexts
def create_conda_context(env: str) -> ExecutionContext:
    """Create a conda execution context."""
    return ExecutionContext("conda", env, "conda")


def create_module_context(modules: List[str]) -> ExecutionContext:
    """Create a module execution context."""
    return ExecutionContext("module", modules, "module")


def create_singularity_context(container: str) -> ExecutionContext:
    """Create a singularity execution context."""
    return ExecutionContext("singularity", container, "singularity")


def create_default_context() -> ExecutionContext:
    """Create a default execution context."""
    return ExecutionContext("default", None, "default")
