"""
Job execution logic for qxub.

This module contains the common job submission and monitoring logic
that was previously duplicated across multiple execution functions.
"""

import base64
import logging
import signal
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click

from .history_manager import history_manager
from .resource_tracker import resource_tracker
from .scheduler import monitor_and_tail, print_status, qsub, start_job_monitoring

# Global variable to track current job for signal handling
_CURRENT_JOB_ID = None  # pylint: disable=invalid-name


def _signal_handler(signum, frame):
    """Handle SIGINT (Ctrl+C) by cleaning up submitted job"""
    # pylint: disable=unused-argument
    if _CURRENT_JOB_ID:
        # Clear the current line completely before printing cleanup messages
        print(" " * 100, end="", flush=True)
        print("🛑 Interrupted! Cleaning up job...")
        from .scheduler import qdel

        success = qdel(_CURRENT_JOB_ID, quiet=False)
        if success:
            print("✅ Job cleanup completed")
        else:
            print(
                f"⚠️  Job cleanup failed - you may need to manually run: qdel {_CURRENT_JOB_ID}"
            )
    print("👋 Goodbye!")
    sys.exit(130)  # Standard exit code for SIGINT


def build_submission_variables(
    command: Tuple[str, ...],
    execdir: str,
    out: Path,
    err: Path,
    quiet: bool,
    pre: Optional[str] = None,
    post: Optional[str] = None,
    **context_vars,
) -> str:
    """
    Build submission variables string for qsub command.

    Args:
        command: Command tuple to execute
        execdir: Execution directory
        out: Output file path
        err: Error file path
        quiet: Whether to run in quiet mode
        pre: Pre-command to run
        post: Post-command to run
        **context_vars: Additional context-specific variables (env, mods, sif, bind, etc.)

    Returns:
        Formatted submission variables string
    """
    # Base64 encode the command to avoid escaping issues
    cmd_str = " ".join(command)
    cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")

    # Build base submission variables
    submission_vars = (
        f'cmd_b64="{cmd_b64}",cwd={execdir},'
        f"out={out},err={err},quiet={str(quiet).lower()}"
    )

    # Add context-specific variables
    for key, value in context_vars.items():
        if value is not None:
            if key in ["env", "sif", "bind"]:
                # String values that need quotes
                submission_vars += f',{key}="{value}"'
            elif key == "mods":
                # Module list needs to be space-separated and quoted
                if isinstance(value, list):
                    mods_str = " ".join(value)
                    submission_vars += f',mods="{mods_str}"'
                else:
                    submission_vars += f',mods="{value}"'

    # Add pre and post commands if provided
    if pre:
        pre_b64 = base64.b64encode(pre.encode("utf-8")).decode("ascii")
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_b64 = base64.b64encode(post.encode("utf-8")).decode("ascii")
        submission_vars += f',post_cmd_b64="{post_b64}"'

    return submission_vars


def submit_and_monitor_job(
    ctx,
    command: Tuple[str, ...],
    template: str,
    context_vars: Optional[Dict] = None,
    pre: Optional[str] = None,
    post: Optional[str] = None,
) -> None:
    """
    Submit a job and monitor its execution.

    This function contains the common logic for job submission and monitoring
    that was previously duplicated across all execution functions.

    Args:
        ctx: Click context object
        command: Command tuple to execute
        template: Path to job script template
        context_vars: Context-specific variables (env, mods, sif, bind, etc.)
        pre: Pre-command to run
        post: Post-command to run
    """
    if context_vars is None:
        context_vars = {}

    # Declare global variables
    global _CURRENT_JOB_ID  # pylint: disable=global-statement

    ctx_obj = ctx.obj
    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Template: %s", template)
    logging.debug("Command: %s", command)
    logging.debug("Context vars: %s", context_vars)
    logging.debug("Pre-commands: %s", pre)
    logging.debug("Post-commands: %s", post)

    # Build submission variables
    submission_vars = build_submission_variables(
        command=command,
        execdir=ctx_obj["execdir"],
        out=out,
        err=err,
        quiet=ctx_obj["quiet"],
        pre=pre,
        post=post,
        **context_vars,
    )

    # Construct qsub command
    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'
    logging.info("Submission command: %s", submission_command)

    # Progress message: Job command constructed
    if not ctx_obj["quiet"]:
        print_status("🔧 Job command constructed", final=True)

    # Handle dry run
    if ctx_obj["dry"]:
        print(f"Dry run - would execute: {submission_command}")
        # Log history even for dry runs
        try:
            history_manager.log_execution(ctx, success=True)
        except Exception as e:
            logging.debug("Failed to log execution history: %s", e)
        return

    # Submit job
    job_id = qsub(submission_command, quiet=ctx_obj["quiet"])

    # Track job ID globally for signal handler
    _CURRENT_JOB_ID = job_id

    # Display job ID to user (unless in quiet mode)
    if not ctx_obj["quiet"]:
        success_msg = f"✅ Job submitted successfully! Job ID: {job_id}"
        print_status(success_msg, final=True)

    # Log execution to history system
    try:
        history_manager.log_execution(ctx, success=True, job_id=job_id)
    except Exception as e:
        logging.debug("Failed to log execution history: %s", e)

    # Log job execution for resource tracking
    try:
        cmd_str = " ".join(command)
        resource_tracker.log_job_resources(
            job_id=job_id,
            resource_data={},  # Will be populated when job completes
            command=cmd_str,
        )
    except Exception as e:
        logging.debug("Failed to log job resources: %s", e)

    # Exit if in quiet mode
    if ctx_obj["quiet"]:
        logging.info("Exiting in quiet mode")
        return

    # Register signal handler for Ctrl+C cleanup
    signal.signal(signal.SIGINT, _signal_handler)

    # Start concurrent monitoring of job and log files
    # Stream log files to STDOUT/STDERR as appropriate
    out.parent.mkdir(parents=True, exist_ok=True)
    err.parent.mkdir(parents=True, exist_ok=True)
    out.touch()
    err.touch()

    # Start job monitoring and get coordinator for signaling
    coordinator, wait_for_completion = start_job_monitoring(
        job_id, out, err, quiet=ctx_obj["quiet"]
    )

    # Signal that submission messages are complete - spinner can start now
    coordinator.signal_submission_complete()

    try:
        exit_status = wait_for_completion()
        # Exit with the job's exit status
        sys.exit(exit_status)
    finally:
        # Clear job ID when monitoring completes (successfully or via interrupt)
        _CURRENT_JOB_ID = None


def validate_execution_context(
    default_flag, conda_env, module_list, container
) -> Tuple[bool, str]:
    """
    Validate execution context parameters.

    Args:
        default_flag: Default execution flag
        conda_env: Conda environment name
        module_list: List of modules to load
        container: Singularity container path

    Returns:
        Tuple of (has_context, context_type)

    Raises:
        click.ClickException: If multiple execution contexts are specified
    """
    execution_contexts = [default_flag, conda_env, module_list, container]
    context_count = sum(bool(x) for x in execution_contexts)

    if context_count > 1:
        raise click.ClickException("Cannot specify multiple execution contexts")

    has_context = context_count > 0

    # Determine context type
    if default_flag:
        return has_context, "default"
    elif conda_env:
        return has_context, "conda"
    elif module_list:
        return has_context, "module"
    elif container:
        return has_context, "singularity"
    else:
        return has_context, "default"
