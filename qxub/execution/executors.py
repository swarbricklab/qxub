"""
Execution context handlers for qxub.

This module contains the simplified execution functions that use the
common job submission logic from the execution module.
"""

from typing import List, Optional, Tuple

import click

from ..core.templates import (
    get_conda_template,
    get_default_template,
    get_module_template,
    get_singularity_template,
)
from .core import submit_and_monitor_job


def _handle_terse_execution(ctx, command, template_path, context_vars, pre, post):
    """Handle terse mode execution for all executor functions."""
    ctx_obj = ctx.obj

    # Handle dry run for terse mode
    if ctx_obj["dry"]:
        click.echo("DRY_RUN")
        return

    # Import required functions
    from pathlib import Path

    from ..core.scheduler import qsub
    from .core import build_submission_variables

    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

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
    submission_command = (
        f'qsub -v {submission_vars} {ctx_obj["options"]} {template_path}'
    )

    # Submit job and emit job ID
    job_id = qsub(submission_command, quiet=True)  # Always quiet for terse
    click.echo(job_id)


def execute_conda(
    ctx,
    command: Tuple[str, ...],
    env: str,
    template: Optional[str] = None,
    pre: Optional[str] = None,
    post: Optional[str] = None,
) -> None:
    """Execute command in conda environment."""
    if not env:
        click.echo(
            "Error: No conda environment available. Either activate a conda environment or use --env to specify one.",
            err=True,
        )
        ctx.exit(2)

    # Use default template if not provided
    template_path = get_conda_template(template)

    # Set up context variables for conda execution
    context_vars = {"env": env}

    # Handle terse mode
    if ctx.obj.get("terse", False):
        return _handle_terse_execution(
            ctx, command, template_path, context_vars, pre, post
        )

    # Submit and monitor the job
    submit_and_monitor_job(ctx, command, template_path, context_vars, pre, post)


def execute_module(
    ctx,
    command: Tuple[str, ...],
    modules: List[str],
    template: Optional[str] = None,
    pre: Optional[str] = None,
    post: Optional[str] = None,
) -> None:
    """Execute command with environment modules."""
    if not modules:
        click.echo("Error: No environment modules specified.", err=True)
        ctx.exit(2)

    # Use default template if not provided
    template_path = get_module_template(template)

    # Set up context variables for module execution
    context_vars = {"mods": modules}

    # Handle terse mode
    if ctx.obj.get("terse", False):
        return _handle_terse_execution(
            ctx, command, template_path, context_vars, pre, post
        )

    # Submit and monitor the job
    submit_and_monitor_job(ctx, command, template_path, context_vars, pre, post)


def execute_singularity(
    ctx,
    command: Tuple[str, ...],
    container: str,
    bind: Optional[str] = None,
    template: Optional[str] = None,
    pre: Optional[str] = None,
    post: Optional[str] = None,
) -> None:
    """Execute command in Singularity container."""
    if not container:
        click.echo("Error: Singularity container path is required.", err=True)
        ctx.exit(2)

    # Use default template if not provided
    template_path = get_singularity_template(template)

    # Set up context variables for singularity execution
    context_vars = {"sif": container}
    if bind:
        context_vars["bind"] = bind

    # Handle terse mode
    if ctx.obj.get("terse", False):
        return _handle_terse_execution(
            ctx, command, template_path, context_vars, pre, post
        )

    # Submit and monitor the job
    submit_and_monitor_job(ctx, command, template_path, context_vars, pre, post)


def execute_default(
    ctx,
    command: Tuple[str, ...],
    template: Optional[str] = None,
    pre: Optional[str] = None,
    post: Optional[str] = None,
) -> None:
    """Execute command using default template (no environment activation)."""
    # Use default template if not provided
    template_path = get_default_template(template)

    # No context variables needed for default execution
    context_vars = {}

    # Submit and monitor the job
    # Handle terse mode
    if ctx.obj.get("terse", False):
        return _handle_terse_execution(
            ctx, command, template_path, context_vars, pre, post
        )

    # Normal mode: use the full monitoring
    submit_and_monitor_job(ctx, command, template_path, context_vars, pre, post)
