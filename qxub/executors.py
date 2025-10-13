"""
Execution context handlers for qxub.

This module contains the simplified execution functions that use the
common job submission logic from the execution module.
"""

from typing import List, Optional, Tuple

import click

from .execution import submit_and_monitor_job, validate_execution_context
from .templates import (
    get_conda_template,
    get_default_template,
    get_module_template,
    get_singularity_template,
)


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
    submit_and_monitor_job(ctx, command, template_path, context_vars, pre, post)
