"""
Job execution logic for qxub.

This module contains the common job submission and monitoring logic
that was previously duplicated across multiple execution functions.
"""

import base64
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import click

from ..core.scheduler import monitor_job_single_thread, print_status, qsub
from ..history import history_manager
from ..resources import resource_tracker


def expand_submission_variables(cmd_str: str) -> str:
    """
    Expand ${var} at submission time, convert ${{var}} to ${var} for execution.
    All other $ usage is preserved unchanged.

    Enhanced with smart quote processing: if command starts and ends with double quotes,
    applies automatic escaping for better readability.

    Args:
        cmd_str: Command string with potential variables

    Returns:
        Command string with submission variables expanded and execution variables prepared

    Examples:
        # Traditional processing (backward compatible)
        expand_submission_variables('echo "${USER} on ${{HOSTNAME}}"')
        # Returns: 'echo "jr9959 on ${HOSTNAME}"'

        # Smart quote processing (enhanced)
        expand_submission_variables('"find /path -exec echo \\"User: ${USER}\\" \\;"')
        # Returns: 'find /path -exec echo "User: jr9959" \\;'

        # Literal dollar preservation
        expand_submission_variables('"awk \\"{print \\\\$1}\\" file"')
        # Returns: 'awk "{print $1}" file'
    """

    # Check for smart quote processing (double-quote wrapped commands)
    if cmd_str.startswith('"') and cmd_str.endswith('"') and len(cmd_str) > 1:
        # Remove outer double quotes
        inner_cmd = cmd_str[1:-1]

        # Process escaped characters first
        inner_cmd = inner_cmd.replace('\\"', '"')  # \" -> "

        # Handle literal dollars: \$ -> temporary placeholder
        inner_cmd = inner_cmd.replace("\\$", "__QXUB_LITERAL_DOLLAR__")

        # Apply variable expansion to the processed command
        processed = _apply_variable_expansion(inner_cmd)

        # Restore literal dollars
        processed = processed.replace("__QXUB_LITERAL_DOLLAR__", "$")

        return processed

    # No smart quotes - use standard processing
    return _apply_variable_expansion(cmd_str)


def _apply_variable_expansion(cmd_str: str) -> str:
    """
    Internal function to apply variable expansion rules.
    Separated for use by both smart quote and standard processing.
    """
    import logging
    import os
    import re

    def replace_submission_var(match):
        """Replace ${var} with environment value."""
        var_name = match.group(1)
        value = os.environ.get(var_name)
        if value is not None:
            return value
        else:
            logging.warning(
                f"Submission variable ${{{var_name}}} not found in environment"
            )
            return match.group(0)  # Keep original ${var}

    def replace_execution_var(match):
        """Convert ${{var}} to ${var} for execution-time expansion."""
        var_name = match.group(1)
        return f"${{{var_name}}}"

    # Process execution variables first to avoid conflicts
    result = re.sub(r"\$\{\{([^}]+)\}\}", replace_execution_var, cmd_str)

    # Then process submission variables
    result = re.sub(r"\$\{([^}]+)\}", replace_submission_var, result)

    return result


def expand_variables_preserving_quotes(cmd_str: str, is_list: bool = False) -> str:
    """
    Expand shell variables while preserving quote structure.

    This function expands variables like $USER and ${HOME} at submission time
    while preserving the exact quoting structure the user intended.
    It does not use shell evaluation to avoid quote mangling.

    Args:
        cmd_str: Command string potentially containing variables

    Returns:
        Command string with variables expanded but quotes preserved

    Examples:
        expand_variables_preserving_quotes('echo "Hello $USER"')
        # Returns: 'echo "Hello jr9959"'

        expand_variables_preserving_quotes('python -c "print(\'$USER\')"')
        # Returns: 'python -c "print(\'jr9959\')"'

        expand_variables_preserving_quotes('echo "Job: \\$PBS_JOBID"')
        # Returns: 'echo "Job: $PBS_JOBID"' (escaped variable preserved)
    """
    # Handle escaped variables first: \$VAR should become $VAR (literal)
    # We need to do this before expanding other variables
    result = cmd_str

    # Replace \$VAR with a temporary placeholder to protect it from expansion
    escaped_vars = []
    escaped_pattern = r"\\(\$\{[^}]+\}|\$[A-Za-z_][A-Za-z0-9_]*)"

    def save_escaped_var(match):
        escaped_vars.append(match.group(1))  # Save the $VAR part
        return f"__ESCAPED_VAR_{len(escaped_vars) - 1}__"

    result = re.sub(escaped_pattern, save_escaped_var, result)

    # Now expand unescaped variables
    def replace_var(match):
        var_name = match.group(1) or match.group(2)  # Handle both ${VAR} and $VAR
        value = os.environ.get(var_name)

        if value is not None:
            return value
        else:
            # Variable not found - keep as-is for potential execution-time expansion
            return match.group(0)

    # Match ${VAR} and $VAR patterns
    pattern = r"\$\{([^}]+)\}|\$([A-Za-z_][A-Za-z0-9_]*)"
    result = re.sub(pattern, replace_var, result)

    # Restore escaped variables (without the escape backslash)
    for i, escaped_var in enumerate(escaped_vars):
        result = result.replace(f"__ESCAPED_VAR_{i}__", escaped_var)

    return result


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
    # Join command and expand variables while preserving quotes
    cmd_str = " ".join(command)
    cmd_expanded = expand_submission_variables(cmd_str)

    # Base64 encode the expanded command
    cmd_b64 = base64.b64encode(cmd_expanded.encode("utf-8")).decode("ascii")

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
        pre_expanded = expand_submission_variables(pre)
        pre_b64 = base64.b64encode(pre_expanded.encode("utf-8")).decode("ascii")
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_expanded = expand_submission_variables(post)
        post_b64 = base64.b64encode(post_expanded.encode("utf-8")).decode("ascii")
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

    # Progress message: Job command constructed (skip for terse mode)
    if not ctx_obj["quiet"] and not ctx_obj.get("terse", False):
        print_status("🔧 Job command constructed", final=True)

    # Handle dry run
    if ctx_obj["dry"]:
        # Handle terse mode - just show "DRY_RUN" for terse
        if ctx_obj.get("terse", False):
            click.echo("DRY_RUN")
            return

        # Show expanded commands for user verification
        cmd_str = " ".join(command)
        cmd_expanded = expand_submission_variables(cmd_str)
        print(f"📝 Command to execute: {cmd_expanded}")

        if pre:
            pre_expanded = expand_submission_variables(pre)
            print(f"📋 Pre-command: {pre_expanded}")
        if post:
            post_expanded = expand_submission_variables(post)
            print(f"📋 Post-command: {post_expanded}")

        # Show full qsub command only if verbose
        verbose_level = ctx_obj.get("verbose", 0)
        if verbose_level > 0:
            print(f"🔧 Full qsub command: {submission_command}")
        else:
            print("Dry run - job would be submitted (use -v to see full qsub command)")

        # Log history even for dry runs
        try:
            history_manager.log_execution(ctx, success=True)
        except Exception as e:
            logging.debug("Failed to log execution history: %s", e)
        return

    # Submit job
    job_id = qsub(submission_command, quiet=ctx_obj["quiet"])

    # Handle terse mode - emit only job ID and return immediately
    if ctx_obj.get("terse", False):
        click.echo(job_id)
        logging.info("Terse mode: emitted job ID %s and exiting", job_id)
        return

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

    # Prepare output files
    out.parent.mkdir(parents=True, exist_ok=True)
    err.parent.mkdir(parents=True, exist_ok=True)
    out.touch()
    err.touch()

    # Use single-thread monitoring for simpler, more reliable operation
    try:
        exit_status = monitor_job_single_thread(
            job_id, out, err, quiet=ctx_obj["quiet"]
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
