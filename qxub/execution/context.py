"""
Unified execution context handling for qxub.

This module consolidates the execution logic that was previously
duplicated across execute_conda, execute_module, execute_singularity,
and execute_default functions.
"""

import base64
import logging
import os
import shlex
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Union

import click

logger = logging.getLogger(__name__)

from .. import backends
from ..backends import SubmitError
from ..config.handler import job_spec
from ..history import history_manager
from ..queue import runs
from ..queue.db import new_virtual_id
from .submit import dispatch, exec_context_summary


def _b64(text: str) -> str:
    return base64.b64encode(text.encode("utf-8")).decode("ascii")


class ExecutionContext:
    """Represents an execution context with its specific configuration."""

    def __init__(
        self,
        context_type: str,
        context_value: Union[str, List[str], dict, None],
        template_type: str,
    ):
        self.context_type = (
            context_type  # 'conda', 'module', 'singularity', 'conda_module', 'default'
        )
        self.context_value = (
            context_value  # env name, module list, container path, dict, or None
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
        elif self.context_type == "conda_module":
            if not self.context_value or not self.context_value.get("env"):
                raise click.ClickException(
                    "Error: No conda environment specified for combined context."
                )
            if not self.context_value.get("modules"):
                raise click.ClickException(
                    "Error: No modules specified for combined context."
                )

    def get_default_template(self) -> str:
        """Get the default template for this execution context."""
        from ..core.templates import get_template

        return get_template(self.template_type)

    def job_env(
        self,
        command: tuple,
        ctx_obj: dict,
        pre: str = None,
        post: str = None,
        bind: str = None,
        user_vars: tuple = None,
    ) -> Dict[str, str]:
        """The variables the job script reads, in the order they are passed.

        Args:
            command: Command tuple to execute
            ctx_obj: Job options (``out``, ``err``, ``execdir``, ``quiet``, ...)
            pre: Pre-command to run before main command
            post: Post-command to run after main command
            bind: Singularity bind mounts
            user_vars: User-defined environment variables (KEY=VALUE tuples)
        """
        from .. import __version__ as qxub_version

        env: Dict[str, str] = {}

        # Context-specific variables
        if self.context_type == "conda":
            env["env"] = self.context_value
        elif self.context_type == "module":
            env["mods"] = " ".join(self.context_value)
        elif self.context_type == "conda_module":
            env["env"] = self.context_value["env"]
            env["mods"] = " ".join(self.context_value["modules"])
        elif self.context_type == "singularity":
            env["sif"] = str(self.context_value)
            if bind:
                env["bind"] = bind

        env.update(
            {
                "cmd_b64": _b64(" ".join(command)),
                # Default to current working directory if execdir is None
                "cwd": str(ctx_obj.get("execdir") or os.getcwd()),
                "out": str(ctx_obj["out"]),
                "err": str(ctx_obj["err"]),
                "quiet": str(bool(ctx_obj.get("quiet"))).lower(),
                "create_execdir": str(ctx_obj.get("create_execdir", False)).lower(),
                # The submitting host's qxub, for the job to report back with
                "qxub_version": qxub_version,
                "qxub_path": shutil.which("qxub") or "unknown",
                "qxub_python": sys.executable,
                "qxub_submit_host": os.uname().nodename,
            }
        )
        if pre:
            env["pre_cmd_b64"] = _b64(pre)
        if post:
            env["post_cmd_b64"] = _b64(post)
        if user_vars:
            # KEY=VALUE pairs joined with semicolons, base64-encoded so no
            # shell or scheduler quoting is needed
            env["user_vars_b64"] = _b64(";".join(user_vars))

        # Slack/Discord notifications, if configured
        if ctx_obj.get("notify") and self._has_push_notifications(ctx_obj):
            env["qxub_notify"] = "true"
            env["qxub_notify_outdir"] = str(Path(ctx_obj.get("out", "")).parent)

        return env

    def _has_push_notifications(self, ctx_obj: dict) -> bool:
        """Check if Slack/Discord notifications are configured."""
        try:
            from ..config import manager as config_mod

            config_manager = config_mod.config_manager
            slack_webhook = config_manager.get_config_value(
                "notifications.slack.webhook_url"
            )
            discord_webhook = config_manager.get_config_value(
                "notifications.discord.webhook_url"
            )
            return bool(slack_webhook or discord_webhook)
        except Exception:
            return False

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
            logger.debug("Context: %s = %s", key, value)

        # Log context-specific information
        if self.context_type == "conda":
            logger.debug("Conda environment: %s", self.context_value)
        elif self.context_type == "module":
            logger.debug("Environment modules: %s", self.context_value)
        elif self.context_type == "conda_module":
            logger.debug(
                "Combined context - Conda: %s, Modules: %s",
                self.context_value["env"],
                self.context_value["modules"],
            )
        elif self.context_type == "singularity":
            logger.debug("Singularity container: %s", self.context_value)
        else:
            logger.debug("Default execution context")

        logger.debug("Jobscript template: %s", template)
        logger.debug("Command: %s", command)
        logger.debug("Pre-commands: %s", pre)
        logger.debug("Post-commands: %s", post)


def _run_record_fields(ctx_obj: dict, command, execution_context) -> dict:
    """What the run file records at submission (see ``qxub.queue.runs``)."""
    from datetime import datetime

    from ..config.handler import as_resource_spec

    cmd_str = " ".join(command)
    return {
        "state": "submitted",
        "submitted_at": datetime.now().isoformat(timespec="seconds"),
        "user": os.environ.get("USER", ""),
        "submit_host": os.uname().nodename,
        "command_b64": base64.b64encode(cmd_str.encode("utf-8")).decode("ascii"),
        "working_dir": ctx_obj.get("execdir") or os.getcwd(),
        "exec_context": {
            "type": execution_context.context_type,
            "value": execution_context.context_value,
        },
        "queue": ctx_obj.get("queue"),
        "job_name": ctx_obj.get("name"),
        "project": ctx_obj.get("project"),
        "resources": as_resource_spec(ctx_obj.get("resources")).as_dict(),
        "tags": list(ctx_obj.get("tags") or []),
        "joblog": ctx_obj.get("joblog"),
        "out": ctx_obj.get("out"),
        "err": ctx_obj.get("err"),
    }


def execute_unified(
    ctx,
    command: tuple,
    execution_context: ExecutionContext,
    template: str = None,
    pre: str = None,
    post: str = None,
    bind: str = None,
    user_vars: tuple = None,
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
        user_vars: User-defined environment variables (KEY=VALUE tuples)
    """
    # Validate execution context
    execution_context.validate()

    # Use default template if not provided
    if not template:
        template = execution_context.get_default_template()

    ctx_obj = ctx.obj

    # Log debug information
    execution_context.log_debug_info(ctx_obj, command, template, pre, post)

    # The variables the job script reads
    job_vars = execution_context.job_env(command, ctx_obj, pre, post, bind, user_vars)

    # The run ID names this run's record file, which the job updates instead
    # of writing the job database (docs/dev/sqlite-robustness.md, Phase 2).
    run_id = new_virtual_id()
    run_path = None
    if not ctx_obj["dry"]:
        run_path = runs.write_submitted(
            run_id,
            _run_record_fields(ctx_obj, command, execution_context),
        )
    backend = backends.get_backend()
    job_vars.update(runs.job_env(run_id, run_path, backend.job_id_variable))
    spec = runs.attach_settings(
        job_spec(ctx_obj, template, job_vars), run_id, write=not ctx_obj["dry"]
    )
    submission_command = shlex.join(backend.submit_command(spec))
    logger.info("Submission command: %s", submission_command)

    # Progress message: Job command constructed (skip for terse mode)
    if not ctx_obj["quiet"] and not ctx_obj.get("terse", False):
        from ..core.scheduler import print_status

        print_status("🔧 Job command constructed", final=True)

    if ctx_obj["dry"]:
        print(f"Dry run - would execute: {submission_command}")
        if spec.env_file:
            print(f"Job settings ({spec.env_file}):")
            print(runs.format_settings(spec.env), end="")
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
            logger.debug("Failed to log execution history: %s", e)
        return

    cmd_str = " ".join(command)
    try:
        job_id, _entry_id = dispatch(
            spec,
            run_id=run_id,
            run_path=run_path,
            command=cmd_str,
            tags=ctx_obj.get("tags") or [],
            working_dir=ctx_obj.get("execdir") or os.getcwd(),
            exec_context=exec_context_summary(
                execution_context.context_type, execution_context.context_value
            ),
            backend=backend,
        )
    except SubmitError as exc:
        # In CLI mode, show the scheduler's message and exit with its code
        click.echo(exc.stderr if exc.stderr else str(exc), err=True)
        sys.exit(exc.exit_code)

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
        logger.debug("Failed to log execution history: %s", e)

    # Handle terse mode - emit real PBS job ID and return immediately (for pipeline use)
    # Virtual IDs will be emitted here once Phase 3 introduces pending-dispatch.
    if ctx_obj.get("terse", False):
        click.echo(job_id)
        logger.info(
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
            run_id=run_id,
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
