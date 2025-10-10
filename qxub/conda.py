"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

import base64
import difflib
import logging

# pylint: disable=duplicate-code
import os
import shutil
import signal
import sys
from pathlib import Path

import click
import pkg_resources

from .cli import qxub
from .resource_tracker import resource_tracker
from .scheduler import get_job_resource_data, monitor_and_tail, print_status, qdel, qsub


class QxubCommand(click.Command):
    """Custom Click command with enhanced error handling for unknown options."""

    def parse_args(self, ctx, args):
        """Override parse_args to provide better error messages for unknown options."""
        try:
            return super().parse_args(ctx, args)
        except click.NoSuchOption as e:
            self.handle_unknown_option_error(ctx, e, args)

    def handle_unknown_option_error(self, ctx, error, args):
        """Provide helpful suggestions for unknown options."""
        unknown_option = error.option_name

        # Get all valid options for this command
        command_options = []
        for param in self.params:
            if isinstance(param, click.Option):
                command_options.extend(param.opts)

        # Get parent qxub options
        parent_options = []
        if ctx.parent and ctx.parent.command:
            for param in ctx.parent.command.params:
                if isinstance(param, click.Option):
                    parent_options.extend(param.opts)

        # Find close matches
        close_matches = difflib.get_close_matches(
            unknown_option, command_options + parent_options, n=3, cutoff=0.6
        )

        # Build helpful error message
        click.echo(f"Error: No such option: {unknown_option}", err=True)

        if close_matches:
            click.echo("\n💡 Did you mean:", err=True)
            for match in close_matches:
                if match in command_options:
                    click.echo(f"   {match}  ({self.name} option)", err=True)
                else:
                    click.echo(
                        f"   {match}  (qxub option - put before '{self.name}')",
                        err=True,
                    )

        # Provide specific guidance for this subcommand
        click.echo(f"\n📖 Option placement for 'qxub {self.name}':", err=True)
        click.echo(
            f"   ✅ qxub --queue normal {self.name} --env myenv command", err=True
        )
        click.echo(
            f"   ❌ qxub {self.name} --queue normal --env myenv command", err=True
        )

        # Check if this looks like a command argument that should use --
        if unknown_option.startswith("-") and self._looks_like_command_argument(
            unknown_option, args
        ):
            click.echo(
                "\n   If this is part of your command, use '--' to separate:", err=True
            )
            click.echo(
                f"   ✅ qxub {self.name} --env myenv -- python script.py {unknown_option}",
                err=True,
            )
            click.echo(
                f"   ❌ qxub {self.name} --env myenv python script.py {unknown_option}",
                err=True,
            )

        # Show relevant help
        click.echo(f"\n🔍 For help: qxub {self.name} --help", err=True)

        ctx.exit(2)

    def _looks_like_command_argument(self, option, args):
        """Check if the unknown option looks like it should be part of a command."""
        # Common command-line options that users might accidentally use
        command_like_options = [
            "-c",
            "-i",
            "-o",
            "-f",
            "-d",
            "-r",
            "-t",
            "-s",
            "-n",
            "-p",
            "-h",
            "--input",
            "--output",
            "--file",
            "--config",
            "--help",
            "--version",
            "--debug",
            "--verbose",
            "--quiet",
            "--force",
            "--recursive",
        ]

        return option in command_like_options


# Global variable to track current job for signal handling
_CURRENT_JOB_ID = None  # pylint: disable=invalid-name


def _signal_handler(signum, frame):
    """Handle SIGINT (Ctrl+C) by cleaning up submitted job"""
    # pylint: disable=unused-argument
    if _CURRENT_JOB_ID:
        # Clear the current line completely before printing cleanup messages
        try:
            terminal_width = shutil.get_terminal_size().columns
            clear_width = terminal_width
        except Exception:
            clear_width = 100

        print("\r" + " " * clear_width + "\r", end="", flush=True)
        print("🛑 Interrupted! Cleaning up job...")
        success = qdel(_CURRENT_JOB_ID, quiet=False)
        if success:
            print("✅ Job cleanup completed")
        else:
            print(
                "⚠️  Job cleanup failed - you may need to manually run: qdel",
                _CURRENT_JOB_ID,
            )
    print("👋 Goodbye!")
    sys.exit(130)  # Standard exit code for SIGINT


def _get_default_template():
    """Get the default PBS template file path."""
    # Try pkg_resources first
    try:
        template_path = pkg_resources.resource_filename(
            __name__, "jobscripts/qconda.pbs"
        )
        if os.path.exists(template_path):
            return template_path
    except (ImportError, AttributeError, FileNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / "jobscripts" / "qconda.pbs"
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate qconda.pbs template file. "
        f"Looked in: {current_dir / 'jobscripts' / 'qconda.pbs'}"
    )


@qxub.command(cls=QxubCommand)
@click.argument("cmd", nargs=-1, required=False)
@click.option(
    "--env",
    default=os.getenv("CONDA_DEFAULT_ENV"),
    help="Conda environment to use (default: active environment)",
)
@click.option(
    "--template",
    default=_get_default_template(),
    help="Jobscript template (optional - for further customization)",
)
@click.option(
    "--pre",
    help="Command to run before the main command (use quotes for commands with options)",
)
@click.option(
    "--post",
    help="Command to run after the main command (only if main command "
    "succeeds, use quotes for commands with options)",
)
@click.pass_context
def conda(
    ctx, cmd, env, template, pre, post
):  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    """
    Constructs and submits a qsub job that will execute the given command
    in the specified conda environment and work directory

    Example:
        qxub --resources mem=50MB conda --env myenv -- python --version

    Here the command "python --version" will run in the "myenv" conda environment
    with 50MB of RAM. See "qxub --help" for other resource options.

    Note the "--" before the command. This is required if command itself also
    contains double-dashes, as in this case ("--version").
    """
    # Validate required parameters
    if not cmd:
        click.echo("Error: Command is required.", err=True)
        ctx.exit(2)

    if not env:
        click.echo(
            "Error: No conda environment available. Either activate a conda environment or use --env to specify one.",
            err=True,
        )
        ctx.exit(2)

    # Get values from context
    ctx_obj = ctx.obj
    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Conda environment: %s", env)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", cmd)
    logging.debug("Pre-commands: %s", pre)
    logging.debug("Post-commands: %s", post)

    # Construct qsub command
    cmd_str = " ".join(cmd)
    # Base64 encode the command to avoid escaping issues
    cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")
    submission_vars = (
        f'env={env},cmd_b64="{cmd_b64}",cwd={ctx_obj["execdir"]},'
        f'out={out},err={err},quiet={str(ctx_obj["quiet"]).lower()}'
    )
    if pre:
        pre_b64 = base64.b64encode(pre.encode("utf-8")).decode("ascii")
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_b64 = base64.b64encode(post.encode("utf-8")).decode("ascii")
        submission_vars += f',post_cmd_b64="{post_b64}"'

    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'

    # Execute the command
    logging.info("Submission command: %s", submission_command)
    if ctx_obj["dry"]:
        print(f"Dry run - would execute: {submission_command}")
        logging.info("Dry run. Exiting")
        sys.exit(0)

    # Register signal handler for Ctrl+C cleanup
    signal.signal(signal.SIGINT, _signal_handler)

    logging.info("Submitting job")
    job_id = qsub(submission_command, quiet=ctx_obj["quiet"])

    # Track job ID globally for signal handler
    global _CURRENT_JOB_ID  # pylint: disable=global-statement
    _CURRENT_JOB_ID = job_id

    # Show success message with job ID
    if not ctx_obj["quiet"]:
        success_msg = f"✓ Job {job_id[:8]} submitted"
        print_status(success_msg, final=False)

    logging.info("Your job has been successfully submitted")
    # Exit if in quiet mode
    if ctx_obj["quiet"]:
        logging.info("Exiting in quiet mode")
        sys.exit(0)

    # Start concurrent monitoring of job and log files
    # Stream log files to STDOUT/STDERR as appropriate
    out.parent.mkdir(parents=True, exist_ok=True)
    err.parent.mkdir(parents=True, exist_ok=True)
    out.touch()
    err.touch()

    # Pass success message to monitor_and_tail for spinner display
    if not ctx_obj["quiet"]:
        success_message = f"🚀 Job {job_id[:8]}"
    else:
        success_message = None

    try:
        job_exit_status = monitor_and_tail(
            job_id, out, err, quiet=ctx_obj["quiet"], success_msg=success_message
        )

        # Capture resource data for resource tracking (only if not a subprocess)
        is_subprocess = os.getenv("QXUB_SUBPROCESS")
        if not is_subprocess:
            try:
                resource_data = get_job_resource_data(job_id)
                if resource_data:
                    # Extract command for logging
                    command = " ".join(sys.argv)
                    resource_tracker.log_job_resources(job_id, resource_data, command)
                    logging.debug("Logged resource data for job %s", job_id)
                else:
                    logging.debug("No resource data available for job %s", job_id)
            except Exception as e:
                logging.debug("Failed to log resource data: %s", e)

        # Exit with the same status as the job
        if job_exit_status != 0:
            logging.info("Job failed with exit status %d", job_exit_status)
        else:
            logging.info("Job completed successfully")

        sys.exit(job_exit_status)

    finally:
        # Clear job ID when monitoring completes (successfully or via interrupt)
        _CURRENT_JOB_ID = None
