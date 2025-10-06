"""
This package provides tools for running qsub jobs on PBS Pro with Singularity containers.
This avoids boilerplate code for setting up container environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

# pylint: disable=duplicate-code
import os
import sys
import signal
import logging
import base64
import shutil
from pathlib import Path
import click
import pkg_resources
from .cli import qxub
from .scheduler import qsub, monitor_and_tail, print_status, qdel

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
        template_path = pkg_resources.resource_filename(__name__, "jobscripts/qsing.pbs")
        if os.path.exists(template_path):
            return template_path
    except (ImportError, AttributeError, FileNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / "jobscripts" / "qsing.pbs"
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate qsing.pbs template file. "
        f"Looked in: {current_dir / 'jobscripts' / 'qsing.pbs'}"
    )


@qxub.command(
    context_settings={
        "ignore_unknown_options": True,
        "allow_extra_args": True,
        "allow_interspersed_args": False,
    }
)
@click.argument("cmd", nargs=-1)
@click.option("--sif", required=True, help="Path to the Singularity .sif container file")
@click.option(
    "--template",
    default=_get_default_template(),
    help="Jobscript template (optional - for further customization)",
)
@click.option(
    "--pre",
    help="Command to run before the singularity command " "(use quotes for commands with options)",
)
@click.option(
    "--post",
    help="Command to run after the singularity command "
    "(only if main command succeeds, use quotes for commands with options)",
)
@click.pass_context
def sing(
    ctx, cmd, sif, template, pre, post
):  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements
    """
    Constructs and submits a qsub job that will execute the given command
    in the specified Singularity container.

    Example:
        qxub --resources mem=50MB sing --sif container.sif -- python script.py

    Here the command 'python script.py' will run inside the Singularity container
    with 50MB of RAM. See 'qxub --help' for other resource options.

    Additional Singularity options can be passed before the '--':
        qxub sing --sif container.sif --bind /data:/data --env VAR=value -- command

    Note the '--' before the command. This is required to separate Singularity
    options from the command to run inside the container.
    """
    # Validate that the .sif file exists
    sif_path = Path(sif)
    if not sif_path.exists():
        raise click.ClickException(f"Singularity container file not found: {sif}")

    # Parse the command arguments to separate singularity options from container command
    # Everything before '--' should be singularity options
    # Everything after '--' should be the container command

    cmd_list = list(cmd)
    sing_options = []
    container_cmd = []

    if "--" in cmd_list:
        separator_index = cmd_list.index("--")
        sing_options = cmd_list[:separator_index]
        container_cmd = cmd_list[separator_index + 1 :]
    else:
        # If no '--' separator, everything is the container command
        container_cmd = cmd_list

    # Get values from context
    ctx_obj = ctx.obj
    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Singularity container: %s", sif)
    logging.debug("Singularity options: %s", sing_options)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Container command: %s", container_cmd)
    logging.debug("Pre-commands: %s", pre)
    logging.debug("Post-commands: %s", post)

    # Construct qsub command
    cmd_str = " ".join(container_cmd)
    sing_options_str = " ".join(sing_options) if sing_options else ""
    # Base64 encode the command to avoid escaping issues
    cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")
    submission_vars = (
        f'sif_file="{sif}",sing_options="{sing_options_str}",'
        f'cmd_b64="{cmd_b64}",cwd={ctx_obj["execdir"]},out={out},err={err},'
        f'quiet={str(ctx_obj["quiet"]).lower()}'
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
        success_msg = f"🚀 Job submitted successfully! Job ID: {job_id}"
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

        # Exit with the same status as the job
        if job_exit_status != 0:
            logging.info("Job failed with exit status %d", job_exit_status)
        else:
            logging.info("Job completed successfully")

        sys.exit(job_exit_status)

    finally:
        # Clear job ID when monitoring completes (successfully or via interrupt)
        _CURRENT_JOB_ID = None
