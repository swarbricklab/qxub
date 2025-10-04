"""
This package provides tools for running qsub jobs on PBS Pro with environment modules.
This avoids boilerplate code for loading modules and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""
# pylint: disable=duplicate-code
import os
import sys
import signal
import logging
import base64
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
        print("\r" + " " * 100 + "\r", end="", flush=True)
        print("🛑 Interrupted! Cleaning up job...")
        success = qdel(_CURRENT_JOB_ID, quiet=False)
        if success:
            print("✅ Job cleanup completed")
        else:
            print("⚠️  Job cleanup failed - you may need to manually run: qdel", _CURRENT_JOB_ID)
    print("👋 Goodbye!")
    sys.exit(130)  # Standard exit code for SIGINT

def _get_default_template():
    """Get the default PBS template file path."""
    # Try pkg_resources first
    try:
        template_path = pkg_resources.resource_filename(__name__, 'jobscripts/qmod.pbs')
        if os.path.exists(template_path):
            return template_path
    except (ImportError, AttributeError, FileNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / 'jobscripts' / 'qmod.pbs'
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate qmod.pbs template file. "
        f"Looked in: {current_dir / 'jobscripts' / 'qmod.pbs'}"
    )

@qxub.command()
@click.argument("cmd", nargs=-1)
@click.option("--mod",
              multiple=True,
              help="Environment module to load (can be used multiple times)")
@click.option("--template",
              default=_get_default_template(),
              help="Jobscript template (optional - for further customization)")
@click.option("--pre",
              help="Command to run before the main command (use quotes for commands with options)")
@click.option("--post",
              help="Command to run after the main command (only if main command "
                   "succeeds, use quotes for commands with options)")
@click.pass_context
def module(ctx, cmd, mod, template, pre, post):  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    """
    Constructs and submits a qsub job that will execute the given command
    with the specified environment modules loaded

    Example:
        qxub --resources mem=50MB module --mod bcftools --mod samtools -- bcftools --version

    Here the command "bcftools --version" will run with the bcftools and samtools
    modules loaded with 50MB of RAM. See "qxub --help" for other resource options.

    Note the "--" before the command. This is required if command itself also
    contains double-dashes, as in this case ("--version").
    """
    # Validate that at least one module is specified
    if not mod:
        raise click.ClickException("At least one module must be specified with --mod")

    # Get values from context
    ctx_obj = ctx.obj
    out = Path(ctx_obj['out'])
    err = Path(ctx_obj['err'])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Modules: %s", mod)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", cmd)
    logging.debug("Pre-commands: %s", pre)
    logging.debug("Post-commands: %s", post)

    # Construct qsub command
    cmd_str = " ".join(cmd)
    modules_str = " ".join(mod)  # Space-separated for the template
    # Base64 encode the command to avoid escaping issues
    cmd_b64 = base64.b64encode(cmd_str.encode('utf-8')).decode('ascii')
    submission_vars = (f'modules="{modules_str}",cmd_b64="{cmd_b64}",'
                       f'cwd={ctx_obj["execdir"]},out={out},err={err},'
                       f'quiet={str(ctx_obj["quiet"]).lower()}')
    if pre:
        pre_b64 = base64.b64encode(pre.encode('utf-8')).decode('ascii')
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_b64 = base64.b64encode(post.encode('utf-8')).decode('ascii')
        submission_vars += f',post_cmd_b64="{post_b64}"'

    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'

    # Execute the command
    logging.info("Submission command: %s", submission_command)
    if ctx_obj['dry']:
        print(f"Dry run - would execute: {submission_command}")
        logging.info("Dry run. Exiting")
        sys.exit(0)

    # Register signal handler for Ctrl+C cleanup
    signal.signal(signal.SIGINT, _signal_handler)

    logging.info("Submitting job")
    job_id = qsub(submission_command, quiet=ctx_obj['quiet'])

    # Track job ID globally for signal handler
    global _CURRENT_JOB_ID  # pylint: disable=global-statement
    _CURRENT_JOB_ID = job_id

    # Show success message with job ID
    if not ctx_obj['quiet']:
        success_msg = f"🚀 Job submitted successfully! Job ID: {job_id}"
        print_status(success_msg, final=False)

    logging.info("Your job has been successfully submitted")
    # Exit if in quiet mode
    if ctx_obj['quiet']:
        logging.info("Exiting in quiet mode")
        sys.exit(0)

    # Start concurrent monitoring of job and log files
    # Stream log files to STDOUT/STDERR as appropriate
    out.parent.mkdir(parents=True, exist_ok=True)
    err.parent.mkdir(parents=True, exist_ok=True)
    out.touch()
    err.touch()

    # Pass success message to monitor_and_tail for spinner display
    if not ctx_obj['quiet']:
        success_message = f"🚀 Job submitted successfully! Job ID: {job_id}"
    else:
        success_message = None

    try:
        monitor_and_tail(job_id, out, err, quiet=ctx_obj['quiet'], success_msg=success_message)
    finally:
        # Clear job ID when monitoring completes (successfully or via interrupt)
        _CURRENT_JOB_ID = None
