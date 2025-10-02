"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""
import os
import sys
import logging
from pathlib import Path
import click
import pkg_resources
from .cli import qxub
from .scheduler import qsub, monitor_and_tail, print_status

def _get_default_template():
    """Get the default PBS template file path."""
    # Try pkg_resources first
    try:
        template_path = pkg_resources.resource_filename(__name__, 'jobscripts/qconda.pbs')
        if os.path.exists(template_path):
            return template_path
    except (ImportError, AttributeError, FileNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / 'jobscripts' / 'qconda.pbs'
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate qconda.pbs template file. "
        f"Looked in: {current_dir / 'jobscripts' / 'qconda.pbs'}"
    )

@qxub.command()
@click.argument("cmd", nargs=-1)
@click.option("--env",
              default=os.getenv('CONDA_DEFAULT_ENV'),
              help="Conda environment to use (default: active environment)")
@click.option("--template",
              default=_get_default_template(),
              help="Jobscript template (optional - for further customization)")
@click.option("--pre",
              help="Command to run before the main command")
@click.option("--post",
              help="Command to run after the main command (only if main command succeeds)")
@click.pass_context
def conda(ctx, cmd, env, template, pre, post):  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
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
    # Get values from context
    ctx_obj = ctx.obj
    out = Path(ctx_obj['out'])
    err = Path(ctx_obj['err'])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Conda environment: %s", env)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", cmd)
    logging.debug("Pre-command: %s", pre)
    logging.debug("Post-command: %s", post)

    # Construct qsub command
    cmd_str = " ".join(cmd)
    submission_vars = f'env={env},cmd="{cmd_str}",cwd={ctx_obj["execdir"]},out={out},err={err}'
    if pre:
        submission_vars += f',pre_cmd="{pre}"'
    if post:
        submission_vars += f',post_cmd="{post}"'

    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'

    # Execute the command
    logging.info("Submission command: %s", submission_command)
    if ctx_obj['dry']:
        logging.info("Dry run. Exiting")
        sys.exit(0)

    logging.info("Submitting job")
    job_id = qsub(submission_command, quiet=ctx_obj['quiet'])

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
    monitor_and_tail(job_id, out, err, quiet=ctx_obj['quiet'], success_msg=success_message)
