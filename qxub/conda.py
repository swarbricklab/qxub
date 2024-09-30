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
from .scheduler import qsub,monitor_and_tail

@qxub.command()
@click.argument("cmd", nargs=-1)
@click.option("--env",
              default=os.getenv('CONDA_DEFAULT_ENV'),
              help="Conda environment to use (default: active environment)")
@click.option("--template",
              default=pkg_resources.resource_filename(__name__, 'jobscripts/qconda.pbs'),
              help="Jobscript template (optional - for further customization)")
@click.pass_context
def conda(ctx, cmd, env, template):
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
    execdir = ctx.obj['execdir']
    options = ctx.obj['options']
    out = Path(ctx.obj['out'])
    err = Path(ctx.obj['err'])
    dry = ctx.obj['dry']
    quiet = ctx.obj['quiet']
    # Log parameters and context
    for key, value in ctx.obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Conda environment: %s", env)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", cmd)
    # Construct qsub command
    full_cmd = " ".join(cmd)
    submission_command = f'qsub -v env={env},cmd="{full_cmd}",cwd={execdir},'
    submission_command+= f'out={out},err={err} {options} {template}'
    # Execute the command
    logging.info("Submission command: %s", submission_command)
    if dry:
        logging.info("Dry run. Exiting")
        sys.exit(0)
    logging.info("Submitting job")
    job_id=qsub(submission_command)
    logging.info("Your job has been successfully submitted")
    # Exit if in quiet mode
    if quiet:
        logging.info("Exiting in quiet mode")
        sys.exit(0)

    # Start concurrent monitoring of job and log files
    # Stream log files to STDOUT/STDERR as appropriate
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(err).parent.mkdir(parents=True, exist_ok=True)
    out.touch()
    err.touch()
    monitor_and_tail(job_id, out, err)
