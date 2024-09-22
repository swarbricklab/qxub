"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""
import os
import logging
import click
import pkg_resources

def setup_logging(verbose):
    """
    Configures the logging level based on the verbosity provided by the user.

    Args:
        verbose (int): The number of '-v' flags used. 
                       - 0: ERROR level (default)
                       - 1: WARNING level
                       - 2: INFO level
                       - 3 or more: DEBUG level

    This function adjusts the logging output to provide more detailed information
    as verbosity increases, allowing users to control the granularity of log messages.
    """
    if verbose == 1:
        logging.basicConfig(level=logging.WARNING)
    elif verbose == 2:
        logging.basicConfig(level=logging.INFO)
    elif verbose >= 3:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.ERROR)

@click.group()
@click.option("-l", "--resources", multiple=True, help="Job resource")
@click.option("-q", "--queue", default="normal", help="Job queue")
@click.option("-N", "--name", default="qconda", help="Job name")
@click.option("-P", "--project", default="a56", help="Project")
@click.option("-v", "--verbose", count=True,
              help="Increase verbosity (use -v, -vv, -vvv for more detail)")
@click.pass_context
def qt(ctx, verbose, **params):
    """Parent command 'qt' with shared options."""
    setup_logging(verbose)
    logging.debug("qsub options: %s", params)
    ctx.ensure_object(dict)
    # Construct the qsub options
    options = "-N {name} -q {queue} -P {project} ".format_map(params)
    options += " ".join([f"-l {resource}" for resource in params['resources']])
    logging.info("Options: %s:, options")
    # Load qsub options into context
    ctx.obj['options'] = options

@qt.command()
@click.argument("cmd", nargs=-1)
@click.option("--env", required=True, help="Conda environment to use")
@click.option("--execdir", default=os.getcwd(), help="Execution directory")
@click.pass_context
def conda(ctx, cmd, env, execdir):
    """Runs a command in a given environment."""
    options = ctx.obj['options']

    logging.info("Conda environment: %s", env)
    logging.info("Execution directory: %s", execdir)
    logging.info("Command: %s", cmd)

    # Get jobscript path
    jobscript = pkg_resources.resource_filename(__name__, 'jobscripts/qconda.pbs')

    full_cmd = " ".join(cmd)
    submission_command = f'qsub -v env={env},cmd="{full_cmd}",cwd={execdir} {options} {jobscript}'

    # Execute the command
    logging.info("Submitting: %s", submission_command)
    os.system(submission_command)
