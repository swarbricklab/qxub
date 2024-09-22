"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""
import os
import logging
import click
import pkg_resources

def setup_logging(verbosity):
    """
    Configures the logging level based on the verbosity provided by the user.

    Args:
        verbosity (int): The number of '-v' flags used. 
                       - 0: ERROR level (default)
                       - 1: WARNING level
                       - 2: INFO level
                       - 3 or more: DEBUG level

    This function adjusts the logging output to provide more detailed information
    as verbosity increases, allowing users to control the granularity of log messages.
    """
    if verbosity == 1:
        logging.basicConfig(level=logging.WARNING)
    elif verbosity == 2:
        logging.basicConfig(level=logging.INFO)
    elif verbosity >= 3:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.ERROR)

@click.group()
@click.option("-l", "--resources", multiple=True, help="Job resource")
@click.option("-q", "--queue", default="normal", help="Job queue (default: normal)")
@click.option("-N", "--name", default="qt", help="Job name (default: qt)")
@click.option("-P", "--project",
              default=os.getenv("PROJECT"),
              help="PBS project code (default: $PROJECT)")
@click.option("-v", "--verbose", count=True,
              help="Increase verbosity (use -v, -vv, -vvv for more detail)")
@click.pass_context
def qt(ctx, verbose, **params):
    """
    Parent command 'qt' with options common to all subcommands.
    PBS Pro options (-l,-q,-N,-q) are passed on as parameters to 'qsub'
    but can be specified in either short form (-l) or long form (--resources)

    Run 'qt {subcommand} --help' for details on each subcommand.

    For example: qt conda --help
    """
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
@click.option("--env",
              default=os.getenv('CONDA_DEFAULT_ENV'),
              help="Conda environment to use (default: active environment)")
@click.option("--execdir",
              default=os.getcwd(),
              help="Execution directory (default: current directory)")
@click.option("--template",
              default=pkg_resources.resource_filename(__name__, 'jobscripts/qconda.pbs'),
              help="Jobscript template (optional - for further customization)")
@click.pass_context
def conda(ctx, cmd, env, execdir,template):
    """
    Constructs and submits a qsub job that will execute the given command
    in the specified conda environment and work directory

    Example:
        qt --resources mem=50MB conda --env myenv -- python --version

    Here the command "python --version" will run in the "myenv" conda environment
    with 50MB of RAM. See "qt --help" for other resource options.

    Note the "--" before the command. This is required if command itself also
    contains double-dashes, as in this case ("--version").
    """
    options = ctx.obj['options']
    # Log parameters and context
    logging.debug("Context: %s", ctx)
    logging.debug("Conda environment: %s", env)
    logging.debug("Execution directory: %s", execdir)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", cmd)
    # Construct qsub command
    full_cmd = " ".join(cmd)
    submission_command = f'qsub -v env={env},cmd="{full_cmd}",cwd={execdir} {options} {template}'
    # Execute the command
    logging.info("Submitting: %s", submission_command)
    os.system(submission_command)
