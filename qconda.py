import click
import pkg_resources
import os

@click.group()
@click.option("-l", "--resources", multiple=True, help="Job resource")
@click.option("-q", "--queue", default="normal", help="Job queue")
@click.option("-N", "--name", default="qconda", help="Job name")
@click.option("-P", "--project", default="a56", help="Project")
@click.option("-v", "--verbose", is_flag=True, default=False, help="Verbose output")
@click.pass_context
def qt(ctx, resources, queue, name, project, verbose):
    """Parent command 'qt' with shared options."""
    ctx.ensure_object(dict)
    # Construct the submission command
    options = f"-N {name} -q {queue} -P {project} "
    for resource in resources:
        options += f"-l {resource} "
    if verbose:
        click.echo(f"Options: {options}")
    ctx.obj['options'] = options
    ctx.obj['verbose'] = verbose

@qt.command()
@click.argument("cmd", nargs=-1)
@click.option("--env", required=True, help="Conda environment to use")
@click.pass_context
def conda(ctx, cmd, env):
    """Runs a command in a given environment."""
    options = ctx.obj['options']
    verbose = ctx.obj['verbose']

    if verbose:
        click.echo(f"Conda environment: {env}")
        click.echo(f'Command: {cmd}')

    # Get jobscript path
    jobscript = pkg_resources.resource_filename(__name__, 'jobscripts/qconda.pbs')
    # Get current working directory
    cwd = os.getcwd()
    
    full_cmd = " ".join(cmd)
    submission_command = f'qsub -v env={env},cmd="{full_cmd}",cwd={cwd} {options} {jobscript}'

    # Execute the command
    click.echo(f'Submitting: {submission_command}')
    os.system(submission_command)