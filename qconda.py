import click
import pkg_resources
import os

@click.command()
@click.argument("cmd", nargs=-1)
@click.option("--env", required=True, help="Conda environment to use")
@click.option("-l","--resources", multiple=True, help="Job resource")
@click.option("-q","--queue", default="normal", help="Job queue")
@click.option("-N","--name", default="qconda", help="Job name")
@click.option("-P","--project", default="a56", help="Project")
@click.option("-v","--verbose", is_flag=True, default=False, help="Verbose output")
def conda_submit(cmd, env, resources, queue, name, project, verbose):
    """Runs a command/script in a given environment with the specified PBS options"""
    if verbose:# Display the received inputs
        click.echo(f"Conda environment: {env}")
        click.echo(f"Resources: {resources}")
        click.echo(f"Queue: {queue}")
        click.echo(f"Name: {name}")
        click.echo(f"Project: {project}")
        click.echo(f'Command: {cmd}')
    
    # Get jobscript path
    jobscript = pkg_resources.resource_filename(__name__, 'jobscripts/qconda.pbs')
    # Get current work directory
    cwd=os.getcwd()

    # Construct the submission command
    full_options=f"-N {name} " + f"-q {queue} " + f"-P {project} "
    for resource in resources:
        full_options+=f"-l {resource}" 
    if verbose:
        click.echo(f"Options: {full_options}")
    full_cmd=" ".join(cmd)
    submission_command = f'qsub -v env={env},cmd="{full_cmd}",cwd={cwd} {jobscript} {full_options}'
    
    # Execute the command
    click.echo(f'Submitting: {submission_command}')
    #os.system(submission_command)

if __name__ == '__main__':
    conda_submit()

