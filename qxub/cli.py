"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

import click

from .alias_cli import alias_cli
from .config_cli import config_cli
from .exec_cli import exec_cli
from .history_cli import history
from .monitor_cli import monitor_cli
from .platform_cli import estimate_cmd, platform_cli, select_queue_cmd, validate_cmd
from .resources_cli import resources
from .run_cli import run_cli


@click.group(invoke_without_command=True)
@click.option(
    "--version",
    is_flag=True,
    help="Show version and exit",
)
@click.pass_context
def qxub(ctx, version):
    """
    qxub - PBS job submission tools

    For job execution, use the exec subcommand:
        qxub exec --env myenv -- python script.py
        qxub exec --mod python3 --mod gcc -- make
        qxub exec --mods python3,gcc -- python script.py
        qxub exec --sif container.sif -- python script.py
        qxub exec --default -- echo "hello world"

    For complex commands with variables:
        qxub exec --env myenv --cmd "python script.py --input ${HOME}/data.txt"
        qxub exec --env myenv --cmd "echo 'Job ${{PBS_JOBID}} on node ${{HOSTNAME}}'"

    Use --queue auto for intelligent queue selection:
        qxub exec --queue auto -l mem=500GB --env myenv -- python big_job.py
    """
    # Handle version flag first
    if version:
        from . import __version__

        click.echo(f"qxub {__version__}")
        ctx.exit()

    # If no subcommand and no version flag, show help
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        ctx.exit(0)


# Register all subcommands
qxub.add_command(exec_cli)
qxub.add_command(config_cli)
qxub.add_command(alias_cli)
qxub.add_command(history)
qxub.add_command(monitor_cli)
qxub.add_command(resources)
qxub.add_command(platform_cli)
qxub.add_command(select_queue_cmd)
qxub.add_command(validate_cmd)
qxub.add_command(estimate_cmd)
qxub.add_command(run_cli)
