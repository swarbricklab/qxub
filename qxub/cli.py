"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""
import os
from datetime import datetime
import logging
from pathlib import Path
import click
from .config import setup_logging

@click.group()
@click.option("--execdir",
              default=os.getcwd(),
              help="Execution directory (default: current directory)")
@click.option("--out",
              default=Path(
                  os.getenv("TMPDIR"), "qt", datetime.now().strftime("%Y%m%d_%H%M%S"),"out"),
              help="STDOUT log file (default: $TMPDIR/qt/timestamp/out)")
@click.option("--err",
              default=Path(
                  os.getenv("TMPDIR"), "qt", datetime.now().strftime("%Y%m%d_%H%M%S"),"err"),
              help="STDERR log file (default: $TMPDIR/qt/timestamp/err)")
@click.option("--joblog", help="PBS Pro job log")
@click.option("--dry", "--dry-run", is_flag=True, default=False,
              help="Generate job submission command but don't submit")
@click.option("--quiet", is_flag=True, default=False,
              help="Display no output")
@click.option("-l", "--resources", multiple=True, help="Job resource")
@click.option("-q", "--queue", default="normal", help="Job queue (default: normal)")
@click.option("-N", "--name", default="qt", help="Job name (default: qt)")
@click.option("-P", "--project",
              default=os.getenv("PROJECT"),
              help="PBS project code (default: $PROJECT)")
@click.option("-v", "--verbose", count=True,
              help="Increase verbosity (use -v, -vv, -vvv for more detail)")
@click.pass_context
def qxub(ctx, execdir, verbose, **params):
    """
    Parent command 'qxub' with options common to all subcommands.
    PBS Pro options (-l,-q,-N,-q) are passed on as parameters to 'qsub'
    but can be specified in either short form (-l) or long form (--resources)

    Run 'qxub {subcommand} --help' for details on each subcommand.

    For example: qxub conda --help
    """
    setup_logging(verbose)
    logging.debug("Execution directory: %s", execdir)
    for key, value in params.items():
        logging.debug("qsub option: %s = %s", key, value)
    ctx.ensure_object(dict)
    # Construct the qsub options
    joblog = params['joblog'] or f"{params['name']}.log"
    options = "-N {name} -q {queue} -P {project} ".format_map(params)
    options += " ".join([f"-l {resource}" for resource in params['resources']])
    options += f" -o {joblog}"
    logging.info("Options: %s", options)
    # Load qsub options into context
    ctx.obj['execdir'] = execdir
    ctx.obj['options'] = options
    ctx.obj['name'] = params['name']
    ctx.obj['out'] = params['out']
    ctx.obj['err'] = params['err']
    ctx.obj['dry'] = params['dry']
    ctx.obj['quiet'] = params['quiet']
