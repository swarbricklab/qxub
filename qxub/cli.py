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
from .config_cli import config_cli
from .alias_cli import alias_cli, alias_test_cli
from .config import setup_logging
from .config_manager import config_manager


def _get_default_output_dir():
    """Get appropriate output directory that works across login and compute nodes."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Prefer shared scratch space over TMPDIR to avoid cross-node access issues
    project = os.getenv("PROJECT", "a56")
    user = os.getenv("USER", "unknown")

    # Try shared scratch first
    scratch_path = Path("/scratch", project, user, "qt", timestamp)
    if scratch_path.parent.parent.exists():  # /scratch/PROJECT/USER exists
        return scratch_path

    # Fallback to current directory
    return Path(os.getcwd(), "qt", timestamp)


def _get_config_default(key: str, fallback=None):
    """Get default value from config or fallback."""
    defaults = config_manager.get_defaults()
    return defaults.get(key, fallback)


def _sanitize_job_name(name: str) -> str:
    """Sanitize job name for PBS compliance.

    PBS job names cannot contain certain characters like /, :, @, etc.
    Replace problematic characters with safe alternatives.
    """
    if not name:
        return name

    # Replace problematic characters with safe alternatives
    sanitized = name.replace("/", "_")  # Paths to underscores
    sanitized = sanitized.replace(":", "_")  # Colons to underscores
    sanitized = sanitized.replace(" ", "_")  # Spaces to underscores
    sanitized = sanitized.replace("@", "_")  # At symbols to underscores

    # Remove any remaining non-alphanumeric characters except hyphens and underscores
    sanitized = "".join(c for c in sanitized if c.isalnum() or c in "-_")

    # Ensure it starts with a letter or number (not special character)
    if sanitized and not sanitized[0].isalnum():
        sanitized = "job_" + sanitized

    # Limit length (PBS has limits on job name length)
    if len(sanitized) > 50:
        sanitized = sanitized[:47] + "..."

    return sanitized or "job"  # Fallback if sanitization results in empty string


def _get_config_default_callable(key: str, fallback_callable):
    """Get default value from config or call fallback function."""

    def _default():
        defaults = config_manager.get_defaults()
        value = defaults.get(key)
        if value is not None:
            # Resolve templates if it's a string with placeholders
            if isinstance(value, str) and "{" in value:
                template_vars = config_manager.get_template_variables()
                return config_manager.resolve_templates(value, template_vars)
            return value
        return fallback_callable()

    return _default


@click.group()
@click.option(
    "--execdir",
    default=os.getcwd(),
    help="Execution directory (default: current directory)",
)
@click.option(
    "--out",
    help="STDOUT log file (default: configured or "
    "/scratch/$PROJECT/$USER/qt/timestamp/out)",
)
@click.option(
    "--err",
    help="STDERR log file (default: configured or "
    "/scratch/$PROJECT/$USER/qt/timestamp/err)",
)
@click.option("--joblog", help="PBS Pro job log (default: configured or {name}.log)")
@click.option(
    "--dry",
    "--dry-run",
    is_flag=True,
    default=False,
    help="Generate job submission command but don't submit",
)
@click.option("--quiet", is_flag=True, default=False, help="Display no output")
@click.option(
    "-l", "--resources", multiple=True, help="Job resource (default: configured)"
)
@click.option("-q", "--queue", help="Job queue (default: configured or normal)")
@click.option("-N", "--name", help="Job name (default: configured or qt)")
@click.option(
    "-P", "--project", help="PBS project code (default: configured or $PROJECT)"
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (use -v, -vv, -vvv for more detail)",
)
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

    # Apply config defaults for any unset parameters
    defaults = config_manager.get_defaults()

    # Apply defaults if not explicitly set by CLI
    if params["out"] is None:
        params["out"] = defaults.get("out", _get_default_output_dir() / "out")
    if params["err"] is None:
        params["err"] = defaults.get("err", _get_default_output_dir() / "err")
    if params["joblog"] is None:
        params["joblog"] = defaults.get("joblog")
    if params["queue"] is None:
        params["queue"] = defaults.get("queue", "normal")
    if params["name"] is None:
        params["name"] = defaults.get("name", "qt")
    if params["project"] is None:
        params["project"] = defaults.get("project", os.getenv("PROJECT"))
    if not params["resources"]:  # Empty tuple means no resources provided
        params["resources"] = defaults.get("resources", [])

    for key, value in params.items():
        logging.debug("qsub option: %s = %s", key, value)
    ctx.ensure_object(dict)

    # Resolve template variables for any config-provided values
    template_vars = config_manager.get_template_variables(
        name=params["name"], project=params["project"], queue=params["queue"]
    )

    # Resolve template strings
    resolved_params = {}
    for key, value in params.items():
        if isinstance(value, str) and "{" in value:
            resolved_params[key] = config_manager.resolve_templates(
                value, template_vars
            )
        else:
            resolved_params[key] = value
    params.update(resolved_params)

    # Sanitize job name for PBS compliance
    if params["name"]:
        params["name"] = _sanitize_job_name(params["name"])

    # Handle joblog default
    joblog = params["joblog"] or f"{params['name']}.log"
    if isinstance(joblog, str) and "{" in joblog:
        joblog = config_manager.resolve_templates(joblog, template_vars)

    # Construct the qsub options
    options = "-N {name} -q {queue} -P {project} ".format_map(params)
    if params.get("resources"):
        options += " ".join([f"-l {resource}" for resource in params["resources"]])
    options += f" -o {joblog}"
    logging.info("Options: %s", options)

    # Load qsub options into context
    ctx.obj["execdir"] = execdir
    ctx.obj["options"] = options
    ctx.obj["name"] = params["name"]
    ctx.obj["out"] = params["out"]
    ctx.obj["err"] = params["err"]
    ctx.obj["dry"] = params["dry"]
    ctx.obj["quiet"] = params["quiet"]


qxub.add_command(config_cli)
qxub.add_command(alias_cli)
qxub.add_command(alias_test_cli)
