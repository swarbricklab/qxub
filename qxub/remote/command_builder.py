"""
Command serialization for remote execution.

This module converts parsed qxub execution context and options back into
CLI arguments suitable for remote execution via SSH.
"""

import shlex
from typing import Any, Dict, List, Tuple

from ..execution_context import ExecutionContext


def build_remote_command(
    execution_context: ExecutionContext,
    options: Dict[str, Any],
    command: List[str],
) -> str:
    """
    Build qxub command string for remote execution.

    Serializes execution context, PBS options, and command into an equivalent
    qxub CLI invocation that can be executed remotely via SSH.

    Args:
        execution_context: Parsed execution context (conda env, modules, container)
        options: Dictionary of parsed CLI options
        command: Command arguments to execute

    Returns:
        Complete qxub command string with proper quoting

    Examples:
        >>> ctx = ExecutionContext('conda', 'pytorch', 'conda')
        >>> opts = {'queue': 'auto', 'mem': '8GB', 'cpus': 4}
        >>> cmd = ['python', 'train.py', '--epochs', '100']
        >>> build_remote_command(ctx, opts, cmd)
        'qxub exec --env pytorch --queue auto --mem 8GB --cpus 4 -- python train.py --epochs 100'
    """
    parts = ["qxub", "exec"]

    # Serialize execution context
    if execution_context.context_type == "conda":
        parts.extend(["--env", shlex.quote(str(execution_context.context_value))])
    elif execution_context.context_type == "module":
        # Module list can be string or list
        if isinstance(execution_context.context_value, list):
            for module in execution_context.context_value:
                parts.extend(["--mod", shlex.quote(str(module))])
        else:
            parts.extend(["--mods", shlex.quote(str(execution_context.context_value))])
    elif execution_context.context_type == "singularity":
        parts.extend(["--sif", shlex.quote(str(execution_context.context_value))])
        # Add bind mounts if specified
        if options.get("bind"):
            parts.extend(["--bind", shlex.quote(str(options["bind"]))])
    elif execution_context.context_type == "default":
        parts.append("--default")

    # Serialize PBS resource options
    _add_pbs_options(parts, options)

    # Serialize workflow-friendly resource options
    _add_workflow_options(parts, options)

    # Serialize job configuration options
    _add_job_options(parts, options)

    # Add command separator and command
    if command:
        parts.append("--")
        parts.extend(shlex.quote(str(arg)) for arg in command)

    return " ".join(parts)


def _add_pbs_options(parts: List[str], options: Dict[str, Any]) -> None:
    """Add PBS-specific options to command parts."""
    # Queue selection
    if options.get("queue"):
        parts.extend(["--queue", shlex.quote(str(options["queue"]))])

    # Resources (PBS format like 'walltime=1:00:00,mem=4GB')
    if options.get("resources"):
        resources = options["resources"]
        if isinstance(resources, (list, tuple)):
            for resource in resources:
                parts.extend(["-l", shlex.quote(str(resource))])
        else:
            parts.extend(["-l", shlex.quote(str(resources))])

    # Array jobs
    if options.get("array"):
        parts.extend(["--array", shlex.quote(str(options["array"]))])


def _add_workflow_options(parts: List[str], options: Dict[str, Any]) -> None:
    """Add workflow-friendly resource options to command parts."""
    # Memory
    if options.get("mem"):
        parts.extend(["--mem", shlex.quote(str(options["mem"]))])

    # Runtime/walltime
    if options.get("runtime"):
        parts.extend(["--runtime", shlex.quote(str(options["runtime"]))])

    # CPUs
    if options.get("cpus"):
        parts.extend(["--cpus", str(options["cpus"])])

    # Disk/jobfs
    if options.get("disk"):
        parts.extend(["--disk", shlex.quote(str(options["disk"]))])

    # Storage volumes
    if options.get("volumes"):
        parts.extend(["--volumes", shlex.quote(str(options["volumes"]))])


def _add_job_options(parts: List[str], options: Dict[str, Any]) -> None:
    """Add job configuration options to command parts."""
    # Job name
    if options.get("name"):
        parts.extend(["-N", shlex.quote(str(options["name"]))])

    # Project code
    if options.get("project"):
        parts.extend(["-P", shlex.quote(str(options["project"]))])

    # Output files
    if options.get("out"):
        parts.extend(["--out", shlex.quote(str(options["out"]))])
    if options.get("err"):
        parts.extend(["--err", shlex.quote(str(options["err"]))])

    # Execution directory
    if options.get("execdir"):
        parts.extend(["--execdir", shlex.quote(str(options["execdir"]))])

    # Email options
    if options.get("email"):
        parts.extend(["--email", shlex.quote(str(options["email"]))])
    if options.get("email_opts"):
        parts.extend(["--email-opts", shlex.quote(str(options["email_opts"]))])

    # Template customization
    if options.get("template"):
        parts.extend(["--template", shlex.quote(str(options["template"]))])

    # Pre/post commands
    if options.get("pre"):
        parts.extend(["--pre", shlex.quote(str(options["pre"]))])
    if options.get("post"):
        parts.extend(["--post", shlex.quote(str(options["post"]))])

    # Execution flags
    if options.get("dry"):
        parts.append("--dry")
    if options.get("quiet"):
        parts.append("--quiet")
    if options.get("terse"):
        parts.append("--terse")

    # Verbosity (count-based, add -v for each level)
    verbose_level = options.get("verbose", 0)
    if verbose_level > 0:
        parts.append("-" + "v" * verbose_level)


def parse_options_from_ctx(ctx_obj: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
    """
    Extract options dict and qsub_options from context object.

    The context object (ctx.obj) contains processed parameters and qsub options.
    This extracts the relevant parts for command serialization.

    Args:
        ctx_obj: Click context object dictionary

    Returns:
        Tuple of (options dict, qsub_options list)
    """
    options = {}
    qsub_options = ctx_obj.get("options", [])

    # Extract standard options from ctx.obj
    option_keys = [
        "queue",
        "name",
        "project",
        "out",
        "err",
        "execdir",
        "email",
        "email_opts",
        "array",
        "mem",
        "runtime",
        "cpus",
        "disk",
        "volumes",
        "template",
        "pre",
        "post",
        "bind",
        "dry",
        "quiet",
        "terse",
        "verbose",
    ]

    for key in option_keys:
        if key in ctx_obj:
            options[key] = ctx_obj[key]

    # Parse PBS resources from qsub_options list
    # qsub_options format: ['-l', 'walltime=1:00:00', '-l', 'mem=4GB', ...]
    resources = []
    i = 0
    while i < len(qsub_options):
        if qsub_options[i] == "-l" and i + 1 < len(qsub_options):
            resources.append(qsub_options[i + 1])
            i += 2
        else:
            i += 1

    if resources:
        options["resources"] = resources

    return options, qsub_options
