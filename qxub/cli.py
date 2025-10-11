"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

import base64
import difflib
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

import click

from .alias_cli import alias_cli
from .config import setup_logging
from .config_cli import config_cli
from .config_manager import config_manager
from .history_cli import history
from .history_manager import history_manager
from .platform_cli import estimate_cmd, platform_cli, select_queue_cmd, validate_cmd
from .resource_tracker import resource_tracker
from .resources_cli import resources
from .scheduler import get_job_resource_data, monitor_and_tail, print_status, qdel, qsub


class QxubGroup(click.Group):
    """Custom Click group with enhanced error handling for unknown options."""

    def get_command(self, ctx, cmd_name):
        """Override get_command to handle execution context."""
        # Check if we have execution context
        execution_options = ["env", "mod", "mods", "sif"]
        has_execution_context = any(ctx.params.get(opt) for opt in execution_options)

        if has_execution_context:
            return None

        # Check if we have protected args (commands after --) without execution context
        # This indicates default execution should be used
        if hasattr(ctx, "protected_args") and ctx.protected_args:
            return None

        return super().get_command(ctx, cmd_name)

    def invoke(self, ctx):
        """Override invoke to handle execution contexts."""
        # Check if we have execution context
        execution_options = ["env", "mod", "mods", "sif"]
        has_execution_context = any(ctx.params.get(opt) for opt in execution_options)

        if (
            has_execution_context or hasattr(ctx, "protected_args")
        ) and ctx.protected_args:
            # We have execution context and protected args (the command)
            # Combine args for the command and invoke the main function directly
            combined_args = list(ctx.protected_args) + ctx.args

            # Temporarily set ctx.args for the main function
            original_args = ctx.args
            ctx.args = combined_args
            try:
                # Extract the parameters that the qxub function expects
                qxub_func = (
                    ctx.command.callback.__wrapped__
                )  # Get the unwrapped function
                return qxub_func(
                    ctx,
                    ctx.params["remote"],
                    ctx.params["config"],
                    ctx.params["execdir"],
                    ctx.params["verbose"],
                    ctx.params["version"],
                    ctx.params["env"],
                    ctx.params["mod"],
                    ctx.params["mods"],
                    ctx.params["sif"],
                    ctx.params["bind"],
                    ctx.params["template"],
                    ctx.params["pre"],
                    ctx.params["post"],
                    **{
                        k: v
                        for k, v in ctx.params.items()
                        if k
                        not in [
                            "remote",
                            "config",
                            "execdir",
                            "verbose",
                            "version",
                            "env",
                            "mod",
                            "mods",
                            "sif",
                            "bind",
                            "template",
                            "pre",
                            "post",
                        ]
                    },
                )
            finally:
                ctx.args = original_args

        return super().invoke(ctx)

    def parse_args(self, ctx, args):
        """Override parse_args to provide better error messages for unknown options."""
        try:
            return super().parse_args(ctx, args)
        except click.NoSuchOption as e:
            self.handle_unknown_option_error(ctx, e, args)

    def handle_unknown_option_error(self, ctx, error, args):
        """Provide helpful suggestions for unknown options."""
        unknown_option = error.option_name

        # Get all valid options for the current command
        all_options = []
        for param in ctx.command.params:
            if isinstance(param, click.Option):
                all_options.extend(param.opts)

        # Also get options from subcommands if we can identify the subcommand
        subcommand_name = None
        subcommand_options = []
        for i, arg in enumerate(args):
            if arg in self.commands:
                subcommand_name = arg
                subcommand = self.commands[arg]
                for param in subcommand.params:
                    if isinstance(param, click.Option):
                        subcommand_options.extend(param.opts)
                break

        # Find close matches using difflib
        close_matches = difflib.get_close_matches(
            unknown_option, all_options + subcommand_options, n=3, cutoff=0.6
        )

        # Build helpful error message
        click.echo(f"Error: No such option: {unknown_option}", err=True)

        if close_matches:
            click.echo("\n💡 Did you mean:", err=True)
            for match in close_matches:
                if match in all_options:
                    click.echo(f"   {match}  (qxub option)", err=True)
                else:
                    click.echo(f"   {match}  ({subcommand_name} option)", err=True)

        # Provide guidance about command vs option separation
        if unknown_option.startswith("-"):
            click.echo(
                "\n📖 Common issue: Mixing command options with qxub options", err=True
            )
            click.echo("   qxub options must come BEFORE the subcommand:", err=True)
            click.echo(
                "   ✅ qxub --queue normal conda --env myenv python script.py", err=True
            )
            click.echo(
                "   ❌ qxub conda --queue normal --env myenv python script.py", err=True
            )

            # Check if this looks like a command argument that should use --
            if self._looks_like_command_argument(unknown_option, args):
                click.echo(
                    "\n   If this is part of your command, use '--' to separate:",
                    err=True,
                )
                click.echo(
                    "   ✅ qxub conda --env myenv -- python script.py -c 'print(\"hello\")'",
                    err=True,
                )
                click.echo(
                    "   ❌ qxub conda --env myenv python script.py -c 'print(\"hello\")'",
                    err=True,
                )

        # Show relevant help
        if subcommand_name:
            click.echo(f"\n🔍 For help: qxub {subcommand_name} --help", err=True)
        else:
            click.echo("\n🔍 For help: qxub --help", err=True)

        ctx.exit(2)

    def _looks_like_command_argument(self, option, args):
        """Check if the unknown option looks like it should be part of a command."""
        # Common command-line options that users might accidentally use
        command_like_options = [
            "-c",
            "-i",
            "-o",
            "-f",
            "-d",
            "-r",
            "-t",
            "-s",
            "-n",
            "-p",
            "-h",
            "--input",
            "--output",
            "--file",
            "--config",
            "--help",
            "--version",
            "--debug",
            "--verbose",
            "--quiet",
            "--force",
            "--recursive",
        ]

        # If it matches common command options, likely a command argument
        if option in command_like_options:
            return True

        # If there are non-option arguments after this, probably a command
        option_index = None
        try:
            option_index = args.index(option)
        except ValueError:
            return False

        # Check if there are command-like things after this option
        for i in range(option_index + 1, len(args)):
            arg = args[i]
            if not arg.startswith("-") and "." in arg:  # Looks like a filename
                return True
            if arg in ["python", "bash", "sh", "perl", "ruby", "node", "java"]:
                return True

        return False


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


def _get_conda_template():
    """Get default conda template path."""
    import pkg_resources

    # Try pkg_resources first
    try:
        template_path = pkg_resources.resource_filename(
            __name__, "jobscripts/qconda.pbs"
        )
        if os.path.exists(template_path):
            return template_path
    except (ImportError, AttributeError, FileNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / "jobscripts" / "qconda.pbs"
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate qconda.pbs template file. "
        f"Looked in: {current_dir / 'jobscripts' / 'qconda.pbs'}"
    )


def _get_module_template():
    """Get default module template path."""
    import pkg_resources

    # Try pkg_resources first
    try:
        template_path = pkg_resources.resource_filename(__name__, "jobscripts/qmod.pbs")
        if os.path.exists(template_path):
            return template_path
    except (ImportError, AttributeError, FileNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / "jobscripts" / "qmod.pbs"
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate qmod.pbs template file. "
        f"Looked in: {current_dir / 'jobscripts' / 'qmod.pbs'}"
    )


def _get_singularity_template():
    """Get default singularity template path."""
    import pkg_resources

    # Try pkg_resources first
    try:
        template_path = pkg_resources.resource_filename(
            __name__, "jobscripts/qsing.pbs"
        )
        if os.path.exists(template_path):
            return template_path
    except (ImportError, AttributeError, FileNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / "jobscripts" / "qsing.pbs"
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate qsing.pbs template file. "
        f"Looked in: {current_dir / 'jobscripts' / 'qsing.pbs'}"
    )


@click.group(cls=QxubGroup, invoke_without_command=True)
@click.option(
    "--remote",
    help="Execute on remote system (defined in ~/.config/qxub/config.yaml)",
)
@click.option(
    "--config",
    help="Alternative configuration file (default: ~/.config/qxub/config.yaml)",
)
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
@click.option(
    "-q",
    "--queue",
    help="Job queue (default: configured or normal, use 'auto' for intelligent selection)",
)
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
@click.option(
    "--version",
    is_flag=True,
    help="Show version and exit",
)
# Execution context options (mutually exclusive)
@click.option("--env", "--conda", help="Conda environment for execution")
@click.option("--mod", multiple=True, help="Environment module to load (repeatable)")
@click.option("--mods", "--modules", help="Comma-separated list of environment modules")
@click.option("--sif", "--sing", "--singularity", help="Singularity container image")
# Additional options from subcommands
@click.option("--bind", help="Singularity bind mounts")
@click.option(
    "--template", help="Jobscript template (optional - for further customization)"
)
@click.option("--pre", help="Command to run before the main command")
@click.option("--post", help="Command to run after the main command")
@click.pass_context
def qxub(
    ctx,
    remote,
    config,
    execdir,
    verbose,
    version,
    env,
    mod,
    mods,
    sif,
    bind,
    template,
    pre,
    post,
    **params,
):
    """
    Submit PBS jobs with conda environments, environment modules, or Singularity containers.

    Management commands:
        qxub config --help
        qxub alias --help
        qxub history --help
        qxub resources --help

    For job execution, use the execution options directly:
        qxub --env myenv -- python script.py
        qxub --mod python3 --mod gcc -- make
        qxub --mods python3,gcc -- python script.py
        qxub --sif container.sif -- python script.py

    Use --queue auto for intelligent queue selection:
        qxub --queue auto -l mem=500GB --env myenv -- python big_job.py
    """
    # Handle version flag first
    if version:
        from . import __version__

        click.echo(f"qxub {__version__}")
        ctx.exit()

    # Handle remote execution early
    if remote:
        return _handle_remote_execution(ctx, remote, config, execdir, **params)

    setup_logging(verbose)
    logging.debug("Execution directory: %s", execdir)
    logging.debug("Remaining args: %s", ctx.args)
    logging.debug("Invoked subcommand: %s", ctx.invoked_subcommand)

    # Get remaining arguments (these would be the command to execute)
    command = tuple(ctx.args) if ctx.args else tuple()

    # Consolidate alternative option names
    conda_env = env  # --env and --conda both map to 'env' parameter

    # Handle module options: --mod (multiple) vs --mods/--modules (comma-separated)
    module_list = None
    if mod:
        module_list = list(mod)  # --mod can be used multiple times
    elif mods:  # --mods and --modules both map to 'mods' parameter
        module_list = [m.strip() for m in mods.split(",")]

    container = sif  # --sif, --sing, and --singularity all map to 'sif' parameter

    # Check if any execution context is specified
    execution_contexts = [conda_env, module_list, container]
    has_execution_context = any(execution_contexts)

    # If a subcommand is being invoked, just set up context and let Click handle it
    if ctx.invoked_subcommand is not None:
        logging.debug("Subcommand '%s' will be invoked", ctx.invoked_subcommand)
        # Fall through to context setup and let Click handle the subcommand
    elif has_execution_context:
        # We have execution context - this is direct execution
        # Validate mutual exclusivity
        if sum(bool(x) for x in execution_contexts) > 1:
            raise click.ClickException("Cannot specify multiple execution contexts")

        # Validate that a command is provided
        if not command:
            click.echo(
                "Error: Command is required when execution context is specified.",
                err=True,
            )
            ctx.exit(2)

        # Clear ctx.args so Click doesn't try to parse them as subcommands
        ctx.args = []
    elif command:
        # Command provided but no execution context - will be handled after context setup
        pass
    else:
        # No command, no execution context, no subcommand - show help
        if ctx.invoked_subcommand is None:
            click.echo(ctx.get_help())
            ctx.exit(0)

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

    # Auto-select queue if needed
    if params["queue"] == "auto":
        try:
            from pathlib import Path

            from .platform import PlatformLoader
            from .resource_utils import parse_memory_size, parse_walltime

            # Check for QXUB_PLATFORM_PATHS environment variable
            platform_paths_env = os.environ.get("QXUB_PLATFORM_PATHS")
            if platform_paths_env:
                search_paths = [Path(p.strip()) for p in platform_paths_env.split(":")]
                loader = PlatformLoader(search_paths=search_paths)
            else:
                loader = PlatformLoader()

            platform_names = loader.list_platforms()

            if not platform_names:
                logging.warning(
                    "No platforms available for auto queue selection, using 'normal'"
                )
                params["queue"] = "normal"
            else:
                # Build requirements from resources
                requirements = {}
                if params.get("resources"):
                    for resource in params["resources"]:
                        if "=" in resource:
                            key, value = resource.split("=", 1)
                            if key == "mem":
                                requirements["memory"] = (
                                    value  # Keep as string for condition parsing
                                )
                            elif key == "walltime":
                                requirements["walltime"] = (
                                    value  # Keep as string for condition parsing
                                )
                            elif key == "ncpus":
                                requirements["cpus"] = int(
                                    value
                                )  # Use "cpus" not "ncpus"
                            elif key in ["ngpus", "gpu"]:
                                gpu_count = int(value) if value.isdigit() else 1
                                requirements["gpu_requested"] = (
                                    gpu_count  # Use "gpu_requested" for conditions
                                )
                                requirements["gpus"] = (
                                    gpu_count  # Keep "gpus" for validation
                                )

                # Try to find best queue from any platform
                best_queue = None
                best_cost = float("inf")

                for platform_name in platform_names:
                    platform = loader.get_platform(platform_name)
                    if not platform:
                        continue

                    try:
                        selected_queue = platform.select_queue(requirements)
                        if selected_queue:
                            # Get estimated cost for comparison
                            queue = platform.get_queue(selected_queue)
                            if queue:
                                # Estimate cost based on cores and walltime
                                cores = requirements.get("cpus", 1)
                                walltime_hours = requirements.get("walltime", 3600)
                                # Convert walltime to hours if it's a string
                                if isinstance(walltime_hours, str):
                                    from .resource_utils import parse_walltime

                                    walltime_hours = (
                                        parse_walltime(walltime_hours) or 1.0
                                    )
                                elif isinstance(walltime_hours, (int, float)):
                                    walltime_hours = (
                                        walltime_hours / 3600.0
                                    )  # Convert seconds to hours

                                cost = queue.estimate_su_cost(cores, walltime_hours)
                                if cost < best_cost:
                                    best_cost = cost
                                    best_queue = selected_queue
                    except Exception as e:
                        logging.debug(
                            f"Failed to select queue from platform {platform.name}: {e}"
                        )
                        continue

                if best_queue:
                    params["queue"] = best_queue
                    logging.info(f"Auto-selected queue: {best_queue}")
                else:
                    logging.warning(
                        "No suitable queue found for requirements, using 'normal'"
                    )
                    params["queue"] = "normal"

        except ImportError:
            logging.warning(
                "Platform system not available for auto queue selection, using 'normal'"
            )
            params["queue"] = "normal"
        except Exception as e:
            logging.warning(f"Auto queue selection failed: {e}, using 'normal'")
            params["queue"] = "normal"

    # Construct the qsub options
    options = "-N {name} -q {queue} -P {project} ".format_map(params)
    if params.get("resources"):
        options += " ".join([f"-l {resource}" for resource in params["resources"]])
    options += f" -o {joblog}"
    logging.info("Options: %s", options)

    # Load qsub options into context
    if ctx.obj is None:
        ctx.obj = {}
    ctx.obj["execdir"] = execdir
    ctx.obj["options"] = options
    ctx.obj["name"] = params["name"]
    ctx.obj["out"] = params["out"]
    ctx.obj["err"] = params["err"]
    ctx.obj["dry"] = params["dry"]
    ctx.obj["quiet"] = params["quiet"]

    # If we reach here and have an execution context, execute the job
    if has_execution_context:
        if conda_env:
            execute_conda(ctx, command, conda_env, template, pre, post)
        elif module_list:
            execute_module(ctx, command, module_list, template, pre, post)
        elif container:
            execute_singularity(ctx, command, container, bind, template, pre, post)
    elif command:
        # No execution context but command provided - use default template
        execute_default(ctx, command, template, pre, post)
    else:
        # No command and no execution context - show help or let Click handle subcommands
        pass


def execute_conda(ctx, command, env, template=None, pre=None, post=None):
    """Execute command in conda environment."""
    if not env:
        click.echo(
            "Error: No conda environment available. Either activate a conda environment or use --env to specify one.",
            err=True,
        )
        ctx.exit(2)

    # Use default template if not provided
    if not template:
        template = _get_conda_template()

    ctx_obj = ctx.obj
    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Conda environment: %s", env)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", command)
    logging.debug("Pre-commands: %s", pre)
    logging.debug("Post-commands: %s", post)

    # Construct qsub command
    cmd_str = " ".join(command)
    # Base64 encode the command to avoid escaping issues
    cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")
    submission_vars = (
        f'env={env},cmd_b64="{cmd_b64}",cwd={ctx_obj["execdir"]},'
        f'out={out},err={err},quiet={str(ctx_obj["quiet"]).lower()}'
    )
    if pre:
        pre_b64 = base64.b64encode(pre.encode("utf-8")).decode("ascii")
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_b64 = base64.b64encode(post.encode("utf-8")).decode("ascii")
        submission_vars += f',post_cmd_b64="{post_b64}"'

    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'
    logging.info("Submission command: %s", submission_command)

    if ctx_obj["dry"]:
        print(f"Dry run - would execute: {submission_command}")
        # Log history even for dry runs
        try:
            history_manager.log_execution(ctx, success=True)
        except Exception as e:
            logging.debug("Failed to log execution history: %s", e)
        return

    # Submit job and handle monitoring
    job_id = qsub(submission_command, quiet=ctx_obj["quiet"])

    # Display job ID to user (unless in quiet mode)
    if not ctx_obj["quiet"]:
        click.echo(f"🚀 Job submitted successfully! Job ID: {job_id}")

    # Log execution to history system
    try:
        history_manager.log_execution(ctx, success=True, job_id=job_id)
    except Exception as e:
        logging.debug("Failed to log execution history: %s", e)

    # Log job execution for resource tracking
    # Note: ResourceTracker expects resource_data which is available after job completion
    # For now, just log the basic command info
    try:
        resource_tracker.log_job_resources(
            job_id=job_id,
            resource_data={},  # Will be populated when job completes
            command=cmd_str,
        )
    except Exception as e:
        logging.debug("Failed to log job resources: %s", e)


def execute_module(ctx, command, modules, template=None, pre=None, post=None):
    """Execute command with environment modules."""
    if not modules:
        click.echo("Error: No environment modules specified.", err=True)
        ctx.exit(2)

    # Use default template if not provided
    if not template:
        template = _get_module_template()

    ctx_obj = ctx.obj
    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Environment modules: %s", modules)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", command)
    logging.debug("Pre-commands: %s", pre)
    logging.debug("Post-commands: %s", post)

    # Construct qsub command
    cmd_str = " ".join(command)
    mods_str = " ".join(modules)
    # Base64 encode the command to avoid escaping issues
    cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")
    submission_vars = (
        f'mods="{mods_str}",cmd_b64="{cmd_b64}",cwd={ctx_obj["execdir"]},'
        f'out={out},err={err},quiet={str(ctx_obj["quiet"]).lower()}'
    )
    if pre:
        pre_b64 = base64.b64encode(pre.encode("utf-8")).decode("ascii")
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_b64 = base64.b64encode(post.encode("utf-8")).decode("ascii")
        submission_vars += f',post_cmd_b64="{post_b64}"'

    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'
    logging.info("Submission command: %s", submission_command)

    if ctx_obj["dry"]:
        print(f"Dry run - would execute: {submission_command}")
        # Log history even for dry runs
        try:
            history_manager.log_execution(ctx, success=True)
        except Exception as e:
            logging.debug("Failed to log execution history: %s", e)
        return

    # Submit job and handle monitoring
    job_id = qsub(submission_command, quiet=ctx_obj["quiet"])

    # Display job ID to user (unless in quiet mode)
    if not ctx_obj["quiet"]:
        click.echo(f"🚀 Job submitted successfully! Job ID: {job_id}")

    # Log execution to history system
    try:
        history_manager.log_execution(ctx, success=True, job_id=job_id)
    except Exception as e:
        logging.debug("Failed to log execution history: %s", e)

    # Log job execution for resource tracking
    try:
        resource_tracker.log_job_resources(
            job_id=job_id,
            resource_data={},  # Will be populated when job completes
            command=cmd_str,
        )
    except Exception as e:
        logging.debug("Failed to log job resources: %s", e)


def execute_singularity(
    ctx, command, container, bind=None, template=None, pre=None, post=None
):
    """Execute command in Singularity container."""
    if not container:
        click.echo("Error: Singularity container path is required.", err=True)
        ctx.exit(2)

    # Use default template if not provided
    if not template:
        template = _get_singularity_template()

    ctx_obj = ctx.obj
    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    # Log parameters and context
    for key, value in ctx_obj.items():
        logging.debug("Context: %s = %s", key, value)
    logging.debug("Singularity container: %s", container)
    logging.debug("Bind mounts: %s", bind)
    logging.debug("Jobscript template: %s", template)
    logging.debug("Command: %s", command)
    logging.debug("Pre-commands: %s", pre)
    logging.debug("Post-commands: %s", post)

    # Construct qsub command
    cmd_str = " ".join(command)
    # Base64 encode the command to avoid escaping issues
    cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")
    submission_vars = (
        f'sif={container},cmd_b64="{cmd_b64}",cwd={ctx_obj["execdir"]},'
        f'out={out},err={err},quiet={str(ctx_obj["quiet"]).lower()}'
    )
    if bind:
        submission_vars += f',bind="{bind}"'
    if pre:
        pre_b64 = base64.b64encode(pre.encode("utf-8")).decode("ascii")
        submission_vars += f',pre_cmd_b64="{pre_b64}"'
    if post:
        post_b64 = base64.b64encode(post.encode("utf-8")).decode("ascii")
        submission_vars += f',post_cmd_b64="{post_b64}"'

    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'
    logging.info("Submission command: %s", submission_command)

    if ctx_obj["dry"]:
        print(f"Dry run - would execute: {submission_command}")
        # Log history even for dry runs
        try:
            history_manager.log_execution(ctx, success=True)
        except Exception as e:
            logging.debug("Failed to log execution history: %s", e)
        return

    # Submit job and handle monitoring
    job_id = qsub(submission_command, quiet=ctx_obj["quiet"])

    # Display job ID to user (unless in quiet mode)
    if not ctx_obj["quiet"]:
        click.echo(f"🚀 Job submitted successfully! Job ID: {job_id}")

    # Log execution to history system
    try:
        history_manager.log_execution(ctx, success=True, job_id=job_id)
    except Exception as e:
        logging.debug("Failed to log execution history: %s", e)

    # Log job execution for resource tracking
    try:
        resource_tracker.log_job_resources(
            job_id=job_id,
            resource_data={},  # Will be populated when job completes
            command=cmd_str,
        )
    except Exception as e:
        logging.debug("Failed to log job resources: %s", e)


def execute_default(ctx, command, template=None, pre=None, post=None):
    """Execute command using default template (no environment activation)."""
    # Use default template if not provided
    if not template:
        template = _get_default_template()

    ctx_obj = ctx.obj
    out = Path(ctx_obj["out"])
    err = Path(ctx_obj["err"])

    # Generate command string
    cmd_str = " ".join(command)

    # Base64 encode the command to safely pass it to the job script
    cmd_b64 = base64.b64encode(cmd_str.encode()).decode()

    # Encode pre and post commands if provided
    pre_cmd_b64 = base64.b64encode(pre.encode()).decode() if pre else ""
    post_cmd_b64 = base64.b64encode(post.encode()).decode() if post else ""

    # Build qsub command
    submission_vars = (
        f'cmd_b64="{cmd_b64}",cwd={ctx_obj["execdir"]},'
        f'out={out},err={err},quiet={str(ctx_obj["quiet"]).lower()}'
    )
    if pre:
        submission_vars += f',pre_cmd_b64="{pre_cmd_b64}"'
    if post:
        submission_vars += f',post_cmd_b64="{post_cmd_b64}"'

    submission_command = f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'
    logging.info("Submission command: %s", submission_command)

    if ctx_obj["dry"]:
        click.echo(f"Dry run - would execute: {submission_command}")
        # Log history even for dry runs
        try:
            history_manager.log_execution(ctx, success=True)
        except Exception as e:
            logging.debug("Failed to log execution history: %s", e)
        return

    # Submit job and handle monitoring
    job_id = qsub(submission_command, quiet=ctx_obj["quiet"])

    # Display job ID to user (unless in quiet mode)
    if not ctx_obj["quiet"]:
        click.echo(f"🚀 Job submitted successfully! Job ID: {job_id}")

    # Log execution to history system
    try:
        history_manager.log_execution(ctx, success=True, job_id=job_id)
    except Exception as e:
        logging.debug("Failed to log execution history: %s", e)

    # Log job execution for resource tracking
    try:
        resource_tracker.log_job_resources(
            job_id=job_id,
            resource_data={},  # Will be populated when job completes
            command=cmd_str,
        )
    except Exception as e:
        logging.debug("Failed to log job resources: %s", e)


def _get_default_template():
    """Get the default job script template."""
    template_path = Path(__file__).parent / "jobscripts" / "qdefault.pbs"
    if not template_path.exists():
        raise FileNotFoundError(
            f"Default template not found: {template_path}. "
            f"Looked in: {Path(__file__).parent / 'jobscripts' / 'qdefault.pbs'}"
        )
    return template_path


def _handle_remote_execution(ctx, remote_name, config_file, execdir, **params):
    """Handle remote execution via SSH."""
    try:
        from .remote_config_loader import ConfigLoadError, get_remote_config
        from .remote_executor import ConnectionError as RemoteConnectionError
        from .remote_executor import RemoteExecutorFactory

        # Get command to execute
        command = ctx.args if ctx.args else []
        if not command:
            click.echo("Error: Command is required for remote execution.", err=True)
            ctx.exit(2)

        # Load remote configuration
        try:
            remote_config = get_remote_config(remote_name, config_file)
        except ConfigLoadError as e:
            click.echo(f"Configuration error: {e}", err=True)
            ctx.exit(1)

        # Create executor
        executor = RemoteExecutorFactory.create(remote_config)

        # Test connection
        if not executor.test_connection():
            click.echo(f"Error: Cannot connect to {remote_config.hostname}", err=True)
            click.echo("Suggestions:", err=True)
            if hasattr(executor, "config") and executor.config.protocol == "ssh":
                click.echo(
                    f"  - Check SSH configuration in {remote_config.config or '~/.ssh/config'}",
                    err=True,
                )
                click.echo(
                    "  - Verify network connectivity and VPN if required", err=True
                )
                click.echo(
                    f"  - Test connection: ssh {remote_config.hostname} echo 'test'",
                    err=True,
                )
            ctx.exit(1)

        # Determine remote working directory
        local_cwd = Path(execdir).resolve()
        explicit_execdir = (
            params.get("execdir") if params.get("execdir") != os.getcwd() else None
        )
        remote_working_dir = remote_config.determine_remote_working_dir(
            local_cwd, explicit_execdir
        )

        # Build remote qxub command
        remote_args = []

        # Add platform file
        remote_args.extend(["--platform-file", remote_config.platform_file])

        # Add execution directory
        remote_args.extend(["--execdir", remote_working_dir])

        # Add other parameters (excluding remote and execdir)
        for key, value in params.items():
            if key in ["remote"] or value is None:
                continue

            # Handle special parameter formatting
            if key == "resources" and value:
                for resource in value:
                    remote_args.extend(["-l", resource])
            elif key == "mod" and value:
                for module in value:
                    remote_args.extend(["--mod", module])
            elif isinstance(value, bool):
                if value:
                    remote_args.append(f'--{key.replace("_", "-")}')
            else:
                remote_args.extend([f'--{key.replace("_", "-")}', str(value)])

        # Add the command to execute
        if command:
            remote_args.append("--")
            remote_args.extend(command)

        # Build final qxub command
        qxub_command = "qxub " + " ".join(
            f"'{arg}'" if " " in arg else arg for arg in remote_args
        )

        # Execute remotely
        exit_code = executor.execute(
            qxub_command, remote_working_dir, stream_output=True
        )
        ctx.exit(exit_code)

    except ImportError as e:
        click.echo(f"Remote execution not available: {e}", err=True)
        click.echo("This feature requires additional dependencies.", err=True)
        ctx.exit(1)
    except RemoteConnectionError as e:
        click.echo(f"Connection error: {e}", err=True)
        for suggestion in e.suggestions:
            click.echo(f"  - {suggestion}", err=True)
        ctx.exit(1)
    except Exception as e:
        click.echo(f"Remote execution failed: {e}", err=True)
        ctx.exit(1)


# CLI Management Commands (these remain as separate commands)
qxub.add_command(config_cli)
qxub.add_command(alias_cli)
qxub.add_command(history)
qxub.add_command(resources)
qxub.add_command(platform_cli)
qxub.add_command(select_queue_cmd)
qxub.add_command(validate_cmd)
qxub.add_command(estimate_cmd)
