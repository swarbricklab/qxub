"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

import base64
import difflib
import logging
import os
import signal
import sys
from datetime import datetime
from pathlib import Path

import click

from .alias_cli import alias_cli
from .config import setup_logging
from .config_cli import config_cli
from .config_manager import config_manager

# Import signal handler from execution module
from .execution import _signal_handler, validate_execution_context
from .executors import (
    execute_conda,
    execute_default,
    execute_module,
    execute_singularity,
)
from .history_cli import history
from .history_manager import history_manager
from .parameters import build_qsub_options, process_parameters
from .platform_cli import estimate_cmd, platform_cli, select_queue_cmd, validate_cmd
from .resource_tracker import resource_tracker
from .resources_cli import resources
from .scheduler import get_job_resource_data, monitor_and_tail, print_status, qdel, qsub


class QxubGroup(click.Group):
    """Custom Click group with enhanced error handling for unknown options."""

    def get_command(self, ctx, cmd_name):
        """Override get_command to handle execution context."""
        # Check if we have execution context FIRST
        execution_options = ["default", "env", "mod", "mods", "sif"]
        has_execution_context = any(ctx.params.get(opt) for opt in execution_options)

        # If we have execution context, don't resolve as subcommand
        if has_execution_context:
            return None

        return super().get_command(ctx, cmd_name)

    def invoke(self, ctx):
        """Override invoke to handle execution contexts."""
        # Check if we have execution context
        execution_options = ["default", "env", "mod", "mods", "sif"]
        has_execution_context = any(ctx.params.get(opt) for opt in execution_options)

        if has_execution_context and ctx.protected_args:
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
                    ctx.params["default"],
                    ctx.params["env"],
                    ctx.params["mod"],
                    ctx.params["mods"],
                    ctx.params["sif"],
                    ctx.params["bind"],
                    ctx.params["template"],
                    ctx.params["pre"],
                    ctx.params["post"],
                    ctx.params["cmd"],
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
                            "default",
                            "env",
                            "mod",
                            "mods",
                            "sif",
                            "bind",
                            "template",
                            "pre",
                            "post",
                            "cmd",
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
@click.option(
    "--default", is_flag=True, help="Default execution (no special environment)"
)
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
@click.option(
    "--cmd",
    help="Command to execute (supports ${var} for submission-time and ${{var}} for execution-time variables)",
)
@click.pass_context
def qxub(
    ctx,
    remote,
    config,
    execdir,
    verbose,
    version,
    default,
    env,
    mod,
    mods,
    sif,
    bind,
    template,
    pre,
    post,
    cmd,
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
        qxub --default -- echo "hello world"
        qxub --env myenv -- python script.py
        qxub --mod python3 --mod gcc -- make
        qxub --mods python3,gcc -- python script.py
        qxub --sif container.sif -- python script.py

    For complex commands with variables, use --cmd:
        qxub --env myenv --cmd "python script.py --input ${HOME}/data.txt"
        qxub --env myenv --cmd "echo 'Job ${{PBS_JOBID}} on node ${{HOSTNAME}}'"

    Use --queue auto for intelligent queue selection:
        qxub --queue auto -l mem=500GB --env myenv -- python big_job.py
        qxub --queue auto -l mem=500GB --default -- my_binary
    """
    # Handle version flag first
    if version:
        from . import __version__

        click.echo(f"qxub {__version__}")
        ctx.exit()

    # Handle remote execution early
    if remote:
        # Gather all parameters for remote execution
        all_params = dict(params)
        all_params.update(
            {
                "verbose": verbose,
                "version": version,
                "env": env,
                "mod": mod,
                "mods": mods,
                "sif": sif,
                "bind": bind,
                "template": template,
                "pre": pre,
                "post": post,
            }
        )
        return _handle_remote_execution(ctx, remote, config, execdir, **all_params)

    setup_logging(verbose)
    logging.debug("Execution directory: %s", execdir)
    logging.debug("Remaining args: %s", ctx.args)
    logging.debug("Invoked subcommand: %s", ctx.invoked_subcommand)

    # Get remaining arguments (these would be the command to execute)
    command = tuple(ctx.args) if ctx.args else tuple()

    # Handle --cmd vs traditional -- syntax
    if cmd and command:
        click.echo(
            "Error: Cannot specify both --cmd and command after --. Use either:\n"
            '  qxub --env base --cmd "command with ${vars}"\n'
            "  qxub --env base -- command args",
            err=True,
        )
        ctx.exit(2)
    elif cmd:
        # Convert --cmd string to tuple for compatibility
        command = (cmd,)

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
    execution_contexts = [default, conda_env, module_list, container]
    has_execution_context, context_type = validate_execution_context(
        default, conda_env, module_list, container
    )

    # If a subcommand is being invoked, just set up context and let Click handle it
    if ctx.invoked_subcommand is not None:
        logging.debug("Subcommand '%s' will be invoked", ctx.invoked_subcommand)
        # Fall through to context setup and let Click handle the subcommand
    elif has_execution_context:
        # We have execution context - this is direct execution
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

    # Process parameters using the new parameter module
    params = process_parameters(params)

    # Build qsub options
    options = build_qsub_options(params)
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
    ctx.obj["verbose"] = verbose

    # If we reach here and have an execution context, execute the job
    if has_execution_context:
        if context_type == "conda":
            execute_conda(ctx, command, conda_env, template, pre, post)
        elif context_type == "module":
            execute_module(ctx, command, module_list, template, pre, post)
        elif context_type == "singularity":
            execute_singularity(ctx, command, container, bind, template, pre, post)
        elif context_type == "default":
            execute_default(ctx, command, template, pre, post)
    elif command:
        # No execution context but command provided - use default template
        execute_default(ctx, command, template, pre, post)
    else:
        # No command and no execution context - show help or let Click handle subcommands
        pass


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

        # Get verbose level for output (early, so we can show info even if connection fails)
        verbose = params.get("verbose", 0)

        # Determine remote working directory early for verbose output
        local_cwd = Path(execdir).resolve()
        explicit_execdir = (
            params.get("execdir") if params.get("execdir") != os.getcwd() else None
        )
        remote_working_dir = remote_config.determine_remote_working_dir(
            local_cwd, explicit_execdir
        )

        # Always show some remote execution info if verbose
        if verbose >= 1:
            click.echo(f"🌐 Remote execution to: {remote_config.url}")
            click.echo(f"📁 Remote working directory: {remote_working_dir}")
            click.echo(f"🐍 Remote conda environment: {remote_config.qxub_env}")
            click.echo(f"📋 Platform file: {remote_config.platform_file}")

            # Show SSH connection details if verbose >= 2
            if verbose >= 2:
                ssh_details = []
                if remote_config.username:
                    ssh_details.append(f"user={remote_config.username}")
                if remote_config.port and remote_config.port != 22:
                    ssh_details.append(f"port={remote_config.port}")
                if remote_config.config:
                    ssh_details.append(f"config={remote_config.config}")
                else:
                    ssh_details.append("config=~/.ssh/config (default)")

                if ssh_details:
                    click.echo(f"🔑 SSH details: {', '.join(ssh_details)}")

        # Build remote qxub command early (for verbose output)
        remote_args = []

        # The platform file will be set via environment variable, not CLI argument
        # This avoids the "--platform-file" invalid option error

        # Add execution directory
        remote_args.extend(["--execdir", remote_working_dir])

        # Add other parameters (excluding remote-specific and config params)
        for key, value in params.items():
            # Skip parameters that are remote-execution specific or shouldn't be forwarded
            if key in ["remote", "config", "execdir"] or value is None:
                continue

            # Handle special parameter formatting
            if key == "resources":
                # Only add resources if there are actual resource specifications
                if value and len(value) > 0:
                    for resource in value:
                        remote_args.extend(["-l", resource])
            elif key == "mod":
                # Only add modules if there are actual module specifications
                if value and len(value) > 0:
                    for module in value:
                        remote_args.extend(["--mod", module])
            elif isinstance(value, bool):
                if value:
                    remote_args.append(f'--{key.replace("_", "-")}')
            elif isinstance(value, (list, tuple)) and len(value) == 0:
                # Skip empty lists/tuples to avoid adding empty parameters
                continue
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

        if verbose >= 2:
            click.echo(f"🚀 Remote command: {qxub_command}")

        # Execute remotely
        is_dry_run = params.get("dry", False)
        if is_dry_run:
            click.echo(f"🧪 Dry run - would execute remotely: {qxub_command}")
            ctx.exit(0)

        # Create executor only if not dry-run
        executor = RemoteExecutorFactory.create(remote_config)

        # Show connection attempt info if verbose
        if verbose >= 1:
            click.echo(f"🔗 Testing SSH connection to {remote_config.hostname}...")
            if verbose >= 2:
                connection_details = []
                if remote_config.username:
                    connection_details.append(f"User: {remote_config.username}")
                if remote_config.port:
                    connection_details.append(f"Port: {remote_config.port}")
                if remote_config.config:
                    connection_details.append(f"Config: {remote_config.config}")
                if connection_details:
                    click.echo(
                        f"   Connection details: {', '.join(connection_details)}"
                    )

        # Test connection
        if not executor.test_connection():
            click.echo(
                f"❌ Error: Cannot connect to {remote_config.hostname}", err=True
            )

            # Show detailed error information if available
            if hasattr(executor, "get_connection_error_details"):
                error_details = executor.get_connection_error_details()
                if error_details:
                    if "returncode" in error_details:
                        click.echo(
                            f"   SSH exit code: {error_details['returncode']}", err=True
                        )
                        if error_details["stderr"]:
                            click.echo(
                                f"   SSH error: {error_details['stderr']}", err=True
                            )
                        if verbose >= 2 and error_details["command"]:
                            click.echo(
                                f"   SSH command: {error_details['command']}", err=True
                            )
                    elif "error" in error_details:
                        click.echo(f"   Error type: {error_details['error']}", err=True)
                        if error_details.get("message"):
                            click.echo(
                                f"   Details: {error_details['message']}", err=True
                            )

            click.echo("💡 Suggestions:", err=True)
            if hasattr(executor, "config") and executor.config.protocol == "ssh":
                click.echo(
                    f"   • Check SSH configuration in {remote_config.config or '~/.ssh/config'}",
                    err=True,
                )
                click.echo(
                    "   • Verify network connectivity and VPN if required", err=True
                )
                click.echo(
                    f"   • Test connection manually: ssh {remote_config.hostname} echo 'test'",
                    err=True,
                )
                if remote_config.username:
                    click.echo(
                        f"   • Test with explicit user: ssh {remote_config.username}@{remote_config.hostname} echo 'test'",
                        err=True,
                    )
                click.echo(
                    "   • Check SSH key permissions: chmod 600 ~/.ssh/id_*", err=True
                )
                click.echo(
                    "   • Add verbose SSH debugging: ssh -vvv hostname", err=True
                )
            ctx.exit(1)

        # Show successful connection if verbose
        if verbose >= 1:
            click.echo(f"✅ SSH connection successful to {remote_config.hostname}")

        # Execute remotely (not dry-run, so actually execute)
        if verbose >= 1:
            click.echo(f"🚀 Executing remote command...")
            if verbose >= 2:
                click.echo(
                    f"   SSH target: {remote_config.username}@{remote_config.hostname}"
                    if remote_config.username
                    else f"   SSH target: {remote_config.hostname}"
                )
                click.echo(f"   Working directory: {remote_working_dir}")

        try:
            exit_code = executor.execute(
                qxub_command, remote_working_dir, stream_output=True, verbose=verbose
            )
            if verbose >= 1:
                click.echo(f"✅ Remote execution completed with exit code: {exit_code}")
            ctx.exit(exit_code)
        except click.exceptions.Exit:
            # Re-raise Click's exit exceptions (this is normal behavior)
            raise
        except Exception as e:
            click.echo(f"Remote execution failed: {e}", err=True)
            ctx.exit(1)

    except ImportError as e:
        click.echo(f"Remote execution not available: {e}", err=True)
        click.echo("This feature requires additional dependencies.", err=True)
        ctx.exit(1)
    except RemoteConnectionError as e:
        click.echo(f"Connection error: {e}", err=True)
        for suggestion in e.suggestions:
            click.echo(f"  - {suggestion}", err=True)
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
