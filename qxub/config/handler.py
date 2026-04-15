"""Configuration handling logic extracted from CLI for better separation of concerns."""

import logging

logger = logging.getLogger(__name__)
import os
from pathlib import Path

import click


def _get_default_output_dir():
    """Get appropriate output directory that works across login and compute nodes."""
    from datetime import datetime

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


def _sanitize_job_name(name):
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

    return sanitized


def _default_job_name(command):
    """Derive a default job name from the command: {cmd}-{date}-{time}."""
    from datetime import datetime

    cmd_word = "qx"
    if command:
        first = command[0] if isinstance(command, (list, tuple)) else command.split()[0]
        cmd_word = Path(first).name or "qx"

    now = datetime.now()
    return f"{cmd_word}-{now.strftime('%Y%m%d')}-{now.strftime('%H%M%S')}"


def apply_config_defaults(params, config_manager):
    """Apply configuration defaults to CLI parameters."""
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
        configured_name = defaults.get("name")
        if configured_name:
            params["name"] = configured_name
        else:
            params["name"] = _default_job_name(params.get("command"))
    if params["project"] is None:
        params["project"] = defaults.get("project", os.getenv("PROJECT"))
    if not params["resources"]:  # Empty tuple means no resources provided
        params["resources"] = defaults.get("resources", [])

    # Email notification defaults
    if params.get("email") == "__DISABLED__":
        # --no-notify was specified, suppress notifications
        params["email"] = None
    elif params.get("email") is None:
        params["email"] = defaults.get("email")
    if params.get("email_opts") is None:
        params["email_opts"] = defaults.get("email_opts", "ae")  # abort + end

    return params


def resolve_template_variables(params, config_manager):
    """Resolve template variables in parameter values."""
    template_vars = config_manager.get_template_variables(
        name=params["name"], project=params["project"], queue=params["queue"]
    )

    # Allow --log-dir CLI flag to override the configured {log_dir} template variable.
    if params.get("log_dir"):
        template_vars["log_dir"] = params["log_dir"]

    # Resolve template strings
    resolved_params = {}
    for key, value in params.items():
        if isinstance(value, str) and "{" in value:
            resolved_params[key] = config_manager.resolve_templates(
                value, template_vars
            )
        else:
            resolved_params[key] = value

    return resolved_params, template_vars


def process_job_configuration(params, config_manager):
    """Process and sanitize job configuration parameters."""
    # Apply defaults
    params = apply_config_defaults(params, config_manager)

    # Resolve template variables
    resolved_params, template_vars = resolve_template_variables(params, config_manager)
    params.update(resolved_params)

    # Sanitize job name for PBS compliance
    if params["name"]:
        params["name"] = _sanitize_job_name(params["name"])

    # Handle joblog default
    joblog = params["joblog"] or f"{params['name']}.log"
    if isinstance(joblog, str) and "{" in joblog:
        joblog = config_manager.resolve_templates(joblog, template_vars)

    # Ensure joblog is always an absolute path so that status checks
    # (which may run from a different CWD) can locate the file.
    if joblog and not os.path.isabs(joblog):
        execdir = params.get("execdir") or os.getcwd()
        joblog = os.path.join(execdir, joblog)

    params["joblog"] = joblog

    # Ensure out/err are always absolute paths.  The PBS jobscripts derive
    # STATUS_DIR from $out *before* ``cd $cwd``, so a relative path would
    # resolve against $HOME (the PBS default CWD) instead of the intended
    # execution directory, causing status file writes to fail after the cd.
    execdir = params.get("execdir") or os.getcwd()
    for key in ("out", "err"):
        val = params.get(key)
        if val and not os.path.isabs(str(val)):
            params[key] = os.path.join(execdir, str(val))

    return params


def select_auto_queue(params):
    """Handle automatic queue selection based on resource requirements."""
    if params["queue"] != "auto":
        return params

    try:
        from pathlib import Path

        from ..platforms import CorePlatformLoader as PlatformLoader
        from ..resources import parse_walltime  # noqa: F811

        # Check for QXUB_PLATFORM_PATHS environment variable
        platform_paths_env = os.environ.get("QXUB_PLATFORM_PATHS")
        if platform_paths_env:
            search_paths = [Path(p.strip()) for p in platform_paths_env.split(":")]
            loader = PlatformLoader(search_paths=search_paths)
        else:
            loader = PlatformLoader()

        platform_names = loader.list_platforms()

        if not platform_names:
            logger.warning(
                "No platforms available for auto queue selection, using 'normal'"
            )
            params["queue"] = "normal"
            return params

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
                        requirements["cpus"] = int(value)  # Use "cpus" not "ncpus"
                    elif key in ["ngpus", "gpu"]:
                        gpu_count = int(value) if value.isdigit() else 1
                        requirements["gpu_requested"] = (
                            gpu_count  # Use "gpu_requested" for conditions
                        )
                        requirements["gpus"] = gpu_count  # Keep "gpus" for validation

        # Add internet connectivity requirement if specified
        if params.get("internet"):
            requirements["internet"] = True

        # Add GPU type for queue selection if specified
        if params.get("gpu_type"):
            gpu_type = params["gpu_type"].lower()
            requirements["gpu_type"] = gpu_type

            # Validate gpu_type against platform definitions
            valid_gpu_types = set()
            for platform_name in platform_names:
                platform = loader.get_platform(platform_name)
                if not platform:
                    continue
                for q in platform.queues.values():
                    if q.gpu_type:
                        valid_gpu_types.add(q.gpu_type.lower())
            if valid_gpu_types and gpu_type not in valid_gpu_types:
                sorted_types = sorted(valid_gpu_types)
                raise click.ClickException(
                    f"Unknown GPU type '{gpu_type}'. "
                    f"Available types: {', '.join(sorted_types)}"
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
                            from ..resources import parse_walltime  # noqa: F811

                            walltime_hours = parse_walltime(walltime_hours) or 1.0
                        elif isinstance(walltime_hours, (int, float)):
                            walltime_hours = (
                                walltime_hours / 3600.0
                            )  # Convert seconds to hours

                        cost = queue.estimate_su_cost(cores, walltime_hours)
                        if cost < best_cost:
                            best_cost = cost
                            best_queue = selected_queue
            except Exception as e:
                logger.debug(
                    f"Failed to select queue from platform {platform.name}: {e}"
                )
                continue

        if best_queue:
            params["queue"] = best_queue
            logger.info("Auto-selected queue: %s", best_queue)
        else:
            logger.warning("No suitable queue found for requirements, using 'normal'")
            params["queue"] = "normal"

    except ImportError:
        logger.warning(
            "Platform system not available for auto queue selection, using 'normal'"
        )
        params["queue"] = "normal"
    except click.ClickException:
        raise
    except Exception as e:
        logger.warning("Auto queue selection failed: %s, using 'normal'", e)
        params["queue"] = "normal"

    return params


def build_qsub_options(params):
    """Build qsub options string from parameters."""
    options = "-N {name} -q {queue} -P {project} ".format_map(params)
    if params.get("resources"):
        options += " ".join([f"-l {resource}" for resource in params["resources"]])
    options += f" -o {params['joblog']}"

    # Email notifications via PBS -M and -m
    if params.get("email"):
        options += f" -M {params['email']}"
        email_opts = params.get("email_opts", "ae")  # default: abort + end
        options += f" -m {email_opts}"

    return options


def adjust_cpus_for_gpus(params):
    """Auto-set ncpus based on cpus_per_gpu when GPUs are requested without explicit CPUs."""
    if params.get("cpus_explicit"):
        return params

    queue_name = params.get("queue")
    if not queue_name or queue_name == "auto":
        return params

    # Find how many GPUs were requested
    ngpus = 0
    resources = list(params.get("resources") or [])
    for r in resources:
        if r.startswith("ngpus="):
            ngpus = int(r.split("=", 1)[1])
            break

    if ngpus <= 0:
        return params

    try:
        from pathlib import Path

        from ..platforms.core import PlatformLoader

        platform_paths_env = os.environ.get("QXUB_PLATFORM_PATHS")
        if platform_paths_env:
            search_paths = [Path(p.strip()) for p in platform_paths_env.split(":")]
            loader = PlatformLoader(search_paths=search_paths)
        else:
            loader = PlatformLoader()

        for platform_name in loader.list_platforms():
            platform = loader.get_platform(platform_name)
            if not platform:
                continue
            queue = platform.get_queue(queue_name)
            if queue and queue.cpus_per_gpu:
                required_cpus = ngpus * queue.cpus_per_gpu
                # Replace or inject ncpus in resources
                has_ncpus = any(r.startswith("ncpus=") for r in resources)
                if has_ncpus:
                    resources = [
                        f"ncpus={required_cpus}" if r.startswith("ncpus=") else r
                        for r in resources
                    ]
                else:
                    resources.append(f"ncpus={required_cpus}")
                params["resources"] = resources
                logger.info(
                    "Auto-set ncpus=%d for %d GPU(s) on queue '%s' "
                    "(%d cpus_per_gpu)",
                    required_cpus,
                    ngpus,
                    queue_name,
                    queue.cpus_per_gpu,
                )
                break
    except Exception as e:
        logger.debug("Could not adjust CPUs for GPU queue: %s", e)

    return params


def adjust_resources_for_queue(params):
    """Adjust resources to respect queue limits for explicitly specified queues.

    Mimics PBS qsub behavior: if user explicitly specifies a queue but doesn't
    explicitly specify ncpus, automatically cap ncpus at the queue's maximum.
    If user explicitly specifies ncpus that exceed the queue max, validation
    will catch it later with a helpful error.
    """
    queue_name = params.get("queue")
    if not queue_name or queue_name == "auto":
        # Auto queue selection handles this differently
        return params

    resources = params.get("resources", [])
    if not resources:
        return params

    # Check if CPUs were explicitly specified by user (tracked in exec_cli.py)
    cpus_explicit = params.get("cpus_explicit", False)

    if cpus_explicit:
        # User explicitly set CPUs - let validation handle any issues
        return params

    # User didn't explicitly set CPUs - check if we need to adjust for queue limits
    try:
        from qxub.platforms import get_current_platform

        # Get the current platform
        platform = get_current_platform()
        if platform:
            queue = platform.get_queue(queue_name)
            if queue and queue.limits.max_cpus:
                # Found the queue - check if we need to adjust CPUs
                current_ncpus = None
                for r in resources:
                    if r.startswith("ncpus="):
                        current_ncpus = int(r.split("=")[1])
                        break

                if current_ncpus and current_ncpus > queue.limits.max_cpus:
                    # Default CPU count exceeds queue limit - adjust it gracefully
                    logger.info(
                        f"Automatically adjusting ncpus from {current_ncpus} to {queue.limits.max_cpus} "
                        f"for queue '{queue_name}' (default exceeds queue maximum)"
                    )
                    # Replace the ncpus value in resources
                    new_resources = [
                        (
                            f"ncpus={queue.limits.max_cpus}"
                            if r.startswith("ncpus=")
                            else r
                        )
                        for r in resources
                    ]
                    params["resources"] = new_resources

    except Exception as e:
        # Platform system not available or other error - don't adjust
        logger.debug(f"Could not adjust resources for queue: {e}")

    return params


def process_job_options(params, config_manager):
    """Complete job option processing pipeline."""
    # Process configuration
    params = process_job_configuration(params, config_manager)

    # Handle auto queue selection
    params = select_auto_queue(params)

    # Auto-set CPUs based on GPU queue's cpus_per_gpu
    params = adjust_cpus_for_gpus(params)

    # Adjust resources for explicitly specified queues
    params = adjust_resources_for_queue(params)

    # Build qsub options
    options = build_qsub_options(params)

    logger.info("Options: %s", options)

    return params, options
