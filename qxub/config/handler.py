"""Configuration handling logic extracted from CLI for better separation of concerns."""

import logging

logger = logging.getLogger(__name__)
import os
import shlex
from dataclasses import replace
from pathlib import Path

import click


def _get_default_output_dir():
    """Get appropriate output directory that works across login and compute nodes.

    Uses ``{scratch}/qt/<timestamp>`` when the current platform defines a
    shared scratch area and it exists, otherwise ``./qt/<timestamp>``.
    """
    from datetime import datetime

    from .site import scratch_dir

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Prefer shared scratch space over TMPDIR to avoid cross-node access issues
    scratch = scratch_dir()
    if scratch is not None and scratch.exists():
        return scratch / "qt" / timestamp

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


def as_resource_spec(resources):
    """``resources`` as a ResourceSpec.

    Accepts a ResourceSpec, or PBS ``-l`` strings as found in config files
    (``defaults.resources``) and older callers.
    """
    from ..backends.pbs import parse_pbs_resources
    from ..resources.spec import ResourceSpec

    if isinstance(resources, ResourceSpec):
        return resources
    return parse_pbs_resources(resources or [])


def _resource_default(config_manager, key):
    """``defaults.<key>`` (or the deprecated root-level ``<key>``)."""
    value = config_manager.get_config_value(f"defaults.{key}")
    if value is not None:
        return value
    value = config_manager.get_config_value(key)
    if value is not None:
        logger.warning(
            "Config key '%s' at root level is deprecated. "
            "Please move to 'defaults.%s' in your config file.",
            key,
            key,
        )
    return value


def resolve_resources(
    requested,
    config_manager,
    *,
    memory=None,
    walltime=None,
    cpus=None,
    disk=None,
    gpus=None,
    volumes=None,
):
    """A job's resources: ``requested`` (a ResourceSpec from ``-l`` or a
    ``resources=`` dict), then qxub's own options, which win, then the
    configured ``defaults.mem``, ``defaults.volumes``, ... for whatever
    neither gave.

    ``qxub exec`` and :func:`qxub.execution.submit_job` both use this, so a
    job submitted from Python gets the same defaults (on Gadi, the
    ``storage=`` its conda environment lives on) as one from the CLI.
    Raises ``ValueError`` for an invalid value.
    """
    from ..resources.spec import ResourceSpec

    def value(given, key, already_requested):
        if given:
            return given
        return None if already_requested else _resource_default(config_manager, key)

    options = ResourceSpec.from_options(
        memory=value(memory, "mem", requested.memory is not None),
        walltime=value(walltime, "runtime", requested.walltime is not None),
        cpus=value(cpus, "cpus", requested.cpus is not None),
        disk=value(disk, "disk", requested.disk is not None),
        gpus=value(gpus, "gpus", requested.gpus is not None),
    )
    volumes = value(volumes, "volumes", "storage" in requested.native)
    if volumes:
        # NCI storage mounts; PBS-only for now, so kept as a native request
        options.native["storage"] = str(volumes)
    return requested.updated(options)


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
        from .site import default_project

        params["project"] = defaults.get("project") or default_project()
    params["resources"] = as_resource_spec(params.get("resources"))
    if not params["resources"]:  # nothing requested: use defaults.resources
        params["resources"] = as_resource_spec(defaults.get("resources", []))

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


def _warn_auto_queue_fallback(reason: str) -> None:
    """Tell the user that queue ``auto`` fell back to ``normal``.

    Printed to stderr so it is seen without ``-v`` and does not mix with
    output scripts read from stdout (e.g. the job ID from ``--terse``).
    """
    message = f"{reason}, using queue 'normal'"
    logger.debug(message)
    click.echo(f"Warning: {message}", err=True)


def _walltime_is_the_only_problem(loader, platform_names, requirements) -> bool:
    """Whether some queue would fit if not for the requested walltime."""
    if not requirements.get("walltime"):
        return False
    without_walltime = {k: v for k, v in requirements.items() if k != "walltime"}
    for name in platform_names:
        platform = loader.get_platform(name)
        if platform and platform.select_queue(without_walltime):
            return True
    return False


def select_auto_queue(params):
    """Handle automatic queue selection based on resource requirements."""
    if params["queue"] != "auto":
        return params

    params["queue_was_auto"] = True

    try:
        from pathlib import Path

        from ..platforms import CorePlatformLoader as PlatformLoader
        from ..resources.spec import format_duration, format_size

        # Check for QXUB_PLATFORM_PATHS environment variable
        platform_paths_env = os.environ.get("QXUB_PLATFORM_PATHS")
        if platform_paths_env:
            search_paths = [Path(p.strip()) for p in platform_paths_env.split(":")]
            loader = PlatformLoader(search_paths=search_paths)
        else:
            loader = PlatformLoader()

        platform_names = loader.list_platforms()

        if not platform_names:
            _warn_auto_queue_fallback("No platforms available for auto queue selection")
            params["queue"] = "normal"
            return params

        # Build requirements from resources. Sizes and durations are passed
        # as strings, as the platform's selection conditions expect.
        requirements = {}
        resources = as_resource_spec(params.get("resources"))
        params["resources"] = resources
        if resources.memory is not None:
            requirements["memory"] = format_size(resources.memory)
        if resources.walltime is not None:
            requirements["walltime"] = format_duration(resources.walltime)
        if resources.cpus is not None:
            requirements["cpus"] = resources.cpus
        if resources.gpus is not None:
            requirements["gpu_requested"] = resources.gpus  # for conditions
            requirements["gpus"] = resources.gpus  # for validation

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
                        walltime_hours = (resources.walltime or 3600) / 3600.0

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
        elif _walltime_is_the_only_problem(loader, platform_names, requirements):
            # Falling back to 'normal' would only get the job rejected by qsub.
            raise click.ClickException(
                f"No queue allows a walltime of {requirements['walltime']} "
                f"with {requirements.get('cpus', 1)} CPU(s). "
                "Request a shorter --runtime, or name a queue with --queue."
            )
        else:
            _warn_auto_queue_fallback("No suitable queue found for requirements")
            params["queue"] = "normal"

    except ImportError:
        _warn_auto_queue_fallback(
            "Platform system not available for auto queue selection"
        )
        params["queue"] = "normal"
    except click.ClickException:
        raise
    except Exception as e:
        _warn_auto_queue_fallback(f"Auto queue selection failed: {e}")
        params["queue"] = "normal"

    return params


def with_status_dir_storage(resources, project):
    """``resources`` plus the storage the job needs to reach qxub's status
    directory, if it is charged to another project (see
    :func:`qxub.config.site.status_dir_storage`)."""
    from .site import status_dir_storage

    needed = status_dir_storage(project)
    if not needed:
        return resources
    volumes = [v for v in resources.native.get("storage", "").split("+") if v]
    if needed in volumes:
        return resources
    return replace(
        resources, native={**resources.native, "storage": "+".join(volumes + [needed])}
    )


def job_spec(params, script="", env=None):
    """The :class:`~qxub.backends.JobSpec` for processed job options."""
    from ..backends import JobSpec
    from .site import default_directives

    return JobSpec(
        script=str(script),
        name=params["name"],
        queue=params.get("queue"),
        account=params.get("project") or None,
        resources=with_status_dir_storage(
            as_resource_spec(params.get("resources")), params.get("project")
        ),
        job_log=str(params["joblog"]) if params.get("joblog") else None,
        mail_to=params.get("email") or None,
        mail_events=params.get("email_opts") or "ae",  # default: abort + end
        env=dict(env or {}),
        # Site-wide scheduler arguments from the platform, e.g. -W umask=0027
        native_args=default_directives(),
    )


def build_qsub_options(params):
    """The scheduler options for ``params``, as one string (for logs).

    Submission uses :func:`job_spec` and the platform's backend; this string
    is only kept for messages and for callers outside qxub.
    """
    from ..backends import get_backend

    argv = get_backend().submit_command(job_spec(params))
    return shlex.join(argv[1:-1])


def adjust_cpus_for_gpus(params):
    """Auto-set ncpus based on cpus_per_gpu when GPUs are requested without explicit CPUs."""
    if params.get("cpus_explicit"):
        return params

    queue_name = params.get("queue")
    if not queue_name or queue_name == "auto":
        return params

    resources = as_resource_spec(params.get("resources"))
    ngpus = resources.gpus or 0
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
                params["resources"] = replace(resources, cpus=required_cpus)
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


def enforce_queue_minimums(params):
    """Override user-specified resources that violate the selected queue's minimums.

    Applies to both auto-selected queues (--queue auto) and explicitly named queues.
    Emits a user-visible warning for each override so the user knows what changed.

    CPU enforcement:
    - If the queue has cpus_per_gpu, CPUs must be a positive multiple of that value.
      Round up non-multiples; raise to the minimum if below it.
    - If the queue has limits.min_cpus, enforce that floor independently.

    Disk (jobfs) enforcement:
    - If the queue has limits.min_local_storage and the user specified jobfs below
      that threshold (or didn't specify it at all), inject/raise the minimum.
    """
    queue_name = params.get("queue")
    if not queue_name or queue_name == "auto":
        return params

    try:
        from pathlib import Path

        from ..platforms.core import PlatformLoader
        from ..resources import parse_memory_size
        from ..resources.spec import format_size

        platform_paths_env = os.environ.get("QXUB_PLATFORM_PATHS")
        if platform_paths_env:
            search_paths = [Path(p.strip()) for p in platform_paths_env.split(":")]
            loader = PlatformLoader(search_paths=search_paths)
        else:
            loader = PlatformLoader()

        queue = None
        for platform_name in loader.list_platforms():
            platform = loader.get_platform(platform_name)
            if platform:
                queue = platform.get_queue(queue_name)
                if queue:
                    break

        if not queue:
            return params

        was_auto = params.get("queue_was_auto", False)
        queue_context = (
            f"auto-selected queue '{queue_name}'"
            if was_auto
            else f"queue '{queue_name}'"
        )

        resources = as_resource_spec(params.get("resources"))

        # --- CPU enforcement ---
        if queue.cpus_per_gpu or queue.limits.min_cpus:
            current_ncpus = resources.cpus

            if current_ncpus is not None:
                adjusted = current_ncpus

                # Enforce minimum from cpus_per_gpu (minimum = 1 * cpus_per_gpu)
                if queue.cpus_per_gpu:
                    min_from_ratio = queue.cpus_per_gpu
                    if adjusted < min_from_ratio:
                        adjusted = min_from_ratio
                    elif adjusted % queue.cpus_per_gpu != 0:
                        adjusted = (
                            (adjusted // queue.cpus_per_gpu) + 1
                        ) * queue.cpus_per_gpu

                # Enforce absolute minimum from limits
                if queue.limits.min_cpus and adjusted < queue.limits.min_cpus:
                    adjusted = queue.limits.min_cpus

                if adjusted != current_ncpus:
                    reason = (
                        f"must be a multiple of {queue.cpus_per_gpu} (cpus_per_gpu)"
                        if queue.cpus_per_gpu
                        and current_ncpus % queue.cpus_per_gpu != 0
                        else f"minimum is {adjusted}"
                    )
                    click.echo(
                        click.style(
                            f"Warning: ncpus={current_ncpus} overridden to ncpus={adjusted} "
                            f"for {queue_context} ({reason})",
                            fg="yellow",
                        ),
                        err=True,
                    )
                    resources = replace(resources, cpus=adjusted)
                    params["resources"] = resources

        # --- Disk (jobfs) enforcement ---
        if queue.limits.min_local_storage:
            min_bytes = parse_memory_size(queue.limits.min_local_storage)
            current_jobfs = resources.disk

            if current_jobfs is not None and min_bytes and current_jobfs < min_bytes:
                click.echo(
                    click.style(
                        f"Warning: jobfs={format_size(current_jobfs)} overridden to "
                        f"jobfs={queue.limits.min_local_storage} "
                        f"for {queue_context} (minimum is {queue.limits.min_local_storage})",
                        fg="yellow",
                    ),
                    err=True,
                )
                params["resources"] = replace(resources, disk=min_bytes)
            elif current_jobfs is None and min_bytes:
                # User didn't specify jobfs at all — inject the minimum silently
                params["resources"] = replace(resources, disk=min_bytes)

    except Exception as e:
        logger.debug("Could not enforce queue minimums: %s", e)

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

    resources = as_resource_spec(params.get("resources"))
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
                current_ncpus = resources.cpus

                if current_ncpus and current_ncpus > queue.limits.max_cpus:
                    # Default CPU count exceeds queue limit - adjust it gracefully
                    logger.info(
                        f"Automatically adjusting ncpus from {current_ncpus} to {queue.limits.max_cpus} "
                        f"for queue '{queue_name}' (default exceeds queue maximum)"
                    )
                    params["resources"] = replace(resources, cpus=queue.limits.max_cpus)

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

    # Enforce queue minimum resource requirements, overriding non-compliant user values
    params = enforce_queue_minimums(params)

    # Adjust resources for explicitly specified queues
    params = adjust_resources_for_queue(params)

    # Build qsub options
    options = build_qsub_options(params)

    logger.info("Options: %s", options)

    return params, options
