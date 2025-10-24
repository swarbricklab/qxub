"""Configuration handling logic extracted from CLI for better separation of concerns."""

import logging
import os
from pathlib import Path


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
        params["name"] = defaults.get("name", "qt")
    if params["project"] is None:
        params["project"] = defaults.get("project", os.getenv("PROJECT"))
    if not params["resources"]:  # Empty tuple means no resources provided
        params["resources"] = defaults.get("resources", [])

    return params


def resolve_template_variables(params, config_manager):
    """Resolve template variables in parameter values."""
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

    params["joblog"] = joblog

    return params


def select_auto_queue(params):
    """Handle automatic queue selection based on resource requirements."""
    if params["queue"] != "auto":
        return params

    try:
        from pathlib import Path

        from ..resources import parse_walltime  # noqa: F811
        from .platform import PlatformLoader

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
                logging.debug(
                    f"Failed to select queue from platform {platform.name}: {e}"
                )
                continue

        if best_queue:
            params["queue"] = best_queue
            logging.info(f"Auto-selected queue: {best_queue}")
        else:
            logging.warning("No suitable queue found for requirements, using 'normal'")
            params["queue"] = "normal"

    except ImportError:
        logging.warning(
            "Platform system not available for auto queue selection, using 'normal'"
        )
        params["queue"] = "normal"
    except Exception as e:
        logging.warning(f"Auto queue selection failed: {e}, using 'normal'")
        params["queue"] = "normal"

    return params


def build_qsub_options(params):
    """Build qsub options string from parameters."""
    options = "-N {name} -q {queue} -P {project} ".format_map(params)
    if params.get("resources"):
        options += " ".join([f"-l {resource}" for resource in params["resources"]])
    options += f" -o {params['joblog']}"
    return options


def process_job_options(params, config_manager):
    """Complete job option processing pipeline."""
    # Process configuration
    params = process_job_configuration(params, config_manager)

    # Handle auto queue selection
    params = select_auto_queue(params)

    # Build qsub options
    options = build_qsub_options(params)

    logging.info("Options: %s", options)

    return params, options
