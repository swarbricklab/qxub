"""
Template management for qxub job scripts.

This module handles finding and loading job script templates for different
execution contexts (conda, modules, singularity, default).
"""

import os
from pathlib import Path
from typing import Dict


def get_template(template_type: str, custom_template: str = None) -> str:
    """
    Get the appropriate job script template path.

    Args:
        template_type: Type of template ('conda', 'module', 'singularity', 'default')
        custom_template: Optional custom template path

    Returns:
        Path to the template file

    Raises:
        FileNotFoundError: If template cannot be found
    """
    if custom_template and custom_template != "None":
        if os.path.exists(custom_template):
            return custom_template
        else:
            raise FileNotFoundError(f"Custom template not found: {custom_template}")

    # Template filename mapping
    template_files = {
        "conda": "qconda.pbs",
        "module": "qmod.pbs",
        "singularity": "qsing.pbs",
        "default": "qdefault.pbs",
    }

    if template_type not in template_files:
        raise ValueError(f"Unknown template type: {template_type}")

    template_filename = template_files[template_type]

    # Try importlib.resources first (modern approach)
    try:
        from importlib import resources

        # For Python 3.9+, use the new API
        if hasattr(resources, "files"):
            jobscripts_path = resources.files("qxub") / "jobscripts" / template_filename
            if jobscripts_path.is_file():
                return str(jobscripts_path)
        else:
            # Fallback for Python 3.8 (deprecated but still works)
            with resources.path("qxub.jobscripts", template_filename) as template_path:
                if template_path.exists():
                    return str(template_path)
    except (ImportError, AttributeError, FileNotFoundError, ModuleNotFoundError):
        pass

    # Fallback to relative path from this module
    current_dir = Path(__file__).parent
    template_path = current_dir / "jobscripts" / template_filename
    if template_path.exists():
        return str(template_path)

    # Last resort - raise an informative error
    raise FileNotFoundError(
        f"Could not locate {template_filename} template file. "
        f"Looked in: {current_dir / 'jobscripts' / template_filename}"
    )


def get_conda_template(custom_template: str = None) -> str:
    """Get conda template path."""
    return get_template("conda", custom_template)


def get_module_template(custom_template: str = None) -> str:
    """Get module template path."""
    return get_template("module", custom_template)


def get_singularity_template(custom_template: str = None) -> str:
    """Get singularity template path."""
    return get_template("singularity", custom_template)


def get_default_template(custom_template: str = None) -> str:
    """Get default template path."""
    return get_template("default", custom_template)


# Legacy function names for backward compatibility
_get_conda_template = get_conda_template
_get_module_template = get_module_template
_get_singularity_template = get_singularity_template
_get_default_template = get_default_template
