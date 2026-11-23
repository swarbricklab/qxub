"""Helpers for building jobs: variable expansion and context validation.

Submission itself is :func:`qxub.execution.submit.dispatch`.
"""

import logging
import os
import re
from typing import Tuple

import click

logger = logging.getLogger(__name__)


def expand_submission_variables(cmd_str: str) -> str:
    """
    Expand ${var} at submission time, convert ${{var}} to ${var} for execution.
    All other $ usage is preserved unchanged.

    Enhanced with smart quote processing: if command starts and ends with double quotes,
    applies automatic escaping for better readability.

    Args:
        cmd_str: Command string with potential variables

    Returns:
        Command string with submission variables expanded and execution variables prepared

    Examples:
        # Traditional processing (backward compatible)
        expand_submission_variables('echo "${USER} on ${{HOSTNAME}}"')
        # Returns: 'echo "jr9959 on ${HOSTNAME}"'

        # Smart quote processing (enhanced)
        expand_submission_variables('"find /path -exec echo \\"User: ${USER}\\" \\;"')
        # Returns: 'find /path -exec echo "User: jr9959" \\;'

        # Literal dollar preservation
        expand_submission_variables('"awk \\"{print \\\\$1}\\" file"')
        # Returns: 'awk "{print $1}" file'
    """

    # Check for smart quote processing (double-quote wrapped commands)
    if cmd_str.startswith('"') and cmd_str.endswith('"') and len(cmd_str) > 1:
        # Remove outer double quotes
        inner_cmd = cmd_str[1:-1]

        # Process escaped characters first
        inner_cmd = inner_cmd.replace('\\"', '"')  # \" -> "

        # Handle literal dollars: \$ -> temporary placeholder
        inner_cmd = inner_cmd.replace("\\$", "__QXUB_LITERAL_DOLLAR__")

        # Apply variable expansion to the processed command
        processed = _apply_variable_expansion(inner_cmd)

        # Restore literal dollars
        processed = processed.replace("__QXUB_LITERAL_DOLLAR__", "$")

        return processed

    # No smart quotes - use standard processing
    return _apply_variable_expansion(cmd_str)


def _apply_variable_expansion(cmd_str: str) -> str:
    """
    Internal function to apply variable expansion rules.
    Separated for use by both smart quote and standard processing.
    """
    import logging
    import os
    import re

    def replace_submission_var(match):
        """Replace ${var} with environment value."""
        var_name = match.group(1)
        value = os.environ.get(var_name)
        if value is not None:
            return value
        else:
            logger.warning(
                f"Submission variable ${{{var_name}}} not found in environment"
            )
            return match.group(0)  # Keep original ${var}

    def replace_execution_var(match):
        """Convert ${{var}} to ${var} for execution-time expansion."""
        var_name = match.group(1)
        return f"${{{var_name}}}"

    # Process execution variables first to avoid conflicts
    result = re.sub(r"\$\{\{([^}]+)\}\}", replace_execution_var, cmd_str)

    # Then process submission variables
    result = re.sub(r"\$\{([^}]+)\}", replace_submission_var, result)

    return result


def expand_variables_preserving_quotes(cmd_str: str, is_list: bool = False) -> str:
    """
    Expand shell variables while preserving quote structure.

    This function expands variables like $USER and ${HOME} at submission time
    while preserving the exact quoting structure the user intended.
    It does not use shell evaluation to avoid quote mangling.

    Args:
        cmd_str: Command string potentially containing variables

    Returns:
        Command string with variables expanded but quotes preserved

    Examples:
        expand_variables_preserving_quotes('echo "Hello $USER"')
        # Returns: 'echo "Hello jr9959"'

        expand_variables_preserving_quotes('python -c "print(\'$USER\')"')
        # Returns: 'python -c "print(\'jr9959\')"'

        expand_variables_preserving_quotes('echo "Job: \\$PBS_JOBID"')
        # Returns: 'echo "Job: $PBS_JOBID"' (escaped variable preserved)
    """
    # Handle escaped variables first: \$VAR should become $VAR (literal)
    # We need to do this before expanding other variables
    result = cmd_str

    # Replace \$VAR with a temporary placeholder to protect it from expansion
    escaped_vars = []
    escaped_pattern = r"\\(\$\{[^}]+\}|\$[A-Za-z_][A-Za-z0-9_]*)"

    def save_escaped_var(match):
        escaped_vars.append(match.group(1))  # Save the $VAR part
        return f"__ESCAPED_VAR_{len(escaped_vars) - 1}__"

    result = re.sub(escaped_pattern, save_escaped_var, result)

    # Now expand unescaped variables
    def replace_var(match):
        var_name = match.group(1) or match.group(2)  # Handle both ${VAR} and $VAR
        value = os.environ.get(var_name)

        if value is not None:
            return value
        else:
            # Variable not found - keep as-is for potential execution-time expansion
            return match.group(0)

    # Match ${VAR} and $VAR patterns
    pattern = r"\$\{([^}]+)\}|\$([A-Za-z_][A-Za-z0-9_]*)"
    result = re.sub(pattern, replace_var, result)

    # Restore escaped variables (without the escape backslash)
    for i, escaped_var in enumerate(escaped_vars):
        result = result.replace(f"__ESCAPED_VAR_{i}__", escaped_var)

    return result


def validate_execution_context(
    default_flag, conda_env, module_list, container
) -> Tuple[bool, str]:
    """
    Validate execution context parameters.

    Args:
        default_flag: Default execution flag
        conda_env: Conda environment name
        module_list: List of modules to load
        container: Singularity container path

    Returns:
        Tuple of (has_context, context_type)

    Raises:
        click.ClickException: If multiple execution contexts are specified
    """
    execution_contexts = [default_flag, conda_env, module_list, container]
    context_count = sum(bool(x) for x in execution_contexts)

    if context_count > 1:
        raise click.ClickException("Cannot specify multiple execution contexts")

    has_context = context_count > 0

    # Determine context type
    if default_flag:
        return has_context, "default"
    elif conda_env:
        return has_context, "conda"
    elif module_list:
        return has_context, "module"
    elif container:
        return has_context, "singularity"
    else:
        return has_context, "default"
