#!/usr/bin/env python3

# Prototype for smart quote processing in --cmd
import os
import re


def process_smart_cmd(cmd_str: str) -> str:
    """
    Enhanced quote processing for --cmd option.

    If command uses double quotes as outer delimiters, apply smart escaping:
    - \" becomes "
    - \$ becomes literal $
    - ${var} still expands at submission time
    - ${{var}} still converts for execution time
    - Single quotes preserved literally
    """

    # Detect double-quote wrapped commands
    if cmd_str.startswith('"') and cmd_str.endswith('"') and len(cmd_str) > 1:
        # Remove outer double quotes
        inner_cmd = cmd_str[1:-1]

        # Process escaped characters
        inner_cmd = inner_cmd.replace('\\"', '"')  # \" -> "

        # Handle literal dollars: \$ -> $
        inner_cmd = inner_cmd.replace("\\$", "__LITERAL_DOLLAR__")

        # Apply existing variable expansion (${var} and ${{var}})
        processed = expand_submission_variables(inner_cmd)

        # Restore literal dollars
        processed = processed.replace("__LITERAL_DOLLAR__", "$")

        return processed

    # No double quotes - use existing processing
    return expand_submission_variables(cmd_str)


def expand_submission_variables(cmd_str: str) -> str:
    """
    Existing variable expansion function (simplified for testing)
    """

    def replace_submission_var(match):
        var_name = match.group(1)
        value = os.environ.get(var_name, f"${{{var_name}}}")  # Keep if not found
        return value

    def replace_execution_var(match):
        var_name = match.group(1)
        return f"${{{var_name}}}"

    # Process execution variables first: ${{var}} -> ${var}
    result = re.sub(r"\$\{\{([^}]+)\}\}", replace_execution_var, cmd_str)

    # Then submission variables: ${var} -> value
    result = re.sub(r"\$\{([^}]+)\}", replace_submission_var, result)

    return result


# Test cases
def test_smart_quotes():
    print("🧪 Testing Smart Quote Processing")
    print()

    test_cases = [
        # Traditional single quotes (backward compatible)
        ('find /path -name "*.txt"', "Single quotes (existing)"),
        # Smart double quotes
        (
            '"find /path -exec sh -c \'echo \\"File: \\$1\\"!\' _ {} \\;"',
            "Smart double quotes",
        ),
        # With variables
        (
            "\"find ${HOME} -exec echo 'User: ${USER}, Job: ${{PBS_JOBID}}' \\;\"",
            "With variables",
        ),
        # Complex find command
        (
            '"find /data -name \\"*.log\\" -exec awk \'{print \\"Line\\" NR \\": \\" \\$0}\' {} \\;"',
            "Complex AWK",
        ),
        # JSON-like output
        (
            '"find /path -exec sh -c \'echo \\"{\\\\"file\\\\": \\\\"\\$1\\\\", \\\\"cost\\\\": \\\\"\\$100\\\\"}\\"!\' _ {} \\;"',
            "JSON output",
        ),
    ]

    for cmd, description in test_cases:
        print(f"Test: {description}")
        print(f"Input:  {cmd}")
        print(f"Output: {process_smart_cmd(cmd)}")
        print()


if __name__ == "__main__":
    test_smart_quotes()
