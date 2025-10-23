"""
Standalone command aliases for qxub.

Provides standalone commands that replicate common bash aliases:
- qx: equivalent to 'qxub exec --'
- qxtat: equivalent to 'qxub status'
- qxet: equivalent to 'qxub config shortcut set'

These are installed as separate executables alongside qxub.
"""

import sys

import click


def qx_main():
    """Standalone 'qx' command equivalent to 'qxub exec --'"""
    from ..exec_cli import exec_cli

    # Import the exec command and invoke it with the provided arguments
    # This mimics the behavior of 'qxub exec --' exactly
    try:
        exec_cli.main(sys.argv[1:], standalone_mode=True)
    except SystemExit as e:
        # Re-raise SystemExit to preserve exit codes
        sys.exit(e.code)


def qxtat_main():
    """Standalone 'qxtat' command equivalent to 'qxub status'"""
    from ..status_cli import status_cli

    # Import the status command group and invoke it with the provided arguments
    try:
        status_cli.main(sys.argv[1:], standalone_mode=True)
    except SystemExit as e:
        # Re-raise SystemExit to preserve exit codes
        sys.exit(e.code)


def qxet_main():
    """Standalone 'qxet' command equivalent to 'qxub config shortcut set'"""
    from ..config_cli import shortcut_config

    # Import the shortcuts set command and invoke it with the provided arguments
    # This gives access to the 'set' subcommand directly
    try:
        # Get the 'set' command from the shortcut config CLI group
        set_command = shortcut_config.commands.get("set")
        if set_command:
            set_command.main(sys.argv[1:], standalone_mode=True)
        else:
            click.echo("Error: shortcut set command not found", err=True)
            sys.exit(1)
    except SystemExit as e:
        # Re-raise SystemExit to preserve exit codes
        sys.exit(e.code)
