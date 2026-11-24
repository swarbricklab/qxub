"""Logging set-up for the qxub CLI."""

import logging


def setup_logging(verbosity: int = None):
    """
    Configures the logging level based on the verbosity provided by the user.

    Args:
        verbosity (int): The number of '-v' flags used, or from config.
                       - 0: ERROR level (default)
                       - 1: WARNING level
                       - 2: INFO level
                       - 3 or more: DEBUG level

    This function adjusts the logging output to provide more detailed information
    as verbosity increases, allowing users to control the granularity of log messages.
    """
    if verbosity is None:
        from .manager import config_manager

        verbosity = int(config_manager.get_config_value("verbosity") or 0)

    # Get the root logger and configure it directly
    root_logger = logging.getLogger()

    # Remove existing handlers to avoid conflicts
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Set the level and format based on verbosity
    if verbosity == 1:
        level = logging.WARNING
        format_str = "%(levelname)s: %(message)s"
    elif verbosity == 2:
        level = logging.INFO
        format_str = "%(levelname)s: %(message)s"
    elif verbosity >= 3:
        level = logging.DEBUG
        format_str = "%(levelname)s:%(name)s: %(message)s"
    else:
        level = logging.ERROR
        format_str = "%(levelname)s: %(message)s"

    # Configure the root logger
    root_logger.setLevel(level)
    handler = logging.StreamHandler()
    handler.setLevel(level)
    formatter = logging.Formatter(format_str)
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)
