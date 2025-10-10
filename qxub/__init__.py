"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

# __init__.py

__version__ = "2.0.0"

from . import config, scheduler
from .cli import qxub
