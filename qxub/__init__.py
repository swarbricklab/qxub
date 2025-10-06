"""
This package provides tools for running qsub jobs on PBS Pro in particular environments.
This avoids boilerplate code for activating environments and switching directories, etc.
In simple cases, the need to create a jobscript can be eliminated entirely.
"""

# __init__.py

__version__ = "1.1.0"

from .cli import qxub
from . import conda
from . import module
from . import sing
from . import scheduler
from . import config
