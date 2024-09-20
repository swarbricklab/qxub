"""
This script configures the installation of the 'qconda' Python package using setuptools.
Defines the package metadata, dependencies, and entry points for the command-line interface (CLI).
The CLI command 'qt' is linked to the 'conda_submit' function, enabling submission of jobs with 
various options such as resource requests, project names, and conda environments.

Run 'pip install -e .' to install the package in editable mode for development purposes.
"""

from setuptools import setup

setup(
    name="qconda",
    version="0.1",
    py_modules=["qconda"],  # Name of the Python file without the extension
    install_requires=["click"],
    package_data={
        '': ['jobscripts/qconda.pbs'],
    },
    entry_points={
        'console_scripts': [
            'qconda = qconda:conda_submit',
        ],
    },
)
