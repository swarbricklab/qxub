"""
This script configures the installation of the 'qxub' Python package using setuptools.
Defines the package metadata, dependencies, and entry points for the command-line interface (CLI).
The CLI command 'qxub' is linked to the 'cli.qxub' function, enabling submission of jobs with 
various options such as resource requests, project names, and conda environments.

Run 'pip install -e .' to install the package in editable mode for development purposes.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="qxub",
    version="0.3.0",
    author="John Reeves",
    author_email="j.reeves@garvan.org.au",
    description="Simplified job submission to HPC",
    long_description = long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/swarbricklab/qsub_tools",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'qxub': ['jobscripts/*.pbs'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'omegaconf',
        'click',
        'tailer',
        'rich'
    ],
    entry_points={
        'console_scripts': ['qxub=qxub.cli:qxub'],
    },
)
