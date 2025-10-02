# qsub_tools

[![Version](https://img.shields.io/badge/version-0.3.0-blue.svg)](https://github.com/swarbricklab/qsub_tools)

This package provides tools for running `qsub` jobs on PBS Pro in particular environments.
This avoids the boilerplate code associated with activating environments and switching work directories that clutters up many jobscripts.
In simple cases, the need to create a jobscript can be eliminated entirely.

## Installation

```bash
# Clone the repository
git clone https://github.com/swarbricklab/qsub_tools.git
cd qsub_tools

# Install in development mode
pip install -e .
```

## Quick Start

```bash
# Run a Python script in a conda environment
qxub conda --env myenv -- python script.py

# Run a command with environment modules
qxub module --mod bcftools --mod samtools -- bcftools --version

# Run a command in a Singularity container
qxub sing --sif /path/to/container.sif -- python analysis.py

# Use dry-run to see what would be executed
qxub --dry-run conda --env myenv -- python script.py
```

## Commands

The main CLI command provided by this tool is `qxub` (short for "extended qsub").
The `qxub` command is a wrapper around `qsub` and so accepts many of the same options, such as `-l mem=16GB`.

> The intention is to support all `qsub` options, but initially only the most common options have been implemented.
> Create a feature request if your favourite option is missing.

In addition to standard `qsub` options, `qxub` also has an `--execdir` option for controlling the work directory where the job is executed.
Unlike `qsub`, `qxub` defaults to executing in the current work directory.

### Common Options

- `--dry-run`: Show the `qsub` command that would be executed without actually submitting it
- `--quiet`: Submit job and exit immediately without monitoring output
- `--execdir`: Set the working directory for job execution (default: current directory)
- `--out`, `--err`: Specify custom paths for STDOUT and STDERR logs

### Output

The STDOUT and STDERR produced by `qxub` jobs is redirected to the following logs:
- `$TMPDIR/qt/timestamp/out`
- `$TMPDIR/qt/timestamp/err`

The location of these logs can be specified via the `--out` and `--err` options.
However, unless `--quiet` is specified `qxub` will monitor these fails and stream their context back to STDOUT and STDERR in the shell where `qxub` was launched.
This makes it possible to view the output from the job in real time, without having to look up the logs (although the logs are still saved).
This also facilitates redirection of `qxub` output via traditional `>` and `2>` operators.

> Note: except for `--quiet` mode, `qxub` does not terminate until the job that it initiates terminates.
> But `qxub` only polls `qstat` every 30 seconds, and so `qxub` may appear to "hang" for up to 30 seconds even after the command has finished.
> In most cases, the time wasted waiting for the next poll will be less than the time required to look up the logs manually.

## Subcommands

`qxub` provides the following subcommands for executing commands or scripts in particular environments:
- `sub`: no special environment, but work directory can controlled and argumemts can be passed into scripts without special gymnastics
- `conda`: execute command or script in the specified conda environment
- `module`: activate specified environment modules before executing command or script
- `sing`: execute command or script in the specified Singularity container
- `config`: view or edit config settings (TODO)

### conda subcommand

Execute commands in a conda environment:

```bash
# Basic usage
qxub conda --env myenv -- python script.py

# With resource specifications
qxub --resources mem=16GB conda --env myenv -- python large_analysis.py

# With pre/post commands
qxub conda --env myenv --pre "echo Starting analysis" --post "echo Analysis complete" -- python script.py
```

### module subcommand

Execute commands with environment modules loaded:

```bash
# Single module
qxub module --mod bcftools -- bcftools --version

# Multiple modules
qxub module --mod bcftools --mod samtools --mod python3 -- python analysis.py

# With pre/post commands
qxub module --mod cuda --pre "nvidia-smi" -- python gpu_training.py
```

### sing subcommand

Execute commands in Singularity containers:

```bash
# Basic usage
qxub sing --sif /path/to/container.sif -- python script.py

# With Singularity options (before --)
qxub sing --sif container.sif --bind /data:/data --env DEBUG=1 -- python script.py

# With pre/post commands
qxub sing --sif container.sif --pre "module load singularity" -- analysis_tool --input data.txt
```

## Signal Handling

All subcommands support graceful interruption with Ctrl+C:
- Automatically deletes submitted jobs with `qdel`
- Provides clear feedback during cleanup
- Ensures no orphaned jobs remain in the queue

## Changelog

### Version 0.3.0
- Added `module` subcommand for environment module support
- Added `sing` subcommand for Singularity container support  
- Added `--pre` and `--post` options for command chaining
- Added Ctrl+C signal handling with automatic job cleanup
- Improved spinner display and output coordination
- Added comprehensive argument passthrough for Singularity options

### Version 0.2.0
- Initial release with `conda` subcommand
- Basic PBS job submission and monitoring functionality
