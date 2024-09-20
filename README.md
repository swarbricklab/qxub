# qsub_tools

This package provides tools for running `qsub` jobs on PBS Pro in particular environments.
This avoids the boilerplate code associated with activating environments and switching work directories that clutters up many jobscripts.
In simple cases, the need to create a jobscript can be eliminated entirely.

## Commands

The main CLI command provided by this tool is `qt` (short for "qsub tools").
The `qt` command is a wrapper around `qsub` and so accepts many of the same options, such as `-l mem=16GB`.

> The intention is to support all `qsub` options, but initially only the most common options have been implemented.
> Create a feature request if your favourite option is missing.

In addition to standard `qsub` options, `qt` also has an `--execdir` option for controlling the work directory where the job is executed.
Unlike `qsub`, `qt` defaults to executing in the current work directory.

Aspirational: unlike `qsub`, `qt` will stream to STDOUT and STDERR, thereby facilitating traditional redirection.

## Subcommands

`qt` provides the following subcommands for executing commands or scripts in particular environments:
- `sub`: no special environment, but work directory can controlled and argumemts can be passed into scripts without special gymnastics
- `conda`: execute command or script in the specified conda environment
- `mod`: activate specified environment modules before executing command or script
- `sing`: execute command or script in the specified singularity environment
- `config`: view or edit config settings

 
