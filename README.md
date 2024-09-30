# qsub_tools

This package provides tools for running `qsub` jobs on PBS Pro in particular environments.
This avoids the boilerplate code associated with activating environments and switching work directories that clutters up many jobscripts.
In simple cases, the need to create a jobscript can be eliminated entirely.

## Commands

The main CLI command provided by this tool is `qxub` (short for "extended qsub").
The `qxub` command is a wrapper around `qsub` and so accepts many of the same options, such as `-l mem=16GB`.

> The intention is to support all `qsub` options, but initially only the most common options have been implemented.
> Create a feature request if your favourite option is missing.

In addition to standard `qsub` options, `qxub` also has an `--execdir` option for controlling the work directory where the job is executed.
Unlike `qsub`, `qxub` defaults to executing in the current work directory.

### Output

The STDOUT and STDERR produced by `qxub` jobs is redirected to the following logs:
- `$TMPDIR/qt/timestamp/out`
- `$TMPDIR/qt/timestamp/err`

TODO: save these logs in `.qt` subdirectory of the current project.

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
- `mod`: activate specified environment modules before executing command or script (TODO)
- `sing`: execute command or script in the specified singularity environment (TODO)
- `config`: view or edit config settings (TODO)
