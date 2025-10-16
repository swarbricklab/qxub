# qxub Tutorial: Complete Guide for Newcomers

Welcome to `qxub`, the modern PBS job submission wrapper that eliminates boilerplate and brings real-time monitoring to HPC environments. This tutorial will take you from basic usage to advanced parallel execution patterns.

## What is qxub?

`qxub` is a sophisticated wrapper around PBS that:
- **Eliminates boilerplate**: No more writing PBS scripts for simple jobs
- **Provides real-time output**: Watch your job output as it runs, just like a local command
- **Intelligent queue selection**: Automatically chooses the best queue for your resources
- **Unified CLI**: All options before `--`, command after - simple and consistent
- **Rich execution contexts**: Run in conda environments, with modules, or in containers
- **Built-in monitoring**: Track multiple jobs and get efficiency reports

## Prerequisites

This tutorial assumes:
- You have access to NCI Gadi and are familiar with basic HPC concepts
- You have activated the `dvc3` conda environment where `qxub` is installed
- The system configuration provides sensible defaults for project `a56`

```bash
# Activate the environment if you haven't already
conda activate dvc3

# Verify qxub is available
qxub --version
```

## Quick Configuration Check

Before starting, let's see what's already configured:

```bash
# See all current configuration with origins
qxub config list --show-origin
```

You should see defaults for project `a56`, standard resource allocations, and several ready-made aliases that we'll explore throughout this tutorial.

## Tutorial Structure

This tutorial is organized into progressive sections:

### Part 1: Foundation
1. **[Basic Usage](01-basic-usage.md)** - Your first qxub commands with real-time output
2. **[Resources and Queues](02-resources-and-queues.md)** - Specifying resources and automatic queue selection
3. **[Debugging and Verbosity](03-debugging-and-verbosity.md)** - Essential debugging tools and self-help

### Part 2: Execution Environments
4. **[Execution Contexts](04-execution-contexts.md)** - Using `--env`, `--mod`, `--mods`, and `--sif`
5. **[Complex Commands](05-complex-commands.md)** - Handling quotes, variables, and the `--cmd` option

### Part 3: Workflow Management
6. **[Job History](06-history.md)** - Tracking and analyzing your job history
7. **[Aliases](07-aliases.md)** - Using and creating command shortcuts
8. **[Configuration](08-configuration.md)** - Understanding config levels and customization

### Part 4: Advanced Usage
9. **[Parallel Execution](09-parallel-execution.md)** - Multiple job patterns and monitoring with `qxub monitor`
10. **[DVC Integration](10-dvc-integration.md)** - Using qxub in data science pipelines
11. **[Remote Execution](11-remote-execution.md)** - Running jobs from your laptop

## Learning Approach

Each section builds on previous knowledge and includes:
- **Hands-on examples** you can run immediately
- **Expected output** so you know what success looks like
- **Common pitfalls** and how to avoid them
- **Best practices** learned from real-world usage

## Getting Help

Throughout your qxub journey:
- Use `qxub --help` and `qxub COMMAND --help` for built-in documentation
- Try `--dry` to see what would happen without submitting jobs
- Use `-v` or `-vv` for more detailed output when debugging
- File bug reports and feature requests - the community appreciates feedback!

## Ready to Start?

Jump into **[Basic Usage](01-basic-usage.md)** to submit your first job and experience the magic of real-time PBS output!

---

*This tutorial assumes qxub v2.3+ with the NCI Gadi platform configuration. Commands and output may vary slightly with different versions or platforms.*
