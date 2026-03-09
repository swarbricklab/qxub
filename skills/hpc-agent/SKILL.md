---
description: Run commands on HPC compute nodes via PBS job scheduler using qxub
---

# HPC Job Execution with qxub

`qxub` submits PBS jobs to HPC compute nodes. Use it when commands need more resources than login nodes provide (heavy computation, large memory, GPUs).
Unlike `qsub`, `qxub` streams output to the current terminal and waits for job completion, so that you can interactively monitor progress and debug issues without SSHing into compute nodes.

**Always activate the environment first:**
```bash
conda activate qxub
```

## Quick Reference

| Command | Alias | Purpose |
|---------|-------|---------|
| `qxub exec` | `qx` | Submit a job and wait for completion |
| `qxub history` | — | View past jobs and their output |
| `qxub status` | `qxtat` | Check job status |

## Running Commands (`qxub exec` / `qx`)

**Syntax:** Options before `--`, command after.

```bash
# Basic: run in conda environment
qx --env myenv -- python script.py

# With resources
qx --env pytorch --mem 32GB --cpus 8 --runtime 2h -- python train.py

# With GPU
qx --env pytorch --mem 32GB --cpus 12 --runtime 4h --resources ngpus=1 -- python train.py

# Using environment modules instead of conda
qx --mod python3 --mod gcc -- python analysis.py

# Preview without submitting
qx --dry --env myenv -- python script.py
```

### Key Options

| Option | Description | Example |
|--------|-------------|---------|
| `--env` | Conda environment | `--env pytorch` |
| `--mod` | Environment module | `--mod python3` |
| `--sif` | Singularity image | `--sif myimage.sif` |
| `--mem` | Memory | `--mem 16GB` |
| `--cpus` | CPU cores | `--cpus 8` |
| `--runtime` | Walltime limit | `--runtime 2h30m` |
| `--queue` | Queue (use `auto` for optimal) | `--queue auto` |
| `--dry` | Preview job script | `--dry` |
| `--quiet` | Suppress progress output | `--quiet` |
| `--terse` | Print job ID only, return immediately | `--terse` |

Run `qxub exec --help` for all options.

## Debugging Jobs (`qxub history`)

After a job completes (or fails), inspect its output:

```bash
# Most recent job info
qxub history latest

# View stdout of most recent job
qxub history out

# View stderr of most recent job
qxub history err

# View PBS log (scheduler info) of most recent job
qxub history log

# Specific job by ID
qxub history out 12345678
qxub history err 12345678
```

Run `qxub history --help` for all subcommands.

## Checking Job Status (`qxub status` / `qxtat`)

```bash
# List recent jobs
qxub status list

# Show specific job
qxub status show 12345678

# Summary of job counts by status
qxub status summary
```

## Typical Agent Workflow

1. **Submit job:**
   ```bash
   qx --env myenv --mem 16GB --cpus 4 --runtime 1h -- python script.py
   ```

2. **If job fails, check errors:**
   ```bash
   qxub history err
   ```

3. **Check stdout for results:**
   ```bash
   qxub history out
   ```

4. **For long jobs, check status:**
   ```bash
   qxub status list
   ```

## Further Documentation

- [docs/examples.md](../../docs/examples.md) — Common usage patterns
- [docs/configuration.md](../../docs/configuration.md) — Config system
- [docs/option-placement.md](../../docs/option-placement.md) — CLI reference
