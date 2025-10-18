# Aliases

Create shortcuts for common workflows.

## Quick Start

```bash
# Create an alias
qxub config alias set gpu --env pytorch --queue gpu -l ngpus=1 -l ncpus=12

# Use the alias
qxub alias gpu -- python train.py

# List all aliases
qxub config alias list
```

## Commands

```bash
qxub config alias set NAME [OPTIONS]    # Create alias
qxub config alias list                  # List aliases
qxub config alias show NAME             # Show alias details
qxub config alias delete NAME           # Delete alias

qxub alias NAME [-- COMMAND]            # Execute alias
```

## Examples

```bash
# GPU training workflow
qxub config alias set train --env pytorch --queue gpu -l ngpus=1 -l ncpus=12 --name "training-{timestamp}"

# Big memory analysis
qxub config alias set bigmem --env pandas -l mem=256GB --queue hugemem

# Quick testing
qxub config alias set test --env base --queue express -l walltime=00:10:00

# Use them
qxub alias train -- python train.py --epochs 100
qxub alias bigmem -- python analyze.py input.csv
qxub alias test -- python test_script.py
```

## Command Override

```bash
# Create alias without command
qxub config alias set gpu --env pytorch --queue gpu -l ngpus=1

# Use with different commands
qxub alias gpu -- python train.py
qxub alias gpu -- python evaluate.py
qxub alias gpu -- jupyter notebook
```

## Option Override

```bash
# Override alias options
qxub alias gpu --queue normal -- python small_job.py  # Use normal queue instead
qxub alias train --name custom-job -- python train.py  # Custom job name
```
