# Examples

## Basic Usage

```bash
# Conda environment
qxub --env pytorch -- python train.py

# Environment modules
qxub --mod python3 --mod gcc -- python analysis.py

# Singularity container
qxub --sif /containers/blast.sif -- blastn -query input.fa -db nt

# Direct submission (no environment)
qxub -- ./my_program
```

## With PBS Options

```bash
# Memory and CPU
qxub -l mem=32GB -l ncpus=8 --env myenv -- python script.py

# GPU job
qxub -l ngpus=1 -l ncpus=12 --queue gpu --env pytorch -- python train.py

# Long-running job
qxub -l walltime=24:00:00 --env myenv -- python long_script.py
```

## Configuration

```bash
# Set defaults once
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# Now just use
qxub --env myenv -- python script.py
```

## Aliases

```bash
# Create shortcuts for common workflows
qxub config alias set train --env pytorch --queue gpu -l ngpus=1 -l ncpus=12
qxub config alias set analyze --env pandas -l mem=64GB

# Use them
qxub alias train -- python train.py
qxub alias analyze -- python analysis.py
```

## Dry Run

```bash
# Preview what qxub would do
qxub --dry --env myenv -- python script.py
```
