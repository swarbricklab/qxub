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

## Complex Commands with Variables

For commands with special characters, quotes, or variables, use `--cmd`:

```bash
# Submission-time variables (expanded when you submit)
qxub --env base --cmd "python script.py --input ${HOME}/data.txt --user ${USER}"

# Execution-time variables (expanded when job runs)
qxub --env base --cmd 'echo "Job ${{PBS_JOBID}} running on node ${{HOSTNAME}}"'

# Mixed variables
qxub --env base --cmd 'python analyze.py --input ${HOME}/data --output ${{TMPDIR}}/results-${USER}/'

# Complex Python commands with quotes
qxub --env base --cmd 'python -c "import sys; print(f\"Args: {sys.argv[1:]}\"); print(\"Job: ${{PBS_JOBID}}\")"'

# AWK/shell commands (literal $ preserved)
qxub --env base --cmd 'awk "{print \$1, \$2}" data.txt | head -10'
```

### Variable Substitution Rules

- `${var}` → Expanded at submission time (your environment)
- `${{var}}` → Converted to `${var}` for execution time (job environment)
- `$other` → Preserved unchanged (literal usage)

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

# See variable expansion
qxub --dry --env myenv --cmd 'echo "User ${USER} job ${{PBS_JOBID}}"'
```
