# qxub

Submit PBS jobs with conda environments, modules, or containers.

## Install
```bash
git clone https://github.com/swarbricklab/qsub_tools.git
cd qsub_tools
pip install -e .
```

## Usage
```bash
# Basic syntax: qxub [options] -- command
qxub --env myenv -- python script.py        # Conda environment
qxub --mod python3 -- python script.py     # Environment module
qxub --sif container.sif -- python script.py  # Singularity container
qxub -- python script.py                    # Direct submission

# Add PBS options before --
qxub -l mem=16GB --queue normal --env myenv -- python script.py

# Cost-optimized auto queue selection (recommended!)
qxub --queue auto -l mem=1200GB --env myenv -- python big_job.py    # → megamem (58% cheaper!)
qxub --queue auto -l ngpus=1 --env pytorch -- python train.py      # → gpuvolta
qxub --queue auto -l ncpus=5000 --env myenv -- python parallel.py  # → normalsr

# Preview without running
qxub --dry --env myenv -- python script.py
```

## Key Rules
- **All qxub options go before `--`**
- **Your command goes after `--`**
- **Use `--dry` to preview**

## Options
| Option | Description | Example |
|--------|-------------|---------|
| `--env` | Conda environment | `--env pytorch` |
| `--mod` | Environment module | `--mod python3` |
| `--sif` | Singularity container | `--sif container.sif` |
| `-l` | PBS resources | `-l mem=16GB,ncpus=8` |
| `--queue` | PBS queue (use `auto` for cost optimization!) | `--queue auto` |
| `--name` | Job name | `--name myjob` |
| `--terse` | Terse output: emit job ID only | `--terse` |

## Parallel Job Execution (v2.3)
```bash
# Submit multiple jobs and monitor them
find -name "*.csv" -exec qxub --terse --env myenv -- process.py {} \; | qxub monitor

# Monitor with live status display (emoji indicators)
qxub monitor 12345.gadi-pbs 12346.gadi-pbs 12347.gadi-pbs

# Summary mode for many jobs
echo -e "job1.gadi-pbs\njob2.gadi-pbs\njob3.gadi-pbs" | qxub monitor --summary

# Clean display without PBS suffix
qxub monitor --suffix .gadi-pbs 12345.gadi-pbs 12346.gadi-pbs

# Quiet monitoring for scripting
echo "12345.gadi-pbs" | qxub monitor --quiet
```

## Configuration (Optional)
```bash
# Set personal defaults to avoid repetition
qxub config set defaults.project "a56"
qxub config set --global defaults.queue "normal"  # explicit user config

# Set monitor defaults
qxub config set monitor.default_suffix ".gadi-pbs"
qxub config set monitor.default_interval 15

# Project-level configuration (v2.3)
qxub config init-project      # Create .qx/ directory
qxub config set --project qxub.defaults.walltime "2:00:00"  # Team defaults (git-tracked)
qxub config set --local qxub.defaults.queue "express"       # Personal overrides (git-ignored)
qxub config set --test qxub.defaults.queue "copyq"          # CI settings (git-tracked)

# Create shortcuts
qxub config alias set gpu --env pytorch --queue gpu -- python train.py
qxub alias gpu  # Run the alias
```

## Help
```bash
qxub --help                    # General help
qxub config --help            # Configuration help
qxub monitor --help           # Monitor multiple jobs
```
