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
# Basic syntax: qxub exec [options] -- command
qxub exec --env myenv -- python script.py        # Conda environment
qxub exec --mod python3 -- python script.py     # Environment module
qxub exec --sif container.sif -- python script.py  # Singularity container
qxub exec -- python script.py                    # Direct submission

# Add PBS options before --
qxub exec --resources mem=16GB --queue normal --env myenv -- python script.py

# Cost-optimized auto queue selection (recommended!)
qxub exec --queue auto --resources mem=1200GB --env myenv -- python big_job.py    # → megamem (58% cheaper!)
qxub exec --queue auto --resources ngpus=1 --env pytorch -- python train.py      # → gpuvolta
qxub exec --queue auto --resources ncpus=5000 --env myenv -- python parallel.py  # → normalsr

# Preview without running
qxub exec --dry --env myenv -- python script.py

# Shortcuts - automatic command detection
qxub exec -- python script.py     # Auto-detects 'python' shortcut (conda: base)
qxub exec -- gcc --version        # Auto-detects 'gcc' shortcut (modules: gcc)

# Explicit shortcut usage
qxub exec --shortcut python -- script.py
```

## Key Rules
- **All qxub exec options go before `--`**
- **Your command goes after `--`**
- **Use `--dry` to preview**
- **Commands starting with known prefixes trigger shortcuts automatically**

## Options
| Option | Description | Example |
|--------|-------------|---------|
| `--env` | Conda environment | `--env pytorch` |
| `--mod` | Environment module | `--mod python3` |
| `--sif` | Singularity container | `--sif container.sif` |
| `--shortcut` | Use predefined shortcut | `--shortcut python` |
| `--cmd` | Complex command (alternative to `--`) | `--cmd "echo \"Hello ${USER}\""` |
| `--resources` | PBS resources | `--resources mem=16GB --resources ncpus=8` |
| `--queue` | PBS queue (use `auto` for cost optimization!) | `--queue auto` |
| `--name` | Job name | `--name myjob` |
| `--project` | PBS project | `--project a56` |
| `--out` | Output file path | `--out /scratch/job.out` |
| `--err` | Error file path | `--err /scratch/job.err` |
| `--quiet` | Suppress progress messages | `--quiet` |
| `--terse` | Terse output: emit job ID only | `--terse` |
| `--remote` | Execute on remote HPC system via SSH | `--remote gadi` |
| `--dry` | Preview without submission | `--dry` |
| `-v/-vv/-vvv` | Verbosity levels | `-vv` |

## Shortcuts System
```bash
# List available shortcuts
qxub config shortcut list

# Create a new shortcut
qxub config shortcut set "python" --env base --description "Python with base conda environment"
qxub config shortcut set "gcc" --mod gcc --description "GCC compiler with modules"

# Use shortcuts (automatic detection)
qxub exec -- python script.py      # Detects 'python' shortcut
qxub exec -- gcc --version         # Detects 'gcc' shortcut

# Explicit shortcut usage
qxub exec --shortcut python -- script.py

# Show shortcut details
qxub config shortcut show python

# Delete a shortcut
qxub config shortcut delete python
```

## Parallel Job Execution & Monitoring
```bash
# Submit multiple jobs and monitor them
find -name "*.csv" -exec qxub exec --terse --env myenv -- process.py {} \; | qxub monitor

# Monitor with live status display (emoji indicators)
qxub monitor 12345.gadi-pbs 12346.gadi-pbs 12347.gadi-pbs

# Summary mode for many jobs
echo -e "job1.gadi-pbs\njob2.gadi-pbs\njob3.gadi-pbs" | qxub monitor --summary

# Clean display without PBS suffix
qxub monitor --suffix .gadi-pbs 12345.gadi-pbs 12346.gadi-pbs

# Quiet monitoring for scripting
echo "12345.gadi-pbs" | qxub monitor --quiet
```

## Job Management
```bash
# View execution history
qxub history executions --limit 10
qxub history latest

# View resource efficiency
qxub resources list --limit 5
qxub resources stats

# Configuration management
qxub config get defaults.project
qxub config set defaults.project "a56"
qxub config files
```

## Configuration (Optional)
```bash
# Set personal defaults to avoid repetition
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# Create shortcuts for common workflows
qxub config shortcut set gpu --env pytorch --description "GPU training with PyTorch"
qxub exec -- gpu train.py    # Use the shortcut

# View configuration
qxub config files
qxub config get defaults
```

## Help
```bash
qxub --help                   # General help
qxub exec --help             # Execution command help
qxub config --help           # Configuration help
qxub history --help          # History management help
qxub resources --help        # Resource tracking help
qxub monitor --help          # Monitor multiple jobs
```
