# qxub

qxub is a sophisticated PBS job submission wrapper designed for HPC environments. It eliminates the boilerplate and complexity of writing PBS job scripts by providing a unified command-line interface for submitting jobs across different execution contexts. Whether you need to run code in conda environments, with environment modules, inside containers, or as direct submissions, qxub handles the PBS integration seamlessly while offering intelligent queue selection, resource management, and comprehensive configuration options.

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

# Built-in command aliases for faster typing
qx dvc push                                      # Short for 'qxub exec -- dvc push'
qxtat 12345.gadi-pbs                             # Short for 'qxub status'
qxet "command" --env base                        # Short for 'qxub config shortcut set "command" --env base'

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
| `--quiet` | Silent monitoring (no progress messages or job ID) | `--quiet` |
| `--terse` | Emit job ID then silent monitoring | `--terse` |
| `--remote` | Execute on remote HPC system via SSH | `--remote gadi` |
| `--dry` | Preview without submission | `--dry` |
| `-v/-vv/-vvv` | Verbosity levels | `-vv` |

## Shortcuts System
```bash
# List available shortcuts
qxub config shortcut list

# Show shortcut source information (system vs user)
qxub config shortcut list --show-origin

# Create user shortcuts (default)
qxub config shortcut set "python" --env base --description "Python with base conda environment"
qxub config shortcut set "gcc" --mod gcc --description "GCC compiler with modules"

# Create system-wide shortcuts (admin/team use)
qxub config shortcut set "dvc data status" --system --env dvc3 --resources mem=64GB,ncpus=16

# Use shortcuts (automatic detection)
qxub exec -- python script.py      # Detects 'python' shortcut
qxub exec -- gcc --version         # Detects 'gcc' shortcut

# Explicit shortcut usage
qxub exec --shortcut python -- script.py

# Show shortcut details with source info
qxub config shortcut show python

# Delete shortcuts
qxub config shortcut delete python              # Delete user shortcut
qxub config shortcut delete "system-wide" --system --yes  # Delete system shortcut

# Rename shortcuts
qxub config shortcut rename "old-name" "new-name"
qxub config shortcut rename "old-name" "new-name" --system  # Rename system shortcut

# Built-in command aliases for faster access
qx --env myenv -- python script.py    # Equivalent to 'qxub exec'
qxtat 12345.gadi-pbs                   # Equivalent to 'qxub status'
qxet "ml-pipeline" --env pytorch       # Equivalent to 'qxub config shortcut set'
```
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

# Create user shortcuts for personal workflows
qxub config shortcut set gpu --env pytorch --description "GPU training with PyTorch"
qxub exec -- gpu train.py    # Use the shortcut

# Create system-wide shortcuts for team use (requires admin permissions)
qxub config shortcut set "dvc repro" --system --env dvc3 --resources mem=32GB --description "DVC data pipeline"

# View configuration and shortcut files
qxub config files
qxub config shortcut files    # Shows both system and user shortcut file locations
qxub config get defaults
```

## Help
```bash
qxub --help                       # General help
qxub exec --help                 # Execution command help
qxub config --help               # Configuration help
qxub config shortcut --help      # Shortcut management help
qxub history --help              # History management help
qxub resources --help            # Resource tracking help
qxub monitor --help              # Monitor multiple jobs

# Built-in alias help
qx --help                        # qxub exec help
qxtat --help                     # qxub status help
qxet --help                      # qxub config shortcut set help
```
