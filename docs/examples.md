# Examples

## Project Setup

```bash
# Initialize project configuration
qxub config init-project

# Set team defaults (tracked in git)
qxub config set --project qxub.defaults.walltime "4:00:00"
qxub config set --project qxub.defaults.queue "normal"

# Personal overrides (ignored by git)
qxub config set --local qxub.defaults.queue "express"

# Now jobs use project defaults automatically
qxub exec --env pytorch -- python train.py
```

## Shortcuts System

```bash
# List available shortcuts
qxub config shortcut list

# Show shortcuts with source information (system vs user)
qxub config shortcut list --show-origin

# Create user shortcuts for personal workflows
qxub config shortcut set "python" --env base --description "Python with base conda environment"
qxub config shortcut set "pytorch" --env pytorch --queue gpuvolta --resources ngpus=1 --description "PyTorch GPU training"
qxub config shortcut set "gcc" --mod gcc --description "GCC compiler with modules"

# Create system-wide shortcuts for team use (requires admin permissions)
qxub config shortcut set "dvc data status" --system --env dvc3 --resources mem=64GB,ncpus=16 --description "DVC data pipeline status check"

# Use shortcuts (automatic detection by command prefix)
qxub exec -- python script.py               # Uses 'python' shortcut
qxub exec -- pytorch train.py               # Uses 'pytorch' shortcut
qxub exec -- gcc -o program source.c        # Uses 'gcc' shortcut

# Explicit shortcut usage
qxub exec --shortcut python -- script.py

# Show shortcut details with source information
qxub config shortcut show python

# Delete shortcuts
qxub config shortcut delete "old-shortcut"
qxub config shortcut delete "system-shortcut" --system --yes   # Delete system shortcut without prompt

# Rename shortcuts
qxub config shortcut rename "old-name" "new-name"
qxub config shortcut rename "old-name" "new-name" --system     # Rename system shortcut

# Built-in command aliases for faster workflow
qx --env pytorch -- train.py                # Short for 'qxub exec'
qxet "ml-training" --env pytorch --resources ngpus=1  # Short for 'qxub config shortcut set'

# Override shortcut settings
qxub exec --shortcut python --name custom-job -- special_script.py
```

## Workflow-Friendly Resource Options

qxub v3.1.0 introduces convenient resource flags that automatically convert to PBS format:

```bash
# New workflow-friendly syntax (recommended for most users)
qxub exec --mem 16GB --cpus 8 --runtime 2h30m --disk 100GB --volumes gdata/a56 --env myenv -- python script.py

# Equivalent traditional PBS syntax
qxub exec --resources mem=16GB,ncpus=8,walltime=2:30:00,jobfs=100GB,storage=gdata/a56 --env myenv -- python script.py

# Mix and match approaches
qxub exec --mem 32GB --cpus 16 --volumes gdata/a56+scratch/a56 --env myenv -- python large_job.py

# Alternative flag names (same functionality)
qxub exec --memory 8GB --threads 4 --time 1h --jobfs 50GB --storage gdata/a56+gdata/px14 --env base -- python analysis.py

# Flexible time formats
qxub exec --runtime 2h30m --env base -- python script.py     # 2 hours 30 minutes
qxub exec --runtime 90m --env base -- python script.py       # 90 minutes
qxub exec --runtime 1:30:00 --env base -- python script.py   # HH:MM:SS format

# Memory format variations
qxub exec --mem 8GB --env base -- python script.py           # Gigabytes
qxub exec --mem 4096MB --env base -- python script.py        # Megabytes

# Quick examples for common scenarios
qxub exec --mem 4GB --cpus 2 --time 30m --env base -- python quick_analysis.py
qxub exec --mem 64GB --cpus 16 --time 4h --disk 200GB --volumes gdata/a56 --env pytorch -- python training.py
qxub exec --mem 1TB --cpus 48 --time 12h --volumes gdata/a56+scratch/a56 --env bioinformatics -- ./genome_assembly.sh

# GPU jobs with workflow-friendly resources
qxub exec --mem 32GB --cpus 12 --time 8h --volumes gdata/a56+gdata/px14 --resources ngpus=1 --env pytorch -- python train_model.py

# Perfect for Snakemake and other workflow engines
snakemake --cluster "qxub exec --mem {resources.mem_gb}GB --cpus {threads} --time {resources.runtime}h --volumes {resources.storage} --"

# Works with all execution contexts
qxub exec --mem 8GB --cpus 4 --mod python3 -- python script.py               # Modules
qxub exec --mem 16GB --cpus 8 --sif container.sif -- ./analysis               # Singularity
qxub exec --mem 4GB --cpus 2 -- ./native_program                              # Direct submission
```

### Config Defaults for Workflow-Friendly Options

Set default values for workflow-friendly options to avoid repeating common resource specifications:

```bash
# Set your typical resource defaults
qxub config set mem "8GB"
qxub config set cpus 4
qxub config set runtime "2h"
qxub config set disk "20GB"
qxub config set volumes "gdata/a56+scratch/a56"

# Now these are automatically applied
qxub exec --env pytorch -- python script.py
# Equivalent to: qxub exec --mem 8GB --cpus 4 --runtime 2h --disk 20GB --volumes gdata/a56+scratch/a56 --env pytorch -- python script.py

# CLI options override config defaults
qxub exec --mem 32GB --env pytorch -- python large_job.py
# Uses 32GB memory but other defaults (cpus=4, runtime=2h, etc.)

# View current config defaults
qxub config list

# Reset to remove all custom defaults
qxub config reset
```

## Basic Execution

```bash
# Conda environment
qxub exec --env pytorch -- python train.py

# Environment modules
qxub exec --mod python3 --mod gcc -- python analysis.py

# Singularity container
qxub exec --sif /containers/blast.sif -- blastn -query input.fa -db nt

# Direct submission (no environment)
qxub exec -- ./my_program

# Shortcuts (automatic detection)
qxub exec -- python train.py     # Auto-detects 'python' shortcut
qxub exec -- gcc --version       # Auto-detects 'gcc' shortcut
```

## Remote Execution (v3.3.0)

Execute jobs on remote HPC systems via SSH with platform definition delegation.

### Setup

Create a configuration file for remote execution:

```yaml
# laptop_config.yaml - for laptop → Gadi execution
defaults:
  platform: nci_gadi
  project: a56
  walltime: "1:00:00"
  mem: 4GB

platforms:
  # Delegated platform - definition resolved on remote
  nci_gadi:
    remote:
      host: gadi                          # SSH hostname from ~/.ssh/config
      working_dir: /scratch/a56/{{user}}  # {{user}} = remote username
      conda_init: |
        eval "$(conda shell.bash hook)"
        conda activate qxub
```

### Template Variables

- `{user}` - Config system username (laptop)
- `{{user}}` - Remote system username (via `$USER`)
- `{project}` - Config system project (laptop `$PROJECT`)
- `{{project}}` - Remote system project (remote `$PROJECT`)

### Basic Remote Execution

```bash
# Simple remote execution (assuming the 'nci_gadi' platform is defined in local config)
qxub exec  --platform nci_gadi --dry -- echo "hello"

# With conda environment (explicitly specifying an alternative environment)
qxub exec --config laptop_config.yaml --platform nci_gadi --env pytorch -- python train.py

# With resource specifications
qxub exec --platform nci_gadi --mem 16GB --cpus 8 -- python script.py
```

### Debugging Remote Execution

```bash
# Dry run with verbose output to see SSH command
qxub exec --platform nci_gadi --dry -vvv -- echo "test"

# Output shows the SSH command that will be executed:
# ssh gadi 'cd /scratch/a56/$USER && eval "$(conda shell.bash hook)" && conda activate qxub && qxub exec --platform nci_gadi ...'
```

### CI/Automation

```bash
# GitHub Actions / CI environments
qxub exec --config ci_config.yaml --platform nci_gadi --walltime "0:05:00" --queue copyq -- echo "CI test"

# Batch processing
for file in *.py; do
  qxub exec --config remote_config.yaml --platform nci_gadi --env myenv -- python "$file"
done
```

## With PBS Options

```bash
# Memory and CPU
qxub exec -l mem=32GB -l ncpus=8 --env myenv -- python script.py

# GPU job (auto-selects cost-effective queue)
qxub exec --queue auto -l ngpus=1 -l ncpus=12 --env pytorch -- python train.py

# Large memory job (auto-selects megamem for best cost)
qxub exec --queue auto -l mem=1500GB --env myenv -- python big_memory_job.py

# Long-running job
qxub exec -l walltime=24:00:00 --env myenv -- python long_script.py
```

## Configuration

```bash
# Set defaults once
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# Now just use
qxub exec --env myenv -- python script.py
```

## Complex Commands with Variables

For commands with special characters, quotes, or variables, use `--cmd`:

### Smart Quote Processing (Recommended)

Use double quotes around `--cmd` for automatic quote handling:

```bash
# Clean, readable syntax with smart quotes
qxub exec --env base --cmd "find /data -name \"*.txt\" -exec echo \"Found: {}\" \;"

# Complex shell commands with nested quotes
qxub exec --env base --cmd "sh -c \"echo \\\"Processing file: \$1\\\"\" arg"

# Mixed variables: submission-time and execution-time
qxub exec --env base --cmd "echo \"User ${USER} running job ${{PBS_JOBID}}\""

# AWK commands with field references preserved
qxub exec --env base --cmd "awk \"{print \\\$1, \\\$2}\" data.txt"

# JSON-like output with complex escaping
qxub exec --env base --cmd "echo \"{\\\"user\\\": \\\"${USER}\\\", \\\"cost\\\": \\\"\$50\\\"}\""
```

### Traditional Syntax (Backward Compatible)

```bash
# Submission-time variables (expanded when you submit)
qxub exec --env base --cmd "python script.py --input ${HOME}/data.txt --user ${USER}"

# Execution-time variables (expanded when job runs)
qxub exec --env base --cmd 'echo "Job ${{PBS_JOBID}} running on node ${{HOSTNAME}}"'

# Mixed variables
qxub exec --env base --cmd 'python analyze.py --input ${HOME}/data --output ${{TMPDIR}}/results-${USER}/'

# Complex Python commands with quotes
qxub exec --env base --cmd 'python -c "import sys; print(f\"Args: {sys.argv[1:]}\"); print(\"Job: ${{PBS_JOBID}}\")"'

# AWK/shell commands (literal $ preserved)
qxub exec --env base --cmd 'awk "{print \$1, \$2}" data.txt | head -10'
```

### Variable Substitution Rules

- `${var}` → Expanded at submission time (your environment)
- `${{var}}` → Converted to `${var}` for execution time (job environment)
- `$other` → Preserved unchanged (literal usage)

### Smart Quote Processing Rules

When using double quotes around `--cmd "..."`:

- `\"` → Literal `"` inside the command
- `\$` → Literal `$` (prevents variable expansion)
- `${var}` → Still expands at submission time
- `${{var}}` → Still converts for execution time
- `'text'` → Single quotes preserved literally

**Example transformations:**
```bash
Input:  --cmd "echo \"Cost: \$100, User: ${USER}\""
Output: echo "Cost: $100, User: jr9959"
```

## Aliases

```bash
# Create shortcuts for common workflows
qxub config alias set train --env pytorch --queue auto -l ngpus=1 -l ncpus=12
qxub config alias set analyze --env pandas -l mem=64GB

# Use them
qxub alias train -- python train.py
qxub alias analyze -- python analysis.py
```

## Parallel Job Execution (v2.3)

Submit multiple jobs and monitor them as a group:

```bash
# Process multiple files in parallel
find data/ -name "*.csv" -exec qxub exec --terse --env pandas -- python process.py {} \; | qxub monitor

# Submit batch of analysis jobs
for sample in sample_*.fastq; do
    qxub exec --terse --env biotools -- analyze.py $sample
done | qxub monitor --quiet

# Monitor specific job IDs
qxub monitor 12345.gadi-pbs 12346.gadi-pbs 12347.gadi-pbs

# DVC pipeline stage with parallel execution
find input/ -name "*.txt" -exec qxub exec --terse --env nlp -- process.py {} \; | qxub monitor
# This blocks until all jobs complete, then DVC can verify outputs exist
```

### Terse Output

The `--terse` flag emits only the job ID for pipeline use:

```bash
# Normal output
qxub exec --env myenv -- python script.py
# ✅ Job submitted successfully! Job ID: 12345.gadi-pbs
# [monitoring and output continues...]

# Terse output
qxub exec --terse --env myenv -- python script.py
# 12345.gadi-pbs
# [returns immediately]
```

### Monitoring Options

#### Interactive Mode (Default)

Live status table with emoji indicators:

```bash
# Interactive mode - shows live status table
qxub monitor 12345.gadi-pbs 12346.gadi-pbs

# Example output:
# ┏━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃ Job ID          ┃ Status     ┃ Job Name                     ┃
# ┡━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
# │ 12345.gadi-pbs  │ 🔄 Running │ process_data_batch1          │
# │ 12346.gadi-pbs  │ ⏳ Queued  │ process_data_batch2          │
# │ 12347.gadi-pbs  │ ✅ Success │ process_data_batch3          │
# └─────────────────┴────────────┴──────────────────────────────┘
# Next update in 28 seconds...
```

#### Summary Mode

Perfect for monitoring many jobs:

```bash
# Summary mode - compact counts
qxub monitor --summary $(cat many_job_ids.txt)

# Example output:
# 📊 Job Summary: 🔄 Running: 5  ⏳ Queued: 12  ✅ Success: 8  ❌ Failed: 1
# Next update in 25 seconds...
```

#### Display Customization

```bash
# Remove PBS suffix for cleaner display
qxub monitor --suffix .gadi-pbs 12345.gadi-pbs 12346.gadi-pbs

# Show only job names
qxub monitor --name-only 12345.gadi-pbs 12346.gadi-pbs

# Show only job IDs (compact)
qxub monitor --job-id-only 12345.gadi-pbs 12346.gadi-pbs

# Custom refresh interval
qxub monitor --interval 60 12345.gadi-pbs  # Check every minute

# Quiet mode - minimal output for scripting
qxub monitor --quiet 12345.gadi-pbs 12346.gadi-pbs
```

## Job History and Log Viewing (v3.2.0)

qxub automatically tracks execution history and provides convenient commands to view job output files.

### View Job Output Files

```bash
# View stdout from most recent job
qxub history out

# View stderr from most recent job
qxub history err

# View PBS log from most recent job
qxub history log

# View specific job output by job ID
qxub history out 153392916.gadi-pbs
qxub history err 153392916.gadi-pbs
qxub history log 153392916.gadi-pbs
```

### History Management

```bash
# List recent executions with details
qxub history executions

# Show details of the most recent execution
qxub history latest

# List computational recipes (reusable job templates)
qxub history recipes

# Clear all history (use with caution)
qxub history clear
```

### Example Output

```bash
$ qxub history out
📄 Stdout from job most recent
(/scratch/a56/user/qxub/test-job_20251026_205314.out):
This is a real job for testing logging
✅ Command completed successfully
🎉 Job completed successfully

$ qxub history executions --limit 3
                                  Recent Executions
┏━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┓
┃ Time           ┃ Recipe   ┃ Command             ┃ Status    ┃ Directory           ┃
┡━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━┩
│ 10-26 20:53:14 │ 603cb1a9 │ qxub exec --mem 1GB │ completed │ /g/data/a56/soft... │
│ 10-26 20:49:57 │ 55dd762f │ qxub exec --dry     │ completed │ /g/data/a56/soft... │
│ 10-26 20:19:49 │ c50c02f2 │ conda exec --env... │ completed │ /a56/user/project   │
└────────────────┴──────────┴─────────────────────┴───────────┴─────────────────────┘
```

**Note:** History tracking works with all execution contexts (conda, modules, containers, default) and includes both dry runs and actual job submissions.

### Real-World Workflows

#### Bioinformatics Pipeline

```bash
# Process samples in parallel
for sample in samples/*.fastq.gz; do
    qxub exec --terse --env biotools -l mem=32GB -- \
        fastqc "$sample" --outdir results/
done | qxub monitor --suffix .gadi-pbs --name-only

# Summary mode for large cohorts
ls cohort/*.bam | xargs -I {} qxub exec --terse --env gatk -- \
    gatk HaplotypeCaller -I {} -O {}.vcf | \
    qxub monitor --summary --interval 120
```

#### Machine Learning Grid Search

```bash
# Hyperparameter sweep
for lr in 0.01 0.001 0.0001; do
    for batch in 32 64 128; do
        qxub exec --terse --env pytorch --queue auto -l ngpus=1 --cmd \
            "python train.py --lr $lr --batch-size $batch --name lr${lr}_batch${batch}"
    done
done | qxub monitor --suffix .gadi-pbs --summary

# Monitor with job names only for clarity
# 📊 Job Summary: 🔄 Running: 2  ⏳ Queued: 5  ✅ Success: 2
# Next update in 28 seconds...
```

#### Data Processing Pipeline

```bash
# Stage 1: Preprocessing (many small jobs)
find raw_data/ -name "*.csv" -exec qxub exec --terse --env pandas -- preprocess.py {} \; | \
    qxub monitor --summary --quiet > stage1_complete.flag

# Stage 2: Analysis (fewer, larger jobs)
if [ -f stage1_complete.flag ]; then
    ls processed_data/*.parquet | xargs -I {} qxub exec --terse --env analytics \
        -l mem=64GB -l ncpus=16 -- analyze.py {} | \
        qxub monitor --name-only
fi
```

## Dry Run

```bash
# Preview what qxub would do
qxub exec --dry --env myenv -- python script.py

# See variable expansion
qxub exec --dry --env myenv --cmd 'echo "User ${USER} job ${{PBS_JOBID}}"'

# Terse dry run
qxub exec --terse --dry --env myenv -- python script.py
# DRY_RUN
```

## Performance and Monitoring

### File-Based Job Monitoring (v3.0+)

Starting with v3.0, qxub uses file-based monitoring for significantly improved performance:

```bash
# File-based monitoring is now the default (15x faster completion detection)
qxub exec --env myenv -- python long_job.py

# Optional: Force legacy qstat polling (not recommended)
QXUB_USE_QSTAT_MONITORING=true qxub exec --env myenv -- python job.py
```

**Performance improvements:**
- **2-second response time** for job completion detection (vs 30-second qstat polling)
- **Reduced PBS scheduler load** by eliminating frequent qstat calls
- **Background resource collection** with 60-second delay for detailed statistics
- **Real-time status files** in `<output_dir>/status/` for external monitoring

**Status file monitoring:**
```bash
# Status files created automatically in job output directory
ls /scratch/$USER/qt/*/status/
# job_started       - Job began execution
# main_started      - Main command started
# pre_exit_code     - Pre-command result
# main_exit_code    - Main command result
# post_exit_code    - Post-command result
# final_exit_code   - Final job result
```

**External monitoring example:**
```bash
# Monitor job progress externally
JOB_DIR="/scratch/a56/user/qt/20251020_194025"
while [ ! -f "$JOB_DIR/status/final_exit_code" ]; do
    if [ -f "$JOB_DIR/status/main_started" ]; then
        echo "Job is running main command..."
    fi
    sleep 5
done
echo "Job completed with exit code: $(head -n1 $JOB_DIR/status/final_exit_code)"
```
