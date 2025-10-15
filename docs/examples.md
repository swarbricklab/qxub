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
qxub --env pytorch -- python train.py
```

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

### Smart Quote Processing (Recommended)

Use double quotes around `--cmd` for automatic quote handling:

```bash
# Clean, readable syntax with smart quotes
qxub --env base --cmd "find /data -name \"*.txt\" -exec echo \"Found: {}\" \;"

# Complex shell commands with nested quotes
qxub --env base --cmd "sh -c \"echo \\\"Processing file: \$1\\\"\" arg"

# Mixed variables: submission-time and execution-time
qxub --env base --cmd "echo \"User ${USER} running job ${{PBS_JOBID}}\""

# AWK commands with field references preserved
qxub --env base --cmd "awk \"{print \\\$1, \\\$2}\" data.txt"

# JSON-like output with complex escaping
qxub --env base --cmd "echo \"{\\\"user\\\": \\\"${USER}\\\", \\\"cost\\\": \\\"\$50\\\"}\""
```

### Traditional Syntax (Backward Compatible)

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
qxub config alias set train --env pytorch --queue gpu -l ngpus=1 -l ncpus=12
qxub config alias set analyze --env pandas -l mem=64GB

# Use them
qxub alias train -- python train.py
qxub alias analyze -- python analysis.py
```

## Parallel Job Execution (v2.3)

Submit multiple jobs and monitor them as a group:

```bash
# Process multiple files in parallel
find data/ -name "*.csv" -exec qxub --terse --env pandas -- python process.py {} \; | qxub monitor

# Submit batch of analysis jobs
for sample in sample_*.fastq; do
    qxub --terse --env biotools -- analyze.py $sample
done | qxub monitor --quiet

# Monitor specific job IDs
qxub monitor 12345.gadi-pbs 12346.gadi-pbs 12347.gadi-pbs

# DVC pipeline stage with parallel execution
find input/ -name "*.txt" -exec qxub --terse --env nlp -- process.py {} \; | qxub monitor
# This blocks until all jobs complete, then DVC can verify outputs exist
```

### Terse Output

The `--terse` flag emits only the job ID for pipeline use:

```bash
# Normal output
qxub --env myenv -- python script.py
# ✅ Job submitted successfully! Job ID: 12345.gadi-pbs
# [monitoring and output continues...]

# Terse output
qxub --terse --env myenv -- python script.py
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

### Real-World Workflows

#### Bioinformatics Pipeline

```bash
# Process samples in parallel
for sample in samples/*.fastq.gz; do
    qxub --terse --env biotools -l mem=32GB -- \
        fastqc "$sample" --outdir results/
done | qxub monitor --suffix .gadi-pbs --name-only

# Summary mode for large cohorts
ls cohort/*.bam | xargs -I {} qxub --terse --env gatk -- \
    gatk HaplotypeCaller -I {} -O {}.vcf | \
    qxub monitor --summary --interval 120
```

#### Machine Learning Grid Search

```bash
# Hyperparameter sweep
for lr in 0.01 0.001 0.0001; do
    for batch in 32 64 128; do
        qxub --terse --env pytorch --queue gpu -l ngpus=1 --cmd \
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
find raw_data/ -name "*.csv" -exec qxub --terse --env pandas -- preprocess.py {} \; | \
    qxub monitor --summary --quiet > stage1_complete.flag

# Stage 2: Analysis (fewer, larger jobs)
if [ -f stage1_complete.flag ]; then
    ls processed_data/*.parquet | xargs -I {} qxub --terse --env analytics \
        -l mem=64GB -l ncpus=16 -- analyze.py {} | \
        qxub monitor --name-only
fi
```

## Dry Run

```bash
# Preview what qxub would do
qxub --dry --env myenv -- python script.py

# See variable expansion
qxub --dry --env myenv --cmd 'echo "User ${USER} job ${{PBS_JOBID}}"'

# Terse dry run
qxub --terse --dry --env myenv -- python script.py
# DRY_RUN
```
