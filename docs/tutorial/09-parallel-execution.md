# Parallel Execution: Advanced Job Patterns and Monitoring

This section covers qxub's advanced parallel execution capabilities. You'll learn how to submit multiple jobs efficiently, monitor them collectively, and coordinate complex workflows using `--terse` mode and `qxub monitor`.

## Understanding Parallel Execution in qxub

qxub supports several patterns for parallel execution:

1. **Sequential submission** - Submit jobs one after another
2. **Batch submission** - Submit many jobs quickly with `--terse`
3. **Coordinated monitoring** - Track multiple jobs with `qxub monitor`
4. **Workflow integration** - Use in scripts and pipelines

The key is the `--terse` flag, which makes qxub output only the job ID, perfect for collecting and monitoring multiple jobs.

## The `--terse` Flag: Quiet Job Submission

### Basic Terse Usage

```bash
# Normal qxub output (verbose)
qxub exec --default -- echo "Hello World"

# Terse output (just job ID)
qxub exec --default --terse -- echo "Hello World"
```

**Terse output:**
```
12345693.gadi-pbs
```

**That's it!** Just the job ID, perfect for capturing in scripts.

### Collecting Job IDs

```bash
# Capture job ID in a variable
JOB_ID=$(qxub exec --terse -- echo "Test job")
echo "Submitted job: $JOB_ID"

# Capture multiple job IDs
JOB_IDS=()
for i in {1..5}; do
    JOB_ID=$(qxub exec --terse py -- python3 -c "
import time
print(f'Processing item {i}')
time.sleep(10)
print(f'Item {i} completed')
")
    JOB_IDS+=($JOB_ID)
    echo "Submitted job $i: $JOB_ID"
done

echo "All jobs submitted: ${JOB_IDS[@]}"
```

## Pattern 1: For Loop Parallel Submission

### Simple Parameter Sweep

```bash
# Submit jobs for different parameters
JOBS=()
for param in 0.1 0.5 1.0 2.0 5.0; do
    JOB_ID=$(qxub exec --terse --name "param-$param" py -- python3 -c "
import numpy as np
import time

param = $param
print(f'Running analysis with parameter: {param}')

# Simulate analysis
np.random.seed(42)
data = np.random.randn(1000)
result = np.mean(data) * param

time.sleep(5)  # Simulate computation time
print(f'Result with param={param}: {result:.4f}')
")
    JOBS+=($JOB_ID)
    echo "Submitted parameter $param: $JOB_ID"
done

echo "Parameter sweep jobs: ${JOBS[@]}"
```

### Processing Multiple Files

```bash
# Create some sample files to process
mkdir -p /tmp/sample_data
for i in {1..5}; do
    echo "Sample data file $i" > /tmp/sample_data/file_$i.txt
done

# Process each file in parallel
JOBS=()
for file in /tmp/sample_data/*.txt; do
    basename=$(basename "$file" .txt)
    JOB_ID=$(qxub exec --terse --name "process-$basename" py -- python3 -c "
import os
import time

filename = '$file'
basename = '$basename'

print(f'Processing {filename}')
time.sleep(3)

# Simulate file processing
with open(filename, 'r') as f:
    content = f.read()

result_file = f'/tmp/results_{basename}.txt'
with open(result_file, 'w') as f:
    f.write(f'Processed: {content}')

print(f'Completed processing {filename}')
print(f'Results saved to {result_file}')
")
    JOBS+=($JOB_ID)
    echo "Processing $basename: $JOB_ID"
done
```

## Pattern 2: While Loop with Conditions

### Processing Based on File Discovery

```bash
# Monitor directory and process new files
process_queue_file="/tmp/process_queue.txt"
echo -e "sample1.dat\nsample2.dat\nsample3.dat" > "$process_queue_file"

JOBS=()
while IFS= read -r filename; do
    [[ -z "$filename" ]] && continue  # Skip empty lines

    JOB_ID=$(qxub exec --terse --name "analyze-$filename" sc -- python3 -c "
import time
import os

filename = '$filename'
print(f'Analyzing {filename}')

# Simulate bioinformatics analysis
time.sleep(8)
print(f'Quality control for {filename}: PASS')
print(f'Analysis for {filename}: 95% confidence')
print(f'Completed analysis of {filename}')
")
    JOBS+=($JOB_ID)
    echo "Analyzing $filename: $JOB_ID"

done < "$process_queue_file"

echo "Analysis jobs: ${JOBS[@]}"
```

## Pattern 3: Find + xargs Pipeline

### Process Files Found by Criteria

```bash
# Create sample directory structure
mkdir -p /tmp/analysis/{raw,processed}
for i in {1..3}; do
    echo "Raw data $i" > "/tmp/analysis/raw/dataset_$i.csv"
done

# Find and process CSV files
find /tmp/analysis/raw -name "*.csv" -print0 | \
xargs -0 -I {} bash -c '
    basename=$(basename "{}" .csv)
    job_id=$(qxub exec --terse --name "process-$basename" py -- python3 -c "
import pandas as pd
import time
import os

input_file = \"{}\"
output_file = \"/tmp/analysis/processed/$basename.processed.csv\"

print(f\"Processing {input_file}\")
time.sleep(4)

# Simulate data processing
with open(input_file, \"r\") as f:
    content = f.read()

with open(output_file, \"w\") as f:
    f.write(f\"Processed: {content}\")

print(f\"Saved results to {output_file}\")
")
    echo "Processing $(basename "{}"): $job_id"
'
```

### Parallel Analysis with Resource Scaling

```bash
# Find large files and process with more resources
JOBS=()
find /scratch/a56/$USER -name "*.fastq" -size +100M 2>/dev/null | head -3 | \
while read fastq_file; do
    basename=$(basename "$fastq_file" .fastq)
    JOB_ID=$(qxub exec --terse --name "bigfile-$basename" \
        --resources mem=32GB,ncpus=8,walltime=4:00:00 \
        --env pysam -- python3 -c "
import time
import os

fastq_file = '$fastq_file'
basename = '$basename'

print(f'Processing large file: {fastq_file}')
print(f'Using 8 CPUs and 32GB RAM')

# Simulate intensive processing
time.sleep(15)
print(f'Quality assessment complete for {basename}')
print(f'Alignment complete for {basename}')
print(f'Variant calling complete for {basename}')
")
    JOBS+=($JOB_ID)
    echo "Large file processing $basename: $JOB_ID"
done
```

## Pattern 4: GNU Parallel Integration

### Using GNU Parallel with qxub

```bash
# Create parameter file
echo -e "0.1\n0.5\n1.0\n2.0\n5.0" > /tmp/parameters.txt

# Use GNU parallel to submit qxub jobs
parallel -j 3 'job_id=$(qxub exec --terse --name "param-{}" py -- python3 -c "
import numpy as np
import time

param = {}
print(f\"Parameter sweep: {param}\")
time.sleep(5)
result = np.random.randn(1000).mean() * param
print(f\"Result: {result:.4f}\")
"); echo "Submitted {}: $job_id"' :::: /tmp/parameters.txt
```

### Parallel with Different Resources

```bash
# Create job specifications
cat > /tmp/job_specs.txt << EOF
small,4GB,2,analysis_small.py
medium,16GB,4,analysis_medium.py
large,32GB,8,analysis_large.py
EOF

# Submit with different resource requirements
parallel --colsep ',' 'job_id=$(qxub exec --terse --name "{1}-job" \
    --resources mem={2},ncpus={3},walltime=2:00:00 py -- python analysis_{1}.py \
    ); echo "Submitted {1}: $job_id"' :::: /tmp/job_specs.txt
```

## Pattern 5: Custom R Script for Job Management

### R Script for Parallel Execution

Create an R script for sophisticated job management:

```bash
cat > /tmp/parallel_analysis.R << 'EOF'
#!/usr/bin/env Rscript
library(tibble)
library(purrr)
library(readr)

# Define analysis parameters
params <- tibble(
    sample_id = paste0("sample_", 1:10),
    replicate = rep(1:2, 5),
    condition = rep(c("control", "treatment"), each = 5)
)

# Function to submit job
submit_job <- function(sample_id, replicate, condition) {
    job_name <- paste0("analysis-", sample_id, "-rep", replicate)

    cmd <- paste0('qxub exec --terse --name "', job_name, '" ',
                  '--env tidyverse --resources mem=8GB,ncpus=2 -- Rscript -e "
                  sample_id <- \\"', sample_id, '\\"
                  replicate <- ', replicate, '
                  condition <- \\"', condition, '\\"

                  cat(\\"Processing\\", sample_id, \\"replicate\\", replicate, \\"condition\\", condition, \\"\\n\\")
                  Sys.sleep(8)

                  # Simulate analysis
                  result <- rnorm(100, mean = ifelse(condition == \\"treatment\\", 1, 0))

                  cat(\\"Mean result:\\", mean(result), \\"\\n\\")
                  cat(\\"Analysis completed for\\", sample_id, \\"\\n\\")
                  "'
    )

    job_id <- system(cmd, intern = TRUE)
    return(job_id)
}

# Submit all jobs
cat("Submitting parallel analysis jobs...\n")
job_ids <- params %>%
    pmap_chr(submit_job)

# Save job IDs
write_lines(job_ids, "/tmp/r_parallel_jobs.txt")
cat("Submitted", length(job_ids), "jobs\n")
cat("Job IDs saved to /tmp/r_parallel_jobs.txt\n")
EOF

# Run the R script
Rscript /tmp/parallel_analysis.R
```

## Monitoring Multiple Jobs with `qxub monitor`

### Basic Monitoring

After submitting multiple jobs, monitor them collectively:

```bash
# Read job IDs from file (created by previous examples)
JOB_IDS=($(cat /tmp/r_parallel_jobs.txt))

# Monitor all jobs
qxub monitor ${JOB_IDS[@]}
```

**Expected output:**
```
📊 Monitoring 10 jobs...

┌─────────────────┬──────────────────────┬─────────┬──────────┬──────────────┐
│ Job ID          │ Name                 │ Status  │ Runtime  │ Queue        │
├─────────────────┼──────────────────────┼─────────┼──────────┼──────────────┤
│ 12345701.gadi   │ analysis-sample_1... │ R       │ 00:02:15 │ normal       │
│ 12345702.gadi   │ analysis-sample_2... │ R       │ 00:02:10 │ normal       │
│ 12345703.gadi   │ analysis-sample_3... │ Q       │ 00:00:30 │ normal       │
│ 12345704.gadi   │ analysis-sample_4... │ R       │ 00:01:45 │ normal       │
│ 12345705.gadi   │ analysis-sample_5... │ Q       │ 00:00:25 │ normal       │
│ 12345706.gadi   │ analysis-sample_6... │ C       │ 00:08:12 │ normal       │ ✅
│ 12345707.gadi   │ analysis-sample_7... │ R       │ 00:03:22 │ normal       │
│ 12345708.gadi   │ analysis-sample_8... │ Q       │ 00:00:15 │ normal       │
│ 12345709.gadi   │ analysis-sample_9... │ R       │ 00:02:58 │ normal       │
│ 12345710.gadi   │ analysis-sample_10.. │ Q       │ 00:00:10 │ normal       │
└─────────────────┴──────────────────────┴─────────┴──────────┴──────────────┘

Status: 5 Running, 4 Queued, 1 Completed, 0 Failed
Progress: ████████░░ 10% (1/10 completed)

Press Ctrl+C to stop monitoring...
```

### Live Monitoring with Updates

```bash
# Monitor with live updates every 10 seconds
qxub monitor --live --refresh 10 ${JOB_IDS[@]}
```

### Summary Monitoring

```bash
# Get summary without continuous updates
qxub monitor --summary ${JOB_IDS[@]}
```

**Expected output:**
```
📈 Job Summary (10 jobs):

Status Distribution:
├── Completed: 8 jobs (80%) ✅
├── Running: 1 job (10%) 🔄
├── Queued: 1 job (10%) ⏳
└── Failed: 0 jobs (0%)

Resource Usage:
├── Total CPU Hours: 2.5 hours
├── Average Memory: 6.2GB
└── Average Walltime: 00:08:15

Efficiency:
├── Walltime Efficiency: 68% (good)
├── Memory Efficiency: 42% (acceptable)
└── Queue Distribution: normal (100%)

💡 Overall: Good resource utilization
```

## Blocking Until All Jobs Complete

### Simple Blocking Pattern

```bash
# Submit jobs and wait for completion
JOBS=()
for i in {1..5}; do
    JOB_ID=$(qxub exec --terse test -- python3 -c "
import time
print(f'Job {i} starting')
time.sleep($((RANDOM % 20 + 10)))  # 10-30 seconds
print(f'Job {i} completed')
")
    JOBS+=($JOB_ID)
done

echo "All jobs submitted: ${JOBS[@]}"

# Monitor until all complete
qxub monitor --wait-for-completion ${JOBS[@]}
echo "All jobs completed!"
```

### Advanced Blocking with Error Handling

```bash
# Submit jobs and handle failures
submit_and_monitor() {
    local job_ids=()

    # Submit jobs
    for param in 1 2 3 4 5; do
        local job_id=$(qxub exec --terse --name "job-$param" py -- python3 -c "
import time
import random
import sys

param = $param
print(f'Processing parameter {param}')
time.sleep(5)

# Simulate occasional failure
if random.random() < 0.2:  # 20% failure rate
    print(f'ERROR: Job {param} failed!')
    sys.exit(1)

print(f'SUCCESS: Job {param} completed')
")
        job_ids+=($job_id)
        echo "Submitted job $param: $job_id"
    done

    echo "Monitoring ${#job_ids[@]} jobs..."

    # Monitor and wait
    if qxub monitor --wait-for-completion --fail-fast "${job_ids[@]}"; then
        echo "✅ All jobs completed successfully"
        return 0
    else
        echo "❌ Some jobs failed"
        return 1
    fi
}

# Run with error handling
if submit_and_monitor; then
    echo "Pipeline completed successfully"
else
    echo "Pipeline failed - check job outputs"
    exit 1
fi
```

## Integration with Job Arrays

### PBS Job Array Alternative

Instead of PBS job arrays, use qxub patterns:

```bash
# Traditional PBS array alternative
submit_parameter_sweep() {
    local params=(0.1 0.5 1.0 2.0 5.0 10.0)
    local jobs=()

    for i in "${!params[@]}"; do
        local param="${params[$i]}"
        local job_id=$(qxub exec --terse --name "sweep-$i-param-$param" py -- python3 -c "
import numpy as np
import time

# Job array simulation
task_id = $i
param = ${param}

print(f'Task {task_id}: Processing parameter {param}')

# Simulate analysis
np.random.seed(task_id)
data = np.random.randn(1000)
result = np.mean(data) * param

time.sleep(10)
print(f'Task {task_id}: Result = {result:.4f}')

# Save result
with open(f'/tmp/result_task_{task_id}.txt', 'w') as f:
    f.write(f'{param},{result:.4f}\\n')
")
        jobs+=($job_id)
        echo "Task $i (param=$param): $job_id"
    done

    echo "Submitted ${#jobs[@]} jobs in parameter sweep"
    qxub monitor --summary "${jobs[@]}"
}

submit_parameter_sweep
```

## Workflow Coordination Patterns

### Sequential Stages with Dependencies

```bash
# Stage 1: Data preprocessing
echo "Stage 1: Data preprocessing"
PREPROCESS_JOBS=()
for dataset in A B C; do
    JOB_ID=$(qxub exec --terse --name "preprocess-$dataset" py -- python3 -c "
import time
dataset = '$dataset'
print(f'Preprocessing dataset {dataset}')
time.sleep(5)
print(f'Dataset {dataset} preprocessing complete')
")
    PREPROCESS_JOBS+=($JOB_ID)
done

echo "Waiting for preprocessing to complete..."
qxub monitor --wait-for-completion "${PREPROCESS_JOBS[@]}"

# Stage 2: Analysis (depends on preprocessing)
echo "Stage 2: Analysis"
ANALYSIS_JOBS=()
for dataset in A B C; do
    JOB_ID=$(qxub exec --terse --name "analyze-$dataset" sc -- python3 -c "
import time
dataset = '$dataset'
print(f'Analyzing preprocessed dataset {dataset}')
time.sleep(8)
print(f'Analysis of dataset {dataset} complete')
")
    ANALYSIS_JOBS+=($JOB_ID)
done

echo "Waiting for analysis to complete..."
qxub monitor --wait-for-completion "${ANALYSIS_JOBS[@]}"

# Stage 3: Summary report
echo "Stage 3: Generating summary report"
SUMMARY_JOB=$(qxub exec --terse --name "summary-report" py -- python3 -c "
import time
print('Generating summary report from all analyses')
time.sleep(3)
print('Summary report generated successfully')
")

qxub monitor --wait-for-completion "$SUMMARY_JOB"
echo "Pipeline completed successfully!"
```

## Key Takeaways

1. **`--terse` for scripting**: Essential for capturing job IDs in parallel workflows
2. **Multiple submission patterns**: for loops, while loops, find+xargs, GNU parallel, R scripts
3. **Centralized monitoring**: `qxub monitor` handles multiple jobs efficiently
4. **Blocking patterns**: Wait for completion with `--wait-for-completion`
5. **Error handling**: Use `--fail-fast` and exit codes for robust workflows

## Next Steps

Now that you understand parallel execution:
- **[DVC Integration](10-dvc-integration.md)** - Use these patterns in data science pipelines
- **[Remote Execution](11-remote-execution.md)** - Run parallel jobs from your laptop

Parallel execution is where qxub really shines for production workflows. These patterns form the foundation for sophisticated HPC pipelines.

---

**💡 Pro Tips:**
- Always use `--terse` when submitting multiple jobs programmatically
- Combine `qxub monitor --wait-for-completion` with job arrays for blocking workflows
- Use meaningful job names (`--name`) to track jobs in complex pipelines
- Monitor resource efficiency across parallel jobs to optimize future runs
- Consider GNU parallel for complex parameter combinations
