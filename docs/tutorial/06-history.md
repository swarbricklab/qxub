# Job History: Tracking and Analyzing Your Work

qxub automatically tracks all your job submissions, creating a comprehensive history that helps you monitor progress, debug issues, and optimize resource usage. This section covers how to view, search, and analyze your job history.

## Understanding qxub History

qxub maintains two types of historical information:
1. **Execution records** - What you ran, when, and with what resources
2. **Resource efficiency** - How well jobs used allocated resources

All history is stored locally and tied to your user account.

## Basic History Commands

### View Recent Jobs

```bash
# Show the 10 most recent jobs
qxub history list
```

**Expected output:**
```
📊 Job History (10 most recent):

┌─────────────────┬──────────────┬─────────────────────┬─────────┬──────────┬─────────────┐
│ Job ID          │ Name         │ Submitted           │ Status  │ Duration │ Efficiency  │
├─────────────────┼──────────────┼─────────────────────┼─────────┼──────────┼─────────────┤
│ 12345691.gadi   │ qx-20241017… │ 2024-10-17 15:05:32 │ ✅ Done │ 00:00:15 │ 25% / 4%    │
│ 12345690.gadi   │ qx-20241017… │ 2024-10-17 15:00:52 │ ✅ Done │ 00:00:25 │ 42% / 12%   │
│ 12345689.gadi   │ qx-20241017… │ 2024-10-17 14:55:22 │ ✅ Done │ 00:00:08 │ 13% / 8%    │
│ 12345688.gadi   │ qx-20241017… │ 2024-10-17 14:50:15 │ ❌ Fail │ 00:00:03 │ 5% / 2%     │
│ 12345687.gadi   │ param-1.0    │ 2024-10-17 14:45:10 │ ✅ Done │ 00:00:22 │ 37% / 15%   │
└─────────────────┴──────────────┴─────────────────────┴─────────┴──────────┴─────────────┘

💡 Efficiency shows: Walltime% / Memory%
💡 Use 'qxub history show <job_id>' for detailed information
```

### Show More Jobs

```bash
# Show last 20 jobs
qxub history list --limit 20

# Show all jobs from today
qxub history list --since today

# Show jobs from the last week
qxub history list --since "1 week ago"
```

## Detailed Job Information

### View Single Job Details

```bash
# Get detailed info about a specific job
qxub history show 12345691.gadi-pbs
```

**Expected output:**
```
📋 Job Details: 12345691.gadi-pbs

Basic Information:
├── Name: qx-20241017-150532
├── Status: Completed ✅
├── Submitted: 2024-10-17 15:05:32
├── Started: 2024-10-17 15:05:45
├── Finished: 2024-10-17 15:06:00
└── Duration: 00:00:15

Resources:
├── Queue: normal
├── Project: a56
├── Requested: mem=4GB, ncpus=1, walltime=2:00:00
├── Used: mem=0.2GB, ncpus=1, walltime=00:00:15
└── Efficiency: 25% walltime, 4% memory

Command:
└── python3 -c "import pandas as pd; print('Hello from history')"

Execution Context:
├── Environment: dvc3 (conda)
├── Working Directory: /g/data/a56/software/qsub_tools
└── Platform: nci_gadi

Files:
├── Output: /scratch/a56/jr9959/qxub/qx-20241017-150532_20241017-150532.out
├── Error: /scratch/a56/jr9959/qxub/qx-20241017-150532_20241017-150532.err
└── Log: /scratch/a56/jr9959/qxub/qx-20241017-150532_20241017-150532.log

Exit Code: 0 (success)

💡 Use 'qxub history efficiency' to see resource usage patterns
```

### View Job Output Files

```bash
# Quickly view output from a historical job
qxub history show 12345691.gadi-pbs --output

# View stderr from a job
qxub history show 12345691.gadi-pbs --error

# View the PBS log
qxub history show 12345691.gadi-pbs --log
```

## Searching and Filtering History

### Filter by Status

```bash
# Show only failed jobs
qxub history list --status failed

# Show only successful jobs
qxub history list --status completed

# Show running jobs (if any)
qxub history list --status running
```

### Filter by Time

```bash
# Jobs from specific date
qxub history list --since "2024-10-17"

# Jobs from last 2 hours
qxub history list --since "2 hours ago"

# Jobs between dates
qxub history list --since "2024-10-16" --until "2024-10-17"
```

### Filter by Command or Environment

```bash
# Search for jobs containing specific text
qxub history search "pandas"

# Jobs that used specific environment
qxub history list --env dvc3

# Jobs that used specific queue
qxub history list --queue express
```

### Filter by Resource Usage

```bash
# Jobs that used more than 8GB memory
qxub history list --min-memory 8GB

# Jobs that ran for more than 1 hour
qxub history list --min-walltime 1:00:00

# Long-running jobs with low efficiency
qxub history list --min-walltime 30:00 --max-efficiency 20%
```

## Resource Efficiency Analysis

### Overall Efficiency Report

```bash
qxub history efficiency
```

**Expected output:**
```
📊 Resource Efficiency Analysis (last 50 jobs):

Walltime Efficiency:
├── Average: 28.5%
├── Median: 22.0%
├── Range: 5% - 85%
└── Under 50%: 42 jobs (84%)

Memory Efficiency:
├── Average: 12.3%
├── Median: 8.5%
├── Range: 2% - 45%
└── Under 25%: 48 jobs (96%)

Queue Usage:
├── normal: 35 jobs (70%)
├── express: 8 jobs (16%)
├── hugemem: 5 jobs (10%)
└── normalsl: 2 jobs (4%)

Resource Recommendations:
├── 🔍 Consider reducing walltime estimates (84% under-utilized)
├── 💾 Consider reducing memory requests (96% under-utilized)
├── ⚡ More jobs could use express queue for faster turnaround
└── 📊 Use 'qxub history optimize' for specific suggestions
```

### Resource Optimization Suggestions

```bash
qxub history optimize
```

**Expected output:**
```
🎯 Resource Optimization Suggestions:

Based on your last 50 jobs:

Memory Optimization:
├── Current typical request: 8GB
├── Actual typical usage: 1.2GB
├── Suggestion: Try --mem 2GB for most jobs
└── Potential savings: ~75% memory allocation

Walltime Optimization:
├── Current typical request: 2:00:00
├── Actual typical usage: 0:15:00
├── Suggestion: Try --walltime 30:00 for quick jobs
└── Potential savings: ~75% walltime allocation

Queue Optimization:
├── 15 jobs could have used express queue (completed <30min)
├── 3 jobs needed more resources than requested
└── Consider using --queue auto for better selection

Most Efficient Recent Jobs:
├── 12345685: 82% walltime, 35% memory (good balance)
├── 12345682: 78% walltime, 42% memory (excellent usage)
└── 12345679: 65% walltime, 38% memory (well optimized)

💡 Use these efficient jobs as templates for similar work
```

## Comparing Jobs

### Compare Resource Usage

```bash
# Compare two specific jobs
qxub history compare 12345690.gadi-pbs 12345691.gadi-pbs
```

**Expected output:**
```
📊 Job Comparison:

                    │ 12345690.gadi-pbs    │ 12345691.gadi-pbs
─────────────────────┼─────────────────────┼──────────────────────
Command             │ python analysis.py   │ python quick_test.py
Environment         │ dvc3                 │ dvc3
Duration            │ 00:00:25             │ 00:00:15
Walltime Requested  │ 2:00:00              │ 2:00:00
Walltime Efficiency │ 42% ⚠️               │ 25% ⚠️
Memory Requested    │ 8GB                  │ 4GB
Memory Used         │ 1.2GB                │ 0.2GB
Memory Efficiency   │ 12% ⚠️               │ 4% ⚠️
Queue              │ normal               │ normal
Status             │ ✅ Completed         │ ✅ Completed

💡 Both jobs are significantly under-utilizing resources
💡 Consider: --mem 2GB --walltime 30:00 for similar future jobs
```

## Exporting History Data

### Export to CSV

```bash
# Export recent history to CSV for analysis
qxub history export --format csv --output my_jobs.csv --limit 100
```

### Export Specific Fields

```bash
# Export only efficiency data
qxub history export --fields job_id,walltime_efficiency,memory_efficiency --format csv
```

## Practical History Usage Patterns

### Debug Failed Jobs

```bash
# Find recent failures
qxub history list --status failed --limit 5

# Get details about a failure
qxub history show 12345688.gadi-pbs

# View error output
qxub history show 12345688.gadi-pbs --error
```

### Track Parameter Sweeps

```bash
# Find all jobs from a parameter sweep
qxub history search "param-"

# Compare efficiency across parameters
qxub history list --pattern "param-*" --format efficiency
```

### Monitor Long-running Analysis

```bash
# Track jobs from a specific project
qxub history list --since "1 week ago" | grep "analysis"

# Find resource-intensive jobs
qxub history list --min-memory 16GB --min-walltime 1:00:00
```

### Optimize Resource Requests

```bash
# Find your most efficient jobs as templates
qxub history list --min-efficiency 50% --limit 10

# Analyze jobs with similar commands
qxub history search "pandas" --format efficiency
```

## History Maintenance

### Clean Old History

```bash
# Remove history older than 30 days
qxub history clean --older-than "30 days"

# Remove only failed job records
qxub history clean --status failed --older-than "7 days"
```

### Backup History

```bash
# Export complete history for backup
qxub history export --all --format json --output backup_$(date +%Y%m%d).json
```

## Key Takeaways

1. **Automatic tracking**: All jobs are automatically recorded
2. **Rich details**: Command, resources, efficiency, and timing all tracked
3. **Powerful filtering**: Search by time, status, resources, and content
4. **Efficiency insights**: Identify over-allocation and optimization opportunities
5. **Debugging support**: Quickly find and analyze failed jobs

## Next Steps

Now that you understand job history:
- **[Aliases](07-aliases.md)** - Save optimized resource combinations as shortcuts
- **[Configuration](08-configuration.md)** - Understand how settings affect all jobs

Job history is invaluable for optimizing your HPC workflows. Use it regularly to improve resource efficiency and debug issues.

---

**💡 Pro Tips:**
- Check `qxub history efficiency` weekly to optimize resource requests
- Use `qxub history show <job_id> --output` to quickly view results
- Export history data for deeper analysis in spreadsheets or Python
- Use failed job analysis to improve error handling in scripts
- Bookmark highly efficient jobs as templates for similar work
