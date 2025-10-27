# Resources and Queues: Smart Allocation and Automatic Selection

Now that you understand qxub basics and debugging, let's explore how to specify resources and leverage qxub's intelligent queue selection. qxub makes resource specification simple while providing smart defaults and automatic optimization.

## Understanding Default Resources

Let's first understand what the system provides by default:

```bash
qxub config get defaults.resources
```

**Output:**
```
📋 Configuration: defaults.resources
└── ['mem=4GB', 'ncpus=1', 'walltime=2:00:00']
```

These defaults work well for:
- Simple scripts and data processing
- Single-threaded applications
- Short-running tasks

## Specifying Custom Resources

qxub provides two ways to specify resources: **workflow-friendly options** (recommended) and traditional PBS format.

### Workflow-Friendly Resource Options (Recommended)

These intuitive options make it easy to specify resources:

```bash
# Memory with --mem (accepts GB, MB, etc.)
qxub exec --dry --mem 8GB -- python analysis.py

# CPUs with --cpus
qxub exec --dry --cpus 4 -- python parallel_analysis.py

# Runtime/walltime with --time (accepts 2h, 30m, 1h30m, or HH:MM:SS)
qxub exec --dry --time 2h -- python long_analysis.py

# Disk/jobfs with --disk
qxub exec --dry --disk 50GB -- python temp_files.py

# Storage volumes with --volumes
qxub exec --dry --volumes "gdata/a56+scratch/a56" -- python data_processing.py
```

### Combining Workflow-Friendly Options

```bash
# Combine multiple resource options
qxub exec --dry --mem 16GB --cpus 4 --time 4h -- python analysis.py
```

**Expected dry run output:**
```
📋 Job Configuration:
├── Resources: mem=16GB, ncpus=4, walltime=4:00:00
├── Queue: normal
...
```

### Traditional PBS Format (Alternative)

You can also use traditional PBS resource format:

```bash
# PBS format with --resources or -l
qxub exec --dry --resources mem=16GB,ncpus=4,walltime=4:00:00 -- python analysis.py

# Or using -l (traditional PBS shorthand)
qxub exec --dry -l mem=16GB -l ncpus=4 -- python analysis.py
```

### Format Examples

#### Memory Formats
```bash
--mem 8GB        # Gigabytes
--mem 2000MB     # Megabytes
--mem 16g        # Case insensitive
```

#### Time Formats
```bash
--time 2h        # 2 hours
--time 30m       # 30 minutes
--time 1h30m     # 1 hour 30 minutes
--time 02:30:00  # 2 hours 30 minutes (HH:MM:SS)
```

#### Disk Formats
```bash
--disk 50GB      # Gigabytes
--disk 1000MB    # Megabytes
```

### Resource Priority

When options are combined, qxub follows this priority (highest to lowest):

1. Workflow-friendly CLI options (`--mem`, `--cpus`, `--time`)
2. Traditional PBS options (`--resources`, `-l`)
3. Shortcut/alias resources
4. Config file defaults

**Example:**
```bash
# Config has cpus: 4, your command has --cpus 8
qxub exec --cpus 8 -- python script.py
# Result: Uses 8 CPUs (CLI overrides config)
```

## Automatic Queue Selection

One of qxub's most powerful features is automatic queue selection with `--queue auto`:

### Basic Automatic Selection

```bash
# Let qxub choose the best queue
qxub exec --dry --default --queue auto --resources mem=8GB,ncpus=2,walltime=1:00:00 -- echo "Auto queue"
```

**Expected output:**
```
🎯 Automatic queue selection:
├── Requested: mem=8GB, ncpus=2, walltime=1:00:00
├── Evaluating queues: normal, express, normalsl, hugemem
├── Best match: normal (fits all constraints, lowest cost)
└── Selected: normal

📋 Job Configuration:
├── Queue: normal
├── Resources: mem=8GB, ncpus=2, walltime=1:00:00
...
```

### Auto Selection for High-Memory Jobs

```bash
# High memory automatically selects appropriate queue
qxub exec --dry --queue auto --resources mem=64GB -- echo "Big memory job"
```

**Expected output:**
```
🎯 Automatic queue selection:
├── Requested: mem=64GB, ncpus=1, walltime=2:00:00
├── Evaluating queues: normal, express, normalsl, hugemem
├── normal: ❌ Memory limit exceeded (max: ~45GB per CPU)
├── express: ❌ Memory limit exceeded
├── normalsl: ❌ Memory limit exceeded
├── hugemem: ✅ Fits all constraints
└── Selected: hugemem

📋 Job Configuration:
├── Queue: hugemem
├── Resources: mem=64GB, ncpus=1, walltime=2:00:00
...
```

### Auto Selection for Many CPUs

```bash
# Many CPUs might trigger normalsl selection
qxub exec --dry --queue auto --resources ncpus=32,walltime=2:00:00 -- echo "Many CPU job"
```

qxub analyzes:
- Available CPUs per queue
- Walltime limits for large jobs
- Cost efficiency (SU billing rates)

## Manual Queue Selection

You can also specify queues directly:

### Standard Queues

```bash
# Specify queue directly
qxub exec --dry --queue express --resources mem=8GB,walltime=30:00 -- python analysis.py

# High-memory jobs typically use hugemem queue
qxub exec --dry --queue hugemem --resources mem=128GB -- python big_analysis.py
```

### Understanding Queue Characteristics

Let's see what queues are available:

```bash
qxub platform show-queues
```

**Expected output:**
```
📋 Available Queues (nci_gadi):

normal:
├── Type: standard
├── Max CPUs: 20,736
├── Max Memory: 190GB
├── Max Walltime: 48:00:00
├── SU Rate: 2.0x
└── Walltime Rules: Shorter for larger jobs

express:
├── Type: standard
├── Max CPUs: 3,168
├── Max Memory: 190GB
├── Max Walltime: 24:00:00
├── SU Rate: 6.0x (premium)
└── Priority: High

normalsl:
├── Type: shared_memory
├── Max CPUs: 192 per node
├── Max Memory: 1.5TB per node
├── Max Walltime: 48:00:00
├── SU Rate: 2.0x
└── Use for: Large parallel jobs on single nodes

hugemem:
├── Type: high_memory
├── Max CPUs: 48
├── Max Memory: 1.5TB
├── Max Walltime: 48:00:00
├── SU Rate: 2.0x
└── Use for: Memory-intensive single-node jobs
```

## Smart Resource Examples

### Data Analysis Job

```bash
# Typical pandas/analysis job with workflow-friendly options
qxub exec --queue auto --mem 16GB --cpus 2 --time 2h -- python3 -c "
import pandas as pd
import numpy as np
print('Running data analysis...')
# Simulate some work
import time; time.sleep(5)
print('Analysis complete!')
"

# Or with traditional PBS format
qxub exec --queue auto --resources mem=16GB,ncpus=2,walltime=2:00:00 -- python analysis.py
```

### Machine Learning Training

```bash
# ML training with multiple cores (workflow-friendly)
qxub exec --queue auto --mem 32GB --cpus 8 --time 6h -- python train_model.py

# Same with PBS format
qxub exec --queue auto --resources mem=32GB,ncpus=8,walltime=6:00:00 -- python train_model.py
```

### Quick Test with Express Queue

```bash
# Fast turnaround for debugging
qxub exec --queue express --time 15m -- python quick_test.py
```

## Resource Validation and Warnings

qxub validates your resource requests:

### Over-allocation Warning

```bash
qxub exec --dry --resources mem=200GB --queue normal -- echo "Too much memory"
```

**Expected warning:**
```
⚠️  Warning: Requested memory (200GB) exceeds normal queue limit (~190GB)
💡 Suggestion: Use --queue auto or --queue hugemem
💡 Alternative: Reduce memory in --resources to 190GB or less
```

### Walltime vs CPU Rules

```bash
qxub exec --dry --resources ncpus=48,walltime 24:00:00 --queue normal -- echo "Large job"
```

**Expected warning:**
```
⚠️  Warning: Large jobs (48+ CPUs) have reduced walltime limits
💡 Max walltime for 48 CPUs: 5:00:00 in normal queue
💡 Suggestion: Reduce walltime in --resources to 5:00:00 or use fewer CPUs
```

## Best Practices

### 1. Start Conservative, Scale Up

```bash
# Start small for testing (workflow-friendly)
qxub exec --dry --mem 4GB --time 30m -- python3 my_script.py

# Scale up after confirming it works
qxub exec --dry --mem 16GB --time 2h -- python3 my_script.py
```

### 2. Use Workflow-Friendly Options

```bash
# Easier to read and write
qxub exec --mem 8GB --cpus 4 --time 2h -- my_analysis.py

# Instead of:
qxub exec --resources mem=8GB,ncpus=4,walltime=2:00:00 -- my_analysis.py
```

### 3. Use Auto Queue Selection

```bash
# Let qxub optimize for you
qxub exec --queue auto --mem 8GB --cpus 4 -- my_analysis.py
```

### 4. Match Resources to Workload

- **I/O intensive**: More memory, fewer CPUs
- **CPU intensive**: More CPUs, standard memory
- **Memory intensive**: High memory, appropriate queue

### 5. Consider Walltime Carefully

```bash
# Better to overestimate slightly than underestimate
qxub exec --time 1h30m -- long_running_task.py  # If you think it takes 1 hour
```

## Monitoring Resource Usage

After jobs complete, qxub shows actual usage:

```bash
qxub exec --mem 8GB --time 1h -- python3 -c "
import time
print('Working...')
time.sleep(10)
print('Done!')
"
```

**End of job output:**
```
🎉 Job completed successfully (exit code: 0)
📊 Walltime used: 00:00:15 / 01:00:00 (25% efficiency)
💾 Memory used: 0.3GB / 8.0GB (4% efficiency)
💡 Suggestion: This job could use --mem 1GB --time 30m
```

## Key Takeaways

1. **Use workflow-friendly options**: `--mem`, `--cpus`, `--time` are easier than PBS format
2. **Smart defaults**: The system provides sensible starting points
3. **Flexible formats**: Accept various formats (GB/MB, 2h/1h30m/HH:MM:SS)
4. **Auto queue selection**: Use `--queue auto` for optimization
5. **Priority matters**: CLI options override shortcuts, aliases, and config defaults
6. **Validation helps**: qxub warns about problematic resource requests
7. **Monitor efficiency**: Learn from actual usage patterns

## Next Steps

Now that you understand resource management:
- **[Execution Contexts](04-execution-contexts.md)** - Using different software environments
- **[Aliases](07-aliases.md)** - Save common resource combinations

Resource specification becomes intuitive with practice. The automatic queue selection feature eliminates much of the guesswork, while the efficiency reporting helps you optimize future jobs.

---

**💡 Pro Tips:**
- **Use workflow-friendly options** (`--mem 8GB --cpus 4`) instead of PBS format for readability
- Use `--queue auto` unless you have specific queue requirements
- Always use `--dry` to verify resource requests before submitting
- Monitor efficiency reports to optimize future resource allocations
- Start with conservative estimates and scale up based on actual usage
- CLI options always override shortcuts, aliases, and config defaults
