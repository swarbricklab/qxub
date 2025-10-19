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

### Memory Specification

qxub accepts various memory formats:

```bash
# Request 8GB of memory (also accepts MB: mem=8000MB)
qxub exec --dry --resources mem=8GB -- python analysis.py
```

### CPU Specification

```bash
# Request 4 CPUs for parallel processing
qxub exec --dry --resources ncpus=4 -- python parallel_analysis.py
```

### Walltime Specification

qxub accepts flexible walltime formats:

```bash
# Walltime formats: MM:SS or H:MM:SS
qxub exec --dry --default --resources walltime=1:30:00 -- python long_analysis.py
```

### Combining Resources

```bash
# A job that needs more resources
qxub exec --dry --default --resources mem=16GB,ncpus=4,walltime=4:00:00 -- echo "Resource-intensive job"
```

**Expected dry run output:**
```
📋 Job Configuration:
├── Resources: mem=16GB, ncpus=4, walltime=4:00:00
├── Queue: normal
...
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
# Typical pandas/analysis job
qxub exec --queue auto --resources mem=16GB,ncpus=2,walltime=2:00:00 -- python3 -c "
import pandas as pd
import numpy as np
print('Running data analysis...')
# Simulate some work
import time; time.sleep(5)
print('Analysis complete!')
"
```

### Machine Learning Training

```bash
# ML training with multiple cores
qxub exec --queue auto --resources mem=32GB,ncpus=8,walltime=6:00:00 -- python train_model.py
```

### Quick Test with Express Queue

```bash
# Fast turnaround for debugging
qxub exec --queue express --resources walltime=15:00 -- python quick_test.py
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
# Start small for testing
qxub exec --dry --resources mem=4GB,walltime=30:00 -- python3 my_script.py

# Scale up after confirming it works
qxub exec --dry --resources mem=16GB,walltime=2:00:00 -- python3 my_script.py
```

### 2. Use Auto Queue Selection

```bash
# Let qxub optimize for you
qxub exec --queue auto --resource mem=8GB,ncpus=4 -- my_analysis.py
```

### 3. Match Resources to Workload

- **I/O intensive**: More memory, fewer CPUs
- **CPU intensive**: More CPUs, standard memory
- **Memory intensive**: High memory, appropriate queue

### 4. Consider Walltime Carefully

```bash
# Better to overestimate slightly than underestimate
qxub exec --resource walltime=1:30:00 -- long_running_task.py  # If you think it takes 1 hour
```

## Monitoring Resource Usage

After jobs complete, qxub shows actual usage:

```bash
qxub exec --resources mem=8GB,walltime=1:00:00 -- python3 -c "
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
💡 Suggestion: This job could use --resources mem=1GB,walltime=30:00
```

## Key Takeaways

1. **Smart defaults**: The system provides sensible starting points
2. **Flexible formats**: Memory and walltime accept various formats
3. **Auto queue selection**: Use `--queue auto` for optimization
4. **Validation helps**: qxub warns about problematic resource requests
5. **Monitor efficiency**: Learn from actual usage patterns

## Next Steps

Now that you understand resource management:
- **[Execution Contexts](04-execution-contexts.md)** - Using different software environments
- **[Aliases](07-aliases.md)** - Save common resource combinations

Resource specification becomes intuitive with practice. The automatic queue selection feature eliminates much of the guesswork, while the efficiency reporting helps you optimize future jobs.

---

**💡 Pro Tips:**
- Use `--queue auto` unless you have specific queue requirements
- Always use `--dry` to verify resource requests before submitting
- Monitor efficiency reports to optimize future resource allocations
- Start with conservative estimates and scale up based on actual usage
