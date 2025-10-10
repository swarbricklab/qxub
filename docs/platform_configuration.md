# Platform Configuration

## Overview

Starting with qxub v2.1, the system includes intelligent queue selection and resource adjustment to help you submit jobs more efficiently. This document explains how these features work and how they can be configured.

## What Are Platforms?

A "platform" in qxub represents a computational environment like NCI Gadi, with its specific queues, resource limits, and scheduling policies. qxub uses platform definitions to:

- **Automatically select the best queue** for your job requirements
- **Suggest resource adjustments** when your requests don't fit queue constraints
- **Prevent common submission errors** before they reach the scheduler

## How Queue Selection Works

When you submit a job, qxub analyzes your resource requirements and automatically chooses the most appropriate queue:

```bash
# GPU job automatically goes to GPU queue
$ qxub --env pytorch -l ngpus=1 -- python train.py
→ Selected: gpuvolta queue

# High memory job goes to high-memory queue
$ qxub --env analysis -l mem=300GB -- python big_analysis.py
→ Selected: hugemem queue

# Regular job uses default queue
$ qxub --env myenv -- python script.py
→ Selected: normal queue
```

## Resource Auto-Adjustment

qxub can help fix common resource specification issues:

### CPU Requirements
Some queues (like GPU queues) have minimum CPU requirements:

```bash
# Original request
$ qxub --env pytorch -l ngpus=1 -l ncpus=4 -- python train.py

# qxub automatically adjusts (if configured)
→ Adjusted CPUs from 4 to 12 (gpuvolta queue minimum)
→ Job submitted with: -l ncpus=12
```

### Memory Suggestions
For memory-intensive jobs, qxub can suggest better queues:

```bash
$ qxub --env analysis -l mem=300GB -- python script.py
→ Selected hugemem queue for 300GB memory request
```

## Configuration Levels

### System Level
Your system administrator configures the platform definitions that describe:
- Available queues and their limits
- Queue selection rules
- Resource adjustment policies

### User Level
You can customize how aggressive qxub is with automation in `~/.config/qxub/config.yaml`:

```yaml
# Conservative approach - just suggest improvements
platform:
  auto_adjust:
    min_cpus: suggest    # Show suggestions, don't auto-change
    memory: suggest      # Suggest better queues
    walltime: user       # Never change walltime

# Aggressive approach - fix issues automatically
platform:
  auto_adjust:
    min_cpus: auto       # Automatically fix CPU requirements
    memory: auto         # Auto-select best queue and adjust memory
    walltime: suggest    # Suggest walltime improvements
```

## Adjustment Policies Explained

Each resource type can have different adjustment behavior:

| Policy | Behavior | Example |
|--------|----------|---------|
| `auto` | Automatically fix conflicts | Changes 4 CPUs → 12 CPUs for GPU jobs |
| `suggest` | Show helpful suggestions | "Consider using hugemem queue for 300GB" |
| `user` | Always use your exact specification | Never changes your request |
| `error` | Fail if there's a conflict | Stops submission with clear error message |

## Common Scenarios

### GPU Jobs
```bash
# What you type
$ qxub --env pytorch -l ngpus=1 -l ncpus=4

# What qxub does (with auto-adjustment)
→ Detects GPU request
→ Selects gpuvolta queue
→ Notices 4 CPUs < 12 minimum
→ Auto-adjusts to 12 CPUs
→ Submits job successfully
```

### Memory-Intensive Jobs
```bash
# What you type
$ qxub --env analysis -l mem=300GB

# What qxub does
→ Detects high memory request
→ Selects hugemem queue (>192GB trigger)
→ Shows: "Selected hugemem queue for high memory job"
→ Submits to appropriate queue
```

### Simple Jobs
```bash
# What you type
$ qxub --env myenv

# What qxub does
→ Uses default queue (normal)
→ Applies default resources
→ No adjustments needed
```

## Overriding Auto-Selection

You can always override qxub's queue selection:

```bash
# Force a specific queue
$ qxub --env myenv -q copyq -- python quick_script.py

# Disable auto-adjustment for one job
$ qxub --env myenv --no-auto-adjust -l ncpus=4 -l ngpus=1
```

## Viewing Platform Information

See what queues and limits are available:

```bash
# List available queues
$ qxub queues

# Show queue details and limits
$ qxub queues --details

# Show platform configuration
$ qxub platform info
```

## Troubleshooting

### "Queue not found" errors
Make sure you're using the correct platform. Check available queues:
```bash
$ qxub queues
```

### Unexpected resource adjustments
Check your auto-adjustment settings:
```bash
$ qxub config show platform.auto_adjust
```

### Jobs going to wrong queue
Review the queue selection rules:
```bash
$ qxub platform info --selection-rules
```

## Advanced Configuration

### Per-Platform Settings
If you work with multiple platforms, you can configure each differently:

```yaml
platforms:
  nci_gadi:
    auto_adjust:
      min_cpus: auto     # NCI has strict GPU→CPU requirements
      memory: suggest    # Just suggest better queues

  university_hpc:
    auto_adjust:
      min_cpus: suggest  # More flexible environment
      memory: auto       # Automatically optimize memory
```

### Disabling Features
Turn off intelligent features entirely:

```yaml
platform:
  auto_selection: false  # Always use default queue
  auto_adjust:
    min_cpus: user      # Never adjust anything
    memory: user
    walltime: user
```

## Getting Help

- `qxub --help` - General help
- `qxub queues --help` - Queue information commands
- `qxub platform --help` - Platform configuration commands
- Check the [full documentation](README.md) for more details

The goal is to make job submission smoother while keeping you in control. If the automation ever gets in your way, you can always override or disable specific features.
