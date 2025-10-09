# Platform Schema Design

## Overview

This document defines the schema for platform definitions in qxub v2.1+. Platforms encapsulate the queue constraints, resource limits, and auto-selection rules for different computational environments.

## Platform Definition Schema

Platforms are defined in YAML files under `~/.config/qxub/platforms/`

### Example: NCI Gadi Platform

```yaml
# ~/.config/qxub/platforms/nci_gadi.yaml
platform:
  name: nci_gadi
  type: pbs_pro
  host: gadi.nci.org.au
  
  # Queue definitions with constraints
  queues:
    normal:
      description: "General purpose compute"
      limits:
        max_cpus: 20736
        min_cpus: 1
        max_memory: "190GB"        # Per-job memory limit
        max_local_storage: "400GB"
        max_jobs_per_user: 300
        max_jobs_queued: 1000
      walltime_rules:
        # Complex NCI walltime rules based on core count
        - cores: "1-672"
          max_walltime: "48:00:00"
        - cores: "720-1440" 
          max_walltime: "24:00:00"
        - cores: "1488-2976"
          max_walltime: "10:00:00"
        - cores: "3024-20736"
          max_walltime: "5:00:00"
      default_walltime: "1:00:00"  # NCI default
      su_rate: 2.0               # SU per core-hour
      constraints:
        - no_gpu: true
      
    hugemem:
      description: "High memory compute"
      limits:
        max_cpus: 192
        min_cpus: 1
        max_memory: "1470GB"       # Per-job memory limit
        max_local_storage: "1400GB"
        min_memory: "193GB"        # Auto-trigger threshold
        max_jobs_per_user: 50
        max_jobs_queued: 1000
      walltime_rules:
        - cores: "1-48"
          max_walltime: "48:00:00"
        - cores: "96"
          max_walltime: "24:00:00"  
        - cores: "144,192"
          max_walltime: "5:00:00"
      default_walltime: "2:00:00"
      su_rate: 3.0
      constraints:
        - no_gpu: true
        
    gpuvolta:
      description: "V100 GPU compute"
      limits:
        max_cpus: 960
        min_cpus: 12              # NCI requirement for GPU jobs
        max_memory: "382GB"       # Per-job memory limit
        max_local_storage: "400GB"
        max_gpus: 4
        min_gpus: 1
        max_jobs_per_user: 50
        max_jobs_queued: 1000
      walltime_rules:
        - cores: "1-96"
          max_walltime: "48:00:00"
        - cores: "144-192"
          max_walltime: "24:00:00"
        - cores: "240-960"
          max_walltime: "5:00:00"
      default_walltime: "1:00:00"
      su_rate: 3.0
      constraints:
        - requires_gpu: true
        - gpu_types: ["V100"]
      resources:
        auto_min_cpus: 12
        
    express:
      description: "High priority compute (6x SU cost)"
      limits:
        max_cpus: 3168
        min_cpus: 1
        max_memory: "190GB"       # Same as normal queue
        max_local_storage: "400GB"
        max_jobs_per_user: 50
        max_jobs_queued: 1000
      walltime_rules:
        - cores: "1-480"
          max_walltime: "24:00:00"
        - cores: "528-3168"
          max_walltime: "5:00:00"
      default_walltime: "1:00:00"
      su_rate: 6.0               # Higher priority = higher cost
      constraints:
        - no_gpu: true
        
    copyq:
      description: "Data transfer queue"
      limits:
        max_cpus: 1
        min_cpus: 1
        max_memory: "190GB"       # Single core with normal memory
        max_local_storage: "200GB"
        max_jobs_per_user: 50
        max_jobs_queued: 1000
      walltime_rules:
        - cores: "1"
          max_walltime: "10:00:00"
      default_walltime: "1:00:00"
      su_rate: 2.0
      internet_connectivity: true  # Only copyq has internet access
      constraints:
        - data_transfer_only: true
        - single_core_only: true

  # Queue selection rules
  auto_selection:
    priority_order: ["normal", "hugemem", "gpuvolta", "copyq"]
    
    rules:
      - if: "gpu_requested > 0"
        then: "gpuvolta"
        
      - if: "memory > 192GB"
        then: "hugemem"
        
      - if: "walltime <= 10:00:00 AND cpus == 1"
        then: "copyq"
        
      - default: "normal"

  # Resource auto-adjustment policies  
  auto_adjust:
    min_cpus: auto      # Automatically adjust to queue minimum
    memory: suggest     # Suggest better queue, don't auto-adjust  
    walltime: user      # Always use user-specified value
```

## Schema Components

### Platform Metadata

```yaml
platform:
  name: string          # Unique platform identifier
  type: string          # Scheduler type (pbs_pro, slurm, k8s, etc.)
  host: string          # Hostname for remote execution
  description: string   # Optional human-readable description
```

### Queue Definition

```yaml
queues:
  <queue_name>:
    description: string
    limits:
      max_cpus: integer           # Maximum CPUs per job
      min_cpus: integer           # Minimum CPUs per job
      max_memory: memory_size     # Maximum memory per job
      max_local_storage: memory_size   # Local scratch storage limit
      max_gpus: integer           # Maximum GPUs per job
      min_gpus: integer           # Minimum GPUs per job
      min_memory: memory_size     # Minimum memory to auto-select this queue
      max_jobs_per_user: integer  # Concurrent jobs limit
      max_jobs_queued: integer    # Total queued jobs limit
    walltime_rules:              # Complex walltime rules based on core count
      - cores: "range"            # e.g., "1-48", "96", "144,192"
        max_walltime: duration    # Maximum walltime for this core range
    default_walltime: duration   # Default walltime if none specified
    su_rate: float              # Service Units per core-hour
    internet_connectivity: boolean  # Whether queue has internet access
    constraints:
      - constraint_expression   # Boolean expressions for validation
    resources:
      auto_min_cpus: integer   # Auto-adjust CPU minimum for this queue
      auto_max_memory: memory_size
```

### Auto-Selection Rules

```yaml
auto_selection:
  priority_order: [queue_names]   # Preference order for tie-breaking
  rules:
    - if: "boolean_expression"     # Resource requirement condition
      then: "queue_name"           # Queue to select
    - default: "queue_name"        # Fallback queue
```

### Auto-Adjustment Policies

```yaml
auto_adjust:
  <resource_type>: policy_value
```

**Policy Values:**
- **`auto`**: Automatically adjust resource to meet queue requirements
- **`suggest`**: Show suggestions but don't change user's request  
- **`user`**: Always use user-specified value, no automation
- **`disabled`**: Don't validate or suggest for this resource
- **`error`**: Fail if resource conflicts with queue constraints

**Example Configurations:**

```yaml
# Conservative platform (minimal automation)
auto_adjust:
  min_cpus: suggest
  memory: suggest  
  walltime: user

# Aggressive platform (maximum automation)  
auto_adjust:
  min_cpus: auto
  memory: auto
  walltime: auto

# NCI-specific (balanced approach)
auto_adjust:
  min_cpus: auto      # Always adjust GPU jobs to 12+ CPUs
  memory: suggest     # Suggest hugemem for >192GB
  walltime: user      # Users must think about time limits
```

## Data Types

### Duration Format
- `"HH:MM:SS"` - PBS walltime format
- `"1h"`, `"30m"`, `"45s"` - Human-readable format (converted internally)

### Memory Size Format
- `"4GB"`, `"512MB"`, `"2TB"` - Standard units
- `"4096"` - Assumed to be MB if no unit specified

### Core Range Format
- `"1-48"` - Continuous range (1 to 48 cores)
- `"144,192"` - Discrete values (exactly 144 or 192 cores)
- `"96"` - Single value (exactly 96 cores)
- `"1-672,720-1440"` - Multiple ranges combined

### Constraint Expressions
- Simple boolean logic: `"gpu_requested > 0"`
- Combined conditions: `"memory > 192GB AND cpus <= 48"`
- Available variables: `memory`, `cpus`, `gpu_requested`, `walltime`, `queue`, `su_cost`

## Validation Rules

1. **Queue names** must be unique within a platform
2. **Limits** must be internally consistent (min <= max)
3. **Auto-selection rules** must have exactly one default
4. **Constraint expressions** must be valid boolean expressions
5. **Resource types** in auto-adjust must correspond to known resource types
6. **Auto-adjustment policies** must be valid policy values
7. **Walltime rules** must cover all possible core counts without gaps
8. **CPU multiples** must be positive integers
9. **SU rates** must be positive numbers
10. **Core ranges** in walltime rules must be valid and non-overlapping

## Behavioral Examples

### SU Cost Estimation

```bash
# qxub can estimate Service Unit costs before submission
$ qxub --env pytorch -l ngpus=1 -l ncpus=12 -l walltime=4:00:00 --estimate-cost

→ Selected queue: gpuvolta
→ Estimated cost: 144 SU (12 cores × 4 hours × 3.0 SU/core-hour)
→ Continue? [y/N]
```

### Walltime Validation and Defaults

```bash
# Request exceeds queue limits
$ qxub --env analysis -l ncpus=1000 -l walltime=48:00:00

→ Selected queue: normal
→ Warning: 48h walltime not available for 1000 cores (max: 10h)
→ Suggested: -l walltime=10:00:00
→ Using default walltime: 1:00:00

# Default walltime applied
$ qxub --env myenv -l ncpus=48

→ Selected queue: normal  
→ Using default walltime: 1:00:00 (normal queue default)
→ Estimated cost: 96 SU (48 cores × 1 hour × 2.0 SU/core-hour)
```

### Auto-Adjustment in Action

```bash
# User requests GPU with insufficient CPUs
$ qxub --env pytorch -l ngpus=1 -l ncpus=4 -- python train.py

# With min_cpus: auto
→ Auto-adjusts to 12 CPUs (gpuvolta minimum)
→ Shows: "Adjusted CPUs from 4 to 12 (gpuvolta queue minimum)"

# With min_cpus: suggest  
→ Keeps 4 CPUs but shows suggestion
→ Shows: "Warning: gpuvolta queue requires minimum 12 CPUs. Consider: -l ncpus=12"

# With min_cpus: error
→ Fails with error
→ Shows: "Error: gpuvolta queue requires minimum 12 CPUs, got 4"
```

### Memory and Queue Selection

```bash
# User requests high memory
$ qxub --env myenv -l mem=300GB -- python big_analysis.py

# With memory: suggest
→ Auto-selects hugemem queue
→ Shows: "Selected hugemem queue for 300GB memory request"
→ Shows: "Estimated cost: 600 SU (1 core × 2 hours × 3.0 SU/core-hour)"

# With memory: auto (hypothetical)
→ Could auto-adjust memory up/down to fit queue constraints
→ Shows: "Adjusted memory from 300GB to 384GB (hugemem queue maximum)"
```

## Extension Points

### Custom Constraint Types
Platforms can define custom constraint validators:

```yaml
constraints:
  - type: "custom_validator"
    validator: "nci_gadi.validate_gpu_memory_ratio"
    params:
      max_ratio: 4
```

### Platform-Specific Resources
Platforms can define additional resource types:

```yaml
platform_resources:
  gpu_type: ["V100", "A100"]
  storage_type: ["scratch", "persistent"]
  network_type: ["infiniband", "ethernet"]
```

## Key Enhancements from NCI Analysis

### Complex Walltime Rules
NCI has sophisticated walltime limits that vary by core count:
- Different maximum walltimes for different core ranges
- Enables better resource utilization and scheduling efficiency
- qxub can validate against these rules and suggest appropriate limits

### Service Unit (SU) Cost Tracking
- Each queue has different SU rates (1.25x to 6x base cost)
- Enables cost estimation before job submission
- Can be combined with actual usage logs for cost analysis
- Helps users make informed decisions about queue selection

### Enhanced Resource Limits
- Local storage limits (important for I/O intensive jobs)
- Job queue limits (both running and queued)
- More granular constraint validation

### Additional Considerations for Implementation

1. **Default Walltime Reporting**: Show users what default will be applied
2. **Cost Estimation**: Calculate and display SU costs before submission
3. **Walltime Validation**: Check against complex core-count-dependent rules
4. **Resource Optimization**: Suggest more cost-effective configurations
5. **Usage Tracking**: Log actual vs estimated costs for analysis
6. **Queue Recommendation**: Suggest cheaper alternatives when appropriate