# User Configuration Schema

## Overview

This document defines the enhanced user configuration schema for qxub v2.1+, which adds platform-aware features and intelligent resource selection capabilities.

## Enhanced Configuration Schema

The user configuration file (`~/.config/qxub/config.yaml`) is extended to support platform awareness and auto-selection preferences.

### Complete Configuration Example

```yaml
# ~/.config/qxub/config.yaml (enhanced for v2.1+)

# Preserve existing v2.0 defaults
defaults:
  project: a56
  queue: "auto"        # New: was "normal" - enables intelligent queue selection
  walltime: "00:30:00"
  mem: "4GB"
  ncpus: 1

  # New platform-aware defaults
  platform: "nci_gadi"  # Default platform for local execution

# Auto-selection preferences
auto_selection:
  queue:
    enabled: true
    strategy: "optimal"    # vs "conservative", "aggressive"

  resources:
    # How to handle resource conflicts
    cpus:
      policy: "auto_adjust"     # vs "warn", "error"
      prefer: "minimum_viable"  # vs "user_specified"

    memory:
      policy: "suggest"         # vs "auto_adjust", "error"

    walltime:
      policy: "warn"           # vs "auto_adjust", "error"

# Override behaviors per queue
queue_overrides:
  gpuvolta:
    auto_adjust_cpus: true    # Always adjust to minimum 12
    confirm_gpu_request: false # Don't prompt for GPU requests

  hugemem:
    confirm_high_memory: true # Prompt for >500GB requests
    suggest_alternatives: true

  copyq:
    max_auto_walltime: "01:00:00"  # Don't auto-select for longer jobs

# Platform profiles (v2.2 feature)
profiles:
  gadi:
    platform: nci_gadi
    defaults:
      project: a56
      queue: auto
    remote:
      host: gadi.nci.org.au
      user: jr9959
      qxub_command: "module load python3; qxub"

  laptop:
    platform: local
    defaults:
      queue: immediate
```

## Configuration Sections

### Defaults Section (Enhanced)

```yaml
defaults:
  # Existing v2.0 options (preserved)
  project: string
  queue: string          # Now supports "auto" for intelligent selection
  walltime: duration
  mem: memory_size
  ncpus: integer

  # New v2.1 options
  platform: string       # Default platform identifier
  auto_adjust: boolean   # Global auto-adjustment toggle
```

### Auto-Selection Section (New)

Controls intelligent queue and resource selection behavior.

```yaml
auto_selection:
  queue:
    enabled: boolean       # Enable/disable queue auto-selection
    strategy: string       # Selection strategy
    fallback: string       # Queue to use if auto-selection fails

  resources:
    <resource_type>:
      policy: string       # How to handle conflicts
      prefer: string       # Preference for conflict resolution
      max_adjustment: string # Maximum allowed auto-adjustment
```

#### Auto-Selection Strategies

- **`optimal`**: Select the most appropriate queue based on resources
- **`conservative`**: Prefer queues with higher limits (safer but potentially slower)
- **`aggressive`**: Prefer specialized queues (faster but potentially more restrictive)

#### Resource Policies

- **`auto_adjust`**: Automatically adjust conflicting resources
- **`suggest`**: Show suggestions but don't auto-adjust
- **`warn`**: Show warnings but proceed with user values
- **`error`**: Fail with error on conflicts

### Queue Overrides Section (New)

Per-queue behavior customization.

```yaml
queue_overrides:
  <queue_name>:
    auto_adjust_cpus: boolean     # Force CPU auto-adjustment for this queue
    auto_adjust_memory: boolean   # Force memory auto-adjustment
    confirm_high_memory: boolean  # Prompt for large memory requests
    confirm_gpu_request: boolean  # Prompt for GPU requests
    suggest_alternatives: boolean # Suggest alternative queues
    max_auto_walltime: duration   # Don't auto-select for longer jobs
```

### Platform Profiles Section (v2.2)

Defines execution profiles for different platforms.

```yaml
profiles:
  <profile_name>:
    platform: string       # Platform identifier
    defaults:              # Override defaults for this profile
      project: string
      queue: string
      # ... other defaults
    remote:                # Remote execution configuration
      host: string         # Remote hostname
      user: string         # Remote username
      port: integer        # SSH port (default: 22)
      qxub_command: string # Command to run qxub on remote
      working_dir: string  # Remote working directory
      sync_files: boolean  # Whether to sync files
```

## Resource Request Schema (Internal)

This represents how qxub internally processes user resource requests.

```yaml
# Internal representation after parsing user input
resource_request:
  # User-specified values
  user_specified:
    queue: string
    cpus: integer
    memory: memory_size
    gpu: integer
    walltime: duration

  # Platform-agnostic core resources
  core_resources:
    compute:
      cpus: integer
      memory: memory_size
      gpu: integer
      walltime: duration

  # Platform-specific resources
  platform_specific:
    queue: string
    project: string
    # ... other PBS-specific options

  # Resolution results
  resolved:
    platform: string
    queue: string           # After auto-selection
    adjusted_cpus: integer  # After auto-adjustment
    adjusted_memory: memory_size
    warnings: [string]      # Warnings to show user
    suggestions: [string]   # Suggestions to show user
    auto_adjustments: [string] # List of auto-adjustments made
```

## Configuration File Evolution

### v2.0 → v2.1 Migration

Existing v2.0 configuration files work unchanged:
- `queue: "normal"` → continues to work exactly as before
- New `queue: "auto"` → enables intelligent selection
- All existing defaults preserved

### v2.1 → v2.2 Migration

Adds platform profiles without breaking existing configuration:
- Local execution behavior unchanged
- New `--profile` option for remote execution
- Profile-specific defaults override global defaults

## Configuration Commands

Extended configuration management commands:

```bash
# View current configuration
qxub config show

# Set auto-selection preferences
qxub config set auto_selection.queue.enabled true
qxub config set auto_selection.queue.strategy optimal

# Set queue overrides
qxub config set queue_overrides.gpuvolta.auto_adjust_cpus true

# Manage profiles (v2.2)
qxub config profile create gadi --platform nci_gadi
qxub config profile set gadi.defaults.project a56
qxub config profile set gadi.remote.host gadi.nci.org.au
```

## Backward Compatibility

### Guaranteed Compatibility

1. **Existing commands work unchanged**
2. **Existing configuration files work unchanged**
3. **Explicit queue specification always honored**
4. **All v2.0 CLI options preserved**

### Opt-in Features

1. **Auto-selection**: Only active when `queue: "auto"` or config enabled
2. **Auto-adjustment**: Only active when policies configured
3. **Platform profiles**: Only active when `--profile` specified

## Validation Rules

1. **Platform references** must exist in platform definitions
2. **Queue overrides** must reference valid queues for the platform
3. **Resource policies** must be valid policy names
4. **Auto-selection strategies** must be recognized strategy names
5. **Profile names** must be unique and valid identifiers

## Future Extensions

### v2.3+ Features

```yaml
# Advanced scheduling preferences
scheduling:
  prefer_idle_queues: boolean
  avoid_peak_hours: boolean
  cost_optimization: boolean

# Team-wide configurations
team:
  shared_defaults: url
  resource_budgets:
    monthly_cpu_hours: integer
    monthly_gpu_hours: integer

# Workflow integration
workflows:
  dvc:
    auto_wrap_stages: boolean
    default_resources: resource_spec
  snakemake:
    profile_template: string
```
