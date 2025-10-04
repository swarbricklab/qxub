# Config and Alias System Design

## Overview

This document outlines the design for qxub's configuration management and alias system, enabling users to define defaults for common options and create ultra-simple shortcuts for complex workflows.

## Design Goals

1. **Ultra-Simple Usage**: Enable patterns like `qxub alias dvc_push` for common workflows
2. **Flexible Configuration**: Support both simple defaults and complex workflow definitions
3. **Clear Precedence**: Unambiguous option resolution with helpful error messages
4. **Institutional Knowledge**: Allow encoding of best practices in shareable configurations

## Configuration System

### File Hierarchy

Configuration follows XDG Base Directory specification with clear precedence:

1. **System Config**: `${XDG_CONFIG_DIRS}/qxub/config.yaml` (e.g., `/etc/xdg/qxub/config.yaml`)
2. **User Config**: `~/.config/qxub/config.yaml`
3. **CLI Arguments**: Highest precedence

**Precedence Rule**: CLI Args > User Config > System Config > Built-in Defaults

### File Format

**Format**: YAML (human-readable, supports complex structures, excellent Python support)

**Structure**:
```yaml
# ~/.config/qxub/config.yaml
defaults:
  # Global qxub options
  name: "qt"
  queue: "normal" 
  project: "a56"
  joblog: "{name}.log"  # Template support
  resources:
    - "mem=4GB"
    - "ncpus=1"
  
  # Output paths (support templates)
  out: "/scratch/{project}/{user}/qt/{timestamp}/out"
  err: "/scratch/{project}/{user}/qt/{timestamp}/err"
  
  # Subcommand-specific defaults
  conda:
    env: "base"
    pre: "echo Starting conda job"
    post: null
  
  module:
    mod: ["python3"]
    pre: null
    post: "echo Module job completed"
  
  sing:
    sif: null  # Must be specified
    bind: ["/scratch", "/g/data"]
    env: []
    pre: null
    post: null

# Alias definitions
aliases:
  # Complete workflow aliases
  dvc_push:
    subcommand: "conda"          # Which subcommand to use
    cmd: "dvc push"              # The actual command to run
    name: "push" 
    queue: "copyq"
    conda:
      env: "dvc3"
  
  train_model:
    subcommand: "conda"
    cmd: "python train.py --epochs 100"
    name: "ml_training"
    resources: ["mem=32GB", "ncpus=8", "ngpus=1"]
    queue: "gpuvolta" 
    conda:
      env: "pytorch"
      pre: "nvidia-smi"
  
  # Resource-only aliases (no subcommand/cmd - used as modifiers)
  gpu:
    resources: ["mem=32GB", "ncpus=8", "ngpus=1"]
    queue: "gpuvolta"
  
  large:
    resources: ["mem=64GB", "ncpus=16"] 
    queue: "express"
  
  # Environment-only aliases
  ml_env:
    subcommand: "conda"
    conda:
      env: "pytorch"
      pre: "nvidia-smi"
```

### Template System

**Supported Variables**:
- `{user}` - Current username
- `{project}` - Project code  
- `{timestamp}` - Current timestamp (YYYYMMDD_HHMMSS)
- `{name}` - Job name
- `{queue}` - Queue name
- `{date}` - Current date (YYYYMMDD)
- `{time}` - Current time (HHMMSS)

**Example Templates**:
```yaml
defaults:
  joblog: "{name}_{date}_{time}.log"
  out: "/scratch/{project}/{user}/{name}/{timestamp}/out"
  err: "/scratch/{project}/{user}/{name}/{timestamp}/err"
```

## Alias System

### Core Concept

Aliases combine:
- **Subcommand**: Which qxub subcommand to use (conda/module/sing)
- **Command**: The actual command to execute
- **Configuration**: All qxub options (resources, environment, etc.)

### Usage Patterns

**Ultra-Simple Execution**:
```bash
qxub alias dvc_push                    # Execute alias as-is
qxub alias train_model                 # Execute alias as-is
```

**Command Appending**:
```bash
qxub alias quick_analysis input2.bam  # Append to existing command
# Becomes: samtools view -c input.bam input2.bam
```

**Option Overrides**:
```bash
qxub alias train_model --queue normal --resources mem=64GB
# Keeps cmd and conda env, changes queue and resources
```

**Complete Command Override**:
```bash
qxub alias train_model --cmd "python validate.py"
# Uses train_model's environment but different command
```

**Global Options**:
```bash
qxub --dry alias dvc_push              # Dry run the alias
qxub --quiet alias train_model         # Execute in quiet mode
```

### Alias Composition (Future Enhancement)

**Combining Aliases**:
```bash
qxub alias gpu+ml_env --cmd "python train.py"    # Combine resource + env
qxub alias train_model+large                     # Use base + resource modifier
```

### Alias Inheritance (Future Enhancement)

```yaml
aliases:
  base-ml:
    resources: ["mem=16GB", "ncpus=4"]
    conda:
      env: "pytorch"
  
  gpu-ml:
    inherits: "base-ml"  # Inherit from base-ml
    resources: ["mem=32GB", "ncpus=8", "ngpus=1"]
    queue: "gpuvolta"
```

## CLI Design

### Command Structure

**Main Commands**:
```bash
qxub [GLOBAL_OPTIONS] SUBCOMMAND [SUBCOMMAND_OPTIONS] [-- COMMAND_ARGS]
qxub [GLOBAL_OPTIONS] alias ALIAS_NAME [OVERRIDE_OPTIONS] [-- COMMAND_ARGS]
qxub config OPERATION [ARGS]
```

### Option Placement Rules

**Critical Rule**: Global options MUST come before subcommand/alias

**Valid Patterns**:
```bash
qxub --dry --name myjob alias dvc_push            # ✅ Global options first
qxub --quiet --resources mem=8GB conda --env dvc3 -- command  # ✅ Clear separation
qxub alias dvc_push --env dvc4 --queue normal     # ✅ Override alias properties
qxub alias dvc_push -- input.txt output.txt       # ✅ Command arguments after --
```

**Invalid Patterns** (with helpful errors):
```bash
qxub alias dvc_push --dry                         # ❌ Error: --dry is global
qxub conda --env dvc3 --quiet -- command          # ❌ Error: --quiet is global
```

**Error Messages**:
```bash
$ qxub alias dvc_push --dry
❌ Error: Global option --dry must come before 'alias'
💡 Did you mean: qxub --dry alias dvc_push
```

### Config Management CLI

```bash
# Config management
qxub config get defaults.name                    # Get single value
qxub config set defaults.name "myjob"           # Set single value
qxub config list                                 # Show all config
qxub config list defaults                       # Show section
qxub config edit                                 # Open in $EDITOR
qxub config reset                                # Reset to defaults
qxub config validate                             # Validate config file

# Config file management
qxub config show-files                          # Show config file locations
qxub config init                                # Create user config template

# Alias management
qxub config alias list                          # List all aliases
qxub config alias show dvc_push                # Show specific alias
qxub config alias test dvc_push                # Test and validate alias
qxub config alias set small --resources mem=4GB,ncpus=1 --queue normal
qxub config alias delete small
```

### Alias CLI

```bash
# Execute aliases
qxub alias dvc_push                             # Execute alias
qxub alias dvc_push input.txt                  # Execute with appended args
qxub --dry alias dvc_push                      # Dry run alias

# Test aliases
qxub alias test dvc_push                       # Detailed test with validation
qxub alias dvc_push --dry                      # Quick dry run (if allowed)

# List aliases  
qxub alias list                                # List all available aliases
```

## Implementation Plan

### Phase 1: Basic Config System
- YAML config loading with XDG paths
- Basic `qxub config` subcommand (get/set/list/edit)
- Simple default overrides for existing options
- Template variable resolution

### Phase 2: Alias System
- Alias definition and storage in config files
- `qxub alias` subcommand with basic execution
- Alias resolution and override logic
- Option placement validation with helpful errors

### Phase 3: Advanced Features
- Alias testing and validation (`qxub alias test`)
- Config file validation
- Alias inheritance
- Alias composition (combining multiple aliases)

## Technical Implementation Notes

### Config Loading Logic
```python
class ConfigManager:
    def __init__(self):
        self.system_config = self._load_system_config()
        self.user_config = self._load_user_config() 
        self.merged_config = self._merge_configs()
    
    def resolve_options(self, cli_args, aliases=None):
        """Resolve final options from config hierarchy + aliases + CLI"""
        # 1. Start with system defaults
        # 2. Apply user config
        # 3. Apply alias(es) in order
        # 4. Apply CLI arguments  
        # 5. Resolve templates (timestamp, user, project, etc.)
        pass
```

### Option Validation
```python
def validate_option_placement(ctx, param, value):
    """Ensure global options come before subcommands"""
    if param.name in GLOBAL_OPTIONS and ctx.info_name != 'qxub':
        raise click.BadParameter(
            f"Global option --{param.name} must come before subcommand.\n"
            f"Try: qxub --{param.name} {' '.join(ctx.parent.params)}"
        )
```

### Alias Resolution
```python
def resolve_alias(alias_name, overrides, global_opts):
    """Resolve alias with overrides and global options"""
    alias_def = config_manager.get_alias(alias_name)
    
    # Apply overrides to alias definition
    resolved = merge_configs(alias_def, overrides)
    
    # Apply global options
    resolved = merge_configs(resolved, global_opts)
    
    # Resolve templates
    resolved = resolve_templates(resolved)
    
    return resolved
```

## Benefits

1. **User Experience**: Ultra-simple `qxub alias dvc_push` for common workflows
2. **Flexibility**: Can still override any aspect when needed
3. **Discoverability**: `qxub alias list` and `qxub config list` show what's available
4. **Institutional Knowledge**: Aliases encode best practices and can be shared
5. **Predictability**: Clear rules about option placement and precedence
6. **Helpful**: Detailed error messages guide users to correct syntax
7. **Maintainable**: YAML config is human-readable and versionable

## Examples

### Simple Daily Workflows
```bash
# Data version control
qxub alias dvc_push
qxub alias dvc_pull  

# Model training
qxub alias train_model
qxub alias train_model --cmd "python train.py --epochs 200"

# Quick analysis
qxub alias analysis input.bam
```

### Override Examples
```bash
# Change resources
qxub alias train_model --resources mem=64GB,ngpus=2

# Change environment  
qxub alias dvc_push --env dvc4

# Dry run
qxub --dry alias train_model

# Quiet mode
qxub --quiet alias dvc_push
```

### Config Examples
```bash
# Set personal defaults
qxub config set defaults.project "xy12"
qxub config set defaults.conda.env "myenv"

# Create new alias
qxub config alias set my_workflow \
  --subcommand conda \
  --cmd "python analysis.py" \
  --env myenv \
  --resources mem=16GB,ncpus=4

# Test alias
qxub alias test my_workflow
```

This design provides a robust foundation for both simple config defaults and sophisticated workflow aliases while maintaining clear, predictable behavior.