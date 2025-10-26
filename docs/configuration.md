# Configuration

Set defaults to avoid repeating common options.

## Quick Setup

```bash
# Set defaults once
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# View configuration
qxub config list
```

## Commands

```bash
qxub config set key value           # Set a value (user config, default)
qxub config set --global key value  # Set in user config (explicit)
qxub config set --system key value  # Set system-wide value
qxub config set --project key value # Set project-level value
qxub config set --local key value   # Set local project value
qxub config set --test key value    # Set test/CI value
qxub config get key                 # Get a value
qxub config list                    # Show all config (merged)
qxub config list --user-only        # Show only user config
qxub config list --system-only      # Show only system config
qxub config list --show-origin      # Show source file for each setting
qxub config edit                    # Edit in $EDITOR
qxub config init-project            # Initialize project config directory
```

## System vs User Configuration

qxub supports both system-wide and user-specific configuration:

### User Configuration (Default)
```bash
# Set user-specific defaults (both forms are equivalent)
qxub config set defaults.project "a56"
qxub config set --global defaults.queue "normal"

# View only user configuration
qxub config list --user-only
```

### System Configuration (Admin)
```bash
# Set system-wide defaults (requires write access to system config)
qxub config set --system defaults.project "a56"
qxub config set --system monitor.default_suffix ".gadi-pbs"

# View only system configuration
qxub config list --system-only
```

**Precedence**: CLI args > Test config > Local config > Project config > User config > System config > Built-in defaults

## Project-Level Configuration

Configure settings specific to a project or repository:

```bash
# Initialize project configuration
qxub config init-project

# Set team-shared defaults (git-tracked)
qxub config set --project qxub.defaults.walltime "2:00:00"
qxub config set --project qxub.defaults.queue "normal"

# Set personal overrides (git-ignored)
qxub config set --local qxub.defaults.queue "express"

# Set CI/testing defaults (git-tracked)
qxub config set --test qxub.defaults.walltime "0:10:00"
qxub config set --test qxub.defaults.queue "copyq"
```

### Project Config Types

- **Project** (`.qx/project.yaml`) - Git-tracked team settings
- **Local** (`.qx/local.yaml`) - Git-ignored personal overrides
- **Test** (`.qx/test.yaml`) - Git-tracked CI/testing settings

Use `qxub config list --show-origin` to see where each setting comes from.

## Defaults

Set defaults for any PBS option:

```bash
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"
qxub config set defaults.resources '["mem=4GB", "ncpus=1"]'
```

## Monitor Configuration

Configure monitoring behavior:

```bash
# Set default suffix removal (cleaner display)
qxub config set monitor.default_suffix ".gadi-pbs"

# Set default refresh interval (seconds)
qxub config set monitor.default_interval 15

# Now monitor commands use these defaults:
qxub monitor 12345.gadi-pbs 12346.gadi-pbs
# Displays as: 12345, 12346 (suffix removed)
# Updates every 15 seconds
```

Monitor options with defaults:
- `monitor.default_suffix` - Suffix to remove from job IDs in display
- `monitor.default_interval` - Refresh interval in seconds (default: 30)

## Aliases

Create shortcuts for common workflows:

```bash
# Create alias
qxub config alias set gpu --env pytorch --queue auto -l ngpus=1 -l ncpus=12

# Use alias
qxub alias gpu -- python train.py

# List aliases
qxub config alias list
```

## Shortcuts

Create command-specific execution profiles with automatic detection:

```bash
# Create user shortcuts (default)
qxub config shortcut set "python" --env base --description "Python with base conda environment"
qxub config shortcut set "dvc" --env dvc3 --resources mem=16GB,ncpus=8

# Create system-wide shortcuts (requires admin permissions)
qxub config shortcut set "dvc data status" --system --env dvc3 --resources mem=64GB,ncpus=16

# List shortcuts with source information
qxub config shortcut list --show-origin

# Use shortcuts (automatic detection)
qxub exec -- python script.py      # Auto-detects 'python' shortcut
qxub exec -- dvc data status       # Auto-detects 'dvc data status' shortcut

# Show shortcut details
qxub config shortcut show "python"

# Delete shortcuts
qxub config shortcut delete "python"                    # Delete user shortcut
qxub config shortcut delete "system-shortcut" --system  # Delete system shortcut

# Built-in command aliases
qx --env myenv -- python script.py     # Short for 'qxub exec'
qxet "ml-pipeline" --env pytorch        # Short for 'qxub config shortcut set'
```

### Shortcut Files
- **System**: `/etc/xdg/qxub/shortcuts.json` - System-wide shortcuts (admin)
- **User**: `~/.config/qxub/shortcuts.json` - Personal shortcuts

**Priority**: User shortcuts override system shortcuts for the same command prefix.

## Template Variables

Use variables in paths and names:

```bash
qxub config set defaults.name "{user}-{timestamp}"
qxub config set defaults.out "/scratch/{user}/logs/{timestamp}/out"
```

Variables: `{user}`, `{project}`, `{timestamp}`

## Configuration Files

- **System**: `/etc/xdg/qxub/config.yaml` - System-wide defaults
- **User**: `~/.config/qxub/config.yaml` - Personal defaults
- **Project**: `.qx/project.yaml` - Team-shared project settings
- **Local**: `.qx/local.yaml` - Personal project overrides (git-ignored)
- **Test**: `.qx/test.yaml` - CI/testing settings

**Priority**: CLI args > Test > Local > Project > User > System > Defaults

## Resource Resolution (v3.2.3+)

Resources are resolved with intelligent per-resource-type precedence to give you fine-grained control:

### Resolution Priority (Highest → Lowest)

1. **CLI Workflow Options** (`--mem`, `--cpus`, `--time`, `--disk`, `--volumes`)
   - Always takes precedence when specified
   - Example: `qxub exec --cpus 48 --mem 128GB -- python script.py`

2. **CLI `--resources` Option**
   - Direct PBS resource specifications
   - Example: `qxub exec --resources ncpus=48,mem=128GB -- python script.py`

3. **Shortcut/Alias Resources**
   - Applied only if not overridden by CLI
   - Shortcuts and aliases are mutually exclusive

4. **Config Workflow Defaults** (per-resource-type)
   - Applied intelligently: each resource type only if not already specified
   - Config precedence: User > System
   - Example: If you specify `--cpus 48` but not `--mem`, config mem defaults still apply

5. **Config `resources` List**
   - Applied only if NO resources specified anywhere else

### Smart Per-Resource Resolution

The key feature of v3.2.3+ is **per-resource-type intelligence**:

```bash
# User config has: cpus: 4, mem: 8GB, runtime: 2h

# Specify only cpus on CLI
qxub exec --cpus 48 -- python script.py
# Result: ncpus=48 (CLI), mem=8GB (config), walltime=2:00:00 (config)

# Mix --resources and workflow options
qxub exec --resources ncpus=48 --mem 32GB -- python script.py
# Result: ncpus=48 (--resources), mem=32GB (--mem override), walltime=2:00:00 (config)

# Shortcut with partial CLI override
# Shortcut 'bigdata': resources=[mem=64GB, ncpus=24]
qxub exec --shortcut bigdata --cpus 48 -- process.py
# Result: mem=64GB (shortcut), ncpus=48 (CLI override), walltime=2:00:00 (config)
```

### Config Merging

Configuration files merge hierarchically:

```yaml
# System config (/etc/xdg/qxub/config.yaml):
defaults:
  cpus: 4
  mem: 8GB
  runtime: 2h
  disk: 15GB

# User config (~/.config/qxub/config.yaml):
cpus: 8        # Overrides system
mem: 16GB      # Overrides system
# runtime and disk inherited from system

# Final merged defaults:
# cpus: 8, mem: 16GB, runtime: 2h, disk: 15GB
```

### Key Behaviors

**Prevents Duplication**: If you specify `--resources ncpus=48`, config defaults won't add `ncpus=4`. But they will still add `mem=8GB` if memory wasn't specified.

**Resource Key Mapping**: PBS keys (`ncpus`, `walltime`, `jobfs`, `storage`) are mapped to workflow keys (`cpus`, `runtime`, `disk`, `volumes`) to detect what's already specified.

**All or Nothing for Config**: The old-style `defaults.resources: ["mem=4GB", "ncpus=1"]` list is only applied if NO resources are specified via CLI, shortcuts, or aliases.
