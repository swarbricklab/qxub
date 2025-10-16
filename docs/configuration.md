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
