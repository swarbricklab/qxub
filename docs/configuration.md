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
qxub config set key value           # Set a value
qxub config get key                 # Get a value
qxub config list                    # Show all config
qxub config edit                    # Edit in $EDITOR
```

## Defaults

Set defaults for any PBS option:

```bash
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"
qxub config set defaults.resources '["mem=4GB", "ncpus=1"]'
```

## Aliases

Create shortcuts for common workflows:

```bash
# Create alias
qxub config alias set gpu --env pytorch --queue gpu -l ngpus=1 -l ncpus=12

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

- User: `~/.config/qxub/config.yaml`
- System: `/etc/xdg/qxub/config.yaml`

Priority: CLI args > User config > System config > Defaults
