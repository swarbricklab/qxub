# Remote Execution

Submit jobs to remote HPC systems via SSH with automatic platform detection.

## Quick Start

```bash
# Simple remote execution - platform auto-detected
qxub --remote gadi --env myenv -- python script.py
```

## Setup

1. **Configure SSH** (in `~/.ssh/config`):
```
Host gadi
    HostName gadi.nci.org.au
    User your_username
    ControlMaster auto
    ControlPersist 10m
```

2. **Configure Remote** (in `~/.config/qxub/config.yaml`):
```yaml
remotes:
  gadi:
    url: ssh://gadi.nci.org.au
    # Optional overrides:
    platform: nci_gadi_custom  # Only if auto-detection insufficient
    conda_env: pytorch         # Only if different from local default
    working_dir: /scratch/a56/{user}/projects  # Custom working directory
```

## Usage

```bash
# Basic remote execution (platform auto-detected)
qxub --remote gadi --env pytorch -- python train.py

# With PBS options (queue selection happens on remote)
qxub --remote gadi --queue auto -l ngpus=1 --env pytorch -- python train.py

# Platform auto-selection works remotely too
qxub --remote gadi --queue auto -l mem=500GB --env myenv -- python big_job.py

# Override platform if needed
qxub --remote gadi --platform custom_gadi --env myenv -- python script.py
```

## How It Works

1. **Local qxub** connects via SSH to remote system
2. **Remote qxub** automatically detects platform based on hostname
3. **Remote qxub** uses its own platform definitions and queue knowledge
4. **Queue selection** happens on remote with current, accurate information
5. **TTY handling** automatically adapts to SSH context:
   - **Local execution**: Interactive spinners and overwriting progress messages
   - **Remote SSH**: Spinner disabled, clean line-by-line progress output

## Benefits

- **Simple configuration**: Only connection details needed
- **Always current**: Remote has latest platform/queue information
- **Self-configuring**: Each HPC system knows its own capabilities
- **Robust**: No platform file synchronization required
- **Adaptive output**: Clean display whether local or remote execution

## Commands

```bash
qxub --remote NAME [options] -- command    # Execute remotely
qxub config get remotes                    # Show remote configs
```

## Requirements

- SSH key authentication
- qxub installed on remote system
- Platform definitions available remotely
