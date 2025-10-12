# Remote Execution

Submit jobs to remote HPC systems via SSH.

## Quick Start

```bash
# Execute on remote system
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
    host: gadi
    conda_base: /apps/Modules/miniconda3
    scratch_base: /scratch/a56/{user}
    platform: nci_gadi
```

## Usage

```bash
# Basic remote execution
qxub --remote gadi --env pytorch -- python train.py

# With PBS options
qxub --remote gadi --queue gpu -l ngpus=1 --env pytorch -- python train.py

# Platform auto-selection works too
qxub --remote gadi --queue auto -l mem=500GB --env myenv -- python big_job.py
```

## Commands

```bash
qxub --remote NAME [options] -- command    # Execute remotely
qxub config get remotes                    # Show remote configs
```

## Requirements

- SSH key authentication
- qxub installed on remote system
- Platform definitions available remotely
