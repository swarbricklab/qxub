# Remote Execution v3.3.0

Execute qxub jobs on remote HPC systems via SSH with platform definition delegation.

## Quick Start

```bash
# Platform delegation - nci_gadi definition resolved on remote
qxub exec --config my_config.yaml --platform nci_gadi --env myenv -- python script.py
```

## Architecture Overview

v3.3.0 remote execution uses **platform definition delegation**:

1. **Local Config**: Provides SSH connection details and working directory
2. **Remote Resolution**: Platform definitions are resolved on the target system
3. **Name Contract**: Platform names must exist on the remote system

## Setup

### 1. SSH Configuration

Configure SSH access in `~/.ssh/config`:
```
Host gadi
    HostName gadi.nci.org.au
    User your_username
    IdentityFile ~/.ssh/id_rsa
    ControlMaster auto
    ControlPersist 10m
```

### 2. Local qxub Configuration

Create a config file with remote platform definitions:

```yaml
# laptop_config.yaml - for laptop → Gadi execution
defaults:
  platform: nci_gadi
  project: a56
  walltime: "1:00:00"
  mem: 4GB
  cpus: 1

platforms:
  # Delegated remote platform - NO definition field
  nci_gadi:
    remote:
      host: gadi                          # SSH hostname from ~/.ssh/config
      working_dir: /scratch/a56/{{user}}  # {{user}} = remote username
      conda_init: |
        eval "$(conda shell.bash hook)"
        conda activate qxub

  # Standard remote platform - WITH definition field
  nci_gadi_explicit:
    definition: https://raw.githubusercontent.com/swarbricklab/qxub-platforms/main/nci_gadi.yaml
    remote:
      host: gadi
      working_dir: /scratch/a56/{{user}}
      conda_init: |
        eval "$(conda shell.bash hook)"
        conda activate qxub
```

### 3. Remote System Setup

Ensure the remote system has:
- qxub installed and available (e.g., in conda environment)
- Platform definition for the platform name you're using
- Working directory exists

## Template Variables

v3.3.0 supports context-aware template variables:

- **`{user}`**: Resolved on **config system** (laptop username)
- **`{{user}}`**: Resolved on **remote system** (remote username via `$USER`)
- **`{project}`**: Resolved on **config system** (laptop `$PROJECT`)
- **`{{project}}`**: Resolved on **remote system** (remote `$PROJECT`)

Example:
```yaml
# Mixed template resolution
working_dir: /home/{user}/scratch/{{user}}/jobs
# Becomes: /home/laptop_user/scratch/$USER/jobs
```

## Execution Modes

### Delegated Execution (Recommended)

Platform name forwarded to remote, definition resolved remotely:

```yaml
platforms:
  nci_gadi:  # Must exist on remote system
    remote:
      host: gadi
      working_dir: /scratch/a56/{{user}}
      # NO definition field = delegation
```

**Flow:**
1. Local: `qxub exec --platform nci_gadi --config laptop_config.yaml ...`
2. SSH: `qxub exec --platform nci_gadi ...` (no --config passed)
3. Remote: Resolves `nci_gadi` from its own config system

### Explicit Definition Execution

Platform definition provided locally, sent to remote:

```yaml
platforms:
  nci_gadi:
    definition: https://raw.githubusercontent.com/.../nci_gadi.yaml # explicitly point to platform definition file
    remote:
      host: gadi
      working_dir: /scratch/a56/{{user}}
```

## Usage Examples

### Basic Remote Execution
```bash
# Basic remote execution with 'nci_gadi' platform defined as remote platform in local config
qxub exec --platform nci_gadi --dry -- echo "hello"

# With conda environment
qxub exec --platform nci_gadi --env pytorch -- python train.py

# With resource specifications
qxub exec  --platform nci_gadi --mem 16GB --cpus 8 --env myenv -- python script.py
```

### CI/Automation
```bash
# GitHub Actions / CI environments
qxub exec --config ci_config.yaml --platform nci_gadi --walltime "0:05:00" --queue copyq -- echo "CI test"
```

### Testing
```bash
# Dry run to verify SSH and command construction
qxub exec --config my_config.yaml --platform nci_gadi --dry -vvv -- echo "test"

# This shows the SSH command that will be executed:
# ssh gadi 'cd /scratch/a56/$USER && eval "$(conda shell.bash hook)" && conda activate qxub && qxub exec --platform nci_gadi ...'
```

## Platform Name Requirements

**Critical**: Platform names must exist on the remote system for delegation to work.

✅ **Correct**: Use real platform names
```yaml
platforms:
  nci_gadi:    # Exists on Gadi
  pawsey:      # Exists on Pawsey
  m3:          # Exists on M3
```

❌ **Incorrect**: Made-up platform names
```yaml
platforms:
  my_custom_gadi:     # Won't exist on remote
  gadi_remote:        # Arbitrary name
  test_platform:      # Made-up name
```

## Output Handling

v3.3.0 provides clean output for remote execution:

- **No TTY allocation**: Prevents staggered output from remote progress indicators
- **Linear output**: Clean, sequential command output
- **CI-friendly**: Perfect for automation and logging
- **Verbose modes**: Use `-v`, `-vv`, `-vvv` for debugging

## Error Handling

### Common Issues

1. **Platform not found on remote**:
   ```
   Error: Platform 'my_platform' not found
   ```
   **Solution**: Use a platform name that exists on the remote system

2. **SSH connection fails**:
   ```
   ssh: Could not resolve hostname gadi: Name or service not known
   ```
   **Solution**: Check `~/.ssh/config` and network connectivity

3. **Working directory doesn't exist**:
   ```
   bash: cd: /scratch/project/user: No such file or directory
   ```
   **Solution**: Ensure working directory exists on remote system

4. **Remote qxub not found**:
   ```
   bash: qxub: command not found
   ```
   **Solution**: Check conda_init activates correct environment

## Requirements

- **Local**: qxub v3.3.0+, SSH client, config file
- **Remote**: qxub installed, platform definitions, working directory
- **Network**: SSH access to remote system
- **Authentication**: SSH key or other non-interactive auth

## Migration from v2.x

v3.3.0 replaces the `--remote` flag with config-based platform definitions:

```bash
# v2.x (deprecated)
qxub --remote gadi --env myenv -- python script.py

# v3.3.0 (new)
qxub exec --config my_config.yaml --platform nci_gadi --env myenv -- python script.py
```

Benefits of v3.3.0 approach:
- **Explicit configuration**: Clear separation of concerns
- **Version control friendly**: Config files can be committed
- **Template variables**: Flexible path resolution
- **Better testing**: Dry-run support for debugging
- **CI integration**: Native support for automation workflows
