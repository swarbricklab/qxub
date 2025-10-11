# qxub Remote Configuration Guide

## Overview

qxub v2.2 introduces remote execution capabilities, allowing you to submit jobs to remote HPC systems from your local machine. This document explains how to configure remote systems for seamless execution.

## Configuration Architecture

qxub uses a clean separation between:
- **Platform definitions**: Describe system capabilities (queues, limits) - stored on the remote system
- **User configuration**: Define how to connect and execute - stored locally (`~/.config/qxub/config.yaml`)
- **SSH configuration**: Handle credentials and connection details - standard `~/.ssh/config`

## User Configuration File

Create `~/.config/qxub/config.yaml` to define your remote systems:

```yaml
# qxub user configuration
remotes:
  nci_gadi:
    # Connection URL (protocol://host:port)
    url: ssh://gadi.nci.org.au

    # Configuration file for connection (optional, smart defaults)
    # For SSH: defaults to ~/.ssh/config
    # For K8s: would default to ~/.kube/config (future)
    config: ~/.ssh/config

    # Conda environment containing qxub on remote system
    qxub_env: qxub-prod
    # OR absolute path to conda environment:
    # qxub_env: /g/data/a56/software/qsub_tools

    # Platform definition file on remote system
    platform_file: /g/data/a56/software/qsub_tools/docs/platforms/nci_gadi.yaml

    # Template for project directories on remote system
    # Supports ${PROJECT} and ${USER} environment variable substitution
    project_root_dir: /scratch/${PROJECT}/${USER}/projects

  other_cluster:
    url: ssh://cluster.university.edu
    qxub_env: /opt/miniconda/envs/qxub
    platform_file: /shared/configs/slurm_cluster.yaml
    project_root_dir: /home/${USER}/projects
```

## Configuration Fields

### `url`
- **Required**: Connection URL in format `protocol://host:port`
- **SSH format**: `ssh://hostname` or `ssh://hostname:port`
- **Future protocols**: `aws://region/resource`, `k8s://cluster/namespace`
- Host should match an entry in your connection config file

### `config`
- **Optional**: Path to protocol-specific configuration file
- **SSH default**: `~/.ssh/config` (authentication, timeouts, keys)
- **Future K8s**: `~/.kube/config` (cluster credentials, contexts)
- **Future AWS**: `~/.aws/config` (regions, profiles)

### `qxub_env`
- **Required**: Conda environment containing qxub on the remote system
- Can be environment name (if in default conda path) or absolute path
- This environment is activated before running remote qxub commands

### `platform_file`
- **Required**: Absolute path to platform definition file on remote system
- Falls back to searching `XDG_CONFIG_DIRS` if file not found at specified path
- Platform files describe system capabilities (queues, limits, etc.)

### `project_root_dir`
- **Optional**: Template for remote working directories
- Supports environment variable substitution: `${PROJECT}`, `${USER}`
- Used to determine default `--execdir` when not explicitly specified
- Default remote execdir: `{project_root_dir}/{local_project_name}`

## Environment Variable Substitution

The following variables are automatically substituted in configuration values:

- `${USER}`: Current username (from `$USER` environment variable)
- `${PROJECT}`: Project identifier (from `$PROJECT` environment variable)

Example:
```yaml
project_root_dir: /scratch/${PROJECT}/${USER}/projects
# Becomes: /scratch/a56/jr9959/projects (if PROJECT=a56, USER=jr9959)
```

## Working Directory Logic

When you run `qxub --remote system_name`, the remote working directory is determined as:

1. **Explicit `--execdir`**: If specified, use that path
2. **Smart default**: `{project_root_dir}/{local_directory_name}`
3. **Fallback**: Remote home directory

Example:
```bash
# Local directory: ~/my-research-project/
# Remote config: project_root_dir: /scratch/a56/jr9959/projects
# Result: Remote execdir = /scratch/a56/jr9959/projects/my-research-project
```

## Example: NCI Gadi Configuration

For NCI Gadi users, a typical configuration would be:

```yaml
remotes:
  nci_gadi:
    url: ssh://gadi.nci.org.au
    qxub_env: qxub-prod
    platform_file: /g/data/a56/software/qsub_tools/docs/platforms/nci_gadi.yaml
    project_root_dir: /scratch/${PROJECT}/${USER}/projects
```

With corresponding SSH configuration in `~/.ssh/config`:
```ssh
Host gadi.nci.org.au
    HostName gadi.nci.org.au
    User jr9959
    IdentityFile ~/.ssh/id_ed25519  # or id_rsa, id_ecdsa, etc.
    ServerAliveInterval 60
    ServerAliveCountMax 3
    ConnectTimeout 30
```

## Platform File Discovery

Platform files are discovered in this order:

1. **Explicit path**: Use `platform_file` from configuration
2. **XDG fallback**: Search standard configuration directories:
   - `/etc/xdg/qxub/platforms/`
   - `~/.config/qxub/platforms/`
   - `~/.local/share/qxub/platforms/`

## Error Handling

qxub provides friendly error messages for common connection issues:

- **SSH connection failed**: Check SSH configuration and network connectivity
- **Conda environment not found**: Verify `qxub_env` path and conda installation
- **Platform file not found**: Check `platform_file` path and permissions
- **qxub not in environment**: Ensure qxub is installed in the specified conda environment

## Security Considerations

- **No credentials in qxub config**: All authentication handled by SSH
- **Use SSH agent**: For seamless key management
- **SSH key forwarding**: Enable if needed for accessing remote git repositories
- **Connection multiplexing**: Configure in `~/.ssh/config` for faster connections

## Troubleshooting

### Test SSH Connection
```bash
ssh your-remote-host echo "Connection successful"
```

### Verify Remote qxub Installation
```bash
ssh your-remote-host "conda activate your-env && qxub --version"
```

### Check Platform File Access
```bash
ssh your-remote-host "cat /path/to/platform.yaml"
```

### Debug Remote Execution
Use the `--verbose` flag to see detailed execution steps:
```bash
qxub --remote nci_gadi --verbose --env pytorch -- python train.py
```

### Test URL Parsing
For troubleshooting connection URL issues:
```bash
# Test SSH connection using URL format
ssh gadi.nci.org.au echo "SSH via hostname works"

# Verify port if specified in URL
ssh -p 2222 hostname echo "SSH via custom port works"
```
