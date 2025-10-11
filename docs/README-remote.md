# qxub v2.2 Remote Execution Documentation

## 📖 Documentation Overview

This directory contains comprehensive documentation for qxub v2.2's remote execution capabilities. The remote execution feature allows you to submit jobs to remote HPC systems directly from your local development environment.

## 📚 Documentation Structure

### Core Documentation

- **[Remote Configuration Guide](remote-configuration.md)** - How to configure qxub for remote systems
- **[Remote Execution Workflow](remote-execution.md)** - Understanding how remote execution works and CLI usage
- **[SSH Configuration Examples](ssh-configuration.md)** - Optimizing SSH settings for qxub remote execution

### Platform Definitions

- **[nci_gadi.yaml](platforms/nci_gadi.yaml)** - NCI Gadi platform definition (example)

## 🚀 Quick Start

### 1. Configure SSH Access

Set up SSH access to your remote HPC system in `~/.ssh/config`:

```ssh
Host gadi.nci.org.au
    HostName gadi.nci.org.au
    User your_username
    IdentityFile ~/.ssh/id_ed25519  # or id_rsa, id_ecdsa, etc.
    ControlMaster auto
    ControlPath ~/.ssh/master-%r@%h:%p
    ControlPersist 10m
```

### 2. Configure qxub Remote

Create `~/.config/qxub/config.yaml`:

```yaml
remotes:
  nci_gadi:
    url: ssh://gadi.nci.org.au
    qxub_env: qxub-prod
    platform_file: /g/data/a56/software/qsub_tools/docs/platforms/nci_gadi.yaml
    project_root_dir: /scratch/${PROJECT}/${USER}/projects
```

### 3. Submit Remote Jobs

```bash
# Submit job to remote system
qxub --remote nci_gadi --env pytorch -- python train.py

# With resource requirements
qxub --remote nci_gadi --cores 24 --memory 48GB --env tensorflow-gpu -- python model.py
```

## 🏗️ Architecture Overview

### Configuration Separation

- **Platform Definitions**: System capabilities (queues, limits) - stored on remote systems
- **User Configuration**: Connection details and preferences - stored locally (`~/.config/qxub/config.yaml`)
- **Connection Configuration**: Authentication and connection settings - standard config files (`~/.ssh/config`)

### Execution Flow

```
Local Machine                    Remote HPC System
     │                                   │
     │ qxub --remote nci_gadi ...        │
     │ ─────────────────────────────────▶│ SSH connection
     │                                   │ conda activate qxub-env
     │                                   │ qxub --platform-file ...
     │                                   │ ─────────────▶ Job submission
     │ ◀─────────────────────────────────│ Output streaming
```

## 🔧 Key Features

### Simplified Configuration
- **Conda-based setup**: Simple environment activation replaces complex setup scripts
- **SSH delegation**: All connection management handled by standard SSH configuration
- **Smart defaults**: Intelligent working directory resolution based on project structure

### Transparent Execution
- **Same CLI**: Use familiar qxub commands with `--remote` flag
- **Real-time output**: Job output streamed back to local terminal
- **Error handling**: Friendly error messages for common connection and configuration issues

### Flexible Platform Management
- **Shared platforms**: Team-wide platform definitions stored on remote systems
- **Individual preferences**: User-specific settings in local configuration
- **Environment variables**: Support for `${PROJECT}` and `${USER}` substitution

## 📁 File Locations

### Local Configuration
```
~/.config/qxub/config.yaml          # User remote configurations
~/.ssh/config                       # SSH connection settings
~/.ssh/id_rsa                       # SSH private key
```

### Remote Platform Definitions
```
/g/data/a56/software/qsub_tools/docs/platforms/nci_gadi.yaml  # NCI Gadi example
/etc/xdg/qxub/platforms/            # System-wide platform definitions
~/.config/qxub/platforms/           # User-specific platform definitions
```

## 🔍 Troubleshooting

### Common Issues

1. **SSH Connection Failed**
   - Check SSH configuration and network connectivity
   - Test: `ssh your-remote-host echo "test"`

2. **Conda Environment Not Found**
   - Verify `qxub_env` setting in configuration
   - Test: `ssh remote-host "conda activate env-name && qxub --version"`

3. **Platform File Not Found**
   - Check `platform_file` path in remote configuration
   - Verify file exists and is readable on remote system

4. **Permission Denied**
   - Check SSH key permissions: `chmod 600 ~/.ssh/id_rsa`
   - Verify remote directory permissions

### Debug Mode

Use verbose mode for detailed execution information:

```bash
qxub --remote nci_gadi --verbose --env pytorch -- python train.py
```

## 🛡️ Security Considerations

- **No credential storage**: qxub never stores SSH credentials
- **SSH agent integration**: Use SSH agent for secure key management
- **Connection multiplexing**: Efficient connection reuse without compromising security
- **Standard SSH practices**: Leverage proven SSH security features

## 🎯 Best Practices

1. **Version Control**: Keep code synchronized between local and remote systems
2. **Shared Storage**: Leverage shared filesystems when available
3. **Resource Planning**: Use appropriate resource requests for efficient queue usage
4. **Connection Optimization**: Configure SSH multiplexing for better performance
5. **Testing**: Validate jobs locally before remote submission

## 🔗 Related Documentation

- [qxub v2.1 Platform System](../platform-guide.md) - Understanding platform abstractions
- [qxub Environment Management](../environment-guide.md) - Working with conda environments
- [qxub Queue Selection](../queue-guide.md) - Understanding automatic queue selection

## 📋 Next Steps

1. Review the [Remote Configuration Guide](remote-configuration.md) for detailed setup instructions
2. Understand the [Remote Execution Workflow](remote-execution.md) for effective usage
3. Optimize your setup with [SSH Configuration Examples](ssh-configuration.md)
4. Start with simple test jobs before running production workloads

---

**qxub v2.2** - Bringing HPC job submission to your local development workflow
