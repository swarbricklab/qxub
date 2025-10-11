# SSH Configuration for qxub Remote Execution

## Overview

qxub delegates all connection management to SSH, allowing you to leverage SSH's robust authentication, connection multiplexing, and security features. This document provides example SSH configurations for optimal qxub remote execution.

## Basic SSH Configuration

Create or edit `~/.ssh/config` to configure remote hosts:

```ssh
# NCI Gadi Configuration
Host gadi.nci.org.au
    HostName gadi.nci.org.au
    User jr9959
    IdentityFile ~/.ssh/id_ed25519  # Modern default (or id_rsa, id_ecdsa)
    Port 22

    # Connection timeouts
    ConnectTimeout 30
    ServerAliveInterval 60
    ServerAliveCountMax 3

    # Compression for better performance over slow connections
    Compression yes

    # Connection multiplexing for faster subsequent connections
    ControlMaster auto
    ControlPath ~/.ssh/master-%r@%h:%p
    ControlPersist 10m
```

## Advanced SSH Configuration

### Connection Multiplexing

Connection multiplexing allows multiple SSH sessions to share a single connection, significantly improving performance for repeated qxub remote executions:

```ssh
Host *.nci.org.au
    ControlMaster auto
    ControlPath ~/.ssh/master-%r@%h:%p
    ControlPersist 10m

    # Reuse existing connections
    ControlMaster auto

    # Keep connections alive for 10 minutes after last use
    ControlPersist 10m

    # Connection health monitoring
    ServerAliveInterval 60
    ServerAliveCountMax 3
```

### Proxy/Jump Hosts

If you need to go through a bastion/jump host:

```ssh
# Jump host configuration
Host bastion.university.edu
    HostName bastion.university.edu
    User myuser
    IdentityFile ~/.ssh/bastion_key

# Target host via jump host
Host cluster.internal
    HostName cluster.internal
    User myuser
    ProxyJump bastion.university.edu
    IdentityFile ~/.ssh/cluster_key
```

Corresponding qxub configuration:
```yaml
remotes:
  university_cluster:
    ssh_host: cluster.internal  # SSH handles the proxy jump
    qxub_env: qxub
    platform_file: /shared/qxub/platforms/slurm.yaml
    project_root_dir: /home/${USER}/projects
```

### SSH Agent Integration

For seamless key management, use SSH agent:

```bash
# Start SSH agent (add to ~/.bashrc or ~/.zshrc)
eval "$(ssh-agent -s)"

# Add your private keys (various key types)
ssh-add ~/.ssh/id_ed25519      # Modern Ed25519 key (recommended)
ssh-add ~/.ssh/id_rsa          # RSA key (legacy but common)
ssh-add ~/.ssh/id_ecdsa        # ECDSA key (good alternative)

# Verify loaded keys
ssh-add -l
```

SSH configuration for agent usage:
```ssh
Host *
    # Use SSH agent for key management
    AddKeysToAgent yes

    # Forward SSH agent to remote hosts (if needed)
    ForwardAgent no  # Enable only if you need git access on remote
```

### Key-specific Configuration

Different keys for different hosts:

```ssh
# NCI Gadi with specific key
Host gadi.nci.org.au
    HostName gadi.nci.org.au
    User jr9959
    IdentityFile ~/.ssh/nci_gadi_ed25519  # Modern Ed25519 key
    IdentitiesOnly yes  # Only use specified key

# University cluster with different key type
Host cluster.university.edu
    HostName cluster.university.edu
    User researcher
    IdentityFile ~/.ssh/university_rsa  # Legacy RSA key if required
    IdentitiesOnly yes

# AWS EC2 instance with ECDSA key
Host ec2-instance
    HostName ec2-12-34-56-78.compute-1.amazonaws.com
    User ec2-user
    IdentityFile ~/.ssh/aws_ecdsa_key  # ECDSA key
    IdentitiesOnly yes
```

### Performance Optimization

```ssh
# Global optimizations
Host *
    # Compression for slower connections
    Compression yes

    # Faster ciphers (if supported by remote)
    Ciphers aes128-gcm@openssh.com,aes256-gcm@openssh.com,chacha20-poly1305@openssh.com

    # Disable features not needed for qxub
    X11Forwarding no
    ForwardAgent no

    # Keep connections alive
    ServerAliveInterval 60
    ServerAliveCountMax 3
```

### Security Hardening

```ssh
# Security-focused configuration
Host production.hpc.edu
    HostName production.hpc.edu
    User scientist

    # Use only specific authentication methods
    PreferredAuthentications publickey
    PasswordAuthentication no
    ChallengeResponseAuthentication no

    # Strict host key checking
    StrictHostKeyChecking yes

    # Specific key only
    IdentityFile ~/.ssh/production_key
    IdentitiesOnly yes

    # No forwarding
    X11Forwarding no
    ForwardAgent no
    AllowTcpForwarding no

    # Connection limits
    ConnectTimeout 10
    ServerAliveInterval 30
    ServerAliveCountMax 2
```

## Multiple Host Patterns

### Wildcard Matching

```ssh
# Apply settings to all NCI hosts
Host *.nci.org.au
    User jr9959
    IdentityFile ~/.ssh/nci_key
    ControlMaster auto
    ControlPath ~/.ssh/master-nci-%r@%h:%p
    ControlPersist 10m

# Apply settings to all university hosts
Host *.university.edu
    User researcher
    IdentityFile ~/.ssh/university_key
    ProxyJump bastion.university.edu

# Apply to all development hosts
Host dev-*
    User developer
    IdentityFile ~/.ssh/dev_key
    StrictHostKeyChecking no  # For dynamic dev environments
```

### Environment-specific Configuration

```ssh
# Production environment
Host prod-*.company.com
    User prod_user
    IdentityFile ~/.ssh/prod_key
    StrictHostKeyChecking yes
    ConnectTimeout 5

# Development environment
Host dev-*.company.com
    User dev_user
    IdentityFile ~/.ssh/dev_key
    StrictHostKeyChecking no
    ConnectTimeout 30

# Testing environment
Host test-*.company.com
    User test_user
    IdentityFile ~/.ssh/test_key
    StrictHostKeyChecking ask
```

## Troubleshooting SSH Issues

### Connection Debugging

Enable verbose SSH output to diagnose connection issues:

```bash
# Test connection with verbose output
ssh -v gadi.nci.org.au echo "Connection test"

# Even more verbose
ssh -vvv gadi.nci.org.au echo "Connection test"
```

### Common SSH Configuration Issues

**Key Permission Problems**
```bash
# Fix key permissions (applies to all key types)
chmod 600 ~/.ssh/id_ed25519 ~/.ssh/id_rsa ~/.ssh/id_ecdsa
chmod 644 ~/.ssh/id_ed25519.pub ~/.ssh/id_rsa.pub ~/.ssh/id_ecdsa.pub
chmod 700 ~/.ssh
chmod 600 ~/.ssh/config
```

**Host Key Verification**
```bash
# Remove old host key if changed
ssh-keygen -R gadi.nci.org.au

# Add new host key
ssh-keyscan -H gadi.nci.org.au >> ~/.ssh/known_hosts
```

**Connection Multiplexing Issues**
```bash
# Check existing control connections
ls ~/.ssh/master-*

# Kill stuck control connections
ssh -O exit gadi.nci.org.au
```

### Testing SSH Configuration

```bash
# Test basic connection
ssh gadi.nci.org.au echo "Basic connection works"

# Test conda activation
ssh gadi.nci.org.au "conda activate qxub-prod && echo 'Conda works'"

# Test qxub availability
ssh gadi.nci.org.au "conda activate qxub-prod && qxub --version"

# Test file access
ssh gadi.nci.org.au "ls /g/data/a56/software/qsub_tools/docs/platforms/"
```

## Integration with qxub

### Recommended Configuration Template

For most qxub users, this SSH configuration provides good performance and security:

```ssh
# Template for qxub remote execution
Host YOUR_HPC_HOST
    HostName your.hpc.system.edu
    User your_username
    IdentityFile ~/.ssh/your_key  # Use ed25519, rsa, or ecdsa as appropriate

    # Performance optimization
    ControlMaster auto
    ControlPath ~/.ssh/master-%r@%h:%p
    ControlPersist 10m
    Compression yes

    # Connection reliability
    ServerAliveInterval 60
    ServerAliveCountMax 3
    ConnectTimeout 30

    # Security
    PasswordAuthentication no
    X11Forwarding no
    ForwardAgent no
```

### Environment Variables in SSH

SSH can use environment variables in configuration:

```ssh
Host gadi.nci.org.au
    HostName gadi.nci.org.au
    User %r  # Uses local username
    IdentityFile ~/.ssh/gadi_%r_key
```

However, for qxub remote configuration, environment variable substitution happens in the qxub config file, not SSH config.

## Security Best Practices

1. **Use SSH Agent**: Avoid storing unencrypted private keys
2. **Specific Keys**: Use different keys for different purposes
3. **Regular Rotation**: Rotate SSH keys periodically
4. **Monitor Access**: Check SSH logs for unauthorized access attempts
5. **Disable Password Auth**: Use key-based authentication only
6. **Connection Limits**: Set appropriate timeouts and retry limits
7. **Host Key Verification**: Always verify host keys on first connection

## Performance Best Practices

1. **Connection Multiplexing**: Reuse connections for better performance
2. **Compression**: Enable for slow networks
3. **Keep-Alive**: Prevent connection drops during long jobs
4. **Fast Ciphers**: Use modern, efficient encryption
5. **Minimal Forwarding**: Disable unnecessary features

This SSH configuration works seamlessly with qxub's remote execution, providing secure, efficient, and reliable connections to your HPC systems.
