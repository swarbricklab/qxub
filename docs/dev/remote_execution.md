# Remote Execution Design

## Overview

This document outlines the design for remote execution capabilities in qxub v2.2, enabling users to submit and monitor jobs on remote platforms from their local machines.

## Architecture

### Client-Server Model

```
┌─────────────────┐    SSH    ┌─────────────────┐    PBS    ┌─────────────────┐
│   Local Client  │ ────────► │  Remote Server  │ ────────► │  Compute Node   │
│   (Laptop)      │           │   (Login Node)  │           │                 │
└─────────────────┘           └─────────────────┘           └─────────────────┘
        │                              │                              │
        │                              │                              │
        ▼                              ▼                              ▼
   Local Logs                   Remote History                  Job Execution
   Exit Codes                   Resource Tracking               Output Files
```

### Execution Flow

1. **Client Process**: User runs qxub on laptop
2. **SSH Connection**: Client connects to remote platform
3. **Remote Invocation**: Client executes qxub on remote login node
4. **Job Submission**: Remote qxub submits job to scheduler
5. **Monitoring**: Remote qxub monitors job, streams output
6. **Output Streaming**: Client receives and displays output
7. **Exit Propagation**: Exit codes propagated: job → remote qxub → client

## Configuration Schema

### Platform Profiles

```yaml
# ~/.config/qxub/config.yaml
profiles:
  gadi:
    platform: nci_gadi           # Platform definition to use
    defaults:                    # Override defaults for this profile
      project: a56
      queue: auto
      walltime: "02:00:00"
    remote:
      host: gadi.nci.org.au      # Remote hostname
      user: jr9959               # Remote username (optional, uses local user)
      port: 22                   # SSH port (optional, default 22)
      qxub_command: "module load python3; qxub"  # Command to run qxub remotely
      working_dir: "/scratch/a56/jr9959"          # Remote working directory
      env_setup: |               # Additional environment setup
        module load python3
        source ~/.bashrc
      sync_files: false          # Whether to sync files (future feature)

  aws_batch:
    platform: aws_batch
    defaults:
      queue: auto
    remote:
      type: cloud               # Cloud execution (future)
      region: us-west-2
```

## Command Syntax

### Profile Selection

```bash
# Execute on specified profile
qxub --profile gadi --env pytorch -- python train.py

# Use profile as default
qxub config set defaults.profile gadi
qxub --env pytorch -- python train.py  # Uses gadi profile

# Override profile settings
qxub --profile gadi --project a99 --env pytorch -- python train.py
```

### Profile Management

```bash
# List available profiles
qxub config profile list

# Create new profile
qxub config profile create myprofile --platform nci_gadi

# Set profile properties
qxub config profile set myprofile.remote.host cluster.university.edu
qxub config profile set myprofile.defaults.project myproject

# Test profile connection
qxub config profile test gadi

# Show profile details
qxub config profile show gadi
```

## Implementation Design

### Core Components

```python
# qxub/remote/executor.py
class RemoteExecutor:
    def __init__(self, profile: Profile):
        self.profile = profile
        self.ssh_client = None

    def execute_remote(self, command_args: List[str]) -> int:
        """Execute qxub command on remote platform"""

    def establish_connection(self) -> bool:
        """Establish SSH connection to remote host"""

    def stream_output(self, remote_process) -> Generator[str, None, None]:
        """Stream output from remote process"""

# qxub/remote/profile.py
class Profile:
    def __init__(self, name: str, config: Dict):
        self.name = name
        self.platform = config['platform']
        self.defaults = config.get('defaults', {})
        self.remote = config.get('remote', {})

    def get_ssh_config(self) -> SSHConfig:
        """Get SSH connection configuration"""

    def get_remote_command(self, local_args: List[str]) -> str:
        """Build remote qxub command from local arguments"""
```

### SSH Connection Management

```python
# qxub/remote/ssh.py
class SSHConnection:
    def __init__(self, config: SSHConfig):
        self.config = config
        self.client = paramiko.SSHClient()

    def connect(self) -> bool:
        """Establish SSH connection using system credentials"""

    def execute_command(self, command: str, stream_output: bool = True):
        """Execute command and optionally stream output"""

    def close(self):
        """Close SSH connection"""
```

### Remote Command Construction

```python
def build_remote_command(profile: Profile, local_args: List[str]) -> str:
    """
    Convert local qxub arguments to remote qxub command

    Example:
    Local:  qxub --profile gadi --env pytorch -- python train.py
    Remote: module load python3; qxub --env pytorch -- python train.py
    """

    # Remove profile-specific arguments
    remote_args = [arg for arg in local_args if not arg.startswith('--profile')]

    # Apply profile defaults
    remote_args = apply_profile_defaults(remote_args, profile.defaults)

    # Build remote command
    setup_command = profile.remote.get('env_setup', '')
    qxub_command = profile.remote.get('qxub_command', 'qxub')

    return f"{setup_command}; {qxub_command} {' '.join(remote_args)}"
```

## Output Streaming

### Real-time Output Display

```python
class OutputStreamer:
    def __init__(self, ssh_session):
        self.session = ssh_session
        self.stdout_thread = None
        self.stderr_thread = None

    def start_streaming(self):
        """Start threads to stream stdout and stderr"""
        self.stdout_thread = threading.Thread(
            target=self._stream_channel,
            args=(self.session.stdout, sys.stdout)
        )
        self.stderr_thread = threading.Thread(
            target=self._stream_channel,
            args=(self.session.stderr, sys.stderr)
        )

    def _stream_channel(self, source, destination):
        """Stream data from source channel to destination"""
        while True:
            data = source.read(1024)
            if not data:
                break
            destination.write(data.decode())
            destination.flush()
```

## Logging and History

### Distributed Logging Strategy

```yaml
# Local execution log
local_history:
  timestamp: "2025-10-09T14:30:00"
  command: "qxub --profile gadi --env pytorch -- python train.py"
  profile: gadi
  remote_host: gadi.nci.org.au
  remote_command: "module load python3; qxub --env pytorch -- python train.py"
  status: success
  exit_code: 0

# Remote execution log (on gadi)
remote_history:
  timestamp: "2025-10-09T14:30:05"
  command: "qxub --env pytorch -- python train.py"
  platform: nci_gadi
  queue: gpuvolta
  job_id: "47291842.gadi-pbs"
  status: success
  exit_code: 0
  resources_used:
    walltime: "01:23:45"
    memory: "16GB"
    cpus: 12
```

### History Synchronization

```python
class DistributedHistory:
    def log_remote_execution(self, profile: Profile, local_command: str,
                           remote_command: str, result: ExecutionResult):
        """Log execution on both local and remote systems"""

        # Log locally
        local_history.log_execution(
            command=local_command,
            profile=profile.name,
            remote_host=profile.remote.host,
            result=result
        )

        # Remote logging handled by remote qxub instance
```

## Error Handling

### Connection Failures

```python
class RemoteExecutionError(Exception):
    pass

class ConnectionError(RemoteExecutionError):
    pass

class AuthenticationError(RemoteExecutionError):
    pass

def handle_connection_errors(func):
    """Decorator for handling SSH connection errors"""
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except paramiko.AuthenticationException:
            raise AuthenticationError("SSH authentication failed")
        except paramiko.SSHException as e:
            raise ConnectionError(f"SSH connection failed: {e}")
        except socket.error as e:
            raise ConnectionError(f"Network error: {e}")
    return wrapper
```

### Graceful Degradation

```bash
# If remote connection fails, offer alternatives
$ qxub --profile gadi --env pytorch -- python train.py
Error: Cannot connect to gadi.nci.org.au
Suggestions:
  1. Check VPN connection to NCI
  2. Verify SSH key authentication: ssh gadi.nci.org.au
  3. Run locally: qxub --env pytorch -- python train.py
  4. Use different profile: qxub --profile local --env pytorch -- python train.py
```

## Security Considerations

### Credential Management

- **SSH Key Authentication**: Rely on system SSH configuration
- **Agent Forwarding**: Support for SSH agent forwarding
- **No Password Storage**: Never store passwords in configuration
- **Connection Validation**: Verify host keys and certificates

### Platform Access Control

```yaml
# Security-conscious profile configuration
profiles:
  gadi:
    remote:
      host: gadi.nci.org.au
      strict_host_key_checking: true    # Verify host keys
      connect_timeout: 30               # Connection timeout
      compression: true                 # Compress SSH traffic
      cipher: "aes256-ctr"             # Specify encryption cipher
```

## Testing Strategy

### Integration Tests

```python
class TestRemoteExecution:
    def test_ssh_connection(self):
        """Test SSH connection establishment"""

    def test_command_execution(self):
        """Test remote command execution"""

    def test_output_streaming(self):
        """Test real-time output streaming"""

    def test_exit_code_propagation(self):
        """Test exit code propagation"""

    def test_connection_failure_handling(self):
        """Test handling of connection failures"""
```

### Mock Remote Environment

```python
class MockRemoteExecutor:
    """Mock remote executor for testing without actual SSH connections"""

    def execute_remote(self, command_args):
        # Simulate remote execution locally
        return subprocess.run(command_args).returncode
```

## Future Enhancements

### File Synchronization

```yaml
profiles:
  gadi:
    remote:
      sync_files:
        enabled: true
        local_dir: "."
        remote_dir: "/scratch/a56/jr9959/work"
        exclude: [".git/", "__pycache__/", "*.pyc"]
        method: "rsync"  # or "scp", "sftp"
```

### Connection Pooling

```python
class ConnectionPool:
    """Maintain persistent SSH connections for repeated use"""

    def get_connection(self, profile: Profile) -> SSHConnection:
        """Get existing or create new connection"""

    def cleanup_idle_connections(self):
        """Close connections idle for too long"""
```

### Advanced Monitoring

```yaml
remote:
  monitoring:
    job_status_interval: 30s      # How often to check job status
    output_buffer_size: 8192      # Output streaming buffer size
    heartbeat_interval: 60s       # Connection heartbeat
    max_connection_retries: 3     # Retry failed connections
```
