# Remote Execution Analysis - v3.3.0

## Current State of Remote Execution

### Overview
Remote execution infrastructure exists but is **incomplete and disconnected** from the main CLI. The following analysis details the current implementation status and identifies gaps.

### What Exists

#### 1. **Remote Package Structure** (`qxub/remote/`)
```
qxub/remote/
├── __init__.py          # Package exports
├── config.py            # RemoteConfig dataclass
├── core.py              # Error classes only
├── executor.py          # SSHRemoteExecutor implementation
└── loader.py            # Configuration loading from YAML
```

#### 2. **RemoteConfig** (`qxub/remote/config.py`)
- URL-based configuration (e.g., `ssh://user@host:port`)
- Protocol validation (SSH only currently supported)
- Fields:
  - `name`: Remote identifier
  - `url`: Connection URL
  - `platform`: Optional platform override (remote will auto-detect if not specified)
  - `conda_env`: Optional conda environment (remote will use default if not specified)
  - `working_dir`: Optional working directory
  - `config`: SSH config file path (defaults to `~/.ssh/config`)
  - `force_tty`: TTY allocation control (None=auto, True=force, False=disable)

#### 3. **SSHRemoteExecutor** (`qxub/remote/executor.py`)
- Implements SSH-based remote command execution
- Features:
  - Real-time output streaming
  - Smart TTY allocation (auto-detects if local terminal is interactive)
  - Conda environment activation
  - Platform override via `QXUB_PLATFORM_OVERRIDE`
  - Connection testing
  - Remote qxub availability testing
  - Error handling with suggestions

#### 4. **Configuration Loader** (`qxub/remote/loader.py`)
- Loads remotes from `~/.config/qxub/config.yaml`
- Expected format:
```yaml
remotes:
  gadi:
    url: ssh://gadi.nci.org.au
    platform: nci_gadi_custom  # Optional
    conda_env: pytorch         # Optional
    working_dir: /scratch/a56/{user}/projects  # Optional
```

#### 5. **Documentation** (`docs/remote-execution.md`)
- User-facing documentation describing intended behavior
- Documents auto-detection, SSH setup, usage patterns
- Claims "simple configuration" but implementation incomplete

#### 6. **Tests** (`tests/test_remote_execution.py`)
- Unit tests for RemoteConfig, loader, and executor
- Test SSH command building, connection testing
- Tests appear comprehensive but isolated from CLI

### What's Missing

#### 1. **CLI Integration** ❌
**CRITICAL**: No `--remote` option in any CLI command!
- Not in `qxub/cli.py` (main entry point)
- Not in `qxub/exec_cli.py` (execution command)
- Remote execution cannot be triggered from command line

#### 2. **Execution Flow** ❌
- No integration between `exec_cli` and `RemoteExecutorFactory`
- No logic to detect `--remote` flag and invoke remote execution
- No command forwarding mechanism

#### 3. **Platform Integration** ⚠️
- `qxub/platform/integration.py` exists but references non-existent `execute_remote_platform` function
- No actual implementation of remote platform execution
- Platform auto-detection mentioned in docs but not implemented

#### 4. **Config Manager Integration** ❌
- `ConfigManager` (`qxub/config/manager.py`) doesn't handle remote configurations
- No `get_remote()` or `list_remotes()` methods in config_manager
- Remote config loading completely separate from main config system

#### 5. **Command Building** ⚠️
- No logic to reconstruct full qxub command for remote execution
- Need to serialize execution context (conda/module/singularity)
- Need to serialize PBS options, working directory, etc.

### Architecture Questions

#### 1. **Remote + Platform Relationship**
Current design appears confused:

**Option A: Remotes ARE platforms** (docs suggest this)
```yaml
# Remote auto-detects its platform
remotes:
  gadi:
    url: ssh://gadi.nci.org.au
    # Platform auto-detected as nci_gadi
```

**Option B: Remotes reference platforms** (code suggests this)
```yaml
remotes:
  gadi:
    url: ssh://gadi.nci.org.au
    platform: nci_gadi  # Explicit reference
```

**Current implementation**: Hybrid - supports both but unclear which is primary

#### 2. **Working Directory Resolution**
Documented behavior:
- Default: Mirror local directory structure at remote base
- Example: `/local/my-project` → `/remote/projects/my-project`

Not implemented in current executor!

#### 3. **Execution Context Forwarding**
How should execution contexts be forwarded?

**Current approach** (in executor):
- Conda: Activate via `conda activate {env}` in SSH command
- Platform: Set `QXUB_PLATFORM_OVERRIDE` environment variable

**Missing**:
- Module execution context (`--mod`, `--mods`)
- Singularity execution context (`--sif`, `--bind`)
- Default execution context
- Pre/post commands
- Custom templates

#### 4. **Queue Selection**
Documentation claims "queue selection happens on remote with current, accurate information"

This implies:
- Remote qxub does the queue selection
- `--queue auto` resolved remotely, not locally
- Local machine just forwards the request

But how are resources forwarded? As-is or transformed?

### Configuration Schema Questions

#### Current Schema
```yaml
remotes:
  name:
    url: ssh://host
    platform: optional_platform_name
    conda_env: optional_env_name
    working_dir: optional_path
    config: optional_ssh_config_path
```

#### Alternative: Profiles-based (from dev docs)
```yaml
profiles:
  gadi:
    platform: nci_gadi
    defaults:
      project: a56
      queue: auto
    remote:
      host: gadi.nci.org.au
      user: jr9959
      qxub_command: "module load python3; qxub"
      working_dir: /scratch/a56/jr9959
```

**Question**: Which schema should we implement?

### Test Coverage

#### What's Tested
- RemoteConfig URL parsing and validation ✅
- Configuration loading from YAML ✅
- SSH command building ✅
- Connection testing ✅

#### What's NOT Tested
- CLI integration ❌
- End-to-end remote execution ❌
- Platform auto-detection ❌
- Working directory resolution ❌
- Execution context forwarding ❌

### Integration Points

Where remote execution needs to hook into existing code:

1. **`qxub/exec_cli.py`**: Add `--remote` option
2. **`qxub/execution/context.py`**: Detect remote execution mode
3. **`qxub/execution/unified.py`**: Fork execution path for remote
4. **`qxub/config/manager.py`**: Add remote config methods
5. **`qxub/platform/*.py`**: Add remote platform resolution

### Key Decisions Needed

1. **Configuration Schema**: Stick with current `remotes` or move to `profiles`?
2. **Platform Relationship**: Auto-detect or explicit reference?
3. **Working Directory**: Smart mirroring or simple passthrough?
4. **Execution Context**: Serialize full context or reconstruct?
5. **Queue Selection**: Local or remote resolution?
6. **CLI Structure**: `--remote NAME` or `qxub remote exec --name NAME`?

### Recommended Next Steps

1. **Design Review**: Answer key decisions above
2. **Schema Finalization**: Pick ONE configuration approach
3. **CLI Integration**: Add `--remote` flag to exec_cli
4. **Execution Flow**: Connect exec_cli → RemoteExecutor
5. **Command Serialization**: Build mechanism to reconstruct remote command
6. **Platform Integration**: Implement auto-detection or explicit mapping
7. **End-to-End Testing**: Create real SSH tests
8. **Documentation Update**: Align docs with actual implementation

### Code Quality Notes

- Remote package is well-structured and isolated ✅
- Good separation of concerns (config/executor/loader) ✅
- Comprehensive error handling ✅
- Missing integration is the main gap ❌
- Documentation ahead of implementation (tech debt) ⚠️

### Complexity Assessment

**Low complexity** (~2-3 days):
- Add `--remote` CLI flag
- Basic command forwarding
- Simple platform mapping

**Medium complexity** (~1 week):
- Full execution context serialization
- Working directory mirroring
- Platform auto-detection
- Config system integration

**High complexity** (~2 weeks):
- Profiles-based architecture
- File synchronization
- Queue selection optimization
- Advanced error recovery

### Current Remote Package Quality: 7/10
- Well-designed but disconnected
- Needs integration, not rewrite
- Documentation overpromises capabilities
