# Platform vs Config: Separation of Concerns

## Core Principle

**Platform Definitions** and **User Configuration** are fundamentally different types of information and must be kept separate.

### Platform Definitions
- **What**: Immutable facts about a platform (queues, limits, billing rates)
- **Who**: Same for all users
- **Where**: Centralized, version-controlled, potentially remote
- **Format**: YAML platform definition files (existing format is good)
- **Examples**: `nci_gadi.yaml`, `pawsey_setonix.yaml`

### User Configuration
- **What**: User choices and preferences
- **Who**: Different per user/admin/system
- **Where**: `~/.config/qxub/config.yaml` (user), `/etc/xdg/qxub/config.yaml` (system)
- **Format**: YAML config files
- **Examples**: Which platform to use, SSH credentials, default resources

## Architecture Design

### Platform Definition Files (Unchanged)

These stay exactly as they are:

```yaml
# nci_gadi.yaml - Platform definition
platform:
  name: nci_gadi
  type: pbs_pro
  host: gadi.nci.org.au
  description: "National Computational Infrastructure - Gadi supercomputer"

  queues:
    normal:
      type: standard
      limits:
        max_cores: 48
        max_memory_gb: 192
        max_walltime: "48:00:00"
      su_billing_rate: 1.0

    gpu:
      type: gpu
      limits:
        max_gpus: 4
        max_walltime: "24:00:00"
      su_billing_rate: 3.0

  auto_selection:
    - condition: "gpus > 0"
      queue: gpu
    - condition: "cores <= 48"
      queue: normal
      is_default: true
```

### User Configuration Files (New Structure)

```yaml
# ~/.config/qxub/config.yaml

# Default platform selection
defaults:
  platform: nci_gadi  # Reference to platform name

# Platform registry - tells qxub where to find platform definitions
platforms:
  nci_gadi:
    # 1. Platform name (used in --platform flag)
    name: nci_gadi

    # 2. Access URL / hostname (for remote execution)
    access:
      url: ssh://gadi.nci.org.au
      # Optional SSH details
      port: 22
      user: jr9959
      ssh_config: ~/.ssh/config  # Use SSH config for connection details

    # 3. Platform definition location
    definition: file:///g/data/a56/software/qxub/platforms/nci_gadi.yaml
    # Or: definition: https://github.com/swarbricklab/qxub-platforms/raw/main/nci_gadi.yaml
    # Or: definition: ssh://gadi-dm.nci.org.au/apps/qxub/platforms/nci_gadi.yaml
    # Or: definition: nci_gadi.yaml  # Defaults to file:// in standard locations

    # Optional: Remote execution settings
    remote:
      working_dir: /scratch/a56/{user}
      conda_init: |
        eval "$(conda shell.bash hook)"

  pawsey_setonix:
    name: pawsey_setonix
    access:
      url: ssh://setonix.pawsey.org.au
    definition: https://github.com/swarbricklab/qxub-platforms/raw/main/pawsey_setonix.yaml
    remote:
      working_dir: /scratch/{project}/{user}

  local_cluster:
    name: local_cluster
    # No access URL = local execution only
    definition: file:///etc/qxub/platforms/local.yaml

# Other config (defaults, aliases, etc.)
defaults:
  project: a56
  walltime: "1:00:00"
  mem: 4GB
```

## Platform Definition Resolution

### URL Schemes Supported

1. **Local file** (default):
   ```yaml
   definition: nci_gadi.yaml
   # Searches in order:
   # 1. ~/.config/qxub/platforms/nci_gadi.yaml
   # 2. /etc/xdg/qxub/platforms/nci_gadi.yaml
   # 3. $QXUB_PLATFORM_PATHS/nci_gadi.yaml
   ```

2. **Explicit local file**:
   ```yaml
   definition: file:///g/data/a56/software/qxub/platforms/nci_gadi.yaml
   ```

3. **HTTPS URL**:
   ```yaml
   definition: https://github.com/swarbricklab/qxub-platforms/raw/main/nci_gadi.yaml
   # Fetched once, cached locally
   ```

4. **SSH URL** (fetch from remote system):
   ```yaml
   definition: ssh://gadi-dm.nci.org.au/apps/qxub/platforms/nci_gadi.yaml
   # Useful for accessing platform files on data mover nodes
   ```

### Caching Strategy

```
~/.cache/qxub/platforms/
  ├── nci_gadi.yaml          # Cached from HTTPS
  ├── nci_gadi.yaml.meta     # Cache metadata (source URL, timestamp, etag)
  └── pawsey_setonix.yaml
```

**Cache invalidation**:
- TTL: 24 hours by default (configurable)
- Force refresh: `qxub platform refresh nci_gadi`
- Clear all: `qxub platform clear-cache`
- Auto-refresh: Check ETag/Last-Modified headers for HTTP URLs

## CLI Usage

### Specifying Platform

```bash
# Use default platform from config
qxub exec --env pytorch -- python train.py

# Explicitly specify platform
qxub exec --platform nci_gadi --env pytorch -- python train.py

# Override at system level
export QXUB_PLATFORM=nci_gadi
qxub exec --env pytorch -- python train.py
```

### Platform Management Commands

```bash
# List available platforms
qxub platform list
# Output:
#   nci_gadi (remote) - National Computational Infrastructure
#   pawsey_setonix (remote) - Pawsey Supercomputing Centre
#   local_cluster (local) - Local PBS cluster

# Show platform details
qxub platform show nci_gadi
# Output: Full platform definition + access config

# Validate platform definition
qxub platform validate nci_gadi

# Refresh cached platform definitions
qxub platform refresh nci_gadi
qxub platform refresh --all

# Test platform connectivity (for remote platforms)
qxub platform test nci_gadi
```

## Execution Flow

### Local Execution
```
User: qxub exec --platform local_cluster --env pytorch -- python train.py
  ↓
Load config: defaults.platform = local_cluster
  ↓
Load platform def: file:///etc/qxub/platforms/local.yaml
  ↓
No access.url → Local execution
  ↓
Queue selection using local platform definition
  ↓
Submit job locally
```

### Remote Execution
```
User: qxub exec --platform nci_gadi --env pytorch -- python train.py
  ↓
Load config: platforms.nci_gadi.access.url = ssh://gadi.nci.org.au
  ↓
Load platform def: Cache or fetch from definition URL
  ↓
access.url present → Remote execution
  ↓
Build SSH command with forwarded arguments
  ↓
SSH to gadi.nci.org.au
  ↓
Remote qxub loads ITS platform definition (from remote's config/defaults)
  ↓
Remote qxub submits job
```

### Key Insight: Platform Definition Used Locally First

For remote execution, the **local** machine loads the platform definition to:
1. Validate resources before SSH connection
2. Provide better error messages locally
3. Determine if platform supports requested features
4. Build appropriate SSH command

Then the **remote** machine also loads the platform definition to:
1. Perform queue auto-selection with current information
2. Submit the job with accurate constraints
3. Use platform-specific job script templates

## Benefits of This Design

### 1. Clean Separation
- Platform facts separate from user choices
- Different maintenance responsibilities
- Different update frequencies

### 2. Flexibility
- Central platform definition repository
- Local overrides possible
- Mixed local/remote platforms

### 3. Consistency
- Same platform definition for all users
- Version-controlled platform definitions
- Easy platform updates (just update URL)

### 4. Scalability
- Add new platforms by adding config entry
- No code changes needed for new platforms
- Community can maintain platform definitions

## Migration Path

### Current State
- Platform definitions in `qxub/platform/definitions/`
- No separate platform registry
- No remote execution

### v3.3.0 Changes
1. Add `platforms:` section to config schema
2. Keep existing platform definitions unchanged
3. Add platform loader with URL support
4. Add caching layer
5. Add `--platform` CLI option
6. Trigger remote execution if `access.url` present

### Backward Compatibility
```yaml
# Old style (still works)
defaults:
  platform: nci_gadi  # Uses built-in definition

# New style (recommended)
platforms:
  nci_gadi:
    name: nci_gadi
    definition: https://github.com/swarbricklab/qxub-platforms/raw/main/nci_gadi.yaml
defaults:
  platform: nci_gadi
```

## Open Questions

1. **Queue selection**: Should local qxub do queue auto-selection before SSH, or always forward to remote?
   - **Proposal**: Local does validation/selection, remote can override if needed

2. **Platform definition versioning**: How to handle breaking changes in platform definitions?
   - **Proposal**: Version in filename? `nci_gadi_v2.yaml`

3. **Platform aliases**: Should `gadi` be a shortcut for `nci_gadi`?
   - **Proposal**: Support aliases in config

4. **SSH key management**: Where to configure SSH keys?
   - **Proposal**: Reference `~/.ssh/config` entries, don't duplicate SSH config

5. **Definition conflicts**: What if local and remote have different platform definitions?
   - **Proposal**: Remote's definition wins for job submission, local used for pre-validation only

## Implementation Checklist

- [ ] Add `platforms:` config schema
- [ ] Implement platform loader with URL support (file://, https://, ssh://)
- [ ] Implement caching layer
- [ ] Add cache management commands
- [ ] Add `--platform` CLI option
- [ ] Integrate with remote executor
- [ ] Update config documentation
- [ ] Create example platform repository
- [ ] Migration guide for existing configs
