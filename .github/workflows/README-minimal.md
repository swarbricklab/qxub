# Minimal Remote Execution Workflows

This directory contains minimal workflow examples demonstrating qxub remote execution from GitHub Actions to NCI Gadi.

## Workflows Overview

### 1. `ultra-minimal-ssh-test.yml` - SSH Connection Verification
**Purpose**: Verify basic SSH connectivity to NCI Gadi
**What it does**:
- Loads SSH keys from GCP Secret Manager (via hpci-scripts setup-hpci action)
- Tests SSH connection to Gadi
- Verifies qxub is available on NCI
- No actual job submission

**Use case**: Quick connectivity test, debugging SSH issues

```bash
# Trigger manually
gh workflow run ultra-minimal-ssh-test.yml
```

### 2. `minimal-remote-test.yml` - Basic qxub Remote Execution
**Purpose**: Minimal demonstration of qxub remote execution
**What it does**:
- Installs qxub on GitHub runner
- Configures qxub for remote platform
- Tests dry-run execution
- Tests with conda environment
- Optional: actual job submission (if `ENABLE_ACTUAL_REMOTE_TEST` is true)

**Use case**: Testing qxub remote execution features, CI validation

```bash
# Trigger manually
gh workflow run minimal-remote-test.yml
```

### 3. `hpci-style-remote.yml` - HPCI Pattern Implementation
**Purpose**: Follow the hpci-scripts pattern with qxub
**What it does**:
- Uses same authentication as hpci-scripts workflows
- Accepts custom commands via workflow inputs
- Submits actual jobs to NCI queues
- Demonstrates conda environment execution

**Use case**: Production-like remote execution, custom command testing

```bash
# Trigger with custom command
gh workflow run hpci-style-remote.yml -f command="ls -la && pwd"

# Or just trigger with default
gh workflow run hpci-style-remote.yml
```

## Architecture

All workflows follow this pattern:

```
GitHub Runner (Self-Hosted)
  ↓
Google Cloud Secret Manager
  ↓ (loads SSH keys)
setup-hpci action
  ↓
qxub (local install)
  ↓ (SSH with keys)
NCI Gadi
  ↓
PBS Job Submission
```

## Key Components

### Authentication Setup
Uses `swarbricklab/hpci-scripts/.github/actions/setup-hpci@main`:
- Authenticates to Google Cloud via Workload Identity
- Downloads SSH keys from Secret Manager
- Sets up keys for SSH access to Gadi

### qxub Configuration
Creates `~/.config/qxub/config.yaml`:
```yaml
platforms:
  gadi:
    remote:
      host: <user>@gadi.nci.org.au
      ssh_key: atlas_key  # Loaded by setup-hpci
      working_dir: /scratch/a56/<user>/ci-jobs
      conda_init: |
        eval "$(conda shell.bash hook)"
        conda activate qxub
```

### Required Secrets
These must be configured in GitHub repository settings:

- `SUBMODULE_TOKEN` - GitHub token for submodule access
- `PROJECT_ID` - GCP project ID
- `PROJECT_NUMBER` - GCP project number
- `REGISTRY_SA` - Service account for registry access
- `USER_NAME` - NCI username for SSH connection

### Optional Variables
- `ENABLE_ACTUAL_REMOTE_TEST` - Set to 'true' to enable real job submission

## Comparison with test-remote-execution-secure.yml

| Feature | Minimal Workflows | test-remote-execution-secure.yml |
|---------|------------------|----------------------------------|
| Lines of code | ~60 | ~160 |
| Test coverage | Basic connectivity | Comprehensive features |
| Job submission | Optional | Multiple scenarios |
| Configuration | Hardcoded | Dynamic with secrets |
| Purpose | Demo/debugging | Full CI testing |

## Usage Tips

1. **Start with ultra-minimal**: Verify SSH connectivity first
2. **Progress to minimal**: Test basic qxub remote execution
3. **Use hpci-style**: For production-like testing with custom commands
4. **Full test suite**: Use test-remote-execution-secure.yml for comprehensive validation

## Troubleshooting

### SSH Connection Fails
```bash
# Run ultra-minimal-ssh-test.yml to isolate SSH issues
# Check that atlas_key was created by setup-hpci action
# Verify USER_NAME secret is correct
```

### qxub Not Found on NCI
```bash
# Ensure qxub conda environment exists on NCI
# Check conda_init path in config
# Verify working_dir is accessible
```

### Job Submission Fails
```bash
# Check queue availability: qxub platform list-queues
# Verify storage access: --storage gdata/a56
# Check project allocation: --project a56
```

## Next Steps

After validating these minimal workflows:
1. Customize configuration for your use case
2. Add project-specific test commands
3. Integrate into your CI/CD pipeline
4. Set up monitoring and notifications

## Related Documentation

- `docs/secure-ci-complete.md` - Full secure CI setup guide
- `docs/remote-execution.md` - qxub remote execution documentation
- `docs/platform_configuration.md` - Platform configuration reference
