# GitHub Actions Workflows

This directory contains CI workflows for qxub testing and release automation.

## Active Workflows

### Code Quality (GitHub Runners)
- **`formatting.yml`** - Black/isort code formatting validation
- **`pylint.yml`** - Python linting checks
- **`release.yml`** - Automated PyPI release on version tags

### qxub Testing (Self-Hosted Runners)
- **`runner-connection-test.yml`** - Progressive connection test (dry-run → live execution)
- **`runner-validation.yml`** - Comprehensive environment and configuration validation

## Self-Hosted Runner Architecture

The qxub testing workflows use self-hosted GitHub Actions runners with pre-installed qxub:

- **Platform**: Self-hosted runners (not GitHub-hosted)
- **qxub**: Pre-installed system-wide (v3.3.0+)
- **Configuration**: `~/.config/qxub/config.yaml` with remote platform configured
- **SSH**: Pre-configured in `~/.ssh/config` and `~/.ssh/known_hosts`
- **Execution Mode**: REMOTE_DELEGATED (runner → SSH → Gadi → PBS)
- **Working Directory**: `/scratch/a56/$USER/ci` on Gadi

### No Authentication Setup Required

Unlike previous experimental workflows, these workflows:
- ❌ Don't use WIF (Workload Identity Federation)
- ❌ Don't generate SSH keys dynamically
- ❌ Don't use custom container images
- ✅ Use pre-configured qxub installation
- ✅ Use pre-configured SSH access
- ✅ Execute qxub commands directly

## Running Workflows

### Manual Trigger
```bash
# Test runner connection and qxub execution
gh workflow run runner-connection-test.yml

# Validate complete runner environment
gh workflow run runner-validation.yml
```

### Automatic Trigger
Both runner workflows automatically trigger on push to `feature/v3.3.0-remote-execution`.

## Documentation

For runner setup and operational details, see:
- [`docs/dev/qxub-runner-strategy.md`](../../docs/dev/qxub-runner-strategy.md) - Runner configuration guide
- [`docs/dev/qxub-ci-strategy.md`](../../docs/dev/qxub-ci-strategy.md) - CI architecture and troubleshooting

## Workflow Evolution

**Previous experimental approaches** (removed Nov 2025):
- WIF authentication workflows
- Container-based execution
- Dynamic SSH key generation
- hpci-scripts integration patterns

**Current stable approach** (Nov 2025+):
- Self-hosted runners with pre-configured qxub
- Direct SSH access (no dynamic auth)
- REMOTE_DELEGATED execution mode
- Minimal workflow complexity
