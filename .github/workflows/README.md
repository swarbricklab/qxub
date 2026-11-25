# GitHub Actions Workflows

This directory contains CI workflows for qxub testing and release automation.

## Branching Model & CI Tiers

Feature/fix branches open PRs **directly against `main`** (see `.github/copilot-instructions.md`). CI is tiered:

| CI tier | Workflows | Trigger | Runner |
|---------|-----------|---------|--------|
| **Lightweight** | `formatting.yml`, `pylint.yml` | Every push & PR | GitHub-hosted |
| **Full** | `runner-connection-test.yml`, `runner-validation.yml` | Push to `main`, or manual dispatch | Self-hosted |
| **Release** | `release.yml` | Push to `main` (version bump) | GitHub-hosted |

**Day-to-day**: PRs into `main` run the GitHub-hosted checks only. The self-hosted runner tests run after merge.
**Release**: merging into `main` creates a GitHub Release when the version is bumped.

The self-hosted runner has SSH access to Gadi, so its workflows must never run on
`pull_request`: a pull request supplies its own copy of the workflow file. To try
the runner tests on a branch, a maintainer can run them with `workflow_dispatch`
(this stops working once the runner group is limited to `refs/heads/main`).

## Active Workflows

### Code Quality (GitHub Runners)
- **`formatting.yml`** - Black/isort code formatting validation
- **`pylint.yml`** - Python linting checks
- **`release.yml`** - GitHub Release on a version bump to `main`. The wheel is built in a read-only job; only the publish job can write.

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
Both runner workflows run on push to `main`.

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
