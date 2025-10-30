# Repository Security Configuration

This document outlines the security configuration requirements for the qxub repository.

## Branch Protection Rules

The following branch protection rules should be configured via GitHub repository settings:

### Main Branch (`main`)
- **Require a pull request before merging**: ✅ Enabled
- **Require approvals**: 1 approval minimum
- **Dismiss stale PR approvals when new commits are pushed**: ✅ Enabled
- **Require review from code owners**: ✅ Enabled
- **Restrict pushes that create files > 100MB**: ✅ Enabled
- **Require status checks to pass before merging**: ✅ Enabled
  - Required checks:
    - `format-check` (Code Formatting)
    - `pylint` (Pylint)
    - `analyze` (CodeQL)
    - `test-remote-execution-secure` (if available)
- **Require branches to be up to date before merging**: ✅ Enabled
- **Require conversation resolution before merging**: ✅ Enabled
- **Restrict pushes from users who can push to matching branches**: ✅ Enabled
- **Do not allow bypassing the above settings**: ✅ Enabled
- **Allow force pushes**: ❌ Disabled
- **Allow deletions**: ❌ Disabled

### Feature Branches (`feature/*`)
- **Require a pull request before merging**: ✅ Enabled
- **Require status checks to pass before merging**: ✅ Enabled
  - Required checks:
    - `format-check` (Code Formatting)
    - `pylint` (Pylint)
    - `analyze` (CodeQL)

## Repository Security Settings

### General Security
- **Private vulnerability reporting**: ✅ Enabled
- **Dependency graph**: ✅ Enabled
- **Dependabot alerts**: ✅ Enabled
- **Dependabot security updates**: ✅ Enabled
- **Dependabot version updates**: ✅ Enabled (configured via `.github/dependabot.yml`)

### Code Security and Analysis
- **Code scanning**: ✅ Enabled (CodeQL)
- **Secret scanning**: ✅ Enabled
- **Secret scanning push protection**: ✅ Enabled
- **Token scanning**: ✅ Enabled

### Actions Security
- **Actions permissions**:
  - Allow select actions and reusable workflows
  - Allow actions created by GitHub: ✅ Enabled
  - Allow verified marketplace actions by verified creators: ✅ Enabled
  - Specified actions: Pin to specific SHA hashes
- **Fork pull request workflows**: Require approval for first-time contributors
- **Fork pull request workflows in private repositories**: Require approval for all outside contributors

## Secrets Management

### Repository Secrets (Legacy - Migrate to GCP Secret Manager)
- `GADI_SSH_PRIVATE_KEY`: SSH private key for Gadi access
- `GADI_SSH_KNOWN_HOSTS`: Known hosts for SSH verification
- `NCI_USERNAME`: NCI username for remote access
- `NCI_PROJECT`: NCI project identifier
- `ENABLE_ACTUAL_REMOTE_TEST`: Boolean flag for real testing

### Google Cloud Secrets (Preferred)
- `PROJECT_ID`: GCP Project ID for Workload Identity
- `PROJECT_NUMBER`: GCP Project Number for Workload Identity
- Secrets stored in GCP Secret Manager:
  - `gadi_ssh_private_key`: SSH private key for Gadi
  - `gadi_ssh_known_hosts`: SSH known hosts

### Environment Variables
- `ENABLE_ACTUAL_REMOTE_TEST`: Repository variable for conditional testing

## Security Monitoring

### Required Status Checks
1. **Code Formatting** (`format-check`): Ensures consistent code style
2. **Pylint** (`pylint`): Static code analysis for Python
3. **CodeQL** (`analyze`): Security vulnerability scanning
4. **Remote Execution Tests** (`test-remote-execution-secure`): Integration testing

### Security Alerts
- **Dependabot alerts**: Monitor for vulnerable dependencies
- **Code scanning alerts**: Monitor for security vulnerabilities
- **Secret scanning alerts**: Monitor for exposed secrets

## Workflow Security

### Required Permissions
All workflows must specify minimal required permissions:
```yaml
permissions:
  contents: read  # Default minimum
  # Add only additional permissions as needed
```

### Action Pinning
All GitHub Actions must be pinned to specific SHA hashes:
```yaml
- uses: actions/checkout@692973e3d937129bcbf40652eb9f2f61becf3332 # v4.1.7
```

### Secret Access
- Use GCP Secret Manager with Workload Identity Federation (preferred)
- Avoid long-lived service account keys
- Minimize secret access to necessary workflows only

## Compliance Checklist

- [ ] Branch protection rules configured for `main` branch
- [ ] Required status checks configured
- [ ] Dependabot configured and enabled
- [ ] CodeQL security scanning enabled
- [ ] Secret scanning enabled with push protection
- [ ] All workflows use minimal permissions
- [ ] All actions pinned to SHA hashes
- [ ] SECURITY.md policy published
- [ ] Workload Identity Federation configured
- [ ] SSH keys migrated to GCP Secret Manager
- [ ] Legacy secrets marked for rotation/removal

## Security Contacts

For security-related configuration questions:
- Repository maintainers via institutional channels
- Security team via appropriate escalation paths

## Review Schedule

This security configuration should be reviewed:
- **Quarterly**: Complete security posture review
- **On incident**: Immediate review after any security event
- **On major changes**: Review when significant features are added
- **Dependency updates**: Review when major dependencies change
