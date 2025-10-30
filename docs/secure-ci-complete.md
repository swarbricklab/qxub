# Secure CI Setup for qxub

Complete guide for setting up secure CI testing using self-hosted GitHub runners with GCP Secret Manager.

## Security Architecture

**Traditional GitHub Secrets:**
- SSH keys stored as GitHub repository secrets
- Keys exposed in runner environment during execution
- Higher risk if runner is compromised

**Self-Hosted + GCP Secrets (recommended):**
- SSH keys stored in GCP Secret Manager
- Workload Identity Federation (no long-lived keys)
- Keys only loaded when needed and cleaned up after use
- Much lower risk exposure

## Current Status

### ✅ Already Implemented
- **Secure workflow**: `.github/workflows/test-remote-execution-secure.yml`
- **Setup script**: `scripts/setup-secure-ci.sh`
- **Self-hosted runner support**: Uses `runs-on: self-hosted`
- **Workload Identity**: `qxub-workload-identity-pool`
- **Service Account**: `qxub-ci-sa@PROJECT_ID.iam.gserviceaccount.com`

### 🔧 Still Needed
1. **GCP Project Setup** - Enable APIs and create workload identity
2. **Create GCP Secrets** - Store SSH keys securely
3. **Deploy Self-Hosted Runner** - Configure runner with GCP access
4. **Test Secure Workflow** - Verify end-to-end functionality

## Quick Setup

### 1. Run the Setup Script
```bash
# This creates GCP secrets and configures permissions
./scripts/setup-secure-ci.sh
```

### 2. Required GCP Secrets
- `gadi_ssh_private_key` - SSH private key for Gadi access
- `gadi_ssh_known_hosts` - Known hosts for gadi.nci.org.au
- `nci_username` - NCI username for testing
- `nci_project` - NCI project code for testing

### 3. Self-Hosted Runner Setup
The runner needs:
- **Python 3.10+** installed
- **gcloud CLI** installed and authenticated
- **Access to GCP Secret Manager**

### 4. Enable Secure Workflow
The secure workflow will automatically run when:
- Self-hosted runners are available
- GCP secrets are configured
- Workload identity is set up

## Security Features

- **No long-lived keys** in GitHub secrets
- **Automatic cleanup** of SSH keys after use
- **Least privilege** service account permissions
- **Workload Identity Federation** for authentication
- **GCP Secret Manager** for secure key storage

## Migration from Current Setup

1. **Keep existing workflow** for GitHub-hosted runners as fallback
2. **Add secure workflow** for enhanced security when self-hosted runners available
3. **Test both workflows** to ensure functionality
4. **Gradually transition** to self-hosted as primary CI method

## Troubleshooting

### Common Issues
- **Permission denied**: Check service account has Secret Manager access
- **Workload identity failed**: Verify pool and provider configuration
- **SSH connection failed**: Check secret content and known_hosts
- **Runner not found**: Ensure self-hosted runner is online and available

### Debug Commands
```bash
# Check GCP authentication
gcloud auth list

# Test secret access
gcloud secrets versions access latest --secret="gadi_ssh_private_key"

# Verify workload identity
gcloud iam workload-identity-pools describe qxub-workload-identity-pool
```

## Files Reference

- **Secure Workflow**: `.github/workflows/test-remote-execution-secure.yml`
- **Setup Script**: `scripts/setup-secure-ci.sh`
- **Current Workflow**: `.github/workflows/test-remote-execution.yml` (fallback)

This setup provides enterprise-grade security for CI testing while maintaining the comprehensive v3.3.0 remote execution test coverage.
