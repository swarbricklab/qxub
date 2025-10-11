# CI Remote Execution Testing Setup

This document explains how to configure GitHub Actions secrets for testing qxub remote execution against NCI Gadi.

## Prerequisites on NCI Gadi

Before setting up CI, run the prerequisites check script on NCI Gadi:

```bash
# SSH to NCI Gadi
ssh gadi.nci.org.au

# Clone the repository or download the script
cd /scratch/$PROJECT/$USER
git clone https://github.com/swarbricklab/qsub_tools.git
cd qsub_tools

# Run the prerequisites check
./scripts/check-nci-prerequisites.sh
```

This script will verify:
- Conda installation and available environments
- qxub availability in different environments
- Platform file locations
- Project directory structure
- SSH connectivity

## Required GitHub Secrets

Configure these secrets in your GitHub repository settings:

### Essential Secrets

1. **`NCI_SSH_PRIVATE_KEY`** (required)
   - Your SSH private key for connecting to NCI Gadi
   - Generate with: `ssh-keygen -t ed25519 -C "ci@github"`
   - Add the public key to `~/.ssh/authorized_keys` on NCI Gadi
   - Copy the private key content (including `-----BEGIN` and `-----END` lines)

2. **`NCI_SSH_HOST_KEY`** (required)
   - NCI Gadi's SSH host key for verification
   - Get with: `ssh-keyscan -t rsa,ed25519 gadi.nci.org.au`
   - Copy the output line(s)

3. **`NCI_USERNAME`** (required)
   - Your NCI username (e.g., `jr9959`)

### Optional Secrets

4. **`NCI_PROJECT`** (optional, default: `a56`)
   - Your NCI project code for scratch directory access

5. **`NCI_QXUB_ENV`** (optional, default: `base`)
   - Conda environment containing qxub on NCI Gadi
   - Use output from prerequisites check script

6. **`NCI_PLATFORM_FILE`** (optional)
   - Path to platform file on NCI Gadi
   - Default: `/g/data/a56/software/qsub_tools/docs/platforms/nci_gadi.yaml`
   - Use a path found by prerequisites check script

7. **`ENABLE_ACTUAL_REMOTE_TEST`** (optional)
   - Set to `'true'` to enable actual job submission tests
   - Leave unset for dry-run only testing

## Setting Up SSH Keys

### 1. Generate SSH Key Pair
```bash
# On your local machine or CI environment
ssh-keygen -t ed25519 -f ~/.ssh/nci_ci_key -C "ci-test@qxub"
```

### 2. Add Public Key to NCI Gadi
```bash
# Copy public key to NCI Gadi
ssh-copy-id -i ~/.ssh/nci_ci_key.pub your-username@gadi.nci.org.au

# Or manually:
cat ~/.ssh/nci_ci_key.pub | ssh your-username@gadi.nci.org.au 'cat >> ~/.ssh/authorized_keys'
```

### 3. Get Host Key
```bash
ssh-keyscan -t rsa,ed25519 gadi.nci.org.au
```

### 4. Add to GitHub Secrets
- Go to your repository → Settings → Secrets and variables → Actions
- Add new repository secrets with the values above

## Testing the Setup

Once secrets are configured, push to the `feature/v2.2-development` branch to trigger the CI workflow:

```bash
git add .
git commit -m "Add CI remote execution testing"
git push origin feature/v2.2-development
```

Monitor the GitHub Actions logs to see:
- SSH connectivity tests
- Remote qxub availability
- Platform file access
- Dry-run remote execution
- (Optional) Actual job submission

## Workflow Behavior

### With SSH Credentials
The workflow will:
1. Set up SSH configuration for NCI Gadi
2. Test SSH connectivity
3. Verify remote qxub installation
4. Check platform file access
5. Test dry-run remote execution
6. (If enabled) Submit actual test job

### Without SSH Credentials
The workflow will:
1. Skip all remote tests
2. Run local qxub functionality tests
3. Display instructions for enabling remote testing

## Security Notes

- SSH private keys are stored as GitHub secrets (encrypted)
- Keys are only used during CI runs
- Consider using dedicated CI keys (not your personal SSH keys)
- Regularly rotate CI SSH keys
- Monitor CI logs for any security issues

## Troubleshooting

### SSH Connection Issues
- Verify SSH key format (include BEGIN/END lines)
- Check NCI Gadi SSH service status
- Verify username and host key accuracy

### qxub Not Found
- Check conda environment name in `NCI_QXUB_ENV`
- Verify qxub installation in that environment
- Run prerequisites check script on NCI Gadi

### Platform File Issues
- Verify path exists and is readable
- Check file permissions
- Try alternative platform file locations

### Project Directory Issues
- Verify PROJECT environment variable
- Check scratch directory permissions
- Ensure sufficient disk space
