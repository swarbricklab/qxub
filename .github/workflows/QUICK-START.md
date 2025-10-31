# Minimal Remote Execution Workflows - Quick Start

## What Was Created

Three new minimal workflows demonstrating qxub remote execution from GitHub Actions to NCI:

| File | Size | Purpose |
|------|------|---------|
| `ultra-minimal-ssh-test.yml` | 1.3K | SSH connectivity verification |
| `minimal-remote-test.yml` | 2.6K | Basic qxub remote execution |
| `hpci-style-remote.yml` | 2.3K | Production-ready pattern ⭐ |

Plus documentation:
- `README-minimal.md` - Detailed guide for all three workflows
- `PATTERN-COMPARISON.md` - Architecture comparison with hpci-scripts

## Quick Test

### 1. Test SSH Connectivity (30 seconds)
```bash
gh workflow run ultra-minimal-ssh-test.yml
```

Expected output:
```
✅ SSH connection successful
✅ qxub available on NCI
```

### 2. Test qxub Remote Execution (2 minutes)
```bash
gh workflow run minimal-remote-test.yml
```

Expected output:
```
✅ Dry-run execution validated
✅ Conda environment configured
✅ Configuration syntax correct
```

### 3. Submit Actual Job (5 minutes) ⭐ RECOMMENDED START
```bash
gh workflow run hpci-style-remote.yml
```

Expected output:
```
✅ Job submitted to NCI
✅ Job completed successfully
✅ Conda environment execution verified
```

## The Key Pattern

All workflows use this pattern learned from hpci-scripts:

```yaml
steps:
  # 1. Get SSH keys from GCP Secret Manager
  - uses: swarbricklab/hpci-scripts/.github/actions/setup-hpci@main
    with:
      submodule_token: ${{ secrets.SUBMODULE_TOKEN }}
      project_id: ${{ secrets.PROJECT_ID }}
      project_number: ${{ secrets.PROJECT_NUMBER }}
      registry_sa: ${{ secrets.REGISTRY_SA }}

  # 2. Install qxub on runner
  - run: pip install -e .

  # 3. Configure qxub for remote platform
  - run: |
      cat > ~/.config/qxub/config.yaml << 'EOF'
      platforms:
        gadi:
          remote:
            host: <user>@gadi.nci.org.au
            ssh_key: atlas_key  # Created by setup-hpci
            working_dir: /scratch/a56/<user>
      EOF

  # 4. Execute commands remotely
  - run: qxub exec --platform gadi -- <command>

  # 5. Clean up (automatic in setup-hpci, but good practice)
  - run: rm -f atlas_key*
```

## Key Differences from hpci-scripts

| Aspect | hpci-scripts | qxub Hybrid |
|--------|--------------|-------------|
| Remote tool | `hpci-exe` (custom binary) | `qxub exec` (standard) |
| Job definition | PBS script files | Command-line arguments |
| Flexibility | Fixed scripts | Any command |
| Configuration | Hardcoded in PBS | Hierarchical config |
| Features | Basic PBS | Full qxub features |

## Authentication (Same as hpci-scripts!)

Both patterns use identical secure authentication:
- ✅ GCP Secret Manager for SSH keys
- ✅ Workload Identity Federation
- ✅ Automatic key cleanup
- ✅ No long-lived secrets in GitHub

## Required GitHub Secrets

Set these in repository settings (same as hpci-scripts):

```bash
SUBMODULE_TOKEN   # GitHub token for submodules
PROJECT_ID        # GCP project ID
PROJECT_NUMBER    # GCP project number
REGISTRY_SA       # Service account email
USER_NAME         # NCI username
```

Optional:
```bash
ENABLE_ACTUAL_REMOTE_TEST=true  # Enable real job submission
```

## Migration from hpci-scripts

If you're using hpci-scripts workflows, migration is straightforward:

### Before (hpci-scripts):
```yaml
- run: |
    hpci-exe \
      --script hpci-scripts/checkout.pbs \
      -c REPO=myrepo,GIT_REF=main
```

### After (qxub):
```yaml
- run: |
    qxub exec --platform gadi --env dvc3 --queue normal -- \
      dvc checkout --relink --force
```

Benefits:
- ✅ No PBS script files to maintain
- ✅ More flexible (any command, any environment)
- ✅ Better error handling
- ✅ Built-in monitoring
- ✅ Same security model

## Recommended Next Steps

1. **Start with hpci-style-remote.yml** - Most production-ready
2. **Customize the command** in workflow dispatch input
3. **Test with your use case** (dvc, snakemake, etc.)
4. **Expand to full test suite** using test-remote-execution-secure.yml patterns

## Troubleshooting

### Workflow not showing up?
- Check you're on branch `feature/v3.3.0-remote-execution`
- Push to trigger `on: push` workflows
- Use `on: workflow_dispatch` for manual testing

### SSH fails?
```bash
# Run ultra-minimal first to isolate issue
gh workflow run ultra-minimal-ssh-test.yml
```

### qxub not found?
```bash
# Verify qxub environment exists on NCI:
ssh gadi.nci.org.au
conda activate qxub
which qxub
```

### Job submission fails?
```bash
# Check in workflow - qxub provides detailed output
# Look for: queue availability, project allocation, storage access
```

## Documentation

- **README-minimal.md** - Comprehensive guide for all workflows
- **PATTERN-COMPARISON.md** - Deep dive into architectural differences
- **docs/secure-ci-complete.md** - Full secure CI setup
- **docs/remote-execution.md** - qxub remote execution documentation

## Success Criteria

After running these workflows, you should see:

✅ SSH connection to NCI working
✅ qxub available in conda environment
✅ Job submission successful
✅ Job completion with correct output
✅ SSH keys automatically cleaned up

## Example Output

```
Run qxub exec --platform gadi --queue copyq --time "0:01:00" -vv -- hostname
[qxub] Platform: gadi
[qxub] Remote host: jr9959@gadi.nci.org.au
[qxub] Connecting via SSH...
[qxub] Creating job script...
[qxub] Submitting to PBS...
[qxub] Job ID: 12345678.gadi-pbs
[qxub] Monitoring job...
[qxub] Job completed successfully
[qxub] Output:
gadi-cpu-v001
✅ Remote execution completed
```

## Ready to Start?

```bash
# Test connectivity
gh workflow run ultra-minimal-ssh-test.yml

# Run minimal test
gh workflow run minimal-remote-test.yml

# Execute custom command
gh workflow run hpci-style-remote.yml -f command="hostname && date"
```

---

**Questions?** See `README-minimal.md` for detailed documentation or `PATTERN-COMPARISON.md` for architecture details.
