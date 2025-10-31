# Remote Execution Pattern Comparison

## Overview

This document compares three approaches to remote job execution from GitHub Actions to NCI:

1. **hpci-scripts pattern** - Custom executable with PBS scripts
2. **qxub remote execution** - Built-in remote execution feature
3. **Hybrid approach** - qxub with hpci-scripts authentication

## Pattern 1: hpci-scripts (Original)

### Architecture
```
GitHub Runner → hpci-exe → SSH → NCI → qsub PBS_SCRIPT
```

### Example from xenium repo
```yaml
steps:
  - uses: swarbricklab/hpci-scripts/.github/actions/setup-hpci@main
  - name: Schedule jobs
    run: |
      hpci-exe \
        --user ${{secrets.USER_NAME}} \
        --host ${{secrets.HOST_NAME}} \
        --publicKey atlas_key.pub \
        --privateKey atlas_key \
        --script hpci-scripts/fetch.pbs \
        --logFile fetch.log \
        -c INSTANCE_DIR=/g/data/a56/dvc/instances,REPO=xenium
```

### Characteristics
- ✅ Centralized PBS scripts in hpci-scripts repo
- ✅ Simple command-line interface
- ✅ Battle-tested in production
- ❌ Requires custom hpci-exe binary
- ❌ PBS scripts must be pre-written
- ❌ Limited to PBS script functionality

## Pattern 2: qxub Remote Execution (Pure)

### Architecture
```
GitHub Runner → qxub → SSH → NCI → qsub (generated script)
```

### Example
```yaml
steps:
  - name: Install qxub
    run: pip install -e .

  - name: Configure remote
    run: |
      cat > ~/.config/qxub/config.yaml << EOF
      platforms:
        gadi:
          remote:
            host: user@gadi.nci.org.au
            ssh_key: ~/.ssh/id_rsa
      EOF

  - name: Execute remotely
    run: qxub exec --platform gadi -- python script.py
```

### Characteristics
- ✅ No custom binaries needed
- ✅ Dynamic job script generation
- ✅ Flexible command-line interface
- ✅ Conda/module/container support built-in
- ❌ Requires qxub installation on runner
- ❌ Configuration must be created
- ❌ SSH key management needed

## Pattern 3: Hybrid (Recommended)

### Architecture
```
GitHub Runner → qxub → SSH (via hpci-scripts keys) → NCI → qsub
```

### Example (from minimal-remote-test.yml)
```yaml
steps:
  - name: Install qxub
    run: pip install -e .

  - name: Set up HPCI authentication
    uses: swarbricklab/hpci-scripts/.github/actions/setup-hpci@main

  - name: Execute via qxub
    run: qxub exec --platform gadi --ssh-key atlas_key -- hostname
```

### Characteristics
- ✅ Leverages hpci-scripts authentication infrastructure
- ✅ Uses qxub's flexible execution model
- ✅ No custom binaries needed
- ✅ Minimal configuration
- ✅ Best of both worlds
- ➖ Requires both repos (but hpci-scripts is standard)

## Feature Comparison

| Feature | hpci-exe | qxub Pure | qxub Hybrid |
|---------|----------|-----------|-------------|
| SSH Key Management | ✅ GCP Secrets | ⚠️ Manual | ✅ GCP Secrets |
| Custom Commands | ❌ PBS only | ✅ Any command | ✅ Any command |
| Conda Environments | ⚠️ In PBS | ✅ Built-in | ✅ Built-in |
| Module Loading | ⚠️ In PBS | ✅ Built-in | ✅ Built-in |
| Container Support | ❌ | ✅ Built-in | ✅ Built-in |
| Queue Selection | ⚠️ In PBS | ✅ Auto/Manual | ✅ Auto/Manual |
| Resource Specs | ⚠️ In PBS | ✅ CLI flags | ✅ CLI flags |
| Job Monitoring | ✅ Via hpci-exe | ✅ Via qxub | ✅ Via qxub |
| Configuration | ❌ Hardcoded | ✅ Hierarchical | ✅ Hierarchical |
| Template Variables | ❌ | ✅ Full support | ✅ Full support |

## Use Case Recommendations

### Use hpci-scripts pattern when:
- You have existing PBS scripts that work
- You need centralized script management
- Multiple repos use the same workflows
- You don't need dynamic commands

### Use qxub pure pattern when:
- You manage your own SSH infrastructure
- You need maximum flexibility
- You don't use hpci-scripts ecosystem
- You want complete control

### Use qxub hybrid pattern when:
- You want modern features + proven auth
- You're already using hpci-scripts for other repos
- You need flexible execution with secure auth
- You want minimal setup complexity ⭐ **RECOMMENDED**

## Migration Path

### From hpci-scripts to qxub hybrid:

1. **Keep existing workflow** as reference
2. **Add qxub installation** step
3. **Replace hpci-exe calls** with qxub exec
4. **Remove PBS script dependency** (optional)
5. **Test side-by-side** before removing old workflow

Example transformation:

**Before (hpci-scripts):**
```yaml
- run: |
    hpci-exe \
      --script hpci-scripts/fetch.pbs \
      -c REPO=myrepo,GIT_REF=main
```

**After (qxub hybrid):**
```yaml
- run: |
    qxub exec --platform gadi --env dvc3 -- \
      git -C /g/data/a56/myrepo fetch
```

### From qxub pure to qxub hybrid:

1. **Add setup-hpci action**
2. **Point ssh_key to atlas_key**
3. **Remove manual SSH setup steps**
4. **Test authentication**

**Before:**
```yaml
platforms:
  gadi:
    remote:
      ssh_key: ~/.ssh/gadi_key
```

**After:**
```yaml
platforms:
  gadi:
    remote:
      ssh_key: atlas_key  # Loaded by setup-hpci
```

## Implementation Examples

See these files for complete working examples:

1. **Ultra-minimal SSH test**: `.github/workflows/ultra-minimal-ssh-test.yml`
   - Just tests connectivity
   - Good for debugging

2. **Minimal remote test**: `.github/workflows/minimal-remote-test.yml`
   - Basic qxub remote execution
   - Dry-run and real submission

3. **HPCI-style remote**: `.github/workflows/hpci-style-remote.yml`
   - Production-ready pattern
   - Custom command support
   - Recommended starting point

4. **Full test suite**: `.github/workflows/test-remote-execution-secure.yml`
   - Comprehensive testing
   - All qxub features
   - CI validation

## Security Comparison

| Aspect | hpci-exe | qxub Hybrid |
|--------|----------|-------------|
| Key Storage | ✅ GCP Secret Manager | ✅ GCP Secret Manager |
| Key Lifecycle | ✅ Auto cleanup | ✅ Auto cleanup |
| Authentication | ✅ Workload Identity | ✅ Workload Identity |
| Least Privilege | ✅ Service Account | ✅ Service Account |
| Audit Trail | ✅ GCP Logging | ✅ GCP Logging |

Both patterns share the same secure authentication infrastructure!

## Conclusion

The **qxub hybrid pattern** provides the best combination of:
- Security (hpci-scripts authentication)
- Flexibility (qxub features)
- Simplicity (minimal configuration)
- Familiarity (for hpci-scripts users)

Start with `hpci-style-remote.yml` and customize for your needs.
