# Testing Remote Execution from Your Laptop

This guide shows how to test qxub remote execution from your laptop using `tests/test_config.yaml`.

## Prerequisites

### 1. SSH Configuration
Configure SSH access to Gadi in your `~/.ssh/config`:

```
Host gadi
  HostName gadi.nci.org.au
  User <your-username>
  IdentityFile ~/.ssh/id_rsa
```

Test SSH access:
```bash
ssh gadi "echo 'SSH works!'"
```

### 2. Remote qxub Installation
Ensure qxub is installed on Gadi in a conda environment:

```bash
ssh gadi "conda activate qxub && qxub --version"
```

If not installed, see the main README for installation instructions.

### 3. Update test_config.yaml
Edit `tests/test_config.yaml` and change these values:
- `defaults.project`: Your NCI project code (currently `a56`)
- `defaults.storage`: Your required storage volumes
- `platforms.gadi_remote.remote.working_dir`: Update project code in path

## Testing from Your Laptop

### 1. Clone the repository
```bash
git clone https://github.com/swarbricklab/qxub.git
cd qxub
git checkout feature/v3.3.0-remote-execution
```

### 2. Install qxub locally (optional but recommended)
```bash
pip install -e .
```

Or use `python -m qxub.cli` instead of `qxub` in commands below.

### 3. Test remote execution (dry run)
```bash
# Dry run - shows what would be executed without submitting
qxub exec --config tests/test_config.yaml \
  --platform gadi_remote \
  --dry \
  -- echo "Hello from laptop"
```

**Expected output**: SSH command with wrapped qxub command

### 4. Test actual remote submission
```bash
# Submit a real job via SSH
qxub exec --config tests/test_config.yaml \
  --platform gadi_remote \
  --env qxub \
  -- qxub --version
```

**Expected output**:
- SSH connection to Gadi
- Remote qxub execution
- PBS job submission on Gadi
- Job ID and output

### 5. Test with your own script
```bash
# Submit your Python script to Gadi from laptop
qxub exec --config tests/test_config.yaml \
  --platform gadi_remote \
  --env pytorch \
  --queue normal \
  --walltime 2:00:00 \
  --mem 8GB \
  --cpus 4 \
  -- python train.py --epochs 10
```

## How It Works

```
┌─────────────┐                              ┌──────────────┐
│   Laptop    │         SSH Connection       │  Gadi (NCI)  │
│             │──────────────────────────────>│              │
│ qxub exec   │                              │  qxub exec   │
│ --config    │   Serialized qxub command    │  --env ...   │
│ --platform  │   + environment setup        │  -- cmd      │
│ gadi_remote │                              │              │
│ -- cmd      │                              │  PBS qsub    │
└─────────────┘                              └──────────────┘
```

1. Your laptop runs `qxub exec --config tests/test_config.yaml --platform gadi_remote`
2. qxub detects `remote:` section in platform config → Remote execution mode
3. qxub serializes the command back to CLI format
4. qxub wraps with environment setup: `cd working_dir && conda_init && export QXUB_PLATFORM=... && qxub exec ...`
5. qxub executes via SSH: `ssh gadi "wrapped command"`
6. Remote qxub on Gadi parses the command and submits to PBS
7. Output streams back to your laptop terminal

## Platform Options

### gadi_remote
Remote execution via SSH (use from laptop/CI)
```bash
qxub exec --config tests/test_config.yaml --platform gadi_remote -- cmd
```

### gadi_local
Direct PBS submission (use when logged into Gadi)
```bash
qxub exec --config tests/test_config.yaml --platform gadi_local -- cmd
```

### gadi_remote_file
Remote execution using local platform definition file (offline mode)
```bash
qxub exec --config tests/test_config.yaml --platform gadi_remote_file -- cmd
```

## Troubleshooting

### SSH Connection Failed
```
Error: ssh: Could not resolve hostname gadi
```
**Solution**: Add Gadi to your `~/.ssh/config` (see Prerequisites)

### Remote qxub Not Found
```
Error: qxub: command not found
```
**Solution**: Update `conda_init` in `test_config.yaml` to activate the correct conda environment

### Permission Denied
```
Error: Permission denied (publickey)
```
**Solution**: Check your SSH key is added to NCI: https://my.nci.org.au/mancini/

### Platform Definition Not Found
```
Error: Failed to load platform definition
```
**Solution**: Switch to `gadi_remote_file` platform which uses local file path, or check internet connection

## CI/CD Integration

### GitHub Actions
```yaml
name: Test on Gadi
on: [push]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Setup SSH
        run: |
          mkdir -p ~/.ssh
          echo "${{ secrets.GADI_SSH_KEY }}" > ~/.ssh/id_rsa
          chmod 600 ~/.ssh/id_rsa
          echo "Host gadi" >> ~/.ssh/config
          echo "  HostName gadi.nci.org.au" >> ~/.ssh/config
          echo "  User ${{ secrets.GADI_USERNAME }}" >> ~/.ssh/config
      - name: Run tests on Gadi
        run: |
          pip install -e .
          qxub exec --config tests/test_config.yaml \
            --platform gadi_remote \
            -- pytest tests/
```

### GitLab CI
```yaml
test_on_gadi:
  script:
    - eval $(ssh-agent -s)
    - echo "$GADI_SSH_KEY" | ssh-add -
    - qxub exec --config tests/test_config.yaml --platform gadi_remote -- pytest
```

## Next Steps

Once remote execution works, you can:
1. Copy `tests/test_config.yaml` to your project as `.qxub/config.yaml`
2. Customize platforms for your workflow
3. Remove `--config` option (qxub will auto-detect project config)
4. See `docs/remote-execution.md` for advanced usage
