# Remote Execution: Running Jobs from Your Laptop

qxub's remote execution capability allows you to submit and monitor HPC jobs directly from your laptop, making development workflows seamless and enabling sophisticated CI/CD pipelines. This section covers setup, usage patterns, and integration possibilities.

## Why Remote Execution?

Remote execution with qxub provides:

- **Local development**: Write and test on your laptop, execute on HPC
- **Seamless workflows**: No need to SSH and upload files manually
- **CI/CD integration**: Automated testing and deployment from GitHub Actions
- **Cross-platform support**: Work from any operating system
- **Real-time monitoring**: Watch HPC jobs from your local terminal

## Setting Up Remote Execution

### Prerequisites

On your laptop, you need:
- SSH access to the HPC system
- qxub installed locally (`pip install qxub`)
- SSH key-based authentication configured

### Basic Remote Configuration

Create a user configuration file on your laptop:

```bash
# On your laptop: ~/.config/qxub/config.yaml
mkdir -p ~/.config/qxub

cat > ~/.config/qxub/config.yaml << 'EOF'
# Remote execution configuration for NCI Gadi

remotes:
  gadi:
    url: ssh://gadi.nci.org.au
    conda_env: dvc3
    platform: nci_gadi
    working_dir: "/scratch/a56/$USER"
    force_tty: false
EOF
```

### SSH Key Setup

Ensure your SSH key is configured:

```bash
# Generate SSH key if you don't have one
ssh-keygen -t rsa -b 4096 -C "your.email@example.com"

# Copy public key to HPC system
ssh-copy-id your_username@gadi.nci.org.au

# Test connection
ssh your_username@gadi.nci.org.au "echo 'Connection successful'"
```

### Verify Remote Setup

```bash
# Test remote connection
qxub exec --remote gadi config list
```

This should show the remote system's configuration, confirming the connection works.

## Basic Remote Execution

### Simple Remote Commands

```bash
# Run a simple command on the remote HPC system
qxub exec --remote gadi --default -- hostname
```

**Expected output:**
```
🔗 Connecting to remote: gadi (gadi.nci.org.au)
🚀 Submitting job...
📋 Job submitted: 12345715.gadi-pbs (qx-20241017-161052)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

gadi-cpu-clx-1234

🎉 Job completed successfully (exit code: 0)
📊 Walltime used: 00:00:05 / 02:00:00
💾 Memory used: 0.1GB / 4.0GB
```

The job runs on the HPC system, but you see real-time output on your laptop!

### Remote Execution with Environments

```bash
# Use conda environment on remote system
qxub exec --remote gadi --env dvc3 -- python3 -c "
import pandas as pd
import sys
print(f'Python version: {sys.version.split()[0]}')
print(f'Pandas version: {pd.__version__}')
print('Running on HPC from laptop!')
"
```

### Remote Execution with Custom Resources

```bash
# Submit resource-intensive job from laptop
qxub exec --remote gadi --mem 16GB --ncpus 4 --walltime 2:00:00 --env sc -- python3 -c "
import scanpy as sc
import multiprocessing
print(f'Single-cell analysis environment ready')
print(f'CPUs available: {multiprocessing.cpu_count()}')
print('High-resource job running remotely!')
"
```

## File Synchronization Patterns

### Automatic File Upload

qxub can automatically sync your local files to the remote system:

```bash
# Sync local directory and run remote job
qxub exec --remote gadi --sync-up ~/my_project /scratch/a56/$USER/my_project --default -- python3 analysis.py
```

### Sync and Download Results

```bash
# Upload code, run analysis, download results
qxub exec --remote gadi \
  --sync-up ~/my_analysis /scratch/a56/$USER/analysis \
  --sync-down /scratch/a56/$USER/analysis/results ~/results \
  --env dvc3 --mem 8GB \
  -- python3 run_analysis.py
```

### Working Directory Patterns

```bash
# Set remote working directory
qxub exec --remote gadi --work-dir /scratch/a56/$USER/project1 --default -- python3 script.py

# Use project-specific remote directory
qxub exec --remote gadi --work-dir '{remote_work}/current_project' --env dvc3 -- python3 analysis.py
```

## Development Workflows

### Local Development → Remote Testing

```bash
# Develop locally, test remotely
cat > test_script.py << 'EOF'
import pandas as pd
import numpy as np

print("Testing data processing pipeline...")
data = pd.DataFrame(np.random.randn(1000, 5))
result = data.describe()
print("Pipeline test successful!")
print(result.head())
EOF

# Test on HPC with minimal resources
qxub exec --remote gadi --sync-up . /scratch/a56/$USER/dev_test test --default -- python3 test_script.py

# Scale up for production
qxub exec --remote gadi --sync-up . /scratch/a56/$USER/production \
  --env dvc3 --mem 16GB --ncpus 4 -- python3 production_script.py
```

### Iterative Analysis Workflow

```bash
# Develop analysis iteratively
for iteration in {1..3}; do
    echo "Running analysis iteration $iteration..."

    # Update analysis parameters
    sed -i "s/n_samples = .*/n_samples = $((iteration * 1000))/" analysis.py

    # Run remotely with results download
    qxub exec --remote gadi \
      --sync-up . /scratch/a56/$USER/iter_$iteration \
      --sync-down /scratch/a56/$USER/iter_$iteration/results ./results_iter_$iteration \
      --env dvc3 -- python3 analysis.py

    echo "Iteration $iteration completed, results downloaded"
done
```

## Parallel Remote Execution

### Remote Parallel Jobs

```bash
# Submit multiple remote jobs in parallel
REMOTE_JOBS=()
for param in 0.1 0.5 1.0 2.0; do
    JOB_ID=$(qxub exec --remote gadi --terse --name "remote-param-$param" \
      --env dvc3 --mem 8GB \
      -- python3 -c "
import time
import numpy as np

param = $param
print(f'Remote parameter analysis: {param}')
time.sleep(10)

result = np.random.randn(1000).mean() * param
print(f'Remote result for param {param}: {result:.4f}')
")
    REMOTE_JOBS+=($JOB_ID)
    echo "Submitted remote job with param $param: $JOB_ID"
done

# Monitor all remote jobs from laptop
echo "Monitoring ${#REMOTE_JOBS[@]} remote jobs..."
qxub exec --remote gadi monitor --wait-for-completion "${REMOTE_JOBS[@]}"
echo "All remote parallel jobs completed!"
```

### Cross-Platform Parallel Analysis

```bash
# Run different parts of analysis on different systems (if you have multiple remotes)
echo "Starting cross-platform analysis..."

# Heavy computation on HPC
HPC_JOB=$(qxub exec --remote gadi --terse --name "heavy-compute" \
  --mem 64GB --ncpus 16 --walltime 4:00:00 \
  -- python3 heavy_computation.py)

# Visualization on system with GPUs (example)
# VIZ_JOB=$(qxub exec --remote gpu_cluster --terse --name "visualization" \
#   --gpu 1 --mem 16GB \
#   -- python3 create_visualizations.py)

echo "HPC computation job: $HPC_JOB"
# echo "GPU visualization job: $VIZ_JOB"

# Monitor both from laptop
qxub exec --remote gadi monitor --wait-for-completion "$HPC_JOB"
echo "Cross-platform analysis completed!"
```

## CI/CD Integration

### GitHub Actions Integration

Create a GitHub Actions workflow that uses qxub for HPC testing:

```yaml
# .github/workflows/hpc-testing.yml
name: HPC Testing

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main ]

jobs:
  hpc-test:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v3

    - name: Setup Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.11'

    - name: Install qxub
      run: pip install qxub

    - name: Setup SSH key
      run: |
        mkdir -p ~/.ssh
        echo "${{ secrets.HPC_SSH_KEY }}" > ~/.ssh/id_rsa
        chmod 600 ~/.ssh/id_rsa
        ssh-keyscan -H gadi.nci.org.au >> ~/.ssh/known_hosts

    - name: Setup qxub remote config
      run: |
        mkdir -p ~/.config/qxub
        cat > ~/.config/qxub/config.yaml << EOF
        remotes:
          gadi:
            url: ssh://${{ secrets.HPC_USERNAME }}@gadi.nci.org.au
            platform: nci_gadi
            config: ~/.ssh/config
        default_remote: gadi
        EOF

    - name: Run HPC tests
      run: |
        # Upload test code and run
        qxub exec --remote gadi --sync-up . /scratch/${{ secrets.HPC_PROJECT }}/${{ secrets.HPC_USERNAME }}/ci_test_${{ github.run_id }} \
          test -- python3 -m pytest tests/ -v

    - name: Run integration test
      run: |
        # Run more comprehensive test with proper resources
        qxub exec --remote gadi --work-dir /scratch/${{ secrets.HPC_PROJECT }}/${{ secrets.HPC_USERNAME }}/ci_test_${{ github.run_id }} \
          --env dvc3 --mem 8GB --ncpus 2 --walltime 30:00 \
          -- python3 integration_test.py

    - name: Cleanup
      if: always()
      run: |
        # Clean up test directory
        qxub exec --remote gadi --default -- rm -rf /scratch/${{ secrets.HPC_PROJECT }}/${{ secrets.HPC_USERNAME }}/ci_test_${{ github.run_id }}
```

### Automated Deployment Pipeline

```yaml
# .github/workflows/deploy-analysis.yml
name: Deploy Analysis Pipeline

on:
  release:
    types: [published]

jobs:
  deploy:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v3

    - name: Setup qxub and deploy
      run: |
        pip install qxub

        # Setup SSH and config (as above)

        # Deploy analysis pipeline
        qxub exec --remote gadi --sync-up . /scratch/a56/production/analysis_v${{ github.event.release.tag_name }} \
          -- echo "Analysis pipeline v${{ github.event.release.tag_name }} deployed"

        # Run validation
        qxub exec --remote gadi --work-dir /scratch/a56/production/analysis_v${{ github.event.release.tag_name }} \
          --env dvc3 -- python3 validate_deployment.py

        # Notify deployment success
        echo "✅ Analysis pipeline deployed to HPC successfully"
```

## Advanced Remote Patterns

### Remote DVC Pipeline Execution

```bash
# Run entire DVC pipeline remotely from laptop
upload_and_run_pipeline() {
    local project_name=$1
    local remote_dir="/scratch/a56/$USER/$project_name"

    echo "Uploading project to HPC..."
    qxub exec --remote gadi --sync-up . "$remote_dir" -- echo "Project uploaded"

    echo "Running DVC pipeline remotely..."
    qxub exec --remote gadi --work-dir "$remote_dir" \
      --env dvc3 --mem 4GB --walltime 4:00:00 \
      -- dvc repro

    echo "Downloading results..."
    qxub exec --remote gadi --sync-down "$remote_dir/results" ./remote_results \
      --sync-down "$remote_dir/metrics" ./remote_metrics \
      -- echo "Results downloaded"

    echo "Remote DVC pipeline completed!"
}

# Usage
upload_and_run_pipeline "my_analysis_$(date +%Y%m%d)"
```

### Remote Jupyter Notebook Execution

```bash
# Convert notebook to script and run remotely
jupyter nbconvert --to script analysis.ipynb

# Run notebook script on HPC
qxub exec --remote gadi \
  --sync-up . /scratch/a56/$USER/notebook_run \
  --sync-down /scratch/a56/$USER/notebook_run/outputs ./notebook_outputs \
  --env jupyterlab --mem 16GB --ncpus 4 \
  -- python3 analysis.py
```

### Remote Resource Optimization

```bash
# Test resource requirements remotely
test_resources() {
    local script=$1
    local test_sizes=("4GB,1" "8GB,2" "16GB,4")

    for size_config in "${test_sizes[@]}"; do
        IFS=',' read -r mem ncpus <<< "$size_config"
        echo "Testing with $mem memory and $ncpus CPUs..."

        start_time=$(date +%s)
        qxub exec --remote gadi --mem "$mem" --ncpus "$ncpus" --walltime 30:00 \
          --env dvc3 -- python3 "$script"
        end_time=$(date +%s)

        duration=$((end_time - start_time))
        echo "Configuration $mem/$ncpus completed in ${duration}s"
    done
}

# Usage
test_resources "performance_test.py"
```

## Troubleshooting Remote Execution

### Connection Issues

```bash
# Test remote connection
qxub exec --remote gadi --dry -- echo "Connection test"

# Verbose connection debugging
qxub exec --remote gadi -vv --dry -- echo "Debug connection"

# Test SSH directly
ssh -v your_username@gadi.nci.org.au "echo 'Direct SSH test'"
```

### File Sync Issues

```bash
# Test file synchronization
qxub exec --remote gadi --sync-up test.txt /scratch/a56/$USER/sync_test/ \
  --dry -- echo "Sync test"

# Check remote file system
qxub exec --remote gadi --default -- ls -la /scratch/a56/$USER/
```

### Performance Optimization

```bash
# Use SSH connection multiplexing for faster subsequent connections
cat >> ~/.ssh/config << 'EOF'
Host gadi gadi.nci.org.au
    ControlMaster auto
    ControlPath ~/.ssh/control-%r@%h:%p
    ControlPersist 10m
EOF
```

## Security Considerations

### SSH Key Management

```bash
# Use specific SSH key for HPC
ssh-keygen -t ed25519 -f ~/.ssh/hpc_key -C "hpc-access"

# Add SSH config for specific key
cat >> ~/.ssh/config << 'EOF'
Host gadi
    HostName gadi.nci.org.au
    User your_username
    IdentityFile ~/.ssh/hpc_key
    ControlMaster auto
    ControlPersist 10m
EOF

# qxub will automatically use SSH config
# No changes needed to qxub config.yaml
```

### Secure File Handling

```bash
# Avoid syncing sensitive files
cat > .qxubignore << 'EOF'
*.key
*.pem
.env
secrets/
passwords.txt
EOF
```

## Key Takeaways

1. **Seamless development**: Code locally, execute on HPC without friction
2. **Real-time monitoring**: See HPC job output in your local terminal
3. **File synchronization**: Automatic upload/download of code and results
4. **CI/CD ready**: Perfect for automated testing and deployment
5. **Cross-platform**: Works from any laptop/desktop operating system

## Summary

Remote execution completes the qxub ecosystem, enabling:
- **Local development** with HPC execution power
- **Automated CI/CD** pipelines for research software
- **Cross-platform workflows** that work everywhere
- **Team collaboration** with consistent HPC access

You now have the complete qxub toolkit for modern HPC workflows!

---

**💡 Pro Tips:**
- Use SSH connection multiplexing for faster remote operations
- Create project-specific remote directories for better organization
- Use `--dry` to test remote commands before executing
- Set up `.qxubignore` to avoid uploading sensitive files
- Consider using remote execution for automated testing in CI/CD pipelines

**🎯 You've completed the full qxub tutorial!** You now understand everything from basic job submission to advanced parallel workflows and remote execution. Start applying these patterns to your own HPC work!
