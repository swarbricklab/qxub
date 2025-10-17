# Debugging and Verbosity: Essential Self-Help Tools

Before diving deeper into qxub's features, let's equip you with the essential debugging tools. These will help you troubleshoot issues independently and understand what qxub is doing behind the scenes.

## The `--dry` Flag: Your Best Friend

The most important debugging tool is `--dry` (or `--dry-run`). It shows you exactly what qxub would do **without actually submitting the job**.

### Basic Dry Run

```bash
qxub --dry -- hostname
```

**Expected output:**
```
🔍 DRY RUN - Would submit the following job:

📋 Job Configuration:
├── Name: qx-20241017-144052
├── Queue: normal
├── Project: a56
├── Resources: mem=4GB, ncpus=1, walltime=2:00:00
├── Output: /scratch/a56/jr9959/qxub/qx-20241017-144052_20241017-144052.out
├── Error: /scratch/a56/jr9959/qxub/qx-20241017-144052_20241017-144052.err
└── Log: /scratch/a56/jr9959/qxub/qx-20241017-144052_20241017-144052.log

📝 PBS Script Preview:
#!/bin/bash
#PBS -N qx-20241017-144052
#PBS -P a56
#PBS -q normal
#PBS -l mem=4GB,ncpus=1,walltime=2:00:00
#PBS -o /scratch/a56/jr9959/qxub/qx-20241017-144052_20241017-144052.out
#PBS -e /scratch/a56/jr9959/qxub/qx-20241017-144052_20241017-144052.err

cd "/g/data/a56/software/qsub_tools"
hostname

🚫 DRY RUN - No job submitted
```

This shows you:
- The exact PBS script that would be generated
- All resource allocations and file paths
- The working directory where the command would run
- **No actual job submission**

### Dry Run with Different Options

Let's see what happens with custom resources:

```bash
qxub --dry --resources mem=8GB,ncpus=2,walltime=30:00 -- python3 -c "print('Hello world')"
```

**Notice the differences:**
```
🔍 DRY RUN - Would submit the following job:

📋 Job Configuration:
├── Name: qx-20241017-144152
├── Queue: normal
├── Project: a56
├── Resources: mem=8GB, ncpus=2, walltime=0:30:00
...

📝 PBS Script Preview:
#!/bin/bash
#PBS -N qx-20241017-144152
#PBS -P a56
#PBS -q normal
#PBS -l mem=8GB,ncpus=2,walltime=0:30:00
...

cd "/g/data/a56/software/qsub_tools"
python3 -c "print('Hello world')"
```

The dry run shows exactly how your resource specifications were interpreted.

## Verbosity Levels: Getting More Information

qxub supports multiple verbosity levels to help you understand what's happening:

### Default Verbosity (Quiet)
```bash
qxub -- echo "Hello"
```
Shows only essential information during job execution.

### Verbose (`-v`)
```bash
qxub -v -- echo "Hello"
```

**Additional information shown:**
```
🔧 Loading configuration from: /g/data/a56/config/xdg/qxub/config.yaml
🎯 Platform detected: nci_gadi
📋 Applying defaults from configuration
🚀 Submitting job with qsub...
📋 Job submitted: 12345684.gadi-pbs (qx-20241017-144252)
🔄 Polling job status every 2 seconds...
⏳ Job queued, waiting for execution...
🎬 Job started, monitoring output streams...
✅ Job started, streaming output...

Hello

🎉 Job completed successfully (exit code: 0)
📊 Resource usage analysis...
📊 Walltime used: 00:00:03 / 02:00:00
💾 Memory used: 0.1GB / 4.0GB
📁 Output files created successfully
```

### Very Verbose (`-vv`)
```bash
qxub -vv -- echo "Hello"
```

**Even more detailed output:**
```
🔧 Loading user configuration from: ~/.config/qxub/config.yaml (not found)
🔧 Loading system configuration from: /g/data/a56/config/xdg/qxub/config.yaml
🔧 Configuration loaded successfully
🎯 Hostname: gadi-login-01
🎯 Platform search paths: ['/g/data/a56/config/xdg/qxub/platforms', ...]
🎯 Platform detected: nci_gadi (from /g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml)
📋 Template variables: {user: jr9959, project: a56, timestamp: 20241017-144352}
📋 Resolving output paths...
📋 Creating job script at: /tmp/qxub_12345685_script.pbs
🚀 PBS command: qsub /tmp/qxub_12345685_script.pbs
📋 Job submitted: 12345685.gadi-pbs (qx-20241017-144352)
🔄 Starting job monitoring thread...
🔄 Starting output streaming threads...
⏳ Job status: Q (queued)
⏳ Job status: R (running)
🎬 Output streaming started
✅ Job started, streaming output...

Hello

🎉 Job status: C (completed)
🎉 Job completed successfully (exit code: 0)
📊 Parsing PBS job log for resource usage...
📊 Walltime used: 00:00:03 / 02:00:00
💾 Memory used: 0.1GB / 4.0GB
🧹 Cleaning up temporary files...
📁 Output files created successfully
```

## Common Debugging Scenarios

### 1. Command Not Found Errors

If you get "command not found" errors:

```bash
qxub --dry -- some_unknown_command
```

The dry run will show you the working directory and environment. You might need:
- A different execution context (`--env`, `--mod`)
- To specify the full path to the command
- To check if the command exists in your PATH

### 2. Resource Issues

If jobs are rejected or queued indefinitely:

```bash
qxub --dry --resources mem=500GB,ncpus=48,walltime=72:00:00 -- echo "test"
```

The dry run might reveal:
- Resources exceed queue limits
- Walltime format issues
- Queue availability problems

### 3. Path and Directory Issues

To see exactly where your command will run:

```bash
qxub --dry -- pwd
```

This shows the working directory in the PBS script. If you need to run somewhere else:

```bash
qxub --dry -- bash -c 'cd /scratch/a56/jr9959 && pwd'
```

### 4. Environment Problems

To debug environment issues:

```bash
qxub -vv --env dvc3 -- which python
```

Very verbose output shows:
- Which conda environment is being activated
- Environment activation commands in the PBS script
- Path resolution

## Understanding Error Messages

qxub provides helpful error messages with suggestions:

### Invalid Resource Format
```bash
qxub --resources mem=4g -- echo "test"  # Wrong format
```

**Error output:**
```
❌ Error: Invalid memory format '4g'
💡 Hint: Use formats like '4GB', '4000MB', or '4096MB'
💡 Examples: --resources mem=8GB,ncpus=2
```

### Unknown Queue
```bash
qxub --queue nonexistent -- echo "test"
```

**Error output:**
```
❌ Error: Queue 'nonexistent' not found
💡 Available queues: normal, express, normalsl, hugemem, gpuvolta
💡 Or use '--queue auto' for automatic selection
```

### Multiple Execution Contexts
```bash
qxub --env dvc3 --mod python3/3.11.7 -- echo "test"
```

**Error output:**
```
❌ Error: Cannot specify multiple execution contexts
💡 Choose one: --env, --mod/--mods, or --sif
💡 Current contexts: conda environment 'dvc3', modules ['python3/3.11.7']
```

## Self-Debugging Workflow

When things go wrong, follow this process:

1. **Start with `--dry`**: See what would happen
   ```bash
   qxub --dry [your options] -- [your command]
   ```

2. **Add verbosity**: Get more details
   ```bash
   qxub -v --dry [your options] -- [your command]
   ```

3. **Check configuration**: Understand your defaults
   ```bash
   qxub config list --show-origin
   ```

4. **Verify the command locally**: Does it work outside qxub?
   ```bash
   [your command]  # Test locally first
   ```

5. **Test with minimal options**: Strip down to basics
   ```bash
   qxub -- [simple version of command]
   ```

## Getting Help

qxub has built-in help at multiple levels:

```bash
qxub --help                    # Main help
qxub config --help             # Config system help
qxub monitor --help            # Monitoring help
qxub --help | grep -A5 resource  # Search help for topics
```

## When to File Bug Reports

File a bug report when:
- qxub crashes with a traceback
- Error messages are unclear or unhelpful
- Expected behavior doesn't match documentation
- `--dry` shows something reasonable but execution fails mysteriously

Include in your bug report:
- The exact command that failed
- Output from `qxub -vv --dry [your command]`
- Your platform (usually `nci_gadi`)
- qxub version (`qxub --version`)

## Key Takeaways

1. **`--dry` is essential**: Always test complex commands first
2. **Verbosity helps**: Use `-v` or `-vv` when troubleshooting
3. **Error messages guide you**: Read them carefully for hints
4. **Start simple**: Test basic versions before adding complexity
5. **Configuration matters**: Check your defaults with `qxub config list`

## Next Steps

With these debugging tools in hand, you're ready to:
- **[Resources and Queues](02-resources-and-queues.md)** - Customize resource allocation
- **[Execution Contexts](04-execution-contexts.md)** - Use different software environments

Remember: when in doubt, `--dry` it out! This simple flag will save you time and help you understand qxub's behavior.

---

**💡 Pro Tips:**
- Always use `--dry` when testing new resource combinations
- Combine `--dry` with `-v` for maximum insight
- Save complex commands as aliases after testing with `--dry`
- The verbose output shows you exactly what PBS script would be generated
