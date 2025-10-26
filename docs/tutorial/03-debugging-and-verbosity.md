# Debugging and Verbosity: Essential Self-Help Tools

Before diving deeper into qxub's features, let's equip you with the essential debugging tools. These will help you troubleshoot issues independently and understand what qxub is doing behind the scenes.

## The `--dry` Flag: Your Best Friend

The most important debugging tool is `--dry` (or `--dry-run`). It shows you exactly what qxub would do **without actually submitting the job**.

### Basic Dry Run

```bash
qxub exec --dry -- hostname
```

**Expected output:**
**Expected dry run output:**
```
� Job command constructed
📝 Command to execute: hostname
Dry run - job would be submitted (use -v to see full qsub command)
```

This shows you:
- The exact PBS script that would be generated
- All resource allocations and file paths
- The working directory where the command would run
- **No actual job submission**

### Dry Run with Different Options

Let's see what happens with custom resources:

```bash
qxub exec --dry --resources mem=8GB,ncpus=2,walltime=30:00 -- python parallel_script.py
```

**Expected output:**
```
� Job command constructed
� Command to execute: python parallel_script.py
Dry run - job would be submitted (use -v to see full qsub command)
```

### Verbose Dry Run

For more detail about the PBS submission, use `-v` with `--dry`:

```bash
qxub -v --dry --default -- echo "Hello"
```

**Expected output:**
```
🔧 Job command constructed
📝 Command to execute: echo Hello
🔧 Full qsub command: qsub -v cmd_b64="ZWNobyBIZWxsbw==",cwd=/g/data/a56/software/qsub_tools,out=/scratch/a56/jr9959/qt/20251018_154522/out,err=/scratch/a56/jr9959/qt/20251018_154522/err,quiet=false -N qt -q normal -P a56  -o qt.log /g/data/a56/software/qsub_tools/qxub/jobscripts/qdefault.pbs
```

This shows the exact `qsub` command that would be executed, including all parameters and file paths.

## Verbosity Levels: Getting More Information

qxub supports multiple verbosity levels to help you understand what's happening:

### Default Verbosity (Quiet)
```bash
qxub exec --default -- echo "Hello"
```
Shows only essential information during job execution.

### Verbose (`-v`)
```bash
qxub exec --default -v -- echo "Hello"
```

**Expected output:**
```
🔧 Job command constructed
✅ Job submitted successfully! Job ID: 152755696.gadi-pbs
Hello
✅ Command completed successfully
🎉 Job completed successfully
```

### Very Verbose (`-vv`)
```bash
qxub -vv --default -- echo "Hello"
```

**Expected output:**
```
INFO: Options: -N qt -q normal -P a56  -o qt.log
INFO: Submission command: qsub -v cmd_b64="ZWNobyBIZWxsbw==",cwd=/g/data/a56/software/qsub_tools,out=/scratch/a56/jr9959/qt/20251018_154333/out,err=/scratch/a56/jr9959/qt/20251018_154333/err,quiet=false -N qt -q normal -P a56  -o qt.log /g/data/a56/software/qsub_tools/qxub/jobscripts/qdefault.pbs
🔧 Job command constructed
✅ Job submitted successfully! Job ID: 152755741.gadi-pbs
Hello
✅ Command completed successfully
🎉 Job completed successfully
INFO: Job 152755741.gadi-pbs completed with status F
INFO: Waiting for PBS cleanup and getting exit status...
INFO: Job 152755741.gadi-pbs exit status: 0
```

Shows detailed internal operations including config loading, platform detection, job monitoring threads, and cleanup steps. Use for troubleshooting qxub itself.

## Common Debugging Scenarios

### 1. Command Not Found Errors

If you get "command not found" errors:

```bash
qxub exec --dry -- some_unknown_command
```

The dry run will show you the working directory and environment. You might need:
- A different execution context (`--env`, `--mod`)
- To specify the full path to the command
- To check if the command exists in your PATH

### 2. Resource Issues

If jobs are rejected or queued indefinitely:

```bash
qxub exec --dry --resources mem=500GB,ncpus=48,walltime=72:00:00 -- echo "test"
```

The dry run might reveal:
- Resources exceed queue limits
- Walltime format issues
- Queue availability problems

### 3. Path and Directory Issues

To see exactly where your command will run:

```bash
qxub exec --dry -- pwd
```

This shows the working directory in the PBS script. If you need to run somewhere else:

```bash
qxub exec --dry -- ./run_in_scratch.sh
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
qxub exec --resources mem=4g -- echo "test"  # Wrong format
```

**Error output:**
```
❌ Error: Invalid memory format '4g'
💡 Hint: Use formats like '4GB', '4000MB', or '4096MB'
💡 Examples: --resources mem=8GB,ncpus=2
```

### Unknown Queue
```bash
qxub exec --queue nonexistent -- echo "test"
```

**Error output:**
```
❌ Error: Queue 'nonexistent' not found
💡 Available queues: normal, express, normalsl, hugemem, gpuvolta
💡 Or use '--queue auto' for automatic selection
```

### Multiple Execution Contexts
```bash
qxub exec --env dvc3 --mod python3/3.11.7 -- echo "test"
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
   qxub exec --dry [your options] -- [your command]
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
   qxub exec --default -- [simple version of command]
   ```

## Getting Help

qxub has built-in help at multiple levels:

```bash
qxub exec --help                    # Main help
qxub config --help             # Config system help
qxub monitor --help            # Monitoring help
qxub exec --help | grep -A5 resource  # Search help for topics
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
- qxub version (`qxub exec --version`)

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
