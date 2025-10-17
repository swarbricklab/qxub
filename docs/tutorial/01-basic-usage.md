# Basic Usage: Your First qxub Commands

Welcome to your first hands-on experience with `qxub`! This section will demonstrate the fundamental difference between traditional PBS and qxub: **real-time output and simplified submission**.

## The qxub Philosophy

Traditional PBS workflow:
1. Write a PBS script
2. Submit with `qsub script.pbs`
3. Wait and hope
4. Check output files later

qxub workflow:
1. Run `qxub --default -- your_command`
2. Watch output in real-time
3. Get results immediately

## Your First qxub Command

Let's start with something simple - a basic system information command:

```bash
qxub --default -- hostname
```

**Expected output:**
```
🚀 Submitting job...
📋 Job submitted: 12345678.gadi-pbs (qx-20241017-143052)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

gadi-login-01

🎉 Job completed successfully (exit code: 0)
📊 Walltime used: 00:00:05 / 02:00:00
💾 Memory used: 0.1GB / 4.0GB
📁 Outputs: /scratch/a56/jr9959/qxub/qx-20241017-143052_20241017-143052.{out,err,log}
```

**What just happened?**
- qxub automatically created a PBS job with sensible defaults
- You saw the job ID and tracking information immediately
- The command output (`gadi-login-01`) streamed in real-time
- Resource usage was reported at the end
- All output files were created with meaningful names

## Real-time Output Magic

Let's try something that takes a bit longer to really see the streaming in action:

```bash
qxub --default -- bash -c 'for i in {1..10}; do echo "Processing step $i"; sleep 2; done'
```

**What you'll see:**
```
🚀 Submitting job...
📋 Job submitted: 12345679.gadi-pbs (qx-20241017-143152)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

Processing step 1
Processing step 2
Processing step 3
...continuing in real-time...
Processing step 10

🎉 Job completed successfully (exit code: 0)
```

Notice how each line appears as it's generated - no waiting for the job to finish!

## Basic Python Example

Let's run a simple Python script:

```bash
qxub -- python3 -c "
import sys
import time
print(f'Python version: {sys.version}')
print('Computing something important...')
time.sleep(3)
print('Result: 42')
"
```

**Expected output:**
```
🚀 Submitting job...
📋 Job submitted: 12345680.gadi-pbs (qx-20241017-143252)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

Python version: 3.11.7 (main, Dec 15 2023, 12:09:04) [GCC 11.2.0] on linux
Computing something important...
Result: 42

🎉 Job completed successfully (exit code: 0)
```

## Understanding the Default Configuration

The commands above used the system defaults. Let's see what those are:

```bash
qxub config get defaults
```

**You should see:**
```
📋 Configuration: defaults
├── project: a56
├── queue: normal
├── resources: ['mem=4GB', 'ncpus=1', 'walltime=2:00:00']
├── name: qx-{timestamp}
├── out: {log_dir}/{name}_{timestamp}.out
├── err: {log_dir}/{name}_{timestamp}.err
└── joblog: {log_dir}/{name}_{timestamp}.log
```

These defaults mean every job gets:
- **4GB memory** and **1 CPU** - perfect for simple tasks
- **2 hours walltime** - plenty for most quick jobs
- **Meaningful names** with timestamps for uniqueness
- **Organized output** in your scratch directory under `qxub/`

## Working with Output and Errors

Let's demonstrate how qxub handles both stdout and stderr:

```bash
qxub -- bash -c 'echo "This goes to stdout"; echo "This goes to stderr" >&2; echo "Back to stdout"'
```

**You'll see both streams in real-time:**
```
🚀 Submitting job...
📋 Job submitted: 12345681.gadi-pbs (qx-20241017-143352)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

This goes to stdout
This goes to stderr
Back to stdout

🎉 Job completed successfully (exit code: 0)
```

The output files are still created separately:
- `.out` file contains stdout
- `.err` file contains stderr
- `.log` file contains the PBS job log

## Handling Job Failures

What happens when a command fails? Let's see:

```bash
qxub -- bash -c 'echo "Starting..."; sleep 2; exit 1'
```

**Expected output:**
```
🚀 Submitting job...
📋 Job submitted: 12345682.gadi-pbs (qx-20241017-143452)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

Starting...

❌ Job failed with exit code: 1
📊 Walltime used: 00:00:03 / 02:00:00
💾 Memory used: 0.1GB / 4.0GB
📁 Outputs: /scratch/a56/jr9959/qxub/qx-20241017-143452_20241017-143452.{out,err,log}
```

qxub clearly reports failures and still provides resource usage information.

## Interrupting Jobs

You can interrupt running jobs with `Ctrl+C`. Let's try:

```bash
qxub -- sleep 30
# Press Ctrl+C after a few seconds
```

**What happens:**
```
🚀 Submitting job...
📋 Job submitted: 12345683.gadi-pbs (qx-20241017-143552)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

^C
⚠️  Interrupted by user
🗑️  Cleaning up job 12345683.gadi-pbs...
✅ Job cancelled successfully
```

qxub automatically cleans up interrupted jobs - no orphaned processes!

## Key Takeaways

1. **Simple syntax**: `qxub -- your_command` - that's it!
2. **Real-time feedback**: See output as it happens, not after
3. **Automatic cleanup**: Sensible defaults and automatic resource management
4. **Error handling**: Clear reporting of both success and failure
5. **Interruption safety**: Ctrl+C cleanly cancels jobs

## Next Steps

Now that you've experienced the qxub basics, you're ready to learn about:
- **[Debugging and Verbosity](03-debugging-and-verbosity.md)** - Essential tools for troubleshooting
- **[Resources and Queues](02-resources-and-queues.md)** - Customizing resources and queue selection

The real-time output feature alone makes qxub a game-changer for interactive HPC work. In the next sections, we'll explore how to customize resources and use different execution environments.

---

**💡 Pro Tips:**
- Use `qxub --help` to see all available options
- The `--` separator is crucial - everything after it is your command
- Output files are always created even when you see real-time output
- Job names include timestamps to avoid conflicts
