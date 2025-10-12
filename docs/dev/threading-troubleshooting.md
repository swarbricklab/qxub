# Threading Troubleshooting Guide

This guide helps developers debug and fix issues in qxub's threading system.

## Quick Diagnostics

### Enable Debug Logging
```bash
export QXUB_LOG_LEVEL=DEBUG
qxub conda --env myenv script.py 2>&1 | grep -E "(Thread|Monitor|Tail|Spinner)"
```

### Check Thread Status
```python
# Add to any thread function for debugging
import threading
logging.debug(f"Thread {threading.current_thread().name} active, "
              f"total threads: {threading.active_count()}")
```

## Common Issues and Solutions

### 1. qxub Hangs on Exit

**Symptoms**: Process doesn't terminate after job completes
**Cause**: Thread not checking `should_shutdown()`

**Debug**:
```bash
# Check if threads are still running
ps aux | grep qxub
# Look for multiple python processes
```

**Fix**: Ensure all thread loops check shutdown condition:
```python
# In every thread loop:
while True:
    if coordinator and coordinator.should_shutdown():
        logging.debug("Thread shutdown requested")
        break
    # ... do work
```

**Prevention**: Always include timeout in thread joins:
```python
thread.join(timeout=5)
if thread.is_alive():
    logging.warning("Thread did not shutdown gracefully")
```

### 2. Spinner Doesn't Clear

**Symptoms**: Spinner continues after output starts
**Cause**: `output_started` not signaled or spinner not checking

**Debug**:
```python
# In tail function, add logging:
if not output_started and coordinator:
    logging.debug("Signaling output_started")
    coordinator.signal_output_started()
    output_started = True
```

**Fix**: Ensure tail threads signal on first output:
```python
# First line of output in tail function:
if not output_started and coordinator:
    coordinator.signal_output_started()
    # Clear spinner line immediately
    print(" " * 120, end="", flush=True)
    coordinator.signal_spinner_cleared()
    output_started = True
```

### 3. Wrong Exit Code Returned

**Symptoms**: qxub always returns 0 even when job fails
**Cause**: Exit status not properly propagated

**Debug**:
```python
# In monitor_qstat function:
logging.debug(f"Job exit status retrieved: {exit_status}")
if coordinator:
    coordinator.job_exit_status = exit_status
    logging.debug(f"Stored in coordinator: {coordinator.job_exit_status}")
```

**Fix**: Ensure exit status flows through:
1. `job_exit_status()` extracts from qstat
2. `wait_for_job_exit_status()` polls until available
3. `monitor_qstat()` stores in coordinator
4. `monitor_and_tail()` returns from coordinator
5. Executor (conda.py/etc) calls `sys.exit()`

### 4. Output Missing or Garbled

**Symptoms**: Job output doesn't appear in terminal
**Cause**: File permissions, paths, or thread exceptions

**Debug**:
```bash
# Check log files exist and are readable
ls -la /path/to/logs/
# Check file permissions
cat /path/to/out.log  # Should show job output
```

**Fix**: Ensure proper file handling:
```python
# In executor functions:
out.parent.mkdir(parents=True, exist_ok=True)
out.touch()  # Create empty file
# Verify file is writable
assert out.exists() and os.access(out, os.R_OK)
```

### 5. Race Conditions

**Symptoms**: Intermittent failures, inconsistent behavior
**Cause**: Threads accessing shared state without synchronization

**Debug**: Add timing logs:
```python
import time
start_time = time.time()
# ... in each critical section:
logging.debug(f"Event at {time.time() - start_time:.2f}s: {event_name}")
```

**Fix**: Use Events for coordination:
```python
# DON'T use shared variables
shared_flag = False  # BAD - race condition

# DO use Events
coordinator.my_event.set()  # GOOD - thread safe
```

## Performance Issues

### Excessive CPU Usage

**Cause**: Tight polling loops without sleep
**Fix**: Add sleep in polling loops:
```python
while self.spinning:
    # ... spinner work
    time.sleep(0.1)  # Don't consume 100% CPU
```

### Memory Leaks

**Cause**: Threads not properly cleaned up
**Debug**:
```python
import gc
import threading

# Before/after thread operations:
print(f"Threads: {threading.active_count()}")
print(f"Objects: {len(gc.get_objects())}")
```

**Fix**: Use daemon threads and proper cleanup:
```python
thread = threading.Thread(target=worker, daemon=True)
thread.start()
# Always join with timeout
thread.join(timeout=5)
```

## Testing Strategies

### Unit Testing Threads

```python
import threading
import time
import unittest

class TestThreading(unittest.TestCase):
    def test_thread_coordination(self):
        coordinator = OutputCoordinator()

        def worker():
            time.sleep(0.1)
            coordinator.signal_output_started()

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()

        # Test event signaling
        self.assertTrue(coordinator.wait_for_output_started(timeout=1))
        thread.join(timeout=1)
        self.assertFalse(thread.is_alive())
```

### Integration Testing

```python
def test_full_flow():
    # Use short job for testing
    result = subprocess.run([
        "qxub", "conda", "--env", "base", "echo", "test"
    ], capture_output=True, text=True, timeout=60)

    assert result.returncode == 0
    assert "test" in result.stdout
```

### Stress Testing

```python
# Test rapid job submission
for i in range(10):
    subprocess.Popen(["qxub", "conda", "--env", "base", "sleep", "1"])
    time.sleep(0.1)
```

## Debugging Tools

### Thread Inspection

```python
def debug_threads():
    """Print current thread status"""
    for thread in threading.enumerate():
        print(f"Thread: {thread.name}")
        print(f"  Alive: {thread.is_alive()}")
        print(f"  Daemon: {thread.daemon}")
        print(f"  Ident: {thread.ident}")
```

### Event Status

```python
def debug_coordinator(coordinator):
    """Print coordinator event states"""
    events = [
        'output_started', 'spinner_cleared', 'job_completed',
        'shutdown_requested', 'eof_detected'
    ]
    for event_name in events:
        event = getattr(coordinator, event_name)
        print(f"{event_name}: {event.is_set()}")
```

### Log Analysis

```bash
# Extract threading-related logs
grep -E "(Thread|Event|Coordinator)" qxub.log | sort

# Timeline analysis
grep "Thread" qxub.log | while read line; do
    echo "$(date): $line"
done
```

## Prevention Best Practices

### 1. Always Use Timeouts

```python
# Good
if event.wait(timeout=30):
    # Handle event
else:
    # Handle timeout

# Good
thread.join(timeout=5)
if thread.is_alive():
    # Handle hung thread
```

### 2. Graceful Error Handling

```python
try:
    # Thread work
    pass
except Exception as e:
    logging.error(f"Thread error: {e}")
    # Signal other threads about the error
    if coordinator:
        coordinator.signal_shutdown()
finally:
    # Always cleanup
    logging.debug("Thread exiting")
```

### 3. Resource Management

```python
# Use context managers for resources
with open(log_file, 'r') as f:
    # File automatically closed
    pass

# Clean up threads
try:
    # Start threads
    pass
finally:
    # Always shutdown threads
    coordinator.signal_shutdown()
    for thread in threads:
        thread.join(timeout=2)
```

### 4. Logging Standards

```python
# Include thread info in logs
import threading
logging.debug(f"[{threading.current_thread().name}] Event occurred")

# Use structured logging
logging.debug("Thread event", extra={
    'thread': threading.current_thread().name,
    'event': 'output_started',
    'timestamp': time.time()
})
```

## Emergency Procedures

### Kill Hung qxub Process

```bash
# Find qxub processes
ps aux | grep qxub

# Gentle termination
pkill -TERM qxub

# Force kill if needed
pkill -KILL qxub

# Check for orphaned jobs
qstat -u $USER
```

### Clean Up Orphaned Jobs

```bash
# List your jobs
qstat -u $USER

# Delete specific job
qdel JOB_ID

# Delete all your jobs (BE CAREFUL!)
qstat -u $USER | awk 'NR>2 {print $1}' | xargs qdel
```

### Reset Log Files

```bash
# If log files are corrupted
rm -rf $TMPDIR/qt/
mkdir -p $TMPDIR/qt/
```

## When to Escalate

Contact the qxub maintainers if you encounter:

1. **Consistent deadlocks** despite following best practices
2. **PBS-specific threading issues** on your cluster
3. **Performance degradation** with increasing job count
4. **Memory leaks** that persist after fixes
5. **Signal handling issues** on specific platforms

Include in your report:
- Debug logs with timestamps
- Thread dumps if available
- Cluster/PBS configuration details
- Reproducible test case
