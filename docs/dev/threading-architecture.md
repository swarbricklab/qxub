# qxub Threading Architecture

qxub uses **single-threaded monitoring** for reliable job tracking with clean output streaming.

## Current Architecture (v3+)

**Single Thread Monitoring**: Uses `monitor_job_single_thread()` which combines:
- Job status polling (via status files created by job scripts)
- Output file streaming (STDOUT/STDERR)
- Progress indication (spinner during job startup)
- Resource collection triggering (background thread after completion)

## Background Processes

1. **Main Thread** - Handles job submission and single-threaded monitoring
2. **Resource Collection Thread** - Daemon thread for post-job resource analysis (60+ second delay)

## Key Implementation Details

**Status File Monitoring**: Job scripts create status files that monitoring watches:
- `job_started_{job_id}` - Written when job begins execution
- `final_exit_code_{job_id}` - Written when job completes with exit status

**Output Streaming**: Tail-follows log files with position tracking for real-time display

**Clean Shutdown**: Handles Ctrl+C with `qdel` cleanup and proper exit codes

## Execution Modes

qxub supports three execution modes that control UI output during monitoring:

| Mode | Job ID Display | Progress Messages | Monitoring | Resource Collection |
|------|---------------|------------------|------------|-------------------|
| **Normal** | ✅ With success message | ✅ Full progress UI | ✅ With output | ✅ Yes |
| **Terse** | ✅ Job ID only | ❌ Silent | ✅ Silent | ✅ Yes |
| **Quiet** | ❌ None | ❌ Silent | ✅ Silent | ✅ Yes |

**Normal Mode**: Interactive use with full progress indication and job output streaming.

**Terse Mode**: Pipeline-friendly - emits job ID immediately then does silent monitoring. Perfect for scripts that need the job ID but want monitoring/resource collection.

**Quiet Mode**: Completely silent operation but still performs full monitoring and resource collection. Useful for automated systems that don't need any output.

## OutputCoordinator

Central hub managing thread lifecycle and events:
- `output_started` - Triggered when output begins, stops spinner
- `job_running` - Set when job transitions from queued to running
- `job_finished` - Set when job completes (success or failure)
- `shutdown_requested` - Graceful shutdown signal for all threads

## Signal Handling

Ctrl+C cleanup:
```python
def _signal_handler(signum, frame):
    global _CURRENT_JOB_ID
    if _CURRENT_JOB_ID:
        subprocess.run(['qdel', _CURRENT_JOB_ID])
    sys.exit(1)
```

## Key Implementation Details

**Thread Safety**: Uses `threading.Event` objects for coordination
**Clean Display**: Carriage returns (`\r`) clear spinner before output
**Exit Code Propagation**: Monitor thread writes job exit status to coordinator
**Resource Cleanup**: Context managers ensure proper thread termination

## Debug Commands

```bash
# Enable threading debug logging
export QXUB_LOG_LEVEL=DEBUG
qxub --env myenv -- python script.py

# Check for hanging processes
ps aux | grep qxub
```

See `qxub/scheduler.py` for implementation details.

The `OutputCoordinator` is the central nervous system that coordinates all threads using Python `threading.Event` objects:

```python
class OutputCoordinator:
    def __init__(self):
        self.output_started = threading.Event()      # Output begins streaming
        self.spinner_cleared = threading.Event()     # Spinner has been cleared
        self.job_completed = threading.Event()       # Job finished (from qstat)
        self.shutdown_requested = threading.Event()  # Ctrl-C pressed
        self.eof_detected = threading.Event()        # Log files reached EOF
        self.submission_complete = threading.Event() # Job submission phase complete
        self.job_running = threading.Event()         # Job status became "R"
        self.job_finished = threading.Event()        # Job status became "F" or "H"
        self.job_exit_status = None                   # Final job exit code
```

**New Event-Driven Events**:
- `job_running` - Signaled when job status changes to "R" (Running)
- `job_finished` - Signaled when job status changes to "F" (Finished) or "H" (Held)
- These events trigger immediate spinner shutdown instead of arbitrary timeouts

**Events are boolean flags that threads can wait for or signal**:
- `wait()` - Block until the event is set
- `set()` - Signal that the event has occurred
- `is_set()` - Check if the event has occurred (non-blocking)

### Coordinated Job Submission

The `start_job_monitoring()` function in `scheduler.py` provides coordinated thread startup:

```python
def start_job_monitoring(job_id, out_file, err_file, quiet=False, success_msg=None):
    """
    Start job monitoring threads and return coordinator for external signaling.

    Returns:
        tuple: (coordinator, monitor_function) where monitor_function() waits for completion
    """
    # Create coordinator for thread synchronization
    coordinator = OutputCoordinator()

    # Create threads for job monitoring and log tailing
    qstat_thread = threading.Thread(
        target=monitor_qstat, args=(job_id, quiet, coordinator, success_msg)
    )
    out_thread = threading.Thread(
        target=tail, args=(out_file, "STDOUT", coordinator), daemon=True
    )
    err_thread = threading.Thread(
        target=tail, args=(err_file, "STDERR", coordinator), daemon=True
    )

    # Set up signal handler for Ctrl-C
    def signal_handler(signum, frame):
        logging.info("Interrupt received, shutting down threads...")
        coordinator.signal_shutdown()

    original_handler = signal.signal(signal.SIGINT, signal_handler)

    # Start all threads
    qstat_thread.start()
    out_thread.start()
    err_thread.start()

    def wait_for_completion():
        """Wait for job monitoring to complete and return exit status"""
        try:
            # Wait for job monitoring to complete or shutdown signal
            qstat_thread.join()
            return coordinator.job_exit_status or 0
        finally:
            # Restore original signal handler
            signal.signal(signal.SIGINT, original_handler)

    return coordinator, wait_for_completion
```

**Coordination Flow**:
1. **Job Construction**: Progress message "🔧 Job command constructed"
2. **Submission**: qsub called, job ID returned
3. **Progress Update**: "✅ Job submitted successfully! Job ID: X"
4. **Thread Creation**: `start_job_monitoring()` creates coordinator and threads
5. **Submission Signal**: `coordinator.signal_submission_complete()`
6. **Monitor Activation**: Monitor thread starts polling qstat after waiting for submission_complete
7. **Spinner Start**: Spinner shows for 10 seconds maximum in monitor thread
8. **Output Detection**: Tail threads detect output and signal `output_started`
9. **Spinner Coordination**: Spinner waits for `output_started`, then clears

### Thread Responsibilities

#### 1. Job Monitor Thread (`monitor_qstat`)
```
Lifecycle: After submission through job completion
Purpose:  Monitor job completion via PBS qstat
Signals:  job_completed, job_exit_status
Waits:    submission_complete, shutdown_requested
```

**Behavior**: Waits for `submission_complete` event before starting qstat polling. This prevents the monitor from starting too early and ensures the job submission process is fully complete before monitoring begins.

- Polls `qstat` every 30 seconds to check job status
- When job completes, waits 5 seconds for PBS cleanup
- Polls for exit status every 5 seconds (up to 60 seconds)
- Stores final exit status in coordinator
- Signals `job_completed` to shut down other threads

#### 2. STDOUT Tail Thread (`tail`)
```
Lifecycle: Until EOF or job completion
Purpose:  Stream job output to terminal STDOUT
Signals:  output_started, eof_detected, spinner_cleared
Waits:    shutdown_requested
```

- Uses `tailer.follow()` to stream log file in real-time
- On first output line, signals `output_started` to clear spinner
- **Enhancement**: Sends carriage return (`\r`) to `/dev/tty` before streaming to clear leftover spinner chars
- Continues until EOF, job completion, or shutdown
- Signals `eof_detected` when log file ends

#### 3. STDERR Tail Thread (`tail`)
```
Lifecycle: Until EOF or job completion
Purpose:  Stream job errors to terminal STDERR
Signals:  output_started, eof_detected, spinner_cleared
Waits:    shutdown_requested
```

- Identical to STDOUT thread but writes to `sys.stderr`
- Both tail threads can signal `output_started` (whichever outputs first)
- **Enhancement**: Also sends carriage return (`\r`) to clear spinner artifacts

#### 4. Spinner Thread (`JobSpinner._spin`)
```
Lifecycle: Until job status changes or output starts
Purpose:  Show minimal progress indicator while waiting
Signals:  spinner_cleared
Waits:    output_started, job_running, job_finished, shutdown_requested
```

- Runs in its own daemon thread created by JobSpinner context manager
- Displays minimal spinner animation (no message by default)
- Created when monitor thread enters `with JobSpinner(...)` context
- Checks for stopping conditions every 0.1 seconds via `should_stop_spinner()`
- Stops automatically when:
  - Output starts (tail threads signal `output_started`)
  - Job starts running (monitor signals `job_running`)
  - Job finishes/fails (monitor signals `job_finished`)
  - Shutdown requested (Ctrl-C)
- Event-driven design eliminates arbitrary timeouts
- **Key Change**: Spinner is NOT used during job submission - only during monitoring phase

## Thread Lifecycle

Here's the complete sequence from job submission to completion:

```
1. Job Construction
   ├── Progress: "🔧 Job command constructed"
   ├── Job script created with execution context
   └── qsub command prepared

2. Job Submission
   ├── qsub executed directly (NO spinner during submission)
   ├── Progress: "✅ Job submitted successfully! Job ID: X" (single message)
   ├── Job ID stored in global _CURRENT_JOB_ID for signal handling
   └── Signal handler registered in execution.py

3. Thread Setup
   ├── start_job_monitoring() called with job_id, out_file, err_file
   ├── OutputCoordinator created
   ├── Monitor thread created (target=monitor_qstat)
   ├── STDOUT tail thread created (target=tail, daemon=True)
   ├── STDERR tail thread created (target=tail, daemon=True)
   ├── Local signal handler registered for thread coordination
   └── All threads started immediately (monitor, stdout, stderr)

4. Monitor Activation
   ├── monitor_qstat waits for submission_complete event
   ├── submission_complete triggered by coordinator.signal_submission_complete()
   ├── Monitor creates JobSpinner context (creates 4th thread - spinner daemon)
   ├── Monitor begins event-driven qstat polling loop
   ├── Spinner animates until job status change or output detected
   ├── Monitor signals job_running when status becomes "R"
   ├── Monitor signals job_finished when status becomes "F" or "H"
   └── STDOUT/STDERR tail threads follow log files throughout

5. Job Starts Running
   ├── Monitor detects status change from qstat polling
   ├── Monitor signals job_running event to coordinator
   ├── Spinner thread detects job_running via should_stop_spinner() and exits
   ├── Monitor displays "\r🚀 Job started running" (carriage return clears spinner)
   ├── Log files get first content (around same time or after)
   ├── Tail thread signals output_started on first line
   ├── Tail thread sends "\r" to /dev/tty to clear any leftover spinner chars
   └── Output streams to terminal in real-time

6. Job Completes
   ├── monitor_qstat detects completion via qstat (status F or H)
   ├── monitor_qstat waits for PBS cleanup (5s)
   ├── monitor_qstat polls for exit status (5s intervals, max 60s)
   ├── monitor_qstat signals job_completed and stores exit status
   └── All threads check should_shutdown() and exit

7. Cleanup
   ├── Tail threads reach EOF and signal eof_detected
   ├── qstat thread completes and wait_for_completion() returns
   ├── coordinator.job_exit_status returned to execution.py
   ├── _CURRENT_JOB_ID cleared
   └── qxub exits with job's exit code via sys.exit()
```

## Signal Flow Diagram

```
   ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
   │  Spinner Thread │    │ STDOUT Thread   │    │ STDERR Thread   │
   │                 │    │                 │    │                 │
   │ Waits for:      │    │ Signals:        │    │ Signals:        │
   │ • output_started│◄───┤ • output_started│    │ • output_started│───┐
   │ • job_running   │    │ • eof_detected  │    │ • eof_detected  │   │
   │ • job_finished  │    │ • spinner_cleared│    │ • spinner_cleared│   │
   │ • shutdown      │    │                 │    │                 │   │
   │                 │    │ Waits for:      │    │ Waits for:      │   │
   │ Signals:        │    │ • shutdown      │    │ • shutdown      │   │
   │ • spinner_cleared│    │                 │    │                 │   │
   └─────────────────┘    └─────────────────┘    └─────────────────┘   │
            │                       │                       │           │
            │                       │                       │           │
            ▼                       ▼                       ▼           │
   ┌─────────────────────────────────────────────────────────────────────┐  │
   │                    OutputCoordinator                                │  │
   │  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐     │  │
   │  │ output_started  │  │ spinner_cleared │  │ job_completed   │     │  │
   │  │ eof_detected    │  │ shutdown_req    │  │ submission_comp │     │  │
   │  │ job_running ⭐ │  │ job_finished ⭐│  │ job_exit_status │     │  │
   │  └─────────────────┘  └─────────────────┘  └─────────────────┘     │  │
   └─────────────────────────────────────────────────────────────────────┘  │
                                      ▲                                      │
                                      │ submission_complete                  │
                                      │ job_running ⭐                     │
                                      │ job_finished ⭐                    │
                              ┌─────────────────┐                           │
                              │ Monitor Thread  │                           │
                              │                 │                           │
                              │ Signals:        │                           │
                              │ • job_completed │                           │
                              │ • job_running ⭐│                           │
                              │ • job_finished ⭐│                           │
                              │ • job_exit_status│                           │
                              │                 │                           │
                              │ Waits for:      │                           │
                              │ • submission_comp│◄─────────────────────────┘
                              │ • shutdown      │
                              └─────────────────┘
                                      ▲
                                      │
                          ┌──────────────────────┐
                          │    Main Thread       │
                          │                      │
                          │ 1. Constructs job    │
                          │ 2. Submits via qsub  │
                          │    (NO spinner)      │
                          │ 3. Signals submission│
                          │    _complete         │
                          │ 4. Waits for monitor │
                          └──────────────────────┘

⭐ = New event-driven enhancements in v2.2
```

## Graceful Shutdown (Ctrl-C Handling)

When a user presses Ctrl-C, qxub performs graceful shutdown through multiple layers:

### Signal Handler in execution.py
The primary signal handler is installed in `execution.py` and tracks the current job globally:

```python
# Global variable to track current job for signal handling
_CURRENT_JOB_ID = None

def _signal_handler(signum, frame):
    """Handle SIGINT (Ctrl+C) by cleaning up submitted job"""
    if _CURRENT_JOB_ID:
        print("🛑 Interrupted! Cleaning up job...")
        from .scheduler import qdel
        success = qdel(_CURRENT_JOB_ID, quiet=False)
        if success:
            print("✅ Job cleanup completed")
        else:
            print(f"⚠️  Job cleanup failed - you may need to manually run: qdel {_CURRENT_JOB_ID}")
    print("👋 Goodbye!")
    sys.exit(130)  # Standard exit code for SIGINT
```

### Thread Coordination Handler
Additionally, `start_job_monitoring()` creates a local signal handler that coordinates thread shutdown:

```python
def signal_handler(signum, frame):
    logging.info("Interrupt received, shutting down threads...")
    coordinator.signal_shutdown()
```

### Shutdown Flow
1. **User presses Ctrl-C**: Both signal handlers may be triggered
2. **Job Cleanup**: Global handler in execution.py attempts to delete the job via `qdel`
3. **Thread Notification**: Local handler signals `coordinator.signal_shutdown()`
4. **Thread Response**: All threads check `should_shutdown()` in their loops
5. **Graceful Exit**: Threads stop processing and exit cleanly
6. **Exit Code**: Returns 130 (standard SIGINT exit code)

## Exit Code Propagation

One of the most important features is that qxub returns the actual job's exit code:

1. **Job Completion**: Monitor detects job finished via qstat
2. **PBS Cleanup**: Waits 5 seconds for PBS to finalize job state
3. **Exit Status Polling**: Queries qstat for Exit_status field every 5 seconds
4. **Status Storage**: Stores exit code in `coordinator.job_exit_status`
5. **Propagation**: Main thread returns this exit code via `sys.exit()`

This allows scripts to check `$?` and respond appropriately to job failures.

## Thread Safety Considerations

### Race Conditions Avoided

1. **Spinner vs Output**: Spinner waits for `output_started` before clearing
2. **Multiple Outputs**: Both STDOUT/STDERR can signal `output_started` safely
3. **Shutdown Coordination**: All threads check the same shutdown conditions
4. **Exit Status**: Only monitor thread writes `job_exit_status`

### Event-Based Synchronization

Using `threading.Event` objects prevents race conditions:
- Events are thread-safe and atomic
- Multiple threads can wait on the same event
- Setting an event is immediate and visible to all waiters

### Daemon Threads

Tail threads are marked as daemon threads:
```python
out_thread = threading.Thread(target=tail, args=(...), daemon=True)
```

This ensures they automatically terminate when the main thread exits, preventing hanging processes.

## Error Handling

### Thread Exceptions

Each thread has its own exception handling:
```python
try:
    # Thread work
    for line in tailer.follow(f):
        # Process line
except Exception as e:
    logging.error(f"Error in tail thread: {e}")
finally:
    logging.debug("Thread exiting")
```

### PBS Command Failures

Monitor thread handles qstat failures gracefully:
- Logs errors but continues monitoring
- Returns sensible defaults (status="C" for completed)
- Attempts fallback parsing for robustness

### File Access Issues

Tail threads handle missing/inaccessible log files:
- Create parent directories if needed
- Touch files to ensure they exist
- Use proper file encoding (UTF-8)

## Recent Improvements (v2.2)

### Event-Driven Spinner Control

The spinner system was enhanced from timeout-based to event-driven control:

**Before**: Spinner ran for fixed 10-second timeout regardless of job status
**After**: Spinner responds immediately to job status changes via events

**New Events**:
- `job_running` - Job status became "R" (Running)
- `job_finished` - Job status became "F" (Finished) or "H" (Held)
- `should_stop_spinner()` - Checks all stopping conditions in one call

**Benefits**:
- Immediate spinner shutdown when job starts or finishes
- No arbitrary delays waiting for timeouts
- More responsive user experience
- Cleaner event coordination

### Message Display Cleanup

Fixed multiple display issues for cleaner user experience:

**Duplicate Message Fix**:
- Removed redundant "Job submitted successfully" print statement
- Single clean submission message displayed

**Spinner Contamination Fix**:
- Removed `JobSpinner` from `qsub()` function during submission phase
- Spinner only runs during monitoring phase, not submission

**Clean Transitions**:
- Added carriage return (`\r`) to "🚀 Job started running" message
- Tail threads send `\r` before streaming output to clear spinner artifacts
- Proper spinner line clearing in context manager

**Result**: Clean, professional display without visual artifacts or duplicate messages

## Development Guidelines

### Adding New Signals

When adding new coordination events:

1. **Add to OutputCoordinator**:
   ```python
   self.my_new_event = threading.Event()
   ```

2. **Add signal method**:
   ```python
   def signal_my_event(self):
       self.my_new_event.set()
   ```

3. **Add wait method** (if needed):
   ```python
   def wait_for_my_event(self, timeout=None):
       return self.my_new_event.wait(timeout)
   ```

4. **Update should_shutdown()** (if relevant):
   ```python
   def should_shutdown(self):
       return (self.shutdown_requested.is_set() or
               self.my_new_event.is_set() or
               # ... other conditions)
   ```

### Thread Communication Patterns

**Don't**: Use shared variables between threads
```python
# BAD - race condition
shared_status = "running"
# Multiple threads modify this
```

**Do**: Use Events and coordinator
```python
# GOOD - thread-safe
coordinator.signal_status_changed()
if coordinator.wait_for_status_change(timeout=1):
    # Handle status change
```

### Critical Coordination Events

The `submission_complete` event is essential for proper timing and is now handled through the execution.py interface:

```python
# In execution.py after qsub and start_job_monitoring
job_id = qsub(submission_command, quiet=ctx_obj["quiet"])
success_msg = f"✅ Job submitted successfully! Job ID: {job_id}"
print_status(success_msg, final=True)

# Start job monitoring and get coordinator
coordinator, wait_for_completion = start_job_monitoring(
    job_id, out, err, quiet=ctx_obj["quiet"]
)

# Signal that submission messages are complete - enables monitor
coordinator.signal_submission_complete()

# Wait for job completion
exit_status = wait_for_completion()
sys.exit(exit_status)
```

In the monitor thread, the coordination works like this:

```python
def monitor_qstat(job_id, quiet=False, coordinator=None, success_msg=None):
    # Wait for job submission to complete before starting spinner
    if coordinator:
        coordinator.wait_for_submission_complete()

    # Start spinner for up to 10 seconds
    with JobSpinner("", show_message=False, quiet=quiet, coordinator=coordinator):
        time.sleep(10)

    # Begin main monitoring loop
    while True:
        # Check status and handle completion...
```

This prevents race conditions where monitoring starts before submission is fully complete.

### Testing Threading Code

Threading code is inherently harder to test. Use these strategies:

1. **Mock the coordinator** in unit tests
2. **Test individual thread functions** separately
3. **Use timeouts** to prevent hanging tests
4. **Test shutdown scenarios** explicitly

```python
def test_thread_shutdown():
    coordinator = OutputCoordinator()
    # Start thread
    thread = threading.Thread(target=my_function, args=(coordinator,))
    thread.start()

    # Signal shutdown
    coordinator.signal_shutdown()

    # Verify graceful exit
    thread.join(timeout=5)
    assert not thread.is_alive()
```

## Performance Characteristics

### Resource Usage

- **Single main thread**: Monitoring, output streaming, and progress indication
- **One background daemon thread**: Resource collection (only after job completion)
- **Low CPU**: Mostly sleep/wait on file I/O and status file polling
- **Memory**: Minimal - only line-by-line file reading with position tracking
- **No network calls**: Status monitoring via local file system (job scripts create status files)

### Latency

- **Output streaming**: Near real-time (0.2 second polling interval)
- **Job status detection**: 0.5 second polling interval for status files
- **Shutdown response**: < 1 second typically
- **Spinner updates**: 10Hz (every 0.1 seconds) during job startup

### Scalability

The current design scales well because:
- Thread count is minimal: 1 main + 1 background daemon (only post-completion)
- No thread pools or dynamic thread creation
- File-based status monitoring eliminates qstat polling overhead
- Clean shutdown prevents resource leaks
- Minimal shared state between threads
- Clean shutdown prevents resource leaks

## Debugging Threading Issues

### Logging

Enable debug logging to trace thread behavior:
```bash
export QXUB_LOG_LEVEL=DEBUG
qxub --env myenv -- script.py
```

Key log messages to watch for:
- `"Thread starting"`
- `"Output started - clearing spinner"`
- `"Thread shutdown requested"`
- `"Thread exiting"`

### Common Issues

1. **Hanging on exit**: Usually means a thread isn't checking `should_shutdown()`
2. **Missing output**: Check file permissions and paths
3. **Spinner not clearing**: Verify `output_started` is being signaled
4. **Wrong exit code**: Monitor thread may not be getting job status

### Manual Thread Inspection

You can inspect thread state in a debugger:
```python
import threading
print(f"Active threads: {threading.active_count()}")
for thread in threading.enumerate():
    print(f"Thread: {thread.name}, alive: {thread.is_alive()}")
```

## Future Improvements

### Potential Enhancements

1. **Configurable poll intervals**: Make 30s monitor interval configurable
2. **Progress indicators**: Show job queue position or estimated time
3. **Multiple job monitoring**: Support monitoring multiple jobs simultaneously
4. **Resume capability**: Reconnect to already-running jobs
5. **Better error recovery**: Retry failed qstat calls with backoff

### Architecture Considerations

The current design is optimized for simplicity and reliability. Any changes should preserve:
- **Clean shutdown behavior**
- **Proper exit code propagation**
- **Real-time output streaming**
- **Low resource usage**

The threading architecture provides a solid foundation that can be extended while maintaining backward compatibility and user experience.
