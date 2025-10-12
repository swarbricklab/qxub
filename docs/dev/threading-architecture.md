# qxub Threading Architecture

qxub uses a sophisticated multi-threaded architecture to provide real-time job monitoring with clean output streaming. This document explains how the threading system works, why it was designed this way, and how to work with it.

## Overview

The threading system coordinates four key components:
1. **Job Status Monitor** - Polls PBS for job completion
2. **STDOUT Tail Thread** - Streams job output to terminal
3. **STDERR Tail Thread** - Streams job errors to terminal
4. **Spinner Thread** - Shows progress indicator until output starts

All threads communicate through a central `OutputCoordinator` that manages their lifecycle and synchronization.

## Why Threading?

Without threading, qxub would need to choose between:
- **Polling-only**: Check job status every 30 seconds, miss real-time output
- **Tailing-only**: Stream output but miss job completion/cleanup

Threading allows qxub to do both simultaneously while providing a smooth user experience.

## Core Components

### OutputCoordinator

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
        self.job_exit_status = None                   # Final job exit code
```

**Events are boolean flags that threads can wait for or signal**:
- `wait()` - Block until the event is set
- `set()` - Signal that the event has occurred
- `is_set()` - Check if the event has occurred (non-blocking)

### Coordinated Job Submission

The `start_job_monitoring()` function in `scheduler.py` provides coordinated thread startup:

```python
def start_job_monitoring(job_id, out_log, err_log):
    """Start job monitoring with proper coordination."""
    coordinator = OutputCoordinator()

    # Start threads but monitor waits for submission_complete
    monitor_thread = start_monitor_thread(job_id, coordinator)
    stdout_thread = start_stdout_tail(out_log, coordinator)
    stderr_thread = start_stderr_tail(err_log, coordinator)
    spinner_thread = start_spinner_thread(job_id, coordinator)

    def wait_for_completion():
        return monitor_thread.join()

    return coordinator, wait_for_completion
```

**Coordination Flow**:
1. **Job Construction**: Progress message "🔧 Job command constructed"
2. **Submission**: qsub called, job ID returned
3. **Progress Update**: "✅ Job submitted successfully! Job ID: X"
4. **Submission Signal**: `coordinator.signal_submission_complete()`
5. **Monitor Activation**: Monitor thread starts polling qstat
6. **Output Detection**: Tail threads detect output and signal `output_started`
7. **Spinner Coordination**: Spinner waits for `output_started`, then clears

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
Signals:  output_started, eof_detected
Waits:    shutdown_requested
```

- Uses `tailer.follow()` to stream log file in real-time
- On first output line, signals `output_started` to clear spinner
- Continues until EOF, job completion, or shutdown
- Signals `eof_detected` when log file ends

#### 3. STDERR Tail Thread (`tail`)
```
Lifecycle: Until EOF or job completion
Purpose:  Stream job errors to terminal STDERR
Signals:  output_started, eof_detected
Waits:    shutdown_requested
```

- Identical to STDOUT thread but writes to `sys.stderr`
- Both tail threads can signal `output_started` (whichever outputs first)

#### 4. Spinner Thread (`JobSpinner._spin`)
```
Lifecycle: Until output starts or job completes
Purpose:  Show progress indicator while waiting
Signals:  spinner_cleared
Waits:    output_started
```

- Displays animated spinner: "🚀 Job 1234abcd - Waiting ⠋"
- Polls for `output_started` every 0.1 seconds
- When output starts, clears spinner line and signals `spinner_cleared`
- Automatically stops when job completes

## Thread Lifecycle

Here's the complete sequence from job submission to completion:

```
1. Job Construction
   ├── Progress: "🔧 Job command constructed"
   ├── Job script created with execution context
   └── qsub command prepared

2. Job Submission
   ├── qsub executed, job ID returned
   ├── Progress: "✅ Job submitted successfully! Job ID: X"
   ├── coordinator.signal_submission_complete() called
   └── All monitoring threads started but monitor waits

3. Monitor Activation
   ├── monitor_qstat waits for submission_complete event
   ├── submission_complete triggered → monitor begins polling
   ├── STDOUT tail thread starts following out.log
   ├── STDERR tail thread starts following err.log
   └── Spinner thread starts: "🚀 Job abc123 - Waiting ⠋"

4. Job Starts Running
   ├── Log files get first content
   ├── Tail thread signals output_started
   ├── Progress: "🚀 Job started running"
   ├── Spinner clears and signals spinner_cleared
   └── Output streams to terminal in real-time

5. Job Completes
   ├── monitor_qstat detects completion via qstat
   ├── monitor_qstat waits for PBS cleanup (5s)
   ├── monitor_qstat polls for exit status (5s intervals)
   ├── monitor_qstat signals job_completed
   └── All threads check should_shutdown() and exit

6. Cleanup
   ├── Tail threads reach EOF and signal eof_detected
   ├── Main thread waits for monitor completion
   ├── coordinator.job_exit_status returned
   └── qxub exits with job's exit code
```

## Signal Flow Diagram

```
   ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
   │  Spinner Thread │    │ STDOUT Thread   │    │ STDERR Thread   │
   │                 │    │                 │    │                 │
   │ Waits for:      │    │ Signals:        │    │ Signals:        │
   │ • output_started│◄───┤ • output_started│    │ • output_started│───┐
   │                 │    │ • eof_detected  │    │ • eof_detected  │   │
   │ Signals:        │    │                 │    │                 │   │
   │ • spinner_cleared   │    │ Waits for:      │    │ Waits for:      │   │
   └─────────────────┘    │ • shutdown      │    │ • shutdown      │   │
            │              └─────────────────┘    └─────────────────┘   │
            │                       │                       │           │
            │                       │                       │           │
            ▼                       ▼                       ▼           │
   ┌─────────────────────────────────────────────────────────────────────┐  │
   │                    OutputCoordinator                                │  │
   │  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐     │  │
   │  │ output_started  │  │ spinner_cleared │  │ job_completed   │     │  │
   │  │ eof_detected    │  │ shutdown_req    │  │ submission_comp │     │  │
   │  │ job_exit_status │  │                 │  │                 │     │  │
   │  └─────────────────┘  └─────────────────┘  └─────────────────┘     │  │
   └─────────────────────────────────────────────────────────────────────┘  │
                                      ▲                                      │
                                      │ submission_complete                  │
                                      │                                      │
                              ┌─────────────────┐                           │
                              │ Monitor Thread  │                           │
                              │                 │                           │
                              │ Signals:        │                           │
                              │ • job_completed │                           │
                              │ • job_exit_status                           │
                              │                 │                           │
                              │ Waits for:      │                           │
                              │ • submission_comp │◄─────────────────────────┘
                              │ • shutdown      │
                              └─────────────────┘
                                      ▲
                                      │
                          ┌──────────────────────┐
                          │    Main Thread       │
                          │                      │
                          │ 1. Constructs job    │
                          │ 2. Submits via qsub  │
                          │ 3. Signals submission│
                          │    _complete         │
                          │ 4. Waits for monitor │
                          └──────────────────────┘
```

## Graceful Shutdown (Ctrl-C Handling)

When a user presses Ctrl-C, qxub performs graceful shutdown:

1. **Signal Handler**: Catches SIGINT and calls `coordinator.signal_shutdown()`
2. **Thread Notification**: All threads check `should_shutdown()` in their loops
3. **Monitor Cleanup**: Monitor thread stops polling and attempts job cleanup
4. **Tail Cleanup**: Tail threads stop following log files
5. **Spinner Cleanup**: Spinner clears and stops animation
6. **Exit Code**: Returns 130 (standard SIGINT exit code)

```python
def signal_handler(signum, frame):
    logging.info("Interrupt received, shutting down threads...")
    coordinator.signal_shutdown()

# In each thread loop:
if coordinator and coordinator.should_shutdown():
    break  # Exit thread gracefully
```

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

The `submission_complete` event is essential for proper timing:

```python
# In execution.py after qsub
job_id = subprocess.run(qsub_cmd, ...)
print_tty(f"✅ Job submitted successfully! Job ID: {job_id}")
coordinator.signal_submission_complete()  # Critical: enables monitor

# In scheduler.py monitor thread
def monitor_qstat(job_id, coordinator):
    # Wait for submission to complete before monitoring
    coordinator.submission_complete.wait()

    while not coordinator.should_shutdown():
        # Now safe to start qstat polling
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

- **4 threads maximum**: Monitor + 2 tails + spinner
- **Low CPU**: Threads mostly sleep/wait on I/O
- **Memory**: Minimal - only line-by-line file reading
- **Network**: Only qstat calls every 30 seconds

### Latency

- **Output streaming**: Near real-time (tailer library)
- **Job completion**: Up to 30 seconds detection lag
- **Shutdown response**: < 1 second typically
- **Spinner updates**: 10Hz (every 0.1 seconds)

### Scalability

The current design scales well because:
- Thread count is fixed regardless of job size
- No thread pools or dynamic thread creation
- Minimal shared state between threads
- Clean shutdown prevents resource leaks

## Debugging Threading Issues

### Logging

Enable debug logging to trace thread behavior:
```bash
export QXUB_LOG_LEVEL=DEBUG
qxub conda --env myenv script.py
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
