# Threading Flow Diagrams

This document provides visual representations of the qxub threading architecture to complement the main [threading-architecture.md](threading-architecture.md) documentation.

## Thread State Diagram

```
                                    qxub starts
                                         │
                                         ▼
                              ┌─────────────────────┐
                              │   Job Submitted     │
                              │                     │
                              │ 1. Create job command
                              │ 2. Submit via qsub  │
                              │    (NO SPINNER)     │
                              │ 3. Call start_job_  │
                              │    monitoring()     │
                              │ 4. Start 3 threads  │
                              └─────────────────────┘
                                         │
                        ┌────────────────┼────────────────┐
                        ▼                ▼                ▼
              ┌─────────────────┐ ┌─────────────┐ ┌─────────────────┐
              │  Monitor Thread │ │ Tail Threads│ │ Monitor enters  │
              │                 │ │             │ │ JobSpinner      │
              │ Waits for       │ │ Follow logs │ │ context manager │
              │ submission_     │ │ stream I/O  │ │                 │
              │ complete, then  │ │ immediately │ │ Creates 4th     │
              │ creates spinner │ │             │ │ thread (daemon) │
              │ thread via      │ │             │ │                 │
              │ context manager │ │             │ │ Event-driven    │
              │ then polls qstat│ │             │ │ spinner waits   │
              │                 │ │             │ │ for events      │
              └─────────────────┘ └─────────────┘ └─────────────────┘
                        │                │                │
                        │                │                │
               ┌────────┼────────────────┼────────────────┘
               │        │                │
               ▼        │                ▼
    ┌─────────────────┐ │     ┌─────────────────────┐
    │ Job Running?    │ │     │ Output detected?    │
    │                 │ │     │                     │
    │ Status: R/Q     │ │     │ First line of logs  │
    │ Continue poll   │ │     │ appears             │
    └─────────────────┘ │     └─────────────────────┘
               │        │                │
               │        │                ▼
               │        │     ┌─────────────────────┐
               │        │     │ Spinner Thread     │
               │        │     │ detects output_    │
               │        │     │ started, clears    │
               │        │     │ itself and exits   │
               │        │     └─────────────────────┘
               │        │                │
               ▼        │                ▼
    ┌─────────────────┐ │     ┌─────────────────────┐
    │ Job Complete?   │ │     │ Stream Real-time    │
    │                 │ │     │                     │
    │ Status: F/C     │ │     │ STDOUT → terminal   │
    │ Get exit code   │ │     │ STDERR → terminal   │
    └─────────────────┘ │     └─────────────────────┘
               │        │                │
               ▼        │                ▼
    ┌─────────────────┐ │     ┌─────────────────────┐
    │ Signal complete │ │     │ Reach EOF or       │
    │                 │ │     │ shutdown signal?    │
    │ job_completed   │ │     │                     │
    │ job_exit_status │ │     │ Signal eof_detected │
    └─────────────────┘ │     └─────────────────────┘
               │        │                │
               └────────┼────────────────┘
                        ▼
                ┌─────────────────┐
                │ All threads     │
                │ check shutdown  │
                │                 │
                │ should_shutdown()│
                │ returns True    │
                └─────────────────┘
                        │
                        ▼
                ┌─────────────────┐
                │ Graceful Exit   │
                │                 │
                │ Return job's    │
                │ exit status     │
                └─────────────────┘
```

## Sequence Diagram

```
User    execution.py  Monitor    STDOUT     STDERR     Spinner      Coordinator
 │        │             │          │          │        Thread            │
 │  cmd   │             │          │          │          │               │
 ├────────►             │          │          │          │               │
 │        │ create      │          │          │          │               │
 │        ├─────────────┼──────────┼──────────┼──────────┼───────────────►
 │        │             │          │          │          │               │
 │        │ start_job_monitoring()  │          │          │               │
 │        ├─────────────►          │          │          │               │
 │        ├────────────────────────►          │          │               │
 │        ├───────────────────────────────────►          │               │
 │        │             │          │          │          │               │
 │        │ signal_submission_complete()       │          │               │
 │        ├─────────────┼──────────┼──────────┼──────────┼───────────────►
 │        │             │          │          │          │               │
 │        │             │ enter    │          │          │               │
 │        │             │ JobSpinner          │          │               │
 │        │             │ context  │          │          │               │
 │        │             ├─ ─ ─ ─ ─ ┤          │          │               │
 │        │             │          │          │          │ create &      │
 │        │             │          │          │          │ start         │
 │        │             │          │          │          ├───────────────►
 │        │             │ sleep    │ follow   │ follow   │ animate       │
 │        │             │ 10s      │ out.log  │ err.log  │ spinner       │
 │        │             ├─ ─ ─ ─ ─ ┤─ ─ ─ ─ ─ ┤─ ─ ─ ─ ─ ┤─ ─ ─ ─ ─ ─ ─ ┤
 │        │             │          │          │          │               │
 │        │             │          │ output!  │          │               │
 │        │             │          ├──────────┼──────────┼───────────────►
 │        │             │          │          │          │ detect &      │
 │        │             │          │          │          │ clear         │
 │        │             │          │          │          ├─ ─ ─ ─ ─ ─ ─ ┤
 │◄───────┼─────────────┼──────────┤          │          │               │
 │ output │             │          │          │          │               │
 │        │             │          │          │ error!   │               │
 │        │             │          │          ├──────────┼───────────────►
 │◄───────┼─────────────┼──────────┼──────────┤          │               │
 │ error  │             │          │          │          │               │
 │        │             │ complete │          │          │               │
 │        │             ├──────────┼──────────┼──────────┼───────────────►
 │        │             │          │ EOF      │ EOF      │               │
 │        │             │          ├──────────┼──────────┼───────────────►
 │        │             │          ├─ ─ ─ ─ ─ ┤─ ─ ─ ─ ─ ┤               │
 │        │             │          │ shutdown │ shutdown │               │
 │        │ wait_for_completion()   │          │          │               │
 │        │ returns exit_status     │          │          │               │
 │◄───────┤             │          │          │          │               │
```

## Event-Driven Flow (v2.2 Enhancement)

```
Event Timeline - Enhanced with Immediate Status Response:

Monitor Thread            Spinner Thread           Tail Threads
      │                        │                        │
      │ qstat polling           │ waiting for            │ following logs
      │ status: Q               │ stop events            │ (empty)
      ▼                        ▼                        ▼
      │ status: R               │ detects job_running    │ first output
      ├─signal_job_running()────►─immediate stop──────────► signal_output_started()
      │                        │ clear & exit           │
      │ "\r🚀 Job started"      │                        │ "\r" clear spinner
      │                        │                        │ stream output
      ▼                        ▼                        ▼
      │ status: F               │ (exited)               │ continuing...
      ├─signal_job_finished()   │                        │
      │ wait for exit status    │                        │ EOF
      ▼                        ▼                        ▼
      │ job_exit_status         │                        │ signal_eof()
      │ signal_job_completed()  │                        │ thread exit
      │ thread exit             │                        │
      ▼                        ▼                        ▼
   [DONE]                   [DONE]                   [DONE]
```

**Key Improvements**:
- **Immediate Response**: Spinner stops as soon as job status changes, not after timeout
- **Clean Transitions**: Carriage returns (`\r`) overwrite spinner characters cleanly
- **Event Coordination**: All threads respond to the same events through OutputCoordinator
- **No Contamination**: Spinner only runs during monitoring, never during submission

## Event Timeline

```
Time →   0s        5s        15s       30s       45s      60s      75s
         │         │         │         │         │         │         │
Monitor: ├─wait────┼─polling─┼─R_detected─signal─┼─polling─┼─F_detected─cleanup─
         │         │         │         │         │         │         │
STDOUT:  ├─waiting─┼─waiting─┼─output──┼─stream──┼─stream──┼─stream──┼─EOF
         │         │         │         │         │         │         │
STDERR:  ├─waiting─┼─waiting─┼─waiting─┼─waiting─┼─errors──┼─stream──┼─EOF
         │         │         │         │         │         │         │
Spinner: ├─animate─┼─animate─┼─STOP────┼─(exit)──┼─────────┼─────────┼──────
         │         │         │ ↑       │         │         │         │
         │         │         │ └─event-driven stop    │         │         │

Legend:
- R_detected: Monitor detects job status "R" (Running)
- F_detected: Monitor detects job status "F" (Finished)
- signal: Monitor signals job_running event to stop spinner immediately
- STOP: Spinner detects event and stops (no timeout needed)
```

**Before v2.2**: Spinner ran for fixed 10s timeout regardless of job status
**After v2.2**: Spinner stops immediately when job status changes or output starts
         │         │         │         │         │         │         │
Spinner: ├─────────┼─active──┼─status──┼─────────┼─────────┼─────────┼─────
         │         │         │ change  │         │         │         │
         │         │         │ detected│         │         │         │
Events:  │         │ submit_ │ job_    │ output_ │ job_    │         │ eof_
         │         │ complete│ running │ started │ complete│         │ detected
         │         │         │ OR      │         │         │         │
         │         │         │ output  │         │         │         │
```

## Control Flow for Different Scenarios

### Scenario 1: Normal Job Completion

```
Job Submit → start_job_monitoring() → Signal submission_complete →
Monitor waits → Create Spinner Thread → Monitor starts qstat polling →
Spinner animates until status change → Job status becomes "R" →
Monitor signals job_running → Spinner exits → Stream Output →
Job Complete → Get Exit Code → wait_for_completion() returns → Exit with Job Status
```

### Scenario 2: User Interruption (Ctrl-C)

```
Job Submit → start_job_monitoring() → Monitor Start → User Ctrl-C →
execution.py Signal Handler → qdel job → coordinator.signal_shutdown() →
All Threads Check → Graceful Stop → Exit 130
```

### Scenario 3: Job Failure

```
Job Submit → Monitor Start → Spinner Show → Output Start →
Clear Spinner → Stream Errors → Job Complete (Failed) →
Get Exit Code (1) → Signal Threads → Cleanup → Exit 1
```

### Scenario 4: Job Never Starts (Queue Wait)

```
Job Submit → start_job_monitoring() → Monitor Start → Create Spinner Thread →
Poll Status & Signal Changes → Status becomes "R" → Signal job_running →
Spinner exits → Continue Polling → Eventually Completes
```

## Thread Communication Matrix

| Thread     | Signals                  | Waits For              | Reads From           |
|------------|--------------------------|------------------------|----------------------|
| Monitor    | job_completed            | submission_complete    | qstat commands       |
|            | job_running              | shutdown_requested     |                      |
|            | job_finished             |                        |                      |
|            | job_exit_status          |                        |                      |
| STDOUT     | output_started           | shutdown_requested     | out.log file         |
|            | eof_detected             |                        |                      |
| STDERR     | output_started           | shutdown_requested     | err.log file         |
|            | eof_detected             |                        |                      |
| Spinner    | spinner_cleared          | job_running            | None (just displays) |
| (daemon)   | (via JobSpinner)         | job_finished           | (created by monitor) |
|            |                          | output_started         |                      |
|            |                          | shutdown_requested     |                      |

## Memory and Resource Usage

```
┌─────────────────┐
│ Main Process    │
│ ┌─────────────┐ │ ≈ 50MB base Python
│ │OutputCoord  │ │ ≈ 1KB (events only)
│ └─────────────┘ │
│ ┌─────────────┐ │ ≈ 8KB each thread
│ │Monitor Thread│ │ (minimal stack)
│ └─────────────┘ │
│ ┌─────────────┐ │
│ │STDOUT Thread│ │
│ └─────────────┘ │
│ ┌─────────────┐ │
│ │STDERR Thread│ │
│ └─────────────┘ │
│ ┌─────────────┐ │ ≈ 8KB (short-lived)
│ │Spinner Thread│ │ (daemon, max 10s)
│ │  (optional) │ │ (created by monitor)
│ └─────────────┘ │
└─────────────────┘
```

Total memory overhead: ~32KB for threading system (negligible)

## Error Propagation Flow

```
                PBS Job Exit Code
                        │
                        ▼
              ┌─────────────────────┐
              │ Monitor detects     │
              │ job completion      │
              │ via qstat           │
              └─────────────────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │ Wait 5s for         │
              │ PBS cleanup         │
              └─────────────────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │ Poll for Exit_status│
              │ field every 5s      │
              │ (max 60s)           │
              └─────────────────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │ Store in            │
              │ coordinator.        │
              │ job_exit_status     │
              └─────────────────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │ wait_for_completion │
              │ returns exit code   │
              └─────────────────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │ sys.exit(exit_code) │
              │ in execution.py     │
              └─────────────────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │ Shell sees correct  │
              │ exit code in $?     │
              └─────────────────────┘
```
