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
                              │ 3. Call start_job_  │
                              │    monitoring()     │
                              │ 4. Start 3 threads  │
                              └─────────────────────┘
                                         │
                        ┌────────────────┼────────────────┐
                        ▼                ▼                ▼
              ┌─────────────────┐ ┌─────────────┐ ┌─────────────────┐
              │  Monitor Thread │ │ Tail Threads│ │ Spinner (in     │
              │                 │ │             │ │ Monitor Thread) │
              │ Waits for       │ │ Follow logs │ │                 │
              │ submission_     │ │ stream I/O  │ │ Shows brief     │
              │ complete, then  │ │ immediately │ │ animation for   │
              │ shows spinner   │ │             │ │ max 10 seconds  │
              │ for 10s, then   │ │             │ │                 │
              │ polls qstat     │ │             │ │                 │
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
               │        │     │ Clear Spinner      │
               │        │     │                     │
               │        │     │ Signal output_started
               │        │     │ Continue streaming  │
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
 │        │             │          │          │        (in Monitor)      │
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
 │        │             │ wait &   │          │          │ show spinner  │
 │        │             │ spinner  │          │          │ (max 10s)     │
 │        │             ├─ ─ ─ ─ ─ ┤          │          ├─ ─ ─ ─ ─ ─ ─ ┤
 │        │             │ poll     │ follow   │ follow   │               │
 │        │             │ qstat    │ out.log  │ err.log  │               │
 │        │             ├─ ─ ─ ─ ─ ┤─ ─ ─ ─ ─ ┤─ ─ ─ ─ ─ ┤               │
 │        │             │          │          │          │               │
 │        │             │          │ output!  │          │               │
 │        │             │          ├──────────┼──────────┼───────────────►
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

## Event Timeline

```
Time →   0s        30s       60s       90s      120s     150s
         │         │         │         │         │         │
         │         │         │         │         │         │
Monitor: ├─polling─┼─polling─┼─polling─┼─complete┼─cleanup─┼─exit
         │         │         │         │         │         │
STDOUT:  ├─waiting─┼─waiting─┼─output──┼─stream──┼─stream──┼─EOF
         │         │         │         │         │         │
STDERR:  ├─waiting─┼─waiting─┼─waiting─┼─waiting─┼─errors──┼─EOF
         │         │         │         │         │         │
Spinner: ├─active──┼─active──┼─clear───┼─────────┼─────────┼─────
         │         │         │         │         │         │
Events:  │         │         │ output_ │ job_    │         │ eof_
         │         │         │ started │ complete│         │ detected
```

## Control Flow for Different Scenarios

### Scenario 1: Normal Job Completion

```
Job Submit → start_job_monitoring() → Signal submission_complete →
Monitor waits → Brief Spinner (10s max) → Start qstat polling →
Output Start → Clear Spinner → Stream Output → Job Complete →
Get Exit Code → wait_for_completion() returns → Exit with Job Status
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
Job Submit → start_job_monitoring() → Monitor Start → Brief Spinner (10s) →
Poll Status (Q) → Continue Polling → Poll Status (Q) → ... →
Eventually Starts or User Interrupts
```

## Thread Communication Matrix

| Thread     | Signals                  | Waits For              | Reads From           |
|------------|--------------------------|------------------------|----------------------|
| Monitor    | job_completed            | submission_complete    | qstat commands       |
|            | job_exit_status          | shutdown_requested     |                      |
| STDOUT     | output_started           | shutdown_requested     | out.log file         |
|            | eof_detected             |                        |                      |
| STDERR     | output_started           | shutdown_requested     | err.log file         |
|            | eof_detected             |                        |                      |
| Spinner    | spinner_cleared          | output_started         | None (just displays) |
| (in Monitor)| (via JobSpinner)        | shutdown_requested     | (max 10s timeout)    |

## Memory and Resource Usage

```
┌─────────────────┐
│ Main Process    │
│ ┌─────────────┐ │ ≈ 50MB base Python
│ │OutputCoord  │ │ ≈ 1KB (events only)
│ └─────────────┘ │
│ ┌─────────────┐ │ ≈ 8KB each thread
│ │Monitor Thread│ │ (minimal stack)
│ │ + Spinner   │ │ (spinner context in same thread)
│ └─────────────┘ │
│ ┌─────────────┐ │
│ │STDOUT Thread│ │
│ └─────────────┘ │
│ ┌─────────────┐ │
│ │STDERR Thread│ │
│ └─────────────┘ │
└─────────────────┘
```

Total memory overhead: ~24KB for threading system (negligible)

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
