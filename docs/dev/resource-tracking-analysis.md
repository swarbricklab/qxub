# PBS Resource Tracking Analysis

## Overview
Analysis of `qstat -fx` output to understand what resource information can be captured for enhanced history logging.

## Resource Categories

### 1. Resource Requests (Resource_List.*)
**What was requested when submitting the job:**

```
Resource_List.jobfs = 104857600b           # Jobfs storage (100MB)
Resource_List.mem = 4294967296b            # Memory (4GB)
Resource_List.mpiprocs = 1                 # MPI processes
Resource_List.ncpus = 1                    # CPU cores
Resource_List.nodect = 1                   # Node count
Resource_List.place = free                 # Placement strategy
Resource_List.select = 1:ncpus=1:mpiprocs=1:mem=4294967296:job_tags=normal:jobfs=104857600
Resource_List.storage = scratch/a56+gdata/a56  # Storage access
Resource_List.walltime = 00:30:00          # Wall clock time limit
Resource_List.wd = 0                       # Working directory flag
```

### 2. Resource Usage (resources_used.*)
**What was actually consumed:**

```
resources_used.cpupercent = 48             # CPU utilization percentage
resources_used.cput = 00:00:02             # CPU time used
resources_used.jobfs = 0b                  # Jobfs storage used
resources_used.mem = 135580kb              # Memory used (132MB)
resources_used.ncpus = 1                   # CPUs used
resources_used.vmem = 135580kb             # Virtual memory used
resources_used.walltime = 00:00:04         # Wall clock time used
```

### 3. Execution Environment
**Where and how the job ran:**

```
exec_host = gadi-cpu-clx-2841/3            # Execution host/core
exec_vnode = (gadi-cpu-clx-2841:ncpus=1:mem=4194304kb:jobfs=102400kb)
queue = normal-exec                        # Queue used
session_id = 2891068                       # Session ID
```

### 4. Timing Information
**Job lifecycle timestamps:**

```
qtime = Sun Oct  5 10:25:10 2025          # Queue time (submitted)
stime = Sun Oct  5 10:25:15 2025          # Start time
etime = Sun Oct  5 10:25:10 2025          # End time
mtime = Sun Oct  5 10:25:30 2025          # Modified time
obittime = Sun Oct  5 10:25:30 2025       # Obituary time
```

### 5. Job Metadata
**Job identification and configuration:**

```
Job_Id = 151612099.gadi-pbs
Job_Name = test-job
Job_Owner = jr9959@10.9.0.35
job_state = F                             # Final state
Exit_status = 0                           # Exit code
project = a56                             # Project allocation
group_list = a56                          # Group membership
```

### 6. File Paths
**Input/output file locations:**

```
Error_Path = gadi.nci.org.au:/g/data/a56/software/qsub_tools/test-job.e151612099
Output_Path = gadi.nci.org.au:/g/data/a56/software/qsub_tools/test-job.log
jobdir = /home/913/jr9959                 # Job home directory
```

## Proposed History Schema Enhancement

### Current executions.yaml Structure
```yaml
- id: "unique-execution-id"
  timestamp: "2025-01-05T10:25:10Z"
  recipe_hash: "abc123..."
  command: ["conda", "run", "-n", "base", "echo", "Hello World"]
  exit_code: 0
  duration: 4.2
  job_id: "151612099"
```

### Enhanced Structure with Resources
```yaml
- id: "unique-execution-id"
  timestamp: "2025-01-05T10:25:10Z"
  recipe_hash: "abc123..."
  command: ["conda", "run", "-n", "base", "echo", "Hello World"]
  exit_code: 0
  duration: 4.2
  job_id: "151612099"

  # NEW: Resource tracking
  resources:
    requested:
      mem: "4294967296b"      # 4GB
      ncpus: 1
      walltime: "00:30:00"
      jobfs: "104857600b"     # 100MB
      storage: "scratch/a56+gdata/a56"

    used:
      mem: "135580kb"         # 132MB actual
      ncpus: 1
      walltime: "00:00:04"
      jobfs: "0b"
      cpupercent: 48
      cput: "00:00:02"
      vmem: "135580kb"

    efficiency:
      mem_efficiency: 3.2%    # 135580kb / 4194304kb * 100
      time_efficiency: 0.2%   # 4s / 1800s * 100
      cpu_efficiency: 48%     # cpupercent

  execution:
    queue: "normal-exec"
    exec_host: "gadi-cpu-clx-2841/3"
    exec_vnode: "gadi-cpu-clx-2841:ncpus=1:mem=4194304kb:jobfs=102400kb"
    project: "a56"
    session_id: 2891068

  timing:
    queued_at: "2025-10-05T10:25:10Z"     # qtime
    started_at: "2025-10-05T10:25:15Z"    # stime
    finished_at: "2025-10-05T10:25:30Z"   # mtime
    queue_wait: 5                         # stime - qtime
    execution_time: 15                    # mtime - stime
```

## Resource Efficiency Calculations

### Memory Efficiency
```python
mem_used_kb = parse_size_to_kb(resources_used.mem)        # 135580kb
mem_requested_kb = parse_size_to_kb(Resource_List.mem)    # 4194304kb
mem_efficiency = (mem_used_kb / mem_requested_kb) * 100   # 3.2%
```

### Time Efficiency
```python
time_used_seconds = parse_time_to_seconds(resources_used.walltime)  # 4s
time_requested_seconds = parse_time_to_seconds(Resource_List.walltime)  # 1800s
time_efficiency = (time_used_seconds / time_requested_seconds) * 100  # 0.2%
```

### CPU Efficiency
```python
cpu_efficiency = int(resources_used.cpupercent)  # 48%
```

## Implementation Priority

### Phase 1: Core Resource Capture
- `resources.requested`: mem, ncpus, walltime, jobfs, storage
- `resources.used`: mem, ncpus, walltime, jobfs, cpupercent, cput, vmem
- `execution`: queue, exec_host, project

### Phase 2: Efficiency Metrics
- Calculate and store efficiency percentages
- Add efficiency analysis to history CLI commands

### Phase 3: Advanced Analytics
- Resource usage trending
- Job efficiency recommendations
- Queue performance analysis

## Data Parsing Considerations

### Size Parsing
- PBS uses various units: `b`, `kb`, `mb`, `gb`, `4294967296b`, `104857600b`
- Need robust size parsing to bytes for comparison

### Time Parsing
- PBS uses `HH:MM:SS` format: `00:00:04`, `00:30:00`
- Convert to seconds for calculations

### Host Parsing
- `exec_host` format: `gadi-cpu-clx-2841/3` (hostname/core)
- `exec_vnode` format: `(hostname:ncpus=1:mem=4194304kb:jobfs=102400kb)`

## Benefits of Enhanced Resource Tracking

1. **Resource Optimization**: Identify over/under-provisioned jobs
2. **Cost Analysis**: Track resource efficiency across projects
3. **Performance Insights**: Understand job resource patterns
4. **Capacity Planning**: Historical resource usage trends
5. **Debugging**: Correlate failures with resource constraints
