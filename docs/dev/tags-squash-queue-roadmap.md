# Tags, Job Squashing, and Queue Management: Roadmap

**Status**: Planning — not yet implemented
**Last updated**: 2026-02-25

This document captures the agreed design for three interconnected features:

1. **Tags** — label jobs for filtering, reporting, and squash matching
2. **Job squashing** — bundle multiple small jobs into fewer, larger PBS jobs
3. **Queue management** — limit concurrent PBS submissions and guarantee dispatch

These features are primarily motivated by Snakemake workflow integration, where a single
workflow may generate hundreds or thousands of jobs that are expensive to schedule
individually.

---

## Implementation Phases

| Phase | Scope | Key deliverables |
|-------|-------|-----------------|
| 1 | Tags | `--tag`/`--tags` on `exec` and `interactive`, stored in DB + history, filterable |
| 2 | Shared DB + status-in-DB | Replace exit-code files with DB writes, virtual job IDs, dispatch-on-status-check |
| 3 | Queue management | Headroom limit, pending entries, `qxub queue` commands |
| 4 | Serial squashing | Squash conditions config, serial mode |
| 5 | Parallel squashing | Parallel mode |

Phases are designed to be independently shippable and useful.

---

## Phase 1: Tags

### CLI

```bash
qxub exec --env myenv --tag rule=align --tag workflow=brca -- bwa mem ref.fa reads.fastq
qxub exec --env myenv --tags "rule=align,workflow=brca,sample=S001" -- bwa mem ...
```

`--tag` accepts a single tag. `--tags` accepts a comma-separated list. Both can be used
together and are merged. A tag without `=` is a plain label; a tag with `=` is a key-value
pair. Both are treated the same in storage and queries.

**Snakemake cluster-generic integration**:

The recommended pattern uses `{rule}` (Snakemake built-in) plus `{resources.workflow}`,
where `workflow` is set once via `--default-resources` in the profile — no per-rule
annotation required:

```bash
# cluster-generic profile config.yaml
submit-cmd: "qxub exec --env {resources.conda_env} [...] --tags \"workflow={resources.workflow},rule={rule}\" -- {exec_job_cmd}"

default-resources:
  workflow: "my_workflow_name"
```

This stamps every job with its rule name and the workflow name automatically.
Additional per-rule tags can be added by defining named resources in individual rules
(e.g. `resources: sample=wildcards.sample`) and extending the `--tags` string.

`--squashable` is a Phase 4 concern and is not part of Phase 1.

### Storage

Tags are stored as a JSON array of strings in a new `tags TEXT` column in
`job_resources` (defaults to `'[]'`). Example stored value:
`'["rule=align", "workflow=brca", "sample=S001"]'`.

The same `tags` field is added to execution records in the history JSON files.

Migration: `_migrate_database_schema` adds the column if absent.

Tags are written to whichever DB path is active (shared if configured, else per-user
fallback — see Phase 2). This means teams using a shared DB get tags in the right place
from Phase 1 onwards.

### Querying

```bash
qxub resources list --tag rule=align
qxub history executions --tag workflow=brca
```

Implemented via SQLite `json_each()`:

```sql
SELECT * FROM job_resources
WHERE EXISTS (
    SELECT 1 FROM json_each(tags)
    WHERE value = 'rule=align'
)
```

For key-prefix queries (e.g. "all executions with any `rule=` tag"):

```sql
WHERE EXISTS (
    SELECT 1 FROM json_each(tags)
    WHERE value LIKE 'rule=%'
)
```

---

## Phase 2: Shared Project DB and Status-in-DB

### Motivation

The current file-based status system (`final_exit_code_{job_id}`, `job_started_{job_id}`)
will generate enormous file counts at workflow scale (Snakemake can submit tens of thousands
of jobs). Writing status to a shared database eliminates this entirely.

### Database layout

A single unified database holds all tables: job resources, queue state, project-level
job tracking, and metadata. The active DB path is resolved at runtime:

1. `shared_db.path` from config (system, user, or project config) — preferred for teams
2. `~/.config/qxub/qxub.db` — per-user fallback when no shared path is configured

This means solo users get a working setup with no configuration, while teams point
everyone at a shared path on scratch.

```yaml
# /etc/xdg/qxub/config.yaml  (system config, sets shared DB for all users on project)
shared_db:
  path: /scratch/a56/qxub/qxub.db
  job_limit: 800  # max concurrent PBS jobs for this project
```

The DB is accessed identically regardless of path — all queries target `shared_db.path`
if set, else the user fallback. Commands like `qxub resources list` just reflect whichever
DB is active; `--project` filtering is only meaningful when using a shared DB that contains
other users' jobs.

> **Migration note**: The existing per-user `resources.db` is renamed `qxub.db` and all
> new tables added to it in the same migration step.

### Job script changes

Job scripts gain a `QXUB_SHARED_DB` environment variable (set by `qxub exec`). The job
script writes to the shared DB on start, completion, and failure instead of creating
status files:

```bash
# Instead of:  touch $STATUS_DIR/job_started_${PBS_JOBID}
python3 -c "from qxub.queue.db import mark_running; mark_running('${PBS_JOBID}')"

# Instead of:  echo $EXIT_CODE > $STATUS_DIR/final_exit_code_${PBS_JOBID}
python3 -c "from qxub.queue.db import mark_complete; mark_complete('${PBS_JOBID}', ${EXIT_CODE})"
```

Stdout (`.out`), stderr (`.err`), and PBS log (`.log`) files are **kept** — those have
value beyond exit status. Only the dedicated exit-code and job-started files are removed.

### Fallback

If `QXUB_SHARED_DB` is unset (old jobs, manual qsub, jobs submitted outside qxub), the
existing file-based status check is retained as a fallback within `job_status_from_files`.

### Virtual job IDs

Jobs enqueued in the qxub queue (pending dispatch) return a virtual job ID immediately:

```
qx-550e8400-e29b-41d4-a716-446655440000
```

`qxub status check qx-...` resolves this:
- If `status = 'pending'` in queue table → return `Q`
- If `status = 'dispatched'` → resolve to PBS job ID, return real PBS status from DB
- If `status = 'cancelled'` → return `F`

Snakemake never sees the difference between a virtual ID and a real PBS ID.

### Dispatch-on-status-check

Snakemake polls `qxub status check <job_id>` every ~10 seconds for every running job. This
is a natural trigger for queue management: on each status-check call, qxub also:

1. Reads current `active_count` from the shared DB (updated when jobs complete)
2. If `active_count < job_limit`: dispatch any pending queue entries for this project
3. For pending squash groups: check if timeout has expired and flush partial batches
4. Returns the requested job status normally

This eliminates the need for a daemon or cron job. As long as any workflow job is being
monitored by Snakemake, the queue will be flushed. `qxub queue flush` provides a manual
override for interactive use.

### Shared DB schema (unified)

```sql
-- /scratch/{project}/qxub/qxub.db  (or ~/.config/qxub/qxub.db for per-user fallback)

CREATE TABLE IF NOT EXISTS queue (
    entry_id        TEXT PRIMARY KEY,   -- 'qx-{uuid}'
    submitted_at    TEXT NOT NULL,
    user            TEXT NOT NULL,
    tags            TEXT DEFAULT '[]',  -- JSON array of tag strings
    squashable      INTEGER DEFAULT 0,
    squash_condition TEXT,              -- matched condition name, NULL if none matched
    squash_batch_id TEXT,               -- entry_id of the batch entry this was absorbed into
    command_b64     TEXT NOT NULL,      -- base64-encoded command
    working_dir     TEXT,
    exec_context    TEXT,               -- JSON: {type, env/mod/sif}
    resources       TEXT,               -- JSON: {mem_mb, cpus, walltime_sec, ...}
    status          TEXT DEFAULT 'pending', -- pending | dispatched | cancelled | completed | failed
    pbs_job_id      TEXT,               -- real PBS job ID once dispatched
    dispatched_at   TEXT,
    completed_at    TEXT,
    exit_code       INTEGER
);

CREATE TABLE IF NOT EXISTS project_jobs (
    pbs_job_id      TEXT PRIMARY KEY,
    entry_id        TEXT,               -- qxub virtual ID, NULL for non-qxub jobs
    user            TEXT NOT NULL,
    submitted_at    TEXT NOT NULL,
    started_at      TEXT,
    completed_at    TEXT,
    status          TEXT DEFAULT 'submitted',  -- submitted | running | completed | failed
    exit_code       INTEGER,
    tags            TEXT DEFAULT '[]',  -- copied from queue entry
    resources_requested TEXT,           -- JSON
    resources_used  TEXT                -- JSON, filled in post-completion
);

CREATE TABLE IF NOT EXISTS meta (
    key             TEXT PRIMARY KEY,
    value           TEXT,
    updated_at      TEXT
);
-- e.g. meta: key='active_count', value='42', updated_at='2026-02-21T10:00:00'
```

---

## Phase 3: Queue Management

### Headroom limit

Configures a project-wide maximum concurrent PBS jobs (default 800 of the 1000 PBS limit).
Jobs submitted beyond this are held in the `queue` table with `status = 'pending'` and
dispatched as headroom becomes available.

`active_count` is maintained in the `meta` table:
- Incremented on each successful `qsub` call
- Decremented when a job writes its completion status to `project_jobs`
- Reconciled by `qxub queue sync` which calls `qstat -P {project}` to recount

Since quotas are per-project, the count covers **all users** on the project. Non-qxub
jobs are invisible to the counter until `qxub queue sync` is run. Configure conservatively
(e.g. 800) to leave headroom for jobs not managed by qxub.

```bash
qxub queue sync        # reconcile active_count with real qstat output
qxub queue sync --dry  # show discrepancy without updating
```

### `qxub queue` commands

```bash
qxub queue list              # pending entries (all users on project)
qxub queue list --user me    # own pending entries only
qxub queue stats             # active_count, pending count, headroom, squash group status
qxub queue flush             # dispatch any pending entries within headroom
qxub queue flush --group <condition-name>  # force-flush a specific squash group
qxub queue cancel <entry_id> # cancel a pending entry
qxub queue sync              # reconcile DB count with qstat
```

The difference between `qxub queue list` and `qxub status list`:

- `qxub status` — PBS-level view: what PBS knows about (queued/running/done in scheduler)
- `qxub queue` — qxub-level view: entries not yet submitted to PBS (pending dispatch)

`qxub status list` will additionally show qxub-pending entries with status `PENDING (qxub)`
so users see the full picture in one command.

---

## Phase 4: Serial Squashing

Serial squashing is implemented before parallel because:
- It's simpler (no thread pool, no resource multiplication)
- It provides immediate value for the most common case: many fast, cheap jobs
- The squash condition config and queue machinery are shared with parallel

### Squash conditions

Defined in config files, following the existing system/user/project hierarchy:

```yaml
# ~/.config/qxub/config.yaml  or  /etc/xdg/qxub/config.yaml
squash_conditions:
  - name: snakemake-fast-align
    match_tags:                    # ALL tags must be present in the job's tags
      - "rule=align"
      - "workflow=brca"
    mode: serial                   # serial | parallel
    max_factor: 50                 # max tasks per PBS job
    timeout_sec: 120               # flush partial batch after this many seconds
    on_failure: worst_exit         # worst_exit | fail_fast
    # For serial mode: walltime = sum(task_walltimes), resources = single task
    # For parallel mode: walltime = max(task_walltimes), resources = N × max(per-task resource)
```

A job is considered squashable if at least one squash condition matches its tags.
`--squashable` is introduced in this phase as a CLI flag (accepting `true`/`false` to
support Snakemake resource interpolation). A job with `--squashable true` but no matching
condition is submitted immediately — no silent failure, no hanging in the queue.

### Walltime aggregation

- **Serial**: `walltime = sum(walltime_i)` — no buffer added; each task's walltime estimate
  is its own responsibility
- **Parallel**: `walltime = max(walltime_i)` — see Phase 5

### on_failure behaviour

| Value | Behaviour |
|-------|-----------|
| `worst_exit` | All tasks run to completion; PBS job exits with the worst (highest) exit code among tasks |
| `fail_fast` | Stop on first task failure; remaining tasks are not run |

No `ignore` mode — `worst_exit` already provides continue-on-failure semantics with a
meaningful exit code for Snakemake.

Per-task exit codes are always recorded in the shared DB regardless of the container job's
final exit code.

### Serial squash job script

A serial squash job script:
1. Reads a task manifest (list of base64-encoded sub-commands, working dirs, contexts)
2. Runs each task sequentially in a bash loop
3. Records per-task exit code to shared DB after each task
4. Exits with `max(exit_codes)` on completion, or immediately on `fail_fast`

Per-task output goes to permanent `.out`/`.err` files (named by virtual job ID so users
can retrieve them with `qxub history out qx-...`).

### Per-task resource tracking

PBS job logs report aggregate resources for the container job, not for individual tasks.
To capture per-task wall time and memory, the squash job script wraps each sub-command
with `/usr/bin/time -v` (GNU time, available on Gadi/Linux):

```bash
/usr/bin/time -v bash -c "$(echo $TASK_CMD_B64 | base64 -d)" 2> "${STATUS_DIR}/time_${ENTRY_ID}.txt"
```

After each task completes, the script parses the `time` output and writes per-task
`resources_used` (wall seconds, max RSS) to the DB alongside the exit code. This is the
only reliable source of per-task memory usage in a squashed job.

### Matching and batching logic (at enqueue time)

```
1. INSERT new entry into queue (status='pending', squashable=True, tags=[...])
2. Acquire SQLite WAL exclusive write lock (busy_timeout=5000ms)
3. Find matching squash condition (first match wins)
4. If no condition matches: dispatch immediately (submit to PBS now)
5. Count pending entries in same condition group
6. If count >= max_factor OR oldest_entry age >= timeout_sec:
     a. Collect up to max_factor pending entries
     b. Build squash job script + task manifest
     c. Submit single PBS job
     d. Update all entries: status='dispatched', pbs_job_id=<batch PBS ID>
7. Else: leave pending (dispatch-on-status-check or manual flush will handle it)
```

---

## Phase 5: Parallel Squashing

Parallel squashing extends serial with:
- Python `concurrent.futures.ThreadPoolExecutor` for task execution
- Resource aggregation: `N × max(resource_i)` for mem and cpus
- Walltime: `max(walltime_i)` across the batch

Resource aggregation uses `max × N` rather than `sum × N` because individual task resource
requests already include headroom. This slightly over-allocates but prevents any task from
hitting memory limits due to variance.

The `resource_aggregation` field in squash conditions is reserved for future use if
`sum × N` is preferable for known-tight requests.

---

## Cross-cutting concerns

### Config hierarchy for squash conditions

Squash conditions obey the same resolution order as other config:

```
CLI args > per-user config > system config > defaults
```

Project-level config (planned) would sit between user and system. This means an
HPC system admin can define standard squash conditions in `/etc/xdg/qxub/config.yaml`
that apply to all users, while users can override or extend them.

### Cross-user resource tracking

When a shared DB is configured, all users on the project write to the same database.
The `project_jobs` and `job_resources` tables therefore contain all users' jobs, enabling
grant reporting and team-wide efficiency analysis:

```bash
qxub resources stats --project   # project-wide summary (shared DB only)
qxub resources list --project    # all users' recent jobs (shared DB only)
qxub resources list              # own jobs only (works with shared or per-user DB)
```

When no shared DB is configured, all commands operate on the per-user fallback DB
(`~/.config/qxub/qxub.db`). The `--project` flag produces a warning in this case.

### `qxub queue sync` and non-qxub jobs

`qxub queue sync` runs `qstat -P {project}` and reconciles `active_count` in the `meta`
table. It also ingests non-qxub jobs into `project_jobs` with `entry_id = NULL`, so
the headroom calculation correctly accounts for jobs submitted outside qxub.

### SQLite concurrency

The shared DB uses WAL mode with `PRAGMA busy_timeout = 5000`. This supports O(10)
simultaneous writers (concurrent `qxub exec` or status-check calls from workflow nodes)
without contention.

---

## Open / deferred decisions

- **Project-level config**: Schema is being designed; squash conditions should work
  naturally once project config is available.
- **Snakemake status script**: The `qxub status check` command will need to be confirmed
  as the recommended pattern in the Snakemake integration docs once Phase 2 is complete.
- **Parallel resource aggregation**: `max × N` is the default. `sum × N` may be added
  as `resource_aggregation: sum_times_n` in a squash condition if use cases emerge.
- **Squash condition templating**: Condition names and tag patterns may support
  template variables (e.g. `{user}`) in a future iteration.
