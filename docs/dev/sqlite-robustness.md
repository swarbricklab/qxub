# Making the job database robust on shared filesystems

Status: **implemented**. All four phases are merged (#132, #133, #135, #136) and released in 3.8.0 and 3.9.0. Revised after the Gadi tests in the #132 review, and to separate live state from resource statistics.

## Why

The job database (`qxub.db`) misbehaves on NCI Gadi when many jobs run at once. We have hit, in order:

| Issue | Symptom | What we did |
|---|---|---|
| #61 | `qxtat check` exits 1 under load, aborting Snakemake | Distinguish DB errors from missing jobs; 30 s timeout in the tracker |
| #69 | SIGBUS kills the Snakemake controller after ~50 jobs update the DB | `PRAGMA mmap_size=0` |
| #75 | WAL's `-shm` file is always mmap'd → SIGBUS; DB corrupted four times | `journal_mode=DELETE` |
| #83 | DB locked → `qxtat check` says "running" forever | Fall back to status files |
| #107 | `exec` jobs missing from the DB since March; `cancel` can't see them | Run files (Phase 2) and per-user DB (Phase 4); to re-test on 3.9.0 |

Each fix treated a symptom. The underlying problem stays: **many processes on many nodes write one SQLite file on a parallel filesystem.** SQLite's own guidance is that this is unsafe. It relies on POSIX byte-range locks being correct and coherent across clients, and on Lustre/GPFS/NFS that depends on mount options (`flock`, `localflock`, `noflock`) and on the filesystem's lock manager under load.

Audit step 3b (keying on `QXUB_RUN_ID`) and moving more state into the DB would both add load. This doc proposes changing the shape of the writes first.

### What the Gadi tests showed (review of #132)

- **`/g/data` and `/scratch` are mounted with `localflock`.** POSIX locks are kept per node, so SQLite processes on different nodes don't see each other's locks at all. There is no mutual exclusion between a job on a compute node and `qxtat check` on a login node, or between two login nodes.
- **Corruption comes from cross-node access, not from contention as such.**
  - 12 writers on one node: 240/240 jobs recorded, `integrity_check` ok.
  - 6 writers on each of two nodes: the DB is corrupt within seconds.
- **The live shared DB is already corrupt.** The damage is in `job_resources`, its seventh corruption. Writes still succeed; full scans fail.

So **no amount of locking, retrying or timeout tuning makes a SQLite file written from more than one node safe here.** Busy timeouts, `BEGIN IMMEDIATE` and retries only coordinate processes on the same node. The fix is to make sure each SQLite file is written from one node only.

## Who writes the DB today

`get_db_path()` resolves to `QXUB_DB_PATH`, then `QXUB_SHARED_DB`, then `shared_db.path`, then `~/.config/qxub/qxub.db`. Submission exports `QXUB_SHARED_DB` to the job, so on the a56 setup every writer below hits `/g/data/a56/config/xdg/qxub/qxub.db`.

| Writer | Where it runs | When | Connections |
|---|---|---|---|
| `create_queue_entry` | login node | each submission | `init_db` (schema DDL, a write lock), then insert |
| job script `update_job_status` → `resource_tracker` | **compute node** | job start, job end | `ResourceTracker()` runs `_init_database` (DDL + migrations), then an update |
| job script → `mark_running` / `mark_complete` | **compute node** | job start, job end | another connection, another write |
| `qxtat check` / `list`, `cancel` | login node | Snakemake polls up to 10/s | mostly reads; some status updates |
| `scheduler` resource collection | login node | after `qx` jobs finish | update |

So **each job opens the shared DB about four times from its compute node, and each open begins with schema DDL that needs a write lock.** With 50–1000 jobs starting and ending around the same time, that is hundreds of near-simultaneous lock requests across nodes, on a filesystem whose lock coherence we don't control. Readers (`qxtat check`) queue behind them in DELETE mode, because a writer excludes all readers.

The compute-node writes also import qxub there, which needs `qxub_python` to be reachable from the node (fixed in #123), and the job's exit path depends on it.

## Two jobs, not one

The shared DB was introduced to collect resource-usage statistics, so that resource requests can be tuned. Over time it also took on live job tracking (`qxtat check`, `cancel`). These two jobs need very different things:

| | Live job state (`qxtat check`, `cancel`) | Resource statistics |
|---|---|---|
| Freshness | seconds | hours or days is fine |
| Writers | every job, on every node | one, occasionally |
| Scope | one run, or one workflow's runs | everything, across users |

Per-run files fit the first column. A database consolidated from those files now and then fits the second. Trying to make one SQLite file do both, with writers on every node, is what corrupts it.

## Principles

1. **Per-run files are the source of truth.** Each run has one small file on the shared filesystem. It is written atomically (temp file + `rename`) and needs no locks, like history (#130). `rename`, `mkdir` and `O_EXCL` creates go through Lustre's metadata server, so they are atomic across nodes, whatever `localflock` does to POSIX locks.
2. **Compute nodes never open SQLite**, and don't need qxub's Python at all.
3. **The shared DB has one writer at a time, and nobody updates it in place.** The consolidator builds a new copy privately and publishes it with `rename`. Readers open it read-only and immutable.
4. **Live state is read from the run files, not the DB.** A DB failure or lag never affects a status check.

## Proposal

### Phase 1: quick hardening (#132, merged)

- **One connection helper** for both `queue/db.py` and `ResourceTracker`, with the same pragmas and timeouts (previously 5 s vs 30 s).
- **Schema init only when `PRAGMA user_version` is behind**, so most connections never need a write lock to start.
- **`BEGIN IMMEDIATE` for writes**, so a writer takes its lock up front. A deferred transaction that upgrades from read to write can get `SQLITE_BUSY` immediately, which `busy_timeout` does not retry.
- **Retry with jittered backoff** on `database is locked`, `database is busy` and `disk I/O error`, then warn once and carry on.
- **Corruption: warn once, leave the file alone.** An earlier version moved a corrupt file aside and started afresh. Under cross-node writers that cascaded: each fresh DB was soon damaged too, and the live file lost nearly all its data.

Under `localflock` this only helps between processes on the same node, such as the per-user DB. **It does not prevent corruption of a DB written from several nodes.**

### Phase 2: one file per run (#133, merged). This is the fix for live state.

- **One file per run, keyed on `QXUB_RUN_ID`:**
  ```
  <spool_dir>/<run_id>.json
  ```
  The file is overwritten by `rename` as the run moves through its states, so each run costs **one inode**. Gadi has inode quotas on `/g/data` and `/scratch`.
  ```
  submitted  # written at submission: command, resources requested, tags, user, PBS job ID, joblog path
  running    # written by the job at start: + host, start time
  finished   # written by the job at end: + exit code, end time
  ```
- **Writing from the job:**
  - The job script merges its fields into the existing record with plain shell.
  - The submission side writes the record in a form the shell can extend without a JSON parser, e.g. one `key=value` per line. Alternatively, the job could write a small sibling file that is folded into the record later.
  - The spool directory is passed to the job explicitly (step 3b), rather than derived from `${out%/*}/status`.
- **No more DB writes from jobs:** the job scripts stop calling `update_job_status`/`mark_running`/`mark_complete`, so there is no `qxub_python` and no SQLite on the compute node.
- **Readers:**
  - `qxtat check` reads the run's file. It falls back to the existing status files and the DB, so in-flight jobs from before the upgrade still work.
  - `cancel` and `list` scan the spool for unfinished runs.
- **Compatibility:** the existing `final_exit_code_*` files stay for one release.

### Phase 3: occasional consolidation into a published snapshot (#135, merged)

- **The consolidator** reads finished run files and builds the statistics DB:
  - It applies them to a **private copy** of `qxub-stats.db` (copy, insert, `PRAGMA integrity_check`).
  - The snapshot lives in a directory of its own (`stats/` beside `shared_db.path`), created group-writable and setgid. Every publishing user needs write access to it, which the config directory beside it shouldn't grant.
  - Then it publishes the copy with `rename`. The DB is a few MB, so copying each round is cheap.
- **Why copy instead of updating in place:** under `localflock`, a *reader* on another node can't see the writer's lock. It would treat the writer's rollback journal as left over from a crash and "roll it back", which corrupts the file even with a single writer.
  - Readers open the published file with `?immutable=1`, so they never touch journals.
  - Readers always see a complete snapshot.
- **One consolidator at a time, without POSIX locks:**
  - The consolidator takes a lease with `mkdir .qxub-stats.db.lease` next to the snapshot, recording a token, host, pid and start time. It sits with the snapshot, not in each user's spool, because every user publishes to the same file.
  - A lease older than a few hours is stale. A process breaks it by renaming it out of the way, and only one rename wins.
  - A process that doesn't get the lease skips this round.
  - Just before publishing, a round checks the lease still holds its token. A round that outlived the stale timeout, and whose lease was taken over, doesn't publish, so it can't drop the other round's runs.
- **Resource usage is collected here:**
  - PBS appends resource usage to the joblog only *after* the job script exits, so the job can never record it itself.
  - The consolidator parses joblogs, as `backfill_resources` does today. It also calls `qstat -fx` while PBS still remembers the job.
  - Batching this is also kinder to PBS than one query per job.
- **Cleanup:**
  - Only run files already in a published snapshot are deleted, and only for finished runs older than a grace period (a few days, so `qxtat check` and `history` still see recent runs).
  - Consolidation is idempotent: a consolidator that dies partway leaves files for the next round.
  - Runs still in progress are never touched.
- **When it runs:**
  - *Opportunistically:* any `qxub` command checks a cheap condition, such as the last publish being more than 6 h ago or more than N finished runs waiting. If so, it starts a detached background consolidator, so the user's command doesn't wait.
  - *Explicitly:* `qxub db consolidate`, for a Snakemake `onsuccess:` hook or by hand.
  - A self-resubmitting `copyq` job could stand in for cron, if the opportunistic trigger proves too irregular.
- **Whose runs:**
  - To start, each user consolidates only their own runs into the shared snapshot. That avoids deleting other users' files, which would need group-writable, setgid spool directories.
  - The lease still keeps it to one writer at a time.
- **Lag:** statistics are as fresh as the last publish. That is fine for tuning allocations, and `qxub db status` says how old the snapshot is.

### Phase 4: retire the shared live DB (#136, merged)

- **Recover once:** `qxub db recover` reads the old shared DB from a private copy (`immutable=1`). If the copy is malformed, `sqlite3 .recover` rebuilds it first. Finished jobs, and jobs with usage, are added to the snapshot through the same lease and publish round as consolidation. Jobs the snapshot already has are kept, so a rerun adds nothing.
  - Rows the corruption garbled are left out: a job ID that isn't a PBS ID, a time that isn't ISO, or text in a numeric column. On Gadi, about 200 of 2600 rows were mixed up with `queue` rows or shifted columns.
  - A job still marked submitted or running that has an exit code lost its last status update. It is stored as completed or failed from the code.
  - It needs a `sqlite3` CLI of 3.29 or later for `.recover`. Gadi's `/usr/bin/sqlite3` is 3.26, so qxub's environments pin `sqlite>=3.29`, and `recover` prefers the `sqlite3` beside its own Python.
- **Stop writing it:** `get_db_path()` no longer returns `shared_db.path`/`QXUB_SHARED_DB`. Those only locate the old file and the snapshot's default `stats/` directory.
- **Read the snapshot:** `qxub resources` reports read the published snapshot, read-only.
- **Keep the per-user DB for now:** `~/.config/qxub/qxub.db` sees one user's login sessions, and Phase 1 covers same-node contention. Retire it later if the run files and the snapshot cover everything it does.

## Alternatives considered

- **WAL mode:** rejected in #75 (SIGBUS from the mmap'd `-shm`).
- **Lock files (`flock`) around every write:** same dependence on cross-node lock semantics that SQLite already has, which `localflock` doesn't provide.
- **A dot-lock (`O_EXCL` lock file) around every SQLite write:** works across nodes without POSIX locks, but needs stale-lock handling when a holder dies, and every compute node would still contend for it. It doesn't remove the coupling Phase 2 removes.
- **A DB server (Postgres, etc.):** none available to users on Gadi, and it would add a service to run.
- **One SQLite file per job, merged on read:** this is the run-file idea with a heavier format. A text file is simpler to write from shell.
- **Per-host SQLite caches of the run files** (an earlier draft of Phase 3): these keep each file to one node's writers, but they are another copy to keep in sync. Live state can be read straight from the run files, and statistics come from the snapshot.
- **Updating the shared stats DB in place from a single consolidator:** unsafe under `localflock`, because readers on other nodes can roll back the writer's journal (see Phase 3).

## How #107 fits

#107 shows submissions not reaching either DB since March. Phase 1's "warn once instead of swallowing" should surface why. Phase 2 gives `cancel` a correct view regardless, from the run files, and since Phase 4 each user's submissions go to their own DB, which only their login sessions write. It needs a re-test on 3.9.0 before closing.

## Open questions

1. ~~What does "flakes out" look like?~~ "disk I/O error" and "malformed" under cross-node writers, from the #132 review.
2. ~~Mount options~~: `/g/data` and `/scratch` are `localflock`. Is the same true on the compute nodes and `gadi-dm-*`? The cross-host results behave as if it is.
3. ~~Is the team-shared view needed?~~ Not for live state. `qxtat list`/`summary`/`cancel` show your own jobs (Phase 4); teammates' finished jobs are in the shared snapshot.
4. ~~Grace period before cleanup~~: 7 days after a run is consolidated (`GRACE_DAYS`).
5. ~~Cross-user consolidation~~: per user. Each user consolidates their own run files into the shared snapshot, so no spool needs to be group-writable; only the snapshot's `stats/` directory is.
6. **Consolidation trigger:** the opportunistic background round after a submission (at most every 6 h), plus `qxub db consolidate` by hand. Revisit a self-resubmitting `copyq` job only if snapshots go stale in practice.
7. **Retire the per-user DB?** It now only holds one user's login-side state. Run files could replace it for `qxtat` and `cancel`; not needed yet.
