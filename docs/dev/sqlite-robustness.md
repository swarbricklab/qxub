# Making the job database robust on shared filesystems

Status: **draft for review**. Phase 1 is in #132. Revised after the Gadi tests in the #132 review.

## Why

The job database (`qxub.db`) misbehaves on NCI Gadi when many jobs run at once. We have hit, in order:

| Issue | Symptom | What we did |
|---|---|---|
| #61 | `qxtat check` exits 1 under load, aborting Snakemake | Distinguish DB errors from missing jobs; 30 s timeout in the tracker |
| #69 | SIGBUS kills the Snakemake controller after ~50 jobs update the DB | `PRAGMA mmap_size=0` |
| #75 | WAL's `-shm` file is always mmap'd → SIGBUS; DB corrupted four times | `journal_mode=DELETE` |
| #83 | DB locked → `qxtat check` says "running" forever | Fall back to status files |
| #107 | `exec` jobs missing from the DB since March; `cancel` can't see them | Open |

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

## Principles

1. **Files are the source of truth; SQLite is an index.** A job's lifecycle is already written to per-job status files (`job_started_*`, `final_exit_code_*`) with plain shell. The DB should be rebuildable from those files plus submission records. Then losing or corrupting it costs a rebuild, not data.
2. **Compute nodes never open SQLite.** Jobs write small files only, the same way #130 made history safe: one file per record, temp file + `rename`, no locks.
3. **Fewer, shorter, idempotent DB writes, from login nodes only.** Status only moves forward (`initiated → dispatched → running → completed/failed`), so applying the same event twice is harmless.
4. **A DB failure never fails a job or a status check.** Readers fall back to the files (as #83 started).

## Proposal

### Phase 1: quick hardening (#132)

- **One connection helper** for both `queue/db.py` and `ResourceTracker`, with the same pragmas and timeouts (previously 5 s vs 30 s).
- **Schema init once per DB file, not per connection:** check `PRAGMA user_version` first and run DDL only when it is behind. Most connections then never need a write lock to start.
- **`BEGIN IMMEDIATE` for writes**, so a writer takes its lock up front. A deferred transaction that upgrades from read to write can get `SQLITE_BUSY` immediately, which `busy_timeout` does not retry.
- **Retry with jittered backoff** on `database is locked`, `database is busy` and `disk I/O error`, then warn once and carry on.
- **Corruption: warn once, leave the file alone.** An earlier version moved a corrupt file aside and started afresh. Under cross-node writers that cascaded: each fresh DB was soon damaged too, and the live file lost nearly all its data. Repair is a separate command, run while nothing writes (Phase 3).

Under `localflock` this only helps between processes on the same node. It reduces lock traffic and turns same-node contention into short waits, but **it does not prevent corruption.** Phase 2 does.

### Phase 2: jobs write events, not rows (with audit step 3b). This is the fix.

- **Spool layout:** each run gets a directory keyed on `QXUB_RUN_ID` under the status directory, which is passed to the job explicitly (step 3b):
  ```
  <status_dir>/<run_id>/started.json    # time, host, PBS job ID
  <status_dir>/<run_id>/finished.json   # exit code, time
  ```
  These are written atomically by shell (`printf … > .tmp && mv .tmp finished.json`). The existing `final_exit_code_*` files can stay for one release for compatibility.
- **No more DB writes from jobs:** the job scripts stop calling `update_job_status`, so there is no `qxub_python` on the compute node and no SQLite there.

### Phase 3: one DB per host, as a cache of the spool

Removing compute-node writes is not enough on its own: Gadi has several login nodes (and `gadi-dm-*` data movers), and `localflock` doesn't coordinate them either. So:

- **Submission records go to the spool too:** `<status_dir>/<run_id>/submitted.json`, written at submission. The spool on the shared filesystem becomes the complete record. It is written once per file, atomically, and needs no locks, like history (#130).
- **SQLite becomes a per-host cache:** e.g. `qxub-<hostname>.db`, written only by processes on that host, so `localflock` is enough. A per-host file can live on the shared filesystem or on local disk.
- **`sync_runs(run_ids)`:** reads the spool for runs the cache doesn't yet show as finished and applies newer state in one short `BEGIN IMMEDIATE` transaction. Called by `qxtat check`/`list`, `cancel`, `qxub history`, and at the end of a foreground `qx`/`qxi`.
- **`qxub db rebuild`:** recreates this host's cache from the spool. A corrupt or missing cache costs a rebuild, not data, and that rebuild is the repair command that Phase 1 leaves out.
- **Reads that need every run** (e.g. `qxtat list` across a workflow) sync first. Scanning the spool for unfinished runs keeps that cheap.

### Phase 4: retire the shared DB

- **Recovery:** once the spool and the per-host caches exist, the team-shared DB on `/g/data` has no writer. Read whatever it still has into the spool or the caches once (`sqlite3 .recover` copes with the existing corruption), then stop using `shared_db.path` and `QXUB_SHARED_DB`.
- **Cross-user `qxtat list` and `cancel`:** these read other users' spool directories, if team members want that.

### Phase 4: where the DB lives

Once compute nodes don't touch the DB, the remaining question is whether a team-shared DB on `/g/data` is worth it. The alternative is a per-user DB with team views built from the spool.

- **Per-user DB:** `~/.config/qxub/qxub.db` on NFS `/home` sees far fewer writers (one user's login sessions).
- **Shared DB:** what does it buy today? Mainly cross-user `qxtat list` and `cancel`. Both could instead read the spool directories.

## Alternatives considered

- **WAL mode:** rejected in #75 (SIGBUS from the mmap'd `-shm`).
- **Lock files (`flock`) around every write:** same dependence on cross-node lock semantics that SQLite already has, which `localflock` doesn't provide.
- **A dot-lock (`O_EXCL` lock file) around every SQLite write:** works across nodes without POSIX locks, but needs stale-lock handling when a holder dies, and every compute node would still contend for it. It doesn't remove the coupling Phase 2 removes.
- **A DB server (Postgres, etc.):** none available to users on Gadi, and it would add a service to run.
- **One SQLite file per job, merged on read:** this is the spool idea with a heavier format. JSON files are simpler to write from shell.

## How #107 fits

#107 shows submissions not reaching either DB since March. Phase 1's "warn once instead of swallowing" and the connection helper should surface why. Phase 3's rebuild gives `cancel` a correct view regardless.

## Open questions

1. ~~What does "flakes out" look like?~~ "disk I/O error" and "malformed" under cross-node writers, from the #132 review.
2. ~~Mount options~~: `/g/data` and `/scratch` are `localflock`. Is the same true on the compute nodes and `gadi-dm-*`? The cross-host results behave as if it is.
3. **Is the team-shared view needed?** Who uses cross-user `qxtat list`/`cancel`? This decides whether the spool must be readable across users.
4. **Spool retention:** how long to keep `<run_id>` directories? Tie it to the status-dir cleanup we already do?
5. **Per-host cache location:** shared filesystem (survives login-node reboots) or local disk (faster; rebuilt after a reboot)?
