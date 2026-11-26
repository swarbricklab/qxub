# PBS coupling audit

**Date:** 2026-11-16 · **Branch audited:** `docs/branching-model` · **Mode:** read-only, no code changed

qxub aims to be a universal executor: PBS Pro/Torque, Slurm, SGE, LSF, local, and eventually Kubernetes/cloud. This audit lists every place where PBS Pro or NCI Gadi assumptions are hard-coded. Four parallel audits covered the whole package; the headline claims were then spot-checked against the source.

| Detail report | Scope | High | Med | Low |
|---|---|---|---|---|
| [A: submission, templates, monitoring](A-submission-templates-monitoring.md) | `core/`, `jobscripts/`, `execution/` | 13 | 20 | 8 |
| [B: resources, platforms, config](B-resources-platforms-config.md) | `resources/`, `platforms/`, `config/`, platform YAML and design docs | 8 | 21 | 8 |
| [C: CLI layer](C-cli.md) | `*_cli.py`, `cli.py`, `notifications/` | 11 | 15 | 17 |
| [D: state, remote, packaging, tests](D-state-remote-packaging-tests.md) | `queue/`, `history/`, `remote/`, packaging, CI, `tests/` | 4 | 16 | 13 |
| **Total (with cross-report overlap)** | | **36** | **72** | **46** |

Several issues appear in more than one report. For example, `PBS_JOBID` in the templates shows up in A and D, and the status directory in A and C. The deduplicated picture is the 10 themes below.

---

## The short version

qxub has **no scheduler seam at all**. The design docs (`docs/dev/multi-platform-architecture.md`, `config_schema.md`) describe a `PlatformResources` ABC, a neutral `core_resources` model and a `scheduler:` platform key. None of these exist in code. Instead:

1. **Resources become PBS `-l` strings at CLI parse time**, before the platform is known. Every later stage edits them by string prefix (`startswith("ncpus=")`).
2. **The submit command is a hand-built `qsub -v … {options} template` f-string**, built on 4 separate code paths. The options string comes from `config/handler.py::build_qsub_options`.
3. **`Platform.type: pbs_pro` is only ever printed.** Nothing dispatches on it. Platform detection is `"gadi" in $HOSTNAME`.
4. **Job identity is the PBS job ID.** Jobscripts key status on `${PBS_JOBID:-unknown}`, the DB primary key is `pbs_job_id`, and the CLI appends `.gadi-pbs` to every ID.
5. **NCI project `a56` is hard-coded in shipped code.** It appears in jobscripts (`#PBS -P a56`, `storage=gdata/a56`, `/g/data/a56/software/qsub_tools`) and as the default `$PROJECT` for status and output dirs.

The good news: several important pieces are already neutral and are the right foundation. These are the execution-context layer (conda/module/singularity), the `qx-<uuid>` virtual ID created *before* submit, the status-file protocol, base64 command transport, `SubmitResult`, the workflow-style CLI flags (`--mem/--cpus/--runtime/--gpus/--disk`), the queue-limit model, and remote delegation to a remote `qxub exec`.

---

## Themes (deduplicated, in priority order)

### 1. No scheduler backend abstraction — **High** · blocks everything

- `core/scheduler.py` mixes TTY spinners, `qsub`/`qdel` wrappers, `qstat -fx` parsing, Gadi joblog parsing, monitoring and a background fork. `core/__init__.py` exports bare `qsub`, `qdel` and `get_job_resource_data`.
- `qsub()` refuses any command not starting with `"qsub"` ([scheduler.py:229](../../../qxub/core/scheduler.py)).
- The CLI imports these PBS helpers directly. `cancel_cli` and `interactive_cli` also shell out to `qstat`/`qsub` themselves.
- There are 4 drifting submission paths: `execution/context.py:260`, `core.py:287`, `executors.py:52` and `submit.py:275`. They differ in var expansion, user vars and notify handling.

### 2. Resources modelled as PBS strings — **High**

- `resources/mappers.py` (`get_pbs_resources`) converts `--mem/--cpus/--runtime/--disk/--gpus` into `mem=`/`ncpus=`/`walltime=`/`jobfs=`/`ngpus=` inside `exec_cli.py:777-805`, `interactive_cli.py:1245-1277` and `config_cli.py:1070-1112`.
- `config/handler.py` rewrites the strings by prefix (auto-queue, GPU→CPU adjustment, queue minimums), then renders `-N -q -P -l -o -M -m` in `build_qsub_options` (`:304-317`).
- **Shortcuts and aliases are persisted as PBS strings**, so users' saved config is not portable. This needs a config migration.
- There are three inconsistent memory parsers (bare number = MB vs bytes; no GiB/Gi distinction). Walltime parsing has no `D-HH:MM:SS`, and a parse failure **silently skips validation**.
- The model has no nodes/tasks concept. Multi-node jobs need a hand-written `-l select=…`.

### 3. Job identity tied to PBS IDs — **High**

- Jobscripts write `job_started_${PBS_JOBID:-unknown}` and similar files. Off PBS, every job collides on `unknown`, and the monitor (`scheduler.py:1003-1030`, no timeout) **waits forever**.
- The `project_jobs` primary key is `pbs_job_id`, with `INSERT OR IGNORE`. Equal Slurm/LSF IDs from two clusters, or reused PIDs, are silently dropped. No row records which backend/cluster owns the job.
- `status_cli.py:175,330` appends `.gadi-pbs` to every ID, so a Slurm job queried by the Snakemake status hook reports `running` forever.
- `history/manager.py:473` matches partial IDs by `split(".")[0]`. `monitor` validation rejects IDs without a `.` (Slurm `123_4`, LSF `123[4]`, k8s names).
- **The fix is already half-built:** key everything on the existing `qx-<uuid>` entry ID, exported to the job as `QXUB_RUN_ID`.

### 4. Hard-coded NCI Gadi / project `a56` values — **High**

- Every `.pbs` template has `#PBS -P a56`, `-l storage=gdata/a56`, `-W umask=0027` and `sys.path.insert(0, '/g/data/a56/software/qsub_tools')`. The last one fails silently everywhere else, so DB status is never updated.
- The status dir is `/scratch/{$PROJECT or "a56"}/$USER/qxub/status` (`scheduler.py:287-299`, `interactive_cli.py:239,355`). The output dir is `/scratch/$PROJECT/$USER/qt/…` (`config/handler.py:12-28`).
- The user-config template (`config/manager.py:533-603`) writes queue `normal`, project `a56`, `gdata/{project}` and binds to `/scratch,/g/data`. The `.qx/test.yaml` template uses `copyq`. `interactive_cli.py:486-492` searches `/g/data/a56/conda`.
- The default queue is the literal `"normal"` (handler, interactive, auto-select fallback). The default account is `$PROJECT`, an NCI login convention.

### 5. Exit status, completion and accounting from PBS/Gadi output — **High/Med**

- Detecting jobs killed by the scheduler (walltime/OOM/cancel) relies **only** on the Gadi epilogue `Exit Status:` block in the joblog (`scheduler.py:375-450`). That block is Gadi-specific, not generic PBS. The fallback globs every `*.log` and matches a fixed-width `"Job Id:             "` string.
- Raw PBS exit codes (negative, 256+signal) go straight into `sys.exit`, which truncates them mod 256 and loses the reason.
- Usage tracking (`resources/tracker.py`, `parser.py`) consumes `qstat -fx` field names and the Gadi epilogue. Its DB columns are `jobfs_*`, `pbs_job_id` and `cpu_percent`.
- Cost estimation is NCI Service Units (`cpus × hours × rate`), which is wrong even on Gadi (NCI bills memory-weighted). It also reads `su_rate` while the YAML sets `su_billing_rate`, so every queue costs the same.

### 6. Job-script templates — **High**

- There are 5 near-duplicate `.pbs` files (~250-320 lines each) that differ only in environment activation. Scheduler directives are baked into the body. *(Done: one `job.sh`; step 6b.)*
- All parameters travel in one comma-separated `qsub -v` string, with unquoted `cwd`/`out`/`err` paths. Singularity bind lists are also comma-separated. Slurm `--export`, LSF and k8s don't share this format.
- Quiet mode relies on PBS `-j oe` spooling into the joblog. The templates assume a `module` command and a module literally named `singularity` (many sites have `apptainer`).
- Notifications submit **another job from inside the job** (`qxub exec --internet … qxub notify send`), which assumes compute nodes can submit and an `internet` queue (Gadi `copyq`) exists.
- Only `jobscripts/*.pbs` is listed in `package_data`/`MANIFEST.in`.

### 7. Platform model and detection — **High/Med**

- `detect_platform()` matches only "gadi" in `$HOSTNAME` and ignores the YAML `host:` field. It otherwise falls back to *whichever platform loaded first* ([core.py:473-502](../../../qxub/platforms/core.py)).
- Auto queue selection searches **every loaded platform**, so it can pick a queue belonging to a different cluster. Lookups match queue names first-wins across platforms.
- Queues are the only scheduling target: there is no partition/QOS/account split (Slurm) and no selectors (k8s). `limits.max_walltime` in YAML is silently ignored, because only Gadi-style per-core-band rules exist.
- The repo's three platform YAMLs spell the scheduler key three ways (`type: pbs_pro`, `type: pbs`, `scheduler: pbs`). There is no schema validation. Two `PlatformLoader` classes exist, and the URL-capable one isn't wired in.

### 8. CLI surface and interactive mode — **High/Med**

- `qxub interactive` is `qsub -I … -- script` (`interactive_cli.py:727-761, 1442, 1493`). Its init script relies on `$PBS_JOBID` and runs `qstat -f` in an EXIT trap. *(Done in 5c: `backend.interactive_command()`, `backend.job_id_variable`, no `qstat -f` in the session.)* The script is written to login-node `$TMPDIR`, which compute nodes may not see. Interactive mode also has its own submission pipeline that bypasses `process_job_options`.
- `cancel` uses `qdel` directly. `cancel all` parses `qstat -u` text and filters on `Q/R/H`. `cancel name` parses `qstat -f`. `--wait` polls `qstat <id>`.
- `resources show` renders raw `qstat -fx` fields. `monitor` has no `--status-dir`.
- PBS-flavoured options and their proposed neutral equivalents are listed in [C-cli.md §"CLI surface"](C-cli.md). The main ones:

  | Current | Proposed |
  |---|---|
  | `-l` | `--native` (PBS-only `-l` alias) |
  | `-P/--project` | `--account` |
  | `-q` | `--queue` (`--partition` alias) |
  | `--jobfs` | `--disk` |
  | `--email-opts abe` | `--mail-on begin,end,fail` |
  | `--storage` | platform capability |
  | `cancel all --status Q/R/H` | `--state pending\|running\|held` |

  `--notify` help says "PBS email" but the flag actually controls Slack/Discord.

### 9. Shared-filesystem and compute-node assumptions — **Med** (blocks k8s/cloud only)

- Monitoring tails `out`/`err` and polls status files on the submit host, so it assumes a shared POSIX FS. Jobs update the submit host's SQLite DB by importing qxub on the compute node.
- This is fine for HPC, but it needs a `JobMonitor`/`StateStore` seam before k8s/cloud. It doesn't block Slurm/SGE/LSF.

### 10. Remote execution, tests and packaging — **High/Med**

- **Remote resubmit safety:** duplicate-submission protection after an SSH drop only recognises PBS-style `1234.server` IDs (`remote/platform_executor.py:22-25`). On Slurm/LSF a dropped connection **resubmits the job**. The regex also matches strings like `2024.log`. Fix it with a `QXUB-SUBMITTED entry_id=…` handshake plus a local `--entry-id` idempotency key. *(Handshake done in step 7b; the idempotency key isn't.)*
- The remote serialiser forwards raw `-l` strings. Remote init forces `eval "$(conda shell.bash hook)" &&`.
- **Tests:** there is no `conftest.py` and no fake scheduler. End-to-end tests assert `"qsub" in stdout`, and some tests `sys.path.insert` `/g/data/a56/...`. One test module imports a non-existent module. CI runs **no pytest at all**; its only integration checks SSH to Gadi.
- **Packaging/docs:** `check-nci-prerequisites.sh`, the Gadi-only CI setup script, the conda recipe ("PBS Pro tool", stale version, imports `qxub.platform`), and the agent skills (`skills/navigate-codebase` says "PBS interaction: always via subprocess.run(['qsub', …])").

---

## Proposed architecture

The four audits independently converged on the same seams:

```python
# qxub/backends/base.py
@dataclass
class ResourceSpec:          # typed, neutral; replaces "-l k=v" strings
    cpus, memory_bytes, walltime_s, gpus, gpu_type, local_disk_bytes,
    nodes, tasks_per_node, internet, extras: dict[backend, dict]   # escape hatch (Gadi storage=)

@dataclass
class JobSpec:               # replaces ctx.obj dict + qsub_options string
    run_id (qx-uuid), name, command, env_context, resources: ResourceSpec,
    queue/target, account, stdout, stderr, scheduler_log, mail, array, env, tags

class SchedulerBackend(Protocol):
    kind: Literal["pbs_pro","torque","slurm","sge","lsf","local","k8s"]
    capabilities: SchedulerCaps        # shared_fs, interactive, arrays, job_log, storage_directive, ...
    job_id_env_var: str                # PBS_JOBID / SLURM_JOB_ID / LSB_JOBID / JOB_ID
    # submission
    def render_directives(spec) -> str           # "#PBS"/"#SBATCH"/"" header for the neutral template body
    def parse_native_resources(args) -> ResourceSpec   # keeps -l back-compat inside PBS backend only
    def submit(spec, script) -> str              # parses its own job-id output
    def interactive_argv(spec, init_script) -> list[str]
    def classify_submit_error(rc, stderr) -> SubmitErrorKind
    def sanitize_job_name(name) -> str
    # identity / query / control
    def normalize_job_id(raw) / display_job_id(id)
    def status(job_id) -> JobState               # neutral enum
    def list_jobs(user, queue, states, name_regex) -> list[JobInfo]
    def cancel(job_id) -> bool
    # completion / accounting
    def final_state(job_id) -> ExitStatus | None # code, signal, reason ∈ {OK,FAILED,WALLTIME,OOM,CANCELLED,NODE_FAIL}
    def job_usage(job_id) -> JobUsage | None     # replaces qstat -fx / joblog parsing

class Platform:              # YAML; `scheduler:` dispatches to BACKENDS registry
    scheduler, detect (hostname regex / env var / binary), queues (with backend selectors),
    default_queue, account_env, paths.{status_dir, scratch, shared_tmp},
    accounting: AccountingModel, policy: QueuePolicy, container_runtime, notify_strategy,
    default_extra_directives   # Gadi: storage=gdata/a56, umask=0027 live here, not in code
```

DB and history changes:
- Add `backend`/`platform`/`cluster` columns.
- Make `entry_id` the primary key and rename `pbs_job_id` to `native_job_id`, with `UNIQUE(backend, cluster, native_job_id)`.
- Add a neutral `JobState` with a CHECK constraint.
- Rename `jobfs_*` to `local_disk_*`.
- Keep compatibility views for one release.

## Suggested order of work

1. **Fix the incidental bugs below.** Several of them break things today.
2. **Test harness first:**
   - add `tests/conftest.py` with a `FakeBackend`, plus PATH shims for `qsub`/`qstat`/`sbatch`/`squeue`;
   - add a hosted `pytest` CI job;
   - capture golden PBS output for current behaviour.

   *Status:* done, except `FakeBackend`, which moves to step 5 so it can implement the real `SchedulerBackend` interface. The fake PBS commands, `sandbox` fixture, golden files and the `Tests` CI workflow are described in `tests/README.md`. The harness found three bugs that have since been fixed: multi-path `--bind` was split by PBS, `cancel` exited 0 on a failed `qdel`, and `--array` was silently dropped (it is now rejected with an error until step 5 adds real array support). Two are still pinned by the golden files: `-P None` when no project is set, and no hugemem auto-selection unless `-q auto`.
3. *(Done. 3a: the `a56`/Gadi values come from `nci_gadi.yaml` and config. 3b: jobs get `QXUB_STATUS_DIR` (`site.default_status_dir()`, the directory `monitor` and `status` read) and name every status file after `QXUB_RUN_ID`, including the per-step `pre/main/post_exit_code` files, which jobs sharing an output directory used to overwrite. Each job also writes `native_<job ID>` with its run ID, so readers given a scheduler ID find its files without the DB (`core.scheduler.status_file_key`). Jobs from older qxub keep the old names and directory. The status directory is under the default project's scratch, so a job charged to another project gets the platform's `paths.scratch_storage` (Gadi: `scratch/<default project>`) added to `storage=`. If the job still can't create it, it falls back to `${out%/*}/status`, and `qxub exec` stops waiting once the joblog shows the job ended. The job DB no longer needs updates from jobs: they write run files (#133).)* **Remove the `a56`/Gadi constants from code** by moving them into `nci_gadi.yaml` (`paths`, `default_extra_directives`, `account_env`). Key status files and DB updates on `QXUB_RUN_ID`, not `PBS_JOBID`. This is a cheap, high-value change, and it also makes qxub usable at *other PBS sites*.
4. *(Done: `qxub/resources/spec.py` has `ResourceSpec` and the single `parse_size`/`parse_duration`. `JobSpec.resources` is a `ResourceSpec`, and `PbsBackend` alone writes `-l`. The golden files changed only in `-l` order (now `ncpus`, `ngpus`, `mem`, `walltime`, `jobfs`, then native requests such as `storage`) and in walltime padding (`1:00:00` is now `01:00:00`). `--cpus` and the other options now replace a matching `-l` request instead of adding a second one. Queues are still looked up by name, not through `(Platform, ResourceSpec)`.)* **Add `ResourceSpec` and `JobSpec`, with one memory parser and one duration parser.** Thread them through `process_job_options`. Show that the PBS output is byte-identical against the golden tests.
5. *(5a done in #65: `qxub/backends/` has `SchedulerBackend` and `PbsBackend` for queries and cancellation, picked by the platform's `type`. Every `qsub`/`qstat`/`qdel` call is in `backends/pbs.py`. 5b done: jobs are described by a `JobSpec` and submitted with `backend.submit()`. `qxub exec` and `submit_job` share one pipeline, `execution.submit.dispatch`, and `executors.py` / `submit_and_monitor_job` are gone. 5c done: `qxub interactive` describes its session as a `JobSpec` and starts it with `backend.interactive_command()` (PBS: `qsub -I … -- script`); the session script reads the job ID from `backend.job_id_variable` and no longer runs `qstat -f` as it ends, since qxub asks the backend for the finished job afterwards. `JobState` and a neutral `JobUsage` are still to do.)* **Add `SchedulerBackend` and `PbsBackend`** wrapping today's `core/scheduler.py`. Add `Platform.scheduler` and the registry. Collapse the 4 submission paths into 1. Route `cancel`/`status`/`monitor`/`resources show`/`interactive` through the backend.
6. *(6a done: the job's variables go in a settings file, `<status dir>/runs/<run ID>.env` (`export NAME='value'` lines, mode 0600), written by `queue.runs.attach_settings`. `JobSpec.env_file` names it and `PbsBackend` passes only `-v QXUB_JOB_FILE=<path>`. Each job script sources it first. If the file can't be written, the job is refused rather than putting the variables on `qsub -v` (other users can read them there, and its quoting is weaker); jobs from older qxub, with no `QXUB_JOB_FILE`, still read `-v`. A custom `--template` that doesn't source `$QXUB_JOB_FILE` is refused under PBS, with the line to add (`. "$QXUB_JOB_FILE"`), unless the user sets `QXUB_ALLOW_QSUB_V_SETTINGS=1` or `defaults.allow_qsub_v_settings: true`, which passes them with `-v` and a notice; Slurm and the local backend pass its variables in the submitting process's environment, so it isn't refused there. Consolidation removes the file with the run file. 6b done: the five templates are replaced by one script, `qxub/jobscripts/job.sh`, which sets up modules, conda and/or a Singularity image from the job's settings. `--pre` and `--post` run in the job's shell in every context, and each step runs as in a plain shell (no `set -e`). The `#PBS -S/-j` lines stay in its header for now: other schedulers ignore them, and a second backend can pass its own directives.)* **Neutral job-script body plus per-backend directive headers.** Pass parameters inside the script or a sidecar file instead of `qsub -v`.
7. *(7a done: schema version 3 adds a `backend` column to `queue` and `job_resources`, backfilled from the ID (`local-…` → `local`, anything else → `pbs`). Submissions record it there and in the run file, and `backends.backend_for_job(job_id, recorded)` picks a job's backend for cancelling (`cancel` uses the recorded one), `qdel`, job-detail lookups and consolidation's batched queries. 7b done: consolidation asks the backend each run file recorded. `backends.live_jobs()` lists the live jobs of the current backend and of every other backend available on this machine (local: its directory exists; PBS: `qstat` is on PATH), so `cancel name` and the stale-row check in `cancel status`/`filter` see local jobs under Gadi and PBS jobs under `local`; `cancel all` does only with `--all-backends`, so `cancel all -f` on the local platform doesn't kill your cluster jobs. A remote qxub asked to (`QXUB_REMOTE_HANDSHAKE=1`) prints `QXUB-SUBMITTED job_id=… backend=…` on stderr once it has submitted; the local side hides the line and won't retry a dropped connection after it (`qxub/remote/handshake.py`), falling back to the PBS-ID match for older remote versions. Jobs get their name as `QXUB_JOB_NAME`, which `job.sh` passes to notifications instead of `PBS_JOBNAME`. Deferred: re-keying `project_jobs`/`job_resources` on the run ID and renaming `queue.pbs_job_id`, which matters only if two backends on one machine hand out the same ID; and a `--entry-id` idempotency key for remote retries.)* 7. **DB/history migration** (`backend` column, `entry_id` key), and the remote `QXUB-SUBMITTED` handshake.
8. *(LocalBackend done: `qxub/backends/local.py`, `type: local` (`docs/platforms/local/local.yaml`). Jobs run as background processes recorded under `$XDG_STATE_HOME/qxub/local/` (the same under every platform); `qx`, `status`, `monitor`, `cancel` and `interactive` work with it, and the tests run real jobs through it. To get there, `job.sh` reads the job ID through `QXUB_JOB_ID_VAR` (the backend's `job_id_variable`) instead of `PBS_JOBID`, `status check` reads a finished run file (which also fixes PBS jobs killed at walltime showing as running without `--log-dir`), `status` appends `.gadi-pbs` only to bare numbers, and `monitor` accepts the backend's IDs. Backends describe jobs with `job_info()`, which returns a `JobInfo`: a neutral `JobState` (queued, running, held, finished, unknown), the exit status, `ended_by` when the scheduler can tell why it ended (PBS: `walltime`, from Exit_status -29; local: `cancelled`), timing, and a `JobUsage` in bytes and seconds. Consolidation, `resources show` and the end of an interactive session read it; `job_details()` remains as the scheduler's own record (`JobInfo.native`), which history stores, and `is_finished()` is gone. Still PBS-shaped: the job log's `Exit Status:` line is how `qx` notices a job that ends early (the local runner writes that line too). Why a job ended is kept: `qxub cancel` writes `ended_by=cancelled` to a running job's run file before cancelling it (and takes it back if the cancel fails), `qx` stops as soon as a cancelled job has recorded its exit and says why a job ended, a walltime kill read from the job log is added to the run file, and the `ended_by` columns (schema 4) reach `status show` and the statistics. The PBS backend reads Exit_status 271 (a `qdel`) as `cancelled` too, so a cancel from outside qxub, or by older qxub, also has a reason; `qx` notices it at its next joblog check, about a minute later. Not yet: a PBS job cancelled while still queued. SlurmBackend done: `qxub/backends/slurm.py`, `type: slurm` (example in `docs/platforms/slurm/slurm.yaml`). It submits with `sbatch --parsable`, passing the job's variables in sbatch's environment with `--export=ALL`; lists with `squeue`; describes jobs from `sacct`, taking `ended_by` from the state (`TIMEOUT`, `CANCELLED`, `OUT_OF_MEMORY`, the last a new `memory` reason) rather than the exit status; cancels with `scancel`; and runs interactive sessions with `srun --pty`. `qx` asks the job's backend how a job ended (`job_end()`; PBS and local read the job log's `Exit Status:` line, Slurm asks `sacct`), so a Slurm job cancelled while queued is noticed too. `cancel all --status Q/R/H` compares each backend's reading of its own states (`live_state()`). With no `--queue`, a job goes to the platform's `default_queue`, else `normal` if the platform has a queue by that name, else the scheduler's own default, instead of always `normal`. `status` adds `.gadi-pbs` to bare numbers only under PBS. Tests use fake Slurm commands (`tests/fakes/slurm`). Test clusters in containers (`ci/`, `tests/cluster`, the `Cluster tests` workflow) run real jobs on Slurm 24.11 and on plain OpenPBS 23.06. The OpenPBS one showed that `qx` waited forever for a job PBS killed or deleted at a PBS without Gadi's job-log epilogue, as at UNSW Katana; the PBS backend now asks `qstat -fx`, once, when the job log comes back without an `Exit Status:` line.)* 8. **Add `SlurmBackend`** as the second implementation to validate the seams, then `LocalBackend` (which also gives tests a real backend). k8s comes after the `JobMonitor`/`StateStore` seams.
9. **CLI renames with deprecated aliases** (`--account`, `--native`, `--disk`, `--mail-on`), a shortcut/alias config migration, and help-text and skills cleanup.

---

## Incidental bugs found (not PBS-specific)

Items marked ✔ were confirmed against the source during synthesis. The rest are as reported by the sub-audits.

| | Location | Bug |
|---|---|---|
| ✔ | `platforms/core.py:824` | `from .config import config_manager` resolves to the non-existent `qxub.platforms.config`, so **`qxub select-queue` fails on every call**. It should be `from ..config import …`. |
| ✔ | `jobscripts/qsing.pbs:11-12,181-203` | The template reads `$sif_file`/`$sing_options`, but Python passes `sif=`/`bind=`. **Singularity jobs run `singularity exec  <cmd>` with no image and no binds.** |
| ✔ | `jobscripts/qmod.pbs:240` | The unload block checks `$modules`, but the variable is `$mods`, so modules are never unloaded. |
| ✔ | `exec_cli.py:723` | `logger` is undefined, so the deprecated root-level config-key path raises `NameError`. |
| ✔ | `interactive_cli.py:1285 vs 1386` | `session_id` is overwritten before `_process_session_completion`, so post-session history and resource updates never find their files. |
| ✔ | `platforms/core.py:448` | `reload_platforms()` calls a non-existent `_load_all_platforms()`. |
| ✔ | `cancel_cli.py:500-511` | `cancel all` skips only 2 header lines and reads state from column 4. PBS Pro's `qstat -u` wide format has a longer header and `S` in a later column, so header lines may be treated as job IDs. **Verify on Gadi before relying on `cancel all`.** |
| | `history/base.py:180`, `history/manager.py:303,345` | `logger` is undefined. |
| | `config/base.py:300` | Imports the non-existent `qxub.config.platforms`. |
| | `config_cli.py:520-546` | `subcommand_def` is undefined for flat aliases. |
| | `alias_cli.py:181-233` | The deprecated `qxub alias` invokes `qxub --env` without `exec` and forwards non-existent options, so it is non-functional. |
| | `interactive_cli.py:217-219` | `_encode_command` uses `base64` without importing it (dead code). |
| | `core/scheduler.py:959` | `return thread` after `os._exit(0)` is unreachable and refers to an undefined name. |
| | `resources/tracker.py:614-637` | `get_resource_trends` calls `dict(row)` without setting `row_factory`. |
| | `notifications/slack.py:194`, `discord.py:112` | The footer suggests a non-existent `qxub history show {job_id} --output` command. |
| | `tests/test_remote_execution_flow.py:14` | Imports the non-existent `qxub.execution_context`. |
| | `recipe/meta.yaml:40` | The test imports the non-existent `qxub.platform`. The recipe version (3.3.1) is stale compared with `setup.py` (3.6.1). |
| | `setup.py:35` | `python_requires>=3.6`, but the code uses `tuple[...]` (3.9+) and `Path \| None` (3.10+). |
| | `execution/mode.py:55,66` | Docstrings use `remote.hostname`, but the validator requires `remote.host`. |
