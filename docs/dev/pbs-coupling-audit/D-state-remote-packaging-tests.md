# Audit D: PBS coupling in job state DB, history, remote execution, packaging and tests

Scope: `qxub/queue/{db.py,__init__.py}`, `qxub/history/*.py`, `qxub/remote/*.py`, packaging/infra (`setup.py`, `pyproject.toml`, `MANIFEST.in`, `Dockerfile`, `recipe/`, `scripts/*.sh`, `skills/`, `environment*.yaml`, `.github/workflows/`), and `tests/`.
This was a read-only audit; no project files were changed. All paths are relative to `/Users/johree/projects/qxub/qxub`. I checked the line numbers against the current files.

Severity: **High** = blocks any non-PBS backend or causes unsafe behaviour on one. **Med** = wrong behaviour, or needs a refactor before a second backend is clean. **Low** = naming, docs or cosmetic issues.

**Totals: 33 findings. High 4, Med 16, Low 13.**

---

## A. Job state database (`qxub/queue/db.py`, `qxub/queue/__init__.py`)

### A1. No backend, scheduler or cluster column. Severity: High
- **Where:** `qxub/queue/db.py:135-187` (`queue` DDL), `196-208` (`project_jobs` DDL), `247-271` (`_QUEUE_NEW_COLUMNS`)
- **Category:** DB schema
- **Problem:** No row records which scheduler or platform issued the job ID. That means status, cancel and resource lookups cannot be sent to the right backend. A single shared DB used by PBS, Slurm and local runs (or by two clusters) cannot be told apart. There is a `hostname` column (163), but it holds the execution host, not the backend identity.
- **Fix:** Add `backend TEXT` (e.g. `pbs`, `slurm`, `lsf`, `sge`, `local`, `k8s`), `platform TEXT` and `cluster TEXT` columns. Add them to `_QUEUE_NEW_COLUMNS` and to a matching `project_jobs` migration. Backfill existing rows with `backend='pbs'`, and make `create_queue_entry()` require `backend`.

### A2. `project_jobs` primary key is the native scheduler ID, and inserts fail silently on collision. Severity: High
- **Where:** `qxub/queue/db.py:197` (`pbs_job_id TEXT PRIMARY KEY`), `453` and `544` (`INSERT OR IGNORE`), and every lookup `WHERE pbs_job_id=?` (`478`, `605`, `612`, `639`, `647`, `702`)
- **Category:** Job-ID handling
- **Problem:** PBS IDs carry a server suffix, so they are mostly globally unique. Slurm and LSF IDs are bare integers per cluster, local PIDs get reused, and K8s names are per namespace. A second job with the same native ID is silently dropped by `INSERT OR IGNORE`, and `mark_running`/`mark_complete` would update every row with that ID. Cross-ref: `qxub/resources/tracker.py:428,701,748-758,814,1004` also writes `queue` keyed on `pbs_job_id`.
- **Fix:** Make `entry_id` (the virtual ID) the primary key of `project_jobs` too. Store the native ID as a non-unique `native_job_id` with a UNIQUE index on `(backend, cluster, native_job_id)`, and route every update through `entry_id`.

### A3. Job-side lifecycle updates are keyed on the native PBS ID, not the qxub virtual ID. Severity: Med
- **Where:** `qxub/queue/db.py:594-651` (`mark_running(pbs_job_id)`, `mark_complete(pbs_job_id, ...)`); called from the job scripts with `${PBS_JOBID:-unknown}` (cross-ref `qxub/jobscripts/qconda.pbs:78-80`, `qconda_mod.pbs:79-81`, `qdefault.pbs:86`)
- **Category:** Job-ID handling
- **Problem:** On Slurm the variable is `$SLURM_JOB_ID`, on LSF `$LSB_JOBID`, on SGE `$JOB_ID`, and local or K8s runs have none. On a non-PBS backend every job reports as `unknown`, which collides on the `WHERE` clause.
- **Fix:** Export `QXUB_ENTRY_ID=qx-...` into the job environment at submit time (it already exists before submit; see `create_queue_entry`, 371-465), and key `mark_running`/`mark_complete` on `entry_id`. Keep a native-ID variant only as a fallback.

### A4. The status model is free text with no enum, and terminal states are binary. Severity: Med
- **Where:** `qxub/queue/db.py:148` (`'initiated'` default), `203` (`'submitted'`), `410`, `455`, `546`, `605`, `611`, `631` (`"completed" if exit_code == 0 else "failed"`)
- **Category:** State enum
- **Problem:** The DB itself stores neutral words, not the PBS letters Q/R/H/E/F (the letters live only in `status_cli.py:351-358` and `monitor_cli.py:134-236`, which other agents cover). But there is no enum or CHECK constraint. Two vocabularies are mixed (`queue` uses `initiated/dispatched/...`, `project_jobs` uses `submitted/...`), and there are no states for cancelled, timeout, out-of-memory, preempted, node-fail or held. Slurm and K8s report those distinctly, and PBS reports them through exit codes above 256.
- **Fix:** Add a neutral `JobState` enum (`PENDING, SUBMITTED, QUEUED, HELD, RUNNING, COMPLETED, FAILED, CANCELLED, TIMEOUT, OOM, PREEMPTED, UNKNOWN`), let each backend map its native states and exit codes onto it, add a CHECK constraint, and make `mark_complete` take a `JobState` rather than inferring one from `exit_code`.

### A5. Resource and metadata columns use PBS/Gadi vocabulary, and there are no GPU columns. Severity: Med
- **Where:** `qxub/queue/db.py:159` and `396` (`joblog_path`, "PBS joblog"), `161-162` and `397` (`queue_name`, "PBS queue name"), `169`, `175` and `182` (`jobfs_*`, a Gadi-specific node-local disk), also in `_QUEUE_NEW_COLUMNS` `254`, `259`, `263`, `268`
- **Category:** DB schema
- **Problem:** Slurm has partitions, not queues. Node-local scratch is `--tmp`, TMPDIR or an ephemeral volume, not jobfs. A separate "joblog" file is an NCI convention. GPU requests and usage, which every modern backend has, are not modelled at all.
- **Fix:** Rename to neutral columns (`queue_or_partition`, `local_disk_requested_mb/used_mb/efficiency`, `scheduler_log_path`), add `gpus_requested`, `gpus_used` and `gpu_type`, and add `backend_extra TEXT` (JSON) for backend-specific fields. Expose the old names as SQL views during migration.

### A6. Cross-node state assumes a shared POSIX filesystem and Python on compute nodes. Severity: Med
- **Where:** `qxub/queue/db.py:50` and `58-60` (`QXUB_SHARED_DB` set by job scripts), `9-10` (docstring), `93-96`
- **Category:** Shared-filesystem assumption
- **Problem:** Running jobs write their own status by opening the submit host's SQLite file and calling `python3 -c "from qxub.queue.db import mark_running"` (cross-ref the `.pbs` templates). This needs the same filesystem path and a qxub-importable Python on every compute node. That is false for K8s and cloud batch, and often false on Slurm clusters with node-local `$HOME` or containers.
- **Fix:** Put a `StateStore` interface in front of the DB. Make job-side reporting optional and backend-driven: either the scheduler is polled (`sacct`, `qstat`, the K8s API), or the job writes a status sentinel to a staged or log location that the submit side ingests.

### A7. Docstrings and comments say PBS throughout, and the WAL docs disagree with the code. Severity: Low
- **Where:** `qxub/queue/db.py:9` (says WAL) versus `94` (`PRAGMA journal_mode=DELETE`); `120-122`, `131`, `161`, `192`, `385-397`, `449`, `469`, `493-503`, `535`, `595`, `625`, `694`; `qxub/queue/__init__.py:6-7`, `16-19`
- **Category:** Naming and docs
- **Problem:** "PBS ID", "after qsub()" and "PBS-level job tracking" are written into the API docs. The WAL claim is also wrong; DELETE mode is probably deliberate, for Lustre/NFS.
- **Fix:** Use neutral wording ("native scheduler job ID", "after backend submit") and correct the journal-mode note.

### A8. The user is read from `$USER`. Severity: Low
- **Where:** `qxub/queue/db.py:405`
- **Category:** Environment assumption
- **Problem:** `USER` is often unset in containers, K8s pods and cron/systemd jobs, so rows get `unknown`.
- **Fix:** Use `getpass.getuser()` and fall back to `os.getuid()`.

**Already well abstracted:** Virtual IDs (`new_virtual_id`/`is_virtual_id`, `356-363`) are created *before* submit, and the native ID is nullable (`149`, `373`). That is the right seam for an opaque `JobRef`. Status words are neutral rather than PBS letters, and `db_path` injection makes the layer easy to test.

---

## B. Job history (`qxub/history/*.py`)

### B1. Job-ID lookup assumes the PBS `<number>.<server>` format. Severity: Med
- **Where:** `qxub/history/manager.py:449-479`; the matching logic is at `473` (`stored_id.split(".")[0] == job_id`)
- **Category:** Job-ID handling
- **Problem:** Partial matching on the text before the first dot only works for PBS IDs like `153392916.gadi-pbs`. Slurm array IDs (`12345_7`), LSF array IDs (`123[4]`) and K8s names are not handled. Lookup also does not filter by backend, so equal native IDs from different backends collide.
- **Fix:** Store `entry_id` and `backend` on each execution record, look records up by `entry_id` first, and delegate partial-ID normalisation to `backend.normalize_job_id()`.

### B2. Execution records have no backend or platform field and store a PBS-shaped resource payload. Severity: Med
- **Where:** `qxub/history/manager.py:243-270` (record body), `252-253` (`job_id`, `joblog`), `256-258` and `324` (`exit_status`), `264-270` and `328-334` (`resources_requested/used/efficiency/execution/timing` from the PBS tracker)
- **Category:** History schema
- **Problem:** History records cannot say which backend ran the job. `joblog` is a PBS/NCI artefact, and the resource keys mirror the PBS joblog/`qstat -f` parser output. The records are not linked to the DB's virtual `entry_id`.
- **Fix:** Add `backend`, `platform` and `entry_id`. Define a neutral `ResourceUsage` dataclass, shared with `queue` (see A5), and rename `joblog` to `scheduler_log`.

### B3. Recipe extraction re-parses `sys.argv`, misreads PBS short flags, and mixes scheduler options into the recipe hash. Severity: Med
- **Where:** `qxub/history/manager.py:59-153` (`_extract_recipe_from_command`), especially `102-119`; the hash is computed at `47-57`
- **Category:** History / replay format
- **Problem:** Given `qxub exec -l walltime=1:00:00 --env x -- cmd`, the parser skips `-l` (120-121) and then treats `walltime=1:00:00` as the start of the target command (116-119). The executor is never detected and the hash is wrong. Long options such as `--queue`, `--project` and `--resources` go into `executor_params` (102-110), so the "what" hash includes scheduler-specific "how" options, and recipes are not portable across backends.
- **Fix:** Build recipes from `ctx.params` and the parsed `ExecutionContext`, not from argv. Keep `executor` (the environment: conda, module, container) separate from `resources` (neutral mem/cpus/runtime/gpus/disk) and `scheduler_passthrough` (raw `-l`/`--slurm-arg`), and hash only executor plus target.

### B4. Raw argv and raw click params are the only replay information. Severity: Low
- **Where:** `qxub/history/manager.py:246` (`command_line`), `250` (`main_options: ctx.params`); `qxub/history/base.py:73`, `81`
- **Category:** History / replay format
- **Problem:** `main_options` contains PBS `-l` tuples, queue and project as typed. Replaying on a different platform would replay PBS syntax.
- **Fix:** Persist the normalised request (see B3), alongside the raw argv kept for display only.

### B5. Docs and examples use Gadi and PBS. Severity: Low
- **Where:** `qxub/history/manager.py:15` ("PBS log"), `450-456`, `494-510` (`153392916.gadi-pbs`, `/scratch/a56/...`); `qxub/history/__init__.py:33-36` (`walltime`/`ncpus` example)
- **Category:** Naming and docs
- **Fix:** Use neutral examples, such as several ID formats and `/path/to/...`.

### B6. `logger` is never defined in the history modules (not a PBS issue). Severity: Low
- **Where:** `qxub/history/base.py:180`; `qxub/history/manager.py:303`, `345`
- **Category:** Bug found while auditing
- **Problem:** The `except` blocks call `logger.debug` without a module-level logger. This raises `NameError` out of the "never break the main command" path. The current callers wrap the call in `try`, which hides it.
- **Fix:** Add `logger = logging.getLogger(__name__)`.

**Already well abstracted:** The executor-type whitelist (conda/module/sing/default, `manager.py:177-185`) describes execution environments, not schedulers. Storage is local YAML under `~/.config/qxub`, which works for any backend.

---

## C. Remote execution (`qxub/remote/*.py`)

### C1. Duplicate-submission protection only recognises PBS-format job IDs. Severity: High
- **Where:** `qxub/remote/platform_executor.py:22-25` (`_PBS_JOB_ID_RE = r"\b(\d{4,})\.[A-Za-z][\w-]*\b"`), `249-251`, `266-278` (hint: `ssh host qstat <id>`), `304`
- **Category:** Remote / job-ID handling
- **Problem:** If SSH drops after submission, the retry is suppressed only when stdout contained a PBS-style `1234.server` ID. Slurm prints `Submitted batch job 12345` (or a bare `12345` with `--parsable`), LSF prints `Job <123> is submitted`, and local or K8s runs print other formats. On those backends a dropped connection causes a **duplicate job submission**. The regex also produces false positives on strings like `2024.log` or `1234.txt`, which suppresses legitimate retries.
- **Fix:** Have the remote qxub print a structured, backend-neutral line immediately after submit, for example `QXUB-SUBMITTED entry_id=qx-… backend=slurm native_id=12345`, and match only that line. Better still, generate the `entry_id` locally and pass it as `--entry-id` (an idempotency key) so the remote can refuse a duplicate. Have the backend supply the "check status" hint command.

### C2. The serialiser forwards PBS resource syntax and qsub-style short flags. Severity: Med
- **Where:** `qxub/remote/command_builder.py:63-64` and `106-123` (`_add_pbs_options`: raw `-l walltime=…`), `157-165` (`-N`, `-P`), `184-185` (`--email-opts`, which uses PBS `-m abe` semantics), `22` (docstring "PBS options")
- **Category:** Remote command format
- **Problem:** The design is sound (see "Already well abstracted" below), but raw `-l` strings are passed through unchanged. If the remote platform is Slurm, PBS resource syntax reaches a backend that cannot interpret it.
- **Fix:** Forward only neutral options (`--mem`, `--cpus`, `--runtime`, `--gpus`, `--disk`, `--name`, `--account`). Send raw `-l` as an explicit `--scheduler-arg` passthrough, tagged with its source backend, and have the remote qxub reject it or warn when its backend differs.

### C3. A conda initialisation is forced on the remote host. Severity: Med
- **Where:** `qxub/remote/platform_executor.py:382-393` (default `eval "$(conda shell.bash hook)"`, chained with `&&` at `405`)
- **Category:** Remote environment assumption
- **Problem:** If the remote host has no conda (pip/venv installs, `module load python`, a container, or a cloud VM), the `&&` chain fails before qxub runs. The code also assumes the remote login shell is bash.
- **Fix:** Make `remote.env_init` an optional list of commands with no default. Add `remote.qxub_path` so qxub can be called by absolute path, and document conda as one option among several.

### C4. The `{project}`/`{{project}}` variables are tied to NCI's `$PROJECT`. Severity: Med
- **Where:** `qxub/remote/platform_executor.py:127-165` (docstring at `133`: "Gadi execution environment"; `155`, `163-164`)
- **Category:** NCI-specific convention
- **Problem:** `$PROJECT` is an NCI login variable. Slurm sites use account names (`$SLURM_ACCOUNT`, or nothing), and cloud backends use neither.
- **Fix:** Support generic substitution (`{env:VAR}` locally, `{{env:VAR}}` remotely), with platform-defined variables such as `account`, and keep `{project}` as a deprecated alias.

### C5. The remote working directory is unquoted, and there is no staging (shared-filesystem assumption). Severity: Med
- **Where:** `qxub/remote/platform_executor.py:380` (`cd {self.working_dir}` without quoting); the class as a whole (`85-125`)
- **Category:** Filesystem assumption
- **Problem:** The remote run assumes inputs already exist on the remote filesystem and that the working directory exists. That holds for an HPC login node but not for cloud or K8s, where inputs need staging. Leaving the path unquoted also breaks paths with spaces and allows shell injection from the config value.
- **Fix:** Apply `shlex.quote` to the path but keep `$VAR` expansion (e.g. `cd "…"`). Add an optional `StagingStrategy` (none / rsync / object store) hook for non-shared-filesystem backends.

### C6. Comments and messages are written for PBS/NCI. Severity: Low
- **Where:** `qxub/remote/platform_executor.py:22-25`, `30-31`, `249`, `272-275`, `304`, `332-335` (comment "CI runners against NCI"); `qxub/remote/command_builder.py:22`, `63`, `105-107`, `271-272`
- **Category:** Naming and docs
- **Fix:** Use neutral wording ("scheduler job ID").

### C7. `parse_options_from_ctx` is unused and parses qsub options. Severity: Low
- **Where:** `qxub/remote/command_builder.py:224-285`
- **Category:** Dead code
- **Problem:** Nothing in the repository calls it. It walks a `qsub_options` list of `-l` strings.
- **Fix:** Delete it.

**Already well abstracted:** Remote execution delegates to a remote `qxub exec` (`command_builder.py:40`, `platform_executor.py:395-402`) rather than building `qsub` locally. The remote side resolves its own platform definition through `QXUB_PLATFORM`, so the transport (SSH options, retries, idle watchdog, `~/.ssh/config` delegation) has nothing to do with the scheduler. Only the job-ID detection (C1) and the option vocabulary (C2) leak PBS details.

---

## D. Packaging, infrastructure, CI and skills

### D1. Package data only ships `*.pbs` templates. Severity: Med
- **Where:** `setup.py:27-29` (`"qxub": ["jobscripts/*.pbs"]`), `MANIFEST.in:2`
- **Category:** Packaging
- **Problem:** Slurm, LSF or SGE templates (`.sbatch`/`.sh`) added under `jobscripts/` would be left out of wheels and sdists.
- **Fix:** Organise templates as `jobscripts/<backend>/*` and use `package_data={"qxub": ["jobscripts/**/*"]}` with `recursive-include qxub/jobscripts *`. Better still, move to `importlib.resources`.

### D2. The conda recipe calls qxub a PBS Pro tool and is stale. Severity: Low
- **Where:** `recipe/meta.yaml:55-64` (summary "via PBS Pro", description "PBS job submission wrapper"), `2` (version 3.3.1 vs `setup.py:18` 3.6.1), `15-19` (no `qxi` entry point), `40` (the test imports `qxub.platform`, which does not exist; only `qxub/platforms` does), `52` and `65-66` (`qsub_tools` URLs); `recipe/README.md` (`/g/data/a56/...` build path)
- **Category:** Packaging
- **Fix:** Use neutral wording ("HPC/batch job executor: PBS, Slurm, …"), generate the version from a single source, and fix the test imports and entry points.

### D3. `python_requires` does not match the syntax the code uses. Severity: Low
- **Where:** `setup.py:35` and `recipe/meta.yaml:23,27` (`>=3.6`), while `qxub/remote/platform_executor.py:425` uses `tuple[int, bool, str]` (3.9+) and `qxub/replay_cli.py:26` uses `Path | None` (3.10+); `environment.dev.yaml:7` pins 3.9; `environment.yaml` is named `qt` and is missing `pyyaml` and `requests`
- **Category:** Packaging hygiene
- **Fix:** Require `>=3.10` and align the environment files with it.

### D4. `scripts/check-nci-prerequisites.sh` is entirely specific to Gadi. Severity: Med
- **Where:** `scripts/check-nci-prerequisites.sh:2-7` (title), `10-18` (conda required), `53-58` (`/g/data/a56/...nci_gadi.yaml`), `73-89` (`/scratch/$PROJECT/$USER`), `64` (GNU-only `stat -c`)
- **Category:** Infrastructure
- **Fix:** Replace it with a `qxub doctor` subcommand that runs generic checks (Python, config, DB writable, SSH) plus backend-supplied checks (`backend.preflight()`: `qstat -B` / `sinfo` / `kubectl auth can-i`).

### D5. `scripts/setup-secure-ci.sh` is specific to Gadi. Severity: Low
- **Where:** `scripts/setup-secure-ci.sh:7` (`nci-automation`), `23-31` (`gadi_private_key`, `ssh-keyscan gadi.nci.org.au`), `37-63`, `79-80` (`NCI_USERNAME`/`NCI_PROJECT`), `72` and `86` (reference `test-remote-execution-secure.yml`, which does not exist)
- **Category:** CI infrastructure
- **Fix:** Parameterise the host and secret names, or move the script to a site-specific repository or directory.

### D6. CI runs no tests, and its only integration checks target Gadi and PBS. Severity: Med
- **Where:** `.github/workflows/pylint.yml`, `formatting.yml` and `release.yml` (none runs `pytest`); `.github/workflows/runner-validation.yml:32` (`Host gadi`), `71-73` (`ssh gadi … qstat -B`); `runner-connection-test.yml:29` (note: "require SSH to NCI Gadi")
- **Category:** CI
- **Problem:** There is no hosted CI that would catch a regression in a second backend. The self-hosted runners only validate access to Gadi.
- **Fix:** Add a `pytest` job on `ubuntu-latest` with a matrix over `backend: [fake, local]`, and later `slurm` through a containerised Slurm image (e.g. `giovtorres/slurm-docker-cluster`) and `pbs` through an OpenPBS container. Keep the Gadi runner as an optional site-specific job.

### D7. The Dockerfile is generic, but it cannot submit to any scheduler and has no test backend. Severity: Low
- **Where:** `Dockerfile:1-44`
- **Category:** Packaging
- **Problem:** The image contains no fake PBS and no PBS assumptions, which is good. It also contains no scheduler client, so it is only useful for a local or K8s backend.
- **Fix:** Document the image as the runtime for `local` and `k8s` backends. Optionally add a `test` stage that puts the fake-scheduler shims from E1 on `PATH`.

### D8. The agent skills describe qxub as PBS-only and teach PBS patterns. Severity: Low
- **Where:** `skills/hpc-agent/SKILL.md:2`, `7-8` ("submits PBS jobs", "Unlike `qsub`"), `35` (`--resources ngpus=1` as the GPU example, even though the neutral `--gpus` exists), `75` ("PBS log"); `skills/navigate-codebase/SKILL.md:45`, `52` ("PBS interaction: Always via `subprocess.run(['qsub', …])`"), `54` (`jobscripts/*.pbs`); `skills/navigate-codebase/modules.json` (8× "PBS", 4× "pbs", 2× "qsub", 2× "QsubError", 2× "gadi")
- **Category:** Docs / agent guidance
- **Problem:** These files steer contributors and agents toward PBS-only code paths.
- **Fix:** Use the neutral `--gpus` in the examples. After the backend seam lands, change navigate-codebase to say "scheduler interaction goes through `Backend.submit()`".

**Already well abstracted:** The Dockerfile (no fake PBS baked in), `pyproject.toml` (black config only), the `install_requires` list (pure Python, no scheduler dependencies), and the CLI entry points.

---

## E. Tests (`tests/`)

### E1. There is no fake-scheduler fixture or backend-parametrised harness, and tests assert on `qsub`. Severity: High
- **Where:** No `conftest.py` exists anywhere in the repository. Nothing mocks `qsub`/`qstat`/`qdel` (grep finds no `patch(...qsub...)`). `tests/end_to_end/test_complete_workflows.py:34`, `55`, `64`, `132`, `155`, `179` assert `"qsub" in result.stdout`. `tests/migration/test_pre_migration_baseline.py:59`, `75` and `tests/migration/test_import_compatibility.py:41-43` treat `qsub`/`qdel` imports as public API.
- **Category:** Test harness
- **Problem:** A second backend would make the end-to-end tests fail. The tests also could not run a non-PBS backend at all, because nothing lets them choose one.
- **Fix:** Add a `tests/conftest.py` with:
  - a `fake_backend` fixture: an in-memory `Backend` that records submissions and returns scripted states;
  - a `fake_scheduler_bin` fixture that puts shim `qsub`/`qstat`/`sbatch`/`squeue` executables on `PATH`, logging JSON so the CLI can be tested end to end;
  - `@pytest.mark.parametrize("backend", ["pbs", "slurm", "local"])`.

  Assert on a neutral dry-run marker (e.g. `"Would submit via <backend>"`) instead of the string `qsub`.

### E2. Tests hard-code Gadi paths and IDs, and one test module is broken. Severity: Med
- **Where:**
  - `tests/test_remote_execution_flow.py:8` and `tests/test_command_serialization.py:7` (`sys.path.insert(0, "/g/data/a56/software/qsub_tools")`).
  - `tests/test_remote_execution_flow.py:14` imports the non-existent `qxub.execution_context`, so the module fails at import.
  - `test_remote_execution_flow.py:47-51`, `64-65` (`gadi`, `/scratch/a56`, `nci_gadi`).
  - `tests/test_queue_db_schema.py:90`, `92`, `102`, `120`, `140-143` (`*.gadi-pbs` IDs).
  - `tests/test_monitor_features.sh:9-30` (`--suffix .gadi-pbs`).
  - `tests/configs/test_remote.yaml:7`, `13-16`, `22` (`file:///g/data/a56/.../nci_gadi.yaml`).
  - `tests/test_remote_execution_integration.sh:21-36`.
  - `tests/test_conda_dry.sh:192` (`/scratch/a56/test`).
- **Category:** Test fixtures
- **Fix:** Drop the `sys.path` hacks; tests should use the installed package. Use `tmp_path` and bundled fixture platform files. Parametrise job-ID formats (`"1234.pbs01"`, `"12345"`, `"12345_3"`, `"123[4]"`, `"qxub-job-abc"`) in the DB and history tests.

### E3. Test fixtures disagree on the platform's scheduler key. Severity: Med
- **Where:** `tests/test_platforms/localhost.yaml:6` (`scheduler: "pbs"`), `tests/test_platforms/test_platform.yaml:2` (`type: pbs`), `tests/test_platform.py:34` (`"type": "pbs"`), `tests/test_enforce_queue_minimums.py:22` (`"type": "pbs_pro"`)
- **Category:** Platform schema / test fixtures
- **Problem:** Three spellings across two keys means the backend selector is not a validated part of the platform schema. Nothing tests how the backend is chosen from it.
- **Fix:** Define a single `scheduler:` field with an enum and validate it. Add a test that a platform with `scheduler: slurm` resolves to the Slurm backend, or fails clearly if that backend is not implemented.

### E4. Most of the suite is shell scripts that print results and need a live PBS/Gadi install. Severity: Low
- **Where:** `tests/README.md` (describes shell scripts as the main suite); `tests/test_queue_cpu_adjustment.sh:53-59` ("PBS qsub will reject"); `test_conda_integration.sh`, `test_unified_cli.sh`, and similar
- **Category:** Test harness
- **Problem:** These scripts cannot run in hosted CI and do not report pass or fail to `pytest`.
- **Fix:** Port the important cases to `pytest` with `subprocess` plus the fake scheduler shims from E1, and mark live-cluster tests `@pytest.mark.site("gadi")` so they are skipped by default.

**Already well abstracted:** `tests/test_queue_db_schema.py` uses `tmp_path` and the `db_path=` override, and does not touch PBS. It can be parametrised for a `backend` column with very little change. `test_command_serialization.py` exercises the neutral options (`--mem`, `--cpus`, `--gpus`, `--runtime`) well.

---

## Proposed seams

1. **Backend identity on every record.** Add `backend`, `platform` and `cluster` columns to `queue` and `project_jobs` through `_QUEUE_NEW_COLUMNS` and a new `_PROJECT_JOBS_NEW_COLUMNS`, backfilling `'pbs'`. Mirror them in history `executions.yaml` records. Replace the `project_jobs` PK with `entry_id` and add `UNIQUE(backend, cluster, native_job_id)`. Rename `pbs_job_id` to `native_job_id`, keeping a compatibility view or alias for one release.

2. **Opaque `JobRef`.**
   ```python
   @dataclass(frozen=True)
   class JobRef:
       entry_id: str              # qx-uuid, qxub-owned, always present
       backend: str               # "pbs" | "slurm" | "local" | "k8s" | ...
       native_id: str | None      # scheduler's ID, opaque string
       cluster: str | None = None
   ```
   Every DB, history and status API takes or returns a `JobRef` or an `entry_id`, never a bare native ID. Each backend provides `parse_submit_output(stdout) -> native_id`, `normalize_job_id(user_input) -> native_id` (replacing the `.gadi-pbs` handling and `split(".")[0]`), and `job_env_var` (`PBS_JOBID` / `SLURM_JOB_ID` / ...).

3. **Neutral `JobState` enum with a per-backend mapping.** Use `PENDING, SUBMITTED, QUEUED, HELD, RUNNING, COMPLETED, FAILED, CANCELLED, TIMEOUT, OOM, PREEMPTED, UNKNOWN`. Each backend implements `map_state(native_state, exit_code) -> JobState`: PBS maps `Q/H/R/E/F` and exit codes above 256, Slurm maps `PD/R/CG/CD/F/TO/OOM/CA/PR`, and K8s maps pod phase and reason. Enforce it with a CHECK constraint in SQLite.

4. **Job-side reporting keyed on `entry_id`.** Inject `QXUB_ENTRY_ID` (and `QXUB_BACKEND`) into every job's environment, and change `mark_running`/`mark_complete` to take `entry_id`. Put a `StateStore` protocol behind it (SQLite on a shared filesystem today; scheduler polling or an HTTP/object-store sink for cloud) so a shared filesystem stops being a hard requirement.

5. **Neutral resource schema.** Define a `ResourceRequest`/`ResourceUsage` dataclass (`mem_mb`, `cpus`, `gpus`, `gpu_type`, `runtime_sec`, `local_disk_mb`, `queue_or_partition`, `account`, plus `extra: dict`). Use it in the DB columns, the history records, and the remote serialiser (forward only neutral options; send raw scheduler arguments through an explicit, backend-tagged passthrough).

6. **Remote submission handshake.** The remote qxub prints `QXUB-SUBMITTED entry_id=… backend=… native_id=…` immediately after submit, and the local side generates the `entry_id` and passes it as an idempotency key (`--entry-id`). This replaces `_PBS_JOB_ID_RE` and makes retries safe on every backend. Also make `env_init` optional instead of the default conda hook.

7. **Test harness.** Add `tests/conftest.py` with:
   - a `FakeBackend` (in-memory, scripted state transitions);
   - a `fake_scheduler_bin` fixture that writes shim `qsub`/`qstat`/`sbatch`/`squeue`/`sacct` executables to `tmp_path/bin`, puts them on `PATH`, and records invocations to JSON;
   - an isolated `QXUB_DB_PATH` and `XDG_CONFIG_HOME` per test;
   - a `backend` parametrisation, and a `@pytest.mark.site(...)` marker for live-cluster tests.

   Add a hosted `pytest` workflow with a backend matrix, plus optional container-based Slurm and OpenPBS jobs.

8. **Packaging.** Use `jobscripts/<backend>/…` with recursive `package_data`, add a `qxub doctor` command with backend-provided preflight checks to replace `check-nci-prerequisites.sh`, and update the recipe, skill descriptions and `python_requires` to match.
