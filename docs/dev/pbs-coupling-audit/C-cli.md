# Audit C — CLI layer: PBS / NCI-Gadi coupling

Scope: `qxub/cli.py`, `exec_cli.py`, `interactive_cli.py`, `status_cli.py`, `monitor_cli.py`, `cancel_cli.py`, `history_cli.py`, `resources_cli.py`, `config_cli.py`, `alias_cli.py`, `replay_cli.py`, `notify_cli.py`, `notifications/*.py`, `__init__.py`. All files were read in full. Line numbers were checked against the current tree. Imports into `core/scheduler.py`, `resources/`, `config/` and `platforms/` were followed only far enough to explain CLI findings.

Severity key:
- **High**: blocks any non-PBS backend.
- **Med**: wrong behaviour on another backend, or needs a refactor.
- **Low**: naming, help text or cosmetic.

## Summary

| Severity | Count |
|---|---|
| High | 11 |
| Med | 15 |
| Low | 17 |
| Incidental (non-PBS) bugs found along the way | 7 |

There is **no scheduler/backend abstraction**: `platforms/core.py` has no scheduler field. Every CLI module that touches a live job imports PBS helpers from `core/scheduler.py` directly (`qdel`, `get_job_resource_data` → `qstat -fx`, `parse_qstat_fx_output`, `job_status_from_files`, `get_default_status_dir`). Some modules also shell out to `qstat`/`qsub` themselves (`cancel_cli`, `interactive_cli`).

---

## 1. Direct scheduler command invocation

### H1 — `qxub interactive` is built on `qsub -I` (interactive_cli.py:727-761, 1306-1314, 1440-1442, 1485-1488, 1493)
- **Category:** Interactive / direct subprocess. **Severity: High.**
- **Problem:** `_build_qsub_command` hard-codes `["qsub", "-I", "-N", name, "-q", queue, "-P", project, "-l", res..., "-l", "storage=..."]`. Line 1442 appends `--` plus a script path, which is PBS Pro syntax. The result is run with `subprocess.run(qsub_cmd)` (1493) or `os.execvp("tmux", [... qsub_cmd_str])` (1485-1488). No other backend can run an interactive session.
- **Fix:** Add `backend.interactive_argv(spec: JobSpec, init_script: Path) -> list[str]`. PBS: `qsub -I ... -- script`. Slurm: `srun --pty ... bash script` or `salloc ... srun --pty`. SGE: `qrsh`. LSF: `bsub -Is`. Local: `bash script`. k8s: `kubectl run -it`/`exec`. Keep the tmux wrapper (1449-1488) generic over the returned argv.

### H2 — `cancel id` calls `qdel` directly (cancel_cli.py:19, 97-98)
- **Category:** Direct subprocess / architecture. **Severity: High.**
- **Problem:** `from .core.scheduler import qdel` and `qdel(job_id, quiet=True)`. Every `cancel` subcommand (`name`, `status`, `filter`, `all`) funnels into this, so cancellation only works on PBS.
- **Fix:** Use `backend.cancel(job_id) -> bool`, or `cancel_many` (Slurm `scancel`, SGE `qdel`, LSF `bkill`, local `kill`, k8s `kubectl delete job`).

### H3 — `cancel all` parses `qstat -u` text output and filters on PBS state letters (cancel_cli.py:471-515)
- **Category:** Direct subprocess / output parsing / job states. **Severity: High.**
- **Problem:** Builds `["qstat", "-u", user, "-q", queue]`. It skips 2 header lines and takes `parts[0]` as the job ID and `parts[4]` as the state. It filters against `--status` with `click.Choice(["Q","R","H"])` (476-481). This won't work on other schedulers.
- **Likely bug on PBS Pro too (unverified):** `qstat -u` switches to the wide format, where `S` is column 9, not 4. That format has a 4-line header, so `[2:]` would also yield the `Job ID ...` and `-----` lines as "jobs".
- **Fix:** Use `backend.list_jobs(user=, queue=, states={JobState.PENDING,...}) -> list[JobInfo]`, with a neutral `--state pending|running|held`.

### H4 — `cancel name` parses `qstat -f` (cancel_cli.py:556-610)
- **Category:** Direct subprocess / output parsing. **Severity: High.**
- **Problem:** Runs `qstat -f` and scans for `Job Id:` and `Job_Name =` lines. These keys are PBS-only. The query also isn't limited to the current user.
- **Fix:** Use `backend.list_jobs(user=current, name=..., name_regex=...)`. Slurm `squeue --me -o "%i %j"`, LSF `bjobs -J`, SGE `qstat -xml`.

### M1 — `--wait` polls `qstat <id>` exit code (cancel_cli.py:613-641)
- **Category:** Direct subprocess / polling. **Severity: Med.**
- **Problem:** Treats a non-zero `qstat job_id` return as "deleted". On a non-PBS host the call raises `FileNotFoundError`, which is caught; the loop breaks and reports a timeout. Slurm `squeue -j` returns 0 with empty output for finished jobs, so the logic would be wrong there even with the command swapped.
- **Fix:** Use `backend.wait_until_gone(job_id, timeout)`, or poll `backend.job_state(job_id)` until terminal/unknown.

### M2 — `cancel filter --job-name` uses `qstat -fx` metadata key `Job_Name` (cancel_cli.py:293, 324-336)
- **Category:** Direct helper call / PBS field names. **Severity: Med.**
- **Problem:** Calls `get_job_resource_data(job_id)` once per DB job, which runs a `qstat -fx` subprocess each time. It then reads `job_data["metadata"]["Job_Name"]`. That's PBS-only, and slow (N subprocesses).
- **Fix:** Store the job name in the tracker DB at submit time and filter in SQL. Fall back to `backend.job_info(id).name` if it's missing.
- **Done (#143):** the name is stored in `job_resources.job_name` at submission and filtered from the DB. `cancel status` and `cancel filter` make one scheduler query, `core.scheduler.list_user_jobs()` (`qstat -w -u USER`), to keep only jobs still queued or running. Rows without a stored name use that listing's (possibly truncated) name. `list_user_jobs()` is the function to move behind `backend.list_jobs()`.

### M3 — `resources show` renders raw `qstat -fx` fields (resources_cli.py:444-521)
- **Category:** Direct helper call / PBS output fields. **Severity: Med.**
- **Problem:** "Get fresh data from PBS" via `get_job_resource_data` (448-452). It displays `exec_host`, `ncpus`, `cpupercent`, `walltime`, `jobfs` and `queue_wait_seconds` from the PBS-parsed dict (462-518).
- **Fix:** Use `backend.job_accounting(job_id) -> JobAccounting` (normalized: `mem_req/used`, `runtime_req/used`, `cpus_req`, `cpu_efficiency`, `local_disk_req/used`, `host`, `queue`, `wait_s`, `exit`). Slurm: `sacct`; LSF: `bacct`; SGE: `qacct`. Rows the backend can't supply should show "n/a".

### M4 — In-job self-accounting uses `$PBS_JOBID` and `qstat -f` (interactive_cli.py:360, 369-370, 283-295)
- **Category:** PBS env vars / direct subprocess. **Severity: Med.**
- **Problem:** The generated script writes `"$PBS_JOBID"` to `job_id` and runs `qstat -f "$PBS_JOBID"` in an EXIT trap. Python then parses the result with `parse_qstat_fx_output`. On other schedulers the job ID is blank and the resources are lost; the session itself still runs.
- **Fix:** Use a backend-provided `job_id_env_var` (`SLURM_JOB_ID`, `LSB_JOBID`, `JOB_ID`, …) or a neutral exported `QXUB_JOB_ID`. Collect accounting after the session ends from Python with `backend.job_accounting()` instead of the in-job `qstat`.

---

## 2. Resource specification and option passthrough

### H5 — Neutral resource options are converted to PBS strings inside the CLI (exec_cli.py:777-805; interactive_cli.py:1245-1277; config_cli.py:1070-1112)
- **Category:** Architecture / resource rendering. **Severity: High.**
- **Problem:** `--mem/--runtime/--cpus/--disk/--volumes/--gpus` are fed to `ResourceMapper` and immediately flattened with `get_pbs_resources()` (`mem=…`, `walltime=H:MM:SS`, `ncpus=`, `jobfs=`, `storage=`, `ngpus=`). Those strings are merged with `-l` (804-805). `interactive_cli` goes further: on parse failure it appends raw `walltime=`/`mem=`/`jobfs=` into `mapper.pbs_resources` (1253, 1260, 1271). Every downstream consumer therefore only ever sees PBS syntax.
- **Fix:** The CLI should build a neutral `ResourceRequest(mem_bytes, runtime_s, cpus, gpus, gpu_type, local_disk_bytes, volumes, extra)`. Only `backend.render_resources(req)` should produce directives (`#PBS -l`, `#SBATCH --mem/--time/-c/--gres/--tmp`, `bsub -R`, k8s `resources.requests`).

### H6 — `-l/--resources` is raw PBS passthrough and is parsed for PBS keys (exec_cli.py:47-52, 681-705, 803-812; interactive_cli.py:805-810, 1277)
- **Category:** CLI option semantics. **Severity: High.**
- **Problem:** `-l` mirrors `qsub -l`. The CLI splits it on `,`/`=` and maps PBS key names (`ncpus`, `walltime`, `jobfs`, `storage`, `ngpus`) to decide which config defaults to skip (691-705). `cpus_explicit` checks `startswith("ncpus=")` (808-812). On another backend these strings are meaningless.
- **Fix:** Parse `-l` only inside the PBS backend: `backend.parse_native_resources(list[str]) -> (ResourceRequest, extras)`. Add a neutral `--native/--scheduler-arg` for backend-specific passthrough, and keep `-l` as a PBS-only alias that warns or errors elsewhere.

### H7 — CLI computes `qsub_options` and hands them down (exec_cli.py:822-871)
- **Category:** Architecture. **Severity: High.**
- **Problem:** The comment says "Extract PBS-specific options for processing". `process_job_options()` returns `qsub_options`, which is stored as `ctx.obj["options"]` for `execute_unified`. The CLI→executor contract is therefore a pre-built qsub option string, not a job spec.
- **Fix:** Pass a typed `JobSpec` (name, queue/partition, account, resources, array, mail, out/err paths, env vars, tags, execdir) in `ctx.obj`, and let the backend render submission args.

### H8 — `config shortcut set` persists PBS-rendered resources (config_cli.py:1067-1116; also 933)
- **Category:** Persisted config / resource rendering. **Severity: High.**
- **Problem:** The workflow-friendly options are run through `ResourceMapper.get_pbs_resources()`. GPUs are written as `f"ngpus={...}"` (1112). The result is stored as `definition["resources"] = ["mem=…","walltime=…"]`. Saved shortcuts are therefore PBS-only and are replayed as `-l` strings on any platform.
- **Fix:** Store neutral keys (`mem: 8GB`, `runtime: 2h`, `cpus: 4`, `gpus: 1`, `disk: …`). Render at submit time. Add a config migration that converts existing `resources:` lists via the PBS parser.

### M5 — `--email/--email-opts/--notify` carry PBS `-M/-m` semantics and contradict each other (exec_cli.py:133-150, 814-821, 835-844; config_cli.py:1001-1002)
- **Category:** CLI option semantics. **Severity: Med.**
- **Problem:** `--email-opts 'abe'` is PBS's `-m abe` letter code. `--notify` help says "Enable PBS email notification", but the value is passed as `params["notify"]` with the comment "Slack/Discord notifications" (844). That's a semantic clash.
- **Fix:** Use `--mail-user ADDR` and `--mail-on begin,end,fail` (a neutral enum that the backend maps to `-m abe`, `--mail-type=BEGIN,END,FAIL`, or `bsub -B -N`). Rename the chat flag to `--chat-notify/--no-chat-notify` and fix the help text.

### M6 — `--volumes/--storage` use the NCI `gdata/a56+scratch/a56` format (exec_cli.py:89-93; interactive_cli.py:851-855, 757-759, 1304; config_cli.py:955-959)
- **Category:** NCI-specific option. **Severity: Med.**
- **Problem:** The help reads "NCI format… Converted to PBS storage=". `interactive_cli` emits `-l storage=…` directly (759). This is a Gadi-only filesystem-mount directive.
- **Fix:** Make it a platform capability (`platform.supports("storage_directive")`). Validate the format in the platform. On backends without it, ignore with a warning or map it to a container/k8s volume mount.

### M7 — `-P/--project` with `$PROJECT` fallback (exec_cli.py:94-96; interactive_cli.py:817-822, 744-751; config_cli.py:962)
- **Category:** CLI option / NCI env var. **Severity: Med.**
- **Problem:** "PBS project code (default: $PROJECT)". `$PROJECT` is an NCI login convention. Slurm uses `--account`, LSF `-P`, and SGE uses `-P` for project and `-A` for account.
- **Fix:** Make `--account` canonical, keeping `-P/--project` as aliases. Take the default from `platform.account_env_var` or from config, not a hard-coded `PROJECT`.

### M8 — Default queue `"normal"` (interactive_cli.py:90-104; exec_cli.py:53 help)
- **Category:** Gadi default. **Severity: Med.**
- **Problem:** Hard-coded Gadi queue name as fallback. Slurm partitions differ, and local/k8s backends have no queue concept.
- **Fix:** Fall back to `platform.default_queue` (or `None` = scheduler default). Describe `--queue` as "queue/partition", with a `--partition` alias.

---

## 3. Job-ID formats and job-state codes

### H9 — `.gadi-pbs` suffix forcibly appended/stripped (status_cli.py:123-126, 174-176, 329-331)
- **Category:** Job-ID format. **Severity: High.**
- **Problem:** `status show` and `status check` (the Snakemake status hook) rewrite any ID to `f"{job_id}.gadi-pbs"` before the file and DB lookup. A Slurm ID `12345` becomes `12345.gadi-pbs`, never matches, and is reported `running` forever (393).
- **Fix:** Use `backend.normalize_job_id(raw)` and `backend.display_job_id(id)`, or store IDs exactly as returned at submit and do no rewriting.

### M9 — `monitor` job-ID validation rejects non-PBS IDs (monitor_cli.py:471-477)
- **Category:** Job-ID format. **Severity: Med.**
- **Problem:** IDs are accepted only if they contain `.` or are all digits. That rejects Slurm array tasks (`12345_7`), LSF arrays (`12345[3]`), k8s/cloud names (`job-abc12`) and qxub virtual IDs.
- **Fix:** Use `backend.is_valid_job_id(s)`, or drop validation and let the lookup report unknown IDs.

### M10 — Virtual-queue entries resolve to a `pbs_job_id` field (status_cli.py:316-324)
- **Category:** Data model naming leaking into the CLI. **Severity: Med** (schema lives in `queue/db.py`; the CLI depends on it).
- **Fix:** Rename the field to `backend_job_id` (plus `backend` name) and read it through a neutral accessor.

### L1 — PBS-style single-letter state model in monitor (monitor_cli.py:77-79, 108-111, 134-143, 159, 221-227)
- **Category:** Job states. **Severity: Low.**
- **Problem:** Internal states are `Q/R/C/F`, mirroring PBS letters. The values come from qxub's own status files, so it works anywhere, but it's PBS-flavoured and can't represent HELD or CANCELLED.
- **Fix:** Use a shared `JobState` enum (PENDING, RUNNING, HELD, COMPLETED, FAILED, CANCELLED, UNKNOWN) from the backend layer.

### L2 — Job-name fallback via `job_id.split(".")[0]` (notify_cli.py:50; notifications/discord.py:74; notifications/slack.py:140)
- **Category:** Job-ID shape. **Severity: Low.** Assumes `NNN.server`.
- **Fix:** Use `backend.display_job_id()`, or just the job ID.

### L3 — Job IDs truncated to 8 chars in tables (resources_cli.py:177, 330)
- **Category:** Job-ID format. **Severity: Low.** Assumes short numeric IDs. UUID/k8s/virtual IDs would collide visually.
- **Fix:** Use `display_job_id()` with sensible width.

### L4 — `.gadi-pbs` examples and test IDs (monitor_cli.py:378, 404-413, 460; notify_cli.py:92)
- **Category:** Help / docs. **Severity: Low.**
- **Fix:** Use neutral examples. Take the `--suffix` default from the backend.

---

## 4. Status, monitor and accounting logic

### H10 — Status directory hard-coded to `/scratch/$PROJECT(a56)/$USER/qxub/status` (interactive_cli.py:239-241, 355; monitor_cli.py:26, 68 via `core/scheduler.get_default_status_dir`; status_cli.py:336 via `job_status_from_files` default)
- **Category:** NCI path. **Severity: High.**
- **Problem:** File-based status is the backbone of `status check`, `monitor` and interactive history. The root is an NCI path with project `a56` as the default. On any other cluster or locally it doesn't exist: the interactive script's `mkdir -p` fails, and monitor/check never see completion.
- **Fix:** Make `status_dir` a platform/config value (e.g. `platform.paths.status_dir`, default `$XDG_STATE_HOME/qxub/status` or a shared-FS path from config). Pass it explicitly from the CLI.

### M11 — `monitor` cannot override the status dir (monitor_cli.py:369-494)
- **Category:** CLI surface. **Severity: Med.**
- **Problem:** `MultiJobMonitor` accepts `status_dir` (43, 68), but the command exposes no `--status-dir`/`--log-dir` (unlike `status check --log-dir`). It always uses the NCI default.
- **Fix:** Add `--status-dir` and resolve its default from platform config.

### M12 — Completion fallback and backfill rely on the PBS joblog (status_cli.py:336-338 → `job_status_from_files` → `_check_pbs_job_log_for_exit`; status_cli.py:396-427 `parse_joblog_resources`)
- **Category:** PBS exit/accounting parsing. **Severity: Med.**
- **Problem:** Detecting jobs killed by the scheduler, and backfilling resource efficiency, both depend on the PBS-written resource-usage block in the joblog. Slurm/LSF don't write it, so killed jobs (OOM/walltime before the qxub trap runs) would stay "running" and efficiency stays empty.
- **Fix:** Add `backend.terminal_status(job_id) -> (state, exit, reason) | None` (Slurm `sacct`, LSF `bhist`/`bjobs -a`, SGE `qacct`) as the second-tier check, and `backend.job_accounting()` for backfill.

### M13 — `history log` and the `joblog` concept (history_cli.py:454-495; alias_cli.py:18; exec_cli.py:832, 873-882)
- **Category:** PBS artefact in UX and data model. **Severity: Med.**
- **Problem:** "PBS log file" is a PBS-specific separate job log. The monitoring config knobs `joblog_check_interval_sec` and `walltime_offset_sec` exist to poll it over NFS. Other backends have no equivalent file, so `history log` would always be empty.
- **Fix:** Rename to "scheduler log" and make it optional (`backend.has_job_log`). On other backends, `history log` can show `backend.job_accounting()` output instead.

---

## 5. Interactive-mode specifics (beyond H1/M4)

### M14 — Init script is written to login-node `$TMPDIR` and passed by path (interactive_cli.py:1424-1442)
- **Category:** Interactive / filesystem assumption. **Severity: Med.**
- **Problem:** `tempfile.gettempdir()/qxub_interactive/init_*.sh` must be readable from the compute node. That's often false where `/tmp` is node-local, and always false for k8s/cloud.
- **Fix:** Let the backend choose the transport: a shared-FS path from platform config, an inline `bash -c` with a base64 payload (the unused `_encode_command` at 217-219 suggests this was intended), or a configmap/`kubectl cp` for k8s.

### L5 — Gadi conda locations baked into the script (interactive_cli.py:486-492)
- **Category:** NCI paths. **Severity: Low.** `"/g/data/a56/conda/miniconda3"` and `"/apps/conda"` are in the search list.
- **Fix:** Use a configurable `defaults.conda_search_paths`, or the platform config.

### L6 — PBS-specific comments and assumptions in the script (interactive_cli.py:394-398, 646-650; 1006-1031 docstring)
- **Category:** Comments / help. **Severity: Low.** "PBS doesn't source .bashrc", "Signal handling challenge in PBS", "survives PBS job ending". The behaviour (bashrc sourcing, `script` pty) is generic and fine on Slurm/SGE.
- **Fix:** Reword as generic batch-scheduler notes.

### L7 — Default runtime `"4:00:00"` in PBS walltime syntax (interactive_cli.py:45, 47)
- **Category:** Format. **Severity: Low.** It parses through the mapper, but the raw fallback at 1253 emits `walltime=4:00:00`.
- **Fix:** Use `"4h"` and the neutral `ResourceRequest`.

### L8 — `{project}` template variable from `$PROJECT` for record_dir (interactive_cli.py:188-191; replay_cli.py:33-36)
- **Category:** NCI env var. **Severity: Low.**
- **Fix:** Source it from the resolved account/platform.

---

## 6. Help text, naming and cosmetic (all Low)

### L9 — "PBS"/"qsub" throughout help, docstrings and labels
- `cli.py:2`, `cli.py:40` ("qxub - PBS job submission tools"), `cli.py:51` (example uses `${{PBS_JOBID}}`), `cli.py:54` (`-l mem=500GB`).
- `__init__.py:2`.
- `exec_cli.py:51, 53, 55, 61, 66, 72, 77, 82, 92, 95, 134, 137, 143, 149, 207`.
- `interactive_cli.py:5, 316, 809, 815, 821, 827, 967-971, 1029-1031, 1235, 1483`.
- `monitor_cli.py:2-6, 32, 49, 103-107, 390-396`.
- `cancel_cli.py:2, 25, 46, 142, 324, 560`.
- `status_cli.py:4, 248, 303, 316-327, 401, 413`.
- `history_cli.py:457, 473, 478, 487, 492, 495`.
- `resources_cli.py:451`.
- `config_cli.py:933, 960-962, 1001-1002, 1019, 1026, 1067, 1106, 1180, 1374-1386` (`shortcut show` prints a "PBS:" label).
- `notify_cli.py:29`; `notifications/__init__.py:5, 19`.
- **Fix:** Use neutral wording ("scheduler", "job", "queue/partition", "account"), and template backend names where relevant.

### L10 — `--internet` help names "copyq on NCI Gadi" (exec_cli.py:164-168; interactive_cli.py:887-891)
- **Fix:** Use generic help text. The mechanism (platform queue selection) is already neutral.

### L11 — `--jobfs` alias of `--disk` (exec_cli.py:74-78; interactive_cli.py:846-850; config_cli.py:950-954; resources_cli.py:291, 503-509 "Jobfs"/"Job Filesystem")
- **Fix:** Use `--disk`/"local scratch" as canonical. Hide `--jobfs` as a deprecated alias.

### L12 — Deprecated `qxub alias` forwards PBS-named options (alias_cli.py:197-202, 18)
- **Problem:** It forwards `--walltime/--ncpus/--jobfs` (options that don't exist on `exec`) and `--joblog`. It also invokes `qxub --env ...` without the `exec` subcommand (181-233), so it's broken regardless.
- **Fix:** Delete the command (it's already deprecated).

### L13 — Notification architecture assumes internet-less compute nodes (notify_cli.py:7-10, 38-39)
- **Problem:** The design is "Job script calls `qxub exec --internet ... qxub notify send`", which suits Gadi. On Slurm, cloud or local runs this extra hop isn't needed.
- **Fix:** Add `platform.compute_has_internet`. When true, send directly from the job epilogue.

### L14 — `--array` grammar documented in PBS form (exec_cli.py:151; config_cli.py:1003)
- **Problem:** `1-100:2` happens to match Slurm and SGE, but the backend must render it (`-J`, `--array`, `-t`, `name[1-100:2]` for LSF).
- **Fix:** Define a neutral grammar. Only backends render it.

### L15 — Exec monitoring knobs are joblog/walltime-specific (exec_cli.py:873-882)
- **Fix:** Namespace them under the PBS backend config, or make them backend-advertised.

### L16 — `resources list --queue` filter and Queue column (resources_cli.py:49-54, 172, 255)
- **Problem:** The concept is mostly portable; the help example "normal, hugemem" is Gadi-specific.
- **Fix:** Neutral example text.

### L17 — Shortcut/alias option merge logic is copy-pasted 3× in exec and 3× in interactive (exec_cli.py:290-536; interactive_cli.py:1036-1160)
- **Category:** Architecture (not PBS per se). Each copy hard-codes the PBS-ish keys (`resources`, `project`, `queue`, `email_opts`). Interactive copies omit keys that exec supports (`project`, `name`, `mem`, …).
- **Fix:** Use a single `JobSpec.merge(shortcut|alias|cli)` helper.

---

## 7. Architecture

### H11 — No backend interface; CLI imports PBS helpers directly
- **Where:**
  - `status_cli.py:20` (`job_status_from_files`)
  - `monitor_cli.py:26` (`get_default_status_dir`, `job_status_from_files`)
  - `cancel_cli.py:19` (`qdel`), `cancel_cli.py:293` (`get_job_resource_data`)
  - `resources_cli.py:448` (`get_job_resource_data`)
  - `interactive_cli.py:287` (`parse_qstat_fx_output`)
  - `cancel_cli.py:13, 493-499, 567-569, 625-627` and `interactive_cli.py:10, 1493` (raw `subprocess` scheduler calls)
- **Category:** Architecture. **Severity: High.**
- **Fix:** Introduce `qxub/backends/base.py::SchedulerBackend` (see Proposed seams), with `PbsBackend` wrapping the existing `core/scheduler.py`. Resolve it via `get_backend(platform_config)` (e.g. a `platform.scheduler: pbs|slurm|sge|lsf|local|k8s` field). CLI modules should import only the backend facade.

### M15 — `interactive` has its own submission pipeline
- **Problem:** `interactive_cli` does its own config defaults, resource mapping and command construction (31-175, 1208-1314). It bypasses `config.handler.process_job_options` and `execution/`, so platform features (auto queue, GPU queue, account defaults, remote delegation) diverge from `exec`.
- **Fix:** Build the same `JobSpec` as `exec` and call `backend.interactive_argv(spec, script)`.

---

## 8. Already reasonably abstracted (keep)
- **Workflow-friendly option names:** `--mem/--memory`, `--runtime/--time`, `--cpus/--threads`, `--gpus`, `--gpu-type`, `--disk` are scheduler-neutral at the surface. Only their rendering (H5) is PBS.
- **Queue auto-selection:** `--queue auto`, `--internet` and `--gpus` delegate to `platforms` (`interactive_cli._get_internet_queue` → `platform.select_queue`; exec sets `queue="auto"` at 856-862). This is the right seam.
- **Remote/delegated execution** (exec_cli.py:587-674): passes CLI options verbatim to a remote qxub that owns its own platform. It's backend-agnostic.
- **Execution contexts:** conda/module/singularity/default are orthogonal to the scheduler and work anywhere.
- **File-based status sentinels:** `job_started_*`, `final_exit_code_*` written by qxub's own job wrapper are scheduler-independent in principle. Only the path (H10) and the PBS-joblog fallback (M12) are coupled.
- **Status vocabulary in the DB and CLI:** `submitted/running/completed/failed/cancelled`, plus `status check` output formats (snakemake/json/exitcode). These are neutral.
- **Other neutral modules:** `history_cli` (apart from `log`), `config_cli` core commands, `replay_cli` (util-linux `script`/`scriptreplay`), and `notifications/*` channels.
- **Tags and vars:** `--tag/--tags/--var/--vars` are neutral.

---

## CLI surface: PBS-flavoured options

| Option (aliases) | Commands | Current PBS/NCI semantics | Proposed neutral name / semantics |
|---|---|---|---|
| `-l`, `--resources` | exec, interactive, `config shortcut set` (`-l`), `config alias set` (`--resources`), `alias` | Raw `qsub -l` strings; parsed for `ncpus/walltime/jobfs/storage/ngpus` | `--native ARG` / `--scheduler-arg` (repeatable, passed verbatim to the active backend). Keep `-l` as a PBS-only alias that warns or errors on other backends. Persist neutral keys, not `-l` strings |
| `-q`, `--queue` | exec, interactive, shortcut set, alias set, `cancel all`, `resources list` | PBS queue; default `normal` | `--queue` (alias `--partition`) = "queue/partition". Default from `platform.default_queue` or none. `auto` stays |
| `-P`, `--project` | exec, interactive, shortcut set, alias | `qsub -P`; default `$PROJECT` | `--account` (aliases `-A`, `-P`, `--project`). Default from platform-configured env var or config |
| `-N`, `--name` | exec, interactive, shortcut set, alias set | `qsub -N` | `--name` (keep `-N`). Backend sanitizes name rules |
| `--mem`, `--memory` | exec, interactive, shortcut set | → `mem=` | Keep. Neutral `ResourceRequest.mem`. Remove "Converted to PBS" from help |
| `--runtime`, `--time` | exec, interactive, shortcut set | → `walltime=H:MM:SS` | Keep. `ResourceRequest.runtime` |
| `--cpus`, `--threads` | exec, interactive, shortcut set | → `ncpus=` | Keep. `ResourceRequest.cpus` (backend chooses `-c` vs `-n` semantics) |
| `--gpus` | exec, shortcut set | → `ngpus=` (also written raw in config_cli:1112) | Keep. `ResourceRequest.gpus` (Slurm `--gres/--gpus`, k8s `nvidia.com/gpu`) |
| `--disk`, `--jobfs` | exec, interactive, shortcut set | → `jobfs=` (Gadi local SSD) | `--disk` / `--local-disk` canonical (Slurm `--tmp`, k8s `ephemeral-storage`). `--jobfs` becomes a hidden deprecated alias |
| `--volumes`, `--storage` | exec, interactive, shortcut set | NCI `storage=gdata/x+scratch/y` | `--volumes` as a platform capability. Format validated per platform; ignored with a warning if unsupported, or mapped to container/k8s mounts |
| `--email` | exec, shortcut set | `qsub -M` | `--mail-user` (keep `--email` alias) |
| `--email-opts` | exec, shortcut set | `qsub -m abe` letters | `--mail-on begin,end,fail` (neutral enum, backend-mapped) |
| `--notify` / `--no-notify` | exec | Help says PBS email; code means Slack/Discord | `--chat-notify/--no-chat-notify` (separate from mail) |
| `--array` | exec, shortcut set | `qsub -J 1-100:2` | `--array START-END[:STEP]` neutral grammar. Backend renders `-J`/`--array`/`-t`/`name[...]` |
| `--internet` | exec, interactive, shortcut set | Selects `copyq` on Gadi | Keep. Help text generic. Platform capability `internet_queue` / `compute_has_internet` |
| `--out` / `--err` / `--log-dir` | exec, alias | `-o/-e`; joblog alongside | Keep. Make the "scheduler log" (`joblog`) optional per backend |
| `--joblog` | alias (deprecated) | PBS joblog | Remove with `qxub alias` |
| `cancel all --status Q/R/H` | cancel all | PBS state letters | `--state pending\|running\|held` (neutral `JobState`) |
| `cancel all --user/-u`, `--queue/-q` | cancel all | `qstat -u`, `-q` | Keep names. Implement via `backend.list_jobs(user=, queue=)` |
| `cancel name --pattern` | cancel name | Matches PBS `Job_Name` from `qstat -f` | Keep. Implement via `backend.list_jobs(name_regex=)` or the DB |
| `cancel filter --job-name` | cancel filter | PBS `Job_Name` via `qstat -fx` | Keep. Use the stored name in the DB |
| `monitor --suffix` | monitor | Strip `.gadi-pbs` | `--strip-suffix`. Default from `backend.display_job_id`. Add `--status-dir` |
| `status check --log-dir` | status check | "PBS job log files" | Keep. Neutral help ("scheduler log dir"). Add `--status-dir` |
| `notify send --job-id` | notify send | "PBS job ID" | Keep. Neutral help |
| `qxub interactive` (whole command) | interactive | `qsub -I` | Keep UX. Backend-provided interactive launcher |

---

## Proposed seams (what the CLI needs from a backend interface)

```python
class JobState(Enum): PENDING, RUNNING, HELD, COMPLETED, FAILED, CANCELLED, UNKNOWN

@dataclass
class ResourceRequest:
    mem: int|None; runtime: int|None; cpus: int|None; gpus: int|None
    gpu_type: str|None; local_disk: int|None; volumes: list[str]
    native: list[str]           # verbatim backend args (old -l)

@dataclass
class JobSpec:
    name, queue, account, resources: ResourceRequest, array, mail_user, mail_on,
    out, err, execdir, env: dict, tags: list, script/command, interactive: bool

class SchedulerBackend(Protocol):
    name: str
    capabilities: set[str]  # interactive, arrays, gpus, storage_directive, local_disk,
                            # mail, job_log, accounting, held_state, internet_queue
    # submission
    def render_resources(self, req: ResourceRequest) -> list[str]
    def parse_native_resources(self, args: list[str]) -> tuple[ResourceRequest, list[str]]
    def submit(self, spec: JobSpec) -> str                      # job id
    def interactive_argv(self, spec: JobSpec, init_script: Path) -> list[str]
    def init_script_transport(self) -> Literal["shared_fs","inline","upload"]
    # identity
    job_id_env_var: str                                         # PBS_JOBID / SLURM_JOB_ID / ...
    def normalize_job_id(self, raw: str) -> str
    def display_job_id(self, job_id: str) -> str
    def is_valid_job_id(self, s: str) -> bool
    # query / control
    def job_state(self, job_id) -> JobState
    def list_jobs(self, *, user=None, queue=None, states=None, name=None, name_regex=None) -> list[JobInfo]
    def cancel(self, job_id) -> bool
    def wait_until_gone(self, job_id, timeout) -> bool
    # completion / accounting
    def terminal_status(self, job_id) -> tuple[JobState, int|None, str|None] | None   # replaces PBS joblog fallback
    def job_accounting(self, job_id) -> JobAccounting | None                           # replaces qstat -fx parsing
    def interpret_exit_code(self, raw: int) -> ExitInfo                                # 271 / 256+sig etc. live here
    def job_log_path(self, spec) -> Path | None                                        # PBS joblog; None elsewhere
```

Platform-level values the CLI currently hard-codes:
- `scheduler` (selects the backend)
- `default_queue`
- `account_env_var` / `default_account`
- `paths.status_dir` (replaces `/scratch/$PROJECT/$USER/qxub/status`)
- `paths.shared_tmp` (for interactive init scripts)
- `conda_search_paths`
- `compute_has_internet` / `internet_queue`
- `supports_storage_directive`

CLI wiring change: every subcommand resolves `backend = get_backend(platform)` once. `exec` and `interactive` build a `JobSpec`. `status`, `monitor`, `cancel` and `resources show` call only the backend facade and the tracker DB, never `core.scheduler` or `subprocess` directly.

---

## Incidental bugs found (not PBS-specific; for the maintainer)
1. **exec_cli.py:721-723:** `logger` is undefined in `exec_cli` (only `import logging` is done), so the deprecated root-level config-key path raises `NameError`.
2. **interactive_cli.py:1285 vs 1386:** `session_id` (the UUID used for the status-dir name in the script) is overwritten with `qxi_{name}_{ts}` before `_process_session_completion(session_id)` at 1506. Post-session history and resource updates therefore never find their files.
3. **interactive_cli.py:217-219:** `_encode_command` uses `base64` without importing it (dead code today).
4. **cancel_cli.py:493-511:** `qstat -u` header parsing and state column index appear wrong even on PBS Pro (see H3). Verify on Gadi.
5. **config_cli.py:520-546:** `subcommand_def` is used at 543 but only defined in the hierarchical branch. A legacy flat alias with `subcommand: sing` raises `NameError`.
6. **alias_cli.py:181-233:** Builds `qxub --env …` without the `exec` subcommand and forwards non-existent options (`--walltime`, `--ncpus`, `--jobfs`), so the command is non-functional.
7. **notifications/slack.py:194-195, discord.py:112-113:** The footer suggests `qxub history show {job_id} --output/--logs`. That command doesn't exist (`history show` takes a recipe hash; the real commands are `history out/err/log JOB_ID`). Also, cancel subcommands `return 0/1` from Click callbacks, which doesn't set the process exit code in standalone mode.
