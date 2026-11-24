# Audit A: PBS coupling in submission, templates, execution and monitoring

**Scope:** `qxub/core/scheduler.py`, `qxub/core/templates.py`, `qxub/core/__init__.py`, `qxub/jobscripts/*.pbs`, `qxub/execution/*.py`. Some lines outside this scope (`config/handler.py`, `resources/parser.py`, `exec_cli.py`, the CLI callers) were checked only to confirm a finding.
**Mode:** Read-only. No project files were modified.
**Line numbers:** Checked against the files on disk on 2026-11-16.

Severity key: **High** = blocks every non-PBS backend (or is Gadi-only and wrong even on other PBS sites). **Med** = wrong behaviour or needs a refactor. **Low** = naming, docs or cosmetic.

## Summary counts

| Severity | Count |
|---|---|
| High | 13 |
| Med | 20 |
| Low | 8 |
| **Total** | **41** |

A further 5 incidental bugs that are not about PBS are listed in section H. They are not counted above.

---

## A. Calls to scheduler commands (submit, cancel, query)

**A1. `core/scheduler.py:215-249` (`qsub()`): submission. High.**
Submission takes a pre-built command *string*. It refuses to run unless the string starts with `qsub` (`:230`), runs it through `shlex.split`, and treats the stdout as the job ID with no parsing or validation (`:249`). No other backend can get past the `cmd[:4] == "qsub"` check.
*Suggest:* `Scheduler.submit(spec: JobSpec, script: Path) -> SubmitHandle`. Each backend builds its own argv (a list, not a string) and parses its own job ID (Slurm: `Submitted batch job N` or `--parsable`; LSF: `Job <N> is submitted`; SGE: `Your job N ...`).

**A2. `core/scheduler.py:240-241`: exit code 166 mapped to 2. Med.**
The comment reads "Map PBS validation errors to exit code 2". 166 is the PBS Pro code for an invalid request, and other schedulers use other codes.
*Suggest:* add `Scheduler.classify_submit_error(returncode, stderr) -> SubmitErrorKind` (VALIDATION / TRANSIENT / AUTH / UNKNOWN), and have the CLI map kinds to exit codes.

**A3. `core/scheduler.py:252-284` (`qdel()`): cancellation. High.**
Cancellation hard-codes `["qdel", job_id]`. It is called on Ctrl-C from `execution/context.py:448-450` and `execution/core.py:375-377`.
*Suggest:* `Scheduler.cancel(job_id) -> bool` (scancel, bkill, qdel for SGE, or a local process kill).

**A4. `core/scheduler.py:453-491` (`get_job_resource_data`): `qstat -fx`. High.**
The function shells out to `qstat -fx`. It does not catch `FileNotFoundError`, so it raises on any host without `qstat`. It is used by `cancel_cli.py:330` and `resources_cli.py:452`.
*Suggest:* `Scheduler.job_info(job_id) -> JobInfo | None` (Slurm: `sacct --json` / `scontrol show job`; LSF: `bjobs -json`; SGE: `qacct -j`). Return a normalised `JobInfo`.

**A5. `core/scheduler.py:494-566` (`parse_qstat_fx_output`): PBS attribute names. Med.**
The parser looks for PBS attribute names: `Exit_status`, `Resource_List.*`, `resources_used.*`, `exec_host`/`exec_vnode`, `qtime/stime/etime/mtime/ctime/obittime`, `Job_Id`, `job_state`. It is also imported directly by `interactive_cli.py:287`.
*Suggest:* move it into a `PBSScheduler` backend. Its output should be a backend-neutral `JobInfo` dataclass (requested, used, timings, hosts, exit), not a dict keyed by PBS names.

**A6. `core/scheduler.py:569-671` (`_enhance_resource_data`): PBS resource names. Med.**
Efficiency maths reads PBS resource names: `mem`, `walltime`, `cpupercent` (`:650`) and the Gadi-specific `jobfs` (`:657-667`). It also calls `parse_exec_host` / `parse_exec_vnode` (`:617-620`), whose docstrings in `resources/parser.py:165-167` use Gadi hostnames (`gadi-cpu-clx-2841/3`).
*Suggest:* compute efficiency from the normalised `JobInfo` (cpu_seconds, mem_bytes, walltime_seconds, scratch_bytes). Each backend fills in the fields it has.

## B. Building the submit command and moving data into the job

**B1. `execution/context.py:260`: hand-built `qsub` command line. High.**
*Done (5b): `PbsBackend.submit(JobSpec)` builds the `qsub` argv; the execution layer no longer builds a command string.*
The live CLI path builds the full command with an f-string: `f'qsub -v {submission_vars} {ctx_obj["options"]} {template}'`. The same pattern appears at `execution/core.py:287`, `execution/executors.py:52-54` and `execution/submit.py:275`.
*Suggest:* replace all four with `scheduler.submit(JobSpec, rendered_script)`. The execution layer should never see argv.

**B2. `qsub -v` as the transport for job parameters. High.** *(Done: a sidecar settings file, `<runs dir>/<run ID>.env`, sourced by the job script; `qsub -v` carries only `QXUB_JOB_FILE`.)*
Sites: `execution/context.py:116, 136-165`, `execution/core.py:204-231`, `execution/submit.py:233-270`.
All job inputs (`cmd_b64`, `cwd`, `out`, `err`, `env`, `mods`, `sif`, `bind`, `pre/post`, `user_vars_b64`, `QXUB_SHARED_DB`, `qxub_notify*`) travel as one comma-separated `-v` list. That format belongs to PBS: Slurm uses `--export=` with different quoting, LSF has no `-v`, and Kubernetes uses pod env. Several values are also unquoted (`cwd=`, `out=`, `err=`, `QXUB_SHARED_DB=`), so paths that contain commas or spaces break. Singularity bind lists also use commas (`a:/a,b:/b`).
*Suggest:* render the parameters into the job script itself (a heredoc or `export` lines inside the script), or write a sidecar `job.env` / `job.json` next to the script. Either removes the dependence on scheduler env-passing and on quoting rules.

**B3. `config/handler.py:304-316` (`build_qsub_options`): opaque PBS flag string. High.**
*Done (5b): `config.handler.job_spec()` builds a `JobSpec`. `build_qsub_options` remains only as a log string derived from the backend.*
This builder is outside scope but is consumed here as `ctx_obj["options"]` / `qsub_options`. It emits `-N -q -P -l k=v -o <joblog> -M -m` as a single string, so the execution layer forwards PBS flags without seeing them.
*Suggest:* pass a structured `JobSpec(name, queue/partition, account, resources, stdout, stderr, mail)` and give each backend a `render_directives(JobSpec)`.

**B4. `execution/submit.py:65-69, 103-108, 177-181, 198-201`: PBS terms in the public API. Med.**
The public API documents `resources` as "same as PBS `-l` keys", `queue` as "PBS queue (default normal)" and `project` as "PBS project code (default `$PROJECT`)". The resources dict is turned straight into `key=value` strings.
*Suggest:* give the API neutral names: `cpus`, `memory`, `walltime`, `gpus`, `scratch`, `queue`/`partition`, `account`. Each backend maps them, and an `extra_directives` escape hatch covers anything else.

**B5. `ncpus` parsed from resource strings. Med.**
Sites: `execution/context.py:27-36` (regex `ncpus=(\d+)`), used at `:318` and `:371`; `execution/submit.py:217, 294, 339`.
The CPU count is pulled out of PBS resource strings. A Slurm-style `--cpus-per-task` would be missed.
*Suggest:* `JobSpec.cpus: int` as the source of truth, with resource strings rendered from it.

**B6. `execution/context.py:326-339`: PBS error handling. Med.**
*Partly done (5b): failures raise `backends.SubmitError` (`QsubError` is an alias); the error kinds are still to do.*
The code catches `QsubError` and calls `sys.exit(qsub_exc.exit_code)`, so the CLI exits with whatever the PBS client returned. The comment says "exit with the PBS exit code".
*Suggest:* use a backend-neutral `SubmissionError(kind, backend_exit_code, stderr)` and map its kinds to CLI exit codes.

**B7. Four overlapping submission paths. Med.**
*Done (5b): one pipeline, `execution.submit.dispatch`. `executors.py` and `core.submit_and_monitor_job` are deleted, and `submit_job` now passes the same job variables as `qxub exec`.*
Sites: `execution/context.py:220` (`execute_unified`, the only one `exec_cli.py:15` uses), `execution/core.py:236` (`submit_and_monitor_job`), `execution/executors.py:21-58` (`_handle_terse_execution`) and `execution/submit.py:57` (`submit_job`).
They have drifted apart:
- `core.build_submission_variables` expands `${var}`; the others do not.
- `submit.py` omits `user_vars` and the notify variables.
- `core.py` exits early in quiet mode; `context.py` does not.

Adding a backend would mean changing all four.
*Suggest:* collapse them into one `JobSpec` builder plus one `submit()` pipeline. Delete `executors.py` / `core.submit_and_monitor_job`, which only the `execution` package's own re-exports reference.

## C. Job script templates (`qxub/jobscripts/*.pbs`)

**C1. Lines 1-7 of all five templates: `#PBS` directives. High.**
Every template starts with:
`#PBS -N <name>`, `#PBS -P a56`, `#PBS -l storage=gdata/a56`, `#PBS -S /bin/bash`, `#PBS -W umask=0027`, `#PBS -j oe`.
Project `a56` and `storage=gdata/a56` are specific to NCI Gadi and to one project, so they are wrong even on other PBS sites. `-l storage` and `-W umask` are Gadi/PBS-only, and `-j oe` changes where stdout/stderr go (see C6).
*Suggest:* drop all scheduler directives from the templates and keep a scheduler-agnostic body. The backend's `JobScriptRenderer` prepends a directive header (`#SBATCH`, `#BSUB`, `#$`, `#PBS`, or none for local/k8s), built from `JobSpec`.

**C2. `PBS_JOBID` as the key for status files. High.** *(Done: status files are named after `QXUB_RUN_ID` and written to `QXUB_STATUS_DIR`; a `native_<job ID>` file maps the scheduler ID to the run ID. Notifications still show the PBS ID.)*
Sites, using `qdefault.pbs` as the example:
- `:56-57` (`job_started_${PBS_JOBID:-unknown}`)
- `:71, 74, 86, 88` (DB updates)
- `:104-105, 197-198, 206, 221-222, 245-246, 256-257` (`final_exit_code_…` / `main_started_…`)
- `:114, 266` (notify)

The same pattern appears throughout `qconda.pbs`, `qconda_mod.pbs`, `qmod.pbs` (`:57-58` and onwards) and `qsing.pbs` (`:58-59` and onwards).
On any non-PBS backend every job writes `job_started_unknown` / `final_exit_code_unknown`. Concurrent jobs then collide, and the monitor, which waits for `job_started_<real id>` (`scheduler.py:994-995`), never sees a start. It waits forever (see E2).
*Suggest:* generate a qxub-owned run ID (the existing virtual ID `qx-<uuid>`) before submission, pass it as `QXUB_RUN_ID`, and key all status files on it. That also removes the race between status files and the job ID. If the scheduler job ID is wanted, the backend can export it as `QXUB_SCHED_JOBID` from `$SLURM_JOB_ID` / `$LSB_JOBID` / `$JOB_ID` / `$PBS_JOBID`.

**C3. `PBS_JOBNAME` in notifications. Low.**
Sites: `qdefault.pbs:116, 268` and the equivalents in the other templates (`qconda.pbs:105, 304`, etc.).
*Suggest:* `QXUB_JOB_NAME`, passed in by qxub.

**C4. Hard-coded NCI path in the DB status helper. High.**
Sites: `qdefault.pbs:67, 82`; `qconda.pbs:58, 74`; `qconda_mod.pbs:59, 75`; `qmod.pbs:28, 43`; `qsing.pbs:29, 44`.
`update_job_status()` runs `sys.path.insert(0, '/g/data/a56/software/qsub_tools')`. This only works on Gadi for project a56. Because of the `except: pass`, failures elsewhere are silent and DB status is never updated.
*Suggest:* call the same `qxub` that submitted (the `qxub_path` that is already passed) through a CLI subcommand such as `qxub _internal mark-running`. Or have it write status files only and let the submit side or monitor update the DB.

**C5. Joblog path derived from the stdout path. Med.**
Sites: `qdefault.pbs:73`, `qconda.pbs:64`, `qconda_mod.pbs:65`, `qmod.pbs:34`, `qsing.pbs:35`.
`joblog = '${out%.out}.log'` assumes the PBS joblog sits next to stdout with a `.log` suffix. That depends on `-o <joblog>` in `build_qsub_options`.
*Suggest:* pass the scheduler-log path explicitly (`QXUB_JOBLOG`), or drop it once resource collection goes through `Scheduler.job_info`.

**C6. Quiet mode depends on PBS output spooling. Med.**
Sites: `qdefault.pbs:42-44, 159-163`; `qconda.pbs:43, 247`; `qconda_mod.pbs:44, 257`; `qmod.pbs:87, 191`; `qsing.pbs:88, 185`.
Quiet mode relies on "output will go to PBS joblog": `-j oe` plus `-o <joblog>` merge the script's stdout and stderr into the joblog, which PBS copies back when the job ends. Slurm writes `slurm-%j.out` live; Kubernetes has pod logs.
*Suggest:* make "where the job's own output goes" a field of `JobSpec` (stdout/stderr/merge). In quiet mode the template can redirect to an explicit log file rather than depend on spooling.

**C7. `core/templates.py:33-39`: template lookup has one axis. Med.**
Templates are keyed by environment context only (`conda`, `module`, `conda_module`, `singularity`, `default`), and the files are named `.pbs`. There is no backend axis, and the five near-copies differ only in how they activate the environment (≈250-320 lines each, mostly duplicated).
*Suggest:* a single `job_body.sh.j2` with environment-activation snippets (`activate_conda`, `load_modules`, `run_container`), plus a per-backend directive header. `get_template(context, backend)`.

**C8. Assumes the Environment Modules / Lmod `module` command. Med.**
Sites: `qsing.pbs:120-122` (`module load singularity`), `qmod.pbs:119-127, 239-248`, `qconda_mod.pbs:130-137`.
These assume `module` exists and that a module named `singularity` exists. On many sites the runtime is `apptainer`, and containers and k8s have no module system.
*Suggest:* make the container runtime (`singularity` / `apptainer` / `docker` / `podman`) and its setup command platform configuration, not template constants.

**C9. Gadi-specific hints in conda error messages. Low.**
Sites: `qconda.pbs:137, 153`; `qconda_mod.pbs:148, 164`.
The messages suggest `-l storage=gdata/conda` and "mounted in your PBS job".
*Suggest:* platform-provided hint text, or generic wording.

**C10. Notifications submit another job from inside the job. Med.**
Sites: `qdefault.pbs:111-119, 263-271` (and the equivalents in every template).
Notifications are sent via `qxub exec --internet ... qxub notify send`, which submits a new job from inside the running job. This assumes compute nodes can submit jobs and that an "internet" queue exists (Gadi `copyq`). Neither holds on many Slurm sites, local mode or Kubernetes.
*Suggest:* add `Platform.notify_strategy` (direct HTTP from the node / submit to a data-mover queue / send from the submit host after the monitor sees completion).

## D. Exit status and job log handling

**D1. `core/scheduler.py:375-403` (`read_exit_from_joblog_file`): exit status from the Gadi epilogue. High.**
The exit status is parsed from `Exit Status:\s*(-?\d+)` in the `.log` file. That block is the **NCI Gadi epilogue** ("Resource Usage on …", see `resources/parser.py:288-300`), not a generic PBS Pro feature. It is the only way qxub detects jobs killed by the scheduler (walltime, OOM, `qdel`).
*Suggest:* `Scheduler.final_state(job_id) -> ExitStatus | None` (sacct State/ExitCode, `bjobs`/`bhist`, `qacct`, k8s pod `terminated.exitCode/reason`). The joblog parser can stay as a PBS-Gadi fallback.

**D2. `core/scheduler.py:406-450` (`_check_pbs_job_log_for_exit`): directory scan with a fixed-width match. Med.**
The function globs every `*.log` in `log_dir`, reads each in full, and matches the fixed-width string `"Job Id:             {job_id}"` (`:434`). That ties it to the Gadi epilogue format, and it reads every file in the directory on each status check (`status_cli.py:336`, `monitor_cli.py:115`).
*Suggest:* the same as D1. Look up by known path or by scheduler query, never by directory scan.

**D3. Exit codes are not normalised. Med.**
Sites: `core/scheduler.py:1262` (returns the raw joblog exit code, which can be negative, e.g. `-29` for walltime), `execution/context.py:444` (`sys.exit(exit_status)`), `execution/core.py:371` (`os._exit`).
PBS uses negative codes and 256+signal; Slurm reports `ExitCode=N:S` and states such as `TIMEOUT` / `OUT_OF_MEMORY`; LSF uses `TERM_*`; k8s uses `reason`. Passing these raw to `sys.exit` truncates them mod 256 (so `-29` becomes `227`) and loses the reason.
*Suggest:* an `ExitStatus(code:int, signal:int|None, reason: Enum{OK, FAILED, WALLTIME, MEMORY, CANCELLED, NODE_FAIL})`, normalised per backend, plus a single `to_cli_exit_code()` policy.

**D4. Walltime gating watches the joblog. Med.**
Sites: `core/scheduler.py:1046-1055` (gating) and `:1150-1158, 1229-1262` (joblog polling inside `stream_job_output_with_status_files`).
The design is reasonable and generic: defer the expensive check until walltime has passed, then rate-limit it. But what it checks is the PBS joblog.
*Suggest:* keep the gating logic and swap the probe for `Scheduler.final_state()`.

**D5. `core/scheduler.py:819-959` (`start_background_resource_collection`): background fork. Med.**
The function forks, sleeps 60 s "for PBS to finalize joblog", then parses the PBS/Gadi joblog (`resources/parser.parse_joblog_resources`).
*Suggest:* `Scheduler.job_info()` with a backend-specific "accounting ready" delay (sacct is usually immediate; PBS `qstat -fx` needs history enabled).

**D6. Joblog path derived with `str.replace`. Low.**
Sites: `core/scheduler.py:756, 1042`.
`str(out_file).replace(".out", ".log")` assumes the PBS joblog naming convention. It also replaces the *first* `.out` anywhere in the path, so a directory such as `/x/.outputs/` breaks it.
*Suggest:* pass the log path explicitly (`joblog_file` already exists), and if it must be derived, use `Path.with_suffix`.

**D7. `core/scheduler.py:302-372` (`job_status_from_files`): PBS state letters. Low.**
The function returns `"Q"/"R"/"C"/"F"`, PBS-style letters, with `"C"` borrowed from Torque.
*Suggest:* `JobState` enum (PENDING, RUNNING, COMPLETED, FAILED, CANCELLED, UNKNOWN), mapped from each backend's codes.

## E. Monitoring and shared-filesystem assumptions

**E1. Monitoring requires a shared filesystem. High.**
Sites: `core/scheduler.py:991-995` (`status_dir = Path(out_file).parent / "status"`), `:1110-1264` (tails `out`/`err` files), `execution/context.py:417-423` (creates the out/err files on the submit host), and the templates' `STATUS_DIR="${out%/*}/status"` (`qdefault.pbs:52`, `qmod.pbs:18`, `qsing.pbs:19`, …).
Monitoring and streaming assume the submit host and the compute node share a POSIX filesystem. That is true on HPC clusters but false for Kubernetes or cloud batch, and in remote-delegated mode.
*Suggest:* a `JobMonitor` / `OutputStream` interface. `SharedFSMonitor` keeps the current behaviour; a k8s backend implements it via pod logs and pod status; local mode uses the process handle.

**E2. `core/scheduler.py:1003-1030`: no timeout or fallback while waiting for the start. Med.**
The loop waits for `job_started_*` with no timeout and never asks the scheduler. A job deleted while queued, rejected by the scheduler, or writing the wrong ID (C2) leaves the client spinning forever. The loop only exits on a status file.
*Suggest:* every N seconds, call `Scheduler.status(job_id)` as a fallback. If it reports terminal without a status file, use `final_state()`.

**E3. `core/scheduler.py:287-299` (`get_default_status_dir`): hard-coded NCI scratch path. High.**
The default is `/scratch/{$PROJECT or "a56"}/{$USER}/qxub/status`, which is NCI-specific and defaults to project `a56`. It is used by `monitor_cli.py:68` and `job_status_from_files`.
*Suggest:* `Platform.paths.status_dir` from config, falling back to `$XDG_STATE_HOME/qxub/status` or the job log directory.

**E4. Tight polling intervals. Low.**
Sites: `core/scheduler.py:1030, 1264`.
The client checks for files every 0.5 s and 0.2 s. That is fine locally but noisy on NFS/Lustre at scale, and remote backends will need different cadences.
*Suggest:* backend or monitor-specific poll intervals.

**E5. Ctrl-C cleanup tells users to run `qdel`. Med.**
Sites: `execution/context.py:445-458`, `execution/core.py:372-385`.
Ctrl-C calls `qdel`, and the failure message tells the user to "manually run: qdel {job_id}".
*Suggest:* `scheduler.cancel()` plus `scheduler.cancel_hint(job_id)` for the message.

## F. Architecture and naming

**F1. No scheduler interface. High.**
Sites: `core/__init__.py:3-25`, `core/scheduler.py` (whole file).
`core` exports bare functions `qsub`, `qdel`, `get_job_resource_data`, `job_status_from_files`. `scheduler.py` mixes TTY UI (spinners, `print_status`), PBS CLI wrappers, PBS output parsing, Gadi joblog parsing, monitoring and background forking. There is no seam to plug in another backend.
*Suggest:*
- Split it into `qxub/scheduler/base.py` (ABC), `pbs.py`, `slurm.py`, `local.py`, …
- Move the UI to `qxub/ui/progress.py` and monitoring to `qxub/monitor/`.
- Select the backend from platform config.

**F2. `QsubError` in the public API. Low.**
Sites: `core/scheduler.py:201-212`; re-exported as public API in `execution/__init__.py:25-26, 94`; `execution/submit.py:123, 127, 313-319`.
*Suggest:* rename to `SubmissionError` and keep `QsubError = SubmissionError` as a deprecated alias.

**F3. `pbs_job_id` field names. Low.**
Sites: `execution/context.py:312, 356`; `execution/submit.py:297, 328`.
These field names, used in queue-DB calls, assume PBS.
*Suggest:* `scheduler_job_id` plus a `backend` column. This needs a DB migration, which is outside this scope.

**F4. PBS and Gadi wording in docs and comments. Low.**
Sites:
- `core/scheduler.py:1-5` ("running qsub jobs on PBS Pro"), `:316, :415` (`"12345.gadi-pbs"`)
- `execution/submit.py:4-21, 38-46, 82-87, 106-108`
- `execution/mode.py:4-6, 33, 41-71` (`nci_gadi`, `/scratch/a56`)
- `execution/context.py:91, 324, 350`

*Suggest:* use neutral wording ("scheduler job ID"), and move examples to per-platform docs.

**F5. `execution/mode.py:15-20, 31-33, 94-98`: `ExecutionMode.LOCAL` naming collision. Med.**
`ExecutionMode.LOCAL` is documented as "direct PBS submission". Here "local" means "submit from this host", which collides with the planned *local execution backend* (no scheduler). The two orthogonal axes, *where qxub runs* (local/SSH) and *what runs the job* (backend), are conflated.
*Suggest:* rename to `Transport.{IN_PROCESS, SSH, SSH_DELEGATED}` and keep `Backend` (pbs/slurm/local/k8s) as a separate platform field.

**F6. No `JobSpec` object. Med.**
Sites: `execution/context.py:82-167`, `execution/submit.py:198-220`.
Job intent is spread across `ctx_obj` (a Click dict) plus pre-rendered PBS strings (`resources` tuple of `"k=v"`, `options` string). `submit.py` even builds a fake Click-style `params` dict to reuse `process_job_options`. Backends need structured intent, not PBS strings.
*Suggest:* add a `JobSpec` dataclass (command, env context, cpus, mem, walltime, gpus, scratch, queue, account, name, stdout/stderr, workdir, env vars, hooks, tags, run_id). The CLI and `submit_job` both build it, and `Scheduler`/`JobScriptRenderer` consume it.

## G. Already well abstracted

- **Environment context is separate from scheduler concerns.** `ExecutionContext` (`execution/context.py:39-217`) already models conda / module / singularity / default without reference to PBS. Moving the scheduler directives out of the template (C1/C7) leaves this layer ready to reuse.
- **Status-file protocol.** `job_started_*`, `main_started_*`, `final_exit_code_*`, `pre/main/post_exit_code` (templates; `scheduler.py:302-372, 1097-1107`) work with any scheduler that has a shared FS, once they are keyed on a qxub run ID instead of `PBS_JOBID` (C2).
- **Base64-encoded command, pre/post hooks and user vars.** These cleanly avoid shell-quoting problems and carry over to any transport.
- **Virtual ID.** The `qx-<uuid>` from `create_queue_entry` before submission (`context.py:308-321`, `submit.py:292-306`) already decouples qxub's identity from the scheduler ID. It is the natural `QXUB_RUN_ID` for C2.
- **`SubmitResult` / `submit_job`.** `submit.py:33-54, 57` is a clean programmatic API: structured result, no `sys.exit`, no monitoring. Its internals only need to swap in `JobSpec` + `Scheduler`.
- **`ExecutionMode` / platform `remote` config.** `mode.py` already has a seam for picking a transport by platform. It needs renaming (F5), not redesign.
- **Template lookup.** `templates.py` supports a custom template override and package-resource lookup, so it can easily gain a backend dimension.
- **Walltime-gated, rate-limited probe.** `scheduler.py:1150-1158, 1229-1236` is a good generic pattern for keeping polling cheap.
- **UI helpers.** `print_status`, `SimpleSpinner`, `JobSpinner` are scheduler-neutral; they just live in the wrong module.

## H. Incidental bugs found while reading (not counted above)

1. **Singularity variables never reach the template.** `qsing.pbs:11-12, 181-203` reads `$sif_file` and `$sing_options`, but Python passes `sif=` and `bind=` (`context.py:129-131`, `core.py:212-214`, `submit.py:249-251`). As read, the singularity path runs `singularity exec  <cmd>` with no image and ignores bind mounts.
2. **Modules are never unloaded.** `qmod.pbs:240-246` checks `$modules`, but the variable passed in is `$mods`, so the unload block never runs.
3. **Dead code after `os._exit`.** `core/scheduler.py:959` has `return thread` after `os._exit(0)`; it is unreachable and refers to an undefined name.
4. **Unused function.** `core/scheduler.py:748-816` (`collect_job_resources_after_completion`) is not referenced anywhere.
5. **Docstring and validator disagree.** The `mode.py` docstring examples use `remote.hostname` (`:55, :66`), but `validate_remote_config` requires `remote.host` (`:118`).

---

## Proposed seams

```python
# qxub/scheduler/base.py
@dataclass
class JobSpec:            # replaces ctx_obj dict + "k=v" resource strings + options string
    run_id: str           # qx-<uuid>, exported to job as QXUB_RUN_ID
    command: list[str]; pre: str|None; post: str|None
    env_context: ExecutionContext
    name: str; queue: str|None; account: str|None
    cpus: int|None; memory_bytes: int|None; walltime_s: int|None
    gpus: int|None; scratch_bytes: int|None
    workdir: Path; stdout: Path; stderr: Path; scheduler_log: Path|None
    env: dict[str,str]; mail: MailSpec|None; tags: list[str]
    extra_directives: list[str]   # backend-specific escape hatch (e.g. PBS storage=)

class Scheduler(ABC):
    name: str
    def render_directives(self, spec: JobSpec) -> str: ...          # "#PBS"/"#SBATCH"/"#BSUB"/"#$"/""
    def submit(self, spec: JobSpec, script: Path) -> str: ...        # returns scheduler job id
    def parse_job_id(self, submit_stdout: str) -> str: ...
    def cancel(self, job_id: str) -> bool: ...
    def status(self, job_id: str) -> JobState: ...                   # PENDING/RUNNING/... enum
    def final_state(self, job_id: str) -> ExitStatus | None: ...     # normalised exit + reason
    def job_info(self, job_id: str) -> JobInfo | None: ...           # requested/used/timings/hosts
    def classify_submit_error(self, rc: int, stderr: str) -> SubmitErrorKind: ...
    def job_id_env_var(self) -> str: ...                             # PBS_JOBID / SLURM_JOB_ID / LSB_JOBID / JOB_ID
    capabilities: SchedulerCaps   # shared_fs, supports_interactive, accounting_delay_s, ...

class JobScriptRenderer:          # header from Scheduler.render_directives + body from env-context snippets
    def render(self, spec: JobSpec, scheduler: Scheduler) -> Path: ...

class JobMonitor(ABC):            # SharedFSMonitor (status files + tail), PodLogMonitor, LocalProcessMonitor
    def wait_started(self, ...) -> None: ...
    def stream_until_done(self, ...) -> ExitStatus: ...

@dataclass
class ExitStatus:                 # ExitStatus normaliser
    code: int; signal: int|None; reason: ExitReason   # OK/FAILED/WALLTIME/MEMORY/CANCELLED/NODE_FAIL
    def to_cli_exit_code(self) -> int: ...

class SubmissionError(Exception): # replaces QsubError (alias kept)
    kind: SubmitErrorKind; backend_exit_code: int; stderr: str
```

Platform configuration would provide: `backend`, `transport` (in-process / SSH), `paths.status_dir`, `container_runtime`, `module_system`, `notify_strategy`, and default `extra_directives` (for example Gadi `storage=gdata/a56` and `umask=0027`).
