# Audit B — PBS coupling in `resources/`, `platforms/`, `config/`

Scope: `qxub/resources/*.py`, `qxub/platforms/*.py`, `qxub/config/*.py`, `docs/platforms/nci_gadi.yaml`, `tests/test_platforms/*.yaml`, and the design docs `docs/dev/{multi-platform-architecture,platform_schema,config_schema}.md`. Imports were followed out of scope only where needed (`exec_cli.py`, `config_cli.py`, `execution/mode.py`). Read-only audit; no project files modified. Line numbers were checked against the files on disk on 2026-11-16.

Severity key: **High** = blocks any non-PBS backend · **Med** = wrong behaviour, or needs a refactor before a second backend can land · **Low** = naming, docs, cosmetic.

**Totals: 8 High, 21 Med, 8 Low** (plus 4 incidental bugs listed separately).

---

## Executive summary

The *platform model* (Queue / QueueLimits / auto-selection rules / internet flag / gpu_type) uses mostly neutral names and is a reasonable starting point. Everything around it is PBS-shaped, though:

1. Resources are stored as PBS `-l key=value` **strings** from the moment the CLI is parsed. `ResourceMapper` creates them, the config handler pipeline string-matches `ncpus=`, `ngpus=` and `jobfs=` to adjust them, and `build_qsub_options` renders `-N -q -P -l -o -M -m`. No typed, scheduler-neutral resource object exists anywhere.
2. `Platform.type` (`pbs_pro`) is parsed but only **printed**. Nothing dispatches on it, there is no backend registry, and the three platform YAMLs in the repo use three different spellings (`pbs_pro`, `pbs`, `scheduler: pbs`).
3. Platform detection is hard-coded: `"gadi" in $HOSTNAME → nci_gadi`. It falls back to "whichever platform loaded first".
4. Usage tracking depends on qstat `-fx` field names and on the **NCI Gadi joblog epilogue** text format, which is Gadi-specific rather than generic PBS. The DB schema bakes in `jobfs_*` and `pbs_job_id`.
5. Defaults and templates hard-code NCI facts: queue `normal`, project `a56`, `/scratch/{project}`, `gdata/{project}`, `/g/data` binds, and `copyq`.

---

## 1. Resource model and rendering

| # | Location | Sev | Finding | Suggested abstraction |
|---|---|---|---|---|
| A1 | `qxub/resources/mappers.py:13-143` (call sites: `exec_cli.py:779-801`, `config_cli.py:1066-1105`, `interactive_cli.py:1246`) | **High** | `ResourceMapper` is the only "convenience layer". It emits PBS strings directly (`mem=`, `walltime=H:MM:SS`, `ncpus=`, `ngpus=`, `jobfs=`, `storage=`) through `get_pbs_resources()` / `get_pbs_string()`. The CLI converts `--mem/--cpus/--runtime/--disk/--gpus` into these strings during argument parsing, before the target platform is known. | Replace it with a `ResourceSpec` builder that has typed fields (`cpus`, `memory_bytes`, `walltime_s`, `gpus`, `gpu_type`, `local_disk_bytes`, `extras`). Render the spec later with a per-backend `ResourceRenderer`. |
| A2 | `qxub/config/handler.py:304-317` `build_qsub_options` | **High** | This is the PBS submit-argument renderer (`-N {name} -q {queue} -P {project} -l … -o {joblog} -M -m`), and it sits in the *config* handler. `process_job_options` (`handler.py:587-609`) always returns its output as `options`. | `SchedulerBackend.render_submit_args(JobSpec) -> list[str]`. A PBS backend keeps this exact output. A Slurm backend emits `--job-name/--partition/--account/--time/--mem/--cpus-per-task/--output/--mail-*`. |
| A3 | `qxub/config/handler.py:198-218, 330-335, 359-367, 436-515, 557-578` | **High** | The whole option pipeline treats `params["resources"]` as a list of PBS `-l` strings. Code finds and rewrites resources by prefix match (`r.startswith("ncpus=")`, `"ngpus="`, `"jobfs="`) and by `key == "mem" / "walltime" / "ncpus" / "ngpus"`. Resources in any other syntax are invisible to auto-selection and minimum enforcement. | Parse the user's native `-l` input into `ResourceSpec` once, at the edge, using the PBS parser. Every later stage works on `ResourceSpec`. Render once, at submit time. |
| A4 | `qxub/resources/mappers.py:120-131`, `:191-195`; `qxub/config/manager.py:560` (`volumes: gdata/{project}`) | Med | `add_storage` / `volumes` pass NCI's `storage=gdata/xx+scratch/xx` through unchanged. This is a Gadi-only PBS resource. | Put it in a backend-namespaced `ResourceSpec.extras["pbs"]["storage"]`. Alternatively, add a `Platform.mounts` capability that only Gadi's renderer understands. |
| A5 | `qxub/resources/mappers.py:109-118` `add_disk` → `jobfs=` | Med | Generic "disk" is rendered as Gadi's `jobfs`. Slurm uses `--tmp`, Kubernetes uses `ephemeral-storage`, SGE and LSF use site-defined names. Note that `QueueLimits` already uses the neutral name `max/min_local_storage`, so the two names are inconsistent. | `ResourceSpec.local_disk_bytes`. The renderer maps it to `jobfs` / `--tmp` / `ephemeral-storage`. Each platform can name its own key in YAML (`resource_names: {local_disk: jobfs}`). |
| A6 | `qxub/resources/mappers.py:197-215` | Med | `map_snakemake_resources` forwards *every* unknown Snakemake key as a PBS custom resource (`key=value`). On another scheduler this either produces invalid arguments or silently drops the key. | Put unknown keys in `extras`. Each backend decides whether to pass, warn or drop them. *(Gone: `map_snakemake_resources` had no callers and was deleted in the dead-code sweep.)* |
| A7 | `qxub/resources/utils.py:15-55`; `qxub/resources/parser.py:15-62`; `qxub/resources/mappers.py:31, 115, 289-299` | Med | There are three inconsistent memory parsers. `utils.parse_memory_size` treats a bare number as **MB** (PBS treats it as bytes; Slurm treats it as MB). `parser.size_to_bytes` accepts only PBS lowercase `kb/mb/gb` and rejects bare numbers. The mapper regex `^\d+[KMGT]?B?$` rejects `1.5GB` and `4GiB`. `standardize_memory_unit` folds `GiB→GB`, which loses the distinction Kubernetes needs (`4G` ≠ `4Gi`). Rendered values like `4GB` are invalid for Slurm, which expects `4G`. | Add one `Quantity` / `parse_bytes()` that returns an int of bytes, with explicit SI/IEC handling. Each backend formats bytes its own way: PBS `4gb`, Slurm `4096M`, Kubernetes `4Gi`. |
| A8 | `qxub/resources/utils.py:72-127`; `qxub/resources/mappers.py:59-63`; `qxub/resources/parser.py:98-151` | Med | Walltime assumes PBS `HH:MM:SS`. A two-part `A:B` value is parsed as HH:MM (`utils.py:99-101`), but Slurm reads `30:00` as MM:SS. `D-HH:MM:SS` (Slurm/SGE) is not parsed: `parse_walltime("7-00:00:00")` returns `None`, and validation then **silently skips** the check (`core.py:168-170`). `_parse_unit_format` does not accept `1h30m` or `2d`, although the mapper does. The mapper's HH:MM:SS regex allows only 1–2 hour digits, so `100:00:00` raises. | Add a `Duration` type holding integer seconds, with one tolerant parser (`D-HH:MM:SS`, `HH:MM:SS`, `1h30m`, `2d`). Formatting belongs to the renderer. Treat a parse failure on a limit as a validation error, never a silent pass. |
| A9 | `qxub/resources/__init__.py:1-38`; `mappers.py:1-7, 14` | Low | The package is documented as "for PBS job scheduling". Public API names include `get_pbs_resources`. | Rename to `to_native()` and keep a deprecated alias. |
| A10 | `qxub/resources/*`, `qxub/platforms/core.py:34-47` | Med | The model has no multi-node or task concept. It has only total `cpus` and total `memory`. Slurm needs `--nodes/--ntasks/--cpus-per-task/--mem` vs `--mem-per-cpu`. PBS needs `select=N:ncpus=…:mem=…`, Kubernetes needs per-pod requests, and LSF needs `-n` / `span[]`. Today users must hand-write `-l select=…`, which is not portable. | Give `ResourceSpec` the fields `nodes`, `tasks_per_node`, `cpus_per_task`, and `mem_per_node` or `mem_total`. Add `QueueLimits.cores_per_node` and `mem_per_node` so renderers can derive chunks. |

## 2. Queue model and auto-selection

| # | Location | Sev | Finding | Suggested abstraction |
|---|---|---|---|---|
| B1 | `qxub/platforms/core.py:194-204, 352-386`; `qxub/platforms/cli.py:52, 135` | **High** | `Platform.type` is a required key (`data["type"]` raises KeyError when it is missing) but is **never dispatched on**. Its only use is `click.echo(f"Type: …")`. There is no backend registry, enum or validation. The repo's YAMLs disagree: `nci_gadi.yaml:7` `pbs_pro`, `test_platform.yaml:2` `pbs`, `localhost.yaml:6` `scheduler: "pbs"`. No extension point exists for a second scheduler. | Add a `Platform.scheduler: SchedulerKind` enum (`pbs_pro`, `torque`, `slurm`, `sge`, `lsf`, `local`, `k8s`) that resolves through a `BACKENDS` registry or entry points to a `SchedulerBackend`. Accept `type` as a deprecated alias. |
| B2 | `qxub/platforms/core.py:34-47, 50-72, 103-108, 163-173, 393-406`; `nci_gadi.yaml:18` etc. | Med | Walltime limits exist **only** as Gadi-style "max walltime by core-count band" (`WalltimeRule`). `QueueLimits` has no `max_walltime`, so the `limits.max_walltime` present on every queue in `nci_gadi.yaml`, and `limits.walltime` in `test_platform.yaml`, is **silently ignored**. When no band matches, walltime is not validated at all. Slurm partitions have a flat `MaxTime` plus QOS `MaxWall`. | Add `QueueLimits.max_walltime` as the primary limit and keep `walltime_rules` as an optional refinement. Also add `QueueLimits.max_nodes`. |
| B3 | `qxub/platforms/core.py:85-101` | Med | `Queue` is the only scheduling target. Slurm separates **partition**, **QOS** and **account**; LSF has queues plus application profiles; Kubernetes has neither (it uses namespace, nodeSelector and priorityClass). `Queue.type` is free text (`standard`, `gpu`, `data`, `compute`) and no code uses it. | Generalise to a `Target` (or keep the name `Queue`) with a backend-specific `selectors` dict (`{partition: gpu, qos: long}` for Slurm, `{queue: gpuvolta}` for PBS, `{nodeSelector: …}` for Kubernetes). Selection logic can stay as it is. |
| B4 | `qxub/config/handler.py:167-301` `select_auto_queue` | **High** | Auto-selection loops over **every loaded platform** (`:249-281`) and takes the "cheapest" queue name from any of them, ignoring which platform or scheduler the job will be submitted to. Every failure path falls back to the literal `"normal"` (`:194, 288, 294, 299`). Requirements are rebuilt from PBS `-l` strings (`:198-218`). | Resolve the target `Platform` once, then call `platform.select_queue(spec)`. Fall back to `platform.default_queue`, or `None` so the backend omits the flag, never `"normal"`. |
| B5 | `qxub/config/handler.py:258-276`; `qxub/platforms/core.py:93, 96, 110-112` | Med | Cost comparison calls `estimate_su_cost`, which uses `su_rate` (default 1.0), not `su_billing_rate`. `nci_gadi.yaml` sets only `su_billing_rate`, so every queue costs the same and the first platform wins. `walltime_hours` defaults to `3600` and is then treated as seconds. There are two overlapping rate fields. | Merge them into one `AccountingModel` (see D-section) with a single rate field. |
| B6 | `qxub/config/handler.py:320-381` `adjust_cpus_for_gpus`, `:384-520` `enforce_queue_minimums` | Med | Both look up the queue **by name across all platforms, first match wins** (`:352-357`, `:416-422`). A name like `normal` can collide between a PBS and a Slurm platform. They rewrite `ncpus=` and `ngpus=` strings and inject `jobfs=` minimums, which is Gadi semantics (`min_local_storage`→`jobfs`). | Operate on `(Platform, ResourceSpec)`. Minimum enforcement becomes a backend-neutral `QueuePolicy.apply(spec)`; only rendering is backend-specific. |
| B7 | `qxub/config/handler.py:523-584` `adjust_resources_for_queue` | Med | Uses `get_current_platform()` (hostname detection), while the sibling functions use "all platforms via `QXUB_PLATFORM_PATHS`". One pipeline therefore resolves the platform three different ways. The docstring says "Mimics PBS qsub behavior". | Resolve the platform once in `process_job_options` and pass it to every step. |
| B8 | `qxub/platforms/core.py:110-112, 551-559, 627-631, 682-697, 773-777`; `qxub/platforms/cli.py:169, 253, 380-445`; `qxub/platforms/__init__.py:8` | Med | Cost is hard-wired to NCI Service Units: the strings "SU", "SU/CPU·hour" and "Service Unit", and the formula `cpus × hours × rate`. The formula is wrong even on Gadi, because NCI charges on `max(ncpus, mem / mem_per_core)`. Other sites bill in core-hours, GPU-hours, fairshare or dollars, and some do not bill at all. | Add `AccountingModel` (`unit_name`, `rate_per_core_hour`, `rate_per_gpu_hour`, `mem_per_core` for memory-weighted charging, or `none`). `estimate(spec, queue) -> Cost(amount, unit)`. |
| B9 | `qxub/platforms/core.py:98, 200-204, 375-376`; `nci_gadi.yaml:294-304`; `platform_schema.md:165-199` | Med | `auto_adjust` is parsed into `Platform.auto_adjust` but never read. The YAML gives it an expression DSL (`round_cores_to_multiple(cores, 48)`), while the schema doc gives it policy values (`auto/suggest/user/error`). `Queue.constraints` is never evaluated. Gadi's "whole-node rounding" rules therefore exist only as dead data. | Choose one schema. Add a small `CoreRounding` policy (`multiple_of`, `allowed_values`), applied by the neutral `QueuePolicy`. |

## 3. Platform detection and loading

| # | Location | Sev | Finding | Suggested abstraction |
|---|---|---|---|---|
| C1 | `qxub/platforms/core.py:473-502` | **High** | `detect_platform()` hard-codes `if "gadi" in $HOSTNAME → "nci_gadi"`. It ignores the YAML `host` field (the comment at `nci_gadi.yaml:8` claims that field drives matching), and `$HOSTNAME` is often not exported to Python. `get_current_platform()` then falls back to `list_platforms()[0]`, which follows glob order, so a Slurm host with several definitions installed can silently get the Gadi platform. `interactive_cli.py:115` and `handler.py:549` depend on it. | Add a `detect:` block per platform YAML (hostname regex, env var such as `SLURM_CLUSTER_NAME` or `PBS_SERVER`, `which sbatch`). Precedence: CLI `--platform` > config `defaults.platform` > detection > error. Never fall back to the first platform. |
| C2 | `qxub/platforms/core.py:273-448` vs `qxub/platforms/loader.py:32-359`; search paths `config/manager.py:710-714`, `platforms/loader.py:198-204`, `config/base.py:38-42` | Med | There are two `PlatformLoader` classes. The core one scans directories and builds objects. The `loader.py` one handles `file://` and `https://` with caching but returns raw dicts and **is not used outside the package**, so a config `platforms.<name>.definition: https://…` URL never reaches the object model. The default search paths differ across the three places: `/etc/qxub/platforms` vs `/etc/xdg/qxub/platforms`. This is not PBS-specific, but it blocks a clean registry. | Use one `PlatformRegistry`: discover with the URL-capable loader, construct with `_create_platform`, and key by name and scheduler. *(The unused `platforms/loader.py` was deleted in the dead-code sweep; it is in git history if a URL-capable registry is built. `qxub.platforms.PlatformLoader` is now the core loader.)* |
| C3 | `tests/test_platforms/localhost.yaml:1-21`; `tests/test_platforms/test_platform.yaml:10-13, 22, 31, 41` | Low | The YAMLs follow different schemas. `localhost.yaml` uses `scheduler:`, a mapping-style `queues`, and top-level `defaults` with `project`, `ncpus` and `mem`; `_create_platform` would raise KeyError on `type` and `host`. `test_platform.yaml` uses `limits.walltime` and `cost_per_hour`, which are silently ignored. The schema is never validated. | Add a versioned JSON-schema or pydantic model with `schema_version`. Reject unknown keys, or warn on them. |

## 4. Accounting, usage tracking and job identity

| # | Location | Sev | Finding | Suggested abstraction |
|---|---|---|---|---|
| D1 | `qxub/config/handler.py:88-89, 306`; `config/manager.py:222-241, 550, 698-700` | **High** | "Project" is rendered as PBS `-P` and defaults to `$PROJECT`, an NCI-set environment variable. The other schedulers differ: Slurm uses `--account`, SGE uses `-P` (a project, with different semantics) or `-A`, LSF uses `-P`, and Kubernetes uses namespace or labels. Some sites need no account. | Add `JobSpec.account`, rendered by the backend, with a configurable env-var source per platform (`account_env: PROJECT` / `SLURM_ACCOUNT`). Keep the `{project}` template variable as an alias of `{account}`. |
| D2 | `qxub/resources/parser.py:288-416`; `qxub/resources/tracker.py:785-825, 827-977` | Med | `parse_joblog_resources` parses NCI Gadi's **custom epilogue** ("Resource Usage on …", "Service Units", "JobFS requested"). It is not generic PBS Pro, Torque output differs, and other schedulers have nothing equivalent. `finalize_job` and `backfill_resources` rely on PBS appending that block to the joblog after the script exits. | Add `UsageCollector` per backend: PBS uses `qstat -fx` or the Gadi epilogue, Slurm uses `sacct --json`, SGE uses `qacct -j`, LSF uses `bhist -l`/`bacct`, Kubernetes uses the metrics API, and local uses `resource.getrusage`. Each returns a neutral `JobUsage`. *(`finalize_job` and `backfill_resources` had no callers and were deleted in the dead-code sweep.)* |
| D3 | `qxub/resources/tracker.py:234-307` `log_job_resources` | **High** | It consumes the qstat `-fx`-derived dict: `resources_requested.{mem,walltime,ncpus,jobfs}`, `resources_used.{mem,walltime,ncpus,jobfs,cpupercent}`, `execution.exec_host_parsed`. It converts values with PBS-only parsers (`size_to_bytes` with lowercase units, `time_to_seconds` that requires `H:M:S`). Slurm `D-HH:MM:SS` or `MaxRSS=1234K` values fail inside a try/except and are stored as NULL. | `log_job_usage(job_id, JobUsage)`, where `JobUsage` holds typed ints (bytes, seconds, cpus, gpus) and the backend fills it. |
| D4 | `qxub/resources/parser.py:154-261`; exported in `resources/__init__.py:42-52, 78-87` | Med | `parse_exec_host` (`host/core`), `parse_exec_vnode` (`(host:ncpus=…:mem=…kb:jobfs=…)`) and `parse_timestamp` (ctime) are PBS formats exposed as the package's public API. | Move them to `backends/pbs/parsing.py`. Keep only neutral helpers in `resources`. |
| D5 | `qxub/resources/tracker.py:121-169` (schema), `:428, 694-704, 747-759, 810-816, 1001-1004` | Med | DB schema and queries are PBS-named: `jobfs_requested_mb`, `jobfs_used_mb`, `jobfs_efficiency`, "Path to PBS joblog", and joins on `queue.pbs_job_id`. `cpu_percent` is PBS's `cpupercent`. | Rename to `scheduler`, `scheduler_job_id` and `local_disk_*` with a migration, and keep compatibility views. |
| D6 | `qxub/resources/tracker.py:76-80, 449-452` | Low | Comments mention `/g/data`, Gadi `/scratch`, and `/g/data/a56/software/...`. The underlying concern (mmap and journal mode on Lustre or GPFS) is generic and correct. | Make the wording neutral. |
| D7 | `qxub/resources/parser.py:1-8, 159-200, 419-443` | Low | Docstrings, examples and the `__main__` block use Gadi hostnames and a PBS-only vocabulary. | Move them with D4. |

## 5. Config keys, defaults and NCI paths

| # | Location | Sev | Finding | Suggested abstraction |
|---|---|---|---|---|
| E1 | `qxub/config/handler.py:12-28` | Med | The default output directory is `/scratch/$PROJECT/$USER/qt/<ts>`, with `PROJECT` defaulting to the literal **`"a56"`**. | Add a platform-provided `paths.scratch` template, falling back to `$XDG_CACHE_HOME/qxub` or the current directory. |
| E2 | `qxub/config/handler.py:80-81` | Med | With no configuration, the queue defaults to `"normal"`. Slurm should omit `--partition` and use the cluster default; Kubernetes has no queues. | Default to `None`, meaning the backend omits the flag, or use `platform.default_queue`. |
| E3 | `qxub/config/manager.py:533-603` `create_user_config_template` | Med | The generated user config is NCI-specific: queue `normal`, project `a56`, `resources: ["mem=4GB","ncpus=1"]`, `/scratch/{project}/…`, `volumes: gdata/{project}`, `sing.bind: [/scratch, /g/data]`, `"Direct PBS execution"`. | Generate the template from the detected platform's YAML, using a `defaults:` and `paths:` block per platform. |
| E4 | `qxub/config/manager.py:516, 523` `init_project_config` | Low | `.qx/test.yaml` defaults to queue `copyq`, a Gadi queue name. | Use `queue: auto` or omit the key. |
| E5 | `qxub/config/handler.py:31-53, 137-139` | Low | `_sanitize_job_name` applies PBS rules ("PBS compliance") and does not enforce PBS's own 15- or 236-character limit. Kubernetes requires lowercase DNS-1123 names of at most 63 characters; Slurm is permissive. | `backend.sanitize_job_name(name)`. |
| E6 | `qxub/config/handler.py:93-100, 311-315` | Low | `email_opts` defaults to `"ae"`, which are PBS `-m` letters. Slurm needs `--mail-type=FAIL,END`. | Use a neutral `NotifyOn = {begin, end, fail}` set and let the backend render it. |
| E7 | `qxub/config/manager.py:222-241`; `qxub/config/handler.py:105-109` | Low | The template variables `{project}` and `{queue}` have PBS-flavoured names. | Add `{account}` and `{partition}` aliases. |
| E8 | `qxub/config/shortcuts.py` (data), `config_cli.py:1066-1105`, `config/manager.py:552` | Med | `ShortcutManager` itself is scheduler-agnostic, but shortcuts and aliases **persist PBS `-l` strings** (`"resources": ["mem=4GB", …]`) because `--mem/--cpus` are converted before saving. Saved shortcuts therefore cannot be reused on another platform. | Store neutral fields (`mem`, `cpus`, `runtime`, `gpus`, `disk`) and allow `native: {pbs: [...], slurm: [...]}` for escape hatches. |

## 6. Incidental bugs found while auditing (not PBS coupling, but they block refactoring)

- `qxub/platforms/core.py:824`: `from .config import config_manager` resolves to the non-existent `qxub.platforms.config` (checked with `importlib.util.find_spec` → `None`). `select_best_queue` raises `ModuleNotFoundError` on every call, so `qxub select-queue` (registered in `qxub/cli.py:82`) is broken. The import should be `from ..config import config_manager`.
- `qxub/platforms/core.py:445-448`: `reload_platforms()` calls the non-existent method `self._load_all_platforms()`.
- `qxub/platforms/core.py:263-270`: `validate_queue_resources` is annotated `-> List[str]` but returns either `ResourceValidationResult` or `list`.
- `qxub/config/base.py:300`: `from .platforms import get_platform` resolves to the non-existent `qxub.config.platforms`. The module is deprecated, but it is still imported by `config/__init__.py:51-55`, which also triggers a `DeprecationWarning` on every import.
- (Minor) `qxub/resources/tracker.py:614-637`: `get_resource_trends` calls `dict(row)` without setting `row_factory`, which raises on tuples.

---

## 7. Already well abstracted (build on these)

- **`evaluate_condition`** (`resources/utils.py:178-267`) uses neutral variables (`cpus`, `memory`, `walltime`, `gpu_requested`, `gpu_type`). Auto-selection rules written in it are portable, apart from the walltime and memory parser issues (A7, A8).
- **`QueueLimits`** names (`max_cpus`, `max_memory`, `max/min_local_storage`, `max/min_gpus`) are neutral, as are `Queue.gpu_type`, `cpus_per_gpu` and `internet_connectivity`, and the scoring in `QueueSelector`. All of these carry over to Slurm partitions.
- **`suggest_resource_adjustment`** and `ResourceValidationResult` hold neutral data (they need the walltime and memory fixes).
- **`platforms/loader.py`** (URL, file and search-path loading with ETag and TTL caching) is fully scheduler-agnostic. It is just not wired into the object model (C2). *(Deleted as unused; see C2.)*
- **`ConfigManager`** (system/user/project/local/test/override hierarchy, OmegaConf merge, template resolution) and **`ShortcutManager`** are agnostic. Only their *defaults and data* are NCI-shaped (E3, E8).
- **`ExecutionMode`** (`execution/mode.py`) separates local, remote and remote-delegated execution cleanly, which is the right place for a `local` backend.
- **Tracker storage**: values are normalised to MB, seconds and percent, and the sqlite settings for Lustre and GPFS (no mmap, `journal_mode=DELETE`) are generic and sound. Only the column names and input parsers are PBS-specific (D3, D5).
- `config/base.py:164-168` `default_resources` already uses neutral keys (`cpus`, `memory`, `walltime`).

---

## 8. Design-vs-code drift

| Design doc says | Code does |
|---|---|
| `multi-platform-architecture.md:42-101`: a `qxub/platforms/base.py` containing a `PlatformResources` ABC (`to_command_args`, `validate`, `estimate_cost`), plus `PBSResources` and `SlurmResources` | None of these exist. No `base.py`, `pbs_pro.py` or `slurm.py`. Resources are plain `-l` strings, and rendering is `build_qsub_options` in the config handler. |
| `:103-162`: a `workflow_adapters/` package with a `WorkflowAdapter.translate_resources(platform, …)` that returns platform-native resources | Workflow mapping lives in `resources/mappers.py` (`map_snakemake/nextflow/cwl`) and always emits PBS. It takes no `platform` argument. |
| `:213-225` Phase 2: `--mem/--cpus/--time/--gpus` convenience flags "auto-translate to platform syntax" | The flags exist (exec_cli), but translate unconditionally to PBS through `ResourceMapper` (A1). |
| `:19-23` "Config-driven platform detection" and "Validation at platform level" | Detection is a hard-coded hostname check (C1). Validation happens in `Queue.validate_resources`, which is fine, but the pipeline resolves the platform inconsistently (B4, B6, B7). |
| `:366-393` platform YAML uses `scheduler: slurm`, `scheduler_config`, `resource_mappings`, `max_walltime: "7-00:00:00"` | The loader reads `type`, not `scheduler`. It ignores `scheduler_config` and `resource_mappings`. `QueueLimits` has no `max_walltime`, and `parse_walltime("7-00:00:00")` returns `None` (A8, B2). |
| `:343-355` config `platform.default`, `remote_platforms.*.scheduler` | The config uses `default_platform` (`manager.py:716-718`), `defaults.platform` (`manager.py:777-785`), and a `platforms.<name>` registry with `definition` and `remote`. There are two different "default platform" keys, and `select_best_queue` reads the one that `exec_cli` does not. No `scheduler` key exists anywhere. |
| `platform_schema.md:7, 51` "loaded via `qxub/platform.py`" | The code lives in `qxub/platforms/core.py`. |
| `platform_schema.md:16-47, 125-152` `queues` is a **mapping** keyed by name; `auto_selection.rules` uses `if/then/default` and `priority_order`; `resources.auto_min_cpus` is nested | The loader expects `queues` as a **list** with `name`, and `auto_selection_rules` with `condition`, `queue` and `is_default`. `auto_min_cpus` is read at the top level of the queue. `priority_order` is not implemented. |
| `platform_schema.md:59-111` (an orphaned fragment) uses `su_rate` | `nci_gadi.yaml` uses `su_billing_rate`. The code has **both** fields with different defaults, and `handler.py` uses the one the YAML never sets (B5). |
| `platform_schema.md:165-199` `auto_adjust: {min_cpus: auto, memory: suggest, …}` policies | `nci_gadi.yaml:295-304` uses an expression DSL, and the code ignores `auto_adjust` entirely (B9). |
| `platform_schema.md:215, 219-220` core ranges `"1-672,720-1440"`; conditions with `AND`; variables `queue` and `su_cost` | `WalltimeRule.matches_core_count` cannot combine `-` and `,`: it splits on `-` into 3 parts and returns `False`, so the rule silently never matches. `evaluate_condition` supports a single comparison only, and only the variables `gpu_requested`, `gpu_type`, `cpus`, `memory` and `walltime`. |
| `platform_schema.md:222-233` validation rules (exactly one default, no walltime gaps, non-overlapping ranges) | No schema validation exists (C3). |
| `platform_schema.md:302-323` extension points (`custom_validator`, `platform_resources`) | Not implemented. `constraints` is parsed and never evaluated. |
| `config_schema.md:7-11` project config at `.qxub/config.yaml` and `.qxub/local.yaml` | The code uses `.qx/project.yaml`, `.qx/local.yaml` and `.qx/test.yaml` (`manager.py:52-102`). |
| `config_schema.md:56` "See `qxub/config_manager.py`" | The code is in `qxub/config/manager.py`. |
| `config_schema.md:101-145` `auto_selection.*` and `queue_overrides.*` sections | Not read anywhere. The code reads `queue_preferences` (`optimization`, `adjustment_policy`, `auto_select`), which is undocumented. |
| `config_schema.md:147-166` `profiles.<name>` with `platform`, `defaults` and `remote` | The code uses `platforms.<name>` with `definition` and `remote` (`manager.py:737-764`). `profiles` is not read. |
| `config_schema.md:168-206` internal `resource_request` with `core_resources` (neutral) and `platform_specific` | This is **the right design, but not implemented**: there is no internal neutral representation. |
| `config_schema.md:15-35` defaults use `ncpus`, `mem` and `walltime` (PBS keys) | The code's convenience defaults use `cpus`, `mem`, `runtime` and `disk` (`manager.py:556-559`). `walltime` and `ncpus` keys are only meaningful inside `resources: [...]` strings. |

The docs describe the intended split: neutral `core_resources` versus `platform_specific`, and a `PlatformResources` ABC per scheduler. The code has none of those seams. The only "platform" layer that exists is queue metadata for selection and validation.

---

## 9. Proposed seams for this area

```python
# qxub/resources/spec.py  — scheduler-neutral, typed
@dataclass
class ResourceSpec:
    cpus: int | None = None              # total, or per task when tasks set
    memory_bytes: int | None = None      # total per job (or per node, see mem_scope)
    mem_scope: Literal["job","node","cpu"] = "job"
    walltime_s: int | None = None
    gpus: int = 0
    gpu_type: str | None = None
    local_disk_bytes: int | None = None  # jobfs / --tmp / ephemeral-storage
    nodes: int | None = None
    tasks_per_node: int | None = None
    internet: bool = False
    extras: dict[str, dict[str, str]] = field(default_factory=dict)  # {"pbs": {"storage": "gdata/a56"}}

# Single tolerant parsers (replace utils/parser/mappers trio)
def parse_bytes(s: str, *, bare_unit: str) -> int: ...      # SI vs IEC explicit
def parse_duration(s: str) -> int: ...                       # D-HH:MM:SS, HH:MM:SS, 1h30m, 2d

@dataclass
class JobSpec:
    name: str
    command: ...
    resources: ResourceSpec
    target: "QueueRef | None"            # queue / partition+qos / none
    account: str | None                  # was "project" (-P)
    stdout: Path; stderr: Path; joblog: Path | None
    notify: set[Literal["begin","end","fail"]]; email: str | None

class ResourceRenderer(Protocol):        # one per backend
    def render(self, spec: ResourceSpec, target: QueueRef | None) -> list[str]: ...
    def parse_native(self, args: list[str]) -> ResourceSpec: ...   # PBS: "-l mem=4GB" -> spec (keeps -l back-compat)

class SchedulerBackend(Protocol):
    kind: SchedulerKind                  # pbs_pro | torque | slurm | sge | lsf | local | k8s
    renderer: ResourceRenderer
    usage: "UsageCollector"
    def render_submit_args(self, job: JobSpec) -> list[str]: ...   # replaces build_qsub_options
    def sanitize_job_name(self, name: str) -> str: ...
    def default_paths(self, account: str | None) -> dict[str, Path]: ...

BACKENDS: dict[SchedulerKind, type[SchedulerBackend]]  # entry-point extensible

@dataclass
class Platform:
    name: str
    scheduler: SchedulerKind             # dispatches to BACKENDS (accept legacy `type`)
    detect: "DetectRule"                 # hostname regex / env var / binary presence
    queues: dict[str, Queue]             # Queue.selectors: backend-specific (partition/qos/queue/nodeSelector)
    default_queue: str | None
    accounting: "AccountingModel"
    policy: "QueuePolicy"                # core rounding, mins, cpus_per_gpu, local-disk minimums
    paths: dict[str, str]                # scratch/output templates (no /scratch, /g/data in code)
    resource_names: dict[str, str]       # e.g. {"local_disk": "jobfs"} for Gadi
    def backend(self) -> SchedulerBackend: ...

@dataclass
class QueueLimits:                       # add: max_walltime, max_nodes, cores_per_node, mem_per_node
    ...

class QueuePolicy:                       # replaces enforce_queue_minimums / adjust_cpus_for_gpus / adjust_resources_for_queue
    def apply(self, spec: ResourceSpec, queue: Queue) -> tuple[ResourceSpec, list[str]]: ...

class AccountingModel:                   # replaces su_rate/su_billing_rate/"SU" strings
    unit: str | None                     # "SU", "core-h", "$", None
    rate_per_core_hour: float | None
    rate_per_gpu_hour: float | None
    mem_per_core_bytes: int | None       # NCI-style memory-weighted charging
    def estimate(self, spec: ResourceSpec, queue: Queue) -> "Cost | None": ...

@dataclass
class JobUsage:                          # neutral input to ResourceTracker
    scheduler_job_id: str; exit_code: int | None
    cpus_req: int | None; cpus_used: float | None
    mem_req_bytes: int | None; mem_used_bytes: int | None
    walltime_req_s: int | None; walltime_used_s: int | None
    local_disk_req_bytes: int | None; local_disk_used_bytes: int | None
    cost: "Cost | None"; host: str | None

class UsageCollector(Protocol):          # PBS: qstat -fx / Gadi joblog; Slurm: sacct; SGE: qacct; LSF: bhist; k8s: metrics
    def collect(self, job_id: str, joblog: Path | None) -> JobUsage | None: ...

class PlatformRegistry:                  # merges core.PlatformLoader + loader.PlatformLoader
    def resolve(self, cli_platform: str | None) -> Platform: ...  # CLI > config > detect > error
```

Suggested order of work:
1. Fix the incidental bugs (section 6).
2. Introduce `ResourceSpec`, `parse_bytes` and `parse_duration`, plus a PBS `parse_native` and `render`, and show that the PBS output is byte-identical.
3. Thread `(Platform, ResourceSpec)` through `process_job_options`, replacing A3, B4, B6 and B7.
4. Add `Platform.scheduler` and the `BACKENDS` registry, and move `build_qsub_options` into `PbsBackend`.
5. Add `AccountingModel` and `UsageCollector`, and rename the DB columns.
6. Move NCI defaults and paths into `nci_gadi.yaml`.
7. Add a Slurm backend as the second implementation to validate the seams.
