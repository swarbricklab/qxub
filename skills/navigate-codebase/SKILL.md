---
name: navigate-codebase
description: >
  Navigate the qxub codebase to find the right module quickly.
  Use this when asked to implement, modify, or debug any qxub/qx/qxet command,
  or when you need to know where a feature lives before reading or editing files.
  Also use this when adding a new module, CLI command, or exported function —
  to look up where it fits and to update modules.json to reflect the change.
---

The structured module map for this codebase lives in `modules.json` alongside
this file. Query it with `jq` before reading source files.

> Note: `jq` is available in the `qxub` conda env and the base conda env on
> this system. If unavailable, fall back to reading `modules.json` directly.

## Useful queries

```sh
# Find the module(s) for a CLI command (e.g. "exec", "config", "status")
jq '.commands[] | select(.command | contains("exec"))' modules.json

# Find a module by keyword in its description or exports
jq '.modules[] | select(.description | ascii_downcase | contains("scheduler"))' modules.json
jq '.modules[] | select(.exports[]? | contains("ResourceMapper"))' modules.json

# List all top-level CLI commands
jq -r '.commands[].command' modules.json

# List all entry-point aliases
jq -r '.commands[] | select(.alias) | "\(.alias) → \(.command)"' modules.json

# Find which module handles shortcuts
jq '.modules[] | select(.exports[]? | contains("ShortcutManager"))' modules.json
```

Run these from the skill directory, or pass the full path to `modules.json`.

## Key structural facts (not in the JSON)

- **CLI entry point**: `qxub/cli.py` defines the main Click group and
  registers all subcommands. Each subcommand delegates to its own `*_cli.py`
  module.
- **Execution flow**: `exec_cli.py` → `execution/mode.py` (local vs remote) →
  `execution/executors.py` (build job script) → `core/scheduler.py` (qsub) →
  `queue/db.py` (record virtual ID) → `core/scheduler.py` (monitor).
- **Configuration precedence**: Override > Test > Local > Project > User >
  System > Defaults. Managed by `config/manager.py`.
- **Standalone aliases**: `qx` → `qxub exec`, `qxet` → `qxub config shortcut set`,
  `qxtat` → `qxub status`, `qxi` → interactive session. Defined in
  `config/aliases.py` and registered as entry points in `setup.py`.
- **PBS interaction**: Always via `subprocess.run(['qsub', …])` /
  `subprocess.run(['qdel', …])` in `core/scheduler.py` — never direct API.
- **Template system**: Job script templates in `qxub/jobscripts/*.pbs`.
  Selected by execution context (conda, module, singularity, default).
- **Lazy imports**: `__init__.py` uses `__getattr__` to defer heavy CLI
  imports until actually invoked.
- **Platform definitions**: YAML files loaded by `platforms/loader.py`.
  Auto-detection uses hostname matching.

## Keeping `modules.json` up to date

When you add a new module, CLI command, or significant exported function,
add or update the relevant entry in `modules.json` as part of the same change.
