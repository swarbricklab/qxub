"""
New dual-log history system for qxub.

This system separates "what" (recipes) from "how/when" (executions):
- history/recipes/<hash>.json: Unique computational recipes indexed by hash
- history/executions/<timestamp>.json: Individual execution records indexed
  by microsecond timestamp

Each record is its own file so concurrent qxub processes cannot lose each
other's records (#129); see ``store.py``. The older ``recipes.yaml`` and
``executions.yaml`` are imported on first use.

The manager provides comprehensive job tracking including:
- Recipe storage and deduplication
- Execution logging with file path tracking
- Job lookup by ID for output file viewing
- Support for all execution contexts (conda, modules, containers, default)

New in v3.2.0:
- File path tracking in execution records (stdout, stderr, PBS log)
- Job ID-based execution lookup for log viewing
- Enhanced CLI commands: qxub history out/err/log [job_id]
"""

import hashlib
import json
import logging
import os
import re
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from .store import RecordStore, ensure_private_dir, make_owner_only

logger = logging.getLogger(__name__)

# Keep this many execution records; prune once there are PRUNE_SLACK more.
MAX_EXECUTIONS = 1000
PRUNE_SLACK = 50
# A record with no job ID yet may still be updated (qxi fills it in when the
# session ends), so pruning leaves it alone for this long.
AWAITING_UPDATE_SECONDS = 24 * 3600


# Values of --var/--vars often carry credentials; history keeps only names.
REDACTED = "***"
_VAR_NAME_RE = re.compile(r"(?:^|,)\s*([^=,\s]+)\s*=")


def _redact_var(item: str) -> str:
    """``NAME=value`` -> ``NAME=***``; a bare ``NAME`` is kept."""
    name, sep, _ = str(item).partition("=")
    return f"{name}={REDACTED}" if sep else name


def _redact_vars(value: str) -> str:
    """``A=1,B=2`` -> ``A=***,B=***``. Values may contain commas, so only
    the names are kept rather than splitting the string."""
    names = _VAR_NAME_RE.findall(str(value))
    return ",".join(f"{n}={REDACTED}" for n in names) if names else REDACTED


def redact_argv(args: List[str]) -> List[str]:
    """``args`` with the values of ``--var``/``--vars`` replaced.

    Arguments after ``--`` belong to the job's command and are kept.
    """
    out: List[str] = []
    redact_next = None
    for i, arg in enumerate(args):
        if redact_next is not None:
            out.append(redact_next(arg))
            redact_next = None
        elif arg == "--":
            out.extend(args[i:])
            break
        elif arg == "--var":
            out.append(arg)
            redact_next = _redact_var
        elif arg == "--vars":
            out.append(arg)
            redact_next = _redact_vars
        elif arg.startswith("--var="):
            out.append("--var=" + _redact_var(arg[len("--var=") :]))
        elif arg.startswith("--vars="):
            out.append("--vars=" + _redact_vars(arg[len("--vars=") :]))
        else:
            out.append(arg)
    return out


def redact_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """Click parameters with the values of ``var``/``vars`` replaced."""
    params = dict(params)
    if params.get("var"):
        params["var"] = [_redact_var(v) for v in params["var"]]
    if params.get("vars"):
        params["vars"] = _redact_vars(params["vars"])
    return params


def _with_imported_baseline(recipe: Dict[str, Any], until: str) -> Dict[str, Any]:
    """Keep a YAML recipe's run count and last use as a baseline.

    ``until`` is the import time as an execution key: executions up to it
    are already in ``run_count``, so ``get_recipes`` only adds later ones.
    """
    metadata = dict(recipe.get("metadata") or {})
    metadata["imported"] = {
        "run_count": metadata.pop("run_count", 0) or 0,
        "last_seen": metadata.pop("last_seen", None),
        "until": until,
    }
    return {**recipe, "metadata": metadata}


def _resources_section(resource_data: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "requested": resource_data.get("resources_requested", {}),
        "used": resource_data.get("resources_used", {}),
        "efficiency": resource_data.get("efficiency", {}),
        "execution_env": resource_data.get("execution", {}),
        "timing": resource_data.get("timing", {}),
    }


class HistoryManager:
    """Manages dual-log history system with recipes and executions."""

    def __init__(self, config_dir: Optional[Path] = None):
        """Initialize history manager with config directory."""
        if config_dir is None:
            # Use XDG config directory
            config_dir = Path.home() / ".config" / "qxub"

        self.config_dir = Path(config_dir)
        # Written by older qxub; imported on first use
        self.recipes_file = self.config_dir / "recipes.yaml"
        self.executions_file = self.config_dir / "executions.yaml"
        self._history_dir = self.config_dir / "history"
        self._recipes = RecordStore(self._history_dir / "recipes")
        self._executions = RecordStore(self._history_dir / "executions")
        self._private = False

        # Ensure config directory exists and is owner-only
        # pylint: disable-next=import-outside-toplevel
        from ..config.files import ensure_private_config_dir

        ensure_private_config_dir(self.config_dir)

    def _ensure_private(self) -> None:
        """Make the history directories owner-only, once per process.

        Records hold full command lines. Directories and records written by
        older qxub with the umask's mode are tightened here.
        """
        if self._private:
            return
        if ensure_private_dir(self._history_dir):
            for legacy in (self.recipes_file, self.executions_file):
                make_owner_only(legacy.with_name(legacy.name + ".imported"))
        for store in (self._recipes, self._executions):
            store._ensure_dir()  # pylint: disable=protected-access
        self._private = True

    def compute_recipe_hash(self, recipe_data: Dict[str, Any]) -> str:
        """Compute stable hash for a recipe based on executor and target."""
        # Only hash the executor and target - ignore metadata and context
        hashable = {
            "executor": recipe_data.get("executor", {}),
            "target": recipe_data.get("target", {}),
        }

        # Canonical JSON for stable hashing
        canonical = json.dumps(hashable, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode()).hexdigest()[:12]

    def _extract_recipe_from_command(self, ctx, command_line: str) -> Dict[str, Any]:
        """Extract recipe definition from Click context and command line."""
        import sys

        recipe = {"executor": {}, "target": {}}

        if len(sys.argv) > 1:
            args = redact_argv(sys.argv[1:])  # Skip script name

            # Handle new v3 CLI structure: qxub exec --env test "command"
            if len(args) > 0 and args[0] == "exec":
                args = args[1:]  # Remove 'exec' subcommand

            # For unified CLI (qxub 2.0/3.0), determine executor type from execution context options
            executor_type = None
            executor_params = {}

            # Check for execution context options
            i = 0
            while i < len(args):
                arg = args[i]

                if arg == "--env" and i + 1 < len(args):
                    # Conda environment specified
                    executor_type = "conda"
                    executor_params["env"] = args[i + 1]
                    i += 2
                elif arg in ["--mod", "--mods"] and i + 1 < len(args):
                    # Module specified (handle both --mod and --mods)
                    executor_type = "module"
                    executor_params["name"] = args[i + 1]
                    i += 2
                elif arg == "--sif" and i + 1 < len(args):
                    # Singularity container specified
                    executor_type = "sing"
                    executor_params["image"] = args[i + 1]
                    i += 2
                elif arg == "--":
                    # Found command separator - everything after is the target command
                    i += 1
                    if i < len(args):
                        recipe["target"] = {"cmd": " ".join(args[i:])}
                    break
                elif (
                    arg.startswith("--")
                    and i + 1 < len(args)
                    and not args[i + 1].startswith("-")
                ):
                    # Other executor options
                    option_name = arg[2:]  # Remove --
                    executor_params[option_name] = args[i + 1]
                    i += 2
                elif arg.startswith("--"):
                    # Flag without value
                    option_name = arg[2:]  # Remove --
                    executor_params[option_name] = True
                    i += 1
                elif not arg.startswith("-"):
                    # This is the start of the target command (no -- separator found)
                    recipe["target"] = {"cmd": " ".join(args[i:])}
                    break
                else:
                    i += 1

            # Set executor type and params
            if executor_type:
                executor_params["type"] = executor_type
                recipe["executor"] = executor_params

            # If no explicit executor type found, check for legacy subcommand format for backward compatibility
            if not executor_type:
                executor_types = ["conda", "module", "sing"]
                for j, arg in enumerate(args):
                    if arg in executor_types:
                        executor_type = arg
                        recipe["executor"] = {"type": executor_type}
                        # Parse remaining arguments for this legacy format
                        k = j + 1
                        while k < len(args):
                            if args[k] == "--":
                                k += 1
                                if k < len(args):
                                    recipe["target"] = {"cmd": " ".join(args[k:])}
                                break
                            elif not args[k].startswith("-"):
                                recipe["target"] = {"cmd": " ".join(args[k:])}
                                break
                            k += 1
                        break

            # If still no executor type found, default to "default" type
            if not executor_type:
                recipe["executor"] = {"type": "default"}

        return recipe

    # -- writing -------------------------------------------------------------

    def log_execution(
        self,
        ctx,
        success: bool = True,
        error: Optional[str] = None,
        job_id: Optional[str] = None,
        resource_data: Optional[Dict] = None,
        file_paths: Optional[Dict[str, str]] = None,
        tags: Optional[List[str]] = None,
    ) -> str:
        """Log a command execution and return the execution timestamp.

        The timestamp is the record's key; pass it to
        ``update_execution_with_resources`` once the job has finished.
        """
        execution_timestamp = str(int(time.time() * 1000000))
        try:
            import sys

            command_line = " ".join(redact_argv(sys.argv))
            recipe_data = self._extract_recipe_from_command(ctx, command_line)

            # Only proceed if we have a valid executor (not config commands, etc.)
            if not recipe_data.get("executor") or recipe_data["executor"].get(
                "type"
            ) not in ["conda", "module", "sing", "default"]:
                return execution_timestamp

            self._ensure_private()
            self._import_legacy_yaml()

            recipe_hash = self.compute_recipe_hash(recipe_data)
            now = datetime.now().isoformat()
            # A recipe's content depends only on its hash, so whichever
            # process writes it first, the result is the same. Run counts and
            # last-seen times are derived from the execution records.
            self._recipes.write_if_absent(
                recipe_hash, {**recipe_data, "metadata": {"first_seen": now}}
            )

            execution_record = {
                "recipe_hash": recipe_hash,
                "context": {
                    "command_line": command_line,
                    "working_directory": os.getcwd(),
                    "user": os.getenv("USER", "unknown"),
                },
                "requested": {
                    "main_options": (
                        redact_params(ctx.params) if ctx and ctx.params else {}
                    )
                },
                "execution": {
                    "job_id": job_id,
                    "joblog": None,  # To be filled in later if available
                    "status": "completed" if success else "failed",
                    "timestamp": now,
                    "exit_code": (
                        resource_data.get("exit_status") if resource_data else None
                    ),
                },
            }
            if resource_data:
                execution_record["resources"] = _resources_section(resource_data)
            if error:
                execution_record["execution"]["error"] = error
            execution_record["tags"] = tags or []
            if file_paths:
                execution_record["files"] = file_paths

            execution_timestamp = self._executions.create_unique(
                int(execution_timestamp), execution_record
            )
            self._executions.prune(
                keep=MAX_EXECUTIONS,
                slack=PRUNE_SLACK,
                protect=self._awaiting_update,
            )
            return execution_timestamp

        except Exception as e:
            # Don't let history logging break the main command
            logger.warning("Failed to log execution to history: %s", e)
            return execution_timestamp

    def _awaiting_update(self, key: str) -> bool:
        """True for a recent record that has no job ID yet.

        The key is the creation time in microseconds, so only records younger
        than ``AWAITING_UPDATE_SECONDS`` are read.
        """
        cutoff = (time.time() - AWAITING_UPDATE_SECONDS) * 1000000
        if int(key) < cutoff:
            return False
        record = self._executions.read(key)
        if record is None:
            return True  # claimed but not written yet: its writer is still running
        return not (record.get("execution") or {}).get("job_id")

    def update_execution_with_resources(
        self, execution_timestamp: str, job_id: str, resource_data: Dict[str, Any]
    ) -> bool:
        """Update an execution record with resource data after job completion."""
        try:
            record = self._executions.read(str(execution_timestamp))
            if record is None:
                logger.warning(
                    "History record %s not found; job %s not recorded",
                    execution_timestamp,
                    job_id,
                )
                return False

            execution = record.setdefault("execution", {})
            execution["job_id"] = job_id
            execution["exit_code"] = resource_data.get("exit_status")
            if resource_data:
                record["resources"] = _resources_section(resource_data)

            self._executions.write(str(execution_timestamp), record)
            return True

        except Exception as e:
            logger.warning("Failed to update execution with resources: %s", e)
            return False

    # -- reading -------------------------------------------------------------

    def _all_executions(self) -> Iterator[Dict[str, Any]]:
        """Execution records, newest first, each with its ``timestamp`` key."""
        self._ensure_private()
        self._import_legacy_yaml()
        for ts, record in self._executions.items_newest_first():
            record["timestamp"] = ts
            yield record

    def get_recipes(self) -> Dict[str, Any]:
        """All recipes, with ``run_count`` and ``last_seen`` filled in.

        Both are derived from the execution records rather than stored as
        counters, which every process would have to rewrite. A recipe
        imported from ``recipes.yaml`` keeps its old totals as a baseline
        (``metadata.imported``); executions logged after the import are
        added to it, and earlier ones are already included in it.

        Executions pruned beyond the newest ``MAX_EXECUTIONS`` stop
        counting, so for recipes created after the import ``run_count`` is
        the number of runs still in history, not a lifetime total.
        """
        self._ensure_private()
        self._import_legacy_yaml()
        recipes = {}
        for recipe_hash in self._recipes.keys():
            recipe = self._recipes.read(recipe_hash)
            if recipe is not None:
                recipes[recipe_hash] = recipe

        def imported(recipe_hash: str) -> Dict[str, Any]:
            recipe = recipes.get(recipe_hash) or {}
            return (recipe.get("metadata") or {}).get("imported") or {}

        usage: Dict[str, Dict[str, Any]] = {}
        for execution in self._all_executions():
            h = execution.get("recipe_hash")
            if int(execution["timestamp"]) <= int(imported(h).get("until", 0)):
                continue  # already counted in the imported baseline
            when = (execution.get("execution") or {}).get("timestamp")
            if h not in usage:
                usage[h] = {"last_seen": when, "run_count": 0}  # newest first
            usage[h]["run_count"] += 1

        for recipe_hash, recipe in recipes.items():
            metadata = recipe.setdefault("metadata", {})
            baseline = imported(recipe_hash)
            used = usage.get(recipe_hash, {})
            metadata["run_count"] = int(baseline.get("run_count") or 0) + used.get(
                "run_count", 0
            )
            seen = [t for t in (baseline.get("last_seen"), used.get("last_seen")) if t]
            metadata["last_seen"] = max(seen) if seen else metadata.get("first_seen")
        return recipes

    def get_executions(
        self,
        limit: int = 10,
        tags: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Get recent executions, optionally filtered by tags."""
        result = []
        for execution in self._all_executions():
            if tags:
                record_tags = list(execution.get("tags") or [])
                if not all(t in record_tags for t in tags):
                    continue
            result.append(execution)
            if len(result) >= limit:
                break
        return result

    def get_recipe_by_hash(self, recipe_hash: str) -> Optional[Dict[str, Any]]:
        """Get a specific recipe by hash."""
        return self.get_recipes().get(recipe_hash)

    def get_executions_for_recipe(
        self, recipe_hash: str, limit: int = 10
    ) -> List[Dict[str, Any]]:
        """Get executions for a specific recipe, most recent first."""
        result = []
        for execution in self._all_executions():
            if execution.get("recipe_hash") == recipe_hash:
                result.append(execution)
                if len(result) >= limit:
                    break
        return result

    def clear_all_history(self):
        """Clear all history (recipes and executions)."""
        self._executions.clear()
        self._recipes.clear()
        for legacy in (self.executions_file, self.recipes_file):
            if legacy.exists():
                legacy.unlink()

    # -- migration from the YAML files ---------------------------------------

    def _import_legacy_yaml(self) -> None:
        """Import ``executions.yaml``/``recipes.yaml`` written by older qxub.

        Records are written only if absent, so running this twice, in two
        processes at once, or after an older qxub re-created the YAML files,
        is harmless. Each file is renamed to ``*.imported`` afterwards.
        """
        for path, store, key in (
            (self.recipes_file, self._recipes, "recipes"),
            (self.executions_file, self._executions, "executions"),
        ):
            if not path.exists():
                continue
            try:
                from omegaconf import OmegaConf

                data = OmegaConf.to_container(OmegaConf.load(path)) or {}
            except Exception as e:  # unreadable: leave it for a human
                logger.warning("Could not import %s: %s", path, e)
                continue
            until = str(int(time.time() * 1000000))
            for record_key, record in (data.get(key) or {}).items():
                if not isinstance(record, dict):
                    continue
                if key == "recipes":
                    record = _with_imported_baseline(record, until)
                store.write_if_absent(str(record_key), record)
            imported_path = path.with_name(path.name + ".imported")
            try:
                os.replace(path, imported_path)
            except FileNotFoundError:
                pass  # another process got there first
            make_owner_only(imported_path)

    def convert_recipe_to_alias(self, recipe_hash: str) -> Optional[Dict[str, Any]]:
        """Convert a recipe to alias format for backward compatibility."""
        recipe = self.get_recipe_by_hash(recipe_hash)
        if not recipe:
            return None

        # Convert to alias format
        alias_def = {}

        # Map executor to subcommand format for compatibility
        if "executor" in recipe:
            alias_def["subcommand"] = recipe["executor"]

        # Map target format
        if "target" in recipe:
            alias_def["target"] = recipe["target"]

        return alias_def

    def get_execution_by_job_id(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get execution record by PBS job ID.

        Supports both full IDs (e.g., "153392916.gadi-pbs") and partial numeric
        IDs (e.g., "153392916"). Partial IDs match the portion before the first dot.

        Returns:
            Execution record dictionary with 'timestamp' field added, or None if not found.
            The record includes 'files' section with paths to out, err, and joblog files.
        """
        for execution in self._all_executions():
            stored_id = (execution.get("execution") or {}).get("job_id") or ""
            if stored_id == job_id or stored_id.split(".")[0] == job_id:
                return execution
        return None

    def get_most_recent_execution(self) -> Optional[Dict[str, Any]]:
        """Get the most recent execution record.

        Returns:
            Most recent execution record dictionary with 'timestamp' field, or None if no executions.
            Equivalent to calling get_executions(limit=1)[0] but more convenient.
        """
        executions = self.get_executions(limit=1)
        if executions:
            return executions[0]
        return None

    def get_execution_file_paths(self, job_id: str = None) -> Optional[Dict[str, str]]:
        """Get file paths (stdout, stderr, PBS log) for a job.

        Args:
            job_id: PBS job ID. If None, returns paths for most recent execution.

        Returns:
            Dictionary with 'out', 'err', and 'joblog' keys containing file paths,
            or None if execution not found. Used by history CLI commands to locate
            job output files for viewing.

        Example:
            >>> manager.get_execution_file_paths("153392916.gadi-pbs")
            {
                'out': '/scratch/ab12/user/qxub/job_20251026_205314.out',
                'err': '/scratch/ab12/user/qxub/job_20251026_205314.err',
                'joblog': '/scratch/ab12/user/qxub/job_20251026_205314.log'
            }
        """
        if job_id:
            execution = self.get_execution_by_job_id(job_id)
        else:
            execution = self.get_most_recent_execution()

        if execution and "files" in execution:
            return execution["files"]
        return None


# Global history manager instance
history_manager = HistoryManager()
