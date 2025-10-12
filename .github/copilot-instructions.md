# qxub Copilot Instructions

## Project Overview

qxub is a sophisticated PBS job submission wrapper for HPC environments that eliminates boilerplate when running jobs in conda environments, with modules, or in containers. The codebase follows a unified CLI architecture with execution contexts, intelligent queue selection, and comprehensive configuration management.

## Key Architecture Patterns

### Unified CLI Structure
- **Core principle**: All options before `--` separator, command after
- **Execution contexts**: `--env` (conda), `--mod`/`--mods` (modules), `--sif` (singularity)
- **Mutual exclusivity**: Only one execution context per command
- **Default execution**: Direct PBS submission when no context specified

```python
# CLI parsing in qxub/cli.py uses custom QxubGroup class
# Execution context detection happens in invoke() method
execution_contexts = [conda_env, module_list, container]
if sum(bool(x) for x in execution_contexts) > 1:
    raise click.ClickException("Cannot specify multiple execution contexts")
```

### Configuration System
- **Hierarchical precedence**: CLI args > User config > System config > Defaults
- **XDG compliance**: Uses `~/.config/qxub/` and `/etc/xdg/qxub/`
- **Template variables**: `{user}`, `{project}`, `{timestamp}` for dynamic substitution
- **Alias system**: Hierarchical structure with main/subcommand/target sections

```python
# Config manager in qxub/config_manager.py
class ConfigManager:
    def _get_user_config_dir(self) -> Path:
        xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
        if xdg_config_home:
            return Path(xdg_config_home) / "qxub"
        return Path.home() / ".config" / "qxub"
```

### Platform System
- **Platform definitions**: YAML files describing HPC system capabilities
- **Queue selection**: `--queue auto` for intelligent selection based on resources
- **Resource validation**: Platform-aware constraint checking
- **Environment discovery**: `QXUB_PLATFORM_PATHS` for custom platform locations

```python
# Platform loader in qxub/platform.py
@dataclass
class Queue:
    name: str
    type: str
    limits: QueueLimits
    su_billing_rate: Optional[float] = None
    walltime_rules: List[WalltimeRule] = field(default_factory=list)
```

## Critical Implementation Details

### Threading Architecture
- **OutputCoordinator**: Central synchronization hub for job monitoring
- **Signal handling**: Ctrl+C cleanup with `_signal_handler()` and global `_CURRENT_JOB_ID`
- **Thread responsibilities**: Monitor thread, STDOUT/STDERR tail, spinner display
- **Graceful shutdown**: Automatic job cleanup via `qdel` on interruption

### Job Script Generation
- **Base64 encoding**: Commands encoded to avoid shell escaping issues
- **Template system**: Configurable job script templates in `qxub/jobscripts/`
- **Environment activation**: Context-specific activation (conda/module/singularity)

```python
# Command encoding pattern in execute_* functions
cmd_b64 = base64.b64encode(cmd_str.encode("utf-8")).decode("ascii")
```

### Resource Management
- **Resource tracking**: SQLite-based efficiency analysis in `qxub/resource_tracker.py`
- **History system**: Dual logging (computational recipes + execution records)
- **Memory/walltime parsing**: Robust parsing with format validation

## Development Workflows

### Environment Setup
- **Virtual environment**: Always activate the qxub virtual environment before development
- **Environment verification**: Check `qxub --version` to ensure correct installation
- **Terminal sessions**: Activate environment in each new terminal session

```bash
# Activate virtual environment (adjust path as needed)
source /path/to/qxub-venv/bin/activate

# Verify qxub is available
qxub --version
which qxub
```

### Code Style
- **Black formatting**: Line length 88, Python 3.10+ target
- **Pre-commit hooks**: Black, isort, trailing whitespace, YAML validation
- **Pylint**: Extensive disable comments for specific patterns

### Testing Strategy
- **Dry-run tests**: `tests/test_conda_dry.sh` for rapid development iteration
- **Integration tests**: `tests/test_conda_integration.sh` with real job submission
- **System config tests**: `tests/test_system_config.sh` for config hierarchy validation
- **Platform tests**: `tests/run_platform_tests.py` for queue selection logic

### Critical Testing Commands
```bash
# Fast development validation
./tests/test_conda_dry.sh

# Test config precedence
./tests/test_realistic_system_config.sh

# Platform system validation
python tests/run_platform_tests.py
```

## Common Patterns

### Error Handling
- **Exit codes**: Validation errors return code 2, execution errors return 1
- **Fuzzy matching**: Typo suggestions for unknown options using `difflib`
- **Context-aware errors**: Different error messages based on execution context

### CLI Command Structure
```python
# Standard execution function pattern
def execute_conda(ctx, command, conda_env, template, pre, post):
    params = ctx.obj  # Access to global parameters
    # Validation, template processing, job submission
```

### Config Management Commands
- `qxub config get/set/list` - Configuration manipulation
- `qxub config alias set/list/show/delete` - Alias management
- `qxub config files` - Show config file locations and status

## File Organization
- **qxub/cli.py**: Main CLI entry point and execution contexts
- **qxub/config_manager.py**: Configuration loading and template resolution
- **qxub/platform*.py**: Platform abstraction and queue selection
- **qxub/*_cli.py**: Modular CLI commands (config, alias, history, platform)
- **qxub/scheduler.py**: PBS interaction and job monitoring
- **docs/dev/**: Technical architecture documentation

## Remote Execution (Upcoming)
- **SSH-based execution**: `qxub --remote REMOTE_NAME` pattern emerging
- **Platform integration**: Remote platforms inherit local platform definitions
- **Config separation**: SSH config vs platform definitions vs user preferences

When working on this codebase, prioritize understanding the execution context flow in `cli.py`, configuration precedence in `config_manager.py`, and platform-aware resource validation patterns.
