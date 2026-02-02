# qxub Developer Documentation

Technical documentation for developers working on qxub v3.0.

## 🎉 v3.0 Migration Complete

The package structure migration is **COMPLETE**. All modules have been successfully reorganized into logical packages. See [migration-complete.md](migration-complete.md) for full details.

## Current Architecture (v3.0)

### Package Structure
```
qxub/
├── *_cli.py (9 files)    # CLI commands at root level
├── config/               # Configuration management
├── core/                 # Core utilities (scheduler, templates, parameters)
├── execution/            # Execution contexts and logic
├── history/              # Job history management
├── platform/             # Platform detection and queue selection
├── remote/               # Remote execution via SSH
├── resources/            # Resource parsing and management
└── jobscripts/           # PBS script templates
```

See: [package-structure-specification.md](package-structure-specification.md)

### Execution Model
qxub uses a simple, single-threaded execution model for job monitoring:
- **Status file monitoring**: Polls job status files for state changes
- **Simple spinner**: Non-blocking spinner during job submission and startup
- **Signal handling**: Clean Ctrl+C handling with automatic job cleanup via `qdel`
- **Real-time output**: Direct streaming of job stdout/stderr once job starts

See: [threading-architecture.md](threading-architecture.md)

### Configuration System
- **Hierarchical precedence**: CLI args > User config > System config > Defaults
- **Template variables**: `{user}`, `{project}`, `{timestamp}` for dynamic substitution
- **Alias system**: Hierarchical structure (main/subcommand/target) for resource and PBS configuration
- **Shortcuts system**: Command-based automation with prefix matching for workflow patterns

See: [config-and-alias-system-design.md](config-and-alias-system-design.md), [config_schema.md](config_schema.md), [shortcuts-system-design.md](shortcuts-system-design.md)

### Platform System
- **Platform definitions**: YAML files describing HPC system capabilities
- **Queue selection**: `--queue auto` for intelligent, cost-optimized selection
- **Resource validation**: Platform-aware constraint checking

See: [platform_schema.md](platform_schema.md)

### CLI Architecture
- **Standard Click interface**: Clean modular design with `qxub exec` subcommand
- **Execution contexts**: `--env`, `--mod/--mods`, `--sif`, `--default` for different environments
- **Management commands**: Preserved as separate subcommands (`config`, `alias`, `history`, etc.)
- **Comprehensive PBS options**: All PBS features accessible via `qxub exec`

### Remote Execution
- **SSH-based**: Submit jobs to remote HPC systems from local machine
- **Platform inheritance**: Remote systems use local platform definitions
- **Output streaming**: Real-time output from remote jobs

See: [remote_execution.md](remote_execution.md)

## Debug Commands
```bash
# Enable debug logging for execution issues
export QXUB_LOG_LEVEL=DEBUG
qxub exec --env myenv -- python script.py

# Test dry-run without job submission
qxub exec --dry-run --env myenv -- python script.py

# Check configuration hierarchy
qxub config files

# Check job status files for debugging
ls -la /scratch/$USER/qxub/status/
```

## Key Files (v3.0 Structure)
- `qxub/cli.py` - Main CLI entry point with Click interface
- `qxub/exec_cli.py` - Comprehensive execution subcommand with PBS options
- `qxub/execution.py` - Unified execution interface
- `qxub/execution/context.py` - Execution context logic
- `qxub/config/manager.py` - Configuration loading and template resolution
- `qxub/platform.py` - Platform abstraction and queue selection
- `qxub/core/scheduler.py` - PBS interaction and job monitoring
- `qxub/remote/executor.py` - SSH-based remote execution

## Development Workflow

```bash
# Run tests after making changes
./tests/test_realistic_system_config.sh
./tests/test_command_edge_cases.sh

# Check for regressions
qxub exec --dry-run --env base -- echo "test"

# Test specific execution contexts
qxub exec --dry-run --mod python3 -- python --version
qxub exec --dry-run --sif container.sif -- whoami
```
