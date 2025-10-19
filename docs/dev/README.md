# qxub Developer Documentation

Technical documentation for developers working on qxub v2.3.2.

## Core Architecture

### Threading System
qxub uses multi-threaded architecture for real-time job monitoring:
- **OutputCoordinator**: Central synchronization hub for job monitoring
- **Signal handling**: Ctrl+C cleanup with `_signal_handler()` and global `_CURRENT_JOB_ID`
- **Thread responsibilities**: Monitor thread, STDOUT/STDERR tail, spinner display
- **Graceful shutdown**: Automatic job cleanup via `qdel` on interruption

See: [threading-architecture.md](threading-architecture.md)

### Configuration System
- **Hierarchical precedence**: CLI args > User config > System config > Defaults
- **Template variables**: `{user}`, `{project}`, `{timestamp}` for dynamic substitution
- **Alias system**: Hierarchical structure with main/subcommand/target sections

See: [config-and-alias-system-design.md](config-and-alias-system-design.md), [config_schema.md](config_schema.md)

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
# Enable debug logging for threading issues
export QXUB_LOG_LEVEL=DEBUG
qxub --env myenv -- python script.py

# Test dry-run without job submission
qxub --dry-run --env myenv -- python script.py

# Check configuration hierarchy
qxub config files
```

## Key Files
- `qxub/cli.py` - Clean main CLI entry point with standard Click interface
- `qxub/exec_cli.py` - Comprehensive execution subcommand with all PBS options
- `qxub/execution_context.py` - Unified execution logic for all environment types
- `qxub/config_manager.py` - Configuration loading and template resolution
- `qxub/platform.py` - Platform abstraction and queue selection
- `qxub/scheduler.py` - PBS interaction and job monitoring
- `qxub/remote_executor.py` - SSH-based remote execution

# Check thread status
ps aux | grep qxub

# View active threads in Python
python -c "import threading; print([t.name for t in threading.enumerate()])"

# Kill hung qxub processes
pkill -TERM qxub
```

### Key Files for Threading
- `qxub/scheduler.py` - Core threading logic, OutputCoordinator
- `qxub/conda.py` - Conda executor with threading integration
- `qxub/module.py` - Module executor with threading integration
- `qxub/sing.py` - Singularity executor with threading integration

### Testing Threading Code
```python
# Always use timeouts in tests
thread.join(timeout=5)
assert not thread.is_alive()

# Test shutdown scenarios
coordinator.signal_shutdown()
assert coordinator.should_shutdown()

# Mock the coordinator for unit tests
from unittest.mock import Mock
coordinator = Mock()
```

## Contributing to Threading System

When modifying the threading system:

1. **Read the architecture docs first** - Understand the event-based coordination
2. **Test shutdown scenarios** - Ensure clean Ctrl-C handling
3. **Use timeouts everywhere** - Prevent hanging threads
4. **Follow the event patterns** - Use Events for coordination, not shared variables
5. **Add debug logging** - Include thread names and timing information
6. **Test with real PBS jobs** - Integration testing is crucial

## Getting Help

If you encounter threading-related issues:

1. **Check the troubleshooting guide** - Most common issues are documented
2. **Enable debug logging** - Use `QXUB_LOG_LEVEL=DEBUG`
3. **Gather thread dumps** - Include in bug reports
4. **Test in isolation** - Reproduce with minimal examples
5. **Document edge cases** - Help improve these docs

For questions or contributions, see the main project README for contact information.
