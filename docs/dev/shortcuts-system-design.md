# qxub Shortcuts System Design (v2.4.0)

## Overview

The shortcuts system is a new feature introduced in qxub v2.4.0 to complement the existing alias system. While aliases handle resource configuration and PBS job settings, shortcuts focus on command automation and execution patterns.

## Design Philosophy

### Separation of Concerns
- **Aliases**: Resource configuration, PBS settings, execution environments (YAML-based)
- **Shortcuts**: Command automation, execution patterns, workflow automation (JSON-based)

### Key Principles
- **Performance**: JSON storage with dictionary-based O(1) lookup optimization
- **Simplicity**: Clear command-based interface without complex hierarchical structure
- **Integration**: Seamless combination with existing alias system via `--alias` flag
- **Discoverability**: Easy-to-use configuration commands similar to alias management

## Architecture

### Core Components

#### ShortcutManager (`qxub/shortcut_manager.py`)
- **Purpose**: Core shortcuts management with optimized lookup
- **Storage**: JSON-based configuration in XDG-compliant directories
- **Features**:
  - Dictionary-based storage for O(command_length) lookup performance
  - XDG Base Directory specification compliance
  - Automatic config directory creation
  - Comprehensive error handling

#### exec_cli (`qxub/exec_cli.py`)
- **Purpose**: Integrated shortcut execution within the main `qxub exec` command
- **Features**:
  - Shortcut execution via `--shortcut <name>` option
  - Automatic shortcut detection based on command prefix matching
  - CLI option precedence (CLI > shortcut > defaults)
  - Full integration with all execution contexts (conda, modules, singularity)

#### Configuration Commands (`qxub/shortcuts_cli.py`)
- **Purpose**: Management interface for shortcuts via `qxub shortcut` subcommand
- **Commands**:
  - `qxub shortcut list` - List all shortcuts with rich table display
  - `qxub shortcut show <name>` - Show shortcut details with context descriptions
  - `qxub shortcut set <name> [OPTIONS]` - Create/update shortcuts with execution context
  - `qxub shortcut delete <name>` - Remove shortcuts
  - `qxub shortcut rename <old> <new>` - Rename shortcuts

### Storage Format

The key innovation is using command strings (or prefixes) as dictionary keys for O(1) lookup performance:

```json
{
  "python analyze.py": {
    "alias": "ml-gpu",
    "qxub_options": ["--queue", "normal", "--name", "analysis"]
  },
  "pytest": {
    "qxub_options": ["--queue", "normal"]
  },
  "snakemake": {
    "alias": "compute-intensive"
  }
}
```

**Key Design**: The dictionary key IS the command (or command prefix). When a user types `qxub exec -- python analyze.py`, the system looks up "python" in the shortcuts dictionary for automatic matching, or the user can explicitly specify `--shortcut python` for named shortcut usage.

### Configuration Hierarchy

1. **CLI arguments** (highest precedence)
2. **Shortcut qxub_options**
3. **Alias configuration** (if specified)
4. **Global defaults** (lowest precedence)

## Implementation Status

### ✅ Completed Components

#### Core ShortcutManager
- [x] XDG-compliant config directory resolution
- [x] JSON-based storage and retrieval
- [x] Dictionary optimization for fast lookup
- [x] CRUD operations (create, read, update, delete)
- [x] Error handling and validation

#### run_cli Command
- [x] Basic shortcut execution
- [x] Command override capability
- [x] Alias integration via `--alias` flag
- [x] CLI option pass-through
- [x] Comprehensive help system

#### Configuration Management
- [x] `qxub config shortcut list` command
- [x] `qxub config shortcut show` command
- [x] `qxub config shortcut set` command with options
- [x] `qxub config shortcut delete` command

#### Version Management
- [x] Version bump to 2.4.0 in `__init__.py` and `setup.py`
- [x] Changelog entry for shortcuts feature

### 🔄 In Progress

#### CLI Integration
- [ ] Fix missing `validate_execution_context` import in `cli.py`
- [x] Integrate shortcuts into `qxub exec` command (completed)
- [ ] Complete installation and testing

### ❌ Pending Implementation

#### Documentation
- [ ] Update `README.md` with shortcuts examples
- [ ] Add shortcuts section to `docs/examples.md`
- [ ] Update `docs/configuration.md` with shortcut config details

#### Testing
- [ ] Unit tests for ShortcutManager
- [ ] Integration tests for run_cli command
- [ ] Configuration command tests
- [ ] End-to-end workflow tests

#### Advanced Features
- [ ] Shortcut discovery and suggestion system
- [ ] Integration with bash completion
- [ ] Shortcut templates and parameterization
- [ ] History integration (create shortcuts from history)

## Usage Examples

### Basic Shortcut Creation and Execution

```bash
# Create shortcuts using command names as keys
qxub shortcut set "python" --env data-science --resources mem=8GB

# Run using automatic shortcut detection
qxub exec -- python analyze.py

# Or explicitly specify the shortcut
qxub exec --shortcut python -- analyze.py

# Prefix matching for convenience
qxub shortcut set "pytest" --env testing

# This matches the "pytest" shortcut
qxub exec -- pytest tests/unit/
```

### Advanced Shortcut with Alias Integration

```bash
# Create shortcut using command prefix as key
qxub config shortcut set "python train_model.py" --alias ml-gpu

# Run with combined shortcut + alias configuration
qxub exec -- python train_model.py

# Override queue setting from shortcut
qxub exec --queue express -- python train_model.py
```

### Shortcut Management

```bash
# List all shortcuts (shows command strings as keys)
qxub config shortcut list

# Show shortcut details for a specific command
qxub config shortcut show "python analyze.py"

# Delete shortcut by command string
qxub config shortcut delete "old-workflow.py"
```

## Technical Challenges Encountered

### CLI Integration Complexity
The qxub CLI uses a custom `QxubGroup` class that overrides Click's command resolution behavior to handle the hybrid execution/management command interface. Adding new commands requires careful integration with this system.

### Import Dependencies
The main CLI module has complex import dependencies due to the unified architecture. Missing imports can cause silent failures where commands aren't registered.

### Performance Considerations
Chose JSON over YAML for shortcuts to avoid parsing overhead during command lookup, especially important for interactive use.

## Future Enhancements

### Shortcut Parameterization
```bash
# Future: Command prefix matching with parameter extension
qxub config shortcut set "python analyze.py" --alias ml-env

# The shortcut would match and user can add parameters
qxub exec -- python analyze.py --input /data/file1.csv --output /results/
```

### Workflow Shortcuts
```bash
# Command pipeline shortcuts
qxub config shortcut set "snakemake --profile pbs" --alias compute-cluster

# Run the workflow shortcut
qxub exec -- snakemake --profile pbs --jobs 20
```

### Integration with External Tools
```bash
# Use command prefixes as keys, not arbitrary names
qxub config shortcut set "snakemake" --alias compute-intensive
qxub config shortcut set "jupyter lab" --alias interactive

# Run these shortcuts with additional parameters
qxub exec -- snakemake --jobs 10 analysis.smk
qxub exec -- jupyter lab --port 8888 --no-browser
```

## Migration Path

### From Aliases to Shortcuts
Users can gradually migrate command-focused aliases to shortcuts:

1. Keep resource/environment aliases as-is
2. Extract command portions into shortcuts
3. Use `--alias` flag to combine both systems

### Backward Compatibility
- All existing alias functionality remains unchanged
- Shortcuts are additive, not replacement
- No breaking changes to existing workflows

## Related Documentation

- `docs/dev/click_argument_parsing_solution.md` - CLI architecture background
- `docs/examples.md` - Will include shortcut examples
- `docs/configuration.md` - Will include shortcut configuration details

---

*Implementation Status: In Progress (v2.4.0 development)*
*Last Updated: October 18, 2025*
