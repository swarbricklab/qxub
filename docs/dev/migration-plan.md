# qxub Package Structure Migration Plan

## Overview

This document outlines the step-by-step migration from the current flat package structure to the organized multi-package architecture. The migration is designed to minimize risk and maintain backwards compatibility throughout the process.

## Current State Analysis

### Current File Count and Organization
```
qxub/ (35+ files in flat structure)
├── CLI Commands: 9 files (*_cli.py)
├── Configuration: 4 files (config*.py)
├── Platform System: 3 files (platform*.py)
├── Execution Logic: 4 files (execution*.py, executors.py)
├── Scheduling: 1 file (scheduler.py)
├── Resources: 4 files (resource*.py)
├── Remote Execution: 4 files (remote*.py)
├── History: 2 files (history*.py)
├── Utilities: 4+ files (templates.py, parameters.py, etc.)
└── Legacy: 1 file (cli_old.py)
```

### Dependency Assessment
Based on import analysis:
- **High Coupling**: CLI modules import from many other modules
- **Medium Coupling**: Execution and platform modules
- **Low Coupling**: Resource utilities, configuration management
- **Circular Dependencies**: Some exist between config, platform, and execution

## Migration Strategy

### Overall Approach
1. **Risk-Based Phases** - Move low-risk packages first
2. **Backwards Compatibility** - Maintain import compatibility during transition
3. **Testing at Each Step** - Comprehensive testing after each phase
4. **Gradual Cleanup** - Remove compatibility shims only after migration complete

### Phase Breakdown

## Phase 1: Foundation Packages (Low Risk)

**Duration**: 1-2 days
**Risk Level**: Low
**Goal**: Establish basic package structure with minimal dependencies

### Phase 1.1: Resources Package
**Files to Move**:
```
qxub/resource_parser.py        → qxub/resources/parser.py
qxub/resource_utils.py         → qxub/resources/utils.py
qxub/resource_tracker.py       → qxub/resources/tracker.py
qxub/resource_mappers.py       → qxub/resources/mappers.py
```

**Rationale**: Resources package has no dependencies, pure utility functions.

**Steps**:
1. Create `qxub/resources/` directory
2. Move files with minimal modifications
3. Update `qxub/resources/__init__.py` with public API
4. Add compatibility imports in `qxub/__init__.py`
5. Run Phase 1 test suite

**Compatibility Layer**:
```python
# qxub/__init__.py - Maintain backwards compatibility
from .resources.parser import parse_memory_size, parse_walltime
from .resources.utils import format_walltime, bytes_to_human
from .resources.tracker import resource_tracker
```

### Phase 1.2: History Package
**Files to Move**:
```
qxub/history.py                → qxub/history/base.py
qxub/history_manager.py        → qxub/history/manager.py
```

**Rationale**: History package has minimal dependencies (only resources utilities).

**Steps**:
1. Create `qxub/history/` directory
2. Move files and update internal imports
3. Update history CLI to use new imports
4. Run Phase 1 test suite

## Phase 2: Configuration Package (Low-Medium Risk)

**Duration**: 2-3 days
**Risk Level**: Low-Medium
**Goal**: Centralize configuration management

### Phase 2.1: Configuration Core
**Files to Move**:
```
qxub/config.py                 → qxub/config/base.py
qxub/config_manager.py         → qxub/config/manager.py
qxub/config_handler.py         → qxub/config/handler.py
qxub/shortcut_manager.py       → qxub/config/shortcuts.py
qxub/standalone_aliases.py     → qxub/config/aliases.py (optional)
```

**Complexity**: Medium - Configuration is imported by many modules

**Steps**:
1. Create `qxub/config/` directory
2. Move core config files first (`config.py`, `config_manager.py`)
3. Update imports throughout codebase
4. Move secondary files (`config_handler.py`, `shortcut_manager.py`)
5. Extensive testing due to wide usage

**Import Updates Required**:
```python
# Many files need updates like:
from .config_manager import config_manager
# becomes:
from .config import config_manager
```

## Phase 3: CLI Package (Medium Risk)

**Duration**: 3-4 days
**Risk Level**: Medium
**Goal**: Organize command-line interface

### Phase 3.1: CLI Commands
**Files to Move**:
```
qxub/exec_cli.py               → qxub/cli/exec_cli.py
qxub/config_cli.py             → qxub/cli/config_cli.py
qxub/alias_cli.py              → qxub/cli/alias_cli.py
qxub/history_cli.py            → qxub/cli/history_cli.py
qxub/monitor_cli.py            → qxub/cli/monitor_cli.py
qxub/cancel_cli.py             → qxub/cli/cancel_cli.py
qxub/status_cli.py             → qxub/cli/status_cli.py
qxub/resources_cli.py          → qxub/cli/resources_cli.py
qxub/platform_cli.py           → qxub/cli/platform_cli.py
```

**Complexity**: Medium - CLI files import from many other modules

**Steps**:
1. Create `qxub/cli/` directory
2. Move CLI files with updated imports
3. Update main `cli.py` to import from new locations
4. Test all CLI commands

**Import Updates Example**:
```python
# Before: qxub/exec_cli.py
from .config_manager import config_manager
from .platform import get_platform

# After: qxub/cli/exec_cli.py
from ..config import config_manager
from ..platforms import get_platform
```

## Phase 4: Platform Package (High Risk)

**Duration**: 4-5 days
**Risk Level**: High
**Goal**: Prepare for multi-platform architecture

### Phase 4.1: Platform Abstraction
**Files to Move**:
```
qxub/platform.py               → qxub/platforms/base.py + loader.py
qxub/platform_integration.py   → qxub/platforms/integration.py
```

**Complexity**: High - Platform code is complex and widely used

**Steps**:
1. **Analyze** `platform.py` to identify logical splits
2. **Split** into `base.py` (classes) and `loader.py` (functions)
3. **Extract** platform detection to `detection.py`
4. **Create** PBS-specific implementation in `pbs_pro.py`
5. **Update** all imports throughout codebase
6. **Extensive testing** - platform logic is critical

**Split Strategy**:
```python
# qxub/platforms/base.py
class Platform:
    """Platform definition and queue management."""

class Queue:
    """Queue definition and validation."""

# qxub/platforms/loader.py
def get_platform(name: str) -> Platform:
    """Platform loading functions."""

def list_platforms() -> List[str]:
    """Platform discovery functions."""

# qxub/platforms/detection.py
def detect_platform() -> Optional[str]:
    """Platform auto-detection logic."""
```

### Phase 4.2: Platform CLI Integration
**Updates Required**:
```python
# qxub/cli/platform_cli.py
from ..platforms import get_platform, list_platforms, detect_platform
```

## Phase 5: Execution and Scheduling (High Risk)

**Duration**: 5-6 days
**Risk Level**: High
**Goal**: Separate execution orchestration from scheduler interaction

### Phase 5.1: Scheduling Package
**Files to Move**:
```
qxub/scheduler.py              → Split into multiple files:
                               → qxub/scheduling/pbs.py (PBS-specific functions)
                               → qxub/scheduling/monitoring.py (monitoring utilities)
                               → qxub/scheduling/scheduler.py (abstract interface)
```

**Complexity**: Very High - `scheduler.py` is 1100+ lines with complex logic

**Split Strategy**:
```python
# qxub/scheduling/pbs.py
def qsub(cmd: str, quiet: bool = False) -> str:
    """PBS job submission."""

def qdel(job_id: str, quiet: bool = False) -> bool:
    """PBS job deletion."""

def job_status(job_id: str) -> str:
    """PBS job status checking."""

# qxub/scheduling/monitoring.py
def monitor_job_single_thread(job_id: str, out_file: Path, err_file: Path):
    """Job monitoring logic."""

def stream_job_output(job_id: str, out_file: Path, err_file: Path):
    """Output streaming logic."""

# qxub/scheduling/scheduler.py (Future)
class SchedulerInterface:
    """Abstract scheduler interface for multi-platform support."""
```

### Phase 5.2: Execution Package
**Files to Move**:
```
qxub/execution_context.py       → qxub/execution/context.py
qxub/execution.py              → qxub/execution/unified.py
qxub/executors.py              → qxub/execution/executors.py
qxub/templates.py              → qxub/execution/templates.py
qxub/parameters.py             → qxub/execution/parameters.py (optional)
```

**Steps**:
1. Move execution context and unified execution logic
2. Update executor implementations
3. Update template management
4. Test all execution contexts (conda, modules, singularity)

## Phase 6: Remote Package (Medium Risk)

**Duration**: 2-3 days
**Risk Level**: Medium
**Goal**: Organize remote execution capabilities

### Phase 6.1: Remote Execution
**Files to Move**:
```
qxub/remote.py                 → qxub/remote/base.py
qxub/remote_config.py          → qxub/remote/config.py
qxub/remote_config_loader.py   → qxub/remote/loader.py
qxub/remote_executor.py        → qxub/remote/executor.py
```

**Steps**:
1. Create `qxub/remote/` directory
2. Move remote execution files
3. Update imports in CLI and execution modules
4. Test remote execution functionality

## Phase 7: Cleanup and Optimization

**Duration**: 2-3 days
**Risk Level**: Low
**Goal**: Final cleanup and optimization

### Phase 7.1: Legacy Removal
1. Remove `cli_old.py` (if no longer needed)
2. Remove compatibility shims from `__init__.py`
3. Clean up any remaining flat imports

### Phase 7.2: Import Optimization
1. Optimize import performance
2. Add lazy imports where appropriate
3. Final circular dependency cleanup

## Risk Mitigation Strategies

### Pre-Migration Preparation
1. **Comprehensive Test Suite** - Ensure 100% test coverage for migration
2. **Import Analysis** - Map all current imports to identify dependencies
3. **Backup Strategy** - Git branching strategy for easy rollback
4. **Documentation** - Document all import changes for reference

### During Migration
1. **One Phase at a Time** - Complete testing before moving to next phase
2. **Compatibility Layers** - Maintain backwards compatibility throughout
3. **Incremental Testing** - Run full test suite after each file move
4. **Rollback Plan** - Be prepared to revert any phase if issues arise

### Post-Migration Validation
1. **Full Regression Testing** - All existing functionality works
2. **Performance Testing** - No performance degradation
3. **Import Testing** - All public APIs still accessible
4. **Documentation Updates** - Update all documentation for new structure

## Testing Strategy Per Phase

### Phase 1 Testing
```bash
# Resources package tests
python -m pytest tests/test_resources/ -v
python -c "from qxub.resources import parse_memory_size; print('✓ Import works')"

# History package tests
python -m pytest tests/test_history/ -v
python -c "from qxub.history import history_manager; print('✓ Import works')"
```

### Phase 2 Testing
```bash
# Configuration tests
python -m pytest tests/test_config/ -v
qxub config get defaults.project  # Test CLI still works
python -c "from qxub.config import config_manager; print('✓ Import works')"
```

### Phase 3 Testing
```bash
# CLI tests
python -m pytest tests/test_cli/ -v
qxub --help                       # Test main CLI
qxub config --help               # Test subcommands
qxub exec --help                 # Test exec command
```

### Phase 4 Testing
```bash
# Platform tests
python -m pytest tests/test_platforms/ -v
qxub platform list               # Test platform CLI
python -c "from qxub.platforms import get_platform; print('✓ Import works')"
```

### Phase 5 Testing
```bash
# Execution and scheduling tests
python -m pytest tests/test_execution/ tests/test_scheduling/ -v
qxub exec --dry --env test -- echo "test"  # Test execution
python -c "from qxub.scheduling import qsub; print('✓ Import works')"
```

### Phase 6 Testing
```bash
# Remote execution tests
python -m pytest tests/test_remote/ -v
python -c "from qxub.remote import RemoteConfig; print('✓ Import works')"
```

### Full Integration Testing
```bash
# Complete regression test
./tests/test_unified_cli.sh        # Full CLI tests
./tests/test_conda_dry.sh          # Conda execution tests
./tests/test_realistic_system_config.sh  # Configuration tests
python tests/run_platform_tests.py # Platform tests
```

## Rollback Procedures

### Git Strategy
```bash
# Before each phase
git checkout -b migration-phase-N
git tag pre-phase-N

# If rollback needed
git checkout main
git reset --hard pre-phase-N
```

### Incremental Rollback
Each phase should be atomic - if issues arise:
1. Identify specific file causing issues
2. Revert only that file move
3. Fix issues before proceeding
4. Re-run phase tests

## Success Criteria

### Phase Completion Criteria
Each phase is complete when:
1. ✅ All files moved successfully
2. ✅ All imports updated correctly
3. ✅ Full test suite passes
4. ✅ No functionality regressions
5. ✅ Performance unchanged
6. ✅ Documentation updated

### Overall Migration Success
Migration is successful when:
1. ✅ All 35+ files organized into logical packages
2. ✅ No breaking changes for users
3. ✅ Clean import hierarchy established
4. ✅ Foundation ready for multi-platform features
5. ✅ Improved code maintainability
6. ✅ Reduced circular dependencies

## Timeline Estimate

### Conservative Timeline (20-25 days)
- **Phase 1**: 2 days (Resources + History)
- **Phase 2**: 3 days (Configuration)
- **Phase 3**: 4 days (CLI)
- **Phase 4**: 5 days (Platform)
- **Phase 5**: 6 days (Execution + Scheduling)
- **Phase 6**: 3 days (Remote)
- **Phase 7**: 2 days (Cleanup)

### Aggressive Timeline (12-15 days)
Possible with dedicated focus and good test coverage, but higher risk.

## Resources Required

### Development Resources
- **Primary Developer**: Full-time focus during migration
- **Testing**: Automated test suite + manual verification
- **Documentation**: Update guides and examples

### Infrastructure
- **Git Branching**: Separate branch for migration work
- **CI/CD**: Ensure tests run on all changes
- **Backup**: Multiple restore points during migration

This migration plan provides a structured approach to reorganizing qxub while minimizing risk and maintaining stability throughout the process.
