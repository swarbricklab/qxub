# qxub Package Structure Specification - FINAL IMPLEMENTATION

## Overview

This document defines the **actual implemented** package structure for qxub v3.0 after completing the migration. This reflects the final reality rather than the original plan, as some decisions changed during implementation.

## Design Principles

### Organizational Principles
1. **Single Responsibility** - Each package has one clear purpose
2. **Dependency Clarity** - Clear import hierarchy with minimal circular dependencies
3. **Future-Proof** - Structure accommodates planned features without refactoring
4. **Logical Grouping** - Related functionality lives together
5. **Public APIs** - Clean interfaces between packages

### Dependency Hierarchy
```
CLI Commands (qxub.*_cli.py + qxub/cli.py)
    ↓
Business Logic (qxub/execution/*, qxub/platform/*)
    ↓
Core Services (qxub/config/*, qxub/core/*)
    ↓
Utilities (qxub/resources/*, qxub/remote/*, qxub/history/*)
```

## Final Package Structure

### Root Level - CLI Commands
```
qxub/
├── __init__.py                    # Main package exports and version
├── cli.py                         # Main CLI entry point and command registration
├── exec_cli.py                   # Main execution command (qxub exec)
├── config_cli.py                # Configuration commands (qxub config)
├── alias_cli.py                  # Alias management (qxub config alias)
├── history_cli.py                # Job history (qxub history)
├── monitor_cli.py                # Job monitoring (qxub monitor)
├── cancel_cli.py                 # Job cancellation (qxub cancel)
├── status_cli.py                 # Job status (qxub status)
├── resources_cli.py              # Resource management (qxub resources)
├── platform_cli.py              # Platform commands (qxub platform)
├── config_handler.py            # Legacy config processing (compatibility)
├── config_manager.py            # Legacy config manager (compatibility)
├── execution_context.py         # Legacy execution context (compatibility)
├── execution.py                  # Legacy execution wrapper (compatibility)
├── platform.py                  # Legacy platform wrapper (compatibility)
├── shortcut_manager.py          # Legacy shortcut management (compatibility)
└── history_manager.py           # Legacy history management (compatibility)
```

**Purpose**: CLI commands at root level for direct access, with legacy compatibility wrappers.

**Public API**:
```python
# qxub/__init__.py
__version__ = "3.0.0"

# Main CLI entry point
from .cli import qxub

# Core public APIs
from .config import config_manager
from .platform import get_platform, select_best_queue
from .execution import execute_unified
from .core.scheduler import qsub, qdel, job_status
```

### Configuration Package (`qxub/config/`)
```
config/
├── __init__.py                   # Config module exports
├── manager.py                    # Main configuration manager
├── handler.py                    # Configuration processing logic
├── base.py                       # Base configuration classes
├── aliases.py                    # Alias management
└── shortcuts.py                  # Shortcut/alias management
```

**Purpose**: Configuration loading, validation, template processing, and hierarchy management.

**Dependencies**:
- Utilities: `resources` (for resource parsing)
- No dependencies on CLI or execution layers

**Public API**:
```python
# qxub/config/__init__.py
from .manager import config_manager, ConfigManager
from .base import ConfigurationError
from .aliases import AliasManager
from .shortcuts import shortcut_manager

__all__ = ['config_manager', 'ConfigManager', 'ConfigurationError', 'AliasManager', 'shortcut_manager']
```

### Core Package (`qxub/core/`)
```
core/
├── __init__.py                   # Core module exports
├── scheduler.py                  # PBS scheduler interface (moved from root)
├── parameters.py                 # Parameter processing utilities (moved from root)
└── templates.py                  # Job script templates (moved from root)
```

**Purpose**: Core scheduling, parameter processing, and template management utilities.

**Dependencies**:
- Utilities: `resources` (for resource parsing)
- No dependencies on higher-level packages

**Public API**:
```python
# qxub/core/__init__.py
from .scheduler import qsub, qdel, job_status, JobSubmitter
from .parameters import process_parameters
from .templates import get_template, list_templates

__all__ = ['qsub', 'qdel', 'job_status', 'JobSubmitter', 'process_parameters', 'get_template', 'list_templates']
```

### Platform Package (`qxub/platform/`)
```
platform/
├── __init__.py                   # Platform module exports
├── core.py                       # Platform core logic and loading
├── cli.py                        # Platform CLI commands
└── integration.py                # Platform integration utilities
```

**Purpose**: Platform abstraction, queue selection, resource validation, and scheduler-agnostic interfaces.

**Dependencies**:
- Core Services: `config` (for platform search paths)
- Utilities: `resources` (for resource parsing and validation)

**Public API**:
```python
# qxub/platform/__init__.py
from .core import Platform, Queue, QueueLimits, PlatformResources
from .core import get_platform, list_platforms, get_current_platform
from .integration import detect_platform

__all__ = [
    'Platform', 'Queue', 'QueueLimits', 'PlatformResources',
    'get_platform', 'list_platforms', 'get_current_platform', 'detect_platform'
]
```
- `PlatformResources`: Abstract resource specification
- `PBSPlatform`: PBS Pro implementation
- `PlatformLoader`: Platform discovery and loading

### Execution Package (`qxub/execution/`)
```
execution/
├── __init__.py                   # Execution module exports
├── context.py                    # Execution contexts (conda, modules, singularity)
├── core.py                       # Core execution logic and job submission
└── executors.py                  # Individual executor implementations
```

**Purpose**: Job execution orchestration, environment setup, and job script generation.

**Dependencies**:
- Core Services: `config`, `core` (scheduler, templates)
- Business Logic: `platform`
- Utilities: `resources`

**Public API**:
```python
# qxub/execution/__init__.py
from .core import execute_unified
from .context import ExecutionContext

__all__ = ['execute_unified', 'ExecutionContext']
```

### History Package (`qxub/history/`)
```
history/
├── __init__.py                   # History module exports
├── base.py                       # Base history classes
└── manager.py                    # History management implementation
```

**Purpose**: Job history tracking, computational recipe storage, and execution records.

**Dependencies**:
- Core Services: `config`
- Utilities: `resources`

**Public API**:
```python
# qxub/history/__init__.py
from .manager import HistoryManager
from .base import HistoryRecord

__all__ = ['HistoryManager', 'HistoryRecord']
```

### Remote Package (`qxub/remote/`)
```
remote/
├── __init__.py                   # Remote module exports
├── config.py                     # Remote configuration management
├── core.py                       # Remote execution core logic
├── executor.py                   # Remote job execution
└── loader.py                     # Remote platform loading
```

**Purpose**: SSH-based remote execution capabilities.

**Dependencies**:
- Core Services: `config`, `core`
- Business Logic: `platform`, `execution`

**Public API**:
```python
# qxub/remote/__init__.py
from .core import RemoteExecutor
from .config import RemoteConfig

__all__ = ['RemoteExecutor', 'RemoteConfig']
```

**Public API**:
```python
# qxub/scheduling/__init__.py
from .pbs import qsub, qdel, job_status, get_job_resource_data
from .monitoring import monitor_job_single_thread, stream_job_output

__all__ = ['qsub', 'qdel', 'job_status', 'get_job_resource_data',
           'monitor_job_single_thread', 'stream_job_output']
```

**Key Functions**:
- `qsub`: Job submission
- `qdel`: Job cancellation
- `job_status`: Status checking
- `monitor_job_single_thread`: Job monitoring

### Resources Package (`qxub/resources/`)
```
resources/
├── __init__.py                   # Resources module exports
├── parser.py                     # Resource parsing utilities
├── utils.py                      # Resource utility functions
├── tracker.py                    # Resource tracking and efficiency analysis
└── mappers.py                    # Workflow engine resource mappers
### Resources Package (`qxub/resources/`)
```
resources/
├── __init__.py                   # Resource module exports
├── parser.py                     # Resource parsing utilities
├── utils.py                      # Resource formatting and utilities
├── tracker.py                    # Resource usage tracking
└── mappers.py                    # Resource mapping and conversion
```

**Purpose**: Resource specification parsing, validation, conversion, and efficiency tracking.

**Dependencies**: None (pure utility functions)

**Public API**:
```python
# qxub/resources/__init__.py
from .parser import parse_memory_size, parse_walltime, size_to_bytes
from .utils import format_walltime, bytes_to_human
from .tracker import resource_tracker
from .mappers import ResourceMapper

__all__ = ['parse_memory_size', 'parse_walltime', 'size_to_bytes',
           'format_walltime', 'bytes_to_human', 'resource_tracker', 'ResourceMapper']
```

### JobScripts Package (`qxub/jobscripts/`)
```
jobscripts/
├── __init__.py                   # JobScripts module exports
├── qdefault.pbs                  # Default PBS job script template
├── qconda.pbs                    # Conda environment PBS template
├── qmod.pbs                      # Environment modules PBS template
└── qsing.pbs                     # Singularity container PBS template
```

**Purpose**: PBS job script templates for different execution contexts.

**Dependencies**: None (template files)

## Migration Completion Status

### ✅ Completed Phases
- **Phase 1**: Resources and History packages ✅
- **Phase 2**: Configuration package ✅
- **Phase 3**: CLI commands (kept at root level) ✅
- **Phase 4**: Platform package ✅
- **Phase 5**: Execution package ✅
- **Phase 6**: Core package (scheduler, parameters, templates) ✅
- **Phase 7**: Remote package ✅
- **Phase 8**: Cleanup and optimization ✅

### Key Architectural Decisions Made During Migration

1. **CLI Commands at Root Level**: Instead of creating a `qxub/cli/` package, CLI commands were kept at the root level (`*_cli.py`) for simpler imports and backwards compatibility.

2. **Core Package Created**: A new `qxub/core/` package was created to house fundamental utilities (scheduler, parameters, templates) that were previously at the root level.

3. **Platform vs Platforms**: The package was named `platform/` (singular) rather than `platforms/` (plural) as originally planned.

4. **Legacy Compatibility Wrappers**: Several root-level files (`config_manager.py`, `execution.py`, etc.) were retained as compatibility wrappers to maintain backwards compatibility.

5. **JobScripts Package**: Job script templates were organized into their own package rather than being embedded in the execution package.

### Final Package Count
- **8 packages**: config, core, execution, history, jobscripts, platform, remote, resources
- **54 Python files** total
- **Root-level CLI commands** for easy access
- **Legacy compatibility wrappers** for smooth migration
```

**Purpose**: Job history tracking, computational recipe storage, and execution records.

**Dependencies**:
- Utilities: `resources` (for efficiency calculations)

**Public API**:
```python
# qxub/history/__init__.py
from .manager import history_manager
from .base import HistoryEntry, ExecutionRecord

__all__ = ['history_manager', 'HistoryEntry', 'ExecutionRecord']
```

### Future Packages

#### Workflow Adapters Package (`qxub/workflow_adapters/`)
```
workflow_adapters/                 # Future: Phase 3
├── __init__.py
├── base.py                       # Base adapter classes
├── snakemake.py                  # Snakemake resource adapter
├── nextflow.py                   # NextFlow resource adapter
└── cwl.py                        # Common Workflow Language adapter
```

**Purpose**: Convert workflow engine resource specifications to platform-native formats.

#### Cloud Platforms (`qxub/platforms/cloud/`)
```
platforms/cloud/                   # Future: Phase 4
├── __init__.py
├── aws_batch.py                  # AWS Batch platform
├── kubernetes.py                 # Kubernetes Jobs platform
└── gcp_batch.py                  # Google Cloud Batch platform
```

**Purpose**: Cloud platform implementations with infrastructure provisioning.

### Templates and Static Files
```
jobscripts/                        # Job script templates (existing)
├── __init__.py
├── qconda.pbs                    # Conda environment template
├── qdefault.pbs                  # Default execution template
├── qmod.pbs                      # Module loading template
└── qsing.pbs                     # Singularity container template
```

**Purpose**: PBS job script templates for different execution contexts.

## Import Guidelines

### Internal Package Imports
```python
# Good: Relative imports within package
from .manager import ConfigManager
from ..resources import parse_memory_size

# Bad: Absolute imports within package
from qxub.config.manager import ConfigManager
```

### Cross-Package Imports
```python
# Good: Import from package public API
from qxub.config import config_manager
from qxub.platforms import get_platform

# Bad: Import internal modules directly
from qxub.config.manager import ConfigManager
from qxub.platforms.loader import PlatformLoader
```

### Circular Dependency Prevention
```python
# Good: Import at function level when needed
def some_function():
    from qxub.execution import execute_unified
    return execute_unified(...)

# Good: Use TYPE_CHECKING for type hints
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from qxub.platforms import Platform
```

## Package Dependencies Matrix

| Package | config | platforms | execution | scheduling | resources | remote | history |
|---------|--------|-----------|-----------|------------|-----------|--------|---------|
| cli     | ✓      | ✓         | ✓         | ✓          | ✓         | ✓      | ✓       |
| config  | -      | -         | -         | -          | ✓         | -      | -       |
| platforms| ✓     | -         | -         | -          | ✓         | -      | -       |
| execution| ✓     | ✓         | -         | ✓          | ✓         | -      | -       |
| scheduling| -    | -         | -         | -          | ✓         | -      | -       |
| resources| -     | -         | -         | -          | -         | -      | -       |
| remote  | ✓      | ✓         | -         | -          | -         | -      | -       |
| history | -      | -         | -         | -          | ✓         | -      | -       |

## Backwards Compatibility

### Migration Period
During migration, maintain backwards compatibility:

```python
# qxub/__init__.py - Compatibility shims
from .config.manager import config_manager
from .scheduling.pbs import qsub, qdel, job_status

# Deprecated: Will be removed in v4.0
import warnings
def deprecated_import():
    warnings.warn("Importing from qxub root is deprecated. Use qxub.config.config_manager",
                  DeprecationWarning, stacklevel=2)
```

### Public API Stability
- Package public APIs (exported in `__init__.py`) are stable
- Internal module organization may change
- Function signatures remain stable across minor versions

## Testing Structure

### Package-Level Tests
```
tests/
├── test_cli/                     # CLI package tests
├── test_config/                  # Configuration tests
├── test_platforms/               # Platform abstraction tests
├── test_execution/               # Execution logic tests
├── test_scheduling/              # Scheduler interaction tests
├── test_resources/               # Resource utilities tests
├── test_remote/                  # Remote execution tests
├── test_history/                 # History management tests
└── integration/                  # Cross-package integration tests
```

### Import Testing
Each package should have tests that verify:
1. Public API imports work correctly
2. No circular dependencies exist
3. Dependency hierarchy is respected

## Future Extensions

### Plugin Architecture Support
The package structure naturally supports plugins:

```python
# Future: Plugin discovery
def discover_workflow_plugins():
    """Discover installed workflow adapter plugins."""
    return find_plugins('qxub.workflow_adapters')

def discover_platform_plugins():
    """Discover installed platform plugins."""
    return find_plugins('qxub.platforms')
```

### Modular Installation
Future versions could support modular installation:

```bash
pip install qxub[core]              # Basic PBS functionality
pip install qxub[remote]            # Add remote execution
pip install qxub[workflows]         # Add workflow engine support
pip install qxub[cloud]             # Add cloud platform support
```

This package structure provides a solid foundation for current functionality while enabling the planned multi-platform and workflow engine capabilities.
