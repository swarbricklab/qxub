# qxub Package Structure Specification

## Overview

This document defines the target package structure for qxub, designed to support current functionality while enabling future multi-platform and workflow engine capabilities.

## Design Principles

### Organizational Principles
1. **Single Responsibility** - Each package has one clear purpose
2. **Dependency Clarity** - Clear import hierarchy with minimal circular dependencies
3. **Future-Proof** - Structure accommodates planned features without refactoring
4. **Logical Grouping** - Related functionality lives together
5. **Public APIs** - Clean interfaces between packages

### Dependency Hierarchy
```
CLI Layer (qxub.cli.*)
    ↓
Business Logic (qxub.execution.*, qxub.platforms.*)
    ↓
Core Services (qxub.config.*, qxub.scheduling.*)
    ↓
Utilities (qxub.resources.*, qxub.remote.*)
```

## Target Package Structure

### Root Level
```
qxub/
├── __init__.py                    # Main package exports and version
├── cli.py                         # Main CLI entry point
└── cli_old.py                     # Legacy CLI (temporary, will be removed)
```

**Purpose**: Package entry points and backwards compatibility during migration.

**Public API**:
```python
# qxub/__init__.py
__version__ = "3.0.0"

# Main CLI entry point
from .cli import qxub

# Core public APIs
from .config import config_manager
from .platforms import get_platform, select_best_queue
from .execution import execute_unified
from .scheduling import qsub, qdel, job_status
```

### CLI Package (`qxub/cli/`)
```
cli/
├── __init__.py                    # CLI module exports
├── exec_cli.py                   # Main execution command (qxub exec)
├── config_cli.py                # Configuration commands (qxub config)
├── alias_cli.py                  # Alias management (qxub alias)
├── history_cli.py                # Job history (qxub history)
├── monitor_cli.py                # Job monitoring (qxub monitor)
├── cancel_cli.py                 # Job cancellation (qxub cancel)
├── status_cli.py                 # Job status (qxub status)
├── resources_cli.py              # Resource management (qxub resources)
└── platform_cli.py               # Platform commands (qxub platform)
```

**Purpose**: All Click-based command-line interface definitions.

**Dependencies**:
- Business Logic: `execution`, `platforms`, `config`
- Core Services: `scheduling`, `resources`

**Public API**:
```python
# qxub/cli/__init__.py
from .exec_cli import exec_cli
from .config_cli import config_cli
from .platform_cli import platform_cli
# ... other CLI commands

__all__ = ['exec_cli', 'config_cli', 'platform_cli', ...]
```

### Configuration Package (`qxub/config/`)
```
config/
├── __init__.py                   # Config module exports
├── manager.py                    # Main configuration manager
├── handler.py                    # Configuration processing logic
├── base.py                       # Base configuration classes
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
from .shortcuts import shortcut_manager

__all__ = ['config_manager', 'ConfigManager', 'ConfigurationError', 'shortcut_manager']
```

**Key Classes**:
- `ConfigManager`: Main configuration interface
- `ConfigurationError`: Configuration-related exceptions
- `shortcut_manager`: Global shortcut management instance

### Platform Package (`qxub/platforms/`)
```
platforms/
├── __init__.py                   # Platform module exports
├── base.py                       # Abstract platform interfaces
├── loader.py                     # Platform discovery and loading
├── pbs_pro.py                    # PBS Pro platform implementation
├── detection.py                  # Platform auto-detection logic
└── integration.py                # Platform integration utilities
```

**Purpose**: Platform abstraction, queue selection, resource validation, and scheduler-agnostic interfaces.

**Dependencies**:
- Core Services: `config` (for platform search paths)
- Utilities: `resources` (for resource parsing and validation)

**Public API**:
```python
# qxub/platforms/__init__.py
from .base import Platform, Queue, QueueLimits, PlatformResources
from .loader import get_platform, list_platforms, get_current_platform
from .detection import detect_platform
from .pbs_pro import PBSPlatform, PBSResources

__all__ = [
    'Platform', 'Queue', 'QueueLimits', 'PlatformResources',
    'get_platform', 'list_platforms', 'get_current_platform', 'detect_platform',
    'PBSPlatform', 'PBSResources'
]
```

**Key Classes**:
- `Platform`: Abstract platform definition
- `PlatformResources`: Abstract resource specification
- `PBSPlatform`: PBS Pro implementation
- `PlatformLoader`: Platform discovery and loading

### Execution Package (`qxub/execution/`)
```
execution/
├── __init__.py                   # Execution module exports
├── context.py                    # Execution contexts (conda, modules, singularity)
├── unified.py                    # Unified execution logic
├── executors.py                  # Individual executor implementations
└── templates.py                  # Job script template management
```

**Purpose**: Job execution orchestration, environment setup, and job script generation.

**Dependencies**:
- Core Services: `config`, `scheduling`
- Business Logic: `platforms`
- Utilities: `resources`

**Public API**:
```python
# qxub/execution/__init__.py
from .unified import execute_unified
from .context import ExecutionContext
from .templates import get_template, list_templates

__all__ = ['execute_unified', 'ExecutionContext', 'get_template', 'list_templates']
```

**Key Classes**:
- `ExecutionContext`: Represents execution environment (conda, module, etc.)
- `execute_unified`: Main execution orchestration function

### Scheduling Package (`qxub/scheduling/`)
```
scheduling/
├── __init__.py                   # Scheduling module exports
├── scheduler.py                  # Abstract scheduler interface
├── pbs.py                        # PBS Pro scheduler implementation
└── monitoring.py                 # Job monitoring utilities
```

**Purpose**: Job scheduler interaction, job submission, status checking, and monitoring.

**Dependencies**:
- Utilities: `resources` (for resource parsing)
- No dependencies on higher-level packages

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

### Remote Package (`qxub/remote/`)
```
remote/
├── __init__.py                   # Remote module exports
├── config.py                     # Remote configuration classes
├── loader.py                     # Remote configuration loading
├── executor.py                   # Remote execution implementation
└── integration.py                # Platform integration for remote execution
```

**Purpose**: SSH-based remote execution, remote platform management, and remote configuration.

**Dependencies**:
- Business Logic: `platforms`
- Core Services: `config`

**Public API**:
```python
# qxub/remote/__init__.py
from .config import RemoteConfig
from .loader import load_remote_configurations, get_remote_config
from .executor import RemoteExecutorFactory

__all__ = ['RemoteConfig', 'load_remote_configurations', 'get_remote_config',
           'RemoteExecutorFactory']
```

### History Package (`qxub/history/`)
```
history/
├── __init__.py                   # History module exports
├── manager.py                    # History management
└── base.py                       # History base classes and data models
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
