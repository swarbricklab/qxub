# qxub Migration Testing Strategy

## Overview

This document defines a comprehensive testing strategy to ensure the package structure migration maintains 100% functionality and introduces no regressions. The testing strategy covers unit tests, integration tests, and end-to-end validation for each migration phase.

## Testing Philosophy

### Core Principles
1. **Safety First** - Catch any breaking changes immediately
2. **Comprehensive Coverage** - Test all functionality, not just happy paths
3. **Automated Validation** - Reduce manual testing burden
4. **Performance Monitoring** - Ensure no performance degradation
5. **User Experience** - Validate that users see no changes

### Testing Pyramid
```
                    E2E Tests (10%)
                Integration Tests (20%)
              Unit Tests (70%)
```

## Test Suite Organization

### Current Test Structure Analysis
```bash
# Analyze existing tests
find tests/ -name "*.py" -o -name "*.sh" | wc -l
find tests/ -name "*.py" -exec grep -l "def test_" {} \; | wc -l
```

**Current Tests**:
- Shell scripts: `tests/test_*.sh` (integration-style tests)
- Python tests: `tests/test_*.py` (unit tests)
- Platform tests: `tests/run_platform_tests.py`

### Target Test Structure
```
tests/
├── unit/                          # Unit tests for individual packages
│   ├── test_config/
│   ├── test_platforms/
│   ├── test_execution/
│   ├── test_scheduling/
│   ├── test_resources/
│   ├── test_remote/
│   ├── test_history/
│   └── test_cli/
├── integration/                   # Cross-package integration tests
│   ├── test_config_platform_integration.py
│   ├── test_execution_scheduling_integration.py
│   └── test_cli_backend_integration.py
├── migration/                     # Migration-specific tests
│   ├── test_import_compatibility.py
│   ├── test_api_stability.py
│   └── test_circular_dependencies.py
├── end_to_end/                    # Complete workflow tests
│   ├── test_conda_execution.py
│   ├── test_platform_detection.py
│   ├── test_remote_execution.py
│   └── test_workflow_compatibility.py
├── performance/                   # Performance regression tests
│   ├── test_import_performance.py
│   ├── test_cli_startup_time.py
│   └── test_job_submission_speed.py
└── legacy/                        # Existing test files (preserved)
    ├── test_conda_dry.sh
    ├── test_unified_cli.sh
    └── run_platform_tests.py
```

## Phase-Specific Test Plans

### Pre-Migration Baseline Tests

#### Baseline Test Suite (`tests/pre_migration_baseline.py`)
```python
"""
Comprehensive baseline test suite to run before migration begins.
Creates a reference point for all functionality.
"""

import pytest
import subprocess
import time
from pathlib import Path

class TestPreMigrationBaseline:
    """Baseline tests to establish reference behavior."""

    def test_all_cli_commands_exist(self):
        """Verify all CLI commands are accessible."""
        commands = [
            "qxub --help",
            "qxub config --help",
            "qxub exec --help",
            "qxub platform --help",
            "qxub history --help",
            "qxub monitor --help",
            "qxub cancel --help"
        ]
        for cmd in commands:
            result = subprocess.run(cmd.split(), capture_output=True, text=True)
            assert result.returncode == 0, f"Command failed: {cmd}"

    def test_import_performance_baseline(self):
        """Measure current import performance as baseline."""
        start_time = time.time()
        import qxub
        import_time = time.time() - start_time

        # Store baseline (should be < 0.5 seconds)
        assert import_time < 0.5, f"Import too slow: {import_time}s"

        # Save baseline for comparison
        with open("tests/.baseline_import_time", "w") as f:
            f.write(str(import_time))

    def test_all_imports_work(self):
        """Verify all current imports work."""
        import_tests = [
            "from qxub.config_manager import config_manager",
            "from qxub.platform import get_platform",
            "from qxub.scheduler import qsub, qdel",
            "from qxub.resource_utils import parse_memory_size",
            "from qxub.execution import execute_unified"
        ]

        for import_test in import_tests:
            exec(import_test)  # Should not raise any exceptions

    def test_dry_run_execution(self):
        """Test basic execution functionality."""
        cmd = ["qxub", "exec", "--dry", "--", "echo", "test"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        assert "qsub" in result.stdout

    def test_configuration_access(self):
        """Test configuration system."""
        cmd = ["qxub", "config", "get", "defaults.project"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        # Should not crash, regardless of output
        assert result.returncode in [0, 1]  # 0 if set, 1 if not set

    def test_platform_detection(self):
        """Test platform detection."""
        cmd = ["qxub", "platform", "list"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
```

### Phase 1 Tests: Resources and History

#### Resources Package Tests (`tests/unit/test_resources/`)
```python
# tests/unit/test_resources/test_parser.py
import pytest
from qxub.resources import parse_memory_size, parse_walltime

class TestResourceParser:
    def test_memory_parsing(self):
        """Test memory size parsing."""
        assert parse_memory_size("4GB") == 4 * 1024**3
        assert parse_memory_size("1000MB") == 1000 * 1024**2

    def test_walltime_parsing(self):
        """Test walltime parsing."""
        assert parse_walltime("2:00:00") == 7200
        assert parse_walltime("1h") == 3600

# tests/unit/test_resources/test_import_compatibility.py
def test_backwards_compatible_imports():
    """Ensure old imports still work during migration."""
    # These should work via compatibility layer
    from qxub.resource_parser import parse_memory_size
    from qxub.resource_utils import format_walltime
    from qxub.resource_tracker import resource_tracker
```

#### History Package Tests (`tests/unit/test_history/`)
```python
# tests/unit/test_history/test_manager.py
import pytest
from qxub.history import history_manager

class TestHistoryManager:
    def test_history_manager_accessible(self):
        """Test history manager is accessible."""
        assert history_manager is not None

    def test_backwards_compatible_imports(self):
        """Test old imports still work."""
        from qxub.history_manager import history_manager as old_import
        assert old_import is not None
```

#### Phase 1 Integration Test
```python
# tests/migration/test_phase1_integration.py
def test_phase1_no_regressions():
    """Comprehensive test that Phase 1 introduces no regressions."""

    # Test all CLI commands still work
    test_all_cli_commands()

    # Test resource functionality
    test_resource_parsing_functionality()

    # Test history functionality
    test_history_functionality()

    # Test performance hasn't degraded
    test_import_performance_unchanged()
```

### Phase 2 Tests: Configuration

#### Configuration Package Tests (`tests/unit/test_config/`)
```python
# tests/unit/test_config/test_manager.py
import pytest
from qxub.config import config_manager

class TestConfigManager:
    def test_config_manager_accessible(self):
        """Test config manager works after move."""
        assert config_manager is not None

    def test_config_loading(self):
        """Test configuration loading."""
        # Should not crash
        config_manager.get_default_platform()

    def test_backwards_compatible_imports(self):
        """Test old imports work via compatibility layer."""
        from qxub.config_manager import config_manager as old_import
        assert old_import is not None
```

#### Configuration Integration Tests
```python
# tests/integration/test_config_integration.py
def test_config_platform_integration():
    """Test config works with platform system."""
    from qxub.config import config_manager
    from qxub.platforms import get_platform

    platform_name = config_manager.get_default_platform()
    if platform_name:
        platform = get_platform(platform_name)
        assert platform is not None
```

### Phase 3 Tests: CLI Package

#### CLI Package Tests (`tests/unit/test_cli/`)
```python
# tests/unit/test_cli/test_exec_cli.py
import subprocess
import pytest

class TestExecCLI:
    def test_exec_command_help(self):
        """Test exec command help works."""
        result = subprocess.run(["qxub", "exec", "--help"],
                              capture_output=True, text=True)
        assert result.returncode == 0
        assert "usage:" in result.stdout.lower()

    def test_exec_dry_run(self):
        """Test exec dry run functionality."""
        cmd = ["qxub", "exec", "--dry", "--", "echo", "test"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0

# tests/unit/test_cli/test_all_commands.py
class TestAllCLICommands:
    """Test all CLI commands work after reorganization."""

    CLI_COMMANDS = [
        "config", "exec", "history", "monitor",
        "cancel", "status", "platform", "alias"
    ]

    @pytest.mark.parametrize("command", CLI_COMMANDS)
    def test_command_help(self, command):
        """Test each CLI command help works."""
        result = subprocess.run(["qxub", command, "--help"],
                              capture_output=True, text=True)
        assert result.returncode == 0
```

### Phase 4 Tests: Platform Package

#### Platform Package Tests (`tests/unit/test_platforms/`)
```python
# tests/unit/test_platforms/test_loader.py
import pytest
from qxub.platforms import get_platform, list_platforms

class TestPlatformLoader:
    def test_platform_loading(self):
        """Test platform loading functionality."""
        platforms = list_platforms()
        assert isinstance(platforms, list)

    def test_platform_detection(self):
        """Test platform detection."""
        from qxub.platforms import detect_platform
        # Should not crash
        detected = detect_platform()

    def test_backwards_compatible_imports(self):
        """Test old platform imports work."""
        from qxub.platform import get_platform as old_import
        assert callable(old_import)
```

#### Platform Integration Tests
```python
# tests/integration/test_platform_integration.py
def test_platform_cli_integration():
    """Test platform CLI works with new package structure."""
    result = subprocess.run(["qxub", "platform", "list"],
                          capture_output=True, text=True)
    assert result.returncode == 0
```

### Phase 5 Tests: Execution and Scheduling

#### Execution Package Tests (`tests/unit/test_execution/`)
```python
# tests/unit/test_execution/test_unified.py
import pytest
from qxub.execution import execute_unified, ExecutionContext

class TestUnifiedExecution:
    def test_execution_context_creation(self):
        """Test execution context creation."""
        ctx = ExecutionContext("conda", "test_env", "conda")
        assert ctx.type == "conda"
        assert ctx.value == "test_env"

    def test_backwards_compatible_imports(self):
        """Test old execution imports work."""
        from qxub.execution import execute_unified as old_import
        assert callable(old_import)
```

#### Scheduling Package Tests (`tests/unit/test_scheduling/`)
```python
# tests/unit/test_scheduling/test_pbs.py
import pytest
from qxub.scheduling import qsub, qdel, job_status

class TestPBSScheduling:
    def test_pbs_functions_accessible(self):
        """Test PBS functions are accessible."""
        assert callable(qsub)
        assert callable(qdel)
        assert callable(job_status)

    def test_backwards_compatible_imports(self):
        """Test old scheduler imports work."""
        from qxub.scheduler import qsub as old_import
        assert callable(old_import)
```

### Phase 6 Tests: Remote Package

#### Remote Package Tests (`tests/unit/test_remote/`)
```python
# tests/unit/test_remote/test_config.py
import pytest
from qxub.remote import RemoteConfig

class TestRemoteConfig:
    def test_remote_config_creation(self):
        """Test remote config creation."""
        config = RemoteConfig(
            name="test",
            host="test.example.com",
            platform="test_platform"
        )
        assert config.name == "test"
```

## Comprehensive Migration Test Suite

### Import Compatibility Test (`tests/migration/test_import_compatibility.py`)
```python
"""
Comprehensive test to ensure all imports work after migration.
This test should pass at every phase of migration.
"""

import pytest

class TestImportCompatibility:
    """Test that all public APIs remain accessible."""

    def test_config_imports(self):
        """Test configuration imports work."""
        # New imports
        from qxub.config import config_manager

        # Old imports (via compatibility layer)
        from qxub.config_manager import config_manager as old_manager

        assert config_manager is old_manager

    def test_platform_imports(self):
        """Test platform imports work."""
        # New imports
        from qxub.platforms import get_platform

        # Old imports (via compatibility layer)
        from qxub.platform import get_platform as old_get_platform

        assert get_platform is old_get_platform

    def test_scheduling_imports(self):
        """Test scheduling imports work."""
        # New imports
        from qxub.scheduling import qsub

        # Old imports (via compatibility layer)
        from qxub.scheduler import qsub as old_qsub

        assert qsub is old_qsub

    def test_resources_imports(self):
        """Test resource imports work."""
        # New imports
        from qxub.resources import parse_memory_size

        # Old imports (via compatibility layer)
        from qxub.resource_parser import parse_memory_size as old_parse

        assert parse_memory_size is old_parse
```

### API Stability Test (`tests/migration/test_api_stability.py`)
```python
"""
Test that public APIs maintain their signatures and behavior.
"""

import inspect
import pytest

class TestAPIStability:
    """Ensure public API signatures don't change."""

    def test_qsub_signature_unchanged(self):
        """Test qsub function signature is unchanged."""
        from qxub.scheduling import qsub

        sig = inspect.signature(qsub)
        params = list(sig.parameters.keys())

        # Expected signature: qsub(cmd, quiet=False)
        assert params == ["cmd", "quiet"]
        assert sig.parameters["quiet"].default is False

    def test_config_manager_interface_unchanged(self):
        """Test config manager interface is unchanged."""
        from qxub.config import config_manager

        # Test expected methods exist
        expected_methods = [
            "get_default_platform",
            "get_queue_preferences",
            "get_platform_search_paths"
        ]

        for method in expected_methods:
            assert hasattr(config_manager, method)
            assert callable(getattr(config_manager, method))
```

### Circular Dependency Test (`tests/migration/test_circular_dependencies.py`)
```python
"""
Test for circular dependencies in the new package structure.
"""

import ast
import pytest
from pathlib import Path

class TestCircularDependencies:
    """Detect circular dependencies between packages."""

    def test_no_circular_imports(self):
        """Test that packages don't have circular imports."""
        # Parse all Python files and build dependency graph
        qxub_path = Path("qxub")
        dependencies = {}

        for py_file in qxub_path.rglob("*.py"):
            if py_file.name == "__init__.py":
                continue

            with open(py_file) as f:
                tree = ast.parse(f.read())

            imports = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        if alias.name.startswith("qxub."):
                            imports.append(alias.name)
                elif isinstance(node, ast.ImportFrom):
                    if node.module and node.module.startswith("qxub."):
                        imports.append(node.module)

            package = str(py_file.relative_to(qxub_path)).replace("/", ".").replace(".py", "")
            dependencies[package] = imports

        # Check for cycles (simplified check)
        # In production, use proper cycle detection algorithm
        for package, deps in dependencies.items():
            for dep in deps:
                if package in dependencies.get(dep, []):
                    pytest.fail(f"Circular dependency detected: {package} <-> {dep}")
```

### Performance Regression Test (`tests/performance/test_performance_regression.py`)
```python
"""
Test that migration doesn't introduce performance regressions.
"""

import time
import pytest

class TestPerformanceRegression:
    """Test performance hasn't degraded after migration."""

    def test_import_performance(self):
        """Test import performance hasn't degraded."""
        # Load baseline
        try:
            with open("tests/.baseline_import_time") as f:
                baseline = float(f.read().strip())
        except FileNotFoundError:
            pytest.skip("No baseline import time available")

        # Measure current import time
        start_time = time.time()
        import qxub
        import_time = time.time() - start_time

        # Allow 20% degradation
        max_allowed = baseline * 1.2
        assert import_time <= max_allowed, \
            f"Import time degraded: {import_time}s > {max_allowed}s (baseline: {baseline}s)"

    def test_cli_startup_time(self):
        """Test CLI startup time hasn't degraded."""
        import subprocess

        start_time = time.time()
        result = subprocess.run(["qxub", "--help"],
                              capture_output=True, text=True)
        startup_time = time.time() - start_time

        assert result.returncode == 0
        assert startup_time < 2.0, f"CLI startup too slow: {startup_time}s"
```

## End-to-End Validation Tests

### Complete Workflow Test (`tests/end_to_end/test_complete_workflows.py`)
```python
"""
End-to-end tests that validate complete qxub workflows work correctly.
"""

import subprocess
import pytest
import tempfile
from pathlib import Path

class TestCompleteWorkflows:
    """Test complete qxub workflows end-to-end."""

    def test_conda_execution_workflow(self):
        """Test complete conda execution workflow."""
        cmd = [
            "qxub", "exec",
            "--dry",
            "--env", "base",
            "--name", "test-job",
            "--",
            "python", "-c", "print('hello world')"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        assert "qsub" in result.stdout
        assert "test-job" in result.stdout

    def test_platform_queue_selection_workflow(self):
        """Test platform detection and queue selection."""
        # Test platform detection
        result = subprocess.run(["qxub", "platform", "list"],
                              capture_output=True, text=True)
        assert result.returncode == 0

        # Test queue selection
        result = subprocess.run([
            "qxub", "select-queue",
            "--cpus", "2",
            "--memory", "4GB",
            "--format", "json"
        ], capture_output=True, text=True)
        assert result.returncode == 0

    def test_configuration_workflow(self):
        """Test configuration management workflow."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_file = Path(tmpdir) / "test_config.yaml"

            # Test config file operations
            result = subprocess.run([
                "qxub", "--config", str(config_file),
                "config", "set", "defaults.project", "test123"
            ], capture_output=True, text=True)
            assert result.returncode == 0

            # Verify setting was saved
            result = subprocess.run([
                "qxub", "--config", str(config_file),
                "config", "get", "defaults.project"
            ], capture_output=True, text=True)
            assert result.returncode == 0
            assert "test123" in result.stdout
```

## Test Execution Strategy

### Automated Test Pipeline
```bash
#!/bin/bash
# tests/run_migration_tests.sh

set -e

echo "🧪 Running Migration Test Suite"

# Phase-specific tests
echo "📦 Testing Phase 1: Resources & History"
python -m pytest tests/unit/test_resources/ tests/unit/test_history/ -v

echo "⚙️  Testing Phase 2: Configuration"
python -m pytest tests/unit/test_config/ -v

echo "🖥️  Testing Phase 3: CLI"
python -m pytest tests/unit/test_cli/ -v

echo "🏗️  Testing Phase 4: Platforms"
python -m pytest tests/unit/test_platforms/ -v

echo "🚀 Testing Phase 5: Execution & Scheduling"
python -m pytest tests/unit/test_execution/ tests/unit/test_scheduling/ -v

echo "🌐 Testing Phase 6: Remote"
python -m pytest tests/unit/test_remote/ -v

# Integration tests
echo "🔗 Testing Integration"
python -m pytest tests/integration/ -v

# Migration-specific tests
echo "📦 Testing Migration Compatibility"
python -m pytest tests/migration/ -v

# Performance tests
echo "⚡ Testing Performance"
python -m pytest tests/performance/ -v

# End-to-end tests
echo "🎯 Testing End-to-End Workflows"
python -m pytest tests/end_to_end/ -v

# Legacy tests (existing test suite)
echo "🔄 Running Legacy Test Suite"
./tests/test_conda_dry.sh
./tests/test_unified_cli.sh
python tests/run_platform_tests.py

echo "✅ All Migration Tests Passed!"
```

### Pre-Commit Hooks
```yaml
# .pre-commit-config.yaml (add to existing)
repos:
  - repo: local
    hooks:
      - id: migration-tests
        name: Migration Tests
        entry: python -m pytest tests/migration/ -x
        language: system
        always_run: true

      - id: import-tests
        name: Import Tests
        entry: python -c "import qxub; from qxub.config import config_manager; print('✓ Imports work')"
        language: system
        always_run: true
```

### CI/CD Integration
```yaml
# .github/workflows/migration-tests.yml
name: Migration Tests

on: [push, pull_request]

jobs:
  migration-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.10

      - name: Install dependencies
        run: |
          pip install -e .
          pip install pytest

      - name: Run migration test suite
        run: |
          ./tests/run_migration_tests.sh
```

## Success Criteria

### Test Completion Criteria
Each phase passes when:
1. ✅ All unit tests pass for affected packages
2. ✅ All integration tests pass
3. ✅ Migration compatibility tests pass
4. ✅ Performance tests show no regression
5. ✅ End-to-end workflows work correctly

### Overall Migration Test Success
Migration testing is successful when:
1. ✅ 100% backwards compatibility maintained
2. ✅ All existing functionality preserved
3. ✅ No performance degradation
4. ✅ Clean import hierarchy established
5. ✅ No circular dependencies
6. ✅ All CLI commands work unchanged
7. ✅ All configuration functionality works
8. ✅ All execution contexts work (conda, modules, singularity)
9. ✅ Platform detection and queue selection work
10. ✅ Remote execution functionality preserved

This comprehensive testing strategy ensures that the migration maintains the stability and reliability that qxub users depend on while enabling the future multi-platform architecture.
