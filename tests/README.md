# qxub Testing Suite

This directory contains comprehensive integration tests and testing utilities for qxub, designed to catch bugs and validate functionality across different configurations and HPC environments.

## Test Scripts

### Integration Tests

### 1. `test_unified_cli.sh` - Unified CLI Tests ⚡

**Purpose**: Validation of the new unified CLI interface (qxub 2.0+).

### 2. `test_conda_dry.sh` - Quick Dry-Run Tests ⚡

**Purpose**: Fast validation of command generation and parsing without submitting actual jobs.

**Usage**:
```bash
./test_conda_dry.sh
```

**Features**:
- ✅ **Safe**: Uses `--dry-run` mode, no jobs submitted
- ⚡ **Fast**: Completes in seconds
- 🔍 **Comprehensive**: Tests all CLI options and configurations
- 🎯 **Development-friendly**: Perfect for rapid iteration

**Test Coverage**:
- All main PBS options (`--name`, `--project`, `--queue`, etc.)
- Resource specifications (`--resources mem=8GB`, etc.)
- Conda-specific options (`--env`, `--pre`, `--post`, `--template`)
- Configuration file integration
- Alias functionality
- Complex command scenarios
- Error condition handling

### 2. `test_conda_integration.sh` - Full Integration Tests 🚀

**Purpose**: Complete end-to-end testing including actual job submission and monitoring.

**Usage**:
```bash
# Standard run (submits real jobs - use carefully!)
./test_conda_integration.sh

# Check what it would do first
head -50 test_conda_integration.sh
```

**⚠️ Warning**: This script submits real jobs to the PBS queue. Use with caution!

**Features**:
- 🎯 **Real-world testing**: Actual job submission and monitoring
- 🔄 **Complete workflows**: Tests full lifecycle including cleanup
- 🛡️ **Safe cleanup**: Handles interruptions and job cleanup
- 📊 **Detailed reporting**: Comprehensive test results

**Test Coverage**:
- Basic conda functionality with job submission
- Config file integration with real jobs
- Alias execution with actual PBS jobs
- Complex command scenarios
- Resource allocation verification
- Error handling and cleanup

### 3. `test_remote_config.sh` - Remote Execution Tests 🌐

**Purpose**: Testing qxub v2.2 remote execution capabilities and `--config` option functionality.

**Usage**:
```bash
./test_remote_config.sh
```

**Features**:
- ✅ **Safe**: Uses dry-run mode for remote testing
- 🔧 **Configuration testing**: Validates `--config` option with alternative config files
- 🌐 **Remote validation**: Tests remote configuration loading and validation
- 🎯 **Error handling**: Verifies proper error messages for connection failures

**Test Coverage**:
- `--config` option recognition and parsing
- Alternative configuration file loading (`tests/test_config.yaml`)
- Remote configuration validation
- SSH connection error handling
- Invalid remote name rejection

**Test Configuration File**: `tests/test_config.yaml` contains example remote configurations:
- `nci_gadi`: NCI Gadi cluster configuration
- `localhost`: Local testing configuration
- `example_cluster`: Generic cluster example
- `custom_port_host`: Custom port example

## Test Configuration

Both test scripts create temporary configuration files to ensure consistent testing:

```yaml
# Example test config created by the scripts
defaults:
  name: "test-job-{timestamp}"
  time: "00:30:00"
  mem: "4GB"
  ncpus: 1
  project: "a56"
  queue: "normal"

aliases:
  quick:
    subcommand: conda
    env: base
    name: "quick-{user}"
    time: "00:15:00"
    mem: "2GB"

  bigmem:
    subcommand: conda
    env: scipy
    name: "bigmem-job"
    time: "01:00:00"
    mem: "32GB"
    ncpus: 4
    queue: "hugemem"
```

## Development Workflow

### 1. **Quick Validation** (Recommended for development)
```bash
# Run dry tests to check basic functionality
./test_conda_dry.sh

# Check specific test failure details
echo $?  # Non-zero indicates failures
run_dry_test "Test description" \
  "qxub --env base -- your-command-here"
### 2. **Full Testing** (Before releases)
```bash
# Run comprehensive tests (careful - submits jobs!)
./test_conda_integration.sh

# Monitor job queue during tests
watch qstat -u $USER
```

### 3. **Debugging Failed Tests**
run_test "Test description" \
  "qxub --env base -- your-command-here" \
  0  # Expected exit code

# Check configuration
qxub config get

# Validate specific aliases
qxub config alias test quick
```

## HPC-Specific Considerations

These tests are designed to catch common issues in HPC environments:

### Configuration Issues
- ✅ XDG config directory handling
- ✅ Template variable resolution (`{user}`, `{project}`, `{timestamp}`)
- ✅ Default value inheritance
- ✅ Config file validation

### Resource Management
- ✅ Memory specification formats (`4GB`, `4096MB`)
- ✅ Time format variations (`01:30:00`, `1:30`, `90`)
- ✅ CPU and GPU resource allocation
- ✅ Queue-specific requirements

### Environment Handling
- ✅ Conda environment activation
- ✅ Module loading integration
- ✅ Path and variable resolution
- ✅ Cross-node compatibility

### Command Execution
- ✅ Complex shell commands with pipes/redirects
- ✅ Quoting and escaping
- ✅ Multi-line commands
- ✅ Pre/post command execution

## Adding New Tests

### For Dry-Run Tests (`test_conda_dry.sh`)
```bash
# Add to appropriate test function
run_dry_test "Test description" \
    "qxub --env base -- your-command-here"

# For error conditions (should fail)
run_dry_test "Error test description" \
    "qxub --env base -- invalid-command" \
    2  # Expected exit code
```

### For Integration Tests (`test_conda_integration.sh`)
```bash
# Add to appropriate test function
run_test "Test description" \
    "qxub --env base -- your-command-here" \
    0  # Expected exit code
```

## Test Output

### Success Output
```
[INFO] Starting qxub --env dry-run integration tests...
[TEST] All main PBS options
Command: qxub --name cli-test --project a56 --queue normal --env base -- echo 'All options test' --dry-run
[PASS] All main PBS options
...
========================================
         DRY RUN TEST RESULTS
========================================
Total Tests: 31
Passed: 31
Failed: 0

🎉 All dry-run tests passed!
```

### Failure Output
```
[FAIL] Missing environment (expected exit code 2, got 0)
...
========================================
         DRY RUN TEST RESULTS
========================================
Total Tests: 31
Passed: 30
Failed: 1

Failed Test Details:
Missing environment: Expected exit 2, got 0
❌ Some tests failed!
```
[INFO] Starting qxub dry-run integration tests...
Command: qxub --name cli-test --project a56 --queue normal --env base -- echo 'All options test' --dry-run

These tests can be integrated into CI/CD pipelines:

```yaml
# Example GitHub Actions integration
- name: Run qxub dry tests
  run: ./test_conda_dry.sh

# For full integration tests (if PBS available)
- name: Run integration tests
  run: ./test_conda_integration.sh
  if: ${{ github.event_name == 'release' }}
```

## Troubleshooting

### Common Issues

1. **Config file conflicts**
   - Tests backup and restore `~/.config/qxub/config.yaml`
   - Manual cleanup: `rm ~/.config/qxub/config.yaml.backup`

2. **PBS environment not available**
   - Dry tests work anywhere
   - Integration tests require PBS environment

3. **Permission issues**
   - Ensure scripts are executable: `chmod +x test_*.sh`
   - Check directory permissions for config creation

4. **Test failures due to HPC quirks**
   - Check module availability (`module avail`)
   - Verify queue names and limits (`qstat -Q`)
   - Confirm project codes (`id`)

### Debug Mode

Enable verbose output for debugging:
```bash
# Add debug logging to test scripts
export QXUB_DEBUG=1
./test_conda_dry.sh

# Or modify tests to use verbose flags
qxub -vvv conda --env base --dry-run echo "Debug"
```
