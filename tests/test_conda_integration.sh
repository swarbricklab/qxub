#!/bin/bash
#
# Integration test for qxub conda subcommand.
#
# Tests all combinations of options and calling patterns for the conda subcommand,
# including config files, aliases, and CLI overrides. Designed to catch HPC-specific
# issues and edge cases in the new config and alias system.
#

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Test results array
declare -a FAILED_TEST_DETAILS

# Helper functions
log_info() {
    echo -e "${YELLOW}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

log_error() {
    echo -e "${RED}[FAIL]${NC} $1"
}

run_test() {
    local test_name="$1"
    local test_cmd="$2"
    local expected_exit_code="${3:-0}"

    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    log_info "Running test: $test_name"
    echo "Command: $test_cmd"

    # Run the command and capture output and exit code
    if output=$(eval "$test_cmd" 2>&1); then
        actual_exit_code=0
    else
        actual_exit_code=$?
    fi

    # Check if exit code matches expected
    if [ "$actual_exit_code" -eq "$expected_exit_code" ]; then
        log_success "$test_name"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        echo "Output: $output"
    else
        log_error "$test_name (expected exit code $expected_exit_code, got $actual_exit_code)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        FAILED_TEST_DETAILS+=("$test_name: Expected exit $expected_exit_code, got $actual_exit_code")
        echo "Output: $output"
    fi
    echo "---"
}

# Setup test environment
setup_test_env() {
    log_info "Setting up test environment..."

    # Backup existing config if it exists
    if [ -f ~/.config/qxub/config.yaml ]; then
        cp ~/.config/qxub/config.yaml ~/.config/qxub/config.yaml.backup
        log_info "Backed up existing config"
    fi

    # Create test config directory
    mkdir -p ~/.config/qxub

    # Create a test config file
    cat > ~/.config/qxub/config.yaml << 'EOF'
defaults:
  name: "test-job"
  project: "a56"
  queue: "normal"
  resources:
    - "walltime=00:30:00"
    - "mem=4GB"
    - "ncpus=1"

aliases:
  test-small:
    main:
      name: "small-test"
      resources:
        - "walltime=00:15:00"
        - "mem=2GB"
    subcommand:
      type: conda
      env: base
    target:
      cmd: echo 'Small test alias'

  test-large:
    main:
      name: "large-test"
      resources:
        - "walltime=01:00:00"
        - "mem=8GB"
        - "ncpus=2"
    subcommand:
      type: conda
      env: myenv
    target:
      cmd: echo 'Large test alias'

  test-gpu:
    main:
      name: "gpu-test"
      queue: "gpuvolta"
      resources:
        - "ngpus=1"
        - "ncpus=12"
    subcommand:
      type: conda
      env: tensorflow
    target:
      cmd: echo 'GPU test alias'
EOF

    log_info "Created test config file"
}

# Cleanup test environment
cleanup_test_env() {
    log_info "Cleaning up test environment..."

    # Restore original config if it existed
    if [ -f ~/.config/qxub/config.yaml.backup ]; then
        mv ~/.config/qxub/config.yaml.backup ~/.config/qxub/config.yaml
        log_info "Restored original config"
    else
        rm -f ~/.config/qxub/config.yaml
        rmdir ~/.config/qxub 2>/dev/null || true
    fi
}

# Test basic conda functionality
test_basic_conda() {
    log_info "=== Testing Basic Conda Functionality ==="

    # Test 1: Basic conda command with minimal options
    run_test "Basic conda command" \
        "qxub conda --env base echo 'Hello World'" \
        0

    # Test 2: Conda with all basic PBS options
    run_test "Conda with full PBS options" \
        "qxub --name test-job --resources walltime=00:15:00 --resources mem=2GB --resources ncpus=1 --project a56 --queue normal conda --env base echo 'Full options'" \
        0

    # Test 3: Conda with storage option (using resources)
    run_test "Conda with storage" \
        "qxub --resources storage='scratch/a56+gdata/a56' conda --env base echo 'Storage test'" \
        0

    # Test 4: Conda with walltime variations
    run_test "Conda with hours:minutes format" \
        "qxub --resources walltime=1:30 conda --env base echo 'Time format test'" \
        0

    # Test 5: Conda with memory variations
    run_test "Conda with MB memory" \
        "qxub --resources mem=2048MB conda --env base echo 'Memory format test'" \
        0
}

# Test config integration
test_config_integration() {
    log_info "=== Testing Config Integration ==="

    # Test 6: Using defaults from config
    run_test "Config defaults" \
        "qxub conda --env base echo 'Config defaults test'" \
        0

    # Test 7: Overriding config defaults
    run_test "Override config defaults" \
        "qxub --name override-test --resources walltime=00:45:00 conda --env base echo 'Override test'" \
        0

    # Test 8: Partial override (some from config, some from CLI)
    run_test "Partial config override" \
        "qxub --resources mem=6GB conda --env base echo 'Partial override test'" \
        0
}

# Test alias functionality
test_alias_functionality() {
    log_info "=== Testing Alias Functionality ==="

    # Test 9: Basic alias execution
    run_test "Basic alias execution" \
        "qxub alias test-small echo 'Alias test'" \
        0

    # Test 10: Alias with CLI overrides
    run_test "Alias with overrides" \
        "qxub alias test-small --resources walltime=00:20:00 --resources mem=3GB echo 'Alias override test'" \
        0

    # Test 11: Large job alias
    run_test "Large job alias" \
        "qxub alias test-large echo 'Large job test'" \
        0

    # Test 12: GPU alias
    run_test "GPU alias" \
        "qxub alias test-gpu echo 'GPU test'" \
        0

    # Test 13: Alias with additional arguments
    run_test "Alias with extra arguments" \
        "qxub alias test-small 'python -c \"print(\\\"Python test\\\")\"'" \
        0
}

# Test edge cases and error conditions
test_edge_cases() {
    log_info "=== Testing Edge Cases ==="

    # Test 14: Missing environment (should fail when no conda env available)
    run_test "Missing environment" \
        "CONDA_DEFAULT_ENV= qxub conda echo 'No env specified'" \
        2

    # Test 15: Invalid environment name
    run_test "Invalid environment" \
        "qxub conda --env nonexistent-env echo 'Invalid env'" \
        0  # Should succeed but might warn

    # Test 16: Invalid time format (should fail)
    run_test "Invalid time format" \
        "qxub --resources walltime=invalid-time conda --env base echo 'Bad time'" \
        2

    # Test 17: Invalid memory format (should fail)
    run_test "Invalid memory format" \
        "qxub --resources mem=invalid-mem conda --env base echo 'Bad memory'" \
        2

    # Test 18: Non-existent alias (should fail)
    run_test "Non-existent alias" \
        "qxub alias nonexistent-alias echo 'Bad alias'" \
        2

    # Test 19: Empty command (should fail)
    run_test "Empty command" \
        "qxub conda --env base" \
        2
}

# Test complex scenarios
test_complex_scenarios() {
    log_info "=== Testing Complex Scenarios ==="

    # Test 20: Long command with pipes and redirects
    run_test "Complex command with pipes" \
        "qxub conda --env base 'echo \"test data\" | wc -l > output.txt && cat output.txt'" \
        0

    # Test 21: Command with quotes and special characters
    run_test "Command with special characters" \
        "qxub conda --env base 'python -c \"print(\\\"Hello, World!\\\")\"'" \
        0

    # Test 22: Multi-line command
    run_test "Multi-line command" \
        "qxub conda --env base 'echo \"Line 1\" && echo \"Line 2\" && echo \"Line 3\"'" \
        0

    # Test 23: Command with environment variables
    run_test "Command with env vars" \
        "qxub conda --env base 'export TEST_VAR=hello && echo \$TEST_VAR'" \
        0

    # Test 24: Resource-intensive simulation
    run_test "Resource intensive job" \
        "qxub --resources walltime=00:05:00 --resources mem=1GB --resources ncpus=2 conda --env base 'python -c \"import time; time.sleep(1); print(\\\"Resource test\\\")\"'" \
        0
}

# Test config file variations
test_config_variations() {
    log_info "=== Testing Config Variations ==="

    # Test 25: Config show command
    run_test "Show config" \
        "qxub config get" \
        0

    # Test 26: Show specific config value
    run_test "Show specific config" \
        "qxub config get defaults.name" \
        0

    # Test 27: List aliases
    run_test "List aliases" \
        "qxub config alias list-aliases" \
        0

    # Test 28: Test alias dry run
    run_test "Alias test (dry run)" \
        "qxub config alias test test-small" \
        0
}

# Print final results
print_results() {
    echo
    echo "========================================"
    echo "           INTEGRATION TEST RESULTS"
    echo "========================================"
    echo "Total Tests: $TOTAL_TESTS"
    echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
    echo -e "Failed: ${RED}$FAILED_TESTS${NC}"

    if [ $FAILED_TESTS -gt 0 ]; then
        echo
        echo "Failed Test Details:"
        printf '%s\n' "${FAILED_TEST_DETAILS[@]}"
        echo
        exit 1
    else
        echo -e "\n${GREEN}🎉 All tests passed!${NC}"
        exit 0
    fi
}

# Main execution
main() {
    log_info "Starting qxub conda integration tests..."

    # Setup
    setup_test_env

    # Run test suites
    test_basic_conda
    test_config_integration
    test_alias_functionality
    test_edge_cases
    test_complex_scenarios
    test_config_variations

    # Cleanup and results
    cleanup_test_env
    print_results
}

# Execute main function
main "$@"
