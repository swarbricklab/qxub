#!/bin/bash
#
# Updated integration test for qxub unified CLI using dry-run mode.
#
# This test validates the new 2.0 unified interface command generation and parsing
# without actually submitting jobs to the PBS queue. Perfect for rapid
# testing during development.
#

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

log_test() {
    echo -e "${BLUE}[TEST]${NC} $1"
}

run_dry_test() {
    local test_name="$1"
    local test_cmd="$2"
    local expected_exit_code="${3:-0}"

    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    log_test "$test_name"

    # Add --dry flag to ensure we don't actually submit, but only for execution commands
    # Don't add --dry to management commands (config, alias, history, resources)
    if [[ "$test_cmd" != *"--dry"* && "$test_cmd" != *"config"* && "$test_cmd" != *"alias"* && "$test_cmd" != *"history"* && "$test_cmd" != *"resources"* ]]; then
        test_cmd="$test_cmd --dry"
    fi

    echo "Command: $test_cmd"

    # Capture both stdout and stderr, and exit code
    local output
    local exit_code

    if output=$(eval "$test_cmd" 2>&1); then
        exit_code=0
    else
        exit_code=$?
    fi

    # Check exit code
    if [ "$exit_code" -eq "$expected_exit_code" ]; then
        log_success "$test_name (exit code: $exit_code)"
        PASSED_TESTS=$((PASSED_TESTS + 1))

        # For successful tests, show key parts of output
        if [ "$exit_code" -eq 0 ]; then
            echo "  Output preview: $(echo "$output" | head -n 1)"
        fi
    else
        log_error "$test_name (expected exit code: $expected_exit_code, got: $exit_code)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        FAILED_TEST_DETAILS+=("$test_name: Expected exit $expected_exit_code, got $exit_code")
        echo "  Output: $output"
    fi

    echo ""
}

# Test execution contexts
test_execution_contexts() {
    log_info "Testing execution contexts..."

    # Conda environment tests
    run_dry_test \
        "Basic conda execution with --env" \
        "qxub --env base -- echo 'Hello conda'"

    run_dry_test \
        "Basic conda execution with --conda" \
        "qxub --conda base -- echo 'Hello conda'"

    # Module tests - single module with --mod
    run_dry_test \
        "Single module with --mod" \
        "qxub --mod python3 -- echo 'Hello module'"

    # Module tests - multiple modules with --mod
    run_dry_test \
        "Multiple modules with --mod" \
        "qxub --mod python3 --mod gcc -- echo 'Hello modules'"

    # Module tests - comma-separated with --mods
    run_dry_test \
        "Comma-separated modules with --mods" \
        "qxub --mods python3,gcc -- echo 'Hello modules'"

    # Module tests - comma-separated with --modules
    run_dry_test \
        "Comma-separated modules with --modules" \
        "qxub --modules python3,gcc -- echo 'Hello modules'"

    # Singularity tests
    run_dry_test \
        "Singularity with --sif" \
        "qxub --sif container.sif -- echo 'Hello singularity'"

    run_dry_test \
        "Singularity with --sing" \
        "qxub --sing container.sif -- echo 'Hello singularity'"

    run_dry_test \
        "Singularity with --singularity" \
        "qxub --singularity container.sif -- echo 'Hello singularity'"

    # Default execution tests (no environment context)
    run_dry_test \
        "Default execution basic command" \
        "qxub -- echo 'Hello default'"

    run_dry_test \
        "Default execution with pre command" \
        "qxub --pre 'echo Starting' -- echo 'Hello default'"

    run_dry_test \
        "Default execution with post command" \
        "qxub --post 'echo Done' -- echo 'Hello default'"

    run_dry_test \
        "Default execution with pre and post" \
        "qxub --pre 'echo Starting' --post 'echo Done' -- echo 'Hello default'"

    run_dry_test \
        "Default execution with PBS resources" \
        "qxub -l walltime=01:00:00 -l mem=8GB -- echo 'Hello default'"
}

# Test mutual exclusivity
test_mutual_exclusivity() {
    log_info "Testing mutual exclusivity validation..."

    run_dry_test \
        "Conda + Module (should fail)" \
        "qxub --env base --mod python3 -- echo 'Should fail'" \
        1

    run_dry_test \
        "Conda + Singularity (should fail)" \
        "qxub --env base --sif container.sif -- echo 'Should fail'" \
        1

    run_dry_test \
        "Module + Singularity (should fail)" \
        "qxub --mod python3 --sif container.sif -- echo 'Should fail'" \
        1

    run_dry_test \
        "All three contexts (should fail)" \
        "qxub --env base --mod python3 --sif container.sif -- echo 'Should fail'" \
        1
}

# Test PBS options integration
test_pbs_options() {
    log_info "Testing PBS options with unified interface..."

    run_dry_test \
        "Conda with queue option" \
        "qxub --env base --queue normal -- echo 'Queue test'"

    run_dry_test \
        "Module with resources" \
        "qxub --mod python3 --resources walltime=01:00:00 --resources mem=8GB -- echo 'Resources test'"

    run_dry_test \
        "Singularity with job name" \
        "qxub --sif container.sif --name test-job -- echo 'Name test'"

    run_dry_test \
        "Conda with project" \
        "qxub --env base --project a56 -- echo 'Project test'"
}

# Test error conditions
test_error_conditions() {
    log_info "Testing error conditions..."

    run_dry_test \
        "No execution context (uses default execution)" \
        "qxub -- echo 'No context'"

    run_dry_test \
        "Execution context without command (should fail)" \
        "qxub --env base" \
        2

    run_dry_test \
        "Empty environment name (uses default execution)" \
        "qxub --env '' -- echo 'Empty env'"
}

# Test complex commands
test_complex_commands() {
    log_info "Testing complex command scenarios..."

    run_dry_test \
        "Command with options after --" \
        "qxub --env base -- python -c 'print(\"hello\")'"

    run_dry_test \
        "Multi-word command" \
        "qxub --mod python3 -- bash -c 'echo \"Complex command\"'"

    run_dry_test \
        "Command with quotes and escaping" \
        "qxub --sif container.sif -- echo 'Quote test: \"hello world\"'"

    run_dry_test \
        "Pipeline command" \
        "qxub --env base -- bash -c 'echo \"test\" | grep test'"
}

# Test alternative option combinations
test_alternative_options() {
    log_info "Testing alternative option name combinations..."

    # Test that alternatives produce same results
    run_dry_test \
        "Alternative conda options equivalent" \
        "qxub --conda myenv -- echo 'Alternative test'"

    run_dry_test \
        "Alternative module options equivalent" \
        "qxub --modules python3,gcc -- echo 'Alternative test'"

    run_dry_test \
        "Alternative singularity options equivalent" \
        "qxub --sing container.sif -- echo 'Alternative test'"
}

# Test management commands still work
test_management_commands() {
    log_info "Testing management commands..."

    run_dry_test \
        "Config help" \
        "qxub config --help"

    run_dry_test \
        "Alias help" \
        "qxub alias --help"

    run_dry_test \
        "History help" \
        "qxub history --help"

    run_dry_test \
        "Resources help" \
        "qxub resources --help"
}

# Test comprehensive scenarios
test_comprehensive_scenarios() {
    log_info "Testing comprehensive real-world scenarios..."

    run_dry_test \
        "Data science workflow" \
        "qxub --conda pytorch --queue gpuvolta --resources ngpus=1 --resources ncpus=12 --name gpu-training -- python train.py --epochs 100"

    run_dry_test \
        "Bioinformatics pipeline" \
        "qxub --modules samtools,bwa --resources walltime=02:00:00 --resources mem=16GB -- bash pipeline.sh"

    run_dry_test \
        "Container workflow" \
        "qxub --singularity /containers/blast.sif --bind /data --name blast-search -- blastn -query input.fa -db nt"
}

# Setup test environment
setup_test_env() {
    log_info "Setting up test environment..."

    # Backup existing config if present
    if [ -f ~/.config/qxub/config.yaml ]; then
        cp ~/.config/qxub/config.yaml ~/.config/qxub/config.yaml.drytest.backup
        log_info "Backed up existing config file"
    fi

    # Ensure config directory exists
    mkdir -p ~/.config/qxub

    # Create minimal test config
    cat > ~/.config/qxub/config.yaml << 'EOF'
defaults:
  name: "test-job"
  time: "00:30:00"
  mem: "4GB"
  ncpus: 1
  project: "a56"
  queue: "normal"
  execdir: "/scratch/{project}/{user}/work"
EOF

    log_info "Created test config file"
}

# Cleanup test environment
cleanup_test_env() {
    log_info "Cleaning up test environment..."

    # Restore original config if it existed
    if [ -f ~/.config/qxub/config.yaml.drytest.backup ]; then
        mv ~/.config/qxub/config.yaml.drytest.backup ~/.config/qxub/config.yaml
        log_info "Restored original config file"
    else
        rm -f ~/.config/qxub/config.yaml
        log_info "Removed test config file"
    fi
}

# Print summary
print_summary() {
    echo ""
    echo "====================================="
    echo "         TEST SUMMARY"
    echo "====================================="
    echo "Total tests: $TOTAL_TESTS"
    echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
    echo -e "Failed: ${RED}$FAILED_TESTS${NC}"

    if [ $FAILED_TESTS -gt 0 ]; then
        echo ""
        echo "FAILED TEST DETAILS:"
        for detail in "${FAILED_TEST_DETAILS[@]}"; do
            echo -e "${RED}  - $detail${NC}"
        done
        echo ""
        echo -e "${RED}❌ Some tests failed${NC}"
        exit 1
    else
        echo ""
        echo -e "${GREEN}✅ All tests passed!${NC}"
        exit 0
    fi
}

# Main execution
main() {
    echo "======================================="
    echo "    qxub 2.0 Unified CLI Test Suite"
    echo "======================================="
    echo ""

    setup_test_env

    # Trap to ensure cleanup happens
    trap cleanup_test_env EXIT

    # Run all test suites
    test_execution_contexts
    test_mutual_exclusivity
    test_pbs_options
    test_error_conditions
    test_complex_commands
    test_alternative_options
    test_management_commands
    test_comprehensive_scenarios

    print_summary
}

# Run main function
main "$@"
