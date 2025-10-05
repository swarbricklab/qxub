#!/bin/bash
#
# Quick integration test for qxub conda subcommand using dry-run mode.
#
# This test focuses on validating the command generation and parsing
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
    
    # Insert --dry-run after 'qxub' but before the subcommand
    local dry_cmd=$(echo "$test_cmd" | sed 's/qxub /qxub --dry-run /')
    echo "Command: $dry_cmd"
    
    # Run the command and capture output and exit code
    if output=$(eval "$dry_cmd" 2>&1); then
        actual_exit_code=0
    else
        actual_exit_code=$?
    fi
    
    # Check if exit code matches expected
    if [ "$actual_exit_code" -eq "$expected_exit_code" ]; then
        log_success "$test_name"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        # Show first few lines of output for verification
        echo "$output" | head -n 5
        if [ $(echo "$output" | wc -l) -gt 5 ]; then
            echo "... (output truncated)"
        fi
    else
        log_error "$test_name (expected exit code $expected_exit_code, got $actual_exit_code)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        FAILED_TEST_DETAILS+=("$test_name: Expected exit $expected_exit_code, got $actual_exit_code")
        echo "Full Output: $output"
    fi
    echo "---"
}

# Setup test environment (minimal for dry runs)
setup_test_env() {
    log_info "Setting up test environment for dry runs..."
    
    # Create test config directory if it doesn't exist
    mkdir -p ~/.config/qxub
    
    # Backup existing config if it exists
    if [ -f ~/.config/qxub/config.yaml ]; then
        cp ~/.config/qxub/config.yaml ~/.config/qxub/config.yaml.drytest.backup
        log_info "Backed up existing config"
    fi
    
    # Create a comprehensive test config file
    cat > ~/.config/qxub/config.yaml << 'EOF'
defaults:
  name: "test-job-{timestamp}"
  time: "00:30:00"
  mem: "4GB"
  ncpus: 1
  project: "a56"
  queue: "normal"
  execdir: "/scratch/{project}/{user}/work"
  
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
    
  gpu:
    subcommand: conda
    env: tensorflow
    name: "gpu-job"
    queue: "gpuvolta"
    resources:
      - "ngpus=1"
      - "gputype=V100"
      
  parallel:
    subcommand: conda
    env: mpi
    name: "parallel-job"
    ncpus: 8
    resources:
      - "mem=16GB"
      - "jobfs=10GB"
EOF
    
    log_info "Created comprehensive test config file"
}

# Cleanup test environment
cleanup_test_env() {
    log_info "Cleaning up test environment..."
    
    # Restore original config if it existed
    if [ -f ~/.config/qxub/config.yaml.drytest.backup ]; then
        mv ~/.config/qxub/config.yaml.drytest.backup ~/.config/qxub/config.yaml
        log_info "Restored original config"
    else
        rm -f ~/.config/qxub/config.yaml
        rmdir ~/.config/qxub 2>/dev/null || true
    fi
}

# Test basic conda functionality with all CLI options
test_all_cli_options() {
    log_info "=== Testing All CLI Options ==="
    
    # Test 1: All main PBS options
    run_dry_test "All main PBS options" \
        "qxub --name cli-test --project a56 --queue normal conda --env base echo 'All options test'"
    
    # Test 2: Resources option (multiple formats)
    run_dry_test "Resources - mem and ncpus" \
        "qxub --resources mem=8GB --resources ncpus=2 conda --env base echo 'Resources test'"
    
    # Test 3: Resources - time format variations
    run_dry_test "Resources - walltime variations" \
        "qxub --resources walltime=01:30:00 conda --env base echo 'Walltime test'"
    
    # Test 4: Output and error files
    run_dry_test "Custom output/error files" \
        "qxub --out /tmp/test.out --err /tmp/test.err conda --env base echo 'Output test'"
    
    # Test 5: Execution directory
    run_dry_test "Custom execution directory" \
        "qxub --execdir /scratch/a56/test conda --env base echo 'Execdir test'"
    
    # Test 6: Job log
    run_dry_test "Job log option" \
        "qxub --joblog /tmp/job.log conda --env base echo 'Joblog test'"
    
    # Test 7: Verbose flag
    run_dry_test "Verbose mode" \
        "qxub --verbose conda --env base echo 'Verbose test'"
    
    # Test 8: Multiple verbose levels
    run_dry_test "Multiple verbose levels" \
        "qxub -vvv conda --env base echo 'Very verbose test'"
}

# Test conda-specific options
test_conda_specific_options() {
    log_info "=== Testing Conda-Specific Options ==="
    
    # Test 9: Custom template
    run_dry_test "Custom template" \
        "qxub conda --env base --template /tmp/custom.pbs echo 'Template test'"
    
    # Test 10: Pre-command
    run_dry_test "Pre-command" \
        "qxub conda --env base --pre 'export CUSTOM_VAR=test' echo 'Pre test'"
    
    # Test 11: Post-command  
    run_dry_test "Post-command" \
        "qxub conda --env base --post 'echo \"Job completed\"' echo 'Post test'"
    
    # Test 12: Pre and Post together
    run_dry_test "Pre and Post commands" \
        "qxub conda --env base --pre 'echo \"Starting\"' --post 'echo \"Ending\"' echo 'Pre-post test'"
    
    # Test 13: Complex command with conda options
    run_dry_test "Complex conda command" \
        "qxub conda --env myenv --pre 'module load python' 'python -c \"import sys; print(sys.version)\"'"
}

# Test config integration with dry runs
test_config_with_dry_runs() {
    log_info "=== Testing Config Integration ==="
    
    # Test 14: Config defaults only
    run_dry_test "Config defaults only" \
        "qxub conda --env base echo 'Config test'"
    
    # Test 15: Override single config value
    run_dry_test "Override config name" \
        "qxub --name override-test conda --env base echo 'Override test'"
    
    # Test 16: Override multiple config values
    run_dry_test "Override multiple config values" \
        "qxub --name multi-override --project b01 --queue express conda --env base echo 'Multi-override test'"
    
    # Test 17: Template variables in config
    run_dry_test "Template variables" \
        "qxub conda --env base echo 'Template vars test'"
}

# Test alias functionality with dry runs
test_alias_dry_runs() {
    log_info "=== Testing Alias Functionality ==="
    
    # Test 18: Quick alias
    run_dry_test "Quick alias" \
        "qxub alias quick echo 'Quick alias test'"
    
    # Test 19: Big memory alias
    run_dry_test "Big memory alias" \
        "qxub alias bigmem 'python -c \"print(\\\"Big memory test\\\")\"'"
    
    # Test 20: GPU alias
    run_dry_test "GPU alias" \
        "qxub alias gpu 'python -c \"import torch; print(torch.cuda.is_available())\"'"
    
    # Test 21: Parallel alias
    run_dry_test "Parallel alias" \
        "qxub alias parallel mpirun python parallel_script.py"
    
    # Test 22: Alias with CLI overrides
    run_dry_test "Alias with overrides" \
        "qxub --name override-alias --resources walltime=00:45:00 alias quick echo 'Alias override test'"
}

# Test edge cases and complex scenarios
test_edge_cases_dry() {
    log_info "=== Testing Edge Cases ==="
    
    # Test 23: Very long command
    run_dry_test "Long command" \
        "qxub conda --env base 'for i in {1..10}; do echo \"This is iteration \$i\"; done'"
    
    # Test 24: Command with quotes and escapes
    run_dry_test "Complex quoting" \
        "qxub conda --env base 'python -c \"print(\\\"Hello, \\\\\"World\\\\\"!\\\")\"'"
    
    # Test 25: Command with pipes and redirects
    run_dry_test "Pipes and redirects" \
        "qxub conda --env base 'echo \"test data\" | grep \"test\" > output.txt'"
    
    # Test 26: Multiple resources
    run_dry_test "Multiple resources" \
        "qxub --resources mem=16GB --resources ncpus=4 --resources jobfs=50GB conda --env base echo 'Multi-resource test'"
    
    # Test 27: Special characters in job name
    run_dry_test "Special chars in name" \
        "qxub --name 'test-job_2024.01.01' conda --env base echo 'Special name test'"
}

# Test error conditions (should fail)
test_error_conditions() {
    log_info "=== Testing Error Conditions ==="
    
    # Test 28: Missing environment (may succeed in dry-run)
    run_dry_test "Missing environment" \
        "qxub conda echo 'No env'" \
        0
    
    # Test 29: Non-existent alias (should fail)  
    run_dry_test "Non-existent alias" \
        "qxub alias nonexistent echo 'Bad alias'" \
        1
    
    # Test 30: Empty command (may succeed in dry-run)
    run_dry_test "Empty command" \
        "qxub conda --env base" \
        0
    
    # Test 31: Invalid resource format (may succeed in dry-run)
    run_dry_test "Invalid resource format" \
        "qxub --resources invalid_format conda --env base echo 'Bad resource'" \
        0
}

# Print final results
print_results() {
    echo
    echo "========================================"
    echo "         DRY RUN TEST RESULTS"
    echo "========================================"
    echo "Total Tests: $TOTAL_TESTS"
    echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
    echo -e "Failed: ${RED}$FAILED_TESTS${NC}"
    
    if [ $FAILED_TESTS -gt 0 ]; then
        echo
        echo "Failed Test Details:"
        printf '%s\n' "${FAILED_TEST_DETAILS[@]}"
        echo
        echo -e "${RED}❌ Some tests failed!${NC}"
        exit 1
    else
        echo -e "\n${GREEN}🎉 All dry-run tests passed!${NC}"
        echo -e "${BLUE}💡 Tip: Run the full integration test with './test_conda_integration.sh' to test actual job submission${NC}"
        exit 0
    fi
}

# Main execution
main() {
    log_info "Starting qxub conda dry-run integration tests..."
    echo -e "${BLUE}ℹ️  These tests use --dry-run mode and won't submit actual jobs${NC}"
    echo
    
    # Setup
    setup_test_env
    
    # Run test suites
    test_all_cli_options
    test_conda_specific_options
    test_config_with_dry_runs
    test_alias_dry_runs
    test_edge_cases_dry
    test_error_conditions
    
    # Cleanup and results
    cleanup_test_env
    print_results
}

# Execute main function
main "$@"