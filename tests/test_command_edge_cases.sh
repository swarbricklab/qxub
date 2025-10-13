#!/bin/bash
#
# Command Handling Edge Cases Test Script
#
# This script tests complex command scenarios to understand current behavior
# and identify potential issues with qxub's base64 encoding approach.
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

# Helper functions
log_info() {
    echo -e "${YELLOW}[INFO]${NC} $1"
}

log_test() {
    echo -e "${BLUE}[TEST]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

log_error() {
    echo -e "${RED}[FAIL]${NC} $1"
}

run_edge_case_test() {
    local test_name="$1"
    local test_cmd="$2"
    local expected_exit_code="${3:-0}"
    local expected_pattern="${4:-}"

    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    log_test "$test_name"

    # Add --dry flag to test command
    local dry_cmd=$(echo "$test_cmd" | sed 's/qxub /qxub --dry /')
    echo "Command: $dry_cmd"

    # Run the command and capture output and exit code
    if output=$(eval "$dry_cmd" 2>&1); then
        actual_exit_code=0
    else
        actual_exit_code=$?
    fi

    # Check exit code
    if [ "$actual_exit_code" -eq "$expected_exit_code" ]; then
        # Check pattern if provided
        if [[ -z "$expected_pattern" ]] || echo "$output" | grep -q "$expected_pattern"; then
            log_success "$test_name"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            log_error "$test_name - Expected pattern '$expected_pattern' not found"
            echo "Output: $output"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        log_error "$test_name - Expected exit code $expected_exit_code, got $actual_exit_code"
        echo "Output: $output"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi

    echo "---"
}

print_results() {
    echo "=================================================="
    echo "           EDGE CASE TEST RESULTS"
    echo "=================================================="
    echo "Total tests: $TOTAL_TESTS"
    echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"
    if [ $FAILED_TESTS -gt 0 ]; then
        echo -e "${RED}Failed: $FAILED_TESTS${NC}"
    else
        echo "Failed: $FAILED_TESTS"
    fi
    echo "=================================================="

    if [ $FAILED_TESTS -eq 0 ]; then
        echo -e "${GREEN}🎉 All edge case tests passed!${NC}"
        exit 0
    else
        echo -e "${RED}❌ Some edge case tests failed.${NC}"
        exit 1
    fi
}

# Test Category 1: Advanced Quoting Edge Cases
test_advanced_quoting() {
    log_info "=== Testing Advanced Quoting Edge Cases ==="

    # Nested quotes
    run_edge_case_test "Nested quotes - single inside double" \
        "qxub --env base -- echo \"outer 'single' quotes\""

    run_edge_case_test "Nested quotes - double inside single" \
        "qxub --env base -- echo 'outer \"double\" quotes'"

    # Python with complex quotes
    run_edge_case_test "Python complex nested quotes" \
        "qxub --env base -- python -c \"print('single inside double: \\\"nested double\\\"')\""

    # Pathological escaping
    run_edge_case_test "Multi-level escaping" \
        "qxub --env base -- python -c \"print(\\\"Hello, \\\\\\\"World\\\\\\\"!\\\")\""

    # Mixed quote types in shell commands
    run_edge_case_test "Mixed quotes in bash command" \
        "qxub --env base -- bash -c 'echo \"level1 '\\''level2'\\'' back to level1\"'"
}

# Test Category 2: Variable Expansion Timing
test_variable_expansion() {
    log_info "=== Testing Variable Expansion Timing ==="

    # Should expand at submission time (current shell)
    run_edge_case_test "Submission-time variable expansion" \
        "qxub --env base -- echo \"Submitted by: \$USER\"" \
        0 "$USER"

    # Should NOT expand at submission time (escaped for execution time)
    run_edge_case_test "Execution-time variable expansion (escaped)" \
        "qxub --env base -- echo \"Job user: \\\$USER\""

    # PBS-specific variables (only available at execution)
    run_edge_case_test "PBS job variables" \
        "qxub --env base -- echo \"PBS job ID: \\\$PBS_JOBID\""

    # Mixed expansion scenario
    run_edge_case_test "Mixed expansion timing" \
        "qxub --env base -- echo \"Job \\\$PBS_JOBID submitted by \$USER on \\\$HOSTNAME\""

    # Parameter expansion
    run_edge_case_test "Parameter expansion" \
        "qxub --env base -- bash -c \"VAR=test; echo \\\"\\\${VAR:-default}\\\"\""
}

# Test Category 3: Shell Metacharacters
test_shell_metacharacters() {
    log_info "=== Testing Shell Metacharacters ==="

    # Command substitution
    run_edge_case_test "Command substitution with \$()" \
        "qxub --env base -- bash -c \"echo \\\"Today is \\\$(date)\\\"\""

    run_edge_case_test "Command substitution with backticks" \
        "qxub --env base -- bash -c \"echo \\\"Today is \\\`date\\\`\\\"\""

    # Nested command substitution
    run_edge_case_test "Nested command substitution" \
        "qxub --env base -- bash -c \"echo \\\"\\\$(basename \\\$(dirname \\\$(pwd)))\\\"\""

    # Here documents
    run_edge_case_test "Here document simulation" \
        "qxub --env base -- bash -c \"cat <<EOF > output.txt
Line 1
Line 2 with \\\$VAR
EOF\""

    # Process substitution (where supported)
    run_edge_case_test "Process substitution" \
        "qxub --env base -- bash -c \"echo \\\"Compare: \\\$(diff <(echo a) <(echo b) || true)\\\"\""

    # Globbing patterns
    run_edge_case_test "Glob patterns" \
        "qxub --env base -- bash -c \"echo file*.txt\""

    run_edge_case_test "Brace expansion" \
        "qxub --env base -- bash -c \"echo file{1,2,3}.txt\""
}

# Test Category 4: Unicode and Encoding
test_unicode_encoding() {
    log_info "=== Testing Unicode and Encoding ==="

    # Basic Unicode
    run_edge_case_test "Unicode characters" \
        "qxub --env base -- echo \"Unicode: 你好世界\""

    # Emoji
    run_edge_case_test "Emoji characters" \
        "qxub --env base -- echo \"Emoji: 🚀🔬📊\""

    # Accented characters
    run_edge_case_test "Accented characters" \
        "qxub --env base -- echo \"Accents: café résumé naïve\""

    # Unicode in Python
    run_edge_case_test "Unicode in Python" \
        "qxub --env base -- python -c \"print('Greek: αβγδε')\""

    # Control characters
    run_edge_case_test "Control characters" \
        "qxub --env base -- printf \"Tab:\\\\t Newline:\\\\n\""

    # ANSI escape sequences
    run_edge_case_test "ANSI color codes" \
        "qxub --env base -- echo -e \"\\\\x1b[31mRed text\\\\x1b[0m\""
}

# Test Category 5: Large and Complex Commands
test_large_commands() {
    log_info "=== Testing Large and Complex Commands ==="

    # Very long command
    LONG_CMD="python -c \"data='x'*1000; print(f'Generated {len(data)} characters')\""
    run_edge_case_test "Large command payload" \
        "qxub --env base -- $LONG_CMD"

    # Complex pipeline
    run_edge_case_test "Complex pipeline" \
        "qxub --env base -- bash -c \"echo 'test data' | tee /tmp/qxub_test_out | wc -l | tee -a /tmp/qxub_test_log\""

    # Multiple commands with error handling
    run_edge_case_test "Command chain with error handling" \
        "qxub --env base -- bash -c \"echo 'step1' && echo 'step2' && echo 'step3' || echo 'error'\""

    # Loop with variables
    run_edge_case_test "Loop with variable expansion" \
        "qxub --env base -- bash -c \"for i in {1..3}; do echo \\\"Iteration \\\$i: \\\$(date)\\\"; done\""
}

# Test Category 6: Cross-Mode Consistency
test_cross_mode_consistency() {
    log_info "=== Testing Cross-Mode Consistency ==="

    local test_cmd="echo 'Mode test: \\\"Hello World\\\"'"

    # Same command across all execution modes
    run_edge_case_test "Default mode quoting" \
        "qxub --default -- $test_cmd"

    run_edge_case_test "Conda mode quoting" \
        "qxub --env base -- $test_cmd"

    run_edge_case_test "Module mode quoting" \
        "qxub --mod python3 -- $test_cmd"

    run_edge_case_test "Singularity mode quoting" \
        "qxub --sif container.sif -- $test_cmd"
}

# Test Category 7: Pre/Post Command Interactions
test_pre_post_interactions() {
    log_info "=== Testing Pre/Post Command Interactions ==="

    # Variable sharing between pre/main/post
    run_edge_case_test "Pre/post variable sharing" \
        "qxub --env base --pre 'export MYVAR=pre_value' --post 'echo Final: \\\$MYVAR' -- echo 'Main: \\\$MYVAR'"

    # Complex pre command with quotes
    run_edge_case_test "Complex pre command" \
        "qxub --env base --pre 'echo \\\"Pre: \\\$(date)\\\"' -- echo 'Main command'"

    # Post command with pipeline
    run_edge_case_test "Complex post command" \
        "qxub --env base --post 'echo \\\"Done\\\" | tee -a /tmp/qxub_completion_log' -- echo 'Main work'"
}

# Main execution
main() {
    echo "=================================================="
    echo "    QXUB COMMAND HANDLING EDGE CASES TEST"
    echo "=================================================="
    echo "Testing complex command scenarios across all execution modes"
    echo "Using --dry flag to test command generation without submission"
    echo "=================================================="

    test_advanced_quoting
    test_variable_expansion
    test_shell_metacharacters
    test_unicode_encoding
    test_large_commands
    test_cross_mode_consistency
    test_pre_post_interactions

    print_results
}

# Execute main function
main "$@"
