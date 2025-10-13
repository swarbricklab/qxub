#!/bin/bash

# Test script for GNU find commands with qxub
# Tests various find patterns that are known to cause shell parsing issues

set -e

echo "🔍 Testing GNU find command patterns with qxub"
echo

# Activate qxub environment
source venv/bin/activate

# Create test directory structure for find commands
TEST_DIR="/tmp/qxub_find_test_$$"
mkdir -p "$TEST_DIR"/{subdir1,subdir2}
touch "$TEST_DIR"/{file1.txt,file2.py,file3.log}
touch "$TEST_DIR"/subdir1/{nested1.txt,nested2.py}
touch "$TEST_DIR"/subdir2/{nested3.txt,nested4.log}

echo "📁 Created test directory: $TEST_DIR"
echo

# Cleanup function
cleanup() {
    echo "🧹 Cleaning up test directory..."
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

# Test counter
test_count=0
pass_count=0
fail_count=0

# Helper function to run test
run_test() {
    local test_name="$1"
    local traditional_cmd="$2"
    local cmd_syntax="$3"

    ((test_count++))
    echo "Test $test_count: $test_name"

    echo "  Traditional syntax:"
    if traditional_result=$(eval "qxub --env base --dry -- $traditional_cmd" 2>&1); then
        echo "    ✅ SUCCESS: Traditional syntax worked"
        echo "    Output: $(echo "$traditional_result" | grep "Command to execute:" | cut -d: -f2-)"
        ((pass_count++))
    else
        echo "    ❌ FAILED: Traditional syntax failed"
        echo "    Error: $traditional_result"
        ((fail_count++))
    fi

    echo "  --cmd syntax:"
    if cmd_result=$(qxub --env base --dry --cmd "$cmd_syntax" 2>&1); then
        echo "    ✅ SUCCESS: --cmd syntax worked"
        echo "    Output: $(echo "$cmd_result" | grep "Command to execute:" | cut -d: -f2-)"
        ((pass_count++))
    else
        echo "    ❌ FAILED: --cmd syntax failed"
        echo "    Error: $cmd_result"
        ((fail_count++))
    fi
    echo
}

# Test 1: Basic find with -exec and \;
run_test "Basic find with -exec and semicolon terminator" \
    "find $TEST_DIR -name '*.txt' -exec echo 'Found: {}' \;" \
    "find $TEST_DIR -name '*.txt' -exec echo 'Found: {}' \;"

# Test 2: Find with -exec and +
run_test "Find with -exec and + terminator" \
    "find $TEST_DIR -name '*.py' -exec echo 'Python file: {}' +" \
    "find $TEST_DIR -name '*.py' -exec echo 'Python file: {}' +"

# Test 3: Find with complex -exec command containing quotes
run_test "Find with complex -exec containing nested quotes" \
    "find $TEST_DIR -type f -exec sh -c 'echo \"Processing: \$1\"' _ {} \;" \
    "find $TEST_DIR -type f -exec sh -c 'echo \"Processing: \$1\"' _ {} \;"

# Test 4: Find with variables in path
run_test "Find with submission-time variables in path" \
    "find \${HOME} -maxdepth 1 -name '.*' -type f" \
    "find \${HOME} -maxdepth 1 -name '.*' -type f"

# Test 5: Find with execution-time variables
run_test "Find with execution-time variables" \
    "find $TEST_DIR -name '*.txt' -exec echo 'Job \${{PBS_JOBID}} found: {}' \;" \
    "find $TEST_DIR -name '*.txt' -exec echo 'Job \${{PBS_JOBID}} found: {}' \;"

# Test 6: Find with -print0 and xargs
run_test "Find with -print0 piped to xargs" \
    "find $TEST_DIR -name '*.txt' -print0 | xargs -0 ls -la" \
    "find $TEST_DIR -name '*.txt' -print0 | xargs -0 ls -la"

# Test 7: Find with regex patterns
run_test "Find with regex patterns" \
    "find $TEST_DIR -regex '.*\.\(txt\|py\)'" \
    "find $TEST_DIR -regex '.*\.\(txt\|py\)'"

# Test 8: Find with -exec and shell metacharacters
run_test "Find with shell metacharacters in -exec" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \$0: \$1 | wc -c' find-script {} \;" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \$0: \$1 | wc -c' find-script {} \;"

# Test 9: Find with -exec and AWK-style variables
run_test "Find with AWK-style field references" \
    "find $TEST_DIR -name '*.txt' -exec awk '{print \"Line \" NR \": \" \$0}' {} \;" \
    "find $TEST_DIR -name '*.txt' -exec awk '{print \"Line \" NR \": \" \$0}' {} \;"

# Test 10: Find with -exec containing dollars and braces
run_test "Find with complex dollar and brace patterns" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"File: \$1, Cost: \\\$100, Fields: \$# args\"' _ {} \;" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"File: \$1, Cost: \\\$100, Fields: \$# args\"' _ {} \;"

# Test 11: Find with glob patterns in -name
run_test "Find with glob patterns in -name" \
    "find $TEST_DIR -name '[fn]*' -type f" \
    "find $TEST_DIR -name '[fn]*' -type f"

# Test 12: Find with -exec containing backticks
run_test "Find with command substitution in -exec" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"Found at: \$(date): \$1\"' _ {} \;" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"Found at: \$(date): \$1\"' _ {} \;"

# Test 13: Find with multiple -exec actions
run_test "Find with multiple -exec actions" \
    "find $TEST_DIR -name '*.txt' -exec echo 'Processing: {}' \; -exec wc -l {} \;" \
    "find $TEST_DIR -name '*.txt' -exec echo 'Processing: {}' \; -exec wc -l {} \;"

# Test 14: Find with -exec and environment variables
run_test "Find mixing submission and execution time variables" \
    "find \${HOME} -maxdepth 1 -name '.*' -exec echo 'User \${USER} found \${{PBS_JOBID}}: {}' \;" \
    "find \${HOME} -maxdepth 1 -name '.*' -exec echo 'User \${USER} found \${{PBS_JOBID}}: {}' \;"

# Test 15: Find with complex shell escaping
run_test "Find with complex shell escaping patterns" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"File: \\\$1\" | sed \"s/.*\\///\"' _ {} \;" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"File: \\\$1\" | sed \"s/.*\\///\"' _ {} \;"

# Test 16: Find with -execdir
run_test "Find with -execdir action" \
    "find $TEST_DIR -name '*.txt' -execdir echo 'In directory: {}' \;" \
    "find $TEST_DIR -name '*.txt' -execdir echo 'In directory: {}' \;"

# Test 17: Find with -delete action
run_test "Find with -delete action (dry run safe)" \
    "find $TEST_DIR -name 'nonexistent*' -delete" \
    "find $TEST_DIR -name 'nonexistent*' -delete"

# Test 18: Find with complex logical expressions
run_test "Find with complex logical expressions" \
    "find $TEST_DIR \\( -name '*.txt' -o -name '*.py' \\) -a -type f" \
    "find $TEST_DIR \\( -name '*.txt' -o -name '*.py' \\) -a -type f"

# Test 19: Find with printf action
run_test "Find with printf action" \
    "find $TEST_DIR -name '*.txt' -printf 'Found: %p, size: %s bytes\\n'" \
    "find $TEST_DIR -name '*.txt' -printf 'Found: %p, size: %s bytes\\n'"

# Test 20: Find with very complex -exec containing JSON-like structure
run_test "Find with JSON-like structure in -exec" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"{\\\"file\\\": \\\"\$1\\\", \\\"user\\\": \\\"\${USER}\\\", \\\"job\\\": \\\"\${{PBS_JOBID}}\\\"}\"' _ {} \;" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'echo \"{\\\"file\\\": \\\"\$1\\\", \\\"user\\\": \\\"\${USER}\\\", \\\"job\\\": \\\"\${{PBS_JOBID}}\\\"}\"' _ {} \;"

# Test 21: Find with -exec and pipes
run_test "Find with -exec containing pipes" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'cat \$1 | head -5 | wc -l' _ {} \;" \
    "find $TEST_DIR -name '*.txt' -exec sh -c 'cat \$1 | head -5 | wc -l' _ {} \;"

# Test 22: Find with special characters in filenames
run_test "Find with special characters and spaces" \
    "find $TEST_DIR -name '*' -exec echo 'File with spaces: \"{}\"' \;" \
    "find $TEST_DIR -name '*' -exec echo 'File with spaces: \"{}\"' \;"

# Summary
echo "📊 Test Summary:"
echo "  Total tests: $((test_count * 2))  # Each test runs both syntaxes"
echo "  Passed: $pass_count"
echo "  Failed: $fail_count"
echo "  Success rate: $(( pass_count * 100 / (test_count * 2) ))%"
echo

echo "🎯 Key Findings:"
if [ $fail_count -eq 0 ]; then
    echo "  ✅ All find command patterns work with both syntaxes!"
else
    echo "  ⚠️  Some patterns failed - check output above for details"
fi

echo
echo "💡 Find Command Patterns Tested:"
echo "  • Basic -exec with \; terminator"
echo "  • -exec with + terminator"
echo "  • Complex nested quotes in -exec"
echo "  • Submission-time variables (\${var})"
echo "  • Execution-time variables (\${{var}})"
echo "  • Pipes with -print0 and xargs"
echo "  • Regex patterns"
echo "  • Shell metacharacters (\$0, \$1, \$#)"
echo "  • AWK field references"
echo "  • Command substitution with backticks"
echo "  • Multiple -exec actions"
echo "  • Complex escaping patterns"
echo "  • -execdir and -delete actions"
echo "  • Logical expressions with parentheses"
echo "  • -printf with format specifiers"
echo "  • JSON-like structures"
echo "  • Pipes in -exec commands"
echo "  • Special characters in paths"

echo
echo "🔍 Usage Patterns:"
echo "  Traditional: qxub --env base -- find /path -name '*.txt' -exec cmd {} \\;"
echo "  --cmd:       qxub --env base --cmd 'find /path -name \"*.txt\" -exec cmd {} \\;'"
echo "  With vars:   qxub --env base --cmd 'find \${HOME} -exec echo \"User \${USER}: {}\" \\;'"
