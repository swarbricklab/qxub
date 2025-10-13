#!/bin/bash

# Test script for smart quote processing in --cmd option
# Validates both backward compatibility and new smart quote features

set -e

echo "🧪 Testing Smart Quote Processing in --cmd Option"
echo

# Activate qxub environment
source venv/bin/activate

# Test counter
test_count=0
pass_count=0

# Helper function to run and validate test
run_test() {
    local test_name="$1"
    local cmd_input="$2"
    local expected_pattern="$3"

    ((test_count++))
    echo "Test $test_count: $test_name"
    echo "  Input: $cmd_input"

    if result=$(qxub --env base --dry --cmd "$cmd_input" 2>&1); then
        output=$(echo "$result" | grep "Command to execute:" | cut -d: -f2- | sed 's/^ *//')
        echo "  Output: $output"

        if echo "$output" | grep -q "$expected_pattern"; then
            echo "  ✅ PASS: Expected pattern found"
            ((pass_count++))
        else
            echo "  ❌ FAIL: Expected pattern '$expected_pattern' not found"
            echo "  Full result: $result"
        fi
    else
        echo "  ❌ FAIL: Command execution failed"
        echo "  Error: $result"
    fi
    echo
}

# Test 1: Backward compatibility - single quotes (no smart processing)
run_test "Backward compatibility with single quotes" \
    'echo "Hello ${USER}"' \
    'echo "Hello jr9959"'

# Test 2: Smart quotes - basic conversion
run_test "Smart quotes basic conversion" \
    '"echo \"Hello ${USER}\""' \
    'echo "Hello jr9959"'

# Test 3: Smart quotes with literal dollars
run_test "Smart quotes with literal dollars" \
    '"echo \"Cost: \$100\""' \
    'echo "Cost: $100"'

# Test 4: Smart quotes with mixed variables
run_test "Mixed submission and execution time variables" \
    '"echo \"User: ${USER}, Job: ${{PBS_JOBID}}\""' \
    'echo "User: jr9959, Job: ${PBS_JOBID}"'

# Test 5: Complex find command (the original problem)
run_test "Complex find command with smart quotes" \
    '"find /tmp -name \"*.txt\" -exec echo \"Found: {}\" \;"' \
    'find /tmp -name "*.txt" -exec echo "Found: {}" \\;'

# Test 6: AWK command with field references
run_test "AWK with field references preserved" \
    '"awk \"{print \\\$1, \\\$2}\" file.txt"' \
    'awk "{print $1, $2}" file.txt'

# Test 7: Shell command with nested quotes
run_test "Shell command with nested quotes" \
    '"sh -c \"echo \\\"User: ${USER}\\\" | wc -c\""' \
    'sh -c "echo \\"User: jr9959\\" | wc -c"'

# Test 8: JSON-like output structure
run_test "JSON-like structure output" \
    '"echo \"{\\\"user\\\": \\\"${USER}\\\", \\\"cost\\\": \\\"\$50\\\"}\"" ' \
    'echo "{\\"user\\": \\"jr9959\\", \\"cost\\": \\"$50\\"}"'

# Test 9: Complex find with sh -c and multiple escaping levels
run_test "Complex find with sh -c command" \
    '"find /tmp -type f -exec sh -c \"echo \\\"Processing: \\\$1 for ${USER}\\\"\" _ {} \;"' \
    'find /tmp -type f -exec sh -c "echo \\"Processing: $1 for jr9959\\"" _ {} \\;'

# Test 10: Pipes and command substitution
run_test "Command with pipes preserved" \
    '"find /tmp -name \"*.log\" | head -5 | wc -l"' \
    'find /tmp -name "*.log" | head -5 | wc -l'

echo "📊 Test Results Summary:"
echo "  Total tests: $test_count"
echo "  Passed: $pass_count"
echo "  Failed: $((test_count - pass_count))"
echo "  Success rate: $(( pass_count * 100 / test_count ))%"
echo

if [ $pass_count -eq $test_count ]; then
    echo "🎉 All tests passed! Smart quote processing is working correctly."
else
    echo "⚠️  Some tests failed. Review output above for details."
fi

echo
echo "💡 Smart Quote Processing Rules:"
echo "  • Use double quotes around --cmd for smart processing"
echo "  • \\\" becomes literal \" inside the command"
echo "  • \\\$ becomes literal \$ (no variable expansion)"
echo "  • \${var} still expands at submission time"
echo "  • \${{var}} still converts for execution time"
echo "  • Single quotes preserved literally inside"
echo "  • Backward compatible with existing single-quote syntax"
echo
echo "🔧 Usage Examples:"
echo "  Before: qxub --cmd 'complex command with ugly nested quotes'"
echo "  After:  qxub --cmd \"clean command with \\\"readable quotes\\\"\""
