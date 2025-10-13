#!/bin/bash

# Test script for --cmd option functionality
# Tests the simplified variable substitution system with ${var} and ${{var}} syntax

set -e

echo "🧪 Testing --cmd option with simplified variable substitution"
echo

# Activate qxub environment
source venv/bin/activate

# Test 1: Basic submission-time variable expansion
echo "Test 1: Submission-time variable expansion"
result=$(qxub --env base --dry --cmd 'echo "Hello ${USER}"' 2>&1)
if echo "$result" | grep -q "Hello jr9959"; then
    echo "✅ PASS: ${USER} expanded to jr9959"
else
    echo "❌ FAIL: ${USER} not expanded correctly"
    echo "Output: $result"
fi
echo

# Test 2: Execution-time variable preparation
echo "Test 2: Execution-time variable preparation"
result=$(qxub --env base --dry --cmd 'echo "Job ${{PBS_JOBID}} running"' 2>&1)
if echo "$result" | grep -q "Job \${PBS_JOBID} running"; then
    echo "✅ PASS: \${{PBS_JOBID}} converted to \${PBS_JOBID}"
else
    echo "❌ FAIL: \${{PBS_JOBID}} not converted correctly"
    echo "Output: $result"
fi
echo

# Test 3: Mixed variables
echo "Test 3: Mixed submission and execution variables"
result=$(qxub --env base --dry --cmd 'echo "User ${USER} running job ${{PBS_JOBID}}"' 2>&1)
if echo "$result" | grep -q "User jr9959 running job \${PBS_JOBID}"; then
    echo "✅ PASS: Mixed variables handled correctly"
else
    echo "❌ FAIL: Mixed variables not handled correctly"
    echo "Output: $result"
fi
echo

# Test 4: Literal dollar signs preserved
echo "Test 4: Literal dollar signs preserved"
result=$(qxub --env base --dry --cmd 'python -c "print(\"Cost: \$100\")"' 2>&1)
if echo "$result" | grep -q '\$100'; then
    echo "✅ PASS: Literal \$100 preserved"
else
    echo "❌ FAIL: Literal dollar signs not preserved"
    echo "Output: $result"
fi
echo

# Test 5: AWK field references preserved
echo "Test 5: AWK field references preserved"
result=$(qxub --env base --dry --cmd 'awk "{print \$1, \$2}" file.txt' 2>&1)
if echo "$result" | grep -q '\$1.*\$2'; then
    echo "✅ PASS: AWK field references preserved"
else
    echo "❌ FAIL: AWK field references not preserved"
    echo "Output: $result"
fi
echo

# Test 6: Complex Python with f-strings
echo "Test 6: Complex Python command with f-strings"
result=$(qxub --env base --dry --cmd 'python -c "import os; print(f\"User: ${USER}, Job: ${{PBS_JOBID}}\")"' 2>&1)
if echo "$result" | grep -q "User: jr9959, Job: \${PBS_JOBID}" && echo "$result" | grep -q "f\\\\\""; then
    echo "✅ PASS: Complex Python with f-strings handled correctly"
else
    echo "❌ FAIL: Complex Python command not handled correctly"
    echo "Output: $result"
fi
echo

# Test 7: Error handling - both --cmd and -- specified
echo "Test 7: Error handling for both --cmd and --"
if result=$(qxub --env base --cmd "echo hello" -- echo world 2>&1); then
    echo "❌ FAIL: Should have exited with error"
    echo "Output: $result"
else
    if echo "$result" | grep -q "Cannot specify both"; then
        echo "✅ PASS: Error correctly raised for both --cmd and --"
    else
        echo "❌ FAIL: Wrong error message"
        echo "Output: $result"
    fi
fi
echo

# Test 8: Non-existent variable warning
echo "Test 8: Non-existent variable handling"
result=$(qxub --env base --dry -v --cmd 'echo "${NONEXISTENT_VAR}"' 2>&1)
if echo "$result" | grep -q "WARNING.*NONEXISTENT_VAR.*not found"; then
    echo "✅ PASS: Warning shown for non-existent variable"
else
    echo "❌ FAIL: No warning for non-existent variable"
    echo "Output: $result"
fi
echo

# Test 9: Verify traditional -- syntax still works
echo "Test 9: Traditional -- syntax compatibility"
result=$(qxub --env base --dry -- echo "Traditional syntax" 2>&1)
if echo "$result" | grep -q "Traditional syntax"; then
    echo "✅ PASS: Traditional -- syntax still works"
else
    echo "❌ FAIL: Traditional -- syntax broken"
    echo "Output: $result"
fi
echo

echo "🎯 All tests completed!"
echo
echo "Variable substitution rules:"
echo "  \${var}     → Expanded at submission time (config/user context)"
echo "  \${{var}}   → Converted to \${var} for execution time (job context)"
echo "  \$other     → Preserved unchanged (literal usage)"
