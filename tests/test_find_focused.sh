#!/bin/bash

# Simplified find command test focused on key edge cases
# This version handles shell quoting issues more carefully

set -e

echo "🔍 Testing critical GNU find command patterns with qxub"
echo

# Create test directory structure
TEST_DIR="/tmp/qxub_find_test_$$"
mkdir -p "$TEST_DIR"/{subdir1,subdir2}
touch "$TEST_DIR"/{file1.txt,file2.py,file3.log}
touch "$TEST_DIR"/subdir1/{nested1.txt,nested2.py}
echo "some content" > "$TEST_DIR/file1.txt"

echo "📁 Created test directory: $TEST_DIR"
echo

# Cleanup function
cleanup() {
    echo "🧹 Cleaning up test directory..."
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

# Test 1: Basic find with -exec and \; (the classic problem case)
echo "Test 1: Basic find with -exec and semicolon terminator"
echo "  Traditional syntax:"
qxub --env base --dry -- find "$TEST_DIR" -name "*.txt" -exec echo "Found:" {} \;
echo "  --cmd syntax:"
qxub --env base --dry --cmd "find $TEST_DIR -name \"*.txt\" -exec echo \"Found:\" {} \;"
echo

# Test 2: Find with complex shell variables and escaping
echo "Test 2: Find with shell variables and complex escaping"
echo "  Traditional syntax (problematic):"
qxub --env base --dry -- find "$TEST_DIR" -name "*.txt" -exec sh -c 'echo "File: $1, User: $USER"' _ {} \;
echo "  --cmd syntax (clean):"
qxub --env base --dry --cmd 'find '"$TEST_DIR"' -name "*.txt" -exec sh -c '"'"'echo "File: $1, User: ${USER}"'"'"' _ {} \;'
echo

# Test 3: Find with qxub variable substitution
echo "Test 3: Find with qxub submission and execution time variables"
echo "  --cmd syntax with mixed variables:"
qxub --env base --dry --cmd 'find ${HOME} -maxdepth 1 -name ".*" -exec echo "User ${USER} found ${{PBS_JOBID}}: {}" \;'
echo

# Test 4: Find with AWK-style field references (should be preserved)
echo "Test 4: Find with AWK field references (literal \$ preservation)"
echo "  --cmd syntax:"
qxub --env base --dry --cmd 'find '"$TEST_DIR"' -name "*.txt" -exec awk '"'"'{print "Line " NR ": " $0}'"'"' {} \;'
echo

# Test 5: Find with JSON-like structure (complex quoting)
echo "Test 5: Find with JSON-like structure in -exec"
echo "  --cmd syntax:"
qxub --env base --dry --cmd 'find '"$TEST_DIR"' -name "*.txt" -exec sh -c '"'"'echo "{\"file\": \"$1\", \"user\": \"${USER}\", \"cost\": \"\$50\"}"'"'"' _ {} \;'
echo

# Test 6: Find with regex patterns and special characters
echo "Test 6: Find with regex patterns"
echo "  Traditional syntax:"
qxub --env base --dry -- find "$TEST_DIR" -regex '.*\.\(txt\|py\)'
echo "  --cmd syntax:"
qxub --env base --dry --cmd 'find '"$TEST_DIR"' -regex '"'"'.*\.\(txt\|py\)'"'"
echo

# Test 7: Find with -print0 and xargs (pipe handling)
echo "Test 7: Find with pipes to xargs"
echo "  --cmd syntax:"
qxub --env base --dry --cmd 'find '"$TEST_DIR"' -name "*.txt" -print0 | xargs -0 wc -l'
echo

# Test 8: Find with glob patterns that should NOT be shell-expanded
echo "Test 8: Find with glob patterns (literal preservation)"
echo "  Traditional syntax (shell expands):"
qxub --env base --dry -- find "$TEST_DIR" -name '[fn]*'
echo "  --cmd syntax (literal preserved):"
qxub --env base --dry --cmd 'find '"$TEST_DIR"' -name '"'"'[fn]*'"'"
echo

# Real execution test
echo "🚀 Real execution test:"
echo "Running a simple find command to verify it actually works..."
qxub --env base --cmd 'find '"$TEST_DIR"' -name "*.txt" -exec echo "Actually found: {}" \;'
echo

echo "🎯 Key Observations:"
echo "  1. Traditional -- syntax: Shell processes quotes and special chars first"
echo "  2. --cmd syntax: Better control over quoting and special characters"
echo "  3. Variable substitution: \${var} for submission, \${{var}} for execution"
echo "  4. AWK patterns: Literal \$ preserved in both syntaxes"
echo "  5. Complex quoting: --cmd syntax handles nested quotes better"
echo
echo "💡 Recommendation: Use --cmd for complex find operations with:"
echo "  • Multiple levels of quotes"
echo "  • Variables that need specific expansion timing"
echo "  • Patterns that shouldn't be shell-expanded"
echo "  • JSON or structured output formatting"
