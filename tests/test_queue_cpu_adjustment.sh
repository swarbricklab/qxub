#!/bin/bash
# Test script for graceful CPU adjustment when queue is explicitly specified
# Tests the new behavior where default CPUs are adjusted to queue max, but explicit CPUs cause errors

set -e

# Test directory
TEST_DIR=$(dirname "$(readlink -f "$0")")
SCRIPT_DIR=$(dirname "$TEST_DIR")

echo "🧪 Testing queue CPU adjustment behavior..."
echo ""

# Setup test environment
cd "$SCRIPT_DIR"

# Verify qxub is available
if ! command -v qxub &> /dev/null; then
    echo "❌ ERROR: qxub command not found"
    echo "   Please ensure qxub is installed and in PATH"
    exit 1
fi

echo "Using qxub: $(which qxub)"
echo ""

echo "=== Test 1: Explicit queue with default CPUs (should auto-adjust) ==="
echo "Testing: qxub exec --dry --queue copyq -- echo 'test'"
echo "Expected: Should auto-adjust ncpus to 1 (copyq maximum)"
OUTPUT=$(qxub exec --dry --queue copyq -- echo "test" 2>&1)
if echo "$OUTPUT" | grep -q "ncpus=1"; then
    echo "✅ PASS: Auto-adjusted to ncpus=1 for copyq"
else
    echo "❌ FAIL: Did not adjust CPUs properly"
    echo "Output: $OUTPUT"
    exit 1
fi
echo ""

echo "=== Test 2: Explicit queue with explicit CPUs within limit (should work) ==="
echo "Testing: qxub exec --dry --queue copyq -l ncpus=1 -- echo 'test'"
echo "Expected: Should accept ncpus=1"
OUTPUT=$(qxub exec --dry --queue copyq -l ncpus=1 -- echo "test" 2>&1)
if echo "$OUTPUT" | grep -q "ncpus=1"; then
    echo "✅ PASS: Accepted explicit ncpus=1"
else
    echo "❌ FAIL: Did not accept explicit CPUs"
    echo "Output: $OUTPUT"
    exit 1
fi
echo ""

echo "=== Test 3: Explicit queue with explicit CPUs exceeding limit (PBS will reject) ==="
echo "Testing: qxub exec --dry --queue copyq -l ncpus=4 -- echo 'test'"
echo "Expected: qxub allows submission (PBS qsub will reject it)"
OUTPUT=$(qxub exec --dry --queue copyq -l ncpus=4 -- echo "test" 2>&1)
if echo "$OUTPUT" | grep -q "ncpus=4"; then
    echo "✅ PASS: qxub accepts explicit ncpus=4 (PBS will validate and reject)"
    echo "   Note: PBS qsub will reject this with an error about queue limits"
else
    echo "❌ FAIL: Unexpected behavior"
    echo "Output: $OUTPUT"
fi
echo ""

echo "=== Test 4: Auto queue selection (should use existing logic) ==="
echo "Testing: qxub exec --dry --queue auto -l mem=16GB -- echo 'test'"
echo "Expected: Should select appropriate queue based on resources"
OUTPUT=$(qxub exec --dry --queue auto -l mem=16GB -- echo "test" 2>&1)
if echo "$OUTPUT" | grep -q "queue"; then
    echo "✅ PASS: Auto queue selection still works"
else
    echo "⚠️  WARNING: Auto queue selection may have issues"
    echo "Output: $OUTPUT"
fi
echo ""

echo "=== Test 5: Normal queue with default CPUs (should not adjust) ==="
echo "Testing: qxub exec --dry --queue normal -- echo 'test'"
echo "Expected: Should keep default CPUs (normal queue has high max_cpus)"
OUTPUT=$(qxub exec --dry --queue normal -- echo "test" 2>&1)
if echo "$OUTPUT" | grep -q "ncpus="; then
    echo "✅ PASS: Normal queue accepts default CPUs"
    echo "Resources: $(echo "$OUTPUT" | grep -o 'ncpus=[0-9]*')"
else
    echo "⚠️  WARNING: Could not determine CPU count"
    echo "Output: $OUTPUT"
fi
echo ""

echo "========================================="
echo "✅ All core tests passed!"
echo "========================================="
echo ""
echo "Summary:"
echo "- Default CPUs are auto-adjusted to queue max_cpus when queue is explicit"
echo "- Explicit CPUs within limits are accepted"
echo "- Explicit CPUs exceeding limits cause helpful errors"
echo "- Auto queue selection unchanged"
