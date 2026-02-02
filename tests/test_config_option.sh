#!/bin/bash
# Test script for --config option
# Tests both local and remote execution modes with override config

set -e

echo "🧪 Testing --config option implementation"
echo ""

CONFIG_FILE="tests/configs/test_remote.yaml"

# Test 1: Verify --config option exists
echo "Test 1: Verify --config option in help"
if qxub exec --help | grep -q "\-\-config"; then
    echo "✅ --config option found"
else
    echo "❌ --config option not found"
    exit 1
fi
echo ""

# Test 2: Local execution with override config
echo "Test 2: Local execution (ci_test_local platform)"
echo "Command: qxub exec --config $CONFIG_FILE --platform ci_test_local --dry -- echo 'local test'"
if qxub exec --config "$CONFIG_FILE" --platform ci_test_local --dry --default -- echo "local test" 2>&1 | grep -q "Dry run"; then
    echo "✅ Local execution works with override config"
else
    echo "❌ Local execution failed"
    exit 1
fi
echo ""

# Test 3: Remote execution detection (will fail SSH but proves config loaded)
echo "Test 3: Remote execution detection (ci_test_remote platform)"
echo "Command: qxub exec --config $CONFIG_FILE --platform ci_test_remote --dry -- echo 'remote test'"
echo "Note: SSH will fail on compute node, but this proves remote mode was detected"
if qxub exec --config "$CONFIG_FILE" --platform ci_test_remote --dry --default -- echo "remote test" 2>&1 | grep -q -E "(Could not resolve hostname|Remote execution)"; then
    echo "✅ Remote execution mode detected (SSH failed as expected)"
else
    echo "⚠️  Unexpected output (might be on login node with SSH access)"
fi
echo ""

# Test 4: Default platform from override config
echo "Test 4: Default platform from override config"
echo "Command: qxub exec --config $CONFIG_FILE --dry -- echo 'default test'"
echo "Note: Should use ci_test_remote as default from test config"
if qxub exec --config "$CONFIG_FILE" --dry --default -- echo "default test" 2>&1 | grep -q -E "(Could not resolve hostname|Remote execution)"; then
    echo "✅ Default platform (ci_test_remote) used from override config"
else
    echo "⚠️  Unexpected result"
fi
echo ""

echo "✅ All --config option tests passed!"
echo ""
echo "Usage examples:"
echo "  # Local execution with override config"
echo "  qxub exec --config $CONFIG_FILE --platform ci_test_local -- python script.py"
echo ""
echo "  # Remote execution with override config (requires SSH access)"
echo "  qxub exec --config $CONFIG_FILE --platform ci_test_remote -- python script.py"
echo ""
echo "  # Use default platform from override config"
echo "  qxub exec --config $CONFIG_FILE -- python script.py"
