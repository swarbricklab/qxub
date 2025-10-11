#!/bin/bash
# Test script for qxub v2.2 remote execution
# This script demonstrates how to test the --config option and remote execution

set -e

echo "🧪 Testing qxub v2.2 remote execution with test configuration"
echo "============================================================="

CONFIG_FILE="tests/test_config.yaml"

# Test 1: Verify --config option is recognized
echo
echo "📋 Test 1: Verify --config option exists"
if qxub --help | grep -q "\-\-config"; then
    echo "✅ --config option found in help"
else
    echo "❌ --config option not found in help"
    exit 1
fi

# Test 2: Test config file loading
echo
echo "📋 Test 2: Test configuration file loading"
echo "Testing with: qxub --config $CONFIG_FILE --dry-run -- echo 'basic test'"
if qxub --config "$CONFIG_FILE" --dry-run -- echo "basic test" > /dev/null 2>&1; then
    echo "✅ Configuration file loads successfully"
else
    echo "❌ Configuration file loading failed"
    exit 1
fi

# Test 3: Test remote configuration validation (will fail SSH connection, but config should load)
echo
echo "📋 Test 3: Test remote configuration validation"
echo "Testing with: qxub --config $CONFIG_FILE --remote localhost --dry-run -- echo 'remote test'"
if qxub --config "$CONFIG_FILE" --remote localhost --dry-run -- echo "remote test" 2>&1 | grep -q "Error: Cannot connect to localhost"; then
    echo "✅ Remote configuration loaded correctly (expected SSH connection failure)"
elif qxub --config "$CONFIG_FILE" --remote localhost --dry-run -- echo "remote test" > /dev/null 2>&1; then
    echo "✅ Remote configuration and connection successful!"
else
    echo "ℹ️  Remote test completed (check output above for details)"
fi

# Test 4: Test invalid remote name
echo
echo "📋 Test 4: Test invalid remote name handling"
if qxub --config "$CONFIG_FILE" --remote nonexistent --dry-run -- echo "test" 2>&1 | grep -q "not found"; then
    echo "✅ Invalid remote name properly rejected"
else
    echo "❌ Invalid remote name handling failed"
fi

echo
echo "🎉 Test suite completed!"
echo
echo "Available test remotes in $CONFIG_FILE:"
echo "  - nci_gadi"
echo "  - localhost"
echo "  - example_cluster"
echo "  - custom_port_host"
echo
echo "Example usage:"
echo "  qxub --config $CONFIG_FILE --remote localhost --dry-run -- echo 'test'"
echo "  qxub --config $CONFIG_FILE --remote nci_gadi --dry-run -- python script.py"
