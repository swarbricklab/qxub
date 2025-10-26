#!/bin/bash

# Test workflow-friendly config defaults functionality
# This test verifies that workflow-friendly options (--mem, --cpus, etc.)
# can have defaults set in config and are properly overridden by CLI args

set -e

echo "🧪 Testing workflow-friendly config defaults..."

#!/bin/bash

# Test workflow-friendly config defaults functionality
# This test verifies that workflow-friendly options (--mem, --cpus, etc.)
# can have defaults set in config and are properly overridden by CLI args

set -e

echo "🧪 Testing workflow-friendly config defaults..."

# Note: This test assumes a clean config or ignores existing values
# qxub doesn't currently have an 'unset' command, so we test functionality as-is

# Test 1: Set workflow-friendly config defaults
echo "🔧 Test 1: Setting workflow-friendly config defaults..."
qxub config set mem "8GB"
qxub config set cpus 4
qxub config set runtime "2h"
qxub config set disk "15GB"
qxub config set volumes "gdata/a56+gdata/px14"

# Test 2: Verify config defaults are used
echo "🔍 Test 2: Checking config defaults are applied..."
OUTPUT=$(qxub exec --dry -- echo "test-config-defaults")
if echo "$OUTPUT" | grep -q "mem=8GB" && \
   echo "$OUTPUT" | grep -q "ncpus=4" && \
   echo "$OUTPUT" | grep -q "walltime=2:00:00" && \
   echo "$OUTPUT" | grep -q "jobfs=15GB" && \
   echo "$OUTPUT" | grep -q "storage=gdata/a56+gdata/px14"; then
    echo "✅ All config defaults are applied correctly"
else
    echo "❌ Config defaults not working correctly"
    echo "$OUTPUT"
    exit 1
fi

# Test 3: CLI options override config defaults
echo "🔍 Test 3: CLI options override config defaults..."
OUTPUT=$(qxub exec --dry --mem 16GB --cpus 8 --runtime 4h -- echo "test-cli-override")
if echo "$OUTPUT" | grep -q "mem=16GB" && \
   echo "$OUTPUT" | grep -q "ncpus=8" && \
   echo "$OUTPUT" | grep -q "walltime=4:00:00" && \
   echo "$OUTPUT" | grep -q "jobfs=15GB" && \
   echo "$OUTPUT" | grep -q "storage=gdata/a56+gdata/px14"; then
    echo "✅ CLI options correctly override config defaults"
else
    echo "❌ CLI override not working correctly"
    echo "$OUTPUT"
    exit 1
fi

# Test 4: Alternative option names work with config
echo "🔍 Test 4: Alternative option names work..."
OUTPUT=$(qxub exec --dry --memory 32GB --threads 16 --time 6h --jobfs 30GB --storage "gdata/xyz" -- echo "test-alternatives")
if echo "$OUTPUT" | grep -q "mem=32GB" && \
   echo "$OUTPUT" | grep -q "ncpus=16" && \
   echo "$OUTPUT" | grep -q "walltime=6:00:00" && \
   echo "$OUTPUT" | grep -q "jobfs=30GB" && \
   echo "$OUTPUT" | grep -q "storage=gdata/xyz"; then
    echo "✅ Alternative option names work correctly"
else
    echo "❌ Alternative option names not working"
    echo "$OUTPUT"
    exit 1
fi

echo "✅ All workflow-friendly config defaults tests passed!"
echo "💡 Note: To clean up, use 'qxub config reset' to remove test values"
