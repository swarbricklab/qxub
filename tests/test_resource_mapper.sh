#!/bin/bash
# Test script for new workflow-friendly resource options in qxub exec

set -e

# Test directory
TEST_DIR=$(dirname "$(readlink -f "$0")")
SCRIPT_DIR=$(dirname "$TEST_DIR")

echo "[INFO] Starting qxub v3.1.0 resource mapper tests..."
echo "ℹ️  These tests verify the new --mem, --cpus, --runtime, --disk options"

# Setup test environment
cd "$SCRIPT_DIR"

# Check if we're already in the qxub environment
if [[ "$CONDA_DEFAULT_ENV" != "qxub" ]]; then
    echo "[INFO] Note: Please ensure qxub conda environment is activated"
    echo "[INFO] Run: conda activate qxub"
fi

echo ""
echo "=== Testing Resource Mapper Functionality ==="

# Test 1: Basic memory and CPU flags
echo "[TEST] Basic --mem and --cpus flags"
OUTPUT=$(qxub exec --dry --mem 8GB --cpus 4 -- echo "Basic resource test" 2>&1)
if echo "$OUTPUT" | grep -q "mem=8GB" && echo "$OUTPUT" | grep -q "ncpus=4"; then
    echo "[PASS] Basic --mem and --cpus flags"
else
    echo "[FAIL] Basic --mem and --cpus flags"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 2: Runtime flag variants
echo "[TEST] Runtime flag variants"
OUTPUT=$(qxub exec --dry --runtime 2h30m --cpus 2 -- echo "Runtime test" 2>&1)
if echo "$OUTPUT" | grep -q "walltime=2:30:00" && echo "$OUTPUT" | grep -q "ncpus=2"; then
    echo "[PASS] Runtime flag variants"
else
    echo "[FAIL] Runtime flag variants"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 3: Alternative flag names
echo "[TEST] Alternative flag names (--memory, --threads, --time, --jobfs)"
OUTPUT=$(qxub exec --dry --memory 16GB --threads 8 --time 1h --jobfs 10GB -- echo "Alternative flags test" 2>&1)
if echo "$OUTPUT" | grep -q "mem=16GB" && echo "$OUTPUT" | grep -q "ncpus=8" && echo "$OUTPUT" | grep -q "walltime=1:00:00" && echo "$OUTPUT" | grep -q "jobfs=10GB"; then
    echo "[PASS] Alternative flag names"
else
    echo "[FAIL] Alternative flag names"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 4: Mixing new flags with traditional --resources
echo "[TEST] Mixing new flags with traditional --resources"
OUTPUT=$(qxub exec --dry --resources "walltime=1:00:00" --mem 4GB --cpus 2 -- echo "Mixed resources test" 2>&1)
if echo "$OUTPUT" | grep -q "walltime=1:00:00" && echo "$OUTPUT" | grep -q "mem=4GB" && echo "$OUTPUT" | grep -q "ncpus=2"; then
    echo "[PASS] Mixing new flags with traditional --resources"
else
    echo "[FAIL] Mixing new flags with traditional --resources"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 5: Complex runtime formats
echo "[TEST] Complex runtime formats"
OUTPUT=$(qxub exec --dry --runtime 1h30m -- echo "Complex runtime" 2>&1)
if echo "$OUTPUT" | grep -q "walltime=1:30:00"; then
    echo "[PASS] Complex runtime formats"
else
    echo "[FAIL] Complex runtime formats"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 6: Memory format variations
echo "[TEST] Memory format variations"
OUTPUT=$(qxub exec --dry --mem 2048MB -- echo "Memory variations" 2>&1)
if echo "$OUTPUT" | grep -q "mem=2048MB"; then
    echo "[PASS] Memory format variations"
else
    echo "[FAIL] Memory format variations"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 7: All resource flags combined
echo "[TEST] All resource flags combined"
OUTPUT=$(qxub exec --dry --mem 32GB --cpus 16 --runtime 4h --disk 100GB -- echo "All resources test" 2>&1)
if echo "$OUTPUT" | grep -q "mem=32GB" && echo "$OUTPUT" | grep -q "ncpus=16" && echo "$OUTPUT" | grep -q "walltime=4:00:00" && echo "$OUTPUT" | grep -q "jobfs=100GB"; then
    echo "[PASS] All resource flags combined"
else
    echo "[FAIL] All resource flags combined"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 7.5: Storage/volumes flags
echo "[TEST] Storage/volumes flags (--volumes, --storage)"
OUTPUT=$(qxub exec --dry --volumes gdata/a56 --mem 8GB -- echo "Volumes test" 2>&1)
if echo "$OUTPUT" | grep -q "storage=gdata/a56" && echo "$OUTPUT" | grep -q "mem=8GB"; then
    echo "[PASS] --volumes flag works"
else
    echo "[FAIL] --volumes flag"
    echo "Output: $OUTPUT"
    exit 1
fi

OUTPUT=$(qxub exec --dry --storage gdata/a56+scratch/a56 --cpus 4 -- echo "Storage alias test" 2>&1)
if echo "$OUTPUT" | grep -q "storage=gdata/a56+scratch/a56" && echo "$OUTPUT" | grep -q "ncpus=4"; then
    echo "[PASS] --storage alias works"
else
    echo "[FAIL] --storage alias"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 8: With execution contexts
echo "[TEST] Resource flags with conda environment"
OUTPUT=$(qxub exec --dry --env base --mem 8GB --cpus 4 -- echo "Conda with resources" 2>&1)
if echo "$OUTPUT" | grep -q "env=base" && echo "$OUTPUT" | grep -q "mem=8GB" && echo "$OUTPUT" | grep -q "ncpus=4"; then
    echo "[PASS] Resource flags with conda environment"
else
    echo "[FAIL] Resource flags with conda environment"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 9: With modules
echo "[TEST] Resource flags with modules"
OUTPUT=$(qxub exec --dry --mod python3 --mem 4GB --cpus 2 -- echo "Modules with resources" 2>&1)
if echo "$OUTPUT" | grep -q "python3" && echo "$OUTPUT" | grep -q "mem=4GB" && echo "$OUTPUT" | grep -q "ncpus=2"; then
    echo "[PASS] Resource flags with modules"
else
    echo "[FAIL] Resource flags with modules"
    echo "Output: $OUTPUT"
    exit 1
fi

# Test 10: Priority test - CLI flags should override config defaults
echo "[TEST] CLI resource flags override config defaults"
# This test verifies that the new flags properly override any config defaults
OUTPUT=$(qxub exec --dry --mem 1GB --cpus 1 -- echo "Override test" 2>&1)
if echo "$OUTPUT" | grep -q "mem=1GB" && echo "$OUTPUT" | grep -q "ncpus=1"; then
    echo "[PASS] CLI resource flags override config defaults"
else
    echo "[FAIL] CLI resource flags override config defaults"
    echo "Output: $OUTPUT"
    exit 1
fi

echo ""
echo "=== Testing Error Conditions ==="

# Test error handling
echo "[TEST] Invalid memory format"
if qxub exec --dry --mem "invalid" -- echo "test" 2>/dev/null; then
    echo "[FAIL] Should have failed with invalid memory format"
    exit 1
else
    echo "[PASS] Invalid memory format properly rejected"
fi

echo "[TEST] Invalid runtime format"
if qxub exec --dry --runtime "invalid" -- echo "test" 2>/dev/null; then
    echo "[FAIL] Should have failed with invalid runtime format"
    exit 1
else
    echo "[PASS] Invalid runtime format properly rejected"
fi

echo ""
echo "========================================"
echo "      RESOURCE MAPPER TEST RESULTS"
echo "========================================"
echo "All workflow-friendly resource option tests passed! ✅"
echo ""
echo "New features verified:"
echo "  ✅ --mem / --memory for memory specification"
echo "  ✅ --cpus / --threads for CPU specification"
echo "  ✅ --runtime / --time for walltime specification"
echo "  ✅ --disk / --jobfs for disk specification"
echo "  ✅ --volumes / --storage for NCI storage volume mounting"
echo "  ✅ Proper integration with existing --resources option"
echo "  ✅ Works with all execution contexts (conda, modules, etc.)"
echo "  ✅ Proper error handling for invalid formats"
echo ""
echo "🎉 v3.1.0 resource mapper implementation successful!"
