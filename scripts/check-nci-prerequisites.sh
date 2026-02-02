#!/bin/bash
# NCI Gadi Prerequisites Check
# Run this script on NCI Gadi to verify prerequisites for qxub remote execution

set -e

echo "🔍 Checking NCI Gadi prerequisites for qxub remote execution..."
echo

# Check conda installation
echo "1. Checking conda installation..."
if command -v conda >/dev/null 2>&1; then
    echo "   ✅ conda found: $(conda --version)"
    conda info --base
else
    echo "   ❌ conda not found in PATH"
    exit 1
fi

# Check available conda environments
echo
echo "2. Checking conda environments..."
conda env list | head -20

# Check if qxub is available in base environment
echo
echo "3. Checking qxub in base environment..."
if conda run -n base qxub --version >/dev/null 2>&1; then
    echo "   ✅ qxub available in base environment: $(conda run -n base qxub --version)"
else
    echo "   ❌ qxub not available in base environment"
    echo "   💡 You may need to install qxub or specify a different environment"
fi

# Check common qxub environments
echo
echo "4. Checking for common qxub environments..."
for env in qxub qxub-prod qxub-dev; do
    if conda env list | grep -q "^$env "; then
        if conda run -n $env qxub --version >/dev/null 2>&1; then
            echo "   ✅ qxub available in $env: $(conda run -n $env qxub --version)"
        else
            echo "   ⚠️  Environment $env exists but qxub not working"
        fi
    else
        echo "   ➖ Environment $env not found"
    fi
done

# Check platform file locations
echo
echo "5. Checking platform file locations..."
platform_paths=(
    "/g/data/a56/software/qsub_tools/docs/platforms/nci_gadi.yaml"
    "/etc/xdg/qxub/platforms/nci_gadi.yaml"
    "~/.config/qxub/platforms/nci_gadi.yaml"
    "~/.local/share/qxub/platforms/nci_gadi.yaml"
)

for path in "${platform_paths[@]}"; do
    expanded_path=$(eval echo "$path")
    if [ -f "$expanded_path" ]; then
        echo "   ✅ Platform file found: $expanded_path"
        echo "      File size: $(stat -c%s "$expanded_path") bytes"
    else
        echo "   ➖ Platform file not found: $expanded_path"
    fi
done

# Check project directories
echo
echo "6. Checking project directory structure..."
if [ -n "$PROJECT" ] && [ -n "$USER" ]; then
    scratch_dir="/scratch/$PROJECT/$USER"
    echo "   PROJECT=$PROJECT, USER=$USER"

    if [ -d "$scratch_dir" ]; then
        echo "   ✅ Scratch directory exists: $scratch_dir"

        # Check if we can create test directory
        test_dir="$scratch_dir/ci-test"
        if mkdir -p "$test_dir" 2>/dev/null; then
            echo "   ✅ Can create test directory: $test_dir"
            rmdir "$test_dir" 2>/dev/null || true
        else
            echo "   ⚠️  Cannot create test directory in $scratch_dir"
        fi
    else
        echo "   ❌ Scratch directory not found: $scratch_dir"
    fi
else
    echo "   ⚠️  PROJECT and/or USER environment variables not set"
    echo "      PROJECT=${PROJECT:-'(not set)'}"
    echo "      USER=${USER:-'(not set)'}"
fi

# Test SSH loopback (for testing remote execution locally)
echo
echo "7. Testing SSH loopback connection..."
if ssh -o ConnectTimeout=5 -o StrictHostKeyChecking=no localhost echo "SSH loopback works" >/dev/null 2>&1; then
    echo "   ✅ SSH loopback connection works"
else
    echo "   ❌ SSH loopback connection failed"
    echo "   💡 This is needed for testing remote execution on the same system"
fi

# Check qxub installation method
echo
echo "8. Checking qxub installation details..."
if command -v qxub >/dev/null 2>&1; then
    qxub_path=$(which qxub)
    echo "   qxub path: $qxub_path"

    if [ -L "$qxub_path" ]; then
        echo "   Real path: $(readlink -f "$qxub_path")"
    fi

    # Check if it's a conda installation
    if echo "$qxub_path" | grep -q conda; then
        echo "   ✅ qxub installed via conda"
    elif echo "$qxub_path" | grep -q pip; then
        echo "   ✅ qxub installed via pip"
    else
        echo "   ℹ️  qxub installation method unclear"
    fi
else
    echo "   ❌ qxub not found in PATH"
fi

echo
echo "🏁 Prerequisites check completed!"
echo
echo "📋 Summary for CI configuration:"
echo "   - Recommended qxub_env: base (or specify a working environment from above)"
echo "   - Platform file: Use a found platform file path from above"
echo "   - Test PROJECT: ${PROJECT:-'(set PROJECT environment variable)'}"
echo "   - Test USER: ${USER:-'(set USER environment variable)'}"
