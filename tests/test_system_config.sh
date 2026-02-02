#!/bin/bash

# Test script for system-level configuration and config precedence
# Tests the hierarchical config system: CLI > User > System > Defaults

set -e  # Exit on error

echo "🧪 Starting system config and precedence tests..."

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to run a test
run_test() {
    local test_name="$1"
    local command="$2"
    local expected_pattern="$3"

    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    echo -e "${YELLOW}[INFO]${NC} Running test: $test_name"
    echo "Command: $command"

    if output=$(eval "$command" 2>&1); then
        if [[ -z "$expected_pattern" ]] || echo "$output" | grep -q "$expected_pattern"; then
            echo -e "${GREEN}[PASS]${NC} $test_name"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            echo -e "${RED}[FAIL]${NC} $test_name - Expected pattern '$expected_pattern' not found"
            echo "Output: $output"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        echo -e "${RED}[FAIL]${NC} $test_name - Command failed"
        echo "Output: $output"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    echo "Output: $output"
    echo "---"
}

# Setup test environment
echo -e "${YELLOW}[INFO]${NC} Setting up test environment..."

# Create temporary directories for testing
TEST_DIR=$(mktemp -d)
SYSTEM_CONFIG_DIR="$TEST_DIR/etc/xdg/qxub"
USER_CONFIG_DIR="$TEST_DIR/home/.config/qxub"

mkdir -p "$SYSTEM_CONFIG_DIR"
mkdir -p "$USER_CONFIG_DIR"

# Backup existing configs if they exist
USER_CONFIG_HOME="${XDG_CONFIG_HOME:-$HOME/.config}/qxub"
if [[ -d "$USER_CONFIG_HOME" ]]; then
    echo -e "${YELLOW}[INFO]${NC} Backing up existing user config"
    cp -r "$USER_CONFIG_HOME" "$TEST_DIR/backup_user_config"
fi

echo -e "${YELLOW}[INFO]${NC} === Testing System Config Discovery ==="

# Test 1: System config file detection
cat > "$SYSTEM_CONFIG_DIR/config.yaml" << 'EOF'
defaults:
  project: "system_project"
  queue: "system_queue"
  resources:
    - "walltime=02:00:00"
    - "mem=4GB"
    - "ncpus=2"

templates:
  system_template: "system_value"
  shared_template: "system_shared"

aliases:
  system_alias:
    main:
      name: "system_job"
      queue: "copyq"
    subcommand:
      type: conda
      env: "system_env"
    target:
      cmd: "echo 'System alias executed'"
EOF

# Test with XDG_CONFIG_DIRS pointing to our test directory
export XDG_CONFIG_DIRS="$TEST_DIR/etc/xdg"
export XDG_CONFIG_HOME="$TEST_DIR/home/.config"

run_test "System config discovery" \
    "qxub config get defaults.project" \
    "system_project"

run_test "System config queue setting" \
    "qxub config get defaults.queue" \
    "system_queue"

run_test "System template variables" \
    "qxub config get templates.system_template" \
    "system_value"

run_test "System alias listing" \
    "qxub config alias list-aliases" \
    "system_alias"

echo -e "${YELLOW}[INFO]${NC} === Testing User Config Override ==="

# Create user config that overrides some system settings
cat > "$USER_CONFIG_DIR/config.yaml" << 'EOF'
defaults:
  project: "user_project"  # Override system
  name: "user_job_{date}"  # New setting
  resources:
    - "walltime=01:00:00"  # Override system
    - "mem=8GB"            # Override system

templates:
  user_template: "user_value"
  shared_template: "user_shared"  # Override system

aliases:
  user_alias:
    main:
      name: "user_job"
      queue: "normal"
    subcommand:
      type: conda
      env: "user_env"
    target:
      cmd: "echo 'User alias executed'"

  system_alias:  # Override system alias
    main:
      name: "overridden_system_job"
      queue: "express"
    subcommand:
      type: conda
      env: "user_override_env"
    target:
      cmd: "echo 'User overridden system alias'"
EOF

# Force reload configs
run_test "User overrides system project" \
    "qxub config get defaults.project" \
    "user_project"

run_test "User inherits system queue (not overridden)" \
    "qxub config get defaults.queue" \
    "system_queue"

run_test "User overrides system template" \
    "qxub config get templates.shared_template" \
    "user_shared"

run_test "User adds new template" \
    "qxub config get templates.user_template" \
    "user_value"

run_test "System template still accessible" \
    "qxub config get templates.system_template" \
    "system_value"

run_test "User and system aliases both available" \
    "qxub config alias list-aliases" \
    "user_alias.*system_alias"

echo -e "${YELLOW}[INFO]${NC} === Testing CLI Override Precedence ==="

# Test CLI arguments override both user and system config
run_test "CLI overrides user project" \
    "qxub --project cli_project config get defaults.project" \
    "cli_project"

run_test "CLI overrides user queue" \
    "qxub --queue cli_queue config get defaults.queue" \
    "cli_queue"

# Test with actual command execution (dry run)
run_test "CLI + User + System precedence (dry run)" \
    "qxub --dry-run --queue cli_queue --project cli_project conda --env test_env echo 'test'" \
    "Project: cli_project.*Queue: cli_queue"

echo -e "${YELLOW}[INFO]${NC} === Testing Alias Precedence ==="

# Test that user alias overrides system alias with same name
run_test "User alias overrides system alias (show)" \
    "qxub config alias show system_alias" \
    "overridden_system_job"

# Test alias execution with overrides
run_test "Alias inherits system config defaults" \
    "qxub --dry-run alias user_alias" \
    "Queue: system_queue"  # Should inherit from system since user didn't override

echo -e "${YELLOW}[INFO]${NC} === Testing Template Variable Precedence ==="

# Create a test that uses template variables from different levels
cat >> "$USER_CONFIG_DIR/config.yaml" << 'EOF'

aliases:
  template_test:
    main:
      name: "{user_template}_{system_template}_{shared_template}"
    subcommand:
      type: conda
      env: "test"
    target:
      cmd: "echo 'template test'"
EOF

run_test "Template variable inheritance in aliases" \
    "qxub --dry-run alias template_test" \
    "user_value.*system_value.*user_shared"

echo -e "${YELLOW}[INFO]${NC} === Testing Multiple System Config Directories ==="

# Test multiple directories in XDG_CONFIG_DIRS
SYSTEM_CONFIG_DIR2="$TEST_DIR/etc/xdg2/qxub"
mkdir -p "$SYSTEM_CONFIG_DIR2"

cat > "$SYSTEM_CONFIG_DIR2/config.yaml" << 'EOF'
defaults:
  queue: "secondary_queue"

templates:
  secondary_template: "secondary_value"

aliases:
  secondary_alias:
    main:
      name: "secondary_job"
    subcommand:
      type: conda
      env: "secondary_env"
    target:
      cmd: "echo 'Secondary system alias'"
EOF

# Test with multiple system config directories (first takes precedence)
export XDG_CONFIG_DIRS="$TEST_DIR/etc/xdg:$TEST_DIR/etc/xdg2"

run_test "First system config takes precedence" \
    "qxub config get defaults.queue" \
    "system_queue"  # Should be from first dir, not secondary_queue

run_test "Secondary system config contributes aliases" \
    "qxub config alias list-aliases" \
    "secondary_alias"

echo -e "${YELLOW}[INFO]${NC} === Testing Error Handling ==="

# Test with invalid system config
INVALID_SYSTEM_DIR="$TEST_DIR/etc/xdg_invalid/qxub"
mkdir -p "$INVALID_SYSTEM_DIR"

cat > "$INVALID_SYSTEM_DIR/config.yaml" << 'EOF'
invalid_yaml: [unclosed bracket
  another_line: with_error
EOF

export XDG_CONFIG_DIRS="$TEST_DIR/etc/xdg_invalid:$TEST_DIR/etc/xdg"

# Should warn about invalid config but continue with valid ones
run_test "Invalid system config handling" \
    "qxub config get defaults.project 2>&1" \
    "Warning.*Failed to load system config"

echo -e "${YELLOW}[INFO]${NC} === Testing Config File Discovery ==="

run_test "Config files discovery" \
    "qxub config files" \
    "system_.*$SYSTEM_CONFIG_DIR/config.yaml"

echo -e "${YELLOW}[INFO]${NC} Cleaning up test environment..."

# Restore original environment
unset XDG_CONFIG_DIRS
unset XDG_CONFIG_HOME

# Restore original user config if it existed
if [[ -d "$TEST_DIR/backup_user_config" ]]; then
    echo -e "${YELLOW}[INFO]${NC} Restoring original user config"
    USER_CONFIG_HOME="${XDG_CONFIG_HOME:-$HOME/.config}/qxub"
    rm -rf "$USER_CONFIG_HOME"
    cp -r "$TEST_DIR/backup_user_config" "$USER_CONFIG_HOME"
fi

# Clean up test directory
rm -rf "$TEST_DIR"

echo "========================================"
echo "          SYSTEM CONFIG TEST RESULTS"
echo "========================================"
echo "Total Tests: $TOTAL_TESTS"
echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed: ${RED}$FAILED_TESTS${NC}"

if [[ $FAILED_TESTS -eq 0 ]]; then
    echo -e "${GREEN}🎉 All system config tests passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ Some tests failed. Please review the output above.${NC}"
    exit 1
fi
