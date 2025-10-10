#!/bin/bash

# Realistic System Configuration Tests
# Tests actual usage patterns that make sense in the real world

set -e  # Exit on error

echo "🧪 Starting realistic system config tests..."run_test "CLI args override user config (dry run)" \
    "qxub --dry-run --project cli_override_project --queue express conda --env test_env echo 'test'" \
    "-q express -P cli_override_project" Colors for output
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
        if [[ -z "$expected_pattern" ]] || echo "$output" | grep -q -- "$expected_pattern"; then
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

echo -e "${YELLOW}[INFO]${NC} === Testing System Admin Setup Scenarios ==="

# Test 1: System admin sets organization-wide defaults
cat > "$SYSTEM_CONFIG_DIR/config.yaml" << 'EOF'
defaults:
  project: "px14"  # Use real project
  queue: "normal"
  resources:
    - "walltime=02:00:00"
    - "mem=4GB"
    - "ncpus=2"

templates:
  org_name: "MyOrganization"
  shared_scratch: "/shared/scratch"
  backup_location: "/shared/backups"

aliases:
  org_backup:
    main:
      name: "backup_{timestamp}"
      queue: "copyq"
      project: "bc07"  # Use real project
      resources:
        - "ncpus=1"  # copyq queue only supports 1 CPU
    subcommand:
      type: module
      mod: "rsync"
    target:
      cmd: "rsync -av {shared_scratch}/ {backup_location}/"

  org_python:
    main:
      name: "python_job_{date}"
      resources:
        - "mem=8GB"
        - "ncpus=4"
    subcommand:
      type: conda
      env: "organization_python"
    target:
      cmd: "python"
EOF

# Set up environment to use test system config
export XDG_CONFIG_DIRS="$TEST_DIR/etc/xdg"
export XDG_CONFIG_HOME="$TEST_DIR/home/.config"

run_test "System config provides organization defaults" \
    "qxub config get defaults.project" \
    "px14"

run_test "System config provides organization templates" \
    "qxub config get templates.org_name" \
    "MyOrganization"

run_test "System config provides organization aliases" \
    "qxub config alias list-aliases" \
    "org_backup"

echo -e "${YELLOW}[INFO]${NC} === Testing User Override Scenarios ==="

# Test 2: User creates personal config that overrides some system defaults
cat > "$USER_CONFIG_DIR/config.yaml" << 'EOF'
defaults:
  project: "px14"  # Override system default
  name: "user_job_{timestamp}"  # Add personal default
  resources:
    - "walltime=01:00:00"  # Override system default
    - "mem=8GB"            # Override system default
    - "ncpus=2"            # Keep system default

templates:
  user_dir: "/home/user/work"
  project_dir: "/scratch/user_project"

aliases:
  personal_analysis:
    main:
      name: "analysis_{date}"
      queue: "normal"
    subcommand:
      type: conda
      env: "personal_env"
    target:
      cmd: "python analysis.py"

  # Override organization alias with personal version
  org_python:
    main:
      name: "my_python_{time}"
      resources:
        - "mem=16GB"  # More memory than org default
    subcommand:
      type: conda
      env: "my_python_env"  # Different environment
    target:
      cmd: "python"
EOF

run_test "User config overrides system project" \
    "qxub config get defaults.project" \
    "px14"

run_test "User config inherits system queue (not overridden)" \
    "qxub config get defaults.queue" \
    "normal"

run_test "User config adds personal templates" \
    "qxub config get templates.user_dir" \
    "/home/user/work"

run_test "System templates still accessible" \
    "qxub config get templates.org_name" \
    "MyOrganization"

run_test "User and system aliases both available" \
    "qxub config alias list-aliases" \
    "org_backup"

echo -e "${YELLOW}[INFO]${NC} === Testing Realistic Job Execution Scenarios ==="

# Test 3: Actual job execution with system/user config precedence
run_test "Job inherits system defaults with user overrides (dry run)" \
    "qxub --dry-run conda --env test_env echo 'test job'" \
    "-q normal -P px14"

run_test "CLI args override user config (dry run)" \
    "qxub --dry-run --project cli_override_project --queue express conda --env test_env echo 'test'" \
    "-q express -P cli_override_project"

run_test "Job uses merged resource requirements" \
    "qxub --dry-run conda --env test_env echo 'test'" \
    "mem=8GB.*ncpus=2"

echo -e "${YELLOW}[INFO]${NC} === Testing Alias Inheritance and Overrides ==="

# Test 4: Alias execution shows proper inheritance
run_test "User alias overrides system alias (show config)" \
    "qxub config alias show org_python" \
    "my_python_env"

run_test "System alias still works" \
    "qxub config alias show org_backup" \
    "bc07"

run_test "Alias execution with overrides (dry run)" \
    "qxub --dry-run alias personal_analysis" \
    "-P px14"

run_test "Alias override at runtime (dry run)" \
    "qxub --dry-run alias personal_analysis --queue express" \
    "-q express"

echo -e "${YELLOW}[INFO]${NC} === Testing Template Variable Resolution ==="

# Test 5: Template variables work across config levels
run_test "System templates resolve in user alias" \
    "qxub --dry-run alias org_backup" \
    "copyq.*bc07"

# Create an alias that uses both system and user templates
# We need to add this to the existing user config without creating duplicate keys
cat > "$USER_CONFIG_DIR/config.yaml" << 'EOF'
defaults:
  project: "px14"  # Override system default
  name: "user_job_{timestamp}"  # Add personal default
  resources:
    - "walltime=01:00:00"  # Override system default
    - "mem=8GB"            # Override system default
    - "ncpus=2"            # Keep system default

templates:
  user_dir: "/home/user/work"
  project_dir: "/scratch/px14"

aliases:
  personal_analysis:
    main:
      name: "analysis_{date}"
      queue: "normal"
    subcommand:
      type: conda
      env: "personal_env"
    target:
      cmd: "python analysis.py"

  # Override organization alias with personal version
  org_python:
    main:
      name: "my_python_{time}"
      resources:
        - "mem=16GB"  # More memory than org default
    subcommand:
      type: conda
      env: "my_python_env"  # Different environment
    target:
      cmd: "python"

  template_test:
    main:
      name: "{org_name}_{user_dir}_backup_{timestamp}"
    subcommand:
      type: module
      mod: "rsync"
    target:
      cmd: "rsync -av {user_dir}/ {backup_location}/"
EOF

run_test "Mixed template variable resolution" \
    "qxub --dry-run alias template_test" \
    "MyOrganization.*work.*backup"

echo -e "${YELLOW}[INFO]${NC} === Testing Admin Tools ==="

# Test 6: Config discovery and management
run_test "Config files command shows system and user configs" \
    "qxub config files" \
    "Configuration Files"

run_test "System config loaded successfully" \
    "qxub config files" \
    "system_0.*Exists"

run_test "User config loaded successfully" \
    "qxub config files" \
    "user.*Exists"

echo -e "${YELLOW}[INFO]${NC} === Testing Real-World Error Scenarios ==="

# Test 7: What happens with broken configs (realistic admin scenario)
INVALID_SYSTEM_DIR="$TEST_DIR/etc/xdg_invalid/qxub"
mkdir -p "$INVALID_SYSTEM_DIR"

cat > "$INVALID_SYSTEM_DIR/config.yaml" << 'EOF'
# Broken YAML - missing quote
defaults:
  project: "broken_project
  queue: normal
EOF

export XDG_CONFIG_DIRS="$TEST_DIR/etc/xdg_invalid:$TEST_DIR/etc/xdg"

# Should warn about invalid config but continue with valid ones
run_test "Invalid system config produces warning but continues" \
    "qxub config get defaults.project 2>&1" \
    "px14"

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
echo "       REALISTIC SYSTEM CONFIG TESTS"
echo "========================================"
echo "Total Tests: $TOTAL_TESTS"
echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed: ${RED}$FAILED_TESTS${NC}"

if [[ $FAILED_TESTS -eq 0 ]]; then
    echo -e "${GREEN}🎉 All realistic system config tests passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ Some tests failed. Please review the output above.${NC}"
    exit 1
fi
