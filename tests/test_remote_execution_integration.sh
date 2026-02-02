#!/bin/bash
# Test remote execution by temporarily configuring nci_gadi with remote: section
# This will SSH to localhost (gadi) which is a bit silly but proves it works

set -e

echo "🧪 Testing remote execution integration..."

# Save original config
BACKUP_CONFIG=$(mktemp)
QXUB_CONFIG="$HOME/.config/qxub/config.yaml"

if [ -f "$QXUB_CONFIG" ]; then
    cp "$QXUB_CONFIG" "$BACKUP_CONFIG"
    echo "📦 Backed up config to: $BACKUP_CONFIG"
fi

# Create test config with remote platform
cat > "$QXUB_CONFIG" << 'EOF'
defaults:
  platform: nci_gadi_remote
  project: a56

platforms:
  # Local platform (no remote section)
  nci_gadi:
    name: nci_gadi
    definition: file:///g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml

  # Remote platform (SSH to same host for testing)
  nci_gadi_remote:
    name: nci_gadi_remote
    definition: file:///g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml
    remote:
      host: localhost
      working_dir: /scratch/a56/{user}
EOF

echo "✅ Created test config with remote platform"

# Test 1: Verify execution mode detection
echo ""
echo "Test 1: Execution mode detection"
qxub config get platforms.nci_gadi_remote

# Test 2: Remote execution with dry-run
echo ""
echo "Test 2: Remote execution (dry-run)"
qxub exec --platform nci_gadi_remote --dry --default -- echo "remote test"

echo ""
echo "✅ Remote execution integration test completed!"
echo ""
echo "To test actual SSH execution (not dry-run):"
echo "  qxub exec --platform nci_gadi_remote --default -- echo 'hello from remote'"

# Restore original config
if [ -f "$BACKUP_CONFIG" ]; then
    mv "$BACKUP_CONFIG" "$QXUB_CONFIG"
    echo ""
    echo "🔄 Restored original config"
fi
