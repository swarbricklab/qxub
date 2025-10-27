# Refactoring Checklist - Preventing Zombie Code

## Problem

During refactoring (especially file moves/reorganizations), old files are often left behind as "zombie code" - unused but not deleted. This happens because:
- Import updates make tests pass
- No verification that old files are unused
- Fear of breaking something
- Documentation doesn't explicitly mandate deletion

## Mandatory Checklist for Code Reorganization

### Before Starting

- [ ] **Document the current state**: List all files that will be moved/renamed
- [ ] **Plan the new structure**: Write down exact new locations
- [ ] **Identify all imports**: Search for all current import statements

### During Refactoring

- [ ] **Use `git mv`** when moving files (not copy + modify)
  ```bash
  git mv old/path/file.py new/path/file.py
  ```

- [ ] **Update ALL imports** in one commit
  ```bash
  # Find all imports of the old location
  git grep "from.*old_module import"
  git grep "import.*old_module"
  ```

- [ ] **Verify old imports are gone**:
  ```bash
  # These should return NOTHING after refactoring:
  git grep "from qxub.config_manager import"
  git grep "import qxub.config_manager"
  ```

### After Refactoring

- [ ] **Check for zombie files**:
  ```bash
  # Look for duplicate modules/packages
  find qxub -name "config*.py" -type f
  find qxub -name "*_manager.py" -type f

  # Check if old files exist alongside new package dirs
  ls -la qxub/config_manager.py 2>/dev/null && echo "⚠️  ZOMBIE FILE!"
  ```

- [ ] **Verify old files are unused**:
  ```bash
  # For each suspected zombie file:
  git grep "from.*zombie_module import"
  git grep "import.*zombie_module"

  # Should return 0 results!
  ```

- [ ] **Delete verified zombie files**:
  ```bash
  git rm qxub/old_file.py
  ```

- [ ] **Update documentation**:
  - Mark old import paths as ❌ DEPRECATED (do not use)
  - Add ✅ NEW import paths
  - Document migration path for external users

- [ ] **Run full test suite**:
  ```bash
  pytest tests/
  python -m qxub.config  # Verify imports work
  ```

- [ ] **Verify package imports**:
  ```python
  # Test that new imports work
  python3 -c "from qxub.config import config_manager; print('✅ New import works')"

  # Test that old imports fail (if intentionally removed)
  python3 -c "from qxub.config_manager import ConfigManager" 2>&1 | grep "ModuleNotFoundError" && echo "✅ Old import correctly fails"
  ```

## Specific Anti-Patterns to Avoid

### ❌ DON'T: Copy and modify
```bash
cp qxub/config_manager.py qxub/config/manager.py
# Edit qxub/config/manager.py
# Update imports
# Commit
# ⚠️ Now you have two files!
```

### ✅ DO: Move and track
```bash
git mv qxub/config_manager.py qxub/config/manager.py
# Edit if needed
# Update imports in same commit
# Commit with clear message: "Move config_manager to config package"
```

### ❌ DON'T: Leave "legacy" files "just in case"
- Either delete them or explicitly maintain them
- "Legacy" without a deprecation timeline = technical debt

### ✅ DO: Be explicit about legacy code
- Add deprecation warnings to the code
- Set a removal date
- Update documentation with migration guide
- Delete after deprecation period

## Recovery: Finding Existing Zombies

Run this audit regularly:

```bash
#!/bin/bash
# audit-zombie-files.sh

echo "=== Checking for zombie files ==="

# Check for duplicate module patterns
for base in config platform execution history resources; do
    if [ -f "qxub/${base}.py" ] && [ -d "qxub/${base}/" ]; then
        echo "⚠️  POTENTIAL ZOMBIE: qxub/${base}.py (package also exists at qxub/${base}/)"
        echo "   Check if it's imported:"
        git grep -q "from.*${base} import" && echo "   ✅ Used" || echo "   ❌ UNUSED - DELETE IT"
    fi

    if [ -f "qxub/${base}_manager.py" ] && [ -d "qxub/${base}/" ]; then
        echo "⚠️  POTENTIAL ZOMBIE: qxub/${base}_manager.py (package exists at qxub/${base}/)"
        echo "   Check if it's imported:"
        git grep -q "from.*${base}_manager import" && echo "   ✅ Used" || echo "   ❌ UNUSED - DELETE IT"
    fi
done

echo ""
echo "=== Checking for standalone vs package conflicts ==="
for standalone in qxub/*_*.py; do
    if [ -f "$standalone" ]; then
        basename=$(basename "$standalone" .py)
        pkg_name=$(echo "$basename" | cut -d_ -f1)
        if [ -d "qxub/${pkg_name}/" ]; then
            echo "⚠️  CHECK: $standalone (package qxub/${pkg_name}/ exists)"
            module_name=$(echo "$basename" | sed 's/_/./g')
            git grep -q "from.*${basename} import" && echo "   ✅ Used" || echo "   ❌ Likely zombie"
        fi
    fi
done
```

## Post-Refactoring Commit Message Template

```
refactor: Move [module] to [new_location] package

BREAKING CHANGE: Old import path removed

Before: from qxub.config_manager import ConfigManager
After:  from qxub.config import config_manager

Changes:
- Move qxub/config_manager.py → qxub/config/manager.py
- Delete qxub/config_manager.py (zombie file)
- Update all imports to new package structure
- Add deprecation notice to __init__.py

Migration: Update imports as shown above
```

## Prevention: Code Review Checklist

Reviewers should verify:
- [ ] If files were moved, old files were deleted
- [ ] No duplicate modules (e.g., both `module.py` and `module/` exist)
- [ ] Import paths use new locations consistently
- [ ] Old import paths don't work (or have deprecation warnings)
- [ ] Documentation updated with new import paths

## Why This Matters

**Impact of zombie files:**
- 🐛 Confusion: "Which one is the real module?"
- 📈 Technical debt: Dead code that needs maintenance
- 🔍 Hard to navigate: Multiple files with similar names
- ⚠️ Merge conflicts: Changes go to wrong file
- 💾 Bloat: Unnecessary files in repository

**In this project:**
- Found 1010 lines of zombie code (config_manager.py + config_handler.py)
- Zero imports of these files
- Caused confusion during development
- Could have led to changes in wrong file

## Remember

**When you move code, DELETE the old location. No exceptions.**

If you're not sure if it's safe:
1. Check imports: `git grep "from.*old_module"`
2. Run tests: `pytest tests/`
3. If both pass with old file gone → DELETE IT

**"Legacy" is not a permanent status. It's a TODO with a deadline.**
