---
description: 'Chat mode for safe code refactoring and reorganization'
tools:
  # Core tools for refactoring
  - search
  - codebase
  - fileSearch
  - readFile
  - editFiles
  - usages
  - problems
  - changes
  - terminal
model: GPT-5
---
# Refactoring Chat Mode

This chat mode enforces strict practices to prevent zombie code, divergent duplicates, and incomplete migrations during code reorganization.

## Core Principle

**"If it's not deleted, it's not migrated."**

When moving code, the old location MUST be deleted in the same commit. No exceptions.

## Mandatory Checklist

Always follow `docs/dev/refactoring-checklist.md`. Key steps:

### Before Starting

1. **Document current state**: List all files being moved/renamed
2. **Plan new structure**: Write down exact new locations
3. **Identify all imports**: Search for current import statements

### During Refactoring

1. **Use `git mv`** when moving files:
   ```bash
   git mv old/path/file.py new/path/file.py
   ```
   Never: `cp old/path/file.py new/path/file.py`

2. **Update ALL imports** in one commit:
   ```bash
   # Find all imports to update
   git grep "from.*old_module import"
   git grep "import.*old_module"
   ```

3. **Verify old imports are gone**:
   ```bash
   # These MUST return nothing after refactoring:
   git grep "from qxub.old_module import"
   ```

### After Refactoring

1. **Check for zombie files**:
   ```bash
   # Look for file/package conflicts
   find qxub -name "module_name.py" -type f
   [ -f "qxub/module.py" ] && [ -d "qxub/module/" ] && echo "⚠️  CONFLICT!"
   ```

2. **Verify old files are unused**:
   ```bash
   # For each suspected zombie:
   git grep "from.*zombie_module import"  # Should return 0 results
   ```

3. **Delete verified zombies**:
   ```bash
   git rm qxub/old_file.py
   ```

4. **Run full test suite**:
   ```bash
   pytest tests/
   # Project-specific tests
   ```

## Anti-Patterns to Prevent

### ❌ NEVER: Copy and modify
```bash
cp qxub/config_manager.py qxub/config/manager.py
# Edit qxub/config/manager.py
# Update imports
# Commit
# ⚠️ Now you have TWO files! (zombie created)
```

### ✅ ALWAYS: Move and track
```bash
git mv qxub/config_manager.py qxub/config/manager.py
# Edit if needed
# Update imports in same commit
git commit -m "Move config_manager to config package"
```

### ❌ NEVER: Leave "legacy" files "just in case"
- Either delete them or add explicit deprecation
- "Legacy" without timeline = technical debt

### ✅ ALWAYS: Be explicit about deprecated code
- Add deprecation warnings in the code
- Set a removal date
- Update documentation with migration guide
- Delete after deprecation period

## Conscious Git Tracking

**ALWAYS stage files individually or by specific subdirectory**, never use blanket `git add .`

### ✅ Good Git Practices
```bash
# Stage specific files
git add qxub/config/manager.py
git add qxub/execution/mode.py

# Stage by logical subdirectory
git add qxub/config/
git add qxub/execution/

# Stage related documentation
git add docs/dev/refactoring-notes.md

# Verify what's staged before committing
git status
git diff --staged
```

### ❌ Bad Git Practices (Causes Zombie Files)
```bash
# DON'T: Blanket add everything
git add .              # ❌ Adds temporary scripts, test files, zombies!
git add -A             # ❌ Same problem
git add qxub/*.py      # ❌ Might add unrelated changes

# DON'T: Add without verifying
git add file.py && git commit  # ❌ No review of what's staged
```

### Conscious Tracking Workflow
```bash
# 1. Check what's changed
git status

# 2. Review each file individually
git diff qxub/config/manager.py

# 3. Stage intentionally
git add qxub/config/manager.py

# 4. Verify staged changes
git diff --staged

# 5. Check for unintended files
git status | grep "modified:" | grep -v "to be committed"

# 6. Commit with clear message
git commit -m "refactor: Move config_manager to package"
```

## File/Package Conflict Detection

Before committing, ALWAYS run this audit:

```bash
#!/bin/bash
echo "=== Checking for file/package conflicts ==="

for base in config platform execution history resources; do
    if [ -f "qxub/${base}.py" ] && [ -d "qxub/${base}/" ]; then
        echo "⚠️  CONFLICT: Both qxub/${base}.py AND qxub/${base}/ exist"
        git grep -q "from.*\.${base} import" && echo "   Package is used" || echo "   ❌ File is zombie"
    fi

    if [ -f "qxub/${base}_manager.py" ] && [ -d "qxub/${base}/" ]; then
        echo "⚠️  POTENTIAL ZOMBIE: qxub/${base}_manager.py"
        git grep -q "from.*${base}_manager import" && echo "   ✅ Used" || echo "   ❌ ZOMBIE - DELETE"
    fi
done

# Check if files have diverged
for base in execution platform; do
    if [ -f "qxub/${base}.py" ] && [ -f "qxub/${base}/core.py" ]; then
        echo ""
        echo "Checking divergence: qxub/${base}.py vs qxub/${base}/core.py"
        if diff -q "qxub/${base}.py" "qxub/${base}/core.py" > /dev/null; then
            echo "   ✅ Files are identical"
        else
            echo "   ⚠️  FILES HAVE DIVERGED - High risk!"
            echo "   Action required: Merge and delete duplicate"
        fi
    fi
done
```

## Import Resolution Verification

When unsure if imports resolve to file or package:

```python
# Check which file is actually loaded
import qxub.execution
print(f"Loaded from: {qxub.execution.__file__}")

# Expected: .../qxub/execution/__init__.py (package)
# Zombie if: .../qxub/execution.py (standalone file)
```

Or use terminal:
```bash
python3 -c "import qxub.execution; print(qxub.execution.__file__)"
```

## Refactoring Commit Message Template

```
refactor: [Brief description of what was moved/renamed]

[BREAKING CHANGE: Old import paths removed (if applicable)]

Before: from qxub.old_module import Class
After:  from qxub.new_location import Class

Changes:
- Move qxub/old.py → qxub/new/location.py (use git mv)
- Delete qxub/old.py (prevent zombie)
- Update all imports to new location
- Add deprecation notice if needed

Files changed:
- Moved: qxub/old.py → qxub/new/location.py
- Updated: [list all files with import changes]
- Deleted: [list any zombie files removed]

Migration: [Instructions for external users if breaking]
Testing: [Confirm tests pass]
```

## When to Use This Mode

Activate this chat mode when:

- Moving files to new locations
- Reorganizing package structure
- Renaming modules or packages
- Creating new package directories from flat files
- Consolidating duplicate code
- Cleaning up technical debt
- Any operation involving `git mv`

## Red Flags to Watch For

🚩 **File and package with same name exist**
- `qxub/module.py` AND `qxub/module/` both exist
- Action: Verify which is canonical, delete the other

🚩 **No git mv in the refactoring**
- Copy operations instead of moves
- Action: Use `git mv` to preserve history

🚩 **Blanket `git add .` used**
- Risk of adding temporary files, zombies, unrelated changes
- Action: Stage files individually

🚩 **"Legacy" or "compatibility" in comments**
- Without deprecation timeline
- Action: Add explicit deprecation or delete immediately

🚩 **Different content in similarly named files**
- `module.py` differs from `module/core.py`
- Action: Files have diverged! Merge and delete duplicate

🚩 **Tests pass but old file still exists**
- Old import path not tested
- Action: Add negative test to ensure old import fails

## Success Criteria

Refactoring is complete ONLY when:

- ✅ Old files deleted (visible in `git status`)
- ✅ Old import paths return zero results in `git grep`
- ✅ No file/package name conflicts exist
- ✅ All tests pass
- ✅ Documentation updated with new paths
- ✅ Commit message clearly describes what was moved
- ✅ No zombie files detected by audit script

## Recovery from Zombie Code

If you discover zombie files:

1. **Identify the canonical version**: Which is actively imported?
   ```bash
   git grep "from.*module import" | wc -l
   ```

2. **Compare for divergence**:
   ```bash
   diff -u old_file.py new_location/file.py
   ```

3. **Merge any unique changes**: Don't lose bug fixes

4. **Delete the zombie**:
   ```bash
   git rm old_file.py
   ```

5. **Document in commit**: Explain which was canonical and why

6. **Add to zombie audit doc**: Learn from the mistake

## Remember

- **Move = `git mv` + update imports + delete old file**
- **Verify = grep for old imports returns nothing**
- **Commit = clear message + conscious staging**
- **Test = full suite passes**

**When in doubt, ask: "Does the old location still exist?" If yes, deletion is incomplete.**
