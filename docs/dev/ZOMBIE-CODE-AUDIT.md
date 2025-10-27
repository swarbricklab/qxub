# Zombie Code Audit - October 27, 2025

## Critical Issue: Divergent Duplicate Code

During v3.3.0 development, I discovered a serious code maintenance issue: **incomplete package migration has created divergent duplicate code**.

## The Problem

When the codebase was migrated from flat structure to packages (v3.0), files were moved but not deleted, creating:

### ✅ DELETED (Confirmed Zombies)
- `qxub/config_manager.py` (765 lines) - Moved to `qxub/config/manager.py`
- `qxub/config_handler.py` (245 lines) - Moved to `qxub/config/handler.py`
- **Total: 1,010 lines of dead code removed**

### ⚠️ DIVERGENT DUPLICATES (Active Maintenance Issue)
- `qxub/execution.py` (441 lines) vs `qxub/execution/core.py` (427 lines)
  - **Status**: Different! Files have diverged
  - **Impact**: Changes may go to wrong file
  - **Risk**: HIGH - Core execution logic

- `qxub/platform.py` (798 lines) vs `qxub/platform/core.py` (796 lines)
  - **Status**: Different! Files have diverged
  - **Impact**: Platform logic inconsistency
  - **Risk**: HIGH - Platform selection logic

### 🤔 INTENTIONAL DUPLICATES? (Need Verification)
- `qxub/execution_context.py` - Package `qxub/execution/` exists
  - **Status**: Need to check if it's a wrapper or duplicate logic
  - **Action Required**: Audit and decide

## Why This Is Dangerous

1. **Maintenance Confusion**: "Which file do I edit?"
2. **Divergent Changes**: Files have different code (already happening!)
3. **Bug Risk**: Fix applied to wrong file
4. **Merge Conflicts**: Changes to both files
5. **Code Review Difficulty**: Reviewers don't know which is "real"

## Evidence of Divergence

```bash
# execution files differ:
$ diff -q qxub/execution.py qxub/execution/core.py
Files qxub/execution.py and qxub/execution/core.py differ

# platform files differ:
$ diff -q qxub/platform.py qxub/platform/core.py
Files qxub/platform.py and qxub/platform/core.py differ

# Changes still being made to BOTH:
$ git log --oneline qxub/execution.py | head -3
b75f647 v3.2.0: Add logging tools (Oct 26, 2025)
05fdddc v3.2.0: Add logging tools (Oct 26, 2025)
2dc8507 feat: add job file logging (Oct 26, 2025)
```

## Root Cause

The v3.0 package migration (commit 411be49 "Complete Phase 2: Configuration package migration") moved files but didn't delete originals. This was documented as "compatibility wrappers" but:

1. **They're not wrappers** - They contain actual implementation
2. **No re-exports** - They don't import from packages
3. **No deprecation warnings** - Silent divergence
4. **Still being edited** - Active development on both copies

## Proper Migration Pattern

### ❌ What Happened (WRONG)
```bash
# Move file
cp qxub/config_manager.py qxub/config/manager.py
# Update imports
sed -i 's/from .config_manager/from .config/' qxub/**/*.py
# Commit
git commit -m "Create config package"
# ⚠️ Old file still exists with implementation
```

### ✅ What Should Happen (CORRECT)

**Option A: Complete Migration (RECOMMENDED)**
```bash
# Move file
git mv qxub/config_manager.py qxub/config/manager.py
# Update imports
sed -i 's/from .config_manager/from .config/' qxub/**/*.py
# Commit
git commit -m "Move config_manager to config package"
# Old file is GONE ✅
```

**Option B: Thin Wrapper (If backwards compatibility needed)**
```python
# qxub/execution.py becomes a thin wrapper:
"""
Legacy import compatibility wrapper.

DEPRECATED: Import from qxub.execution package instead.
This file will be removed in qxub v4.0.
"""
import warnings

warnings.warn(
    "qxub.execution module is deprecated. Use 'from qxub.execution import ...' instead.",
    DeprecationWarning,
    stacklevel=2
)

# Re-export from package
from .execution.core import *  # noqa
from .execution.context import *  # noqa
from .execution.executors import *  # noqa
from .execution.mode import *  # noqa
```

## Recovery Plan

### Phase 1: Immediate (Before v3.3.0)

1. **Audit execution files**:
   - Compare `qxub/execution.py` vs `qxub/execution/core.py`
   - Determine which has the latest changes
   - Merge any unique fixes

2. **Audit platform files**:
   - Compare `qxub/platform.py` vs `qxub/platform/core.py`
   - Ensure platform registry changes are in right place
   - Merge any divergent fixes

3. **Delete or convert**:
   - Either delete standalone files completely
   - Or convert to thin deprecation wrappers
   - **Do NOT leave as-is**

4. **Update documentation**:
   - Mark old imports as deprecated
   - Add migration guide
   - Update package structure docs

### Phase 2: Verification

1. **Verify imports**:
   ```bash
   # Should return NOTHING:
   git grep "import qxub.execution$"  # direct file import
   git grep "import qxub.platform$"   # direct file import
   ```

2. **Run tests**:
   ```bash
   pytest tests/
   ./tests/test_conda_dry.sh
   ```

3. **Check for other zombies**:
   ```bash
   bash /tmp/audit-zombies.sh  # Run audit script
   ```

### Phase 3: Prevention

1. **Add pre-commit hook** to detect duplicate files
2. **Update refactoring-checklist.md** (done ✅)
3. **Code review checklist** includes zombie audit
4. **CI check** for file/package name conflicts

## Lessons Learned

### Why This Kept Happening

1. **Pattern Blindness**: I see a file and use it without checking for package
2. **No Verification**: Don't check if imports resolve to file or package
3. **Fear of Breaking**: Rather leave old file "just in case"
4. **Documentation Ambiguity**: "Legacy" doesn't mean "delete immediately"
5. **No Automated Checks**: Nothing catches duplicate file/package names

### How to Stop It

**When you see a module name that could be either file or package:**

```python
# ⚠️ AMBIGUOUS - Could be EITHER:
from qxub.execution import something
# Could be: qxub/execution.py
# Could be: qxub/execution/__init__.py

# ✅ ALWAYS VERIFY:
import qxub.execution
print(qxub.execution.__file__)  # Shows which file is actually loaded
```

**Before using a module, run:**
```bash
# Check for file/package conflict
if [ -f "qxub/execution.py" ] && [ -d "qxub/execution/" ]; then
    echo "⚠️  CONFLICT: Both file and package exist!"
    echo "   Run: python3 -c 'import qxub.execution; print(qxub.execution.__file__)'"
fi
```

**When migrating code:**
- Use `git mv` (tracks moves)
- Delete old files in SAME commit as import updates
- Add test to verify old imports fail
- Update documentation immediately

## Action Items

- [ ] Audit and merge `execution.py` vs `execution/core.py`
- [ ] Audit and merge `platform.py` vs `platform/core.py`
- [ ] Delete or wrap standalone files
- [ ] Add pre-commit hook for duplicate detection
- [ ] Update CI to check for file/package conflicts
- [ ] Add deprecation warnings if keeping wrappers
- [ ] Set removal date for deprecated files

## Commit Message Template

```
refactor: Remove divergent duplicate execution/platform files

BREAKING CHANGE: Standalone qxub/execution.py and qxub/platform.py removed

These files were left behind during v3.0 package migration and have
since diverged from the package implementations, creating maintenance
confusion and bug risk.

Changes:
- Merge unique changes from standalone files into packages
- Delete qxub/execution.py (use qxub.execution package)
- Delete qxub/platform.py (use qxub.platform package)
- Verify all imports use package (not standalone file)

Migration: Already using packages via existing imports
No external API changes required.
```

---

**Remember**: When you find duplicate code, immediately:
1. ✅ Determine which is canonical
2. ✅ Merge any unique changes
3. ✅ Delete the duplicate
4. ✅ Verify tests pass
5. ✅ Commit with clear message

**"If it's not deleted, it's not migrated."**
