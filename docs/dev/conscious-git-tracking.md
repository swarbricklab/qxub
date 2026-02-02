# Conscious Git Tracking

**Never use `git add .` or `git add -A`**

These create zombie files, commit temp scripts, and bundle unrelated changes.

## The Rule

**Stage files individually or by specific subdirectory.**

```bash
# ✅ Good
git add qxub/config/manager.py
git add qxub/execution/
git add docs/dev/v3.3.0-*.md

# ❌ Bad
git add .
git add -A
```

## Workflow

```bash
git status                       # See what changed
git add qxub/config/manager.py   # Stage deliberately
git diff --staged --name-only    # Verify
git commit -m "feat: message"    # Commit
```

## Before Every Commit

```bash
# Verify what's staged
git diff --staged --name-only

# Check for mistakes (temp files, zombies)
git diff --staged --name-only | grep -E "(tmp|temp|\.bak|debug_)"
```
