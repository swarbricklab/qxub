# Enhanced --dry Output Implementation Plan

## Current Discovery

✅ **qxub already does submission-time variable expansion!**

Testing revealed:
- `qxub --dry --env base -- echo "Hello $USER"`
- Base64 decodes to: `echo Hello jr9959`
- Variables expand when user expects them to

## The Real Issue

Current `--dry` shows:
```
Dry run - would execute: qsub -v cmd_b64="ZWNobyBIZWxsbyBqcjk5NTksIHRlc3RpbmcgdmFyaWFibGVz",...
```

Users can't easily see what command will actually execute without manual base64 decoding.

## Proposed Enhancement

Enhance `--dry` output to show both the decoded command AND the qsub command:

```
🔧 Job command constructed
📝 Command to execute: echo Hello jr9959, testing variables
📋 Pre-command: export CUDA_VISIBLE_DEVICES=0
📋 Post-command: echo "Job completed"
Dry run - would execute: qsub -v cmd_b64="...",pre_cmd_b64="...",post_cmd_b64="..." ...
```

## Implementation

### Step 1: Add command decoding utility
```python
def decode_command_for_display(cmd_b64: str, pre_b64: str = None, post_b64: str = None):
    """Decode base64 commands for user-friendly display."""
    import base64

    cmd = base64.b64decode(cmd_b64).decode('utf-8') if cmd_b64 else None
    pre = base64.b64decode(pre_b64).decode('utf-8') if pre_b64 else None
    post = base64.b64decode(post_b64).decode('utf-8') if post_b64 else None

    return cmd, pre, post
```

### Step 2: Enhance dry run output in execution.py
```python
# In submit_and_monitor_job()
if ctx_obj["dry"]:
    # Show decoded commands for user verification
    cmd_str = " ".join(command)

    print(f"📝 Command to execute: {cmd_str}")

    if pre:
        print(f"📋 Pre-command: {pre}")
    if post:
        print(f"📋 Post-command: {post}")

    print(f"Dry run - would execute: {submission_command}")
    return
```

### Step 3: Consider verbose levels
```python
if ctx_obj["dry"]:
    # Always show decoded command
    cmd_str = " ".join(command)
    print(f"📝 Command to execute: {cmd_str}")

    # Show pre/post if present
    if pre:
        print(f"📋 Pre-command: {pre}")
    if post:
        print(f"📋 Post-command: {post}")

    # Show full qsub command only if verbose
    if ctx_obj.get("verbose", 0) > 0:
        print(f"🔧 Full qsub command: {submission_command}")
    else:
        print("Dry run - job would be submitted (use -v to see full qsub command)")
    return
```

## Benefits

1. **User-friendly**: Shows exactly what will execute
2. **Debugging**: Still shows full qsub command when needed
3. **No new flags**: Enhances existing `--dry` behavior
4. **Backward compatible**: Doesn't break existing scripts
5. **Covers all modes**: Works for conda, modules, singularity, default

## Testing

Add to edge cases test:
```bash
# Test enhanced dry run output
run_edge_case_test "Enhanced dry run shows decoded command" \
    "qxub --dry --env base -- echo 'Hello \$USER'" \
    0 "Command to execute: echo 'Hello \$USER'"

run_edge_case_test "Enhanced dry run shows pre/post commands" \
    "qxub --dry --env base --pre 'export VAR=test' --post 'echo done' -- echo 'main'" \
    0 "Pre-command: export VAR=test"
```

## Documentation Update

Update docs to show new dry run output:
```markdown
### Preview Commands with --dry

```bash
qxub --dry --env myenv -- python script.py --input $HOME/data.txt
# Output:
# 📝 Command to execute: python script.py --input /home/user/data.txt
# Dry run - job would be submitted
```

The `--dry` flag shows exactly what command will execute after variable expansion.
```

This enhancement solves the transparency issue without requiring new flags and leverages the already-correct variable expansion behavior.
