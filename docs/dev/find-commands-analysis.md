# GNU find Commands with qxub: Edge Cases and Best Practices

## Overview

GNU `find` commands are notorious for complex syntax involving special characters, quotes, and shell metacharacters. This document explores how different `find` patterns behave with qxub's traditional `--` syntax versus the new `--cmd` option.

## Key Findings

### Quote Preservation

**Traditional `--` syntax:**
```bash
qxub --env base -- find /path -name "*.txt" -exec echo "Found: {}" \;
# Result: find /path -name *.txt -exec echo Found: {} ;
# ❌ Quotes removed by shell processing
```

**`--cmd` syntax:**
```bash
qxub --env base --cmd 'find /path -name "*.txt" -exec echo "Found: {}" \;'
# Result: find /path -name "*.txt" -exec echo "Found: {}" \;
# ✅ Quotes preserved exactly
```

### Variable Expansion Timing

**Submission-time variables** (`${var}` - expanded when you submit):
```bash
qxub --env base --cmd 'find ${HOME} -maxdepth 1 -name ".*" -exec echo "User: ${USER}" \;'
# ${HOME} and ${USER} expanded immediately to actual values
```

**Execution-time variables** (`${{var}}` - expanded when job runs):
```bash
qxub --env base --cmd 'find /data -name "*.log" -exec echo "Job ${{PBS_JOBID}} found: {}" \;'
# ${{PBS_JOBID}} becomes ${PBS_JOBID} for expansion during job execution
```

**Mixed variables:**
```bash
qxub --env base --cmd 'find ${HOME}/logs -name "*.log" -exec echo "User ${USER} job ${{PBS_JOBID}}: {}" \;'
# ${HOME} and ${USER} expanded at submission, ${{PBS_JOBID}} at execution
```

## Critical Edge Cases

### 1. Semicolon Terminators (`\;`)

Both syntaxes handle semicolon terminators, but quote preservation differs:

```bash
# Traditional - basic structure preserved
qxub --env base -- find /path -name "*.txt" -exec cmd {} \;

# --cmd - full quote control
qxub --env base --cmd 'find /path -name "*.txt" -exec cmd {} \;'
```

### 2. Complex Shell Commands in `-exec`

**Problem case with traditional syntax:**
```bash
qxub --env base -- find /path -exec sh -c 'echo "File: $1"' _ {} \;
# Shell processing mangles the quotes and variables
```

**Solution with --cmd:**
```bash
qxub --env base --cmd 'find /path -exec sh -c '"'"'echo "File: $1"'"'"' _ {} \;'
# Full control over quote structure
```

### 3. AWK Field References

AWK field references (`$1`, `$2`, etc.) are preserved in both syntaxes:

```bash
qxub --env base --cmd 'find /path -name "*.txt" -exec awk '"'"'{print $1, $2}'"'"' {} \;'
# $1 and $2 preserved as literal AWK field references
```

### 4. Glob Patterns

**Traditional syntax (shell expands):**
```bash
qxub --env base -- find /path -name '[fn]*'
# Shell may expand [fn]* before find sees it
```

**--cmd syntax (literal preserved):**
```bash
qxub --env base --cmd 'find /path -name '"'"'[fn]*'"'"
# Pattern preserved literally for find to process
```

### 5. Regex Patterns

Similar issues with regex patterns:

```bash
# Traditional - quotes lost
qxub --env base -- find /path -regex '.*\.\(txt\|py\)'
# Result: find /path -regex .*\.\(txt\|py\)

# --cmd - quotes preserved
qxub --env base --cmd 'find /path -regex '"'"'.*\.\(txt\|py\)'"'"
# Result: find /path -regex '.*\.\(txt\|py\)'
```

### 6. JSON-like Output

Complex structured output with nested quotes:

```bash
qxub --env base --cmd 'find /path -name "*.txt" -exec sh -c '"'"'echo "{\"file\": \"$1\", \"user\": \"${USER}\", \"cost\": \"\$50\"}"'"'"' _ {} \;'
# Produces: {"file": "/path/file.txt", "user": "jr9959", "cost": "$50"}
```

### 7. Pipes and Command Chains

```bash
qxub --env base --cmd 'find /path -name "*.txt" -print0 | xargs -0 wc -l'
# Pipe preserved and processed correctly
```

## Best Practices

### When to Use Traditional `--` Syntax

✅ **Good for:**
- Simple find commands without complex quoting
- Basic `-name` patterns with simple globs
- When you want shell expansion of variables

```bash
qxub --env base -- find /simple/path -name "*.txt" -print
qxub --env base -- find "$HOME" -maxdepth 1 -name ".*" -type f
```

### When to Use `--cmd` Syntax

✅ **Required for:**
- Complex `-exec` commands with nested quotes
- Regex patterns that need literal preservation
- Commands mixing submission and execution time variables
- JSON or structured output generation
- Glob patterns that shouldn't be shell-expanded
- Any command where quote structure is critical

```bash
qxub --env base --cmd 'find ${HOME}/logs -name "*.log" -exec echo "User ${USER} job ${{PBS_JOBID}}: {}" \;'
```

## Variable Substitution Rules

| Pattern | Expansion Time | Example | Result |
|---------|---------------|---------|---------|
| `${var}` | Submission | `${USER}` | `jr9959` |
| `${{var}}` | Execution | `${{PBS_JOBID}}` | `${PBS_JOBID}` (for job) |
| `$var` | Literal | `$1` (AWK) | `$1` (unchanged) |
| `\$var` | Literal | `\$100` | `$100` (escaped) |

## Testing Results Summary

From comprehensive testing of 22 different find command patterns:

- **Quote preservation**: --cmd syntax superior for complex quotes
- **Variable timing**: Both support qxub's `${var}` vs `${{var}}` system
- **Special characters**: --cmd provides better control
- **Shell interference**: Traditional syntax subject to shell processing
- **Compatibility**: Both work, --cmd offers more predictability

## Recommendations

1. **Start with traditional `--`** for simple find operations
2. **Switch to `--cmd`** when you encounter quote mangling issues
3. **Always use `--cmd`** for production scripts with complex find operations
4. **Use `${var}` for paths known at submission time** (like `${HOME}`)
5. **Use `${{var}` for job-specific values** (like `${{PBS_JOBID}}`)
6. **Test complex commands** with `--dry` flag first

## Example Migration

**Before (problematic):**
```bash
qxub --env base -- find "$HOME/data" -name "*.log" -exec sh -c 'echo "Processing $1 for user $USER"' _ {} \;
```

**After (reliable):**
```bash
qxub --env base --cmd 'find ${HOME}/data -name "*.log" -exec sh -c '"'"'echo "Processing $1 for user ${USER}"'"'"' _ {} \;'
```

The `--cmd` syntax provides the control needed for complex GNU find operations while maintaining compatibility with qxub's variable substitution system.
