# Smart Quote Processing in qxub

## Overview

qxub's `--cmd` option supports **Smart Quote Processing** to make complex commands with nested quotes much more readable and maintainable. This feature automatically handles quote escaping when you use double quotes around your `--cmd` argument.

## The Problem

Complex commands often require nested quotes, leading to ugly escaping patterns:

```bash
# Traditional approach - hard to read and error-prone
qxub --env base --cmd 'find /data -exec sh -c '"'"'echo "Found: $1"'"'"' _ {} \;'
```

## The Solution

Smart Quote Processing allows you to write clean, readable commands:

```bash
# Smart quote approach - clean and readable
qxub --env base --cmd "find /data -exec sh -c \"echo \\\"Found: \$1\\\"\" _ {} \;"
```

## How It Works

### Detection

Smart Quote Processing activates when your `--cmd` argument:
1. Starts with a double quote (`"`)
2. Ends with a double quote (`"`)
3. Has content between the quotes

```bash
qxub --cmd "your command here"  # ← Smart processing
qxub --cmd 'your command here'  # ← Traditional processing
```

### Transformation Rules

| Input Pattern | Output | Purpose |
|--------------|--------|---------|
| `\"` | `"` | Literal double quote |
| `\$` | `$` | Literal dollar (no expansion) |
| `${var}` | `value` | Submission-time expansion |
| `${{var}}` | `${var}` | Execution-time preparation |
| `'text'` | `'text'` | Single quotes preserved |

## Examples

### Basic Quote Handling

```bash
# Input
qxub --cmd "echo \"Hello World\""

# Output
echo "Hello World"
```

### Literal Dollar Signs

```bash
# Input
qxub --cmd "echo \"Cost: \$100\""

# Output
echo "Cost: $100"
```

### Mixed Variables

```bash
# Input
qxub --cmd "echo \"User: ${USER}, Job: ${{PBS_JOBID}}\""

# Output (assuming USER=jr9959)
echo "User: jr9959, Job: ${PBS_JOBID}"
```

### Complex Find Commands

```bash
# Input
qxub --cmd "find /data -name \"*.log\" -exec sh -c \"echo \\\"Processing: \$1\\\"\" _ {} \;"

# Output
find /data -name "*.log" -exec sh -c "echo \"Processing: $1\"" _ {} \;
```

### AWK Field References

```bash
# Input
qxub --cmd "awk \"{print \\\$1, \\\$2}\" data.txt"

# Output
awk "{print $1, $2}" data.txt
```

### JSON Output Generation

```bash
# Input
qxub --cmd "echo \"{\\\"file\\\": \\\"\$1\\\", \\\"user\\\": \\\"${USER}\\\"}\""

# Output (assuming USER=jr9959)
echo "{\"file\": \"$1\", \"user\": \"jr9959\"}"
```

## Comparison: Before vs After

### GNU find with -exec

**Before (ugly nested quotes):**
```bash
qxub --cmd 'find /path -exec sh -c '"'"'echo "Processing $1 for user $USER"'"'"' _ {} \;'
```

**After (clean smart quotes):**
```bash
qxub --cmd "find /path -exec sh -c \"echo \\\"Processing \$1 for user ${USER}\\\"\" _ {} \;"
```

### Python Commands with f-strings

**Before (quote mangling issues):**
```bash
qxub --cmd 'python -c "print(f\"User: {os.environ.get(\"USER\")}, Args: {len(sys.argv)}\")"'
# Often breaks due to shell processing
```

**After (protected from shell):**
```bash
qxub --cmd "python -c \"print(f\\\"User: {os.environ.get(\\\"USER\\\")}, Args: {len(sys.argv)}\\\")\""
```

## Backward Compatibility

Smart Quote Processing is **completely backward compatible**:

- Single-quoted commands work exactly as before
- Existing scripts and aliases continue to work unchanged
- No breaking changes to existing functionality

```bash
# These continue to work exactly as before
qxub --cmd 'echo "Hello ${USER}"'
qxub --cmd 'python script.py --arg value'
```

## Best Practices

### When to Use Smart Quotes

✅ **Use Smart Quotes For:**
- Commands with nested quotes (find -exec, sh -c, etc.)
- JSON or structured output generation
- Complex shell commands with multiple escaping levels
- AWK/sed commands with field references
- Any command where readability matters

### When to Use Traditional Quotes

✅ **Use Traditional Quotes For:**
- Simple commands without nested quotes
- Existing working commands (no need to change)
- Shell patterns you want expanded before qxub sees them

### Escaping Guidelines

**In Smart Quote mode:**
```bash
# For literal quotes inside the command
\"     # Becomes "

# For literal dollars (no variable expansion)
\$     # Becomes $

# For qxub variable expansion (still works normally)
${var} # Expands at submission time
${{var}} # Expands at execution time
```

## Troubleshooting

### Common Mistakes

**Wrong: Forgetting to escape internal quotes**
```bash
qxub --cmd "echo "Hello World""  # Syntax error
```

**Right: Properly escaped quotes**
```bash
qxub --cmd "echo \"Hello World\""  # Works correctly
```

**Wrong: Using backslashes unnecessarily**
```bash
qxub --cmd "echo \\"Hello World\\""  # Double-escaped
```

**Right: Minimal necessary escaping**
```bash
qxub --cmd "echo \"Hello World\""  # Clean and correct
```

### Debugging

Use `--dry` to see how your command is processed:

```bash
qxub --env base --dry --cmd "your command here"
```

This shows exactly what command qxub will execute, helping you verify the quote processing worked as expected.

## Migration Guide

### Updating Existing Commands

1. **Identify complex commands** with nested quotes
2. **Wrap in double quotes** instead of single quotes
3. **Replace internal double quotes** with `\"`
4. **Escape literal dollars** with `\$` if needed
5. **Test with --dry** to verify processing

### Example Migration

```bash
# Before: Ugly but functional
qxub --cmd 'find /data -exec sh -c '"'"'echo "Found: $1"'"'"' _ {} \;'

# After: Clean and readable
qxub --cmd "find /data -exec sh -c \"echo \\\"Found: \$1\\\"\" _ {} \;"
```

Smart Quote Processing makes qxub's `--cmd` option much more user-friendly for complex shell operations while maintaining full backward compatibility and the powerful variable substitution system.
