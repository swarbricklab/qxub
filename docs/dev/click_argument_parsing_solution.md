# Click Argument Parsing Solution for Unified CLI Interface

## Problem Statement

During the qxub 2.0 migration, we needed to transform the CLI from a subcommand-based architecture to a unified interface while preserving management commands. This created a complex Click framework challenge:

**Before (1.x):**
```bash
qxub conda --env myenv python script.py
qxub module --mods python3,gcc make
qxub sing --sif container.sif python script.py
qxub config --help  # Management command
```

**After (2.x):**
```bash
qxub --env myenv -- python script.py
qxub --mods python3,gcc -- make
qxub --sif container.sif -- python script.py
qxub --default -- echo "hello world"  # Default execution (no special environment)
qxub config --help  # Management command (unchanged)
```

## The Challenge

Click's argument parsing was designed for traditional subcommand patterns. When we tried to implement both execution options (`--env`, `--mod`, `--sif`, `--default`) and management subcommands in the same group, we encountered:

1. **Argument Capture Conflicts**: Using `@click.argument("command", nargs=-1)` captured everything, preventing Click from recognizing subcommands
2. **Subcommand Resolution Issues**: Click tried to resolve execution commands (like `echo`) as subcommands, causing "No such command" errors
3. **Protected Args vs Args**: The `--` separator correctly put commands in `ctx.protected_args`, but Click still attempted subcommand resolution

## Technical Details

### The Core Issue

```python
# This was the problem:
qxub --env base -- echo "test"
```

Click's parsing flow:
1. Parse options: `--env base` ✅
2. Encounter `--`: Put remaining in `protected_args: ['echo']` and `args: ['test']` ✅
3. Try to resolve `'echo'` as subcommand ❌ **FAILS HERE**

### The Solution: Custom Click Group

We implemented a custom `QxubGroup` class that overrides Click's command resolution behavior:

```python
class QxubGroup(click.Group):
    """Custom Click group with enhanced error handling for unknown options."""

    def get_command(self, ctx, cmd_name):
        """Override get_command to handle execution context."""
        # Check if we have execution context
        execution_options = ['env', 'mod', 'mods', 'sif']
        has_execution_context = any(ctx.params.get(opt) for opt in execution_options)

        if has_execution_context:
            return None  # Don't resolve as subcommand

        return super().get_command(ctx, cmd_name)

    def invoke(self, ctx):
        """Override invoke to handle execution contexts."""
        execution_options = ['env', 'mod', 'mods', 'sif']
        has_execution_context = any(ctx.params.get(opt) for opt in execution_options)

        if has_execution_context and ctx.protected_args:
            # Combine protected_args and args for execution
            combined_args = list(ctx.protected_args) + ctx.args

            # Temporarily set ctx.args for the main function
            original_args = ctx.args
            ctx.args = combined_args
            try:
                # Call the main qxub function directly
                qxub_func = ctx.command.callback.__wrapped__
                return qxub_func(
                    ctx,
                    ctx.params['execdir'],
                    ctx.params['verbose'],
                    # ... other parameters
                )
            finally:
                ctx.args = original_args

        return super().invoke(ctx)
```

### Key Insights

1. **`get_command()` Override**: When execution context is detected, return `None` to prevent subcommand resolution
2. **`invoke()` Override**: Manually combine `protected_args` and `args` to reconstruct the full command
3. **Direct Function Call**: Bypass Click's subcommand machinery and call the main function directly
4. **Parameter Extraction**: Carefully extract and pass parameters to avoid Click's `@pass_context` conflicts

## Argument Flow Visualization

### Execution Command Flow
```
qxub --env base -- echo "test"
    ↓
parse_args(): ['--env', 'base', '--', 'echo', 'test']
    ↓
Result: ctx.params['env']='base', ctx.protected_args=['echo'], ctx.args=['test']
    ↓
get_command('echo'): has_execution_context=True → return None
    ↓
invoke(): combine args=['echo', 'test'] → call qxub() directly
```

### Management Command Flow
```
qxub config --help
    ↓
parse_args(): ['config', '--help']
    ↓
Result: ctx.protected_args=['config'], ctx.args=['--help']
    ↓
get_command('config'): has_execution_context=False → return ConfigGroup
    ↓
invoke(): standard Click subcommand handling
```

## Alternative Approaches Considered

### 1. Separate Entry Points
**Approach**: Create separate entry points for execution vs management
**Rejected**: Would break backward compatibility and user experience

### 2. Context Detection in Main Function
**Approach**: Let Click parse normally, detect contexts in main function
**Problem**: Click errors occur before main function is reached

### 3. Argument Preprocessing
**Approach**: Preprocess `sys.argv` before Click sees it
**Problem**: Too fragile and would break Click's built-in features

## Implementation Benefits

1. **Backward Compatibility**: Management commands work unchanged
2. **Clean User Experience**: Unified interface as intended
3. **Robust Error Handling**: Click's error system still works for most cases
4. **Maintainable**: Solution is contained in the custom group class

## Lessons Learned

1. **Click is Powerful but Opinionated**: Great for traditional patterns, requires creativity for novel interfaces
2. **Method Override Strategy**: Sometimes the best approach is surgical overrides rather than fighting the framework
3. **Protected Args**: The `--` separator is crucial for complex argument parsing
4. **Context State Management**: Careful handling of Click's context state is essential

## Future Considerations

- Monitor Click library updates for potential conflicts
- Consider this pattern for other complex CLI interfaces
- Document any edge cases discovered during usage
- Possible upstream contribution to Click for better hybrid command support

## Related Files

- `qxub/cli.py`: Main implementation
- `docs/dev/migration_roadmap.md`: Overall migration context
- `tests/`: Test coverage for argument parsing scenarios

---
*This solution was developed during the qxub 2.0 migration (October 2024) to resolve Click framework limitations with hybrid command interfaces.*

## v2.2 Enhancement: Explicit Default Execution

In v2.2 (October 2025), the solution was further simplified by requiring an explicit `--default` flag for default execution:

```bash
# v2.0-2.1 (ambiguous)
qxub -- echo "hello"  # Could be confusing

# v2.2+ (explicit)
qxub --default -- echo "hello"  # Clear intent
```

**Benefits:**
- **Eliminates ambiguity**: No more protected_args fallback logic needed
- **Cleaner code**: Simplified `get_command()` and `invoke()` methods
- **Better UX**: Users must be explicit about execution intent
- **Consistent syntax**: All execution requires a flag (`--env`, `--mod`, `--sif`, `--default`)

This change makes the CLI more predictable and the codebase more maintainable.
