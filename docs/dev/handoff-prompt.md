# Handoff Prompt for CLI Refactoring Completion

## Context
You are continuing CLI refactoring work on the qxub project at `/g/data/a56/software/qsub_tools` on NCI Gadi. The previous agent completed major consolidation work but hit file corruption issues with VS Code's `create_file` tool.

## Current State
- **Branch**: `feature/v3-cli-refactoring`
- **Location**: `/g/data/a56/software/qsub_tools`
- **Python Environment**: `source venv/bin/activate` (already set up)

## What's Been Completed ✅
1. **Template consolidation** → `qxub/templates.py` (working)
2. **Execution consolidation** → `qxub/unified_execution.py` (working)
3. **Config extraction** → `qxub/config_handler.py` (working)
4. **Clean CLI created** → `qxub/cli_new.py` (ready to deploy)

## What Needs To Be Done 🎯

### Immediate Task: Recreate exec_cli.py
The execution subcommand file was corrupted. You need to recreate it using **terminal commands only** (avoid `create_file`):

```bash
# Use cat with heredoc to create the file safely
cat > qxub/exec_cli.py << 'EOF'
[content from the specification below]
EOF
```

### Required exec_cli.py Content
Create a comprehensive Click command with:
- All PBS options: `-l`, `-q`, `-N`, `-P`, `-v`, `--dry`, `--quiet`, etc.
- All execution contexts: `--env`, `--mod`, `--mods`, `--sif`, `--default`
- Additional options: `--bind`, `--template`, `--pre`, `--post`, `--cmd`
- Command argument handling for both `-- command args` and `--cmd "string"` syntax
- Uses `ConfigHandler` and `execute_job` from our consolidated modules

### Final Steps After Creation
1. **Replace CLI**: `mv qxub/cli.py qxub/cli_old.py && mv qxub/cli_new.py qxub/cli.py`
2. **Test Interface**: `qxub --help` should show clean Click interface with exec command
3. **Test Execution**: `qxub exec --help` should show all execution options
4. **Smoke Tests**: Try basic execution commands with different contexts

## Critical Notes
- **DO NOT use `create_file` tool** - it corrupts files in this VS Code session
- **Use terminal commands** for any file creation (`cat`, `echo`, etc.)
- The consolidated modules (`templates.py`, `unified_execution.py`, `config_handler.py`) are working and tested
- Current `cli.py` is the original 737-line complex version with `QxubGroup`

## Success Criteria
- Clean `qxub --help` output showing standard Click interface
- `qxub exec` command available with comprehensive options
- All management commands still work (`config`, `alias`, `history`, etc.)
- No file corruption or syntax errors

## Reference Files
- See `docs/dev/cli-refactoring-progress.md` for detailed progress summary
- Original copilot instructions in `.github/copilot-instructions.md`
- Working consolidated modules already in place

Start by creating `exec_cli.py` with terminal commands, then proceed with CLI replacement and testing.
