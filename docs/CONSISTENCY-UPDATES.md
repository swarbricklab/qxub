# Documentation Consistency Updates

## Summary of Changes Made

Updated qxub v2.2 documentation to ensure consistency across all files regarding remote execution terminology and architecture.

## Key Terminology Changes

### CLI Flag Consistency
- **OLD**: `--profile` (from early design)
- **NEW**: `--remote` (final design)

### Architecture Terminology
- **OLD**: "Platform profiles" combining host + platform + credentials
- **NEW**: "Remote configurations" for connection details, separate from platform definitions

### Environment Setup
- **OLD**: Complex "remote_setup_commands" with module loading
- **NEW**: Simple "qxub_env" conda environment activation

## Files Updated

### ✅ `/DEVELOPMENT.md`
- Updated usage examples from `--profile` to `--remote`
- Updated architecture description to reflect simplified design
- Removed references to "distributed logging" in favor of simplified approach
- Updated feature list to reflect conda-based setup

### ✅ `/docs/dev/v2_roadmap.md`
- Changed example from `--profile gadi` to `--remote nci_gadi`

### ✅ `/docs/dev/config_schema.md`
- Updated migration notes to reference `--remote` instead of `--profile`
- Updated validation rules to reference "remote names" instead of "profile names"
- Updated opt-in features section

### ✅ `/docs/dev/remote_execution.md`
- Updated code examples from `--profile` to `--remote`
- Changed setup commands from module loading to conda activation
- Updated platform specification to use `--platform-file`
- Updated error messages and suggestions

## Verified Consistent Documentation

### ✅ New Documentation (Already Consistent)
- `/docs/remote-configuration.md` - Uses `--remote` throughout
- `/docs/remote-execution.md` - Uses `--remote` throughout
- `/docs/ssh-configuration.md` - No CLI flags, SSH-focused
- `/docs/README-remote.md` - Uses `--remote` throughout
- `/docs/platforms/nci_gadi.yaml` - Clean platform definition

## Architecture Consistency Verified

### Configuration Separation
✅ **Platform Definitions**: Describe system capabilities (queues, limits)
- Location: Remote systems (e.g., `/g/data/a56/.../platforms/nci_gadi.yaml`)
- Content: Universal system capabilities, no user-specific data

✅ **User Configuration**: Connection and execution preferences
- Location: Local machine (`~/.config/qxub/config.yaml`)
- Content: SSH hosts, conda environments, working directory templates

✅ **SSH Configuration**: Authentication and connection details
- Location: Local machine (`~/.ssh/config`)
- Content: Credentials, timeouts, connection optimization

### CLI Consistency
✅ **Remote Execution**: `qxub --remote REMOTE_NAME [options] -- command`
✅ **Local Execution**: `qxub --platform PLATFORM_NAME [options] -- command`
✅ **Platform Files**: `qxub --platform-file /path/to/platform.yaml [options] -- command`

## Key Design Principles Maintained

1. **🔒 No Credential Storage**: qxub never stores SSH credentials
2. **🐍 Conda-Only Setup**: Simple environment activation, no complex scripts
3. **📁 Clean Separation**: Connection details vs. system capabilities
4. **🔧 SSH Delegation**: All connection management via standard SSH tools
5. **📍 Explicit Paths**: Clear platform file locations, no magic discovery

## Verification Commands

To verify consistency, these commands should return no results:

```bash
# Check for old --profile usage in new docs
grep -r "--profile" docs/remote-*.md docs/README-remote.md docs/platforms/

# Check for old terminology in main development docs
grep -r "platform profile" DEVELOPMENT.md docs/dev/

# Check for complex setup commands in new architecture
grep -r "module load\|remote_setup_commands" docs/remote-*.md
```

All documentation now consistently reflects the simplified v2.2 remote execution architecture with `--remote` CLI flag and conda-based environment setup.
