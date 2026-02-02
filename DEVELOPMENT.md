# qxub v2.3 Development Branch

## 🚀 Welcome to v2.3 Development!

This branch is for developing **Phase 2.3: Enhanced Configuration & Pipeline Support** features.

### 🎯 Phase 2.3 Goals: Project Configuration & Pipeline Integration

The goal of v2.3 is to enhance qxub with project-level configuration management and improved pipeline integration capabilities.

### � Core Features Implemented

- **✅ Project-level configuration** with `.qx/` directory structure
- **✅ Multi-tier config hierarchy** (system → user → project → local → test)
- **✅ Git integration** for portable team configurations
- **✅ Terse output mode** for pipeline-friendly job submission
- **✅ Enhanced monitoring** with rich terminal UI and emoji status indicators
- **✅ Origin tracking** for configuration debugging

### 📋 Key Features in v2.3

1. **Project Configuration System**
   - **`.qx/` directory structure** for project-specific settings
   - **Three config types**: project (git-tracked), local (git-ignored), test (CI)
   - **Intelligent discovery** using git/dvc markers and project roots
   - **Configuration precedence**: CLI args > test > local > project > user > system

2. **Pipeline Integration**
   - **`--terse` flag** for machine-readable job ID output
   - **`qxub monitor` command** with rich terminal UI
   - **Parallel job submission** workflows
   - **Exit code propagation** for proper error handling

3. **Enhanced Configuration Management**
   - **Origin tracking** with `--show-origin` flag
   - **Multi-level set commands** (`--project`, `--local`, `--test`)
   - **Project initialization** with `qxub config init-project`
   - **Git integration** for portable team configurations

### 💡 Usage Examples

```bash
# Initialize project configuration
qxub config init-project

# Set team defaults (git-tracked)
qxub config set --project qxub.defaults.walltime "4:00:00"
qxub config set --project qxub.defaults.queue "normal"

# Personal overrides (git-ignored)
qxub config set --local qxub.defaults.queue "express"

# Pipeline integration
find *.csv -exec qxub --terse --env myenv -- process.py {} \; | qxub monitor
```

### 🛠️ Technical Architecture

Building on the solid foundation of v2.1's platform abstraction and v2.2's remote execution:

- **Project discovery**: Automatic detection of project roots using `.qx`, `.git`, `.dvc` markers
- **Configuration hierarchy**: Five-tier precedence system with clear override rules
- **Git integration**: Proper .gitignore handling for local vs shared configurations
- **Rich UI**: Terminal-based configuration display with origin tracking
- **Pipeline support**: Machine-readable output and comprehensive monitoring

### 📈 Success Criteria

- ✅ Team-shared project configurations that work across all developers
- ✅ Git-ignored personal overrides for individual workflow preferences
- ✅ CI/testing configurations that override other settings appropriately
- ✅ Clear visibility into configuration sources and precedence
- ✅ Pipeline-friendly job submission and monitoring workflows
- ✅ Backward compatibility with existing user/system configurations

---

**Configuration management made simple and powerful!** 🚀
