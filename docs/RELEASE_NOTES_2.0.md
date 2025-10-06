# qxub 2.0.0 Release Notes

## 🎉 Major Release: Unified CLI Architecture

**Release Date**: October 2025  
**Version**: 2.0.0  
**Breaking Changes**: Yes - see migration guide below

## Overview

qxub 2.0 represents a major architectural improvement that eliminates subcommands in favor of a unified CLI with execution context options. This change significantly improves user experience by providing a cleaner, more intuitive interface while maintaining all existing functionality.

## 🚀 Major Changes

### ✨ **Unified CLI Interface**
Replaced complex subcommand structure with intuitive execution context options:

```bash
# Before (qxub 1.x)
qxub conda --env myenv python script.py
qxub module --mod python3,gcc python script.py  
qxub sing --sif container.sif python script.py

# After (qxub 2.0)
qxub --env myenv -- python script.py
qxub --mod python3 --mod gcc -- python script.py
qxub --sif container.sif -- python script.py
```

### 🎯 **Multiple Option Names**
Added alternative option names for improved usability:

- **Conda**: `--env` or `--conda`
- **Modules**: `--mod` (repeatable), `--mods`, or `--modules` (comma-separated)
- **Singularity**: `--sif`, `--sing`, or `--singularity`

### ⚡ **Default Execution Context**
New capability for environment-agnostic job submission:

```bash
# Direct PBS submission without environment activation
qxub -- echo "Hello PBS"
qxub --pre "module load python3" -- python script.py
```

### 🔒 **Improved Validation**
- Mutual exclusivity enforcement between execution contexts
- Clear error messages for invalid option combinations
- Better handling of edge cases

### 📊 **Enhanced History System**
- Updated recipe extraction for unified CLI format
- Maintains compatibility with existing history data
- Improved execution tracking across all contexts

## 🔧 **Technical Improvements**

### **Performance**
- No regression in job submission speed
- Optimized option parsing and validation
- Efficient mutual exclusivity checks

### **Testing**
- Comprehensive test suite: 39/39 tests passing
- Real-world scenario validation
- Edge case coverage for all execution contexts

### **Documentation**
- Complete overhaul of all examples
- Updated alias system documentation
- Improved help system organization

## 📋 **Migration Guide**

### **For Individual Users**

#### **Command Syntax Changes**
All subcommands have been removed. Update your commands as follows:

```bash
# Conda Environment Execution
# Old: qxub conda --env myenv python script.py
# New: qxub --env myenv -- python script.py
# Alt: qxub --conda myenv -- python script.py

# Module Environment Execution  
# Old: qxub module --mod python3,gcc python script.py
# New: qxub --mods python3,gcc -- python script.py
# Alt: qxub --mod python3 --mod gcc -- python script.py

# Singularity Container Execution
# Old: qxub sing --sif container.sif python script.py  
# New: qxub --sif container.sif -- python script.py
# Alt: qxub --sing container.sif -- python script.py

# All PBS options remain the same
qxub --env myenv --queue normal --resources mem=8GB -- python script.py
```

#### **Important Notes**
- **`--` separator required**: Separates qxub options from your command options
- **All PBS options unchanged**: `--queue`, `--resources`, `--name`, etc. work exactly the same
- **Management commands unchanged**: `qxub config`, `qxub alias`, `qxub history` work as before

### **For Scripts and Automation**
Search and replace patterns for automated migration:

```bash
# Basic replacements
sed -i 's/qxub conda --env \([^ ]*\)/qxub --env \1 --/g' your_script.sh
sed -i 's/qxub module --mod \([^ ]*\)/qxub --mod \1 --/g' your_script.sh  
sed -i 's/qxub sing --sif \([^ ]*\)/qxub --sif \1 --/g' your_script.sh
```

### **Alias System**
The alias system has been updated to work with the unified CLI. Existing aliases will be automatically converted to the new format when first accessed.

## 🆕 **New Features**

### **Default Execution Template**
Execute commands directly with PBS without any environment activation:

```bash
# Simple command execution
qxub -- echo "Hello from PBS"

# With pre/post commands
qxub --pre "module load python3" --post "echo Done" -- python script.py

# With full PBS options
qxub --queue gpuvolta --resources ngpus=1 -- python gpu_script.py
```

### **Enhanced Module Support**
Multiple ways to specify environment modules:

```bash
# Repeatable single options
qxub --mod python3 --mod gcc --mod cmake -- make

# Comma-separated list
qxub --mods python3,gcc,cmake -- make
qxub --modules python3,gcc,cmake -- make  # alternative name
```

### **Flexible Container Options**
Multiple names for Singularity containers:

```bash
qxub --sif container.sif -- command
qxub --sing container.sif -- command  
qxub --singularity container.sif -- command
```

## 🔄 **Backward Compatibility**

### **What's Preserved**
- All PBS options work exactly the same
- Configuration files unchanged
- History data compatibility maintained
- Alias functionality enhanced (not broken)
- Template system unchanged

### **What's Changed**
- **CLI syntax**: Subcommands removed, execution contexts now options
- **Alias structure**: Simplified flat structure (automatically converted)
- **Command separation**: `--` separator now required

## 🚨 **Breaking Changes**

### **Removed Features**
- **Subcommands**: `qxub conda`, `qxub module`, `qxub sing` no longer exist
- **Complex option placement rules**: Simplified with `--` separator

### **Required Updates**
1. **All command invocations**: Must update to new syntax
2. **Scripts**: Need find/replace for subcommand usage
3. **Documentation**: Update any internal docs with examples

## 🐛 **Bug Fixes**
- Fixed option placement confusion with clear `--` separator
- Resolved mutual exclusivity validation issues
- Improved error messages for invalid combinations
- Enhanced resource specification validation

## 📈 **Performance**
- **No regression** in job submission speed
- **Optimized** option parsing and validation
- **Reduced memory usage** in CLI processing
- **Faster** error detection and reporting

## 🔧 **Development**
- **Comprehensive testing**: 39/39 tests passing
- **Code quality**: Black formatting applied throughout
- **Documentation**: Complete examples update
- **Type safety**: Improved type hints and validation

## 🎯 **User Benefits**

### **Simplified Interface**
- **No more subcommand confusion**: Clear execution context options
- **Intuitive syntax**: `--env` clearly indicates conda execution
- **Flexible options**: Multiple names for the same functionality

### **Enhanced Functionality** 
- **Default execution**: New capability for environment-agnostic jobs
- **Better validation**: Clear error messages for invalid combinations
- **Improved help**: Reorganized and more accessible

### **Maintained Power**
- **All existing features**: Every capability preserved and enhanced
- **Same performance**: No regression in job submission speed
- **Enhanced history**: Better tracking with new command structure

## 📚 **Documentation**

### **Updated Resources**
- [Main README](../README.md) - Complete examples update
- [Aliases Guide](docs/aliases.md) - New structure documentation
- [Configuration Guide](docs/configuration.md) - Updated for 2.0
- [Migration Roadmap](docs/dev/migration_roadmap.md) - Complete implementation status

### **Getting Help**
- `qxub --help` - Updated for unified interface
- `qxub config --help` - Configuration management
- `qxub alias --help` - Alias system help
- `qxub history --help` - Execution history help

## 🙏 **Acknowledgments**

Special thanks to our three active users for their patience during this major architectural improvement. The feedback and testing during development were invaluable.

## 🔮 **What's Next**

### **Future Enhancements**
- Enhanced default execution context detection
- Additional container runtime support
- Improved resource optimization features
- Extended history and analytics capabilities

### **Community**
- Continue gathering feedback for minor improvements
- Enhance documentation based on user experience
- Explore additional productivity features

---

## 📞 **Support**

For questions, issues, or feedback about qxub 2.0:

- **GitHub Issues**: [qsub_tools repository](https://github.com/swarbricklab/qsub_tools)
- **Direct Support**: Contact the development team
- **Documentation**: Check the updated docs in this repository

---

**qxub 2.0 represents a significant step forward in job submission simplicity and power. We're excited to see how it improves your computational workflows!** 🚀