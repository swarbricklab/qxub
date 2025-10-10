# qsub_tools

[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)](https://github.com/swarbricklab/qsub_tools)

A powerful wrapper around PBS `qsub` that eliminates boilerplate when running jobs in conda environments, with modules, or in containers.

## Quick Start

### Installation

```bash
git clone https://github.com/swarbricklab/qsub_tools.git
cd qsub_tools
pip install -e .
```

### Basic Usage

```bash
# Run in conda environment
qxub --env myenv -- python script.py
qxub --conda myenv -- python script.py  # alternative option name

# Run with environment modules
qxub --mod samtools -- samtools --version
qxub --mods "python3,gcc" -- python analysis.py  # comma-separated
qxub --mod python3 --mod gcc -- python analysis.py  # repeatable

# Run in Singularity container
qxub --sif container.sif -- python script.py
qxub --singularity container.sif -- python script.py  # alternative

# Run without specific environment (direct PBS submission)
qxub -- echo "Hello world"
qxub --pre "echo Starting" --post "echo Done" -- python script.py

# Add PBS options (all options go before --)
qxub -l mem=16GB --queue normal --env myenv -- python script.py

# Preview commands without running
qxub --dry --env myenv -- python script.py
```

## Advanced Features

For more complex workflows, qxub provides configuration and alias systems:

```bash
# Set up defaults (optional)
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# Create workflow aliases (optional)
qxub config alias set dvc_push --env dvc3 --queue copyq --cmd "dvc push"

# Execute aliases
qxub alias dvc_push
```

**📖 See detailed guides:**
- **[Configuration Guide](docs/configuration.md)** - Set up defaults and template variables
- **[Alias Guide](docs/aliases.md)** - Create reusable workflow shortcuts

## Important Notes

> ⚠️ **Key things to remember:**
>
> 1. **Option placement**: All qxub options must come **before** the `--` separator
> 2. **Command separation**: Use `--` to separate qxub options from your command and its arguments
> 3. **GPU queues**: GPU queues like `gpuvolta` require both GPU and CPU specifications: `-l ngpus=1,ncpus=12`
> 4. **Alternative options**: Use `--env`/`--conda`, `--mod`/`--mods`/`--modules`, `--sif`/`--sing`/`--singularity` based on preference

## Execution Contexts

| Context | Options | Description |
|---------|---------|-------------|
| **Conda** | `--env`, `--conda` | Execute in conda environment |
| **Modules** | `--mod` (repeatable), `--mods`, `--modules` | Execute with environment modules |
| **Singularity** | `--sif`, `--sing`, `--singularity` | Execute in Singularity container |
| **Default** | *(no context option)* | Direct PBS submission without environment |

## Common Options

| Option | Description | Example |
|--------|-------------|---------|
| `--dry` | Preview command without executing | `qxub --dry --env myenv -- script.py` |
| `--quiet` | Submit and exit (no monitoring) | `qxub --quiet --env myenv -- script.py` |
| `-l` | PBS resource requirements | `-l mem=16GB -l ncpus=4` |
| `--queue` | PBS queue name or 'auto' for intelligent selection | `--queue normal` or `--queue auto` |
| `--name` | Job name | `--name analysis_job` |
| `--pre` | Command to run before main command | `--pre "echo Starting"` |
| `--post` | Command to run after main command | `--post "echo Done"` |

## Management Commands

| Command | Description |
|---------|-------------|
| `qxub config` | Manage configuration (see [docs](docs/configuration.md)) |
| `qxub alias` | Execute workflow aliases (see [docs](docs/aliases.md)) |
| `qxub history` | View execution history |
| `qxub resources` | View resource usage |

## Intelligent Queue Selection

Use `--queue auto` to enable intelligent queue selection based on your resource requirements:

```bash
# Auto-select queue based on resources
qxub --queue auto -l mem=2GB -l ncpus=2 --env myenv -- python small_job.py
# → Automatically selects 'small' queue for lightweight jobs

qxub --queue auto -l mem=500GB --env myenv -- python big_memory_job.py
# → Automatically selects 'hugemem' queue for high memory requirements

qxub --queue auto -l ngpus=1 -l ncpus=12 --env pytorch -- python train.py
# → Automatically selects 'gpu' queue for GPU workloads
```

### Auto-Selection Rules

When `--queue auto` is used, qxub intelligently selects the best queue based on:

- **Memory requirements**: Jobs requiring >100GB automatically use high-memory queues
- **GPU requirements**: Jobs with GPU specifications automatically use GPU queues
- **CPU requirements**: Small jobs (≤2 CPUs, ≤2GB) use express/small queues when available
- **Cost optimization**: Among valid queues, selects the most cost-effective option
- **Fallback**: Uses 'normal' queue if no specific rules match

> **Note**: Auto-queue selection requires platform definitions. Use `QXUB_PLATFORM_PATHS` environment variable to specify custom platform definitions.

## Output and Monitoring

- **Real-time monitoring**: By default, `qxub` monitors job output and streams it to your terminal
- **Log files**: STDOUT/STDERR saved to configurable paths (default: `$TMPDIR/qt/{timestamp}/`)
- **Quiet mode**: Use `--quiet` to submit and exit immediately
- **Signal handling**: Ctrl+C automatically cleans up submitted jobs

## Getting Help

```bash
qxub --help                    # General help
qxub config --help            # Configuration help
qxub alias --help             # Alias management help
qxub history --help           # History management help
qxub resources --help         # Resource tracking help
```

## Examples

### Real-world Usage Patterns

```bash
# Data science workflow
qxub --conda pytorch --queue gpuvolta -l ngpus=1 -l ncpus=12 --name gpu-training -- python train.py --epochs 100

# Bioinformatics pipeline
qxub --modules samtools,bwa -l walltime=02:00:00 -l mem=16GB -- bash pipeline.sh

# Container workflow
qxub --singularity /containers/blast.sif --bind /data --name blast-search -- blastn -query input.fa -db nt

# Simple script with setup and cleanup
qxub --pre "mkdir -p /tmp/work" --post "rm -rf /tmp/work" -- python analyze_data.py

# Direct PBS submission (no environment)
qxub -l walltime=01:00:00 -l mem=8GB -- ./my_compiled_program
```

---

## Changelog

### Version 2.0.0 - Unified CLI Architecture

**🎯 Major Features:**
- **Unified Command Interface**: Eliminated subcommands in favor of execution context options (`--env`, `--mod`, `--sif`)
- **Alternative Option Names**: Multiple option styles for better usability (`--env`/`--conda`, `--sif`/`--sing`/`--singularity`)
- **Default Execution Template**: Direct PBS job submission without environment activation
- **Enhanced Option Handling**: All options before `--` separator for unambiguous parsing
- **Comprehensive Testing**: 39 automated tests ensuring reliability

**🔧 Technical Improvements:**
- **Simplified Architecture**: No more subcommand complexity
- **Mutual Exclusivity Validation**: Prevents conflicting execution contexts
- **Black Code Formatting**: Consistent code style across entire codebase
- **Context Initialization Fixes**: Proper handling of execution contexts

**⚠️ Breaking Changes:**
- **CLI Syntax**: Changed from `qxub conda --env myenv -- cmd` to `qxub --env myenv -- cmd`
- **Option Placement**: All qxub options must come before `--` separator
- **Subcommand Removal**: `conda`, `module`, and `sing` subcommands no longer exist

**✨ New Syntax Examples:**
```bash
# Before (1.x)                        # After (2.0)
qxub conda --env myenv -- python script.py    →    qxub --env myenv -- python script.py
qxub module --mod python3 -- python script.py    →    qxub --mod python3 -- python script.py
qxub sing --sif container.sif -- python script.py    →    qxub --sif container.sif -- python script.py
```
### Version 2.1.0 - Platform Abstraction and Intelligent Queue Selection

**🎯 Major Features:**
- **Intelligent Queue Selection**: Use `--queue auto` for automatic queue selection based on resource requirements
- **Platform Abstraction System**: YAML-based platform definitions with queue limits, validation, and cost estimation
- **Environment Variable Support**: `QXUB_PLATFORM_PATHS` for custom platform definition discovery
- **Resource Validation**: Automatic validation of resource requests against platform queue limits
- **Cost Optimization**: Auto-selection prioritizes cost-effective queues among valid options

**🔧 Platform System:**
- **Platform Loader**: Automatic discovery and loading of platform YAML definitions
- **Queue Management**: Comprehensive queue definitions with limits, walltime rules, and constraints
- **Auto-Selection Engine**: Rules-based selection with memory, GPU, and CPU requirement handling
- **Resource Utilities**: Parse and compare memory sizes, walltime formats, and resource conditions

**🧪 Testing Infrastructure:**
- **Comprehensive Test Suite**: 31 test cases with 100% pass rate covering all platform functionality
- **Integration Testing**: Real-world auto-queue selection scenarios with multiple platform types
- **Debug Scripts**: Extensive debugging and validation tools for platform system development

### Version 1.1.0 - Enhanced CLI, Resource Tracking, and System Improvements

**🎯 Major Features:**
- **Enhanced CLI Error Handling**: Intelligent error messages with fuzzy matching for typos and contextual guidance for option placement
- **Resource Efficiency Tracking**: SQLite-based system for analyzing job performance with efficiency metrics and trend analysis
- **Comprehensive Threading Improvements**: Fixed infinite loop bugs, proper exit code propagation, and enhanced job monitoring
- **Dual-Log History System**: Separate computational recipes from execution records for better analysis and alias conversion
- **System Configuration Support**: XDG-compliant system-wide configurations for organizational defaults

**🔧 Technical Enhancements:**
- **Custom Click Classes**: QxubGroup and QxubCommand with intelligent error handling and suggestion system
- **Resource Tracking Integration**: Automatic efficiency calculations integrated with job completion workflow
- **Enhanced Documentation**: Comprehensive threading architecture docs and development guides
- **Robust Testing**: System configuration tests and realistic usage scenario validation

### Version 1.0.1 - Documentation and Linting Improvements

**🔧 Improvements:**
- **Migrated from Pylint to Black**: Faster, simpler code formatting with industry standard tools
- **Organized Documentation**: Restructured docs into dedicated guides (configuration, aliases, option placement)
- **Comprehensive Integration Tests**: Added 28 test cases covering all CLI patterns and edge cases
- **Enhanced Option Placement Guide**: Clear rules and examples for CLI syntax
- **CI Optimization**: Faster GitHub Actions workflow with Black formatting checks

### Version 1.0.0 - Major Release: Configuration and Alias System

**🎉 Major Features:**
- **Configuration Management**: XDG-compliant config system with hierarchical precedence
- **Workflow Aliases**: Ultra-simple execution with `qxub alias name` syntax
- **Template Variables**: Dynamic path resolution with `{user}`, `{project}`, `{timestamp}`, etc.
- **Enhanced Module Options**: New `--mods` for comma-separated module lists alongside `--mod`
- **Hierarchical Alias Structure**: Separates main qxub options from subcommand options for proper CLI syntax

**🔧 Configuration System:**
- `qxub config init` - Create configuration template
- `qxub config get/set/list` - Manage configuration values
- `qxub config edit` - Edit config in $EDITOR
- `qxub config validate` - Validate configuration files
- Template variables in all config values
- Hierarchical config: CLI args > User config > System config > Defaults

**🚀 Alias System:**
- `qxub config alias set/list/show/delete` - Manage workflow aliases
- `qxub alias name` - Execute aliases with ultra-simple syntax
- Option overrides: `qxub alias name --queue normal`
- Command appending: `qxub alias name input.txt output.txt`
- Complete command override: `qxub alias name --cmd "new command"`
- Support for all subcommands (conda, module, sing)
- New hierarchical structure prevents CLI option placement errors

**💡 Enhanced CLI:**
- All defaults now configurable (project, queue, resources, paths, etc.)
- Automatic template resolution for paths and job names
- Improved module handling: `--mod` (single) and `--mods` (comma-separated)
- **Consistent error handling**: All validation errors return exit code 2
- **Improved environment handling**: `conda` subcommand uses active environment by default, provides clear error when none available
- **PBS error mapping**: Invalid resource formats now consistently return exit code 2
- Better error messages and validation
- Comprehensive help text and examples

**🐛 Critical Bug Fixes:**
- Fixed option placement issues: main options (--name, --queue, --resources) must precede subcommands
- Fixed GPU queue validation: gpuvolta queue now properly requires GPU and CPU specifications
- Fixed null resource handling to prevent TypeErrors
- Fixed alias option override logic to prevent empty values from overriding alias defaults
- Fixed exit code consistency for script integration

**⚠️ Breaking Changes:**
- Alias configuration structure changed to hierarchical format (main/subcommand/target sections)
- Option placement now strictly enforced (main options before subcommand)
- Some error exit codes changed from 1 to 2 for consistency

### Version 0.3.0
- Added `module` subcommand for environment module support
- Added `sing` subcommand for Singularity container support
- Added `--pre` and `--post` options for command chaining
- Added Ctrl+C signal handling with automatic job cleanup
- Improved spinner display and output coordination
- Added comprehensive argument passthrough for Singularity options

### Version 0.2.0
- Initial release with `conda` subcommand
- Basic PBS job submission and monitoring functionality
