# qsub_tools

[![Version](https://img.shields.io/badge/version-1.0.1-blue.svg)](https://github.com/swarbricklab/qsub_tools)

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
# Run in conda environment (uses active environment by default)
qxub conda -- python script.py
qxub conda --env myenv -- python script.py

# Run with environment modules  
qxub module --mod samtools -- samtools --version
qxub module --mods "python3,gcc" -- python analysis.py

# Run in Singularity container
qxub sing --sif container.sif -- python script.py

# Add PBS options (main options go before subcommand)
qxub --resources mem=16GB --queue normal conda --env myenv -- python script.py

# Preview commands without running
qxub --dry-run conda --env myenv -- python script.py
```

## Advanced Features

For more complex workflows, qxub provides configuration and alias systems:

```bash
# Set up defaults (optional)
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# Create workflow aliases (optional)
qxub config alias set dvc_push --subcommand conda --cmd "dvc push" --env dvc3 --queue copyq

# Execute aliases
qxub alias dvc_push
```

**📖 See detailed guides:**
- **[Configuration Guide](docs/configuration.md)** - Set up defaults and template variables
- **[Alias Guide](docs/aliases.md)** - Create reusable workflow shortcuts

## Important Notes

> ⚠️ **Key things to remember:**
>
> 1. **Option placement**: Main options (`--name`, `--queue`, `--resources`) must come **before** the subcommand
> 2. **GPU queues**: GPU queues like `gpuvolta` require both GPU and CPU specifications: `--resources ngpus=1,ncpus=12`
> 3. **Global flags**: Options like `--dry-run` must come immediately after `qxub`

## Subcommands

| Command | Description |
|---------|-------------|
| `conda` | Execute in conda environment |
| `module` | Execute with environment modules |
| `sing` | Execute in Singularity container |
| `config` | Manage configuration (see [docs](docs/configuration.md)) |
| `alias` | Execute workflow aliases (see [docs](docs/aliases.md)) |

## Common Options

| Option | Description | Example |
|--------|-------------|---------|
| `--dry-run` | Preview command without executing | `qxub --dry-run conda --env myenv script.py` |
| `--quiet` | Submit and exit (no monitoring) | `qxub --quiet conda --env myenv script.py` |
| `--resources` | PBS resource requirements | `--resources mem=16GB,ncpus=4` |
| `--queue` | PBS queue name | `--queue normal` |
| `--name` | Job name | `--name analysis_job` |

## Output and Monitoring

- **Real-time monitoring**: By default, `qxub` monitors job output and streams it to your terminal
- **Log files**: STDOUT/STDERR saved to configurable paths (default: `$TMPDIR/qt/{timestamp}/`)
- **Quiet mode**: Use `--quiet` to submit and exit immediately
- **Signal handling**: Ctrl+C automatically cleans up submitted jobs

## Getting Help

```bash
qxub --help                    # General help
qxub conda --help             # Subcommand help
qxub config --help            # Configuration help
```

## Version 1.0.0 Highlights

**🎉 Major Features:**
- **Hierarchical configuration** with defaults, template variables, and CLI overrides
- **Workflow aliases** for ultra-simple execution (`qxub alias name`)
- **Multi-environment support** - conda, modules, Singularity containers
- **Smart error handling** with consistent exit codes
- **Template variables** for dynamic paths and job names

**🔧 Bug Fixes:**
- Fixed option placement validation (main options before subcommands)
- Fixed GPU queue resource requirements
- Fixed alias option override logic
- Consistent exit code 2 for all validation errors

See full [CHANGELOG](#changelog) below for complete details.

---

## Changelog

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