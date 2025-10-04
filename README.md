# qsub_tools

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/swarbricklab/qsub_tools)

This package provides tools for running `qsub` jobs on PBS Pro in particular environments.
This avoids the boilerplate code associated with activating environments and switching work directories that clutters up many jobscripts.
In simple cases, the need to create a jobscript can be eliminated entirely.

## New in Version 1.0: Configuration and Aliases

Version 1.0 introduces a powerful configuration and alias system that allows you to:
- **Configure defaults** for common options (project, queue, resources, etc.)
- **Create workflow aliases** for ultra-simple execution (`qxub alias dvc_push`)
- **Template-based paths** with variables like `{user}`, `{project}`, `{timestamp}`
- **Hierarchical config** (system → user → CLI arguments)

### Quick Examples

```bash
# Set up configuration defaults
qxub config init                          # Create config template
qxub config set defaults.project "xy12"  # Set default project
qxub config set defaults.queue "normal"  # Set default queue

# Create workflow aliases
qxub config alias set dvc_push --subcommand conda --cmd "dvc push" --env dvc3 --queue copyq
qxub config alias set analysis --subcommand module --mods "samtools,python3" --cmd "samtools view -c"

# Execute aliases with ultra-simple syntax
qxub alias dvc_push                       # Push data with DVC
qxub alias analysis input.bam            # Analyze BAM file
qxub alias dvc_push --queue normal        # Override queue for this run
```

## Installation

```bash
# Clone the repository
git clone https://github.com/swarbricklab/qsub_tools.git
cd qsub_tools

# Install in development mode
pip install -e .
```

## Quick Start

```bash
# Run a Python script in a conda environment
qxub conda --env myenv -- python script.py

# Run a command with environment modules
qxub module --mod bcftools --mod samtools -- bcftools --version

# Run a command in a Singularity container
qxub sing --sif /path/to/container.sif -- python analysis.py

# Use dry-run to see what would be executed
qxub --dry-run conda --env myenv -- python script.py
```

## Commands

The main CLI command provided by this tool is `qxub` (short for "extended qsub").
The `qxub` command is a wrapper around `qsub` and so accepts many of the same options, such as `-l mem=16GB`.

> The intention is to support all `qsub` options, but initially only the most common options have been implemented.
> Create a feature request if your favourite option is missing.

In addition to standard `qsub` options, `qxub` also has an `--execdir` option for controlling the work directory where the job is executed.
Unlike `qsub`, `qxub` defaults to executing in the current work directory.

### Common Options

- `--dry-run`: Show the `qsub` command that would be executed without actually submitting it
- `--quiet`: Submit job and exit immediately without monitoring output
- `--execdir`: Set the working directory for job execution (default: current directory)
- `--out`, `--err`: Specify custom paths for STDOUT and STDERR logs

## Configuration System

`qxub` 1.0 introduces a comprehensive configuration system that eliminates repetitive command-line options and enables powerful workflow aliases.

### Configuration Files

Configuration follows the XDG Base Directory specification:
- **User config**: `~/.config/qxub/config.yaml`
- **System config**: `/etc/xdg/qxub/config.yaml` (optional)
- **Precedence**: CLI arguments > User config > System config > Built-in defaults

### Getting Started with Configuration

```bash
# Create a configuration template
qxub config init

# View current configuration
qxub config list

# Set defaults
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"
qxub config set defaults.resources '["mem=4GB", "ncpus=1"]'

# Edit configuration directly
qxub config edit
```

### Template Variables

Configuration values support template variables that are resolved at runtime:

- `{user}` - Current username
- `{project}` - Project code
- `{name}` - Job name
- `{queue}` - Queue name
- `{timestamp}` - Current timestamp (YYYYMMDD_HHMMSS)
- `{date}` - Current date (YYYYMMDD)
- `{time}` - Current time (HHMMSS)

**Example Configuration:**
```yaml
defaults:
  name: "qt"
  queue: "normal"
  project: "a56"
  joblog: "{name}_{date}_{time}.log"
  resources: ["mem=4GB", "ncpus=1"]
  out: "/scratch/{project}/{user}/qt/{timestamp}/out"
  err: "/scratch/{project}/{user}/qt/{timestamp}/err"
```

### Configuration Management Commands

```bash
# View configuration
qxub config list                     # Show all config
qxub config list defaults           # Show defaults section
qxub config get defaults.project    # Get specific value

# Modify configuration
qxub config set defaults.name "myjob"       # Set single value
qxub config edit                             # Open in $EDITOR

# Validate and manage
qxub config validate                 # Check config files
qxub config show-files              # Show config file locations
qxub config reset                   # Reset user config
```

## Alias System

Aliases combine subcommands, commands, and configuration into reusable workflows with ultra-simple execution syntax.

### Creating Aliases

```bash
# Conda-based aliases
qxub config alias set dvc_push \
  --subcommand conda \
  --cmd "dvc push" \
  --env dvc3 \
  --queue copyq \
  --name push

# Module-based aliases
qxub config alias set analysis \
  --subcommand module \
  --cmd "samtools view -c" \
  --mod samtools \
  --name analysis

# Multiple modules
qxub config alias set pipeline \
  --subcommand module \
  --cmd "python analysis.py" \
  --mods "python3,samtools,gcc" \
  --name pipeline

# Singularity aliases
qxub config alias set container_job \
  --subcommand sing \
  --cmd "python analysis.py" \
  --sif "/path/to/container.sif" \
  --name container_analysis
```

### Executing Aliases

```bash
# Basic execution
qxub alias dvc_push                  # Execute as defined
qxub alias analysis input.bam       # Execute with command arguments

# Override options
qxub alias dvc_push --queue normal   # Override queue
qxub alias analysis --mod python3   # Override module
qxub alias pipeline --mods "R,python3"  # Override modules

# Global options
qxub --dry alias dvc_push           # Dry run
qxub --quiet alias analysis         # Quiet mode
```

### Alias Management

```bash
# List and inspect aliases
qxub config alias list              # List all aliases
qxub config alias show dvc_push     # Show alias details
qxub config alias test dvc_push     # Test and validate alias

# Modify aliases
qxub config alias set existing_alias --queue newqueue  # Update alias
qxub config alias delete old_alias  # Remove alias
```

### Advanced Alias Features

**Command Appending:**
```bash
# If alias defines: cmd: "samtools view -c"
qxub alias analysis input.bam output.bam
# Executes: samtools view -c input.bam output.bam
```

**Complete Command Override:**
```bash
qxub alias dvc_push --cmd "dvc pull"
# Uses dvc_push environment/settings but runs dvc pull instead
```

**Template Resolution in Aliases:**
```yaml
aliases:
  backup:
    subcommand: conda
    cmd: "rsync -av data/ /backup/{user}/{date}/"
    env: "tools"
    name: "backup_{timestamp}"
```

### Output

The STDOUT and STDERR produced by `qxub` jobs is redirected to the following logs:
- `$TMPDIR/qt/timestamp/out`
- `$TMPDIR/qt/timestamp/err`

The location of these logs can be specified via the `--out` and `--err` options.
However, unless `--quiet` is specified `qxub` will monitor these fails and stream their context back to STDOUT and STDERR in the shell where `qxub` was launched.
This makes it possible to view the output from the job in real time, without having to look up the logs (although the logs are still saved).
This also facilitates redirection of `qxub` output via traditional `>` and `2>` operators.

> Note: except for `--quiet` mode, `qxub` does not terminate until the job that it initiates terminates.
> But `qxub` only polls `qstat` every 30 seconds, and so `qxub` may appear to "hang" for up to 30 seconds even after the command has finished.
> In most cases, the time wasted waiting for the next poll will be less than the time required to look up the logs manually.

## Subcommands

`qxub` provides the following subcommands for executing commands or scripts in particular environments:
- `conda`: execute command or script in the specified conda environment
- `module`: activate specified environment modules before executing command or script
- `sing`: execute command or script in the specified Singularity container
- `config`: manage configuration settings and aliases
- `alias`: execute workflow aliases

### conda subcommand

Execute commands in a conda environment:

```bash
# Basic usage
qxub conda --env myenv -- python script.py

# With resource specifications
qxub --resources mem=16GB conda --env myenv -- python large_analysis.py

# With pre/post commands
qxub conda --env myenv --pre "echo Starting analysis" --post "echo Analysis complete" -- python script.py
```

### module subcommand

Execute commands with environment modules loaded:

```bash
# Single module
qxub module --mod bcftools -- bcftools --version

# Multiple modules (using comma-separated list)
qxub module --mods "bcftools,samtools,python3" -- python analysis.py

# Or multiple --mod options (still supported)
qxub module --mod bcftools --mod samtools --mod python3 -- python analysis.py

# With pre/post commands
qxub module --mod cuda --pre "nvidia-smi" -- python gpu_training.py
```

### sing subcommand

Execute commands in Singularity containers:

```bash
# Basic usage
qxub sing --sif /path/to/container.sif -- python script.py

# With Singularity options (before --)
qxub sing --sif container.sif --bind /data:/data --env DEBUG=1 -- python script.py

# With pre/post commands
qxub sing --sif container.sif --pre "module load singularity" -- analysis_tool --input data.txt
```

### config subcommand

Manage configuration and aliases:

```bash
# Configuration management
qxub config init                     # Create config template
qxub config list                     # Show all configuration
qxub config set defaults.project "a56"  # Set configuration value
qxub config edit                     # Edit config in $EDITOR

# Alias management  
qxub config alias list              # List all aliases
qxub config alias set myalias --subcommand conda --env myenv --cmd "python script.py"
qxub config alias show myalias      # Show alias details
qxub config alias delete myalias    # Remove alias
```

### alias subcommand

Execute workflow aliases:

```bash
# Basic execution
qxub alias myalias                   # Execute alias as defined
qxub alias myalias input.txt         # Execute with additional arguments

# Override alias options
qxub alias myalias --queue normal    # Override queue
qxub alias myalias --env newenv      # Override environment

# Global options apply
qxub --dry alias myalias            # Dry run the alias
qxub --quiet alias myalias          # Execute in quiet mode
```

## Signal Handling

All subcommands support graceful interruption with Ctrl+C:
- Automatically deletes submitted jobs with `qdel`
- Provides clear feedback during cleanup
- Ensures no orphaned jobs remain in the queue

## Common Workflows

Here are some practical examples of how to use `qxub` with the configuration and alias system:

### Data Science Pipeline

```bash
# Set up environment
qxub config set defaults.project "ds01"
qxub config set defaults.queue "normal"
qxub config set defaults.out "/scratch/ds01/{user}/logs/{timestamp}/out"

# Create analysis aliases
qxub config alias set preprocess \
  --subcommand conda \
  --env datasci \
  --cmd "python preprocess.py" \
  --name "preprocess_{date}"

qxub config alias set train_model \
  --subcommand conda \
  --env pytorch \
  --cmd "python train.py" \
  --resources "mem=32GB,ncpus=8,ngpus=1" \
  --queue gpuvolta \
  --name "training_{timestamp}"

qxub config alias set evaluate \
  --subcommand conda \
  --env datasci \
  --cmd "python evaluate.py" \
  --name "eval_{date}"

# Execute pipeline
qxub alias preprocess data/raw/
qxub alias train_model --cmd "python train.py --epochs 100"
qxub alias evaluate models/latest/
```

### Bioinformatics Analysis

```bash
# Set up bioinformatics environment
qxub config set defaults.project "bio01"
qxub config set defaults.resources '["mem=16GB", "ncpus=4"]'

# Create analysis aliases
qxub config alias set qc_check \
  --subcommand module \
  --mods "fastqc,multiqc" \
  --cmd "fastqc" \
  --name "qc_{date}"

qxub config alias set align \
  --subcommand module \
  --mods "bwa,samtools" \
  --cmd "bwa mem ref.fa" \
  --resources "mem=32GB,ncpus=16" \
  --name "align_{timestamp}"

qxub config alias set variant_call \
  --subcommand sing \
  --sif "/containers/gatk.sif" \
  --cmd "gatk HaplotypeCaller" \
  --name "variants_{date}"

# Execute analysis
qxub alias qc_check reads.fastq
qxub alias align reads.fastq
qxub alias variant_call -I aligned.bam -R ref.fa
```

### Data Management

```bash
# Set up data management aliases
qxub config alias set dvc_push \
  --subcommand conda \
  --env dvc3 \
  --cmd "dvc push" \
  --queue copyq \
  --name "push_{time}"

qxub config alias set backup \
  --subcommand module \
  --mod rsync \
  --cmd "rsync -av" \
  --queue copyq \
  --name "backup_{date}"

qxub config alias set sync_remote \
  --subcommand conda \
  --env tools \
  --cmd "rclone sync" \
  --queue copyq \
  --name "sync_{timestamp}"

# Execute data operations
qxub alias dvc_push
qxub alias backup data/ /backup/{user}/
qxub alias sync_remote local/ remote:bucket/
```

### Container-based Workflows

```bash
# Set up container aliases
qxub config alias set rnaseq \
  --subcommand sing \
  --sif "/containers/rnaseq.sif" \
  --bind "/data:/data,/scratch:/scratch" \
  --env "THREADS=8" \
  --cmd "nextflow run nf-core/rnaseq" \
  --resources "mem=64GB,ncpus=16" \
  --queue express

qxub config alias set jupyter \
  --subcommand sing \
  --sif "/containers/jupyter.sif" \
  --bind "/home:/home,/data:/data" \
  --cmd "jupyter lab --no-browser --ip=0.0.0.0" \
  --name "jupyter_{user}"

# Execute container workflows
qxub alias rnaseq --input data/reads/ --outdir results/
qxub alias jupyter --port 8888
```

## Changelog

### Version 1.0.0 - Major Release: Configuration and Alias System

**🎉 Major Features:**
- **Configuration Management**: XDG-compliant config system with hierarchical precedence
- **Workflow Aliases**: Ultra-simple execution with `qxub alias name` syntax
- **Template Variables**: Dynamic path resolution with `{user}`, `{project}`, `{timestamp}`, etc.
- **Enhanced Module Options**: New `--mods` for comma-separated module lists alongside `--mod`

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

**💡 Enhanced CLI:**
- All defaults now configurable (project, queue, resources, paths, etc.)
- Automatic template resolution for paths and job names
- Improved module handling: `--mod` (single) and `--mods` (comma-separated)
- Better error messages and validation
- Comprehensive help text and examples

**📋 Example Workflows:**
```bash
# Set up common defaults
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# Create reusable aliases
qxub config alias set dvc_push --subcommand conda --cmd "dvc push" --env dvc3 --queue copyq
qxub config alias set analysis --subcommand module --mods "samtools,python3" --cmd "samtools view -c"

# Execute with simple syntax
qxub alias dvc_push
qxub alias analysis input.bam
qxub alias dvc_push --queue normal  # Override queue
```

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
