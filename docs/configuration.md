# Configuration Guide

`qxub` 2.0 introduces a comprehensive configuration system that eliminates repetitive command-line options and enables powerful workflow management.

## Quick Start

The configuration system lets you set defaults for commonly used options like project codes, queues, and resource requirements. Instead of typing `--project a56 --queue normal --resources mem=4GB` every time, you can set these as defaults and just run your commands.

### Basic Setup

```bash
# Create a configuration template
qxub config init

# Set common defaults
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"
qxub config set defaults.resources '["mem=4GB", "ncpus=1"]'

# View your configuration
qxub config list
```

### Everyday Usage

Once configured, your commands become much simpler:

```bash
# Before configuration:
qxub --project a56 --queue normal -l mem=4GB --env myenv -- python script.py

# After configuration:
qxub --env myenv -- python script.py  # Uses your defaults automatically
```

## Configuration Management

Managing your configuration is straightforward with these commands:

### Viewing Configuration

```bash
# Show all configuration
qxub config list

# Show specific sections
qxub config list defaults
qxub config list aliases

# Get specific values
qxub config get defaults.project
qxub config get defaults.resources
```

### Modifying Configuration

```bash
# Set single values
qxub config set defaults.name "myjob"
qxub config set defaults.project "xy12"

# Set complex values (use quotes for YAML)
qxub config set defaults.resources '["mem=8GB", "ncpus=2"]'

# Edit interactively
qxub config edit  # Opens in $EDITOR
```

### Validation and Management

```bash
# Validate configuration files
qxub config validate

# Show configuration file locations
qxub config show-files

# Reset user configuration (prompts for confirmation)
qxub config reset
```

## How Configuration Works

Configuration follows the XDG Base Directory specification with hierarchical precedence:
- **CLI arguments** (highest priority)
- **User config**: `~/.config/qxub/config.yaml`
- **System config**: `/etc/xdg/qxub/config.yaml` (optional)
- **Built-in defaults** (lowest priority)

This means you can set broad defaults in your user config, but still override them on the command line when needed.

## Configuration Structure

Your configuration file is written in YAML format and contains defaults for common options:

```yaml
defaults:
  name: "qt"
  queue: "normal"
  project: "a56"
  joblog: "{name}_{date}_{time}.log"
  resources: 
    - "mem=4GB"
    - "ncpus=1"
  out: "/scratch/{project}/{user}/qt/{timestamp}/out"
  err: "/scratch/{project}/{user}/qt/{timestamp}/err"

aliases:
  # See aliases.md for detailed alias configuration
```

### Available Default Options

| Option | Type | Description | Example |
|--------|------|-------------|---------|
| `name` | string | Default job name | `"analysis"` |
| `queue` | string | Default PBS queue | `"normal"` |
| `project` | string | Default project code | `"a56"` |
| `resources` | list | Default resource requirements | `["mem=4GB", "ncpus=1"]` |
| `joblog` | string | Job log file path | `"job_{timestamp}.log"` |
| `out` | string | STDOUT log path | `"/logs/{name}_out"` |
| `err` | string | STDERR log path | `"/logs/{name}_err"` |

## Template Variables

Configuration values support dynamic template variables resolved at runtime:

| Variable | Description | Example |
|----------|-------------|---------|
| `{user}` | Current username | `jr9959` |
| `{project}` | Project code | `a56` |
| `{name}` | Job name | `analysis` |
| `{queue}` | Queue name | `normal` |
| `{timestamp}` | Current timestamp | `20231005_143022` |
| `{date}` | Current date | `20231005` |
| `{time}` | Current time | `143022` |

### Template Examples

```yaml
defaults:
  name: "job_{date}"
  out: "/scratch/{project}/{user}/logs/{timestamp}/out"
  err: "/scratch/{project}/{user}/logs/{timestamp}/err"
  joblog: "/home/{user}/logs/{name}_{time}.log"
```

## Configuration Management Commands

### Viewing Configuration

```bash
# Show all configuration
qxub config list

# Show specific sections
qxub config list defaults
qxub config list aliases

# Get specific values
qxub config get defaults.project
qxub config get defaults.resources
```

### Modifying Configuration

```bash
# Set single values
qxub config set defaults.name "myjob"
qxub config set defaults.project "xy12"

# Set complex values (use quotes for YAML)
qxub config set defaults.resources '["mem=8GB", "ncpus=2"]'

# Edit interactively
qxub config edit  # Opens in $EDITOR
```

### Validation and Management

```bash
# Validate configuration files
qxub config validate

# Show configuration file locations
qxub config show-files

# Reset user configuration (prompts for confirmation)
qxub config reset
```

## Advanced Configuration

### Environment-Specific Configurations

You can create different configurations for different environments:

```bash
# Development configuration
qxub config set defaults.queue "express"
qxub config set defaults.resources '["mem=2GB", "ncpus=1"]'

# Production configuration  
qxub config set defaults.queue "normal"
qxub config set defaults.resources '["mem=16GB", "ncpus=4"]'
```

### Project-Specific Settings

Using template variables for project-specific paths:

```yaml
defaults:
  project: "bio01"
  out: "/scratch/{project}/{user}/analysis/{timestamp}/out"
  err: "/scratch/{project}/{user}/analysis/{timestamp}/err"
  joblog: "/home/{user}/projects/{project}/logs/{name}_{date}.log"
```

### Resource Templates

Create resource templates for different workload types:

```yaml
defaults:
  # Light workloads
  resources: ["mem=4GB", "ncpus=1", "walltime=01:00:00"]

# Override for specific workflows via aliases or CLI
# Heavy: --resources "mem=32GB,ncpus=8,walltime=06:00:00"
# GPU: --resources "mem=16GB,ncpus=12,ngpus=1,walltime=02:00:00" --queue gpuvolta
```

## Configuration Examples

### Data Science Workflow

```yaml
defaults:
  name: "ds_{date}"
  project: "ds01"
  queue: "normal"
  resources: ["mem=8GB", "ncpus=2", "walltime=02:00:00"]
  out: "/scratch/ds01/{user}/experiments/{timestamp}/out"
  err: "/scratch/ds01/{user}/experiments/{timestamp}/err"
  joblog: "/home/{user}/ds_logs/{name}_{time}.log"
```

### Bioinformatics Pipeline

```yaml
defaults:
  name: "bio_{timestamp}"
  project: "bio01"  
  queue: "normal"
  resources: ["mem=16GB", "ncpus=4", "walltime=04:00:00"]
  out: "/scratch/bio01/{user}/analysis/{date}/out"
  err: "/scratch/bio01/{user}/analysis/{date}/err"
```

### GPU Computing

```yaml
defaults:
  name: "gpu_{date}"
  project: "ml01"
  queue: "gpuvolta"
  resources: ["mem=32GB", "ncpus=12", "ngpus=1", "walltime=08:00:00"]
  out: "/scratch/ml01/{user}/training/{timestamp}/out"
  err: "/scratch/ml01/{user}/training/{timestamp}/err"
```

## Troubleshooting

### Common Issues

**Invalid YAML syntax:**
```bash
qxub config validate  # Check for syntax errors
```

**Template variables not resolving:**
- Ensure variables are spelled correctly: `{user}`, not `{username}`
- Check that referenced config values exist (e.g., `{project}` requires `defaults.project`)

**Permission errors:**
- Check that output directories exist and are writable
- Verify project codes are valid for your account

**Resource validation:**
- Use `qxub --dry-run` to preview generated commands
- Check PBS queue limits and requirements

### Getting Help

```bash
qxub config --help        # General config help
qxub config set --help    # Help for specific commands
qxub --dry-run <command>  # Preview what would be executed
```