# qsub_tools Documentation

Welcome to the complete documentation for qsub_tools 2.0! 

## Quick Navigation

### User Documentation
- **[Configuration Guide](configuration.md)** - Complete guide to setting up defaults, template variables, and config management
- **[Alias Guide](aliases.md)** - Create and manage workflow aliases for ultra-simple execution
- **[Option Placement Guide](option-placement.md)** - Understand command-line option ordering rules
- **[README](../README.md)** - Quick start and essential information

### Developer Documentation
- **[Developer Guide](dev/README.md)** - Technical documentation for developers
- **[Threading Architecture](dev/threading-architecture.md)** - Multi-threaded job monitoring system
- **[Threading Troubleshooting](dev/threading-troubleshooting.md)** - Debug threading issues

## Documentation Overview

### Getting Started
1. **Read the [README](../README.md)** for installation and basic usage
2. **Set up [Configuration](configuration.md)** to eliminate repetitive options
3. **Create [Aliases](aliases.md)** for your common workflows
4. **Learn [Option Placement](option-placement.md)** to avoid command-line errors

### Configuration System
The configuration system allows you to:
- Set defaults for project, queue, resources, and paths
- Use template variables like `{user}`, `{project}`, `{timestamp}`
- Manage hierarchical configuration (CLI > User > System > Defaults)
- Validate and troubleshoot configuration files

**Key Commands:**
```bash
qxub config init                    # Create initial configuration
qxub config set defaults.project "a56"  # Set defaults
qxub config list                    # View current configuration
qxub config edit                    # Edit configuration file
```

### Alias System
Aliases enable ultra-simple workflow execution by combining:
- PBS job options (queue, resources, name)
- Environment setup (conda, modules, containers)  
- Commands to execute

**Key Commands:**
```bash
qxub config alias set myalias --env myenv --cmd "python script.py"
qxub alias myalias                  # Execute alias
qxub alias myalias --queue normal   # Execute with overrides
qxub config alias list              # List all aliases
```

## Common Workflows

### Data Science
```bash
# Set up environment
qxub config set defaults.project "ds01"

# Create analysis pipeline using new unified CLI format
qxub config alias set preprocess --env datasci --cmd "python preprocess.py"
qxub config alias set train --env pytorch --cmd "python train.py" --queue gpuvolta --resources "ngpus=1,ncpus=12"
qxub config alias set evaluate --env datasci --cmd "python evaluate.py"

# Execute pipeline
qxub alias preprocess data/
qxub alias train
qxub alias evaluate
```

### Bioinformatics
```bash
# Set up tools using unified CLI format
qxub config alias set qc --mods "fastqc,multiqc" --cmd "fastqc"
qxub config alias set align --mods "bwa,samtools" --cmd "bwa mem ref.fa"
qxub config alias set variants --sif "/containers/gatk.sif" --cmd "gatk HaplotypeCaller"

# Execute analysis
qxub alias qc reads.fastq.gz
qxub alias align reads_1.fastq.gz reads_2.fastq.gz  
qxub alias variants -I aligned.bam -R ref.fa
```

### Data Management
```bash
# Set up data operations using unified CLI format
qxub config alias set dvc_push --env dvc3 --cmd "dvc push" --queue copyq
qxub config alias set backup --mod rsync --cmd "rsync -av" --queue copyq
qxub config alias set sync_cloud --env tools --cmd "rclone sync"

# Execute operations
qxub alias dvc_push
qxub alias backup /data/ /backup/
qxub alias sync_cloud local/ remote:bucket/
```

## Advanced Topics

### Template Variables
Use dynamic variables in configuration and aliases:
- `{user}` - Current username
- `{project}` - Project code from config
- `{timestamp}` - Current timestamp (YYYYMMDD_HHMMSS)
- `{date}` - Current date (YYYYMMDD)
- `{time}` - Current time (HHMMSS)

```yaml
defaults:
  name: "job_{date}"
  out: "/scratch/{project}/{user}/logs/{timestamp}/out"
  joblog: "/home/{user}/logs/{name}_{time}.log"
```

### GPU Computing
For GPU workloads, always specify both GPU and CPU requirements:

```bash
qxub config alias set gpu_job \
  --subcommand conda \
  --env pytorch \
  --cmd "python train.py" \
  --queue gpuvolta \
  --resources "ngpus=1,ncpus=12,mem=32GB"
```

### Error Handling
All validation errors return exit code 2 for consistent script integration:
- Missing required parameters
- Invalid resource formats  
- Non-existent aliases or environments
- Option placement errors

## Troubleshooting

### Common Issues

**"No such option" errors:**
- Check option placement: main options (--name, --queue, --resources) before subcommand
- Use `qxub --help` and `qxub <subcommand> --help` for valid options

**GPU queue failures:**
- Ensure both GPU and CPU resources specified: `--resources ngpus=1,ncpus=12`
- Check queue availability and limits

**Alias not found:**
- Use `qxub config alias list` to see available aliases
- Check spelling and case sensitivity

**Template variable errors:**
- Ensure referenced variables exist in configuration
- Use `qxub config list` to check available values

### Getting Help

```bash
qxub --help                    # General help
qxub <subcommand> --help       # Subcommand-specific help
qxub config --help            # Configuration management
qxub config alias --help      # Alias management
qxub --dry-run <command>      # Preview execution
```

## File Locations

**Configuration files:**
- User config: `~/.config/qxub/config.yaml`
- System config: `/etc/xdg/qxub/config.yaml` (optional)

**Default log locations:**
- STDOUT: `$TMPDIR/qt/{timestamp}/out`
- STDERR: `$TMPDIR/qt/{timestamp}/err`
- Job logs: Configurable via `defaults.joblog`

## Contributing

Found a bug or want to contribute? Visit the [GitHub repository](https://github.com/swarbricklab/qsub_tools) to:
- Report issues
- Submit feature requests
- Contribute code improvements
- Add documentation

---

*For the most up-to-date information, always refer to the command-line help: `qxub --help`*