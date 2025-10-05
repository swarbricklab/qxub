# Alias Guide

Aliases are the most powerful feature of qxub 1.0, enabling ultra-simple execution of complex workflows. Instead of typing long commands repeatedly, create reusable aliases that encapsulate your entire workflow configuration.

## Quick Start

Think of aliases as shortcuts for your common workflows. Instead of remembering and typing complex commands with many options, you create an alias once and then just run it by name.

```bash
# Create a simple alias for a common task
qxub config alias set dvc_push \
  --subcommand conda \
  --cmd "dvc push" \
  --env dvc3 \
  --queue copyq

# Now just execute it
qxub alias dvc_push

# Override options when needed
qxub alias dvc_push --queue normal
```

## Managing Aliases

Working with aliases is straightforward using these commands:

### Listing and Inspection

```bash
# List all aliases
qxub config alias list

# Show detailed alias information
qxub config alias show myalias

# Test alias configuration (dry run)
qxub config alias test myalias
```

### Creating and Modifying Aliases

```bash
# Create a new alias
qxub config alias set myalias \
  --subcommand conda \
  --cmd "python script.py" \
  --env myenv \
  --name "job_{date}"

# Update existing alias
qxub config alias set existing_alias --queue newqueue --resources "mem=8GB"

# Delete alias
qxub config alias delete old_alias
```

### Executing Aliases

```bash
# Execute as defined
qxub alias myalias

# Execute with additional command arguments
qxub alias analysis input.txt output.txt

# Use global qxub options (these must come before 'alias')
qxub --dry-run alias myalias     # Preview execution
qxub --quiet alias myalias       # Run in quiet mode

# Override alias options (these come after alias name)
qxub alias myalias --queue normal --resources mem=16GB
qxub alias myalias --env different_env --mods "python3,gcc"
qxub alias myalias --cmd "python different_script.py"
```

> 📖 **Note**: The `alias` subcommand is special - it accepts override options after the alias name. 
> Global options like `--dry-run` still come before `alias`.
> See the **[Option Placement Guide](option-placement.md)** for complete details.

## Hierarchical Alias Structure

Version 1.0 introduces a hierarchical structure that prevents common CLI syntax errors by separating main qxub options from subcommand-specific options:

```yaml
aliases:
  my_workflow:
    main:                    # Options for qxub command level
      name: "workflow_{date}"
      queue: "normal"
      resources:
        - "mem=16GB"
        - "ncpus=4"
    subcommand:             # Subcommand-specific options
      type: conda           # or 'module', 'sing'
      env: "myenv"
    target:                 # Command to execute
      cmd: "python analysis.py"
```

This structure ensures proper option placement and prevents the common error of putting main options after subcommands.

## Creating Aliases

### Conda-Based Aliases

```bash
# Basic conda alias
qxub config alias set analysis \
  --subcommand conda \
  --cmd "python analyze.py" \
  --env myenv \
  --name "analysis_{date}"

# GPU-enabled conda alias (note proper resource specification)
qxub config alias set gpu_training \
  --subcommand conda \
  --cmd "python train.py" \
  --env pytorch \
  --queue gpuvolta \
  --resources "ngpus=1,ncpus=12,mem=32GB" \
  --name "training_{timestamp}"

# Data processing alias
qxub config alias set preprocess \
  --subcommand conda \
  --cmd "python preprocess.py" \
  --env datasci \
  --resources "mem=8GB,ncpus=2" \
  --name "preprocess_{date}"
```

### Module-Based Aliases

```bash
# Single module
qxub config alias set samtools_analysis \
  --subcommand module \
  --cmd "samtools view -c" \
  --mod samtools \
  --name "sam_analysis_{timestamp}"

# Multiple modules (comma-separated)
qxub config alias set bioinformatics \
  --subcommand module \
  --cmd "python pipeline.py" \
  --mods "python3,samtools,gcc" \
  --resources "mem=16GB,ncpus=4" \
  --name "bio_pipeline_{date}"

# R analysis
qxub config alias set r_analysis \
  --subcommand module \
  --cmd "Rscript analysis.R" \
  --mods "R,gcc" \
  --resources "mem=12GB,ncpus=2" \
  --name "r_analysis_{time}"
```

### Singularity Container Aliases

```bash
# Basic container alias
qxub config alias set container_analysis \
  --subcommand sing \
  --cmd "python analysis.py" \
  --sif "/containers/analysis.sif" \
  --bind "/data:/data" \
  --name "container_{timestamp}"

# Jupyter notebook in container
qxub config alias set jupyter \
  --subcommand sing \
  --cmd "jupyter lab --no-browser --ip=0.0.0.0" \
  --sif "/containers/jupyter.sif" \
  --bind "/home:/home,/data:/data" \
  --name "jupyter_{user}_{time}"

# Complex container workflow
qxub config alias set rnaseq \
  --subcommand sing \
  --cmd "nextflow run nf-core/rnaseq" \
  --sif "/containers/nextflow.sif" \
  --bind "/data:/data,/scratch:/scratch" \
  --env "THREADS=8" \
  --resources "mem=64GB,ncpus=16" \
  --queue express \
  --name "rnaseq_{date}"
```

### Command Appending

If an alias defines a base command, you can append additional arguments:

```bash
# Alias defines: cmd: "samtools view -c"
qxub alias count_reads input.bam        # Executes: samtools view -c input.bam
qxub alias count_reads *.bam           # Executes: samtools view -c *.bam

# Alias defines: cmd: "python analysis.py"  
qxub alias analyze --input data.csv    # Executes: python analysis.py --input data.csv
```

## Advanced Alias Features

### Template Variables in Aliases

Aliases support all template variables available in configuration:

```yaml
aliases:
  backup:
    main:
      name: "backup_{timestamp}"
      queue: "copyq"
    subcommand:
      type: conda
      env: "tools"
    target:
      cmd: "rsync -av data/ /backup/{user}/{date}/"

  analysis:
    main:
      name: "{user}_analysis_{date}"
      project: "{project}"  # Uses project from config
      out: "/scratch/{project}/{user}/analysis/{timestamp}/out"
    subcommand:
      type: conda
      env: "analysis"
    target:
      cmd: "python analyze.py --output results_{timestamp}.csv"
```

### Conditional Resource Allocation

Create aliases for different resource requirements:

```bash
# Light processing
qxub config alias set light_task \
  --subcommand conda \
  --cmd "python light_script.py" \
  --env myenv \
  --resources "mem=2GB,ncpus=1" \
  --name "light_{time}"

# Heavy processing  
qxub config alias set heavy_task \
  --subcommand conda \
  --cmd "python heavy_script.py" \
  --env myenv \
  --resources "mem=32GB,ncpus=8" \
  --name "heavy_{timestamp}"

# GPU processing
qxub config alias set gpu_task \
  --subcommand conda \
  --cmd "python gpu_script.py" \
  --env pytorch \
  --resources "mem=16GB,ncpus=12,ngpus=1" \
  --queue gpuvolta \
  --name "gpu_{date}"
```

## Real-World Examples

### Data Science Workflow

```bash
# Set up data science aliases
qxub config alias set preprocess \
  --subcommand conda \
  --env datasci \
  --cmd "python preprocess.py" \
  --resources "mem=8GB,ncpus=2" \
  --name "preprocess_{date}"

qxub config alias set train \
  --subcommand conda \
  --env pytorch \
  --cmd "python train.py" \
  --resources "mem=32GB,ncpus=12,ngpus=1" \
  --queue gpuvolta \
  --name "training_{timestamp}"

qxub config alias set evaluate \
  --subcommand conda \
  --env datasci \
  --cmd "python evaluate.py" \
  --resources "mem=4GB,ncpus=1" \
  --name "eval_{date}"

# Execute pipeline
qxub alias preprocess data/raw/
qxub alias train --cmd "python train.py --epochs 100 --lr 0.001"
qxub alias evaluate models/best_model.pt
```

### Bioinformatics Pipeline

```bash
# Quality control
qxub config alias set qc \
  --subcommand module \
  --cmd "fastqc" \
  --mods "fastqc,multiqc" \
  --resources "mem=4GB,ncpus=2" \
  --name "qc_{date}"

# Alignment
qxub config alias set align \
  --subcommand module \
  --cmd "bwa mem ref.fa" \
  --mods "bwa,samtools" \
  --resources "mem=32GB,ncpus=16" \
  --name "align_{timestamp}"

# Variant calling
qxub config alias set variants \
  --subcommand sing \
  --cmd "gatk HaplotypeCaller" \
  --sif "/containers/gatk.sif" \
  --bind "/data:/data" \
  --resources "mem=16GB,ncpus=4" \
  --name "variants_{date}"

# Execute analysis
qxub alias qc reads.fastq.gz
qxub alias align reads_1.fastq.gz reads_2.fastq.gz
qxub alias variants -I aligned.bam -R ref.fa -O variants.vcf
```

### Data Management

```bash
# DVC operations
qxub config alias set dvc_push \
  --subcommand conda \
  --env dvc3 \
  --cmd "dvc push" \
  --queue copyq \
  --name "push_{time}"

qxub config alias set dvc_pull \
  --subcommand conda \
  --env dvc3 \
  --cmd "dvc pull" \
  --queue copyq \
  --name "pull_{time}"

# Backup operations
qxub config alias set backup \
  --subcommand module \
  --mod rsync \
  --cmd "rsync -av --progress" \
  --queue copyq \
  --name "backup_{date}"

# Cloud sync
qxub config alias set sync_cloud \
  --subcommand conda \
  --env tools \
  --cmd "rclone sync --progress" \
  --queue copyq \
  --name "sync_{timestamp}"

# Execute data operations
qxub alias dvc_push
qxub alias backup /home/user/data/ /backup/user/
qxub alias sync_cloud local/ remote:bucket/data/
```

## Important Notes

### GPU Queue Requirements

When creating GPU aliases, always specify both GPU and CPU requirements:

```bash
# ✅ Correct: Specify both GPU and CPU
qxub config alias set gpu_job \
  --subcommand conda \
  --env pytorch \
  --cmd "python train.py" \
  --queue gpuvolta \
  --resources "ngpus=1,ncpus=12,mem=32GB"

# ❌ Incorrect: Missing CPU requirement
qxub config alias set gpu_job \
  --queue gpuvolta \
  --resources "ngpus=1"  # Will fail - gpuvolta requires minimum 12 CPUs
```

### Option Placement

The hierarchical structure ensures options go to the right place automatically, but when overriding, remember the rules:

```bash
# ✅ Correct: Main options can be overridden
qxub alias myalias --queue normal --resources "mem=8GB"

# ✅ Correct: Subcommand options can be overridden  
qxub alias myalias --env newenv --mods "python3,gcc"

# ✅ Correct: Command can be completely overridden
qxub alias myalias --cmd "python different_script.py"
```

## Troubleshooting

### Common Issues

**Alias not found:**
```bash
qxub config alias list  # Check available aliases
```

**Wrong option placement:**
- The hierarchical structure should prevent this, but if you see "No such option" errors, check that you're using valid override options

**Resource validation:**
```bash
qxub config alias test myalias  # Test alias configuration
qxub --dry-run alias myalias    # Preview execution
```

**Template variable errors:**
- Ensure all referenced variables exist in your configuration
- Use `qxub config list` to check available template values

### Getting Help

```bash
qxub config alias --help       # Alias management help
qxub alias --help             # Alias execution help
qxub config alias test myalias # Validate specific alias
```