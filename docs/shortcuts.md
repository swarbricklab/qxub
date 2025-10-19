# Shortcuts System

Shortcuts provide a fast way to execute common commands with predefined execution contexts and PBS settings. They use command-prefix matching to automatically detect and apply the right environment settings.

## Quick Start

```bash
# List available shortcuts
qxub shortcut list

# Use a shortcut (automatic detection)
qxub exec -- python script.py    # Auto-detects 'python' shortcut

# Explicit shortcut usage
qxub exec --shortcut python -- script.py
```

## Creating Shortcuts

```bash
# Basic shortcut with conda environment
qxub shortcut set "python" --env base --description "Python with base conda environment"

# GPU training shortcut
qxub shortcut set "pytorch" \
  --env pytorch \
  --queue gpuvolta \
  --resource ngpus=1 \
  --resource ncpus=12 \
  --description "PyTorch GPU training"

# Compiler shortcut with modules
qxub shortcut set "gcc" --mod gcc --description "GCC compiler with modules"

# Complex analysis pipeline
qxub shortcut set "bigdata" \
  --env analytics \
  --queue hugemem \
  --resource mem=500GB \
  --resource ncpus=16 \
  --description "Big data analysis pipeline"
```

## Managing Shortcuts

```bash
# List all shortcuts in a table
qxub shortcut list

# Show details of a specific shortcut
qxub shortcut show python

# Rename a shortcut
qxub shortcut rename python py3

# Delete a shortcut
qxub shortcut delete python

# Refresh shortcuts cache
qxub shortcut refresh
```

## How Shortcuts Work

### Automatic Detection

When you run `qxub exec -- COMMAND`, the system:

1. Extracts the first word of COMMAND as the prefix
2. Checks if a shortcut exists for that prefix
3. If found, applies the shortcut's settings automatically
4. Executes the command with the shortcut's context

```bash
# These commands are equivalent:
qxub exec -- python script.py                    # Auto-detection
qxub exec --shortcut python -- script.py         # Explicit
```

### Command Prefix Matching

Shortcuts match based on the first word of your command:

```bash
# Command: "python script.py" → prefix: "python"
qxub exec -- python script.py

# Command: "gcc -o program source.c" → prefix: "gcc"
qxub exec -- gcc -o program source.c

# Command: "pytorch-train --epochs 100" → prefix: "pytorch-train"
qxub exec -- pytorch-train --epochs 100
```

### Overriding Shortcut Settings

You can override any shortcut setting from the command line:

```bash
# Use python shortcut but override job name and queue
qxub exec --shortcut python --job-name my-analysis --queue express -- script.py

# Use pytorch shortcut but with different memory
qxub exec --shortcut pytorch --resource mem=64GB -- train.py
```

## Storage and File Locations

Shortcuts are stored in JSON format in XDG-compliant directories:

- **User shortcuts**: `~/.config/qxub/shortcuts.json`
- **System shortcuts**: `/etc/xdg/qxub/shortcuts.json`

```bash
# View shortcuts file locations
qxub shortcut files
```

## Example Workflows

### Data Science Pipeline

```bash
# Create shortcuts for different stages
qxub shortcut set "preprocess" --env pandas --resource mem=16GB
qxub shortcut set "train" --env pytorch --queue gpuvolta --resource ngpus=1
qxub shortcut set "evaluate" --env scipy --resource mem=8GB

# Use them in your pipeline
qxub exec -- preprocess raw_data.csv
qxub exec -- train model.py --epochs 100
qxub exec -- evaluate results.pkl
```

### Bioinformatics Analysis

```bash
# Create domain-specific shortcuts
qxub shortcut set "blast" --sif /containers/blast.sif --resource mem=32GB
qxub shortcut set "gatk" --env biotools --queue normal --resource mem=64GB
qxub shortcut set "samtools" --mod samtools --resource mem=16GB

# Run analysis steps
qxub exec -- blast -query input.fa -db nt
qxub exec -- gatk HaplotypeCaller -I sample.bam
qxub exec -- samtools view -b sample.sam
```

### Compilation and Development

```bash
# Compiler shortcuts for different languages
qxub shortcut set "gcc" --mod gcc --description "GNU C Compiler"
qxub shortcut set "nvcc" --mod cuda --queue gpuvolta --description "NVIDIA CUDA Compiler"
qxub shortcut set "mpicc" --mod openmpi --resource ncpus=8 --description "MPI C Compiler"

# Build projects
qxub exec -- gcc -O3 -o program source.c
qxub exec -- nvcc -o cuda_program cuda_source.cu
qxub exec -- mpicc -o mpi_program mpi_source.c
```

## Advanced Usage

### Dynamic Command Construction

Shortcuts preserve the full command after the prefix:

```bash
# Shortcut: python → conda env: base
qxub exec -- python -c "import sys; print(sys.version)"
qxub exec -- python script.py arg1 arg2 --verbose
qxub exec -- python -m pip install package
```

### Integration with PBS Options

Shortcuts work seamlessly with PBS resource specifications:

```bash
# Apply shortcut context + additional PBS options
qxub exec --shortcut python -l walltime=02:00:00 --job-name long-job -- train.py

# Override shortcut queue with auto-selection
qxub exec --shortcut pytorch --queue auto -- distributed_training.py
```

### Backup and Sharing

```bash
# Copy shortcuts to share with team members
cp ~/.config/qxub/shortcuts.json team_shortcuts.json

# Install team shortcuts
cp team_shortcuts.json ~/.config/qxub/shortcuts.json
qxub shortcut refresh
```

## Migration from Aliases

If you have existing aliases in your qxub configuration, you can create equivalent shortcuts:

```bash
# Old alias: qxub config alias set gpu --env pytorch --queue gpu
# New shortcut equivalent:
qxub shortcut set "gpu" --env pytorch --queue gpu --description "GPU training alias"

# Usage remains similar:
# Old: qxub alias gpu python train.py
# New: qxub exec -- gpu python train.py
```

## Troubleshooting

### Shortcut Not Found
```bash
# Check if shortcut exists
qxub shortcut list | grep python

# Create if missing
qxub shortcut set "python" --env base
```

### Command Not Matching
```bash
# Check the exact command prefix
qxub exec --dry -- your-command-here  # Shows what would be executed

# Show shortcut details to verify settings
qxub shortcut show your-prefix
```

### Cache Issues
```bash
# Refresh shortcuts cache if changes aren't appearing
qxub shortcut refresh
```
