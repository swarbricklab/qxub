# Shortcuts System

Shortcuts provide a fast way to execute common commands with predefined execution contexts and PBS settings. They use command-prefix matching to automatically detect and apply the right environment settings. Shortcuts can be created at user level (personal) or system level (team-wide).

## Quick Start

```bash
# List available shortcuts
qxub config shortcut list

# Show shortcuts with source information (system vs user)
qxub config shortcut list --show-origin

# Use a shortcut (automatic detection)
qxub exec -- python script.py    # Auto-detects 'python' shortcut

# Explicit shortcut usage
qxub exec --shortcut python -- script.py

# Built-in command aliases for faster access
qx --env myenv -- python script.py     # Short for 'qxub exec'
qxet "my-shortcut" --env base           # Short for 'qxub config shortcut set'
```

## Built-in Command Aliases

For convenience, qxub provides built-in aliases for common commands:

- **`qx`** - Short for `qxub exec` (fastest way to run commands)
- **`qxtat`** - Short for `qxub status` (check job status)
- **`qxet`** - Short for `qxub config shortcut set` (create shortcuts quickly)

```bash
# These are equivalent:
qxub exec --env pytorch -- python train.py
qx --env pytorch -- python train.py

# Quick shortcut creation:
qxub config shortcut set "python" --env base
qxet "python" --env base

# Check job status:
qxub status
qxtat
```

## Creating Shortcuts

### User Shortcuts (Default)

```bash
# Basic shortcut with conda environment
qxub config shortcut set "python" --env base --description "Python with base conda environment"

# GPU training shortcut
qxub config shortcut set "pytorch" \
  --env pytorch \
  --queue gpuvolta \
  --resources ngpus=1,ncpus=12 \
  --description "PyTorch GPU training"

# Compiler shortcut with modules
qxub config shortcut set "gcc" --mod gcc --description "GCC compiler with modules"

# Complex analysis pipeline
qxub config shortcut set "bigdata" \
  --env analytics \
  --queue hugemem \
  --resources mem=500GB,ncpus=16 \
  --description "Big data analysis pipeline"
```

### System-Wide Shortcuts (Admin)

```bash
# Create system-wide shortcuts for team use (requires admin permissions)
qxub config shortcut set "dvc data status" \
  --system \
  --env dvc3 \
  --resources mem=64GB,ncpus=16,storage=gdata/a56+gdata/px14+gdata/bc07 \
  --description "DVC data pipeline status check"

# System-wide GPU pipeline
qxub config shortcut set "team-ml-pipeline" \
  --system \
  --env pytorch \
  --queue gpuvolta \
  --resources ngpus=2,ncpus=24,mem=128GB \
  --description "Team ML training pipeline"
```

## Managing Shortcuts

```bash
# List all shortcuts in a table
qxub config shortcut list

# Show shortcuts with source information
qxub config shortcut list --show-origin

# Show details of a specific shortcut
qxub config shortcut show python

# Rename a shortcut (user shortcut)
qxub config shortcut rename python py3

# Rename a system shortcut (requires admin permissions)
qxub config shortcut rename python py3 --system

# Delete a user shortcut
qxub config shortcut delete python

# Delete a system shortcut (requires admin permissions)
qxub config shortcut delete python --system

# Delete with confirmation prompt bypass
qxub config shortcut delete python --yes
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
qxub exec --shortcut pytorch --resources mem=64GB -- train.py
```

## Storage and File Locations

Shortcuts are stored in JSON format in XDG-compliant directories:

- **User shortcuts**: `~/.config/qxub/shortcuts.json`
- **System shortcuts**: `/etc/xdg/qxub/shortcuts.json`

System shortcuts have precedence over user shortcuts when both exist for the same prefix.

```bash
# View shortcuts with their source (system vs user)
qxub config shortcut list --show-origin
```

## Example Workflows

### Data Science Pipeline

```bash
# Create shortcuts for different stages
qxub config shortcut set "preprocess" --env pandas --resources mem=16GB
qxub config shortcut set "train" --env pytorch --queue gpuvolta --resources ngpus=1
qxub config shortcut set "evaluate" --env scipy --resources mem=8GB

# Use them in your pipeline
qxub exec -- preprocess raw_data.csv
qxub exec -- train model.py --epochs 100
qxub exec -- evaluate results.pkl
```

### Bioinformatics Analysis

```bash
# Create domain-specific shortcuts
qxub config shortcut set "blast" --sif /containers/blast.sif --resources mem=32GB
qxub config shortcut set "gatk" --env biotools --queue normal --resources mem=64GB
qxub config shortcut set "samtools" --mod samtools --resources mem=16GB

# Run analysis steps
qxub exec -- blast -query input.fa -db nt
qxub exec -- gatk HaplotypeCaller -I sample.bam
qxub exec -- samtools view -b sample.sam
```

### Compilation and Development

```bash
# Compiler shortcuts for different languages
qxub config shortcut set "gcc" --mod gcc --description "GNU C Compiler"
qxub config shortcut set "nvcc" --mod cuda --queue gpuvolta --description "NVIDIA CUDA Compiler"
qxub config shortcut set "mpicc" --mod openmpi --resources ncpus=8 --description "MPI C Compiler"

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
qxub config shortcut set "gpu" --env pytorch --queue gpu --description "GPU training alias"

# Usage remains similar:
# Old: qxub alias gpu python train.py
# New: qxub exec -- gpu python train.py
```

## Troubleshooting

### Shortcut Not Found
```bash
# Check if shortcut exists
qxub config shortcut list | grep python

# Create if missing
qxub config shortcut set "python" --env base
```

### Command Not Matching
```bash
# Check the exact command prefix
qxub exec --dry -- your-command-here  # Shows what would be executed

# Show shortcut details to verify settings
qxub config shortcut show your-prefix
```

### Permission Issues
```bash
# System shortcuts require appropriate permissions
sudo qxub config shortcut set "team-python" --system --env base

# Check source of shortcuts to understand precedence
qxub config shortcut list --show-origin
```
