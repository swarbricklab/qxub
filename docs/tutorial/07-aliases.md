# Shortcuts and Aliases: Two Powerful Systems for Command Optimization

qxub provides two complementary systems for optimizing your workflow:

1. **Shortcuts** - Modern command-prefix based system with automatic detection
2. **Aliases** - Traditional configuration-based system for predefined workflows

This section covers both systems, when to use each, and how they work together.

## Understanding the Two Systems

### Shortcuts (New System)
- **Automatic detection**: `qxub exec -- python script.py` automatically uses 'python' shortcut
- **Command-prefix based**: Matches the first word of your command
- **JSON storage**: Stored in `~/.config/qxub/shortcuts.json`
- **Modern interface**: Managed with `qxub shortcut` commands

### Aliases (Legacy System)
- **Explicit invocation**: `qxub alias quick -- python script.py`
- **Configuration-based**: Defined in YAML config files
- **Hierarchical**: Support for complex nested configurations
- **Legacy interface**: Managed with `qxub config alias` commands

## Using Shortcuts (Recommended for New Workflows)

### View Available Shortcuts

```bash
qxub shortcut list
```

**Expected output:**
```
                             Available Shortcuts
┏━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Command Prefix ┃ Context              ┃ Command   ┃ Description            ┃
┡━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━┩
│ python         │ conda: base          │ (dynamic) │ Python with base conda │
│                │                      │           │ environment            │
│ echo           │ default              │ (dynamic) │ Simple echo test       │
│ gcc            │ default              │ (dynamic) │ GCC compiler           │
└────────────────┴──────────────────────┴───────────┴────────────────────────┘
```

### Using Shortcuts Automatically

Shortcuts work seamlessly - just run your command normally:

```bash
# These automatically detect and use shortcuts:
qxub exec -- python script.py        # Uses 'python' shortcut (conda: base)
qxub exec -- echo "Hello world"      # Uses 'echo' shortcut (default execution)
qxub exec -- gcc -o program source.c # Uses 'gcc' shortcut (default execution)

# See what shortcut would be used:
qxub exec --dry -- python script.py
```

### Creating Shortcuts

```bash
# Create a shortcut for PyTorch development
qxub shortcut set "pytorch" --env pytorch --queue gpuvolta --resource ngpus=1 --description "PyTorch GPU training"

# Create a shortcut for data analysis
qxub shortcut set "pandas" --env analytics --resource mem=16GB --description "Pandas data analysis"

# Use them automatically:
qxub exec -- pytorch train.py
qxub exec -- pandas analyze.py
```

## Using Aliases (Legacy System)

### View Available Aliases

```bash
qxub alias list
```

**Expected output:**
```
Available aliases:
  • quick: conda (env: base) - (requires command args)
  • bigmem: conda (env: scipy) - (requires command args)
  • gpu: conda (env: tensorflow) - (requires command args)
  • parallel: conda (env: mpi) - (requires command args)
```

### Using Aliases Explicitly

Aliases require explicit invocation with the `qxub alias` command:

```bash
# Use aliases with explicit commands:
qxub alias quick -- python script.py
qxub alias bigmem -- python memory_intensive.py
qxub alias gpu -- python train_model.py
qxub alias parallel -- mpirun -n 8 parallel_app
```

### Combining Aliases with Additional Options

You can override alias settings:

```bash
# Use alias but override specific settings
qxub alias quick --queue express -- python test.py
qxub alias bigmem --resource mem=64GB -- python bigger_job.py
qxub alias gpu --job-name my-training -- python train.py
```

**Expected output:**
```
🔧 Using alias 'py': --env dvc3 --resources mem=8GB,ncpus=2,walltime=1:00:00
🚀 Submitting job...
📋 Job submitted: 12345692.gadi-pbs (qx-20241017-151052)
...
```

## Creating Your Own Aliases (Legacy System)

### User-Level Aliases

Create aliases in your personal configuration:

```bash
# Create a user alias for machine learning
qxub config set aliases.ml.env "pytorch"
qxub config set aliases.ml.queue "gpuvolta"
qxub config set aliases.ml.resources '["mem=32GB", "ncpus=12", "ngpus=1"]'

# Create an alias for quick data exploration
qxub config set aliases.explore.env "base"
qxub config set aliases.explore.resources '["mem=4GB", "ncpus=1", "walltime=30:00"]'
```

### Using Custom Aliases

```bash
# Use your machine learning alias
qxub alias ml -- python train_model.py

# Use your exploration alias
qxub alias explore -- python explore_data.py
```

## When to Use Each System

### Use Shortcuts When:
- ✅ You want **automatic detection** based on command names
- ✅ You prefer **modern JSON storage** and management
- ✅ You want **command-prefix based workflows** (e.g., all `python` commands use the same context)
- ✅ You're starting new workflows and want the **latest features**

### Use Aliases When:
- ✅ You have **complex hierarchical configurations**
- ✅ You want **explicit control** over when optimizations are applied
- ✅ You're working with **existing alias-based workflows**
- ✅ You need **configuration-file based management**

## Advanced Workflows Combining Both Systems

You can use both systems together for maximum flexibility:

```bash
# Create shortcuts for automatic detection
qxub shortcut set "python" --env base
qxub shortcut set "jupyter" --env datascience --resource mem=8GB

# Create aliases for specific workflows
qxub config set aliases.gpu-train.env "pytorch"
qxub config set aliases.gpu-train.queue "gpuvolta"
qxub config set aliases.gpu-train.resources '["ngpus=1", "ncpus=12"]'

# Use shortcuts automatically:
qxub exec -- python analyze.py           # Auto-detects python shortcut
qxub exec -- jupyter notebook           # Auto-detects jupyter shortcut

# Use aliases explicitly for specific workflows:
qxub alias gpu-train -- python train_model.py
```

## Best Practices

### For Shortcuts:
- Use **descriptive command prefixes** that match your actual commands
- Create shortcuts for **commonly used command patterns**
- Keep shortcuts **simple and focused** on execution context

### For Aliases:
- Use aliases for **complex multi-step workflows**
- Leverage **configuration hierarchy** for team vs personal settings
- Use **meaningful alias names** that describe the workflow purpose

### General Guidelines:
- **Start with shortcuts** for new workflows - they're more intuitive
- **Keep aliases** for existing complex configurations
- **Document your shortcuts and aliases** for team members
- **Use overrides sparingly** - if you override often, create a new shortcut/alias

```bash
# Create different variations of data analysis
qxub config alias set analysis small --env dvc3 --resources mem=4GB,ncpus=2,walltime=1:00:00
qxub config alias set analysis large --env dvc3 --resources mem=64GB,ncpus=8,walltime=8:00:00 --queue hugemem
```

### Using Sub-aliases

```bash
# Small analysis
qxub analysis small -- python quick_analysis.py

# Large analysis (automatically uses hugemem queue)
qxub analysis large -- python big_analysis.py
```

### View Hierarchical Aliases

```bash
qxub config alias show analysis
```

**Expected output:**
```
📋 Alias: analysis
├── small: --env dvc3 --resources mem=4GB,ncpus=2,walltime=1:00:00
└── large: --env dvc3 --resources mem=64GB,ncpus=8,walltime=8:00:00 --queue hugemem
```

## Project-Level Aliases

For team projects, you can create project-specific aliases:

### Create Project Config Directory

```bash
# In your project directory
mkdir -p .qxub
```

### Create Project Aliases

```bash
# Create project-specific alias (in .qxub/config.yaml)
cat > .qxub/config.yaml << EOF
aliases:
  preprocess:
    main:
      env: dvc3
      resources: ["mem=16GB", "ncpus=4", "walltime=2:00:00"]

  train:
    quick:
      env: dvc3
      resources: ["mem=8GB", "ncpus=2", "walltime=1:00:00"]
    full:
      env: dvc3
      resources: ["mem=32GB", "ncpus=8", "walltime=8:00:00"]

  evaluate:
    main:
      env: dvc3
      resources: ["mem=4GB", "ncpus=1", "walltime=30:00"]
EOF
```

### Using Project Aliases

```bash
# These aliases are only available when running from this project directory
qxub preprocess -- python3 preprocess_data.py
qxub train quick -- python3 train.py --epochs 10
qxub train full -- python3 train.py --epochs 100
qxub evaluate -- python3 evaluate_model.py
```

## Alias Configuration Hierarchy

qxub resolves aliases from multiple sources, with this precedence:

1. **Command-line options** (highest priority)
2. **Project-level aliases** (`.qxub/config.yaml`)
3. **User-level aliases** (`~/.config/qxub/config.yaml`)
4. **System-level aliases** (`/g/data/a56/config/xdg/qxub/config.yaml`)

### Viewing Alias Origins

```bash
qxub config alias list --show-origin
```

**Expected output:**
```
📋 Available Aliases (with origins):

py (Python Data Science):
└── main: --env dvc3 --mem 8GB --ncpus 2 --walltime 1:00:00
    (/g/data/a56/config/xdg/qxub/config.yaml)

ml (Machine Learning):
└── main: --env dvc3 --mem 32GB --ncpus 8 --walltime 4:00:00
    (~/.config/qxub/config.yaml)

preprocess (Data Preprocessing):
└── main: --env dvc3 --mem 16GB --ncpus 4 --walltime 2:00:00
    (.qxub/config.yaml)
```

## Advanced Alias Patterns

### Environment-Specific Aliases

```bash
# Create aliases for different R environments
qxub config alias set seurat3 main --env seurat3 --mem 16GB --ncpus 4 --walltime 3:00:00
qxub config alias set seurat5 main --env seurat_5.1 --mem 16GB --ncpus 4 --walltime 3:00:00

# Use for different analyses
qxub seurat3 -- Rscript legacy_analysis.R
qxub seurat5 -- Rscript modern_analysis.R
```

### Resource-Scaled Aliases

```bash
# Create aliases for different data scales
qxub config alias set genomics small --env pysam --mem 8GB --ncpus 2 --walltime 2:00:00
qxub config alias set genomics medium --env pysam --mem 32GB --ncpus 8 --walltime 4:00:00
qxub config alias set genomics large --env pysam --mem 128GB --ncpus 16 --walltime 8:00:00 --queue hugemem

# Scale based on data size
qxub genomics small -- python3 analyze_sample.py
qxub genomics large -- python3 analyze_genome.py
```

### Pipeline Stage Aliases

```bash
# Create aliases for different pipeline stages
qxub config alias set pipeline qc --env dvc3 --mem 4GB --ncpus 2 --walltime 1:00:00
qxub config alias set pipeline align --env pysam --mem 16GB --ncpus 8 --walltime 4:00:00
qxub config alias set pipeline variant --env pysam --mem 32GB --ncpus 4 --walltime 6:00:00

# Run pipeline stages
qxub pipeline qc -- python3 quality_control.py
qxub pipeline align -- python3 alignment.py
qxub pipeline variant -- python3 variant_calling.py
```

## Managing Aliases

### List All Aliases

```bash
# Show all aliases from all sources
qxub config alias list --all

# Show only user aliases
qxub config alias list --user

# Show only system aliases
qxub config alias list --system
```

### Delete Aliases

```bash
# Delete a user alias
qxub config alias delete ml

# Delete a specific sub-alias
qxub config alias delete analysis large
```

### Modify Existing Aliases

```bash
# Update an existing alias
qxub config alias set ml main --env dvc3 --mem 48GB --ncpus 12 --walltime 6:00:00

# Add a new sub-alias to existing hierarchy
qxub config alias set analysis xlarge --mem 128GB --ncpus 16 --walltime 12:00:00 --queue hugemem
```

## Debugging Aliases

### Test Aliases with --dry

Always test aliases before using them:

```bash
qxub --dry py -- python3 test.py
```

**Expected output:**
```
🔧 Using alias 'py': --env dvc3 --mem 8GB --ncpus 2 --walltime 1:00:00
🔍 DRY RUN - Would submit the following job:

📋 Job Configuration:
├── Environment: dvc3 (conda)
├── Resources: mem=8GB, ncpus=2, walltime=1:00:00
...
```

### Verbose Alias Resolution

```bash
qxub -v py -- python3 test.py
```

This shows exactly how the alias was resolved and applied.

## Best Practices

### 1. Start with System Aliases

Use the built-in aliases (`py`, `r`, `sc`, etc.) before creating custom ones.

### 2. Create Meaningful Names

```bash
# ✅ Good: Descriptive names
qxub config alias set singlecell main --env sc --mem 16GB
qxub config alias set bioinformatics main --env pysam --mem 8GB

# ❌ Bad: Cryptic names
qxub config alias set x main --env sc --mem 16GB
qxub config alias set tmp main --env pysam --mem 8GB
```

### 3. Use Hierarchical Organization

```bash
# ✅ Good: Organized hierarchy
qxub config alias set ml train --env dvc3 --mem 32GB
qxub config alias set ml evaluate --env dvc3 --mem 8GB
qxub config alias set ml infer --env dvc3 --mem 4GB

# ❌ Bad: Flat namespace
qxub config alias set ml_train main --env dvc3 --mem 32GB
qxub config alias set ml_eval main --env dvc3 --mem 8GB
qxub config alias set ml_infer main --env dvc3 --mem 4GB
```

### 4. Share Project Aliases

Include `.qxub/config.yaml` in your project repository so team members can use the same aliases.

## Key Takeaways

1. **Two systems available**: Shortcuts (modern, automatic) and Aliases (legacy, explicit)
2. **Shortcuts for new workflows**: Use `qxub exec --` with automatic command detection
3. **Aliases for complex configurations**: Use `qxub alias` for explicit workflow management
4. **Both can coexist**: Use shortcuts for daily commands, aliases for special workflows
5. **Start simple**: Begin with built-in shortcuts and aliases before creating custom ones

## Migration Strategy

If you're transitioning from aliases to shortcuts:

```bash
# Old alias approach:
qxub alias python -- python script.py

# New shortcut approach:
qxub shortcut set "python" --env base
qxub exec -- python script.py  # Automatic detection
```

## Next Steps

Now that you understand both systems:
- **[Configuration](08-configuration.md)** - Understand the full configuration system
- **[Parallel Execution](09-parallel-execution.md)** - Use shortcuts and aliases in parallel job patterns

Both shortcuts and aliases are essential for efficient qxub usage. Shortcuts provide modern automatic detection, while aliases offer explicit control for complex workflows.

---

**💡 Pro Tips:**
- **Shortcuts**: Use `qxub shortcut list` to see available shortcuts
- **Aliases**: Use `qxub alias list` to see available aliases
- **Testing**: Use `--dry` with both systems to preview behavior
- **Overrides**: Both systems support command-line overrides
- **Documentation**: Use `qxub shortcut show` and `qxub config alias show` for details
- **Team sharing**: Shortcuts use JSON files, aliases use YAML config files
