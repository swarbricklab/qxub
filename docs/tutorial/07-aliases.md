# Aliases: Command Shortcuts and Workflow Optimization

Aliases are one of qxub's most powerful features, allowing you to save complex command patterns as simple shortcuts. This section covers using the built-in aliases, creating your own, and managing alias hierarchies.

## Understanding qxub Aliases

Aliases in qxub are pre-configured combinations of:
- Execution environments (`--env`, `--mod`, etc.)
- Resource specifications (`--resources mem=8GB,ncpus=2,walltime=1:00:00`)
- Queue selections
- Any other qxub options

They're organized hierarchically and can be overridden at different configuration levels.

## Using Built-in Aliases

The system configuration already includes several ready-made aliases. Let's explore them:

### View Available Aliases

```bash
qxub config alias list
```

**Expected output:**
```
📋 Available Aliases:

py (Python Data Science):
└── main: --env dvc3 --resources mem=8GB,ncpus=2,walltime=1:00:00

r (R Analysis):
└── main: --env tidyverse --resources mem=8GB,ncpus=2,walltime=1:00:00

sc (Single-cell Analysis):
└── main: --env sc --resources mem=16GB,ncpus=4,walltime=2:00:00

bio (Bioinformatics):
└── main: --env pysam --resources mem=8GB,ncpus=4,walltime=2:00:00

bigmem (High Memory Jobs):
└── main: --resources mem=64GB,ncpus=8,walltime=4:00:00 --queue hugemem

test (Quick Testing):
└── main: --resources mem=2GB,ncpus=1,walltime=0:15:00 --queue express

parallel (Parallel Processing):
└── main: --resources mem=4GB,ncpus=8,walltime=2:00:00
```

### Using Simple Aliases

Replace long resource specifications with short aliases:

```bash
# Instead of this long command:
qxub --env dvc3 --resources mem=8GB,ncpus=2,walltime=1:00:00 -- python3 my_analysis.py

# Use this simple alias:
qxub py -- python3 my_analysis.py
```

**Expected output:**
```
🔧 Using alias 'py': --env dvc3 --resources mem=8GB,ncpus=2,walltime=1:00:00
🚀 Submitting job...
📋 Job submitted: 12345692.gadi-pbs (qx-20241017-151052)
...
```

### More Alias Examples

```bash
# R analysis with appropriate resources
qxub r -- Rscript my_analysis.R

# Single-cell analysis with high memory
qxub sc -- python3 scanpy_analysis.py

# Quick test with express queue
qxub test -- python test_script.py

# High-memory job automatically uses hugemem queue
qxub bigmem -- python3 memory_intensive.py

# Parallel processing with 8 CPUs
qxub parallel -- python3 parallel_analysis.py
```

### Combining Aliases with Additional Options

You can override or extend alias settings:

```bash
# Use py alias but with more memory
qxub py --mem 16GB -- python3 big_analysis.py

# Use test alias but with different walltime
qxub test --walltime 30:00 -- python3 longer_test.py

# Use sc alias but with different queue
qxub sc --queue normal -- python3 sc_analysis.py
```

**The alias provides the base configuration, and your additional options override specific settings.**

## Creating Your Own Aliases

### User-Level Aliases

Create aliases in your personal configuration:

```bash
# Create a user alias for machine learning
qxub config alias set ml --env dvc3 --resources mem=32GB,ncpus=8,walltime=4:00:00

# Create an alias for quick data exploration
qxub config alias set explore --env dvc3 --resources mem=4GB,ncpus=1,walltime=30:00
```

### View Your Custom Aliases

```bash
qxub config alias show ml
```

**Expected output:**
```
📋 Alias: ml
└── main: --env dvc3 --mem 32GB --ncpus 8 --walltime 4:00:00

🔍 Origin: User configuration (~/.config/qxub/config.yaml)
```

### Using Custom Aliases

```bash
# Use your machine learning alias
qxub ml -- python3 train_model.py

# Use your exploration alias
qxub explore -- python3 explore_data.py
```

## Hierarchical Aliases

Aliases can have multiple sub-commands for different variations:

### Creating Sub-aliases

```bash
# Create different variations of data analysis
qxub config alias set analysis small --env dvc3 --mem 4GB --ncpus 2 --walltime 1:00:00
qxub config alias set analysis medium --env dvc3 --mem 16GB --ncpus 4 --walltime 4:00:00
qxub config alias set analysis large --env dvc3 --mem 64GB --ncpus 8 --walltime 8:00:00 --queue hugemem
```

### Using Sub-aliases

```bash
# Small analysis
qxub analysis small -- python3 quick_analysis.py

# Medium analysis
qxub analysis medium -- python3 standard_analysis.py

# Large analysis (automatically uses hugemem queue)
qxub analysis large -- python3 big_analysis.py
```

### View Hierarchical Aliases

```bash
qxub config alias show analysis
```

**Expected output:**
```
📋 Alias: analysis
├── small: --env dvc3 --mem 4GB --ncpus 2 --walltime 1:00:00
├── medium: --env dvc3 --mem 16GB --ncpus 4 --walltime 4:00:00
└── large: --env dvc3 --mem 64GB --ncpus 8 --walltime 8:00:00 --queue hugemem
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

1. **Built-in aliases**: Start with system-provided shortcuts (`py`, `r`, `sc`, etc.)
2. **Custom aliases**: Create shortcuts for your common workflows
3. **Hierarchical organization**: Use sub-aliases for variations
4. **Configuration hierarchy**: Project > User > System precedence
5. **Always test**: Use `--dry` to verify alias behavior

## Next Steps

Now that you understand aliases:
- **[Configuration](08-configuration.md)** - Understand the full configuration system
- **[Parallel Execution](09-parallel-execution.md)** - Use aliases in parallel job patterns

Aliases are essential for efficient qxub usage. They eliminate repetitive typing and ensure consistent resource allocation across your workflows.

---

**💡 Pro Tips:**
- Use `qxub config alias list --show-origin` to understand where aliases come from
- Create project-specific aliases in `.qxub/config.yaml` for team sharing
- Test new aliases with `--dry` before relying on them
- Use hierarchical aliases (e.g., `analysis small`, `analysis large`) for related workflows
- Override alias settings when needed: `qxub py --mem 16GB -- script.py`
