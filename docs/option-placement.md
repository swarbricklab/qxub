# Option Placement Guide

Understanding where to place options in `qxub` commands is crucial for successful usage. This guide explains the rules and provides clear examples.

## The Golden Rule

**Main options must come BEFORE subcommands, subcommand options come AFTER subcommands.**

```bash
qxub [MAIN_OPTIONS] SUBCOMMAND [SUBCOMMAND_OPTIONS] [ARGUMENTS]
```

## Main Options vs Subcommand Options

### Main Options (for qxub itself)
These control PBS job submission and must come **before** any subcommand:

| Option | Description | Example |
|--------|-------------|---------|
| `--name` | Job name | `--name analysis_job` |
| `--queue` | PBS queue | `--queue normal` |
| `--project` | Project code | `--project a56` |
| `--resources` | PBS resources | `--resources mem=16GB,ncpus=4` |
| `--out` | STDOUT log path | `--out /logs/output.log` |
| `--err` | STDERR log path | `--err /logs/error.log` |
| `--joblog` | Job log path | `--joblog /logs/job.log` |
| `--execdir` | Working directory | `--execdir /scratch/work/` |
| `--dry-run` | Preview command | `--dry-run` |
| `--quiet` | Submit and exit | `--quiet` |

### Subcommand Options (specific to each subcommand)
These configure the execution environment and come **after** the subcommand:

#### conda subcommand options:
- `--env` - Conda environment name
- `--pre` - Command to run before main command  
- `--post` - Command to run after main command

#### module subcommand options:
- `--mod` - Single module to load
- `--mods` - Comma-separated modules to load
- `--pre` - Command to run before main command
- `--post` - Command to run after main command

#### sing subcommand options:
- `--sif` - Singularity image file
- `--bind` - Bind mount paths
- `--env` - Environment variables
- `--pre` - Command to run before main command
- `--post` - Command to run after main command

## Correct Examples

### Basic Commands

```bash
# ✅ Correct: Main options before subcommand
qxub --name myjob --queue normal conda --env myenv -- python script.py
qxub --resources mem=16GB module --mod samtools -- samtools --version
qxub --project a56 --queue gpuvolta --resources ngpus=1,ncpus=12 conda --env pytorch -- python train.py

# ✅ Correct: Global flags immediately after qxub
qxub --dry-run conda --env myenv -- python script.py
qxub --quiet --name background_job conda --env myenv -- python script.py
```

### Alias Execution

**The `alias` subcommand is special** - it accepts both global qxub options and override options:

```bash
# ✅ Correct: Global options before alias subcommand
qxub --dry-run alias myalias
qxub --quiet alias myalias

# ✅ Correct: Override options after alias name  
qxub alias dvc_push --queue normal
qxub alias analysis --env different_env --resources mem=16GB
qxub alias gpu_job --queue gpuvolta --cmd "python train.py --epochs 100"

# ✅ Correct: Combination of both
qxub --dry-run alias myalias --queue express
```

The alias subcommand accepts override options because it needs to modify the original alias definition with your changes.

### Complex Commands

```bash
# ✅ Correct: Full option placement
qxub --name "analysis_$(date +%Y%m%d)" --queue normal --resources mem=16GB,ncpus=4 \
     conda --env analysis --pre "echo Starting" --post "echo Done" -- \
     python analyze.py --input data.csv --output results.csv

# ✅ Correct: Singularity with multiple options
qxub --project bio01 --queue express --resources mem=32GB,ncpus=8 \
     sing --sif /containers/analysis.sif --bind /data:/data --env DEBUG=1 -- \
     python pipeline.py
```

## Common Mistakes

### Wrong Option Placement

```bash
# ❌ Wrong: Main options after subcommand
qxub conda --env myenv --name myjob --queue normal -- python script.py
#    ^^^^^ conda subcommand    ^^^^^^^^^^^^^^^^^^^ main options (wrong place!)

# ❌ Wrong: Global flags at the end  
qxub conda --env myenv -- python script.py --dry-run
#    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ command  ^^^^^^^^^ global flag (wrong place!)

# ❌ Wrong: Global flags after alias options
qxub alias myalias --queue normal --dry-run
#    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ alias execution ^^^^^^^^^ global flag (wrong place!)
```

### Corrected Versions

```bash
# ✅ Correct: Fixed option placement
qxub --name myjob --queue normal conda --env myenv -- python script.py
qxub --dry-run conda --env myenv -- python script.py
qxub --dry-run alias myalias --queue normal
```

## Special Case: Alias Subcommand

The `alias` subcommand is unique because it accepts override options after the alias name. This allows you to modify any aspect of the original alias definition:

```bash
# ✅ Correct: Override any alias option
qxub alias myalias --queue normal        # Override queue
qxub alias myalias --resources mem=16GB  # Override resources  
qxub alias myalias --env different_env   # Override environment
qxub alias myalias --cmd "new command"   # Override command completely

# ✅ Correct: Multiple overrides
qxub alias gpu_training --queue express --resources mem=32GB --cmd "python train.py --debug"

# ✅ Correct: Global options still come first
qxub --dry-run alias myalias --queue normal
qxub --quiet alias background_job --resources mem=8GB
```

### Why Alias is Special

Unlike other subcommands that have fixed option sets, aliases need to accept any possible override because:

1. **Aliases can contain any combination of options** from their original definition
2. **Users need flexibility** to override any part of an alias at runtime  
3. **Aliases are shortcuts** - if you couldn't override options, you'd need separate aliases for every variation

### Alias Override Categories

The alias subcommand accepts these override types:

**Main qxub options:** `--name`, `--queue`, `--project`, `--resources`, `--out`, `--err`, `--joblog`

**Subcommand options:** `--env`, `--mod`, `--mods`, `--sif`, `--bind`, `--pre`, `--post`

**Command options:** `--cmd` (completely replace the command)

All of these come **after** the alias name in the command.

## Why This Matters

The option placement rules exist because:

1. **qxub needs main options first** to configure the PBS job submission
2. **Subcommands need their options** to configure the execution environment
3. **Command arguments come last** so they can be properly passed through

When options are in the wrong place:
- qxub may not recognize them (`No such option` errors)
- Options may be ignored silently
- The generated PBS script may be incorrect

## Option Placement by Use Case

### Data Science

```bash
# Basic analysis
qxub --name analysis --project ds01 conda --env datasci -- python analyze.py

# GPU training  
qxub --name training --queue gpuvolta --resources ngpus=1,ncpus=12,mem=32GB \
     conda --env pytorch -- python train.py

# With alias override
qxub --queue express alias gpu_training
```

### Bioinformatics

```bash
# Quality control
qxub --name qc --resources mem=8GB,ncpus=2 \
     module --mods "fastqc,multiqc" -- fastqc reads.fastq.gz

# Alignment
qxub --name alignment --resources mem=32GB,ncpus=16 \
     module --mods "bwa,samtools" -- bwa mem ref.fa reads1.fq reads2.fq

# Container analysis
qxub --name variants --project bio01 \
     sing --sif /containers/gatk.sif --bind /data:/data -- \
     gatk HaplotypeCaller -I input.bam -R ref.fa -O variants.vcf
```

### Data Management

```bash
# Data transfer
qxub --name backup --queue copyq \
     module --mod rsync -- rsync -av /data/ /backup/

# Cloud sync with alias
qxub --queue copyq alias cloud_sync

# DVC operations
qxub --name dvc_push --queue copyq \
     conda --env dvc3 -- dvc push
```

## Quick Reference

### Template for New Commands

```bash
qxub [--dry-run] [--quiet] \
     [--name NAME] [--queue QUEUE] [--project PROJECT] \
     [--resources RESOURCES] [--out OUT] [--err ERR] \
     SUBCOMMAND [SUBCOMMAND_OPTIONS] -- COMMAND [ARGS...]
```

### Template for Alias Execution

```bash
qxub [--dry-run] [--quiet] \
     [--name NAME] [--queue QUEUE] [--project PROJECT] \
     [--resources RESOURCES] \
     alias ALIAS_NAME [SUBCOMMAND_OPTION_OVERRIDES] [EXTRA_ARGS...]
```

## Troubleshooting

### "No such option" Error

If you see this error, check option placement:

```bash
# Error: qxub conda --env myenv --name myjob
# Solution: Move --name before conda
qxub --name myjob conda --env myenv
```

### Options Being Ignored

If options seem to be ignored:

1. Check they're in the right position (main vs subcommand)
2. Use `--dry-run` to see the generated command
3. Verify option names are correct (`--resources` not `--resource`)

### Getting Help

```bash
qxub --help                    # See all main options
qxub SUBCOMMAND --help         # See subcommand-specific options
qxub --dry-run [command]       # Preview what will be executed
```

Remember: **Main options before subcommands, subcommand options after subcommands!**