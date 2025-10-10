# Option Placement Guide

Understanding where to place options in `qxub` commands is crucial for successful usage. With qxub 2.0's unified CLI, option placement is much simpler than before.

## The Golden Rule

**All qxub options must come BEFORE the `--` separator, all target command options come AFTER.**

```bash
qxub [OPTIONS] -- [COMMAND] [COMMAND_ARGUMENTS]
```

## Option Categories

### Execution Context Options (mutually exclusive)
Choose ONE execution context for your job:

| Option | Alternative Names | Description | Example |
|--------|-------------------|-------------|---------|
| `--env` | `--conda` | Conda environment | `--env myenv` |
| `--mod` | | Single module (repeatable) | `--mod python3 --mod gcc` |
| `--mods` | `--modules` | Comma-separated modules | `--mods python3,gcc` |
| `--sif` | `--sing`, `--singularity` | Singularity container | `--sif container.sif` |

### PBS Job Options
These control PBS job submission:

| Option | Description | Example |
|--------|-------------|---------|
| `--name` | Job name | `--name analysis_job` |
| `--queue` | PBS queue | `--queue normal` |
| `--project` | Project code | `--project a56` |
| `-l` | PBS resources (repeatable) | `-l mem=16GB -l ncpus=4` |
| `--out` | STDOUT log path | `--out /logs/output.log` |
| `--err` | STDERR log path | `--err /logs/error.log` |
| `--joblog` | Job log path | `--joblog /logs/job.log` |
| `--execdir` | Working directory | `--execdir /scratch/work/` |

### Control Options
These control qxub behavior:

| Option | Description | Example |
|--------|-------------|---------|
| `--dry` | Preview command only | `--dry` |
| `--quiet` | Submit and exit | `--quiet` |
| `--verbose` | Increase verbosity | `-v`, `-vv`, `-vvv` |

### Pre/Post Command Options
These add setup and cleanup commands:

| Option | Description | Example |
|--------|-------------|---------|
| `--pre` | Command before main command | `--pre "echo Starting"` |
| `--post` | Command after main command | `--post "echo Done"` |

## Correct Examples

### Basic Commands

```bash
# ✅ Correct: Conda environment execution
qxub --env myenv -- python script.py
qxub --conda myenv --queue normal --name myjob -- python script.py

# ✅ Correct: Module execution
qxub --mod samtools -- samtools --version
qxub --mods python3,gcc -l mem=16GB -- python analysis.py
qxub --mod python3 --mod gcc -- python script.py

# ✅ Correct: Singularity execution
qxub --sif container.sif -- python script.py
qxub --singularity /containers/blast.sif --bind /data --project bio01 -- blastn -query input.fa

# ✅ Correct: Default execution (no environment)
qxub -- echo "Hello world"
qxub -l mem=8GB -l walltime=01:00:00 -- ./my_program

# ✅ Correct: Control options
qxub --dry --env myenv -- python script.py
qxub --quiet --name background_job --env myenv -- python script.py
```

### Complex Commands

```bash
# ✅ Correct: Full option placement with pre/post commands
qxub --name "analysis_$(date +%Y%m%d)" --queue normal -l mem=16GB -l ncpus=4 \
     --env analysis --pre "echo Starting" --post "echo Done" -- \
     python analyze.py --input data.csv --output results.csv

# ✅ Correct: GPU job with multiple resources
qxub --project a56 --queue gpuvolta -l ngpus=1 -l ncpus=12 -l mem=32GB \
     --conda pytorch -- python train.py --epochs 100

# ✅ Correct: Container with binds and environment variables
qxub --project bio01 --queue express -l mem=32GB -l ncpus=8 \
     --sif /containers/analysis.sif --bind /data:/data -- \
     bash -c 'export DEBUG=1 && python pipeline.py'
```

### Alias Execution

```bash
# ✅ Correct: Basic alias execution
qxub alias myalias

# ✅ Correct: Global options before alias
qxub --dry alias myalias
qxub --quiet alias myalias

# ✅ Correct: Override options after alias name
qxub alias dvc_push --queue normal
qxub alias analysis --env different_env -l mem=16GB
qxub alias gpu_job --queue gpuvolta --cmd "python train.py --epochs 100"

# ✅ Correct: Combination of both
qxub --dry alias myalias --queue express
```

## Common Mistakes

### Wrong Option Placement

```bash
# ❌ Wrong: Execution context option after command
qxub -- python script.py --env myenv
#    ^^^^^^^^^^^^^^^^^^^^^^ command   ^^^^^^^^^^^^ qxub option (wrong place!)

# ❌ Wrong: PBS options after --
qxub --env myenv -- python script.py --queue normal
#    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ correct  ^^^^^^^^^^^^^^ qxub option (wrong place!)

# ❌ Wrong: Control flags mixed with command arguments
qxub --env myenv -- python script.py --dry
#    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ command  ^^^^^ qxub option (wrong place!)
```

### Corrected Versions

```bash
# ✅ Correct: All qxub options before --
qxub --env myenv -- python script.py
qxub --env myenv --queue normal -- python script.py
qxub --dry --env myenv -- python script.py
```

### Mutually Exclusive Execution Contexts

```bash
# ❌ Wrong: Multiple execution contexts
qxub --env myenv --mod python3 -- python script.py
qxub --conda pytorch --sif container.sif -- python script.py

# ✅ Correct: Single execution context
qxub --env myenv -- python script.py
qxub --mod python3 -- python script.py
qxub --sif container.sif -- python script.py
```

## Alternative Option Names

qxub 2.0 provides multiple option names for better usability:

```bash
# All equivalent conda options:
qxub --env myenv -- python script.py
qxub --conda myenv -- python script.py

# All equivalent module options:
qxub --mod python3 -- python script.py        # Single module
qxub --mods python3,gcc -- python script.py   # Comma-separated
qxub --modules python3,gcc -- python script.py # Alternative name

# All equivalent Singularity options:
qxub --sif container.sif -- python script.py
qxub --sing container.sif -- python script.py
qxub --singularity container.sif -- python script.py
```

## Special Case: Default Execution

When no execution context is specified, qxub uses default execution (direct PBS submission):

```bash
# ✅ Correct: Default execution patterns
qxub -- echo "Hello world"
qxub -l walltime=02:00:00 -- ./compiled_program
qxub --pre "mkdir -p /tmp/work" --post "rm -rf /tmp/work" -- python script.py
qxub --name data-processing -l mem=32GB -- ./process_data.sh
```

This is useful for:
- Running compiled programs that don't need specific environments
- Simple shell commands and scripts
- Jobs where you want to use system-installed tools
- Custom setups that don't fit conda/module/container patterns

### Why Alias is Special

## Why This Matters

The `--` separator rules exist because:

1. **Clear separation**: qxub options vs target command options
2. **Unambiguous parsing**: No confusion about which options belong where
3. **Consistent behavior**: Same rules apply regardless of execution context
4. **Prevents conflicts**: Target command options can't interfere with qxub

When options are in the wrong place:
- qxub may not recognize them (`No such option` errors)
- Options may be passed to the wrong component
- The generated PBS script may be incorrect

## Option Placement by Use Case

### Data Science

```bash
# Basic analysis
qxub --name analysis --project ds01 --env datasci -- python analyze.py

# GPU training
qxub --name training --queue gpuvolta -l ngpus=1 -l ncpus=12 -l mem=32GB \
     --conda pytorch -- python train.py --epochs 100 --batch-size 64

# With alias override
qxub --queue express alias gpu_training
```

### Bioinformatics

```bash
# Quality control
qxub --name qc -l mem=8GB -l ncpus=2 \
     --mods fastqc,multiqc -- fastqc reads.fastq.gz

# Alignment
qxub --name alignment -l mem=32GB -l ncpus=16 \
     --mods bwa,samtools -- bwa mem ref.fa reads1.fq reads2.fq

# Container analysis
qxub --name variants --project bio01 \
     --sif /containers/gatk.sif --bind /data:/data -- \
     gatk HaplotypeCaller -I input.bam -R ref.fa -O variants.vcf
```

### System Administration

```bash
# Data transfer
qxub --name backup --queue copyq -- rsync -av /data/ /backup/

# System maintenance
qxub -l walltime=02:00:00 --name cleanup \
     --pre "echo Starting cleanup" --post "echo Cleanup complete" -- \
     bash cleanup_script.sh

# Monitoring
qxub --quiet --name monitoring -- ./monitor_system.py --daemon
```

## Quick Reference

| Syntax | Description |
|--------|-------------|
| `qxub [OPTIONS] -- [COMMAND]` | Basic pattern |
| `qxub --env NAME -- [COMMAND]` | Conda execution |
| `qxub --mod MODULE -- [COMMAND]` | Single module |
| `qxub --mods MOD1,MOD2 -- [COMMAND]` | Multiple modules |
| `qxub --sif IMAGE -- [COMMAND]` | Singularity container |
| `qxub -- [COMMAND]` | Default execution |

**Remember**: All qxub options before `--`, all command options after `--`

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
