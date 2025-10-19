# Option Placement

qxub exec supports two comma## Wrong Examples

```bash
# ❌ qxub options after --
qxub exec -- python script.py --env myenv

# ❌ Both --cmd and -- together
qxub exec --env myenv --cmd "echo hello" -- echo world

# ❌ Missing -- (old CLI syntax)
qxub --env myenv python script.py

# ❌ Mixed placement
qxub --env myenv python script.py --queue normal
```n methods:

## Method 1: Traditional `--` Separator

```bash
qxub exec [OPTIONS] -- [YOUR_COMMAND]
```

## Method 2: `--cmd` Option

```bash
qxub exec [OPTIONS] --cmd "YOUR_COMMAND"
```

Use `--cmd` for complex commands with variables, quotes, or special characters.

## Correct Examples

```bash
# ✅ Traditional separator - simple commands
qxub exec --env myenv --queue normal -l mem=16GB -- python script.py

# ✅ --cmd option - complex commands with variables
qxub exec --env myenv --cmd "python script.py --input ${HOME}/data.txt"
qxub exec --env myenv --cmd 'echo "Job ${{PBS_JOBID}} on ${{HOSTNAME}}"'

# ✅ Smart quotes - clean handling of nested quotes
qxub exec --env myenv --cmd "find /data -exec echo \"Found: {}\" \;"
qxub exec --env myenv --cmd "sh -c \"echo \\\"User: ${USER}\\\"\""

# ✅ Execution context options
qxub exec --env pytorch -- python train.py
qxub exec --mod python3 -- python analysis.py
qxub exec --sif container.sif -- blastn -query input.fa

# ✅ Shortcuts - automatic detection
qxub exec -- python train.py    # Auto-detects 'python' shortcut
qxub exec -- gcc --version      # Auto-detects 'gcc' shortcut

# ✅ Explicit shortcuts
qxub exec --shortcut python -- script.py

# ✅ PBS options
qxub exec -l mem=32GB -l ncpus=8 --queue auto -- python script.py
```

## Wrong Examples

```bash
# ❌ qxub options after --
qxub exec -- python script.py --env myenv

# ❌ Both --cmd and -- together
qxub exec --env myenv --cmd "echo hello" -- echo world

# ❌ Missing --
qxub exec --env myenv python script.py

# ❌ Mixed placement
qxub exec --env myenv python script.py --queue normal
```

## Why This Matters

The `--` separator tells qxub where its options end and your command begins. This prevents conflicts between qxub options and your command's arguments.

## Command Arguments

Everything after `--` is passed directly to your command:

```bash
qxub exec --env myenv -- python script.py --input data.csv --output results.json
```

Here `--input` and `--output` belong to your Python script, not qxub.
