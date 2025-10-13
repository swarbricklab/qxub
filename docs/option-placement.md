# Option Placement

qxub supports two command specification methods:

## Method 1: Traditional `--` Separator

```bash
qxub [QXUB_OPTIONS] -- [YOUR_COMMAND]
```

## Method 2: `--cmd` Option

```bash
qxub [QXUB_OPTIONS] --cmd "YOUR_COMMAND"
```

Use `--cmd` for complex commands with variables, quotes, or special characters.

## Correct Examples

```bash
# ✅ Traditional separator - simple commands
qxub --env myenv --queue normal -l mem=16GB -- python script.py

# ✅ --cmd option - complex commands with variables
qxub --env myenv --cmd "python script.py --input ${HOME}/data.txt"
qxub --env myenv --cmd 'echo "Job ${{PBS_JOBID}} on ${{HOSTNAME}}"'

# ✅ Execution context options
qxub --env pytorch -- python train.py
qxub --mod python3 -- python analysis.py
qxub --sif container.sif -- blastn -query input.fa

# ✅ PBS options
qxub -l mem=32GB -l ncpus=8 --queue gpu -- python script.py
```

## Wrong Examples

```bash
# ❌ qxub options after --
qxub -- python script.py --env myenv

# ❌ Both --cmd and -- together
qxub --env myenv --cmd "echo hello" -- echo world

# ❌ Missing --
qxub --env myenv python script.py

# ❌ Mixed placement
qxub --env myenv python script.py --queue normal
```

## Why This Matters

The `--` separator tells qxub where its options end and your command begins. This prevents conflicts between qxub options and your command's arguments.

## Command Arguments

Everything after `--` is passed directly to your command:

```bash
qxub --env myenv -- python script.py --input data.csv --output results.json
```

Here `--input` and `--output` belong to your Python script, not qxub.
