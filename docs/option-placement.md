# Option Placement

qxub uses `--` to separate qxub options from your command.

## The Rule

```bash
qxub [QXUB_OPTIONS] -- [YOUR_COMMAND]
```

## Correct Examples

```bash
# ✅ All qxub options before --
qxub --env myenv --queue normal -l mem=16GB -- python script.py

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
