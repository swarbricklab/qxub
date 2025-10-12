# qxub

Submit PBS jobs with conda environments, modules, or containers.

## Install
```bash
git clone https://github.com/swarbricklab/qsub_tools.git
cd qsub_tools
pip install -e .
```

## Usage
```bash
# Basic syntax: qxub [options] -- command
qxub --env myenv -- python script.py        # Conda environment
qxub --mod python3 -- python script.py     # Environment module
qxub --sif container.sif -- python script.py  # Singularity container
qxub -- python script.py                    # Direct submission

# Add PBS options before --
qxub -l mem=16GB --queue normal --env myenv -- python script.py

# Preview without running
qxub --dry --env myenv -- python script.py
```

## Key Rules
- **All qxub options go before `--`**
- **Your command goes after `--`**
- **Use `--dry` to preview**

## Options
| Option | Description | Example |
|--------|-------------|---------|
| `--env` | Conda environment | `--env pytorch` |
| `--mod` | Environment module | `--mod python3` |
| `--sif` | Singularity container | `--sif container.sif` |
| `-l` | PBS resources | `-l mem=16GB` |
| `--queue` | PBS queue | `--queue normal` |
| `--name` | Job name | `--name myjob` |

## Configuration (Optional)
```bash
# Set defaults to avoid repetition
qxub config set defaults.project "a56"
qxub config set defaults.queue "normal"

# Create shortcuts
qxub config alias set gpu --env pytorch --queue gpu -- python train.py
qxub alias gpu  # Run the alias
```

## Help
```bash
qxub --help                    # General help
qxub config --help            # Configuration help
```
