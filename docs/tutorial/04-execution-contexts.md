# Execution Contexts: Running in Different Software Environments

One of qxub's most powerful features is the ability to easily run jobs in different software environments. Whether you need specific Python packages, R libraries, or specialized tools, qxub handles the environment setup automatically.

## Understanding Execution Contexts

qxub supports four execution contexts (but only one per job):

- **Default**: Use the login node environment
- **`--env`**: Activate a conda environment
- **`--mod` / `--mods`**: Load environment modules
- **`--sif`**: Run in a Singularity container

**Important**: You can only specify **one** execution context per job. Mixing them will result in an error.

## Conda Environments with `--env`

### Basic Conda Environment Usage

Let's run a Python job in the `dvc3` environment:

```bash
qxub exec --dry --env dvc3 -- python check_environment.py
```

**Expected dry run output:**
```
� Job command constructed
�📝 Command to execute: dvc doctor
Dry run - job would be submitted (use -v to see full qsub command)
```

To see more details about what PBS script would be generated, use `-v`:

```bash
qxub -v --dry --env dvc3 -- dvc doctor
```

**Verbose dry run output:**
```
🔧 Job command constructed
📝 Command to execute: dvc doctor
🔧 Full qsub command: qsub -v cmd_b64="ZHZjIGRvY3Rvcg==",cwd=/g/data/a56/software/qsub_tools,out=/scratch/a56/jr9959/qt/20251018_181306/out,err=/scratch/a56/jr9959/qt/20251018_181306/err,quiet=false,env="dvc3" -N qt -q normal -P a56 -o qt.log /g/data/a56/software/qsub_tools/qxub/jobscripts/qconda.pbs
```

Notice how qxub automatically uses the `qconda.pbs` job script template and passes the conda environment name (`env="dvc3"`) to activate the specified environment.

### Running the Conda Job

```bash
qxub exec --env dvc3 -- dvc doctor
```

**Expected output:**
```
� Job command constructed
✅ Job submitted successfully! Job ID: 152762000.gadi-pbs
🚀 Job started running
DVC version: 3.63.0 (conda)
---------------------------
Platform: Python 3.10.0 on Linux-4.18.0-553.62.1.el8.nci.x86_64-x86_64-with-glibc2.28
Subprojects:
        dvc_data = 3.16.4
        dvc_objects = 5.1.0
        dvc_render = 1.0.1
        dvc_task = 0.3.0
        scmrepo = 3.5.2
Supports:
        gdrive (pydrive2 = 1.21.3),
        gs (gcsfs = 2024.2.0),
        http (aiohttp = 3.9.1, aiohttp-retry = 2.8.3),
        https (aiohttp = 3.9.1, aiohttp-retry = 2.8.3),
        s3 (s3fs = 2024.2.0, boto3 = 1.34.34),
        ssh (sshfs = 2023.10.0)
Config:
        Global: /home/913/jr9959/.config/dvc
        System: /g/data/a56/config/xdg/dvc
✅ Command completed successfully
🎉 Job completed successfully
```

### Data Science Example with Resources

```bash
# Combine conda environment with resource requirements
qxub exec --env dvc3 --resources mem=8GB,ncpus=2 -- python data_analysis.py
```

### Other Environment Examples

```bash
# R analysis environment
qxub exec --env tidyverse -- Rscript analysis.R

# Single-cell analysis with high memory
qxub exec --env sc --resources mem=16GB -- python scanpy_analysis.py
```

## Environment Modules with `--mod` and `--mods`

### Single Module (`--mod`)

Load a single environment module:

```bash
qxub exec --dry --mod python3/3.11.7 -- python --version
```

**Expected dry run output:**
```
🔧 Job command constructed
📝 Command to execute: python --version
Dry run - job would be submitted (use -v to see full qsub command)
```

To see more details, use `-v`:

```bash
qxub -v --dry --mod python3/3.11.7 -- python --version
```

**Verbose dry run output:**
```
� Job command constructed
📝 Command to execute: python --version
🔧 Full qsub command: qsub -v cmd_b64="cHl0aG9uIC0tdmVyc2lvbg==",cwd=/g/data/a56/software/qsub_tools,out=/scratch/a56/jr9959/qt/20251018_181557/out,err=/scratch/a56/jr9959/qt/20251018_181557/err,quiet=false,mods="python3/3.11.7" -N qt -q normal -P a56 -o qt.log /g/data/a56/software/qsub_tools/qxub/jobscripts/qmod.pbs
```

Notice how qxub uses the `qmod.pbs` job script template and passes the module list (`mods="python3/3.11.7"`) to load the specified modules.

### Multiple Modules (`--mods`)

Load multiple modules at once:

```bash
qxub exec --dry --mods python3/3.11.7,gcc/11.1.0 -- python --version
```

**Expected dry run output:**
```
� Job command constructed
📝 Command to execute: python --version
Dry run - job would be submitted (use -v to see full qsub command)
```

To see the full command with multiple modules:

```bash
qxub -v --dry --mods python3/3.11.7,gcc/11.1.0 -- python --version
```

**Verbose dry run output:**
```
🔧 Job command constructed
📝 Command to execute: python --version
🔧 Full qsub command: qsub -v cmd_b64="cHl0aG9uIC0tdmVyc2lvbg==",cwd=/g/data/a56/software/qsub_tools,out=/scratch/a56/jr9959/qt/20251018_181632/out,err=/scratch/a56/jr9959/qt/20251018_181632/err,quiet=false,mods="python3/3.11.7 gcc/11.1.0" -N qt -q normal -P a56 -o qt.log /g/data/a56/software/qsub_tools/qxub/jobscripts/qmod.pbs
```

Notice how the modules are passed as a space-separated list: `mods="python3/3.11.7 gcc/11.1.0"`.

### Real Module Example

```bash
# Load Python module and run script
qxub exec --mod python3/3.11.7 -- python computation.py
```

### Combining Modules for Compilation

```bash
# Example: compiling with GCC (if you had source code)
qxub exec --mods gcc/11.1.0,python3/3.11.7 -- ./compile_project.sh
```

## Singularity Containers with `--sif`

qxub supports Singularity containers for reproducible environments:

```bash
qxub exec --dry --sif /g/data/a56/containers/example.sif -- echo "Hello from container"
```

**Expected dry run output:**
```
� Job command constructed
📝 Command to execute: echo Hello from container
Dry run - job would be submitted (use -v to see full qsub command)
```

To see the full command with Singularity:

```bash
qxub -v --dry --sif /g/data/a56/containers/example.sif -- echo "Hello from container"
```

**Verbose dry run output:**
```
🔧 Job command constructed
📝 Command to execute: echo Hello from container
🔧 Full qsub command: qsub -v cmd_b64="ZWNobyBIZWxsbyBmcm9tIGNvbnRhaW5lcg==",cwd=/g/data/a56/software/qsub_tools,out=/scratch/a56/jr9959/qt/20251018_181659/out,err=/scratch/a56/jr9959/qt/20251018_181659/err,quiet=false,sif="/g/data/a56/containers/example.sif" -N qt -q normal -P a56 -o qt.log /g/data/a56/software/qsub_tools/qxub/jobscripts/qsing.pbs
```

Notice how qxub uses the `qsing.pbs` job script template and passes the container path (`sif="/g/data/a56/containers/example.sif"`) to run commands inside the container.

## Execution Context Rules and Errors

### Cannot Mix Contexts

```bash
# This will fail - mixing conda and modules
qxub exec --env dvc3 --mod python3/3.11.7 -- echo "This fails"
```

**Expected error:**
```
❌ Error: Cannot specify multiple execution contexts
💡 Choose one: --env, --mod/--mods, or --sif
💡 Current contexts: conda environment 'dvc3', modules ['python3/3.11.7']
```

### Valid Single Context Examples

```bash
# ✅ Good: Conda only
qxub exec --env dvc3 -- python3 script.py

# ✅ Good: Single module only
qxub exec --mod python3/3.11.7 -- python3 script.py

# ✅ Good: Multiple modules only
qxub exec --mods python3/3.11.7,gcc/11.1.0 -- python3 script.py

# ✅ Good: Container only
qxub exec --sif container.sif -- python3 script.py

# ✅ Good: No execution context (uses login environment)
qxub exec --default -- python3 script.py
```

## Choosing the Right Execution Context

### Use `--env` when:
- You need specific Python/R packages
- Working with data science or bioinformatics
- You want a complete, curated software stack
- You need consistent package versions

### Use `--mod`/`--mods` when:
- You need system-optimized software
- Working with compiled applications
- You need specific compiler versions
- You want minimal overhead

### Use `--sif` when:
- You need perfect reproducibility
- Working with complex software stacks
- Collaborating across different systems
- Using pre-built scientific containers

### Use default (no context) when:
- Running simple shell commands
- Using software available in base environment
- Testing basic functionality

## Debugging Execution Contexts

### Check Available Environments

```bash
# See conda environments
conda env list

# See available modules
module avail

# Check current environment
qxub exec --env dvc3 -- python check_environment.py
```

### Troubleshooting Environment Issues

Use `--dry` and `-v` to debug environment problems:

```bash
qxub exec --dry -v --env nonexistent -- python test.py
```

This will show you exactly how qxub tries to activate the environment and where it might fail.

## Practical Examples

### Bioinformatics Pipeline

```bash
# Using pysam for genomics work
qxub exec --env pysam --resources mem=8GB -- python3 -c "
import pysam
import sys
print(f'Pysam version: {pysam.__version__}')
print(f'Python: {sys.version.split()[0]}')
print('Genomics environment ready')
"
```

### Machine Learning with Resources

```bash
# ML training in sc environment
qxub exec --env sc --resources mem=32GB,ncpus=4,walltime=2:00:00 -- python3 -c "
import scanpy as sc
import pandas as pd
import numpy as np
print('ML environment loaded')
print(f'Available CPUs for parallel processing: 4')
print('Training could start here...')
"
```

### Scientific Computing

```bash
# Numerical computation with modules
qxub exec --mod python3/3.11.7 --resources mem=16GB,ncpus=8 -- python3 -c "
import multiprocessing
print(f'CPUs available: {multiprocessing.cpu_count()}')
print('Scientific computing environment ready')
"
```

## Key Takeaways

1. **One context per job**: Cannot mix `--env`, `--mod`, `--mods`, and `--sif`
2. **Environment isolation**: Each context provides a clean, predictable environment
3. **Automatic setup**: qxub handles all activation/loading automatically
4. **Resource compatibility**: All execution contexts work with custom resources
5. **Debugging support**: Use `--dry` to see exactly how environments are activated

## Next Steps

Now that you understand execution contexts:
- **[Complex Commands](05-complex-commands.md)** - Handle quotes, variables, and multi-line commands
- **[Aliases](07-aliases.md)** - Save common environment + resource combinations

Execution contexts are fundamental to productive HPC work. They ensure your jobs run in the exact software environment you need, every time.

---

**💡 Pro Tips:**
- Use `--dry` to verify environment activation before submitting
- Combine execution contexts with appropriate resources (`--resources mem=8GB,ncpus=2`)
- Keep commonly used environment+resource combinations as aliases
- Test environment availability locally before using in jobs
