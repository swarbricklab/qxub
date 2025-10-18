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
qxub --dry --env dvc3 -- python check_environment.py
```

**Expected dry run output:**
```
📝 PBS Script Preview:
#!/bin/bash
#PBS -N qx-20241017-144052
...

cd "/g/data/a56/software/qsub_tools"

# Activate conda environment
source /g/data/a56/conda/miniconda3/etc/profile.d/conda.sh
conda activate dvc3

python3 -c "
import sys
print(f'Python: {sys.version}')
import pandas as pd
print(f'Pandas version: {pd.__version__}')
"
```

Notice how qxub automatically:
- Sources the conda initialization script
- Activates the specified environment
- Runs your command in that environment

### Running the Conda Job

```bash
qxub --env dvc3 -- python3 -c "
import sys, pandas as pd, numpy as np
print(f'Python: {sys.version.split()[0]}')
print(f'Pandas: {pd.__version__}')
print(f'NumPy: {np.__version__}')
print('All imports successful!')
"
```

**Expected output:**
```
🚀 Submitting job...
📋 Job submitted: 12345689.gadi-pbs (qx-20241017-145052)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

Python: 3.11.7
Pandas: 2.1.4
NumPy: 1.24.3
All imports successful!

🎉 Job completed successfully (exit code: 0)
```

### Data Science Example with Resources

```bash
# Pandas analysis in a specialized environment
qxub --env dvc3 --resources mem=8GB,ncpus=2 -- python3 -c "
import pandas as pd
import numpy as np
print('Creating sample dataset...')
df = pd.DataFrame(np.random.randn(1000, 4), columns=['A', 'B', 'C', 'D'])
print(f'Dataset shape: {df.shape}')
print('Computing statistics...')
print(df.describe())
print('Analysis complete!')
"
```

### R Environment Example

```bash
# Using the tidyverse environment for R
qxub --env tidyverse -- Rscript -e "
library(tidyverse)
print('R version:')
print(R.version.string)
print('Tidyverse loaded successfully')
data <- mtcars %>%
  group_by(cyl) %>%
  summarise(mean_mpg = mean(mpg))
print(data)
"
```

### Single-Cell Analysis Environment

```bash
# Single-cell analysis with specialized packages
qxub --env sc --resources mem=16GB -- python3 -c "
import scanpy as sc
import pandas as pd
print('Scanpy version:', sc.__version__)
print('Single-cell analysis environment ready')
# Could load and analyze data here
print('Environment validation complete')
"
```

## Environment Modules with `--mod` and `--mods`

### Single Module (`--mod`)

Load a single environment module:

```bash
qxub --dry --mod python3/3.11.7 -- python3 -c "
import sys
print(f'Python: {sys.version}')
"
```

**Expected dry run output:**
```
📝 PBS Script Preview:
#!/bin/bash
#PBS -N qx-20241017-145152
...

cd "/g/data/a56/software/qsub_tools"

# Load environment modules
module load python3/3.11.7

python3 -c "
import sys
print(f'Python: {sys.version}')
"
```

### Multiple Modules (`--mods`)

Load multiple modules at once:

```bash
qxub --dry --mods python3/3.11.7,gcc/11.1.0 -- ./compile_project.sh
```

**Expected dry run output:**
```
📝 PBS Script Preview:
...
# Load environment modules
module load python3/3.11.7 gcc/11.1.0
...
```

### Real Module Example

```bash
# Load Python and run a simple computation
qxub --mod python3/3.11.7 -- python3 -c "
import sys
import math
print(f'Python {sys.version.split()[0]} loaded via module')
print('Computing square roots...')
for i in [1, 4, 9, 16, 25]:
    print(f'sqrt({i}) = {math.sqrt(i)}')
print('Computation complete!')
"
```

### Combining Modules for Compilation

```bash
# Example: compiling with GCC (if you had source code)
qxub --mods gcc/11.1.0,python3/3.11.7 -- ./compile_project.sh
```

## Singularity Containers with `--sif`

qxub supports Singularity containers for reproducible environments:

```bash
# Example (if you had a container)
qxub --dry --sif /path/to/container.sif -- python analysis.py
```

**Expected dry run output:**
```
📝 PBS Script Preview:
...
# Run in Singularity container
singularity exec /path/to/container.sif python analysis.py
```

**Note**: Container examples are not included here as they require specific container files, but the syntax follows the same pattern as other execution contexts.

## Execution Context Rules and Errors

### Cannot Mix Contexts

```bash
# This will fail - mixing conda and modules
qxub --env dvc3 --mod python3/3.11.7 -- echo "This fails"
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
qxub --env dvc3 -- python3 script.py

# ✅ Good: Single module only
qxub --mod python3/3.11.7 -- python3 script.py

# ✅ Good: Multiple modules only
qxub --mods python3/3.11.7,gcc/11.1.0 -- python3 script.py

# ✅ Good: Container only
qxub --sif container.sif -- python3 script.py

# ✅ Good: No execution context (uses login environment)
qxub --default -- python3 script.py
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
qxub --env dvc3 -- python check_environment.py
```

### Troubleshooting Environment Issues

Use `--dry` and `-v` to debug environment problems:

```bash
qxub --dry -v --env nonexistent -- python test.py
```

This will show you exactly how qxub tries to activate the environment and where it might fail.

## Practical Examples

### Bioinformatics Pipeline

```bash
# Using pysam for genomics work
qxub --env pysam --resources mem=8GB -- python3 -c "
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
qxub --env sc --resources mem=32GB,ncpus=4,walltime=2:00:00 -- python3 -c "
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
qxub --mod python3/3.11.7 --resources mem=16GB,ncpus=8 -- python3 -c "
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
