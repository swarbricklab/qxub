# Complex Commands: Handling Quotes, Variables, and Multi-line Scripts

Now that you're comfortable with basic qxub usage and execution contexts, it's time to tackle more complex commands. This section covers shell quoting, variable substitution, multi-line scripts, and the powerful `--cmd` option.

## Understanding Command Processing

qxub processes your command in several steps:
1. Parse options before `--`
2. Pass everything after `--` to the PBS job
3. Encode the command to avoid shell escaping issues
4. Execute in the specified environment

This means you need to be careful with quotes, variables, and complex shell constructs.

## Basic Quoting Patterns

### Simple Quotes

For commands with spaces or special characters:

```bash
# Single quotes preserve everything literally
qxub -- echo 'Hello World! $USER will not be expanded'
```

**Expected output:**
```
Hello World! $USER will not be expanded
```

```bash
# Double quotes allow variable expansion
qxub -- echo "Hello World! Current user is $USER"
```

**Expected output:**
```
Hello World! Current user is jr9959
```

### Shell Variables

Variables are expanded by the shell when you type the command:

```bash
# Variable expanded on the login node (when you type it)
export MY_VAR="test123"
qxub -- echo "My variable: $MY_VAR"
```

**Expected output:**
```
My variable: test123
```

But be careful with variables that should be expanded on the compute node:

```bash
# This expands $HOSTNAME on the login node (probably not what you want)
qxub -- echo "Running on: $HOSTNAME"

# Better: escape the $ to expand on compute node
qxub -- echo "Running on: \$HOSTNAME"
```

**Expected output:**
```
Running on: gadi-cpu-clx-XXXX
```

## Multi-line Commands and Scripts

### Using Bash -c for Complex Logic

```bash
qxub -- bash -c '
echo "Starting complex script..."
for i in {1..5}; do
    echo "Processing item $i"
    sleep 1
done
echo "Script completed successfully"
'
```

**Expected output:**
```
🚀 Submitting job...
📋 Job submitted: 12345690.gadi-pbs (qx-20241017-150052)
⏳ Job queued, waiting for execution...
✅ Job started, streaming output...

Starting complex script...
Processing item 1
Processing item 2
Processing item 3
Processing item 4
Processing item 5
Script completed successfully

🎉 Job completed successfully (exit code: 0)
```

### Python Multi-line Scripts

```bash
qxub --env dvc3 -- python3 -c '
import pandas as pd
import numpy as np

print("Creating sample data...")
data = {
    "name": ["Alice", "Bob", "Charlie"],
    "age": [25, 30, 35],
    "score": [85.5, 92.0, 78.5]
}

df = pd.DataFrame(data)
print("Data created:")
print(df)

print(f"Average age: {df.age.mean():.1f}")
print(f"Average score: {df.score.mean():.1f}")
'
```

### R Multi-line Scripts

```bash
qxub --env tidyverse -- Rscript -e '
library(tidyverse)

cat("Creating sample data...\n")
df <- data.frame(
  name = c("Alice", "Bob", "Charlie"),
  age = c(25, 30, 35),
  score = c(85.5, 92.0, 78.5)
)

print(df)

cat("Summary statistics:\n")
df %>%
  summarise(
    avg_age = mean(age),
    avg_score = mean(score)
  ) %>%
  print()
'
```

## The `--cmd` Option: Advanced Command Handling

The `--cmd` option provides an alternative way to specify commands, especially useful for complex cases:

### Basic `--cmd` Usage

```bash
# Instead of using --
qxub --cmd 'echo "Hello from --cmd"'
```

This is equivalent to:
```bash
qxub -- echo "Hello from --cmd"
```

### When to Use `--cmd`

The `--cmd` option is particularly useful when:
- Your command contains options that look like qxub options
- You need to read commands from files or variables
- You're building commands programmatically

### Reading Commands from Files

Create a simple script file:
```bash
echo 'echo "Script from file"
date
hostname' > /tmp/my_script.sh
```

Then run it:
```bash
qxub --cmd "bash /tmp/my_script.sh"
```

### Complex Quoting with `--cmd`

```bash
# Handling nested quotes
qxub --cmd 'python3 -c "
import subprocess
result = subprocess.run([\"echo\", \"nested quotes work\"], capture_output=True, text=True)
print(f\"Output: {result.stdout.strip()}\")
"'
```

## Variable Substitution Patterns

### Environment Variables in Jobs

```bash
# Variables expanded on compute node
qxub -- bash -c '
echo "Job running on: $HOSTNAME"
echo "Working directory: $PWD"
echo "User: $USER"
echo "Temporary directory: $TMPDIR"
echo "Number of CPUs: $PBS_NCPUS"
'
```

### Passing Local Variables to Jobs

```bash
# Method 1: Export and reference
export DATA_DIR="/scratch/a56/$USER/data"
qxub -- bash -c "echo 'Data directory: $DATA_DIR'"

# Method 2: Inline definition
qxub -- bash -c 'DATA_DIR="/scratch/a56/$USER/data"; echo "Data directory: $DATA_DIR"'
```

### Template Variable Expansion

qxub supports template variables in job names and paths:

```bash
# These are expanded by qxub, not the shell
qxub --name "analysis-{timestamp}" --out "/scratch/a56/{user}/results-{timestamp}.out" -- echo "Template test"
```

## Practical Complex Examples

### Data Processing Pipeline

```bash
qxub --env dvc3 --mem 8GB --walltime 1:00:00 -- bash -c '
echo "Starting data processing pipeline..."

# Set up directories
WORK_DIR="/scratch/a56/$USER/analysis_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

echo "Working in: $WORK_DIR"

# Create sample data
python3 -c "
import pandas as pd
import numpy as np

# Generate sample dataset
np.random.seed(42)
data = pd.DataFrame({
    \"id\": range(1000),
    \"value\": np.random.randn(1000),
    \"category\": np.random.choice([\"A\", \"B\", \"C\"], 1000)
})

data.to_csv(\"sample_data.csv\", index=False)
print(f\"Created dataset with {len(data)} rows\")
print(data.head())
"

echo "Data processing completed successfully"
echo "Results saved in: $WORK_DIR"
'
```

### Bioinformatics Example

```bash
qxub --env pysam --mem 16GB --ncpus 4 -- bash -c '
echo "Bioinformatics analysis starting..."

# Set analysis parameters
SAMPLE_ID="sample_001"
OUTPUT_DIR="/scratch/a56/$USER/analysis/$SAMPLE_ID"
mkdir -p "$OUTPUT_DIR"

echo "Sample ID: $SAMPLE_ID"
echo "Output directory: $OUTPUT_DIR"
echo "CPUs available: $PBS_NCPUS"

# Simulate analysis steps
python3 -c "
import pysam
import time

print(f\"Pysam version: {pysam.__version__}\")
print(\"Step 1: Quality control...\")
time.sleep(2)

print(\"Step 2: Alignment...\")
time.sleep(3)

print(\"Step 3: Variant calling...\")
time.sleep(2)

print(\"Analysis pipeline completed!\")
"
'
```

### Parameter Sweep

```bash
# Running analysis with different parameters
for param in 0.1 0.5 1.0; do
    qxub --name "param-$param" --env dvc3 -- python3 -c "
import numpy as np
import time

param = $param
print(f'Running analysis with parameter: {param}')

# Simulate computation
np.random.seed(42)
result = np.random.randn(100).mean() * param
time.sleep(2)

print(f'Result with param={param}: {result:.4f}')
print('Analysis completed')
"
done
```

## Debugging Complex Commands

### Use `--dry` with Complex Commands

Always test complex commands with `--dry` first:

```bash
qxub --dry --env dvc3 -- bash -c '
echo "Complex script..."
for i in {1..3}; do
    echo "Step $i"
done
'
```

This shows you exactly how your command will be encoded and executed.

### Common Pitfalls and Solutions

#### Problem: Variable Expansion Timing
```bash
# ❌ Wrong: $USER expanded on login node
qxub -- echo "User: $USER"

# ✅ Right: \$USER expanded on compute node
qxub -- echo "User: \$USER"
```

#### Problem: Nested Quotes
```bash
# ❌ Wrong: Quote confusion
qxub -- python3 -c "print("Hello")"

# ✅ Right: Escape inner quotes
qxub -- python3 -c "print(\"Hello\")"

# ✅ Alternative: Use single quotes outside
qxub -- python3 -c 'print("Hello")'
```

#### Problem: Special Characters
```bash
# ❌ Wrong: Semicolon interpreted by login shell
qxub -- echo "First"; echo "Second"

# ✅ Right: Wrap in bash -c
qxub -- bash -c 'echo "First"; echo "Second"'
```

## File-based Scripts

For very complex scripts, consider using files:

### Create Script File
```bash
cat > /tmp/analysis_script.py << 'EOF'
#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

print("Starting analysis script...")
print(f"Python version: {sys.version.split()[0]}")

# Your complex analysis here
data = pd.DataFrame({"x": np.random.randn(100)})
print(f"Generated {len(data)} data points")
print(f"Mean: {data.x.mean():.3f}")
print(f"Std: {data.x.std():.3f}")

print("Analysis completed successfully")
EOF

chmod +x /tmp/analysis_script.py
```

### Run Script File
```bash
qxub --env dvc3 --mem 8GB -- python3 /tmp/analysis_script.py
```

## Key Takeaways

1. **Quote carefully**: Understand when variables are expanded (login vs compute node)
2. **Use `bash -c`**: For multi-line scripts and complex shell logic
3. **Test with `--dry`**: Always verify complex commands before submission
4. **Consider `--cmd`**: Alternative syntax for complex cases
5. **Use files for complexity**: Very complex scripts are better as separate files

## Next Steps

Now that you can handle complex commands:
- **[Job History](06-history.md)** - Track and analyze your job executions
- **[Aliases](07-aliases.md)** - Save complex command patterns for reuse

Complex command handling is crucial for real-world HPC work. The patterns shown here will handle most scenarios you'll encounter.

---

**💡 Pro Tips:**
- Always use `--dry` to test complex commands before submission
- Use single quotes to prevent unwanted variable expansion
- Wrap complex logic in `bash -c '...'`
- Consider creating reusable script files for very complex operations
- The `--cmd` option provides an alternative when `--` syntax becomes unwieldy
