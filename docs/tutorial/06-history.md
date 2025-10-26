# Job History: Tracking and Analyzing Your Work

qxub automatically tracks all your job submissions, creating a comprehensive history that helps you monitor progress, debug issues, and optimize resource usage. This section covers how to view, search, and analyze your job history.

**New in v3.2.0:** Job output viewing commands (`qxub history out/err/log`) make it easy to quickly check job results and debug issues.

## Understanding qxub History

qxub maintains two types of historical information:
1. **Execution records** - What you ran, when, and with what resources (now including file paths)
2. **Computational recipes** - Reusable job templates indexed by unique hash

All history is stored locally and tied to your user account.

## Basic History Commands

### View Recent Jobs

```bash
# Show the 10 most recent jobs
qxub history executions

# Show specific number of recent jobs
qxub history executions --limit 20
```

**Expected output:**
```
                                  Recent Executions
┏━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┓
┃ Time           ┃ Recipe   ┃ Command             ┃ Status    ┃ Directory           ┃
┡━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━┩
│ 10-26 21:01:17 │ 54d59eab │ conda/envs/usr/jr9… │ completed │ g/data/a56/softwar… │
│                │          │ exec --dry --env    │           │                     │
│                │          │ base --...          │           │                     │
│ 10-26 20:53:14 │ 603cb1a9 │ conda/envs/usr/jr9… │ completed │ g/data/a56/softwar… │
│                │          │ exec --mem 1GB      │           │                     │
│                │          │ --name te...        │           │                     │
│ 10-26 20:49:57 │ 55dd762f │ conda/envs/usr/jr9… │ completed │ g/data/a56/softwar… │
│                │          │ exec --dry --mem    │           │                     │
│                │          │ 1GB --n...          │           │                     │
└────────────────┴──────────┴─────────────────────┴───────────┴─────────────────────┘
```

### View Latest Job Details

```bash
# Show detailed information about the most recent job
qxub history latest
```

### Browse Computational Recipes

```bash
# Show all saved computational recipes (reusable job templates)
qxub history recipes

# Show executions for a specific recipe
qxub history runs <recipe_hash>
```

## Detailed Job Information

### View Latest Job Details

```bash
# Get detailed info about the most recent job
qxub history latest
```

**Expected output:**
```
📋 Latest Execution

Recipe: 54d59eab4954
Command: /g/data/a56/conda/envs/usr/jr9959/qxub/bin/qxub exec --env dvc3 --mem 4GB
-- python analysis.py
Directory: /g/data/a56/user/project
Status: completed
Time: 2025-10-26T20:53:14.825717

🏗️ Recipe Structure:
  1 executor:
  2   env: dvc3
  3   mem: 4GB
  4   type: conda
  5 target:
  6   cmd: python analysis.py
  7
```

### View Recipe Details

```bash
# Show detailed information about a computational recipe
qxub history show <recipe_hash>
```

This shows the full recipe definition that can be reused or converted to an alias.

### View Job Output Files (v3.2.0)

qxub now provides convenient commands to quickly view job output files:

```bash
# View stdout from most recent job
qxub history out

# View stderr from most recent job
qxub history err

# View PBS log from most recent job
qxub history log

# View output from specific job by ID
qxub history out 12345691.gadi-pbs
qxub history err 12345691.gadi-pbs
qxub history log 12345691.gadi-pbs
```

**Expected output:**
```bash
$ qxub history out
📄 Stdout from job most recent
(/scratch/a56/jr9959/qxub/analysis_20251026_205314.out):
Processing dataset...
Analysis complete: 1,234 records processed
Results saved to output.csv
✅ Command completed successfully
🎉 Job completed successfully

$ qxub history err 12345691.gadi-pbs
📄 Stderr file is empty:
/scratch/a56/jr9959/qxub/analysis_20251026_205314.err

$ qxub history log
📄 PBS log from job most recent
(/scratch/a56/jr9959/qxub/analysis_20251026_205314.log):
📁 Execution directory: /g/data/a56/user/project
🔧 Conda environment: dvc3
💻 Command: python analysis.py
📄 STDOUT file: /scratch/a56/jr9959/qxub/analysis_20251026_205314.out
📄 STDERR file: /scratch/a56/jr9959/qxub/analysis_20251026_205314.err
🚀 Job started at 2025-10-26T20:53:28+11:00
💻 Executing main command...
...
✅ Job completed successfully
```

**Key Features:**
- **Smart defaults**: Commands without job ID show most recent job output
- **Syntax highlighting**: Output is formatted for easy reading
- **File status**: Clear indication if files are empty or missing
- **Rich formatting**: Headers show job ID and file paths for context

## Working with History and Recipes

### Understanding the History System

qxub tracks two types of data:
1. **Computational Recipes** - Unique job templates indexed by hash (what you want to run)
2. **Execution Records** - Individual job submissions with timestamps (when you ran it)

This separation allows qxub to:
- Deduplicate similar jobs into reusable recipes
- Track multiple executions of the same recipe
- Provide templates for frequently used job patterns

### Exploring Recipes

```bash
# List all computational recipes
qxub history recipes

# View details of a specific recipe
qxub history show <recipe_hash>

# See all executions of a specific recipe
qxub history runs <recipe_hash>

# Convert a recipe to an alias for easy reuse
qxub history to-alias <recipe_hash> <alias_name>
```

### Managing History

```bash
# Clear all history (use with caution!)
qxub history clear

# This removes both recipes and execution records
# Useful for testing or starting fresh
```

## Practical History Usage Patterns

### Debugging Workflow Issues

When something goes wrong with a job, history helps you quickly investigate:

```bash
# 1. Check recent executions to identify the problem job
qxub history executions --limit 10

# 2. View the specific execution details
qxub history latest  # or qxub history show <recipe_hash>

# 3. Check the job output for errors
qxub history err     # stderr from most recent job
qxub history out     # stdout from most recent job
qxub history log     # PBS log with full job details

# 4. If you know a specific job ID that failed
qxub history err 12345.gadi-pbs
```

### Reusing Successful Job Patterns

```bash
# 1. Find a job pattern that worked well
qxub history recipes

# 2. View the details of a successful recipe
qxub history show 54d59eab4954

# 3. See how many times you've used this pattern
qxub history runs 54d59eab4954

# 4. Convert it to an alias for easy reuse
qxub history to-alias 54d59eab4954 "data-analysis"

# 5. Now use the alias
qxub alias data-analysis --mem 8GB -- python new_dataset.py
```

### Tracking Job Progress Over Time

```bash
# Check your recent work
qxub history executions --limit 20

# Look at a specific recipe's execution history
qxub history runs cf7cc954  # shows all times this pattern was used

# View detailed information about successful patterns
qxub history show cf7cc954  # see the full recipe structure
```

## Working with Job Output Files

The new v3.2.0 logging features make it easy to check job results:

### Quick Output Checks

```bash
# Check if your most recent job produced the expected output
qxub history out | grep "Analysis complete"

# Check for any errors in the most recent job
qxub history err

# View the PBS log to understand resource usage
qxub history log | grep "Memory Used"
```

### Debugging Failed Jobs

```bash
# When a job fails, follow this pattern:
# 1. Check recent executions to identify the failed job
qxub history executions --limit 5

# 2. Look at the error output
qxub history err  # most recent job
# or for a specific job:
qxub history err 12345.gadi-pbs

# 3. Check the PBS log for system-level issues
qxub history log | tail -20
```

### Comparing Different Runs

```bash
# Compare output from different executions of similar jobs
qxub history out 12345.gadi-pbs > run1_output.txt
qxub history out 12346.gadi-pbs > run2_output.txt
diff run1_output.txt run2_output.txt
```

## Recipe Management Best Practices

### Understanding Recipe Hashes

Each unique combination of execution context and options creates a recipe:

```bash
# These create different recipes:
qxub exec --env dvc3 --mem 4GB -- python script.py      # Recipe A
qxub exec --env dvc3 --mem 8GB -- python script.py      # Recipe B
qxub exec --env base --mem 4GB -- python script.py      # Recipe C
```

### Converting Recipes to Aliases

When you find a recipe you use frequently:

```bash
# 1. Find the recipe hash in your executions
qxub history executions

# 2. Convert it to an alias
qxub history to-alias 54d59eab4954 "my-analysis"

# 3. Use the alias with modifications
qxub alias my-analysis --name "analysis-v2" -- python script.py --version 2
```

## History Maintenance

### Clear History

```bash
# Clear all history (recipes and executions)
# Use with caution - this removes all tracked data
qxub history clear

# This is useful when:
# - Starting fresh with a new project
# - Testing and development
# - Cleaning up after many experimental runs
```

**Note:** There's currently no selective deletion - clearing removes both recipes and execution records. Consider this before clearing if you have valuable job templates.

## Available History Commands Summary

Here's a complete reference of the history commands available in qxub:

### Core Commands
- `qxub history executions [--limit N]` - List recent job executions
- `qxub history latest` - Show details of most recent execution
- `qxub history recipes` - List all computational recipes
- `qxub history show <recipe_hash>` - Show recipe details
- `qxub history runs <recipe_hash>` - Show executions for a recipe

### New in v3.2.0: Log Viewing
- `qxub history out [job_id]` - View job stdout (defaults to most recent)
- `qxub history err [job_id]` - View job stderr (defaults to most recent)
- `qxub history log [job_id]` - View PBS log (defaults to most recent)

### Management Commands
- `qxub history clear` - Clear all history data
- `qxub history to-alias <recipe_hash> <alias_name>` - Convert recipe to alias

## Key Takeaways

1. **Automatic tracking**: All jobs are automatically recorded with file paths
2. **Recipe deduplication**: Similar jobs are grouped into reusable recipes
3. **Easy log viewing**: New v3.2.0 commands make checking job output effortless
4. **Template creation**: Convert successful recipes to aliases for reuse
5. **Debugging support**: Quickly access job output files for troubleshooting

## Next Steps

Now that you understand job history:
- **[Aliases](07-aliases.md)** - Save optimized resource combinations as shortcuts
- **[Configuration](08-configuration.md)** - Understand how settings affect all jobs

Job history is invaluable for optimizing your HPC workflows. Use it regularly to identify successful patterns and debug issues.

---

**💡 Pro Tips:**
- Use `qxub history out` after every job to quickly check results
- Convert frequently-used recipes to aliases with `qxub history to-alias`
- Check `qxub history err` first when debugging failed jobs
- Use `qxub history latest` to understand the structure of your last job
- The recipe hash in executions can be used with `qxub history show` for details
