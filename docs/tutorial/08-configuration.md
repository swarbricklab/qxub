# Configuration: Understanding Config Levels and Customization

qxub uses a sophisticated hierarchical configuration system that provides sensible defaults while allowing customization at multiple levels. This section explains how configuration works and when you might need to customize it.

## Configuration Hierarchy

qxub loads configuration from multiple sources with this precedence (highest to lowest):

1. **Command-line options** - Highest priority, overrides everything
2. **Project configuration** - `.qxub/config.yaml` in current directory
3. **User configuration** - `~/.config/qxub/config.yaml`
4. **System configuration** - `/g/data/a56/config/xdg/qxub/config.yaml`

**For most users, the system configuration provides everything needed.**

## Viewing Current Configuration

### See All Configuration with Origins

```bash
qxub config list --show-origin
```

**Expected output:**
```
📋 Configuration (with origin)
├── defaults:
│   ├── project: a56
│   │   (/g/data/a56/config/xdg/qxub/config.yaml)
│   ├── queue: normal
│   │   (/g/data/a56/config/xdg/qxub/config.yaml)
│   ├── resources: ['mem=4GB', 'ncpus=1', 'walltime=2:00:00']
│   │   (/g/data/a56/config/xdg/qxub/config.yaml)
│   └── name: qx-{timestamp}
│       (/g/data/a56/config/xdg/qxub/config.yaml)
├── templates:
│   ├── scratch_base: /scratch/{project}/{user}
│   │   (/g/data/a56/config/xdg/qxub/config.yaml)
│   └── log_dir: {scratch_base}/qxub
│       (/g/data/a56/config/xdg/qxub/config.yaml)
└── aliases:
    ├── py: {...}
    │   (/g/data/a56/config/xdg/qxub/config.yaml)
    └── r: {...}
        (/g/data/a56/config/xdg/qxub/config.yaml)
```

This shows where each setting comes from, helping you understand what can be overridden and where.

### View Specific Configuration Sections

```bash
# View just defaults
qxub config get defaults

# View template variables
qxub config get templates

# View platform settings
qxub config get platform_search_paths
```

## Understanding the System Configuration

Most users work exclusively with the system configuration provided for project a56. Let's understand what it provides:

### Default Resource Allocation

```bash
qxub config get defaults.resources
```

**Output:**
```
📋 Configuration: defaults.resources
└── ['mem=4GB', 'ncpus=1', 'walltime=2:00:00']
```

These defaults work well for:
- Simple data processing scripts
- Quick analysis tasks
- Single-threaded applications

### Template Variables

```bash
qxub config get templates
```

**Output:**
```
📋 Configuration: templates
├── scratch_base: /scratch/{project}/{user}
├── log_dir: {scratch_base}/qxub
└── job_prefix: a56-{user}
```

These variables are automatically expanded:
- `{project}` → `a56`
- `{user}` → Your username (e.g., `jr9959`)
- `{timestamp}` → Current timestamp (e.g., `20241017-151052`)

### Platform Configuration

```bash
qxub config get platform_search_paths
```

**Output:**
```
📋 Configuration: platform_search_paths
└── ['/g/data/a56/config/xdg/qxub/platforms', '/etc/qxub/platforms', '/scratch/{project}/{user}/qxub/platforms']
```

This defines where qxub looks for platform definitions (like `nci_gadi.yaml`).

## When You Need User Configuration

**Most users never need personal configuration.** However, you might create user config for:

### Personal Default Preferences

If you consistently want different defaults:

```bash
# Create user config directory
mkdir -p ~/.config/qxub

# Set personal defaults
cat > ~/.config/qxub/config.yaml << EOF
defaults:
  resources:
    - mem=8GB        # You prefer more memory by default
    - ncpus=2        # You often use parallel processing
    - walltime=4:00:00  # You run longer jobs
  queue: auto        # You prefer automatic queue selection

# Personal aliases
aliases:
  myanalysis:
    main:
      env: dvc3
      resources: ["mem=16GB", "ncpus=4", "walltime=3:00:00"]
EOF
```

### Personal Template Customization

```bash
# Add to your user config
cat >> ~/.config/qxub/config.yaml << EOF
templates:
  my_data_dir: /scratch/a56/{user}/my_projects
  my_results_dir: {my_data_dir}/results
EOF
```

### Viewing User Configuration Impact

```bash
# After creating user config, see the changes
qxub config list --show-origin | grep "~/.config"
```

## Project-Level Configuration

Project configuration is useful when working on team projects that need consistent settings:

### Creating Project Configuration

```bash
# In your project directory
mkdir -p .qxub

cat > .qxub/config.yaml << EOF
# Project-specific configuration for MyProject

defaults:
  name: "myproject-{timestamp}"
  resources:
    - mem=16GB       # This project needs more memory
    - ncpus=4        # Parallel processing required
    - walltime=6:00:00  # Long-running analyses

templates:
  project_data: /scratch/a56/{user}/myproject/data
  project_results: /scratch/a56/{user}/myproject/results

aliases:
  preprocess:
    main:
      env: dvc3
      resources: ["mem=8GB", "ncpus=2", "walltime=2:00:00"]
      out: "{project_results}/preprocess_{timestamp}.out"

  analyze:
    main:
      env: sc
      resources: ["mem=32GB", "ncpus=8", "walltime=8:00:00"]
      out: "{project_results}/analyze_{timestamp}.out"
EOF
```

### Using Project Configuration

```bash
# When you run qxub from this project directory:
qxub preprocess -- python3 preprocess_data.py
qxub analyze -- python3 run_analysis.py

# The project config automatically applies
qxub config get defaults.resources
# Shows: ['mem=16GB', 'ncpus=4', 'walltime=6:00:00']
```

### Sharing Project Configuration

Include `.qxub/config.yaml` in your project repository:

```bash
# Add to version control
git add .qxub/config.yaml
git commit -m "Add qxub project configuration"
```

Now team members get consistent qxub settings when working on the project.

## Configuration File Locations

### Check Where Files Are Located

```bash
qxub config files
```

**Expected output:**
```
📁 Configuration File Locations:

System Configuration:
├── Path: /g/data/a56/config/xdg/qxub/config.yaml
├── Status: ✅ Found and loaded
└── Purpose: Team a56 defaults and aliases

User Configuration:
├── Path: ~/.config/qxub/config.yaml
├── Status: 📂 Not found (optional)
└── Purpose: Personal preferences and aliases

Project Configuration:
├── Path: .qxub/config.yaml
├── Status: 📂 Not found in current directory
└── Purpose: Project-specific settings

Platform Definitions:
├── System: /g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml ✅
├── User: ~/.config/qxub/platforms/ (not found)
└── Project: .qxub/platforms/ (not found)
```

## Advanced Configuration Features

### Template Variable Expansion

Templates support recursive expansion:

```yaml
templates:
  base_dir: /scratch/{project}/{user}
  data_dir: {base_dir}/data
  results_dir: {base_dir}/results
  log_dir: {base_dir}/logs

defaults:
  out: "{log_dir}/{name}_{timestamp}.out"
  err: "{log_dir}/{name}_{timestamp}.err"
```

### Conditional Configuration

While not directly supported, you can use different project directories for different configurations:

```bash
# Different project configurations
mkdir -p ~/projects/quick_analysis/.qxub
mkdir -p ~/projects/deep_analysis/.qxub

# Quick analysis config
cat > ~/projects/quick_analysis/.qxub/config.yaml << EOF
defaults:
  resources: ["mem=4GB", "ncpus=1", "walltime=30:00"]
  queue: express
EOF

# Deep analysis config
cat > ~/projects/deep_analysis/.qxub/config.yaml << EOF
defaults:
  resources: ["mem=32GB", "ncpus=8", "walltime=12:00:00"]
  queue: normal
EOF
```

## Configuration Validation

### Test Configuration Changes

```bash
# Test how your config changes affect job submission
qxub exec --dry -- echo "Testing config"
```

### Verify Template Expansion

```bash
# See how templates are expanded
qxub config get templates --expand
```

## Managing Configuration

### Reset to System Defaults

```bash
# Remove user configuration to use only system defaults
rm -f ~/.config/qxub/config.yaml

# Verify reset worked
qxub config list --show-origin | grep -v "system"
```

### Backup Configuration

```bash
# Backup your user configuration
cp ~/.config/qxub/config.yaml ~/.config/qxub/config.yaml.backup

# Backup project configuration
cp .qxub/config.yaml .qxub/config.yaml.backup
```

## Common Configuration Patterns

### Research Group Setup

For a research group, you might use system config for:

```yaml
# System config (/g/data/a56/config/xdg/qxub/config.yaml)
defaults:
  project: a56
  queue: normal

templates:
  group_data: /g/data/a56/data
  group_scratch: /scratch/a56/{user}

aliases:
  genomics:
    main:
      env: pysam
      resources: ["mem=16GB", "ncpus=8", "walltime=4:00:00"]
```

### Personal Productivity Setup

Individual users might add personal config for:

```yaml
# User config (~/.config/qxub/config.yaml)
defaults:
  queue: auto  # Always use automatic queue selection

aliases:
  quick:
    main:
      resources: ["mem=2GB", "ncpus=1", "walltime=15:00"]
      queue: express

  notebook:
    main:
      env: jupyterlab
      resources: ["mem=8GB", "ncpus=2", "walltime=2:00:00"]
```

## Troubleshooting Configuration

### Configuration Not Loading

```bash
# Check file permissions
ls -la ~/.config/qxub/config.yaml

# Validate YAML syntax
python3 -c "import yaml; yaml.safe_load(open('~/.config/qxub/config.yaml'))"
```

### Unexpected Settings

```bash
# See exactly what's being applied and from where
qxub config list --show-origin

# Test with verbose output
qxub -vv --dry -- echo "test"
```

### Template Variable Issues

```bash
# Check template expansion
qxub config get templates

# Test template in dry run
qxub exec --dry --name "test-{timestamp}" -- echo "test"
```

## Key Takeaways

1. **System config is sufficient**: Most users only need the provided team configuration
2. **Hierarchical precedence**: Command line > Project > User > System
3. **Project config for teams**: Share consistent settings via `.qxub/config.yaml`
4. **User config for preferences**: Personal defaults and aliases
5. **Templates enable flexibility**: Dynamic path and name generation

## Next Steps

Now that you understand configuration:
- **[Parallel Execution](09-parallel-execution.md)** - Advanced parallel job patterns
- **[DVC Integration](10-dvc-integration.md)** - Using qxub in data science pipelines

Configuration understanding helps you customize qxub to your workflow while maintaining team consistency.

---

**💡 Pro Tips:**
- Start with system defaults - they're designed for common use cases
- Use `qxub config list --show-origin` to understand setting sources
- Create project configs for team consistency, user configs for personal preferences
- Test configuration changes with `--dry` before relying on them
- Include `.qxub/config.yaml` in version control for team projects
