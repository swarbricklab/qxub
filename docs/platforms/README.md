# Platform Definitions

This directory contains complete platform definitions for various HPC systems. These serve as reference implementations and can be installed in system-level configuration.

## Available Platforms

### NCI Gadi (`nci_gadi.yaml`)
Complete queue definitions for the National Computational Infrastructure Gadi supercomputer, including:

- **Standard queues**: normal, express
- **Memory queues**: hugemem, megamem
- **GPU queues**: gpuvolta, dgxa100
- **Architecture-specific**: normalbw, normalsl, normalsr (Broadwell, Skylake, Sapphire Rapids)
- **Special purpose**: copyq (data transfer with internet)

Based on official NCI documentation: https://opus.nci.org.au/spaces/Help/pages/90308823/Queue+Limits

### This machine, no scheduler (`local/local.yaml`)
`type: local` runs each job as a background process on the current machine
(`qxub.backends.local.LocalBackend`): `qx` streams and reports it as usual,
and `qxub status`, `monitor`, `cancel` and `interactive` work too. Resource
requests are not enforced and every job starts at once. Select it with:

```bash
export QXUB_PLATFORM_FILE=/path/to/qxub/docs/platforms/local/local.yaml
```

It lives in its own directory so that a search path pointing at
`docs/platforms` loads only real clusters.

### A Slurm cluster (`slurm/slurm.yaml`)
An example for `type: slurm` (`qxub.backends.slurm.SlurmBackend`). Copy it
and replace its partitions and limits with your cluster's. qxub's queues are
partitions and `--project/-P` is `--account`. Jobs are submitted with
`sbatch`, listed with `squeue`, looked up with `sacct` and cancelled with
`scancel`; interactive sessions use `srun --pty`. `-l name=value` requests
other than qxub's own become `--name=value` (`-l constraint=skylake`).

With no `--queue`, a job goes to the platform's `default_queue`, else to a
queue named `normal` if the platform has one, else to Slurm's default
partition. The status directory must be visible from the compute nodes:
set `paths.scratch`, or rely on `~/.local/state/qxub` where home is shared.

## Installation

### System-Level (Recommended)
```bash
# Install for all users (requires admin)
sudo mkdir -p /etc/qxub/platforms
sudo cp nci_gadi.yaml /etc/qxub/platforms/

# Verify installation
qxub platform list
```

### User-Level
```bash
# Install for current user only
mkdir -p ~/.config/qxub/platforms
cp nci_gadi.yaml ~/.config/qxub/platforms/

# Verify installation
qxub platform list
```

## Usage

Once installed, qxub will automatically:

1. **Select appropriate queues** based on resource requirements
2. **Validate resource requests** against queue limits
3. **Estimate Service Unit costs** before submission
4. **Apply default walltimes** per queue
5. **Suggest optimizations** for better resource utilization

### Examples

```bash
# GPU job automatically uses gpuvolta queue
qxub --env pytorch -l ngpus=1 -l ncpus=12 script.py

# High memory job uses hugemem queue
qxub --env analysis -l mem=500GB script.py

# Internet access uses copyq queue
qxub --env download --internet wget https://example.com/data.zip

# Cost estimation
qxub --env myenv -l ncpus=48 -l walltime=4:00:00 --estimate-cost script.py
```

## Platform Structure

Each platform definition includes:

- **Queue definitions** with accurate limits and constraints
- **Complex walltime rules** based on core count
- **Service Unit rates** for cost estimation
- **Auto-selection rules** for intelligent queue choice
- **Auto-adjustment policies** for resource optimization
- **Site conventions** that qxub would otherwise have to guess:

  ```yaml
  platform:
    host: gadi                 # matched case-insensitively against the hostname
    account_env: PROJECT       # env var holding the default project (-P)
    paths:
      scratch: "/scratch/{project}/{user}"   # default output and status files
    default_directives:        # added to every qsub
      - "-W umask=0027"
  ```

  qxub picks the platform whose `host` appears in the current hostname.
  Without a platform, qxub sends no `-P` unless one is configured, adds no
  directives, writes default output to `./qt/` and keeps status files in
  `$XDG_STATE_HOME/qxub/status`.

  Storage that your group always needs (for example `gdata/<project>`) is a
  group choice, not a platform one, so set it in config:
  `qxub config set defaults.volumes gdata/<project>`.

## Contributing

To add support for a new platform:

1. Create a new YAML file following the schema in `../dev/platform_schema.md`
2. Include comprehensive queue definitions with accurate limits
3. Add intelligent auto-selection rules
4. Test with realistic job scenarios
5. Document any platform-specific features

## Validation

Platform definitions are validated against the schema when loaded. Common issues:

- **Overlapping walltime rules**: Core ranges must not overlap
- **Inconsistent limits**: min_cpus must be ≤ max_cpus
- **Missing defaults**: Auto-selection rules must include a default
- **Invalid expressions**: Constraint expressions must be valid boolean logic

Use `qxub platform validate <platform>` to check your definitions.

## Updates

Platform definitions should be updated when:

- Queue limits change
- New queues are added/removed
- Service Unit rates change
- Scheduling policies are modified

Check the official documentation for your platform periodically for updates.
