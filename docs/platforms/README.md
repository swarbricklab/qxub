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