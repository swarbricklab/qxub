# Multi-Platform Resource Architecture Design

## Overview

This document outlines the architectural design for extending qxub to support multiple job schedulers (PBS Pro, Slurm, SGE) and workflow engines (Snakemake, NextFlow, CWL) while maintaining backwards compatibility and enabling future cloud/container platform integration.

## Design Philosophy

### Core Principles

1. **No Premature Abstraction**: Start with native platform resources, add universal layer only if proven valuable
2. **Platform Fidelity**: Never lose platform-specific capabilities through abstraction
3. **Natural Evolution**: Architecture supports adding universal resources later without breaking changes
4. **Immediate Validation**: Catch impossible jobs before submission with helpful feedback
5. **Future Vision**: Design accommodates cloud platforms, job squashing, and infrastructure provisioning

### Key Architectural Decisions

- **Native Resources First**: Platforms use their natural resource syntax (e.g., PBS `-l mem=4GB`)
- **Adapter Pattern**: Workflow engines translate to current platform via adapters
- **Layered Approach**: Optional convenience layer can be added without changing core
- **Config-Driven Platform Detection**: Platform determined by config or CLI override
- **Validation at Platform Level**: Each platform validates its own resource constraints

## Current State Analysis

### Existing qxub Architecture
- **PBS Pro Only**: Current implementation is PBS-specific
- **Resource Handling**: Uses `-l/--resources` options that mirror `qsub` directly
- **Platform System**: YAML-based platform definitions exist but are PBS-centric
- **Configuration**: Hierarchical config system already supports platform detection

### Migration Strategy
- **Preserve Existing Interface**: All current PBS functionality continues to work
- **Extract Scheduler Interface**: Abstract current PBS functions behind interface
- **Maintain Backwards Compatibility**: Existing users see no changes

## Target Architecture

### Phase 1: Native Resources + Workflow Adapters

#### Platform Resource Interface
```python
# qxub/platforms/base.py
class PlatformResources:
    """Base class for platform-specific resource specifications."""

    @abstractmethod
    def to_command_args(self) -> List[str]:
        """Convert to scheduler command-line arguments."""
        pass

    @abstractmethod
    def validate(self, platform_limits: PlatformLimits) -> ValidationResult:
        """Validate against platform constraints."""
        pass

    @abstractmethod
    def estimate_cost(self) -> Optional[float]:
        """Estimate resource cost (SU, credits, etc.) if applicable."""
        pass

# qxub/platforms/pbs_pro.py
class PBSResources(PlatformResources):
    """PBS Pro native resource specification."""

    def __init__(self, **resources):
        # Direct PBS resource dict: mem=4GB, ncpus=2, walltime=2:00:00, etc.
        self.resources = resources

    def to_command_args(self) -> List[str]:
        """Convert to qsub -l arguments."""
        return [f"-l {k}={v}" for k, v in self.resources.items()]

    def validate(self, platform_limits: PlatformLimits) -> ValidationResult:
        """Validate against PBS platform constraints."""
        # Check memory limits, CPU limits, walltime limits per queue
        # Return ValidationResult with errors/warnings/suggestions
        pass

# qxub/platforms/slurm.py
class SlurmResources(PlatformResources):
    """Slurm native resource specification."""

    def __init__(self, mem=None, cpus=None, time=None, gres=None, **kwargs):
        # Native Slurm parameters
        self.mem = mem          # --mem=4G
        self.cpus = cpus        # --cpus-per-task=2
        self.time = time        # --time=02:00:00
        self.gres = gres        # --gres=gpu:2
        self.extras = kwargs    # Other Slurm options

    def to_command_args(self) -> List[str]:
        """Convert to sbatch arguments."""
        args = []
        if self.mem: args.append(f"--mem={self.mem}")
        if self.cpus: args.append(f"--cpus-per-task={self.cpus}")
        if self.time: args.append(f"--time={self.time}")
        if self.gres: args.append(f"--gres={self.gres}")
        return args
```

#### Workflow Engine Adapters
```python
# qxub/workflow_adapters/base.py
class WorkflowAdapter:
    """Base class for workflow engine resource translation."""

    @abstractmethod
    def translate_resources(self, platform: str, **workflow_resources) -> PlatformResources:
        """Translate workflow-specific resources to platform-native format."""
        pass

    @abstractmethod
    def supported_platforms(self) -> List[str]:
        """Return list of platforms this adapter supports."""
        pass

# qxub/workflow_adapters/snakemake.py
class SnakemakeAdapter(WorkflowAdapter):
    """Snakemake resource adapter."""

    def translate_resources(self, platform: str, **kwargs) -> PlatformResources:
        """
        Translate Snakemake resources to platform format.

        Snakemake resources: mem_mb, runtime (minutes), threads, disk_mb
        """
        if platform == "pbs":
            pbs_resources = {}

            if "mem_mb" in kwargs:
                pbs_resources["mem"] = f"{kwargs['mem_mb']}MB"

            if "runtime" in kwargs:
                # Convert minutes to HH:MM:SS
                minutes = kwargs["runtime"]
                hours = minutes // 60
                mins = minutes % 60
                pbs_resources["walltime"] = f"{hours:02d}:{mins:02d}:00"

            if "threads" in kwargs:
                pbs_resources["ncpus"] = kwargs["threads"]

            if "disk_mb" in kwargs:
                pbs_resources["jobfs"] = f"{kwargs['disk_mb']}MB"

            return PBSResources(**pbs_resources)

        elif platform == "slurm":
            return SlurmResources(
                mem=f"{kwargs.get('mem_mb', 1000)}M",
                cpus=kwargs.get("threads", 1),
                time=f"{kwargs.get('runtime', 60)}:00"  # minutes to MM:SS
            )

        else:
            raise UnsupportedPlatformError(f"Platform {platform} not supported")

    def supported_platforms(self) -> List[str]:
        return ["pbs", "slurm"]
```

#### Enhanced CLI Interface
```python
# qxub/exec_cli.py modifications
@click.command(name="exec")
# Existing PBS options (preserve backwards compatibility)
@click.option("-l", "--resources", multiple=True,
              help="Platform-native resource specification (e.g., 'walltime=1:00:00,mem=4GB')")

# New workflow engine options
@click.option("--workflow-engine", type=click.Choice(["snakemake", "nextflow", "cwl"]),
              help="Workflow engine resource format")
@click.option("--workflow-resources", multiple=True,
              help="Workflow-engine specific resources (e.g., 'mem_mb=4000,runtime=120')")

# Platform override
@click.option("--platform",
              help="Override platform detection (for remote execution)")

def exec_command(ctx, command, resources, workflow_engine, workflow_resources, platform, **options):
    """Enhanced exec command supporting multiple resource formats."""

    # Determine target platform
    target_platform = platform or config_manager.get_platform()

    # Process resources based on input format
    if workflow_engine and workflow_resources:
        # Use workflow adapter
        adapter = get_workflow_adapter(workflow_engine)
        workflow_kwargs = parse_workflow_resources(workflow_resources)
        platform_resources = adapter.translate_resources(target_platform, **workflow_kwargs)
    else:
        # Use native platform resources (existing behavior)
        platform_resources = create_native_resources(target_platform, resources)

    # Validate resources against platform
    platform_def = get_platform_definition(target_platform)
    validation = platform_resources.validate(platform_def.limits)

    if validation.errors:
        # Show validation errors with suggestions
        show_validation_errors(validation)
        ctx.exit(2)  # Validation error exit code

    if validation.warnings:
        show_validation_warnings(validation)

    # Continue with existing execution logic...
```

### Phase 2: Optional Convenience Layer

Add common shortcuts that work across platforms:

```python
# Optional convenience options (added later)
@click.option("--mem", help="Memory requirement (auto-translates to platform syntax)")
@click.option("--cpus", type=int, help="CPU cores (auto-translates to platform syntax)")
@click.option("--time", help="Walltime/runtime (auto-translates to platform syntax)")
@click.option("--gpus", type=int, help="GPU count (auto-translates to platform syntax)")
```

These would internally create the appropriate `PlatformResources` object for the current platform.

### Phase 3: Advanced Features (Future)

#### Job Squashing Framework
```python
# qxub/job_squashing.py
@dataclass
class JobBatch:
    """Collection of jobs that can be executed together."""
    jobs: List[Job]
    shared_resources: PlatformResources
    execution_strategy: str  # "parallel", "sequential", "array"

    def can_squash_with(self, other_job: Job) -> bool:
        """Check if another job can be squashed with this batch."""
        return (
            self.shared_resources.compatible_with(other_job.resources) and
            self.execution_strategy == other_job.preferred_execution_strategy
        )

class JobSquasher:
    """Intelligently batch compatible jobs for efficiency."""

    def squash_jobs(self, jobs: List[Job]) -> List[JobBatch]:
        """Group jobs into efficient batches."""
        # Algorithm:
        # 1. Group by squash_key (user-defined batching hint)
        # 2. Group by resource compatibility
        # 3. Group by execution context (conda env, modules, etc.)
        # 4. Respect walltime constraints
        pass
```

#### Cloud Platform Integration
```python
# qxub/platforms/aws_batch.py
class AWSBatchResources(PlatformResources):
    """AWS Batch resource specification."""

    def __init__(self, vcpus=None, memory=None, job_queue=None, **kwargs):
        self.vcpus = vcpus
        self.memory = memory      # In MB
        self.job_queue = job_queue
        self.extras = kwargs

    def provision_infrastructure(self) -> bool:
        """Create compute environment if needed."""
        # Integration with Terraform/CloudFormation
        pass

# qxub/platforms/kubernetes.py
class KubernetesResources(PlatformResources):
    """Kubernetes Job resource specification."""

    def __init__(self, memory=None, cpu=None, gpu=None, **kwargs):
        self.memory = memory     # "4Gi"
        self.cpu = cpu          # "2000m"
        self.gpu = gpu          # {"nvidia.com/gpu": 2}
        self.extras = kwargs
```

## Implementation Roadmap

### Immediate (Phase 1)

1. **Extract Platform Interface**
   - Create `qxub/platforms/base.py` with `PlatformResources` base class
   - Implement `PBSResources` by refactoring existing PBS code
   - Add validation framework with helpful error messages

2. **Create Workflow Adapter System**
   - Implement `qxub/workflow_adapters/base.py`
   - Create `SnakemakeAdapter` as proof of concept
   - Add adapter registry and discovery mechanism

3. **Enhance CLI Integration**
   - Add `--workflow-engine` and `--workflow-resources` options
   - Implement platform detection from config
   - Add comprehensive validation with suggestions

4. **Add Slurm Support**
   - Implement `SlurmResources` class
   - Extend `SnakemakeAdapter` to support Slurm
   - Test cross-platform workflow compatibility

### Near-term (Phase 2)

5. **Add Convenience Layer**
   - Implement `--mem`, `--cpus`, `--time` shortcuts
   - Create auto-translation to current platform
   - Maintain backwards compatibility

6. **Enhanced Validation**
   - Platform-aware resource constraint checking
   - Cost estimation and optimization suggestions
   - Queue selection based on resource requirements

### Long-term (Phase 3+)

7. **Job Squashing**
   - Implement job batching algorithms
   - Create user-facing job tracking that hides batching
   - Integration with monitoring and status systems

8. **Cloud Platform Support**
   - Kubernetes Jobs integration
   - AWS Batch platform implementation
   - Infrastructure provisioning framework

9. **Multi-Cloud Orchestration**
   - Cost-based platform selection
   - Auto-scaling and resource optimization
   - Hybrid cloud/HPC workflows

## Configuration Schema

### Platform Configuration Extension
```yaml
# ~/.config/qxub/config.yaml
platform:
  default: "nci_gadi"  # Default platform for local execution

  # Platform definitions can specify scheduler type
  remote_platforms:
    cluster2:
      platform: "slurm_cluster"
      scheduler: "slurm"
    aws_dev:
      platform: "aws_batch"
      scheduler: "aws_batch"

# Workflow engine preferences
workflow_engines:
  snakemake:
    default_platform: "nci_gadi"
    resource_multipliers:
      memory: 1.1  # Add 10% memory buffer
      walltime: 1.2  # Add 20% time buffer
```

### Platform Definition Schema Extension
```yaml
# docs/platforms/slurm_example.yaml
platform:
  name: slurm_cluster
  scheduler: slurm
  host: "cluster.example.edu"

  # Scheduler-specific configuration
  scheduler_config:
    default_partition: "normal"
    account_required: true
    time_format: "HH:MM:SS"

  # Resource mappings for convenience layer
  resource_mappings:
    memory: "--mem"
    cpus: "--cpus-per-task"
    walltime: "--time"
    gpus: "--gres=gpu"

  queues:
    - name: normal
      limits:
        max_cpus: 128
        max_memory: "512GB"
        max_walltime: "7-00:00:00"
```

## Testing Strategy

### Unit Tests
- Test each `PlatformResources` implementation
- Test workflow adapters with various input combinations
- Test validation logic with edge cases

### Integration Tests
- Test cross-platform workflow execution
- Test resource translation accuracy
- Test validation error messages and suggestions

### End-to-End Tests
- Real Snakemake workflows on PBS and Slurm
- Resource constraint validation on real platforms
- Cost estimation accuracy

## Migration Plan

### Backwards Compatibility
- All existing qxub commands continue to work unchanged
- Existing PBS `-l` options preserved exactly
- Configuration files remain compatible

### User Migration
- Document workflow adapter usage patterns
- Provide migration examples for common use cases
- Create platform-specific quick-start guides

### Developer Migration
- Clear interfaces for adding new platforms
- Template implementations for common schedulers
- Plugin architecture documentation

## Success Metrics

### Phase 1 Success Criteria
- [ ] Snakemake workflows work on both PBS and Slurm
- [ ] Resource validation catches impossible jobs before submission
- [ ] Zero breaking changes for existing PBS users
- [ ] Clear error messages with actionable suggestions

### Long-term Vision Success
- [ ] Support for 5+ job schedulers (PBS, Slurm, SGE, K8s, AWS Batch)
- [ ] Support for 3+ workflow engines (Snakemake, NextFlow, CWL)
- [ ] Job squashing reduces scheduler load by 10x for array-like workloads
- [ ] Cloud integration enables cost-optimized hybrid workflows

This architecture provides a clear path from the current PBS-only implementation to a universal job dispatch system while maintaining qxub's core value proposition of eliminating boilerplate and providing intelligent defaults.
