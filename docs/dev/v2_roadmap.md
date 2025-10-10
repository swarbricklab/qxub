# qxub v2.1+ Development Roadmap

## Overview

This document outlines the development roadmap for qxub v2.1 and beyond, focusing on intelligent queue selection and platform abstraction.

## Vision Statement

Transform qxub from a PBS-specific tool into a universal computational job orchestrator while maintaining simplicity and backward compatibility. Enable seamless execution across platforms (HPC, cloud, local) with intelligent resource management and remote execution capabilities.

## Development Phases

### Phase 2.1: Intelligent Queue Selection 🎯 **Priority**
**Goal**: Smart queue selection and resource validation for NCI Gadi

#### Core Features
- **Queue constraint validation** against platform definitions
- **Automatic queue selection** based on resource requirements
- **Resource auto-adjustment** to meet queue minimums
- **Enhanced error messages** with suggestions
- **Platform schema framework** for future extension

#### Implementation Tasks
1. **Platform schema parser** (`qxub/platforms/`)
2. **NCI Gadi platform plugin** with queue definitions
3. **Queue validation engine**
4. **Auto-selection algorithm**
5. **Resource adjustment policies**
6. **Enhanced CLI feedback**

#### Success Criteria
- Users can specify `--queue auto` for intelligent selection
- GPU requests automatically select `gpuvolta` and adjust CPU minimum
- Memory >192GB requests automatically suggest `hugemem` queue
- Clear warnings/suggestions for constraint violations
- 100% backward compatibility with existing commands

#### Timeline: 4-6 weeks

---

### Phase 2.2: Remote Execution 🌐
**Goal**: Execute jobs on remote platforms from local machines

#### Core Features
- **Client-server architecture** for remote job submission
- **SSH-based remote execution** with credential delegation
- **Platform profiles** combining host + platform + credentials
- **Distributed logging** (client + server)
- **Exit code propagation** through the execution chain

#### Implementation Tasks
1. **Platform profiles** with remote host configuration
2. **SSH execution backend**
3. **Remote qxub invocation** and monitoring
4. **Output streaming** from remote to local
5. **Distributed history logging**

#### Usage Example
```bash
qxub --profile gadi --env dvc3 -- dvc --version
```

#### Success Criteria
- Execute jobs on NCI from laptop with same syntax
- Real-time output streaming to local terminal
- Proper exit code handling
- Execution logged both locally and remotely

#### Timeline: 6-8 weeks

---

### Phase 3.0: Multi-Platform Orchestration 🚀
**Goal**: Universal job orchestration across platforms

#### Core Features
- **Plugin architecture** for multiple platforms
- **Abstract resource specifications** (CPU, memory, GPU)
- **Cross-platform resource translation**
- **Cloud platform support** (AWS Batch, GCP, K8s)
- **Workflow integration** (DVC, Snakemake profiles)

#### Implementation Strategy
- **Core resource abstraction** layer
- **Platform-specific translators**
- **Unified configuration** system
- **Advanced scheduling** and cost optimization

#### Timeline: 3-6 months (major version)

## Technical Architecture

### Platform Plugin System

```python
# qxub/platforms/base.py
class Platform(ABC):
    @abstractmethod
    def validate_resources(self, request: ResourceRequest) -> ValidationResult:
        """Validate resource request against platform constraints"""

    @abstractmethod
    def suggest_queue(self, request: ResourceRequest) -> Optional[str]:
        """Suggest appropriate queue for resource request"""

    @abstractmethod
    def adjust_resources(self, request: ResourceRequest, queue: str) -> ResourceRequest:
        """Auto-adjust resources to meet queue requirements"""

# qxub/platforms/nci_gadi.py
class NCIGadiPlatform(Platform):
    def __init__(self, config_path: Path):
        self.queues = self._load_queue_definitions(config_path)

    def validate_resources(self, request):
        # Validate against NCI queue constraints

    def suggest_queue(self, request):
        # Apply NCI-specific queue selection rules
```

### Resource Request Processing

```python
# Flow: User Input → Validation → Auto-Selection → Auto-Adjustment → Execution

def process_resource_request(user_args, config):
    # 1. Parse user resource specification
    request = parse_user_resources(user_args)

    # 2. Load platform and validate
    platform = load_platform(config.platform)
    validation = platform.validate_resources(request)

    # 3. Auto-select queue if needed
    if request.queue == "auto":
        request.queue = platform.suggest_queue(request)

    # 4. Auto-adjust resources if configured
    if config.auto_adjust_enabled:
        request = platform.adjust_resources(request, request.queue)

    # 5. Final validation and warnings
    final_validation = platform.validate_resources(request)

    return request, final_validation
```

## Integration Points

### DVC Integration
```yaml
# dvc.yaml
stages:
  train:
    cmd: python train.py
    qxub:
      env: pytorch
      resources:
        memory: 32GB
        gpu: 1
        walltime: 4h
```

### Snakemake Integration
```python
# Snakemake profile using qxub
def qxub_submit(job, **kwargs):
    return qxub_platform.submit_job(
        command=job.format_command(),
        resources=job.resources,
        env=job.conda_env
    )
```

## Configuration Management

### Platform Definitions
- **Location**: `~/.config/qxub/platforms/`
- **Format**: YAML with schema validation
- **Versioning**: Platform schema version tracking
- **Updates**: Automatic/manual platform definition updates

### User Configuration Evolution
```yaml
# v2.0 (current)
defaults:
  queue: normal

# v2.1 (intelligent)
defaults:
  queue: auto
auto_selection:
  enabled: true

# v2.2 (remote)
profiles:
  gadi:
    platform: nci_gadi
    remote:
      host: gadi.nci.org.au
```

## Backward Compatibility Strategy

### Guaranteed Compatibility
1. **All v2.0 commands work unchanged**
2. **Existing configuration files work unchanged**
3. **Explicit options always override auto-selection**
4. **Platform abstraction is invisible when not used**

### Migration Path
- **v2.0 → v2.1**: Zero-effort migration, new features opt-in
- **v2.1 → v2.2**: Add profiles for remote execution
- **v2.2 → v3.0**: Enhanced platform support, optional migration

## Success Metrics

### v2.1 Success
- [ ] 100% NCI queue constraints modeled accurately
- [ ] Auto-selection reduces user queue selection errors by >80%
- [ ] Resource suggestions accepted by users >60% of time
- [ ] Zero regression in existing functionality
- [ ] Positive feedback from 3 active users

### v2.2 Success
- [ ] Remote execution works transparently from laptop
- [ ] <5 second latency for remote job submission
- [ ] Exit codes propagated correctly 100% of time
- [ ] Distributed logging captures all relevant events

### v3.0 Success
- [ ] Support for 3+ platforms (NCI + 2 cloud providers)
- [ ] DVC/Snakemake integration working
- [ ] Community adoption beyond initial 3 users
- [ ] Performance equivalent to platform-native tools

## Risk Mitigation

### Technical Risks
- **Complexity creep**: Maintain simple CLI surface
- **Performance regression**: Benchmark all changes
- **Platform divergence**: Keep abstraction layer thin

### User Experience Risks
- **Configuration complexity**: Provide smart defaults
- **Backward compatibility**: Maintain comprehensive test suite
- **Learning curve**: Extensive documentation and examples

### Development Risks
- **Scope creep**: Stick to phased delivery
- **Resource constraints**: Focus on high-impact features first
- **Platform changes**: Design for platform evolution

## Next Steps

### Immediate (Week 1-2)
1. **Finalize platform schema** based on discussion
2. **Validate NCI queue constraints** against documentation
3. **Create basic platform parser** and NCI plugin
4. **Begin queue validation implementation**

### Short-term (Month 1)
1. **Complete v2.1 core features**
2. **Test with existing users**
3. **Refine based on feedback**
4. **Plan v2.2 detailed design**

### Medium-term (Month 2-3)
1. **Implement remote execution**
2. **Design multi-platform architecture**
3. **Begin v3.0 planning**
4. **Community feedback and adoption**
