# qxub 2.0 Migration Roadmap

## Overview

qxub 2.0 will eliminate the `conda`, `module`, and `sing` subcommands in favor of a unified CLI with execution context options. This represents a major architectural simplification that will significantly improve user experience.

## Current vs Proposed Architecture

### Current (1.1.x)
```bash
qxub conda --env myenv python script.py
qxub module --mod python3,gcc python script.py  
qxub sing --sif container.sif python script.py
```

### Proposed (2.0)
```bash
qxub --env myenv -- python script.py          # or --conda myenv
qxub --mod python3,gcc -- python script.py   # or --module python3,gcc
qxub --sif container.sif -- python script.py # or --sing / --singularity
```

**Alternative Options**: For better usability, multiple option names are supported:
- `--env` / `--conda`: Conda environment execution
- `--mod` / `--module`: Environment module execution  
- `--sif` / `--sing` / `--singularity`: Singularity container execution

**Note**: The `--` separator is required to clearly separate qxub options from target command options.

## Breaking Changes

### 1. Subcommand Removal
- **Impact**: All existing scripts using `qxub conda`, `qxub module`, `qxub sing` will break
- **Migration**: Users must remove subcommands and move options to main command level
- **Timeline**: 2.0.0 release (clean break, no backward compatibility)

### 2. Alias Structure Changes
- **Current**: Hierarchical structure with `subcommand.type`
- **Proposed**: Flat structure with execution context options
- **Migration**: Automatic alias conversion tool needed

### 3. Option Placement with `--` Separator
- **Current**: Complex rules about main vs subcommand options
- **Proposed**: All qxub options before `--`, all target command options after `--`
- **Benefit**: Unambiguous parsing, eliminates option placement confusion entirely

## Implementation Phases

### Phase 1: Foundation (2.0.0-alpha.1)
**Goal**: Core CLI restructuring

#### CLI Architecture Changes
- [ ] Remove subcommand definitions from `cli.py`
- [ ] Add alternative option names to main command level:
  - [ ] `--env` / `--conda` for conda environments
  - [ ] `--mod` / `--module` / `--mods` for environment modules
  - [ ] `--sif` / `--sing` / `--singularity` for containers
- [ ] Implement option consolidation logic for alternatives
- [ ] Add mutual exclusivity validation between execution contexts
- [ ] Update help system for unified command structure with alternative options

#### Core Execution Logic
- [ ] Modify `scheduler.py` to detect execution context from main options
- [ ] Update job script generation to work without subcommand context
- [ ] Ensure resource tracking continues to work with new structure

#### Testing Infrastructure
- [ ] Update all test cases to use new syntax
- [ ] Validate that existing functionality remains intact
- [ ] Test mutual exclusivity enforcement

### Phase 2: Alias System Migration (2.0.0-alpha.2)
**Goal**: Seamless alias transition

#### Alias Structure Updates
- [ ] Design new flat alias structure
- [ ] Implement automatic conversion from hierarchical to flat format
- [ ] Update alias creation/editing commands
- [ ] Maintain backward compatibility for existing alias files during transition

#### Migration Tools
- [ ] Create `qxub migrate aliases` command
- [ ] Automatic detection and conversion of old alias format
- [ ] Backup mechanism for existing alias configurations
- [ ] Validation tool to ensure converted aliases work correctly

#### Template System Updates
- [ ] Ensure template variables work with new structure
- [ ] Update template resolution logic
- [ ] Test complex template scenarios

### Phase 3: Documentation and UX (2.0.0-beta.1)
**Goal**: Complete user experience overhaul

#### Documentation Overhaul
- [ ] Rewrite all CLI examples to use new syntax
- [ ] Update alias documentation with new structure
- [ ] Create migration guide for existing users
- [ ] Update troubleshooting guides

#### Enhanced Error Handling
- [ ] Simplify error handling (no more subcommand confusion)
- [ ] Improve mutual exclusivity error messages
- [ ] Add suggestions for common migration mistakes
- [ ] Update fuzzy matching for new option structure

#### Help System Improvements
- [ ] Reorganize help output for unified command
- [ ] Group related options together logically
- [ ] Improve examples and usage patterns
- [ ] Add quick migration tips to help text

### Phase 4: Advanced Features (2.0.0-beta.2)
**Goal**: Enhanced functionality and performance

#### Default Execution Context
- [ ] Implement intelligent default detection
- [ ] Add configuration for preferred execution context
- [ ] Support environment variable hints
- [ ] Graceful handling when no context specified

#### History System Updates
- [ ] Update recipe hashing for new command structure
- [ ] Ensure history compatibility across versions
- [ ] Update conversion tools for history migration
- [ ] Test recipe-to-alias conversion with new format

#### System Integration
- [ ] Update system configuration examples
- [ ] Ensure organizational defaults work with new structure
- [ ] Test hierarchical config inheritance
- [ ] Validate XDG compliance maintained

### Phase 5: Finalization (2.0.0-rc.1)
**Goal**: Production readiness

#### Performance Optimization
- [ ] Profile new command parsing performance
- [ ] Optimize option validation logic
- [ ] Ensure no regression in job submission speed
- [ ] Memory usage analysis

#### Edge Case Handling
- [ ] Test complex option combinations
- [ ] Validate resource specification edge cases
- [ ] Test with various PBS environments
- [ ] Ensure compatibility with different shell environments

#### Final Documentation
- [ ] Complete migration guide
- [ ] Update all examples in README
- [ ] Create video/visual guides for migration
- [ ] Prepare release notes

## Migration Strategy

### For Users

#### Automatic Migration Tools
```bash
# Check current usage patterns
qxub migrate analyze

# Convert aliases automatically
qxub migrate aliases --backup

# Validate migrated configuration
qxub migrate validate

# Get migration suggestions for scripts
qxub migrate scripts /path/to/scripts/
```

#### Manual Migration Examples
```bash
# Before (1.x)
qxub conda --env myenv --queue normal python script.py

# After (2.0) - Multiple option styles supported
qxub --env myenv --queue normal -- python script.py
qxub --conda myenv --queue normal -- python script.py  # Alternative

# Before (1.x)
qxub module --mods python3,gcc --resources mem=8GB make

# After (2.0) - Multiple option styles supported
qxub --mods python3,gcc --resources mem=8GB -- make
qxub --mod python3,gcc --resources mem=8GB -- make      # Alternative
qxub --module python3,gcc --resources mem=8GB -- make   # Alternative

# Before (1.x)
qxub sing --sif /path/to/container.sif --bind /data python

# After (2.0) - Multiple option styles supported
qxub --sif /path/to/container.sif --bind /data -- python
qxub --sing /path/to/container.sif --bind /data -- python        # Alternative
qxub --singularity /path/to/container.sif --bind /data -- python # Alternative

# Examples with target command options
# Before (1.x)
qxub conda --env dvc3 dvc repro --dry-run

# After (2.0) - Natural, intuitive syntax
qxub --conda dvc3 -- dvc repro --dry-run

# Before (1.x) 
qxub module --mods python3 python train.py --epochs 100 --lr 0.001

# After (2.0) - Clear, explicit option names
qxub --module python3 -- python train.py --epochs 100 --lr 0.001

# Real-world examples showing natural feel
qxub --conda pytorch -- python train.py --gpu
qxub --module samtools,bwa -- bwa mem -t 16 ref.fa reads.fastq
qxub --singularity /containers/blast.sif -- blastn -query input.fa
```

### For System Administrators

#### Organization-wide Migration
- [ ] Communication plan for breaking changes
- [ ] Centralized script migration tools
- [ ] System configuration updates
- [ ] User training materials

#### Rollout Strategy
- [ ] Alpha testing with developer team
- [ ] Beta testing with power users
- [ ] Gradual production rollout
- [ ] Fallback plan to 1.x if needed

## Technical Implementation Details

### CLI Structure Changes

#### New Main Command Options
```python
@click.command()
@click.option('--env', '--conda', help='Conda environment for execution')
@click.option('--mod', '--module', '--mods', help='Environment modules to load')
@click.option('--sif', '--sing', '--singularity', help='Singularity container image')
@click.option('--bind', help='Singularity bind mounts')
# ... all existing PBS options remain at main level
@click.argument('command', nargs=-1)
def qxub(env, conda, mod, module, mods, sif, sing, singularity, bind, command, **kwargs):
    # Consolidate alternative option names
    conda_env = env or conda
    modules = mod or module or mods
    container = sif or sing or singularity
    
    # Validate mutual exclusivity
    execution_contexts = [conda_env, modules, container]
    if sum(bool(x) for x in execution_contexts) > 1:
        raise click.ClickException("Cannot specify multiple execution contexts")
    
    # Detect execution context and delegate
    if conda_env:
        execute_conda(conda_env, command, **kwargs)
    elif modules:
        execute_module(modules, command, **kwargs)
    elif container:
        execute_singularity(container, bind, command, **kwargs)
    else:
        # Default behavior - needs definition
        handle_no_execution_context(command, **kwargs)
```

#### Alias Structure Evolution
```yaml
# Current (1.x)
aliases:
  my_job:
    main:
      queue: normal
      resources: ["mem=8GB"]
    subcommand:
      type: conda
      env: myenv
    target:
      cmd: "python script.py"

# New (2.0)
aliases:
  my_job:
    queue: normal
    resources: ["mem=8GB"]
    env: myenv
    cmd: "python script.py"
```

### Alternative Option Names
Providing multiple option names for the same functionality improves usability:

#### Benefits
- **Intuitive**: `--conda` clearly indicates conda execution
- **Explicit**: `--module` is more descriptive than `--mod`
- **Familiar**: `--singularity` matches the tool name
- **Flexible**: Users can choose their preferred style
- **Migration-friendly**: Easier transition from subcommand names

#### Implementation
Click's multiple option names feature handles this elegantly:
```python
@click.option('--env', '--conda', help='Conda environment')
@click.option('--mod', '--module', '--mods', help='Environment modules')
@click.option('--sif', '--sing', '--singularity', help='Singularity container')
```

### Validation Logic
- **Mutual Exclusivity**: Ensure only one execution context is specified
- **Alternative Consolidation**: Merge alternative option names into single variables
- **Context Requirements**: Validate required options for each context
- **Resource Compatibility**: Ensure resources make sense for execution context

### Backward Compatibility Strategy
- **None for CLI**: Clean break in 2.0
- **Alias Migration**: Automatic conversion with backup
- **Configuration**: Maintain existing config structure where possible
- **History**: Maintain compatibility with existing history data

## Testing Strategy

### Automated Testing
- [ ] Unit tests for new CLI parsing
- [ ] Integration tests for all execution contexts
- [ ] Migration tool testing
- [ ] Performance regression testing

### User Acceptance Testing
- [ ] Power user testing with real workflows
- [ ] Documentation walkthrough testing
- [ ] Migration tool usability testing
- [ ] Error message clarity testing

### Compatibility Testing
- [ ] Various PBS environments
- [ ] Different shell environments
- [ ] System configuration scenarios
- [ ] Complex alias patterns

## Risk Mitigation

### Major Risks
1. **User Adoption Resistance**: Breaking changes always cause friction
2. **Migration Complexity**: Complex aliases might not migrate cleanly
3. **Script Incompatibility**: Widespread script breakage
4. **Performance Regression**: New parsing logic slower than subcommands

### Mitigation Strategies
1. **Communication**: Early and frequent communication about benefits
2. **Tools**: Robust automatic migration tools
3. **Support**: Comprehensive migration guides and support
4. **Testing**: Extensive testing before release

### Rollback Plan
- Maintain 1.x branch for critical fixes
- Clear downgrade path for organizations
- Documentation for reverting migrations
- Support timeline for 1.x (6 months minimum)

## Success Metrics

### User Experience
- [ ] Reduced CLI error rates
- [ ] Shorter command lengths
- [ ] Faster new user onboarding
- [ ] Positive user feedback

### Technical Metrics
- [ ] No performance regression
- [ ] Successful alias migration rate >95%
- [ ] Documentation completeness score
- [ ] Test coverage maintenance

### Adoption Metrics
- [ ] Migration tool usage rates
- [ ] Support ticket volume (should decrease)
- [ ] Feature request patterns
- [ ] Community engagement

## Timeline

### Q4 2025
- **2.0.0-alpha.1**: Core CLI restructuring (Week 1-2)
- **2.0.0-alpha.2**: Alias migration tools (Week 3-4)

### Q1 2026
- **2.0.0-beta.1**: Documentation and UX (Week 1-2)
- **2.0.0-beta.2**: Advanced features (Week 3-4)

### Q2 2026
- **2.0.0-rc.1**: Release candidate (Week 1-2)
- **2.0.0**: Final release (Week 3-4)

## Communication Plan

### Developer Community
- [ ] GitHub issue for discussion
- [ ] RFC process for major decisions
- [ ] Regular progress updates
- [ ] Beta testing coordination

### End Users
- [ ] Migration announcement
- [ ] Preview releases with feedback collection
- [ ] Training materials and webinars
- [ ] Gradual rollout communication

### System Administrators
- [ ] Advanced notice for infrastructure planning
- [ ] Migration tools and documentation
- [ ] Support for organizational rollouts
- [ ] Coordination with IT teams

---

## Next Steps

1. **Validate Approach**: Get feedback on this roadmap
2. **Start Phase 1**: Begin core CLI restructuring
3. **Create Migration Tools**: Build automatic conversion utilities
4. **User Testing**: Early testing with willing power users
5. **Iterative Refinement**: Adjust based on testing feedback

This migration represents a significant improvement in qxub's usability and maintainability, eliminating a major source of user confusion while preserving all existing functionality.