# GitHub Issue: Platform-aware environment discovery

**Title**: Platform-aware environment discovery (modules, conda, singularity)
**Labels**: enhancement, v2.1, platform-abstraction
**Priority**: Low

## Overview

Add platform-specific environment discovery capabilities to qxub to enable intelligent module, conda, and singularity image management across different HPC systems.

## Proposed Platform Properties

Add these three properties to platform definitions:

1. **`modulepaths`** - List of paths to search for environment modules
2. **`condarc`** - List of conda configuration files to source  
3. **`simg`** - List of directories containing Singularity images

## Example Configuration

```yaml
platform:
  name: nci_gadi
  modulepaths:
    - "/apps/Modules/modulefiles"
    - "/g/data/hh5/public/modules"
    - "/opt/nci/modulefiles"
  
  condarc:
    - "/g/data/hh5/public/apps/miniconda3/etc/conda/condarc"
    - "/apps/conda/etc/conda/condarc"
    - "~/.condarc"
  
  simg:
    - "/g/data/hh5/public/apps/singularity/containers"
    - "/apps/singularity/containers"
    - "/scratch/$PROJECT/singularity"
```

## Proposed Functionality

### Module Discovery (`--mod`)
- Search `modulepaths` to validate modules exist before job submission
- Provide tab completion for available modules
- Generate appropriate `module load` commands
- Enable `qxub modules list` command

### Conda Environment Discovery (`--env`)  
- Source `condarc` files to discover available environments
- Find environments across different conda installations
- Validate environments exist before submission
- Enable `qxub envs list` command

### Singularity Image Discovery (`--sif`)
- Search `simg` directories for .sif/.simg/.sqsh files
- Support both full paths and name-based resolution
- Provide tab completion for available containers
- Generate appropriate `singularity exec` commands

## Benefits

1. **Environment Validation**: Check environments exist before submitting jobs
2. **Tab Completion**: Auto-complete module/environment/container names
3. **Cross-Platform Discovery**: Find environments across different installations  
4. **User Experience**: No need to remember exact paths or names
5. **Error Prevention**: Catch missing environments at submission time

## Implementation Notes

- This would be part of the v2.1+ platform abstraction system
- Should integrate with the existing execution context options (`--env`, `--mod`, `--sif`)
- Could provide fallback behavior when platform definitions don't include these paths
- Tab completion would require shell integration

## Related

- Platform schema design: docs/dev/platform_schema.md
- v2.1 roadmap: docs/dev/v2_roadmap.md

## Priority

Low priority - this is a nice-to-have enhancement that would improve user experience but isn't critical for core functionality.

---

**Note**: This issue should be created in the GitHub repository when the MCP server is available.