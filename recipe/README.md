# qxub Conda Recipe

This directory contains the conda recipe for building qxub as a conda package.

## Building the Package

### Prerequisites

- conda-build installed: `conda install conda-build`
- conda-forge channel enabled (recommended)

### Build Locally

```bash
# From the repository root
conda build recipe/

# Or specify the recipe directory explicitly
conda build /g/data/a56/software/qsub_tools/recipe
```

### Build Output

The built package will be in your conda-bld directory:
- Linux: `~/miniconda3/conda-bld/noarch/qxub-3.3.1-py_0.tar.bz2`
- Location varies by platform and conda installation

### Installing from Local Build

```bash
# Install from local build
conda install --use-local qxub

# Or specify the full path
conda install ~/miniconda3/conda-bld/noarch/qxub-3.3.1-py_0.tar.bz2
```

## Testing the Package

After building, conda will automatically run the tests defined in `meta.yaml`:
- Import tests for all major modules
- Command availability tests for all CLI entry points

To test manually after installation:

```bash
# Activate environment and test commands
conda activate test_env
qxub --version
qxub --help
qxub exec --help
qx --help
```

## Publishing to Conda-Forge

To publish qxub to conda-forge:

1. Fork https://github.com/conda-forge/staged-recipes
2. Add this recipe to `recipes/qxub/`
3. Submit a pull request
4. Address feedback from conda-forge reviewers
5. Once merged, a feedstock will be created at https://github.com/conda-forge/qxub-feedstock

## Publishing to a Custom Channel

```bash
# Build the package
conda build recipe/

# Upload to anaconda.org
anaconda upload ~/miniconda3/conda-bld/noarch/qxub-3.3.1-py_0.tar.bz2

# Users can then install with:
# conda install -c YOUR_CHANNEL qxub
```

## Recipe Maintenance

When updating the version:
1. Update `version` in `recipe/meta.yaml`
2. Update `version` in `qxub/__init__.py`
3. Update `version` in `setup.py`
4. Rebuild and test

## Notes

- This is a `noarch: python` package (platform-independent)
- Requires Python >=3.6
- All dependencies are available on conda-forge
- Entry points include: qxub, qx, qxtat, qxet
