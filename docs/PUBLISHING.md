# Publishing qxub - Complete Guide

This guide covers all publishing channels for qxub when you're ready to release.

## Table of Contents

1. [Pre-Release Checklist](#pre-release-checklist)
2. [GitHub Release](#github-release)
3. [PyPI (Python Package Index)](#pypi-python-package-index)
4. [Conda-Forge](#conda-forge)
5. [Version Management](#version-management)
6. [Post-Release Tasks](#post-release-tasks)

---

## Pre-Release Checklist

Before publishing any release:

- [ ] All tests pass locally and in CI
- [ ] Version bumped in **both** `qxub/__init__.py` and `setup.py`
- [ ] CHANGELOG or release notes updated
- [ ] Documentation is current
- [ ] Branch merged to `main` (or release branch is ready)
- [ ] No outstanding critical issues
- [ ] Git tag created with version (e.g., `v3.3.1`)

---

## GitHub Release

GitHub releases provide downloadable source archives and release notes.

### 1. Create and Push Git Tag

```bash
# Ensure you're on the release branch or main
git checkout main  # or release/v3.3.1

# Create annotated tag
git tag -a v3.3.1 -m "Release v3.3.1: Conda recipe and polish"

# Push tag to GitHub
git push origin v3.3.1
```

### 2. Create Release on GitHub

```bash
gh release create v3.3.1 \
  --title "v3.3.1 - Conda recipe and polish" \
  --notes "See CHANGELOG.md for details"
```

**Release Notes Template:**

```markdown
# v3.3.1 - Conda Recipe and Polish

## 🎁 New Features

- **Conda Recipe**: Official conda recipe for easy installation via conda/mamba
- **MIT License**: Added explicit license file

## 📦 Installation

### pip
\`\`\`bash
pip install git+https://github.com/swarbricklab/qxub.git@v3.3.1
\`\`\`

### conda (local build)
\`\`\`bash
git clone https://github.com/swarbricklab/qxub.git
cd qxub
git checkout v3.3.1
conda build recipe/
conda install --use-local qxub
\`\`\`

## 📝 Changes

- Created comprehensive conda recipe in `recipe/` directory
- Added MIT LICENSE file
- Recipe includes all CLI entry points and comprehensive tests

## 🔗 Links

- [Full Changelog](https://github.com/swarbricklab/qxub/compare/v3.3.0...v3.3.1)
- [Documentation](https://github.com/swarbricklab/qxub/tree/v3.3.1/docs)
```

---

## PyPI (Python Package Index)

PyPI allows users to install with `pip install qxub`.

### 1. Prerequisites

```bash
# Install build and publish tools
pip install build twine

# Get PyPI credentials ready
# Register at https://pypi.org if you haven't
# Set up API token at https://pypi.org/manage/account/token/
```

### 2. Build Distribution Packages

```bash
# From repository root
python -m build

# This creates:
# dist/qxub-3.3.1-py3-none-any.whl
# dist/qxub-3.3.1.tar.gz
```

### 3. Test on TestPyPI (Optional but Recommended)

```bash
# Upload to test instance
twine upload --repository testpypi dist/*

# Test installation
pip install --index-url https://test.pypi.org/simple/ qxub
```

### 4. Upload to PyPI

```bash
# Upload to production PyPI
twine upload dist/*

# Enter your API token when prompted
# Or use token via environment:
# export TWINE_USERNAME=__token__
# export TWINE_PASSWORD=pypi-...your-token...
```

### 5. Verify Installation

```bash
pip install qxub==3.3.1
qxub --version
```

### Using GitHub Actions (Automated)

Create `.github/workflows/publish-pypi.yml`:

```yaml
name: Publish to PyPI

on:
  release:
    types: [published]

jobs:
  publish:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.10'

      - name: Install build dependencies
        run: |
          python -m pip install --upgrade pip
          pip install build twine

      - name: Build package
        run: python -m build

      - name: Publish to PyPI
        env:
          TWINE_USERNAME: __token__
          TWINE_PASSWORD: ${{ secrets.PYPI_API_TOKEN }}
        run: twine upload dist/*
```

**Setup:**
1. Get PyPI API token from https://pypi.org/manage/account/token/
2. Add as GitHub secret: Settings → Secrets → Actions → `PYPI_API_TOKEN`
3. Push a release tag - workflow runs automatically

---

## Conda-Forge

Conda-Forge is the community-driven conda channel. Publishing here makes `conda install qxub` work.

### Method 1: Submit to Conda-Forge (Recommended for Public Release)

**Step 1: Fork staged-recipes**

```bash
git clone https://github.com/YOUR_USERNAME/staged-recipes.git
cd staged-recipes
```

**Step 2: Add recipe**

```bash
# Create recipe directory
mkdir recipes/qxub

# Copy your recipe files
cp /path/to/qxub/recipe/meta.yaml recipes/qxub/
cp /path/to/qxub/recipe/build.sh recipes/qxub/
cp /path/to/qxub/recipe/bld.bat recipes/qxub/

# Important: Update source in meta.yaml to use GitHub release
```

**Step 3: Update meta.yaml for conda-forge**

Change the source section from local path to GitHub:

```yaml
source:
  url: https://github.com/swarbricklab/qxub/archive/refs/tags/v{{ version }}.tar.gz
  sha256: <calculate with: curl -sL URL | sha256sum>
```

**Step 4: Submit Pull Request**

```bash
git checkout -b add-qxub
git add recipes/qxub/
git commit -m "Add qxub recipe"
git push origin add-qxub
```

Then open PR at https://github.com/conda-forge/staged-recipes

**Step 5: Review Process**

- Conda-forge bot will run tests
- Maintainers will review
- Address any feedback
- Once merged, feedstock created automatically

**Step 6: Maintain Feedstock**

After merge, you'll have access to https://github.com/conda-forge/qxub-feedstock

For future releases:
- Bot will auto-detect new GitHub releases
- Or manually update version and sha256 in feedstock

### Method 2: Custom Anaconda Channel

For internal/team distribution before public release:

**Step 1: Create Anaconda.org account**

Sign up at https://anaconda.org

**Step 2: Build package**

```bash
conda build recipe/
```

**Step 3: Upload**

```bash
# Login
anaconda login

# Upload
anaconda upload ~/miniconda3/conda-bld/noarch/qxub-3.3.1-py_0.tar.bz2

# Or specify channel/label
anaconda upload -u CHANNEL_NAME --label dev ~/miniconda3/conda-bld/noarch/qxub-3.3.1-py_0.tar.bz2
```

**Step 4: Users install with**

```bash
conda install -c YOUR_USERNAME qxub

# Or from specific label
conda install -c YOUR_USERNAME/label/dev qxub
```

### Method 3: Private/Team Channel

For internal use without public anaconda.org:

**Step 1: Set up local channel**

```bash
# Create channel directory
mkdir -p /shared/conda-channel/noarch

# Copy built package
cp ~/miniconda3/conda-bld/noarch/qxub-3.3.1-py_0.tar.bz2 \
   /shared/conda-channel/noarch/

# Index the channel
conda index /shared/conda-channel
```

**Step 2: Team members add channel**

```bash
conda config --add channels file:///shared/conda-channel
conda install qxub
```

---

## Version Management

### Version Scheme

qxub follows semantic versioning: `MAJOR.MINOR.PATCH`

- **MAJOR**: Breaking changes, major rewrites
- **MINOR**: New features, backward compatible
- **PATCH**: Bug fixes, small improvements

### Version Bump Workflow

**For each release:**

1. **Update version in both files:**
   - `qxub/__init__.py`: `__version__ = "3.3.1"`
   - `setup.py`: `version="3.3.1"`

2. **Update recipe if publishing to conda:**
   - `recipe/meta.yaml`: `{% set version = "3.3.1" %}`

3. **Commit version bump:**
   ```bash
   git add qxub/__init__.py setup.py recipe/meta.yaml
   git commit -m "Bump version to 3.3.1"
   ```

4. **Tag the release:**
   ```bash
   git tag -a v3.3.1 -m "Release v3.3.1"
   ```

### Automated Version Management (Optional)

Use `bump2version` for consistency:

```bash
pip install bump2version

# Configure .bumpversion.cfg
[bumpversion]
current_version = 3.3.1
commit = True
tag = True

[bumpversion:file:qxub/__init__.py]
search = __version__ = "{current_version}"
replace = __version__ = "{new_version}"

[bumpversion:file:setup.py]
search = version="{current_version}"
replace = version="{new_version}"

[bumpversion:file:recipe/meta.yaml]
search = {{% set version = "{current_version}" %}}
replace = {{% set version = "{new_version}" %}}

# Then bump with one command:
bump2version patch  # 3.3.1 → 3.3.2
bump2version minor  # 3.3.1 → 3.4.0
bump2version major  # 3.3.1 → 4.0.0
```

---

## Post-Release Tasks

After publishing:

- [ ] Announce release (Slack, email, etc.)
- [ ] Update documentation sites if separate
- [ ] Close milestone in GitHub (if used)
- [ ] Move issues to next milestone
- [ ] Update badges in README if needed (version, downloads, etc.)
- [ ] Monitor for installation issues
- [ ] Consider blog post for major releases

### Example Announcement

```markdown
🎉 qxub v3.3.1 released!

New in this release:
- Official conda recipe for easy installation
- MIT license added

Install:
- pip: `pip install qxub` (once on PyPI)
- conda: `conda install -c conda-forge qxub` (once on conda-forge)
- git: `pip install git+https://github.com/swarbricklab/qxub.git@v3.3.1`

Full release notes: https://github.com/swarbricklab/qxub/releases/tag/v3.3.1
```

---

## Quick Reference: Complete Release Workflow

```bash
# 1. Prepare release
git checkout -b release/v3.3.1
# ... make changes ...
git commit -m "Feature X"

# 2. Bump version (in both files)
# Edit qxub/__init__.py and setup.py
git add qxub/__init__.py setup.py
git commit -m "Bump version to 3.3.1"

# 3. Merge to main
git checkout main
git merge release/v3.3.1
git push origin main

# 4. Tag release
git tag -a v3.3.1 -m "Release v3.3.1"
git push origin v3.3.1

# 5. Create GitHub release (web UI or gh CLI)
gh release create v3.3.1 --title "v3.3.1" --notes "Release notes here"

# 6. Publish to PyPI (if ready)
python -m build
twine upload dist/*

# 7. Submit to conda-forge (if ready)
# Fork staged-recipes, add recipe, submit PR

# Done! 🎉
```

---

## Troubleshooting

### PyPI upload fails with "File already exists"

You can't re-upload the same version. Options:
- Bump to patch version (3.3.1 → 3.3.2)
- Use post-release version (3.3.1.post1)
- Delete from TestPyPI and retry there

### Conda build fails

Check:
- All dependencies available in conda-forge
- Python version compatibility
- Recipe syntax (Jinja2 templates)

### GitHub release not triggering CI

Check:
- Workflow file in `.github/workflows/`
- `on: release: types: [published]` in workflow
- Secrets configured (PYPI_API_TOKEN, etc.)

---

## Additional Resources

- **PyPI Publishing Guide**: https://packaging.python.org/tutorials/packaging-projects/
- **Conda-Forge Guide**: https://conda-forge.org/docs/maintainer/adding_pkgs.html
- **Semantic Versioning**: https://semver.org/
- **GitHub Releases**: https://docs.github.com/en/repositories/releasing-projects-on-github
