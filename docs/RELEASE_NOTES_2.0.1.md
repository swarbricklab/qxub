# qxub 2.0.1 - Documentation and Automation Improvements

## 🔧 **Infrastructure Improvements**

### 🤖 **Automated Release Workflow**
- Added GitHub Actions workflow for automated releases
- Automatic wheel building and GitHub release creation
- Smart version bump detection
- Prevents duplicate releases

### 📚 **Documentation Enhancements**
- Improved development documentation organization
- Enhanced README and configuration guides
- Better examples and usage instructions

## 🛠️ **Technical Details**

### **Release Automation**
- Triggers on version bumps in `setup.py`
- Builds Python wheels automatically
- Creates GitHub releases with release notes
- Ready for future PyPI integration

### **Documentation Structure**
- Consolidated development docs in `docs/dev/`
- Improved navigation and cross-references
- Better examples for common use cases

## 📦 **Installation**

```bash
pip install qxub==2.0.1
# or download the wheel from this release
```

## 🎯 **For Developers**

This release primarily improves the development and release process:
- Streamlined release workflow
- Better documentation for contributors
- Enhanced automation for future releases

## 🔄 **No Breaking Changes**

This is a patch release with no breaking changes. All qxub 2.0.0 functionality remains unchanged.

---

**Note**: This release tests our new automated release workflow. Future releases will benefit from this improved automation.