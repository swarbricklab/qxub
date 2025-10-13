# qxub v2.2 Development Branch

## 🚀 Welcome to v2.2 Development!

This branch is for developing **Phase 2.2: Remote Execution** features as outlined in the [v2 roadmap](docs/dev/v2_roadmap.md).

### 🎯 Phase 2.2 Goals: Remote Execution

The goal of v2.2 is to enable executing jobs on remote platforms from local machines, transforming qxub into a distributed job submission system.

### 🌐 Core Features Planned

- **✅ Client-server architecture** for remote job submission
- **✅ SSH-based remote execution** with credential delegation
- **✅ Remote configurations** combining connection details and preferences
- **✅ Conda environment setup** for simplified remote qxub activation
- **✅ Real-time output streaming** from remote to local terminal

### 📋 Key Tasks for v2.2

1. **Remote user configuration** in `~/.config/qxub/config.yaml`
2. **SSH execution backend** for secure remote job submission
3. **Conda environment activation** for remote qxub setup
4. **Output streaming** from remote to local terminal
5. **Smart working directory resolution** based on project structure

### 💡 Usage Vision

```bash
# Execute jobs on NCI from local laptop
qxub --remote nci_gadi --env dvc3 -- dvc push

# Auto-select queue on remote platform
qxub --remote nci_gadi --queue auto -l mem=500GB --env pytorch -- python train.py

# Stream output from remote job to local terminal in real-time
qxub --remote cluster --env myenv -- python long_running_job.py
```

### 🛠️ Technical Architecture

Building on the solid foundation of v2.1's platform abstraction:

- **Remote configurations**: User-specific settings for connecting to remote systems
- **SSH backend**: Secure execution channel with credential delegation to SSH
- **Conda environments**: Simple remote setup via environment activation
- **Smart working directories**: Project-based remote path resolution
- **Real-time streaming**: Output forwarding from remote to local terminal

### 📈 Success Criteria

- Execute jobs on remote HPC from local machine with same qxub syntax
- Real-time output streaming to local terminal
- Proper exit code handling and error propagation
- Clean separation of platform capabilities and user connection settings
- Seamless integration with existing v2.1 platform system

---

**Ready to build the future of distributed HPC job submission!** 🚀

For detailed implementation plans, see [docs/dev/v2_roadmap.md](docs/dev/v2_roadmap.md).
