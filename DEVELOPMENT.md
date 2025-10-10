# qxub v2.2 Development Branch

## 🚀 Welcome to v2.2 Development!

This branch is for developing **Phase 2.2: Remote Execution** features as outlined in the [v2 roadmap](docs/dev/v2_roadmap.md).

### 🎯 Phase 2.2 Goals: Remote Execution

The goal of v2.2 is to enable executing jobs on remote platforms from local machines, transforming qxub into a distributed job submission system.

### 🌐 Core Features Planned

- **✅ Client-server architecture** for remote job submission
- **✅ SSH-based remote execution** with credential delegation
- **✅ Platform profiles** combining host + platform + credentials
- **✅ Distributed logging** (client + server)
- **✅ Exit code propagation** through the execution chain

### 📋 Key Tasks for v2.2

1. **Platform profiles** with remote host configuration
2. **SSH execution backend** for secure remote job submission
3. **Remote qxub invocation** and monitoring
4. **Output streaming** from remote to local terminal
5. **Distributed history logging** across client and server

### 💡 Usage Vision

```bash
# Execute jobs on NCI from local laptop
qxub --profile gadi --env dvc3 -- dvc push

# Auto-select queue on remote platform
qxub --profile gadi --queue auto -l mem=500GB --env pytorch -- python train.py

# Stream output from remote job to local terminal in real-time
qxub --profile cluster --conda myenv -- python long_running_job.py
```

### 🛠️ Technical Architecture

Building on the solid foundation of v2.1's platform abstraction:

- **Platform profiles**: Extend platform definitions with remote host configuration
- **SSH backend**: Secure execution channel with credential management
- **Distributed monitoring**: Real-time output streaming and status updates
- **History integration**: Unified logging across local and remote executions

### 📈 Success Criteria

- Execute jobs on remote HPC from local machine with same qxub syntax
- Real-time output streaming to local terminal
- Proper exit code handling and error propagation
- Unified execution history (local + remote)
- Seamless integration with existing v2.1 platform system

---

**Ready to build the future of distributed HPC job submission!** 🚀

For detailed implementation plans, see [docs/dev/v2_roadmap.md](docs/dev/v2_roadmap.md).
