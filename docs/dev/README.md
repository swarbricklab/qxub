# qxub Developer Documentation

Welcome to the qxub developer documentation! This section provides detailed technical information for developers working on or extending qxub.

## Threading System Documentation

qxub uses a sophisticated multi-threaded architecture to provide real-time job monitoring with clean output streaming. Understanding this system is crucial for debugging issues and making modifications.

### [📚 Threading Architecture](threading-architecture.md)
**Complete technical guide** to qxub's threading system:
- **OutputCoordinator**: Central thread synchronization hub
- **Thread Responsibilities**: Monitor, STDOUT/STDERR tail, spinner
- **Signal Flow**: Event-based coordination mechanisms
- **Exit Code Propagation**: How job exit codes flow through threads
- **Graceful Shutdown**: Ctrl-C handling and resource cleanup
- **Performance**: Resource usage and scalability characteristics

### [📊 Threading Diagrams](threading-diagrams.md)
**Visual representations** of thread interactions:
- **State Diagrams**: Thread lifecycle and state transitions
- **Sequence Diagrams**: Message flow between components
- **Event Timeline**: Timing of signals and coordination
- **Control Flow**: Different execution scenarios (success, failure, interruption)
- **Memory Layout**: Resource usage patterns

### [🔧 Threading Troubleshooting](threading-troubleshooting.md)
**Practical debugging guide** for threading issues:
- **Common Problems**: Hanging processes, wrong exit codes, missing output
- **Diagnostic Tools**: Debug logging, thread inspection, performance analysis
- **Testing Strategies**: Unit testing, integration testing, stress testing
- **Emergency Procedures**: Kill hung processes, clean up orphaned jobs
- **Prevention**: Best practices for thread-safe development

## System Design Documentation

### [⚙️ Config and Alias System Design](config-and-alias-system-design.md)
**Architecture documentation** for configuration and alias systems:
- **Configuration Hierarchy**: File precedence and inheritance
- **Alias System**: Hierarchical structure and execution flow
- **Template Variables**: Dynamic value substitution
- **Validation**: Option validation and error handling

## Quick Reference

### Essential Debugging Commands
```bash
# Enable debug logging for threading issues
export QXUB_LOG_LEVEL=DEBUG
qxub conda --env myenv script.py

# Check thread status
ps aux | grep qxub

# View active threads in Python
python -c "import threading; print([t.name for t in threading.enumerate()])"

# Kill hung qxub processes
pkill -TERM qxub
```

### Key Files for Threading
- `qxub/scheduler.py` - Core threading logic, OutputCoordinator
- `qxub/conda.py` - Conda executor with threading integration
- `qxub/module.py` - Module executor with threading integration
- `qxub/sing.py` - Singularity executor with threading integration

### Testing Threading Code
```python
# Always use timeouts in tests
thread.join(timeout=5)
assert not thread.is_alive()

# Test shutdown scenarios
coordinator.signal_shutdown()
assert coordinator.should_shutdown()

# Mock the coordinator for unit tests
from unittest.mock import Mock
coordinator = Mock()
```

## Contributing to Threading System

When modifying the threading system:

1. **Read the architecture docs first** - Understand the event-based coordination
2. **Test shutdown scenarios** - Ensure clean Ctrl-C handling
3. **Use timeouts everywhere** - Prevent hanging threads
4. **Follow the event patterns** - Use Events for coordination, not shared variables
5. **Add debug logging** - Include thread names and timing information
6. **Test with real PBS jobs** - Integration testing is crucial

## Getting Help

If you encounter threading-related issues:

1. **Check the troubleshooting guide** - Most common issues are documented
2. **Enable debug logging** - Use `QXUB_LOG_LEVEL=DEBUG`
3. **Gather thread dumps** - Include in bug reports
4. **Test in isolation** - Reproduce with minimal examples
5. **Document edge cases** - Help improve these docs

For questions or contributions, see the main project README for contact information.
