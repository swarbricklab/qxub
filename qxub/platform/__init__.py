"""
Platform abstraction package for qxub.

This package provides platform-aware functionality for HPC job submission:
- Platform definition loading and validation
- Queue constraint checking and resource validation
- Intelligent queue selection based on resource requirements
- Service Unit cost estimation and resource optimization
- CLI commands for platform management

Main classes and functions:
- Platform: Core platform abstraction
- Queue: Queue definition with limits and rules
- QueueLimits: Resource constraints for queues
- QueueSelectionResult: Result from queue selection
- PlatformLoader: Platform discovery and loading
- get_platform: Get specific platform by name
- list_platforms: List all available platforms
- detect_platform: Auto-detect current platform
- get_current_platform: Get platform for current environment
- select_best_queue: Intelligent queue selection
"""

# CLI commands for platform management
from .cli import estimate_cmd, platform_cli, select_queue_cmd, validate_cmd

# Core platform classes
from .core import (
    AutoSelectionRule,
    Platform,
    PlatformLoader,
    Queue,
    QueueLimits,
    QueueSelectionResult,
    QueueSelector,
    ResourceValidationResult,
    WalltimeRule,
    detect_platform,
    get_current_platform,
    get_platform,
    get_platform_loader,
    list_platforms,
    select_best_queue,
)

# Platform integration utilities
# from .integration import load_platform_with_remote

__all__ = [
    # Core classes
    "AutoSelectionRule",
    "Platform",
    "PlatformLoader",
    "Queue",
    "QueueLimits",
    "QueueSelectionResult",
    "QueueSelector",
    "ResourceValidationResult",
    "WalltimeRule",
    # Core functions
    "detect_platform",
    "get_current_platform",
    "get_platform",
    "get_platform_loader",
    "list_platforms",
    "select_best_queue",
    # CLI commands
    "estimate_cmd",
    "platform_cli",
    "select_queue_cmd",
    "validate_cmd",
    # Integration
    # "load_platform_with_remote",
]
