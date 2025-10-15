"""
Platform abstraction for qxub - queue definitions, resource validation, and intelligent selection.

This module provides platform-aware functionality for HPC job submission, including:
- Platform definition loading and validation
- Queue constraint checking and resource validation
- Intelligent queue selection based on resource requirements
- Service Unit cost estimation and resource optimization
- Walltime validation with complex core-count dependent rules

Platform definitions are loaded from YAML files in standard locations:
- System-level: /etc/qxub/platforms/*.yaml
- User-level: ~/.config/qxub/platforms/*.yaml
"""

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

from .resource_utils import (
    evaluate_condition,
    format_walltime,
    parse_memory_size,
    parse_walltime,
    suggest_resource_adjustment,
)

logger = logging.getLogger(__name__)


@dataclass
class QueueLimits:
    """Resource limits and constraints for a queue."""

    max_cpus: Optional[int] = None
    min_cpus: Optional[int] = None
    max_memory: Optional[str] = None
    max_local_storage: Optional[str] = None
    max_gpus: Optional[int] = None
    min_gpus: Optional[int] = None
    min_memory: Optional[str] = None  # Trigger threshold for auto-selection
    max_jobs_per_user: Optional[int] = None
    max_jobs_queued: Optional[int] = None


@dataclass
class WalltimeRule:
    """Walltime limit based on core count range."""

    cores: str  # e.g., "1-48", "96", "144,192"
    max_walltime: str  # e.g., "48:00:00"

    def matches_core_count(self, core_count: int) -> bool:
        """Check if core count matches this rule."""
        if "-" in self.cores:
            # Range like "1-48" or "720-1440"
            parts = self.cores.split("-")
            if len(parts) == 2:
                start, end = int(parts[0]), int(parts[1])
                return start <= core_count <= end
        elif "," in self.cores:
            # Discrete values like "144,192"
            values = [int(x.strip()) for x in self.cores.split(",")]
            return core_count in values
        else:
            # Single value like "96"
            return core_count == int(self.cores)
        return False


@dataclass
class ResourceValidationResult:
    """Result of resource validation against queue limits."""

    is_valid: bool = True
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)


@dataclass
class Queue:
    """Queue definition with limits, scheduling, and cost information."""

    name: str
    type: str
    limits: QueueLimits
    priority: str = "normal"  # high/normal/low
    su_billing_rate: Optional[float] = None
    walltime_rules: List[WalltimeRule] = field(default_factory=list)
    default_walltime: str = "1:00:00"
    su_rate: float = 1.0
    internet_connectivity: bool = False
    constraints: List[str] = field(default_factory=list)
    auto_min_cpus: Optional[int] = None

    def get_max_walltime(self, core_count: int) -> Optional[str]:
        """Get maximum walltime for given core count."""
        for rule in self.walltime_rules:
            if rule.matches_core_count(core_count):
                return rule.max_walltime
        return None

    def estimate_su_cost(self, cores: int, walltime_hours: float) -> float:
        """Estimate Service Unit cost for this queue."""
        return cores * walltime_hours * self.su_rate

    def validate_resources(self, resources: Dict[str, Any]) -> ResourceValidationResult:
        """Validate resource request against queue limits."""
        result = ResourceValidationResult()

        # CPU validation
        cpus = resources.get("cpus", 1)
        if self.limits.max_cpus and cpus > self.limits.max_cpus:
            result.errors.append(
                f"CPU count {cpus} exceeds queue maximum {self.limits.max_cpus}"
            )
        if self.limits.min_cpus and cpus < self.limits.min_cpus:
            result.errors.append(
                f"CPU count {cpus} below queue minimum {self.limits.min_cpus}"
            )

        # Memory validation
        memory_str = resources.get("memory")
        if memory_str and self.limits.max_memory:
            req_memory = parse_memory_size(memory_str)
            max_memory = parse_memory_size(self.limits.max_memory)
            if req_memory and max_memory and req_memory > max_memory:
                result.errors.append(
                    f"Memory {memory_str} exceeds queue maximum {self.limits.max_memory}"
                )

        if memory_str and self.limits.min_memory:
            req_memory = parse_memory_size(memory_str)
            min_memory = parse_memory_size(self.limits.min_memory)
            if req_memory and min_memory and req_memory < min_memory:
                result.warnings.append(
                    f"Memory {memory_str} below recommended minimum {self.limits.min_memory}"
                )

        # GPU validation
        gpus = resources.get("gpus", 0)
        if gpus > 0:
            if self.limits.max_gpus and gpus > self.limits.max_gpus:
                result.errors.append(
                    f"GPU count {gpus} exceeds queue maximum {self.limits.max_gpus}"
                )
            if self.limits.min_gpus and gpus < self.limits.min_gpus:
                result.errors.append(
                    f"GPU count {gpus} below queue minimum {self.limits.min_gpus}"
                )
        elif self.limits.min_gpus and self.limits.min_gpus > 0:
            result.warnings.append(
                f"This queue requires at least {self.limits.min_gpus} GPUs"
            )

        # Walltime validation
        walltime_str = resources.get("walltime")
        if walltime_str:
            max_walltime_str = self.get_max_walltime(cpus)
            if max_walltime_str:
                req_walltime = parse_walltime(walltime_str)
                max_walltime = parse_walltime(max_walltime_str)
                if req_walltime and max_walltime and req_walltime > max_walltime:
                    result.errors.append(
                        f"Walltime {walltime_str} exceeds queue maximum {max_walltime_str} for {cpus} CPUs"
                    )

        # Auto-CPU adjustment suggestion
        if self.auto_min_cpus and cpus < self.auto_min_cpus:
            result.suggestions.append(
                f"Consider increasing CPU count to {self.auto_min_cpus} for better efficiency on this queue"
            )

        result.is_valid = len(result.errors) == 0
        return result


@dataclass
class AutoSelectionRule:
    """Rule for automatic queue selection."""

    condition: str  # Boolean expression like "gpu_requested > 0"
    queue: str  # Queue name to select
    is_default: bool = False


@dataclass
class Platform:
    """Complete platform definition with queues and selection rules."""

    name: str
    type: str
    host: str
    description: str = ""
    queues: Dict[str, Queue] = field(default_factory=dict)
    auto_selection_rules: List[AutoSelectionRule] = field(default_factory=list)
    auto_adjust: Dict[str, str] = field(default_factory=dict)

    def get_queue(self, name: str) -> Optional[Queue]:
        """Get queue by name."""
        return self.queues.get(name)

    def list_queues(self) -> List[str]:
        """List all available queue names."""
        return list(self.queues.keys())

    def select_queue(self, resources: Dict[str, Any]) -> Optional[str]:
        """Select best queue for given resource requirements."""
        # Apply auto-selection rules in order
        for rule in self.auto_selection_rules:
            if rule.is_default:
                continue

            # Evaluate rule condition
            if evaluate_condition(rule.condition, resources):
                return rule.queue

        # Return default queue
        for rule in self.auto_selection_rules:
            if rule.is_default:
                return rule.queue

        # Fallback to first available queue
        queues = list(self.queues.keys())
        return queues[0] if queues else None

    def validate_queue_resources(
        self, queue_name: str, resources: Dict[str, Any]
    ) -> List[str]:
        """Validate resources against specific queue."""
        queue = self.get_queue(queue_name)
        if not queue:
            return [f"Queue '{queue_name}' not found"]
        return queue.validate_resources(resources)


class PlatformLoader:
    """Loads and manages platform definitions from YAML files."""

    def __init__(self, search_paths: Optional[List[Path]] = None):
        if search_paths:
            self.search_paths = search_paths
        else:
            # Import here to avoid circular import
            from .config_manager import config_manager

            self.search_paths = config_manager.get_platform_search_paths()

        self.platforms: Dict[str, Platform] = {}
        self._load_platforms()

    def _load_platforms(self):
        """Load all platform definitions from search paths."""
        # Check for specific platform file environment variable first
        specific_platform_file = os.getenv("QXUB_PLATFORM_FILE")
        if specific_platform_file:
            platform_path = Path(specific_platform_file)
            if platform_path.exists():
                logger.debug(
                    f"Loading platform file from QXUB_PLATFORM_FILE: {platform_path}"
                )
                try:
                    self._load_platform_file(platform_path)
                    return  # Only load the specific file, don't load others
                except Exception as e:
                    logger.error(
                        f"Failed to load platform file from QXUB_PLATFORM_FILE {platform_path}: {e}"
                    )
            else:
                logger.warning(
                    f"Platform file specified in QXUB_PLATFORM_FILE does not exist: {specific_platform_file}"
                )

        # Default behavior: load from search paths
        for search_path in self.search_paths:
            if not search_path.exists():
                logger.debug(f"Platform search path does not exist: {search_path}")
                continue

            for yaml_file in search_path.glob("*.yaml"):
                try:
                    self._load_platform_file(yaml_file)
                except Exception as e:
                    logger.error(f"Failed to load platform file {yaml_file}: {e}")

    def _load_platform_file(self, yaml_file: Path):
        """Load a single platform definition file."""
        logger.debug(f"Loading platform file: {yaml_file}")

        with open(yaml_file, "r") as f:
            data = yaml.safe_load(f)

        if not data:
            logger.warning(f"Empty platform file: {yaml_file}")
            return

        # Handle both single platform and multiple platforms in one file
        if "platform" in data:
            # Single platform format
            platform_data = data["platform"]
            platform = self._create_platform(platform_data)
            self.platforms[platform.name] = platform
            logger.info(f"Loaded platform: {platform.name}")
        elif "platforms" in data:
            # Multiple platforms format
            for platform_data in data["platforms"]:
                platform = self._create_platform(platform_data)
                self.platforms[platform.name] = platform
                logger.info(f"Loaded platform: {platform.name}")
        else:
            # Assume the whole file is a platform definition
            platform = self._create_platform(data)
            self.platforms[platform.name] = platform
            logger.info(f"Loaded platform: {platform.name}")

    def _create_platform(self, data: Dict[str, Any]) -> Platform:
        """Create Platform object from YAML data."""
        platform_name = data["name"]
        platform_type = data["type"]
        platform_host = data["host"]
        platform_description = data.get("description", "")

        # Create queues
        queues = {}
        for queue_data in data.get("queues", []):
            queue = self._create_queue(queue_data)
            queues[queue.name] = queue

        # Create auto-selection rules
        auto_selection_rules = []
        for rule_data in data.get("auto_selection_rules", []):
            rule = AutoSelectionRule(
                condition=rule_data["condition"],
                queue=rule_data["queue"],
                is_default=rule_data.get("is_default", False),
            )
            auto_selection_rules.append(rule)

        # Auto-adjust policies
        auto_adjust = data.get("auto_adjust", {})

        return Platform(
            name=platform_name,
            type=platform_type,
            host=platform_host,
            description=platform_description,
            queues=queues,
            auto_selection_rules=auto_selection_rules,
            auto_adjust=auto_adjust,
        )

    def _create_queue(self, data: Dict[str, Any]) -> Queue:
        """Create Queue object from YAML data."""
        queue_name = data["name"]
        queue_type = data["type"]

        # Create queue limits
        limits_data = data.get("limits", {})
        limits = QueueLimits(
            max_cpus=limits_data.get("max_cpus"),
            min_cpus=limits_data.get("min_cpus"),
            max_memory=limits_data.get("max_memory"),
            max_local_storage=limits_data.get("max_local_storage"),
            max_gpus=limits_data.get("max_gpus"),
            min_gpus=limits_data.get("min_gpus"),
            min_memory=limits_data.get("min_memory"),
            max_jobs_per_user=limits_data.get("max_jobs_per_user"),
            max_jobs_queued=limits_data.get("max_jobs_queued"),
        )

        # Create walltime rules
        walltime_rules = []
        for rule_data in data.get("walltime_rules", []):
            rule = WalltimeRule(
                cores=rule_data["cores"], max_walltime=rule_data["max_walltime"]
            )
            walltime_rules.append(rule)

        return Queue(
            name=queue_name,
            type=queue_type,
            limits=limits,
            priority=data.get("priority", "normal"),
            su_billing_rate=data.get("su_billing_rate"),
            walltime_rules=walltime_rules,
            default_walltime=data.get("default_walltime", "1:00:00"),
            su_rate=data.get("su_rate", 1.0),
            internet_connectivity=data.get("internet_connectivity", False),
            constraints=data.get("constraints", []),
            auto_min_cpus=data.get("auto_min_cpus"),
        )

    def get_platform(self, name: str) -> Optional[Platform]:
        """Get platform by name."""
        return self.platforms.get(name)

    def list_platforms(self) -> List[str]:
        """List all available platform names."""
        return list(self.platforms.keys())

    def reload(self):
        """Reload all platform definitions."""
        self.platforms.clear()
        self._load_platforms()

    def get_platform(self, name: str) -> Optional[Platform]:
        """Get platform by name."""
        return self.platforms.get(name)

    def list_platforms(self) -> List[str]:
        """List all available platform names."""
        return list(self.platforms.keys())

    def reload_platforms(self):
        """Reload all platform definitions."""
        self.platforms.clear()
        self._load_all_platforms()


# Global platform loader instance
_platform_loader: Optional[PlatformLoader] = None


def get_platform_loader() -> PlatformLoader:
    """Get global platform loader instance."""
    global _platform_loader
    if _platform_loader is None:
        _platform_loader = PlatformLoader()
    return _platform_loader


def get_platform(name: str) -> Optional[Platform]:
    """Get platform by name."""
    return get_platform_loader().get_platform(name)


def list_platforms() -> List[str]:
    """List all available platforms."""
    return get_platform_loader().list_platforms()


def detect_platform() -> Optional[str]:
    """Detect current platform based on hostname or environment."""
    # Check for platform override from remote execution
    platform_override = os.getenv("QXUB_PLATFORM_OVERRIDE")
    if platform_override:
        return platform_override

    hostname = os.getenv("HOSTNAME", "")

    # NCI Gadi detection
    if "gadi" in hostname.lower():
        return "nci_gadi"

    # Add other platform detection logic here

    return None


def get_current_platform() -> Optional[Platform]:
    """Get current platform based on detection or configuration."""
    platform_name = detect_platform()
    if platform_name:
        return get_platform(platform_name)

    # Fallback to first available platform
    platforms = list_platforms()
    if platforms:
        return get_platform(platforms[0])

    return None


@dataclass
class QueueSelectionResult:
    """Result of queue selection with recommendations and alternatives."""

    best_queue: Optional[str] = None
    valid_queues: List[str] = field(default_factory=list)
    invalid_queues: Dict[str, str] = field(default_factory=dict)  # queue: reason
    resource_adjustments: Optional[Dict[str, Any]] = None
    estimated_cost: Optional[float] = None
    warnings: List[str] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)


class QueueSelector:
    """Intelligent queue selection engine for resource optimization."""

    def __init__(self, platform: Platform):
        self.platform = platform

    def select_queue(
        self, resources: Dict[str, Any], preferences: Optional[Dict[str, Any]] = None
    ) -> QueueSelectionResult:
        """
        Select the best queue for given resource requirements.

        Args:
            resources: Resource requirements (cpus, memory, walltime, gpus, etc.)
            preferences: Selection preferences (priority, optimization, policy)

        Returns:
            QueueSelectionResult with selection details and recommendations
        """
        result = QueueSelectionResult()
        preferences = preferences or {}

        # First, check platform auto-selection rules
        auto_selected_queue = self.platform.select_queue(resources)
        if auto_selected_queue and auto_selected_queue in self.platform.queues:
            # Validate the auto-selected queue
            queue = self.platform.queues[auto_selected_queue]
            validation_result = queue.validate_resources(resources)

            if validation_result.is_valid:
                result.best_queue = auto_selected_queue
                result.valid_queues = [auto_selected_queue]

                # Estimate cost
                if queue.su_billing_rate:
                    cpus = resources.get("cpus", 1)
                    walltime_hours = (
                        parse_walltime(resources.get("walltime", "1:00:00")) or 1.0
                    )
                    result.estimated_cost = (
                        queue.su_billing_rate * cpus * walltime_hours
                    )

                result.suggestions.append(f"Auto-selected based on platform rules")
                return result

        # Validate and collect eligible queues
        valid_queues = []
        invalid_queues = {}

        for queue_name, queue in self.platform.queues.items():
            validation_result = queue.validate_resources(resources)

            if validation_result.is_valid:
                valid_queues.append((queue_name, queue))
            else:
                invalid_queues[queue_name] = "; ".join(validation_result.errors)

        result.valid_queues = [name for name, _ in valid_queues]
        result.invalid_queues = invalid_queues

        if not valid_queues:
            result.warnings.append("No queues satisfy the resource requirements")
            # Try to suggest adjustments for best candidate
            best_candidate = self._find_best_candidate(resources)
            if best_candidate:
                adjustments = suggest_resource_adjustment(
                    resources,
                    best_candidate[1].limits.__dict__,
                    preferences.get("adjustment_policy", "suggest"),
                )
                if adjustments:
                    result.resource_adjustments = adjustments
                    result.suggestions.append(
                        f"Consider adjusting resources for queue '{best_candidate[0]}'"
                    )
            return result

        # Score and rank valid queues
        scored_queues = []
        for queue_name, queue in valid_queues:
            score = self._score_queue(queue, resources, preferences)
            scored_queues.append((score, queue_name, queue))

        # Sort by score (highest first)
        scored_queues.sort(reverse=True)

        # Select best queue
        best_score, best_name, best_queue = scored_queues[0]
        result.best_queue = best_name

        # Estimate cost for best queue
        if best_queue.su_billing_rate:
            cpus = resources.get("cpus", 1)
            walltime_hours = parse_walltime(resources.get("walltime", "1:00:00")) or 1.0
            result.estimated_cost = best_queue.su_billing_rate * cpus * walltime_hours

        # Add suggestions for optimization
        self._add_optimization_suggestions(
            result, resources, scored_queues, preferences
        )

        return result

    def _find_best_candidate(
        self, resources: Dict[str, Any]
    ) -> Optional[Tuple[str, "Queue"]]:
        """Find the queue that's closest to satisfying requirements."""
        candidates = []

        for queue_name, queue in self.platform.queues.items():
            # Count how many constraints are violated
            violations = 0
            validation_result = queue.validate_resources(resources)
            violations = len(validation_result.errors)

            candidates.append((violations, queue_name, queue))

        if candidates:
            candidates.sort()  # Sort by fewest violations
            return (candidates[0][1], candidates[0][2])

        return None

    def _score_queue(
        self, queue: "Queue", resources: Dict[str, Any], preferences: Dict[str, Any]
    ) -> float:
        """Score a queue based on suitability for given resources and preferences."""
        score = 0.0

        optimization = preferences.get("optimization", "balanced")
        cpus = resources.get("cpus", 1)
        walltime_hours = parse_walltime(resources.get("walltime", "1:00:00")) or 1.0

        # Base score for valid queue
        score += 50.0

        # Priority bonus
        if queue.priority == "high":
            score += 30.0
        elif queue.priority == "normal":
            score += 20.0
        elif queue.priority == "low":
            score += 10.0

        # Optimization-based scoring
        if optimization == "cost":
            # Prefer queues with lower billing rate
            if queue.su_billing_rate:
                score += max(0, 20.0 - queue.su_billing_rate)
        elif optimization == "speed":
            # Prefer high-priority queues and shorter wait times
            if queue.priority == "high":
                score += 25.0
            # Prefer queues with lower utilization (if available)
        elif optimization == "balanced":
            # Balance cost and speed
            if queue.su_billing_rate:
                cost_factor = max(0, 10.0 - queue.su_billing_rate / 2)
                score += cost_factor
            if queue.priority == "high":
                score += 15.0

        # Resource efficiency scoring
        if queue.limits.max_cpus:
            cpu_efficiency = min(1.0, cpus / queue.limits.max_cpus)
            if cpu_efficiency > 0.8:
                score += 10.0  # Good utilization
            elif cpu_efficiency < 0.1:
                score -= 5.0  # Poor utilization

        # Memory efficiency
        memory_str = resources.get("memory")
        if memory_str and queue.limits.max_memory:
            req_memory = parse_memory_size(memory_str)
            max_memory = parse_memory_size(queue.limits.max_memory)
            if req_memory and max_memory:
                memory_efficiency = req_memory / max_memory
                if memory_efficiency > 0.8:
                    score += 5.0
                elif memory_efficiency < 0.1:
                    score -= 2.0

        # Walltime fit scoring
        if queue.walltime_rules:
            best_walltime_limit = None
            for rule in queue.walltime_rules:
                if rule.matches_core_count(cpus):
                    rule_limit = parse_walltime(rule.max_walltime)
                    if rule_limit and rule_limit >= walltime_hours:
                        if (
                            best_walltime_limit is None
                            or rule_limit < best_walltime_limit
                        ):
                            best_walltime_limit = rule_limit

            if best_walltime_limit:
                # Prefer tight fit to walltime limits
                walltime_efficiency = walltime_hours / best_walltime_limit
                if 0.7 <= walltime_efficiency <= 1.0:
                    score += 15.0  # Good walltime fit
                elif walltime_efficiency < 0.2:
                    score -= 5.0  # Poor walltime utilization

        # GPU efficiency
        gpus_requested = resources.get("gpus", 0)
        if gpus_requested > 0:
            if queue.limits.max_gpus and queue.limits.max_gpus >= gpus_requested:
                score += 20.0  # GPU queue available
            else:
                score -= 10.0  # GPU queue not suitable
        elif queue.limits.max_gpus and queue.limits.max_gpus > 0:
            score -= 5.0  # Penalize GPU queue for non-GPU jobs

        return score

    def _add_optimization_suggestions(
        self,
        result: QueueSelectionResult,
        resources: Dict[str, Any],
        scored_queues: List[Tuple[float, str, "Queue"]],
        preferences: Dict[str, Any],
    ):
        """Add optimization suggestions to the result."""
        optimization = preferences.get("optimization", "balanced")

        if len(scored_queues) > 1:
            best_score, best_name, best_queue = scored_queues[0]
            second_score, second_name, second_queue = scored_queues[1]

            # Suggest alternatives if they're close in score
            if abs(best_score - second_score) < 10.0:
                if (
                    optimization == "cost"
                    and second_queue.su_billing_rate
                    and best_queue.su_billing_rate
                ):
                    if second_queue.su_billing_rate < best_queue.su_billing_rate:
                        result.suggestions.append(
                            f"Consider queue '{second_name}' for lower cost "
                            f"({second_queue.su_billing_rate} vs {best_queue.su_billing_rate} SU/CPU·hour)"
                        )
                elif optimization == "speed" and second_queue.priority == "high":
                    result.suggestions.append(
                        f"Consider queue '{second_name}' for faster scheduling (high priority)"
                    )

        # Suggest resource optimizations
        cpus = resources.get("cpus", 1)
        walltime_str = resources.get("walltime", "1:00:00")

        if best_queue.limits.max_cpus and cpus < best_queue.limits.max_cpus / 2:
            result.suggestions.append(
                f"Consider increasing CPU count to utilize queue '{best_name}' more efficiently "
                f"(current: {cpus}, max: {best_queue.limits.max_cpus})"
            )

        # Walltime optimization suggestions
        walltime_hours = parse_walltime(walltime_str)
        if walltime_hours and best_queue.walltime_rules:
            for rule in best_queue.walltime_rules:
                if rule.matches_core_count(cpus):
                    max_walltime_hours = parse_walltime(rule.max_walltime)
                    if max_walltime_hours and walltime_hours < max_walltime_hours / 2:
                        result.suggestions.append(
                            f"You could request up to {rule.max_walltime} walltime "
                            f"for {cpus} CPUs on queue '{best_name}'"
                        )
                    break


def select_best_queue(
    resources: Dict[str, Any],
    platform_name: Optional[str] = None,
    preferences: Optional[Dict[str, Any]] = None,
) -> QueueSelectionResult:
    """
    High-level function to select the best queue for given resources.

    Args:
        resources: Resource requirements
        platform_name: Platform name (auto-detected if None)
        preferences: Selection preferences (merged with config)

    Returns:
        QueueSelectionResult with selection details
    """
    # Import here to avoid circular import
    from .config_manager import config_manager

    # Get platform
    if platform_name:
        platform = get_platform(platform_name)
    else:
        platform_name = config_manager.get_default_platform()
        if platform_name:
            platform = get_platform(platform_name)
        else:
            platform = get_current_platform()

    if not platform:
        result = QueueSelectionResult()
        result.warnings.append("No platform available or detected")
        return result

    # Merge preferences with config
    config_prefs = config_manager.get_queue_preferences()
    final_preferences = {**config_prefs}
    if preferences:
        final_preferences.update(preferences)

    # Select queue
    selector = QueueSelector(platform)
    return selector.select_queue(resources, final_preferences)
