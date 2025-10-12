"""
Platform management CLI commands for qxub.

Provides commands for platform discovery, queue selection, resource validation,
and cost estimation using the qxub platform abstraction system.
"""

import sys
from pathlib import Path
from typing import Any, Dict, Optional

import click

from .config import get_config, setup_logging
from .platform import (
    Platform,
    QueueSelectionResult,
    get_current_platform,
    get_platform,
    list_platforms,
    select_best_queue,
)
from .resource_utils import format_walltime, parse_memory_size, parse_walltime


@click.group(name="platform")
def platform_cli():
    """Platform and queue management commands."""
    pass


@platform_cli.command("list")
@click.option("-v", "--verbose", count=True, help="Show detailed platform information")
def list_platforms_cmd(verbose):
    """List available platforms."""
    platforms = list_platforms()

    if not platforms:
        click.echo("No platforms found.")
        click.echo("Check platform search paths with: qxub config show")
        return

    if verbose == 0:
        # Simple list
        for platform_name in sorted(platforms):
            click.echo(platform_name)
    else:
        # Detailed information
        for platform_name in sorted(platforms):
            platform = get_platform(platform_name)
            if platform:
                click.echo(f"\n{platform_name}:")
                click.echo(f"  Type: {platform.type}")
                click.echo(f"  Host: {platform.host}")
                if platform.description:
                    click.echo(f"  Description: {platform.description}")
                click.echo(f"  Queues: {len(platform.queues)} available")

                if verbose >= 2:
                    for queue_name, queue in platform.queues.items():
                        click.echo(f"    {queue_name}: {queue.type} queue")


@platform_cli.command("info")
@click.argument("platform_name", required=False)
def platform_info(platform_name):
    """Show detailed information about a platform."""
    if not platform_name:
        # Show current platform
        platform = get_current_platform()
        if not platform:
            click.echo("No platform detected or configured.")
            click.echo("Use 'qxub platform list' to see available platforms.")
            return
        platform_name = platform.name
    else:
        platform = get_platform(platform_name)
        if not platform:
            click.echo(f"Platform '{platform_name}' not found.")
            click.echo("Use 'qxub platform list' to see available platforms.")
            return

    click.echo(f"Platform: {platform.name}")
    click.echo(f"Type: {platform.type}")
    click.echo(f"Host: {platform.host}")
    if platform.description:
        click.echo(f"Description: {platform.description}")

    click.echo(f"\nQueues ({len(platform.queues)}):")
    for queue_name, queue in platform.queues.items():
        click.echo(f"\n  {queue_name} ({queue.type}):")
        click.echo(f"    Priority: {queue.priority}")

        # Show limits
        limits = queue.limits
        if limits.max_cpus:
            cpu_range = f"{limits.min_cpus or 1}-{limits.max_cpus}"
            click.echo(f"    CPUs: {cpu_range}")

        if limits.max_memory:
            mem_info = f"up to {limits.max_memory}"
            if limits.min_memory:
                mem_info = f"{limits.min_memory} - {limits.max_memory}"
            click.echo(f"    Memory: {mem_info}")

        if limits.max_gpus:
            gpu_range = f"{limits.min_gpus or 0}-{limits.max_gpus}"
            click.echo(f"    GPUs: {gpu_range}")

        # Show walltime rules
        if queue.walltime_rules:
            click.echo(f"    Walltime limits:")
            for rule in queue.walltime_rules:
                click.echo(f"      {rule.cores} cores: {rule.max_walltime}")

        # Show billing
        if queue.su_billing_rate:
            click.echo(f"    Billing: {queue.su_billing_rate} SU/CPU·hour")

    # Show auto-selection rules
    if platform.auto_selection_rules:
        click.echo(f"\nAuto-selection rules:")
        for rule in platform.auto_selection_rules:
            if rule.is_default:
                click.echo(f"  Default: {rule.queue}")
            else:
                click.echo(f"  If {rule.condition}: {rule.queue}")


@click.command("select-queue")
@click.option("--cpus", type=int, help="Number of CPUs required")
@click.option("--memory", help="Memory required (e.g., 8GB, 512MB)")
@click.option("--walltime", help="Walltime required (e.g., 2:00:00, 1h)")
@click.option("--gpus", type=int, help="Number of GPUs required")
@click.option("--platform", help="Platform name (auto-detected if not specified)")
@click.option(
    "--optimization",
    type=click.Choice(["cost", "speed", "balanced"]),
    help="Optimization preference",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["human", "json", "queue-name"]),
    default="human",
    help="Output format",
)
def select_queue_cmd(
    cpus, memory, walltime, gpus, platform, optimization, output_format
):
    """Select the best queue for given resource requirements."""
    # Build resource requirements
    resources = {}
    if cpus:
        resources["cpus"] = cpus
    if memory:
        resources["memory"] = memory
    if walltime:
        resources["walltime"] = walltime
    if gpus:
        resources["gpus"] = gpus

    # Set defaults if nothing specified
    if not resources:
        resources = {"cpus": 1, "memory": "4GB", "walltime": "1:00:00"}
        click.echo(
            "No resources specified, using defaults: 1 CPU, 4GB memory, 1 hour",
            err=True,
        )

    # Build preferences
    preferences = {}
    if optimization:
        preferences["optimization"] = optimization

    # Select queue
    result = select_best_queue(resources, platform, preferences)

    if output_format == "queue-name":
        if result.best_queue:
            click.echo(result.best_queue)
        else:
            sys.exit(1)
    elif output_format == "json":
        import json

        output = {
            "best_queue": result.best_queue,
            "valid_queues": result.valid_queues,
            "invalid_queues": result.invalid_queues,
            "estimated_cost": result.estimated_cost,
            "warnings": result.warnings,
            "suggestions": result.suggestions,
            "resource_adjustments": result.resource_adjustments,
        }
        click.echo(json.dumps(output, indent=2))
    else:  # human format
        if result.best_queue:
            click.echo(f"✅ Recommended queue: {result.best_queue}")

            if result.estimated_cost:
                click.echo(f"💰 Estimated cost: {result.estimated_cost:.2f} SU")

            if result.valid_queues and len(result.valid_queues) > 1:
                other_queues = [
                    q for q in result.valid_queues if q != result.best_queue
                ]
                click.echo(f"📋 Other valid queues: {', '.join(other_queues)}")
        else:
            click.echo("❌ No suitable queue found")

        # Show warnings
        for warning in result.warnings:
            click.echo(f"⚠️  {warning}")

        # Show suggestions
        for suggestion in result.suggestions:
            click.echo(f"💡 {suggestion}")

        # Show resource adjustments
        if result.resource_adjustments:
            click.echo("🔧 Suggested resource adjustments:")
            for key, value in result.resource_adjustments.items():
                if not key.endswith("_reason"):
                    reason_key = f"{key}_reason"
                    reason = result.resource_adjustments.get(reason_key, "")
                    if reason:
                        click.echo(f"   {key}: {value} ({reason})")
                    else:
                        click.echo(f"   {key}: {value}")

        # Show invalid queues
        if result.invalid_queues:
            click.echo("❌ Invalid queues:")
            for queue_name, reason in result.invalid_queues.items():
                click.echo(f"   {queue_name}: {reason}")

        if not result.best_queue:
            sys.exit(1)


@click.command("validate")
@click.option("--cpus", type=int, help="Number of CPUs to validate")
@click.option("--memory", help="Memory to validate (e.g., 8GB, 512MB)")
@click.option("--walltime", help="Walltime to validate (e.g., 2:00:00, 1h)")
@click.option("--gpus", type=int, help="Number of GPUs to validate")
@click.option("--queue", help="Specific queue to validate against")
@click.option("--platform", help="Platform name (auto-detected if not specified)")
def validate_cmd(cpus, memory, walltime, gpus, queue, platform):
    """Validate resource requirements against platform constraints."""
    # Build resource requirements
    resources = {}
    if cpus:
        resources["cpus"] = cpus
    if memory:
        resources["memory"] = memory
    if walltime:
        resources["walltime"] = walltime
    if gpus:
        resources["gpus"] = gpus

    if not resources:
        click.echo("No resources specified to validate.")
        return

    # Get platform
    if platform:
        plat = get_platform(platform)
        if not plat:
            click.echo(f"Platform '{platform}' not found.")
            sys.exit(1)
    else:
        plat = get_current_platform()
        if not plat:
            click.echo("No platform detected. Specify --platform.")
            sys.exit(1)

    click.echo(f"Validating resources against platform: {plat.name}")
    click.echo(f"Resources: {resources}")

    if queue:
        # Validate against specific queue
        if queue not in plat.queues:
            click.echo(f"❌ Queue '{queue}' not found on platform '{plat.name}'")
            sys.exit(1)

        queue_obj = plat.queues[queue]
        result = queue_obj.validate_resources(resources)

        click.echo(f"\nQueue: {queue}")
        if result.is_valid:
            click.echo("✅ Resources are valid for this queue")
        else:
            click.echo("❌ Resources are invalid for this queue")
            for error in result.errors:
                click.echo(f"   Error: {error}")

        for warning in result.warnings:
            click.echo(f"⚠️  Warning: {warning}")

        for suggestion in result.suggestions:
            click.echo(f"💡 Suggestion: {suggestion}")
    else:
        # Validate against all queues
        click.echo("\nValidation results by queue:")
        valid_queues = []
        invalid_queues = []

        for queue_name, queue_obj in plat.queues.items():
            result = queue_obj.validate_resources(resources)
            if result.is_valid:
                click.echo(f"✅ {queue_name}: Valid")
                valid_queues.append(queue_name)
            else:
                click.echo(f"❌ {queue_name}: Invalid")
                for error in result.errors:
                    click.echo(f"   {error}")
                invalid_queues.append(queue_name)

        click.echo(f"\nSummary:")
        click.echo(
            f"Valid queues: {', '.join(valid_queues) if valid_queues else 'None'}"
        )
        click.echo(
            f"Invalid queues: {', '.join(invalid_queues) if invalid_queues else 'None'}"
        )


@click.command("estimate")
@click.option("--cpus", type=int, help="Number of CPUs")
@click.option("--walltime", help="Walltime (e.g., 2:00:00, 1h)")
@click.option("--queue", help="Specific queue to estimate cost for")
@click.option("--platform", help="Platform name (auto-detected if not specified)")
def estimate_cmd(cpus, walltime, queue, platform):
    """Estimate Service Unit cost for resource requirements."""
    cpus = cpus or 1
    walltime_str = walltime or "1:00:00"

    # Parse walltime
    walltime_hours = parse_walltime(walltime_str)
    if not walltime_hours:
        click.echo(f"Invalid walltime format: {walltime_str}")
        sys.exit(1)

    # Get platform
    if platform:
        plat = get_platform(platform)
        if not plat:
            click.echo(f"Platform '{platform}' not found.")
            sys.exit(1)
    else:
        plat = get_current_platform()
        if not plat:
            click.echo("No platform detected. Specify --platform.")
            sys.exit(1)

    click.echo(f"Cost estimation for platform: {plat.name}")
    click.echo(f"Resources: {cpus} CPUs × {format_walltime(walltime_hours)}")

    if queue:
        # Estimate for specific queue
        if queue not in plat.queues:
            click.echo(f"❌ Queue '{queue}' not found")
            sys.exit(1)

        queue_obj = plat.queues[queue]
        if queue_obj.su_billing_rate:
            cost = cpus * walltime_hours * queue_obj.su_billing_rate
            click.echo(f"\nQueue: {queue}")
            click.echo(f"Rate: {queue_obj.su_billing_rate} SU/CPU·hour")
            click.echo(f"💰 Estimated cost: {cost:.2f} SU")
        else:
            click.echo(f"❌ No billing rate available for queue '{queue}'")
    else:
        # Estimate for all queues
        click.echo("\nCost estimates by queue:")

        estimates = []
        for queue_name, queue_obj in plat.queues.items():
            if queue_obj.su_billing_rate:
                cost = cpus * walltime_hours * queue_obj.su_billing_rate
                estimates.append((cost, queue_name, queue_obj.su_billing_rate))
                click.echo(
                    f"  {queue_name}: {cost:.2f} SU (rate: {queue_obj.su_billing_rate})"
                )
            else:
                click.echo(f"  {queue_name}: No billing rate available")

        if estimates:
            estimates.sort()  # Sort by cost
            cheapest_cost, cheapest_queue, cheapest_rate = estimates[0]
            click.echo(
                f"\n💰 Cheapest option: {cheapest_queue} ({cheapest_cost:.2f} SU)"
            )


# Export for registration
__all__ = ["platform_cli", "select_queue_cmd", "validate_cmd", "estimate_cmd"]
