"""
Notify CLI for qxub - Send job completion notifications.

Provides the `qxub notify` command for sending notifications via Slack,
Discord, or other channels when jobs complete.

Architecture:
- Job script calls `qxub notify --job-id X --exit-code N` at completion
- This submits a tiny job to a data-mover queue (has internet access)
- The data-mover job calls `qxub notify --send ...` to actually send

This two-stage approach is needed because compute nodes typically lack
internet access on HPC systems.
"""

import json
import logging
import os
import sys
from pathlib import Path
from typing import Optional

import click

logger = logging.getLogger(__name__)


@click.group(name="notify")
def notify_cli():
    """Send job completion notifications (Slack, Discord, etc.)."""
    pass


@notify_cli.command(name="queue")
@click.option("--job-id", required=True, help="PBS job ID of the completed job")
@click.option("--exit-code", type=int, required=True, help="Exit code of the job")
@click.option("--job-name", default=None, help="Job name (for notification message)")
@click.option(
    "--output-dir",
    default=None,
    help="Directory containing job output files (for including tail in notification)",
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without submitting"
)
def queue_notification(job_id, exit_code, job_name, output_dir, dry_run):
    """Queue a notification job on a data-mover node with internet access.

    This is called from job scripts at completion. It submits a small job
    to an internet-capable queue that will send the actual notification.
    """
    from .config import manager as config_mod
    from .execution import submit_job

    config_manager = config_mod.config_manager

    # Check if notifications are configured
    slack_webhook = config_manager.get_config_value("notifications.slack.webhook_url")
    slack_user_id = config_manager.get_config_value("notifications.slack.user_id")
    discord_webhook = config_manager.get_config_value(
        "notifications.discord.webhook_url"
    )

    if not any([slack_webhook, slack_user_id, discord_webhook]):
        logger.debug("No notification channels configured, skipping")
        return

    # Build the notification command
    notify_cmd = [
        "qxub",
        "notify",
        "send",
        "--job-id",
        job_id,
        "--exit-code",
        str(exit_code),
    ]

    if job_name:
        notify_cmd.extend(["--job-name", job_name])

    if output_dir:
        notify_cmd.extend(["--output-dir", output_dir])

    # Get project from config for the notification job
    project = config_manager.get_config_value("defaults.project") or os.getenv(
        "PROJECT"
    )

    # Get the internet-capable queue (default: copyq for NCI Gadi)
    notify_queue = config_manager.get_config_value("notifications.queue") or "copyq"

    if dry_run:
        click.echo(f"Would submit notification job: {' '.join(notify_cmd)}")
        click.echo(f"  Project: {project}")
        click.echo(f"  Queue: {notify_queue}")
        return

    try:
        # Submit a tiny job to send the notification on internet-capable queue
        result = submit_job(
            command=tuple(notify_cmd),
            resources={"walltime": "00:05:00", "mem": "1GB", "ncpus": 1},
            name=f"notify-{job_id.split('.')[0]}",
            project=project,
            queue=notify_queue,
        )
        logger.info("Queued notification job %s for job %s", result.job_id, job_id)
    except Exception as e:
        # Don't fail the main job if notification fails
        logger.warning("Failed to queue notification job: %s", e)


@notify_cli.command(name="send")
@click.option("--job-id", required=True, help="PBS job ID of the completed job")
@click.option("--exit-code", type=int, required=True, help="Exit code of the job")
@click.option("--job-name", default=None, help="Job name")
@click.option(
    "--output-dir", default=None, help="Directory containing job output files"
)
def send_notification(job_id, exit_code, job_name, output_dir):
    """Send notification to configured channels (Slack, Discord, etc.).

    This is called by the data-mover job to actually send the notification.
    It runs on a node with internet access.
    """
    from .config import manager as config_mod
    from .notifications import send_all_notifications

    config_manager = config_mod.config_manager

    # Build notification context
    context = {
        "job_id": job_id,
        "exit_code": exit_code,
        "job_name": job_name or job_id.split(".")[0],
        "status": "SUCCESS" if exit_code == 0 else "FAILED",
        "output_tail": None,
    }

    # Try to get last 10 lines of output
    if output_dir:
        out_files = list(Path(output_dir).glob("*.out"))
        if out_files:
            try:
                with open(out_files[0]) as f:
                    lines = f.readlines()
                    context["output_tail"] = "".join(lines[-10:])
            except Exception as e:
                logger.debug("Could not read output file: %s", e)

    # Send to all configured channels
    try:
        send_all_notifications(context, config_manager)
        click.echo(f"✅ Notification sent for job {job_id}")
    except Exception as e:
        click.echo(f"❌ Failed to send notification: {e}", err=True)
        sys.exit(1)


@notify_cli.command(name="test")
@click.argument(
    "channel", type=click.Choice(["slack", "discord", "all"]), default="all"
)
def test_notification(channel):
    """Send a test notification to verify configuration.

    Example:
        qxub notify test slack
        qxub notify test all
    """
    from .config import manager as config_mod
    from .notifications import send_all_notifications, send_discord, send_slack

    config_manager = config_mod.config_manager

    test_context = {
        "job_id": "test-12345.gadi-pbs",
        "exit_code": 0,
        "job_name": "qxub-notify-test",
        "status": "SUCCESS",
        "output_tail": "This is a test notification from qxub.\nConfiguration is working correctly!",
    }

    channels_sent = []

    if channel in ("slack", "all"):
        webhook = config_manager.get_config_value("notifications.slack.webhook_url")
        if webhook:
            try:
                send_slack(test_context, config_manager)
                channels_sent.append("slack")
                click.echo("✅ Slack test notification sent")
            except Exception as e:
                click.echo(f"❌ Slack failed: {e}", err=True)
        else:
            click.echo("⚠️  Slack not configured (set notifications.slack.webhook_url)")

    if channel in ("discord", "all"):
        webhook = config_manager.get_config_value("notifications.discord.webhook_url")
        if webhook:
            try:
                send_discord(test_context, config_manager)
                channels_sent.append("discord")
                click.echo("✅ Discord test notification sent")
            except Exception as e:
                click.echo(f"❌ Discord failed: {e}", err=True)
        else:
            click.echo(
                "⚠️  Discord not configured (set notifications.discord.webhook_url)"
            )

    if not channels_sent:
        click.echo("\nNo notification channels are configured.")
        click.echo("Configure with:")
        click.echo(
            '  qxub config set notifications.slack.webhook_url "https://hooks.slack.com/..."'
        )
        click.echo(
            '  qxub config set notifications.discord.webhook_url "https://discord.com/api/webhooks/..."'
        )
        sys.exit(1)
