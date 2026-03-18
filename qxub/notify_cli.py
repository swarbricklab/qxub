"""
Notify CLI for qxub - Send job completion notifications.

Provides the `qxub notify` command for sending notifications via Slack,
Discord, or other channels when jobs complete.

Architecture:
- Job script calls: qxub exec --internet --terse --default -- qxub notify send ...
- The notification runs on an internet-capable queue (e.g., copyq)
- This approach reuses qxub exec's queue selection instead of duplicating logic
"""

import logging
import sys
from pathlib import Path

import click

logger = logging.getLogger(__name__)


@click.group(name="notify")
def notify_cli():
    """Send job completion notifications (Slack, Discord, etc.)."""
    pass


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
