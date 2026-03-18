"""
Notification channels for qxub job completion.

This module provides functions to send notifications via various channels
(Slack, Discord, etc.) when PBS jobs complete.
"""

import logging
from typing import Any, Dict, List

logger = logging.getLogger(__name__)


def send_all_notifications(context: Dict[str, Any], config_manager) -> List[str]:
    """Send notifications to all configured channels.

    Args:
        context: Dictionary with job information:
            - job_id: PBS job ID
            - exit_code: Job exit code
            - job_name: Job name
            - status: "SUCCESS" or "FAILED"
            - output_tail: Last N lines of output (optional)
        config_manager: qxub ConfigManager instance

    Returns:
        List of channel names that were successfully notified.

    Raises:
        Exception: If all configured channels fail.
    """
    channels_sent = []
    errors = []

    # Check Slack
    slack_webhook = config_manager.get_config_value("notifications.slack.webhook_url")
    if slack_webhook:
        try:
            send_slack(context, config_manager)
            channels_sent.append("slack")
        except Exception as e:
            logger.error("Slack notification failed: %s", e)
            errors.append(f"Slack: {e}")

    # Check Discord
    discord_webhook = config_manager.get_config_value(
        "notifications.discord.webhook_url"
    )
    if discord_webhook:
        try:
            send_discord(context, config_manager)
            channels_sent.append("discord")
        except Exception as e:
            logger.error("Discord notification failed: %s", e)
            errors.append(f"Discord: {e}")

    if not channels_sent and errors:
        raise Exception(f"All notification channels failed: {'; '.join(errors)}")

    return channels_sent


def send_slack(context: Dict[str, Any], config_manager) -> None:
    """Send notification to Slack.

    Args:
        context: Job context dictionary.
        config_manager: qxub ConfigManager instance.

    Raises:
        Exception: If sending fails.
    """
    from .slack import SlackNotifier

    notifier = SlackNotifier(config_manager)
    notifier.send(context)


def send_discord(context: Dict[str, Any], config_manager) -> None:
    """Send notification to Discord.

    Args:
        context: Job context dictionary.
        config_manager: qxub ConfigManager instance.

    Raises:
        Exception: If sending fails.
    """
    from .discord import DiscordNotifier

    notifier = DiscordNotifier(config_manager)
    notifier.send(context)
