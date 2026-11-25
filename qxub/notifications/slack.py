"""
Slack notification channel for qxub.

Supports two modes:
1. Incoming Webhook: Posts to a channel via webhook URL
2. Bot Token + User ID: Sends DM to a specific user

Configuration:
    notifications:
      slack:
        webhook_url: "https://hooks.slack.com/services/..."  # Option 1
        # OR
        bot_token: "${QXUB_SLACK_BOT_TOKEN}"  # Option 2
        user_id: "U12345ABCDE"
"""

import json
import logging
import os
from typing import Any, Dict, Optional

from .webhook import NotificationError, check_webhook_url, post_json

logger = logging.getLogger(__name__)


class SlackNotifier:
    """Send notifications to Slack."""

    def __init__(self, config_manager):
        self.config_manager = config_manager
        self.webhook_url = self._resolve_env(
            config_manager.get_config_value("notifications.slack.webhook_url")
        )
        self.bot_token = self._resolve_env(
            config_manager.get_config_value("notifications.slack.bot_token")
        )
        self.user_id = config_manager.get_config_value("notifications.slack.user_id")

    def _resolve_env(self, value: Optional[str]) -> Optional[str]:
        """Resolve environment variable references like ${VAR}."""
        if not value:
            return value
        if value.startswith("${") and value.endswith("}"):
            env_var = value[2:-1]
            return os.environ.get(env_var)
        return value

    def send(self, context: Dict[str, Any]) -> None:
        """Send notification to Slack.

        Args:
            context: Job context with job_id, exit_code, status, etc.

        Raises:
            Exception: If sending fails.
        """
        if self.webhook_url:
            self._send_webhook(context)
        elif self.bot_token and self.user_id:
            self._send_dm(context)
        else:
            raise NotificationError(
                "Slack not configured: set webhook_url or (bot_token + user_id)"
            )

    def _send_webhook(self, context: Dict[str, Any]) -> None:
        """Send via incoming webhook."""
        check_webhook_url(self.webhook_url, "notifications.slack.webhook_url")
        post_json(self.webhook_url, self._build_message(context), "Slack webhook")

    def _send_dm(self, context: Dict[str, Any]) -> None:
        """Send DM via Bot API."""
        headers = {"Authorization": f"Bearer {self.bot_token}"}

        # First, open a conversation with the user
        result = self._api_call(
            "https://slack.com/api/conversations.open",
            {"users": self.user_id},
            headers,
        )
        if not result.get("ok"):
            raise NotificationError(f"Failed to open DM: {result.get('error')}")
        channel_id = result["channel"]["id"]

        # Send the message
        payload = self._build_message(context)
        payload["channel"] = channel_id
        result = self._api_call(
            "https://slack.com/api/chat.postMessage", payload, headers
        )
        if not result.get("ok"):
            raise NotificationError(f"Failed to send message: {result.get('error')}")

    @staticmethod
    def _api_call(url: str, payload: Dict[str, Any], headers) -> Dict[str, Any]:
        """POST to the Slack API; errors never quote the bot token."""
        body = post_json(url, payload, "Slack API", headers=headers)
        try:
            return json.loads(body.decode("utf-8"))
        except ValueError:
            raise NotificationError("Slack API returned an unreadable reply") from None

    def _build_message(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """Build Slack message payload using Block Kit."""
        job_id = context.get("job_id", "unknown")
        job_name = context.get("job_name", job_id.split(".")[0])
        status = context.get("status", "UNKNOWN")
        exit_code = context.get("exit_code", "?")
        output_tail = context.get("output_tail")

        # Status emoji and color
        if status == "SUCCESS":
            emoji = "✅"
            color = "good"
        else:
            emoji = "❌"
            color = "danger"

        # Build blocks
        blocks = [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": f"{emoji} Job {status}: {job_name}",
                    "emoji": True,
                },
            },
            {
                "type": "section",
                "fields": [
                    {"type": "mrkdwn", "text": f"*Job ID:*\n`{job_id}`"},
                    {"type": "mrkdwn", "text": f"*Exit Code:*\n`{exit_code}`"},
                ],
            },
        ]

        # Add output tail if available
        if output_tail:
            blocks.append({"type": "divider"})
            blocks.append(
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": f"*Last 10 lines of output:*\n```{output_tail[:2000]}```",
                    },
                }
            )

        # Add helpful commands
        blocks.append({"type": "divider"})
        blocks.append(
            {
                "type": "context",
                "elements": [
                    {
                        "type": "mrkdwn",
                        "text": (
                            f"View output: `qxub history out {job_id}`\n"
                            f"View logs: `qxub history log {job_id}`"
                        ),
                    }
                ],
            }
        )

        return {
            "text": f"{emoji} Job {status}: {job_name} ({job_id})",  # Fallback
            "blocks": blocks,
        }
