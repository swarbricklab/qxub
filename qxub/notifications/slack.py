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
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

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
            raise ValueError(
                "Slack not configured: set webhook_url or (bot_token + user_id)"
            )

    def _send_webhook(self, context: Dict[str, Any]) -> None:
        """Send via incoming webhook."""
        payload = self._build_message(context)

        data = json.dumps(payload).encode("utf-8")
        request = Request(
            self.webhook_url,
            data=data,
            headers={"Content-Type": "application/json"},
        )

        try:
            with urlopen(request, timeout=30) as response:
                if response.status != 200:
                    raise Exception(f"Slack webhook returned {response.status}")
        except HTTPError as e:
            raise Exception(f"Slack webhook failed: {e.code} {e.reason}")
        except URLError as e:
            raise Exception(f"Slack webhook connection failed: {e.reason}")

    def _send_dm(self, context: Dict[str, Any]) -> None:
        """Send DM via Bot API."""
        # First, open a conversation with the user
        conv_url = "https://slack.com/api/conversations.open"
        conv_data = json.dumps({"users": self.user_id}).encode("utf-8")
        conv_request = Request(
            conv_url,
            data=conv_data,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {self.bot_token}",
            },
        )

        try:
            with urlopen(conv_request, timeout=30) as response:
                result = json.loads(response.read().decode("utf-8"))
                if not result.get("ok"):
                    raise Exception(f"Failed to open DM: {result.get('error')}")
                channel_id = result["channel"]["id"]
        except HTTPError as e:
            raise Exception(f"Slack API failed: {e.code} {e.reason}")
        except URLError as e:
            raise Exception(f"Slack API connection failed: {e.reason}")

        # Send the message
        payload = self._build_message(context)
        payload["channel"] = channel_id

        msg_url = "https://slack.com/api/chat.postMessage"
        msg_data = json.dumps(payload).encode("utf-8")
        msg_request = Request(
            msg_url,
            data=msg_data,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {self.bot_token}",
            },
        )

        try:
            with urlopen(msg_request, timeout=30) as response:
                result = json.loads(response.read().decode("utf-8"))
                if not result.get("ok"):
                    raise Exception(f"Failed to send message: {result.get('error')}")
        except HTTPError as e:
            raise Exception(f"Slack API failed: {e.code} {e.reason}")
        except URLError as e:
            raise Exception(f"Slack API connection failed: {e.reason}")

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
                            f"View output: `qxub history show {job_id} --output`\n"
                            f"View logs: `qxub history show {job_id} --logs`"
                        ),
                    }
                ],
            }
        )

        return {
            "text": f"{emoji} Job {status}: {job_name} ({job_id})",  # Fallback
            "blocks": blocks,
        }
