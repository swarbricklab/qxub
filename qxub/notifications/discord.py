"""
Discord notification channel for qxub.

Uses Discord webhook to post notifications.

Configuration:
    notifications:
      discord:
        webhook_url: "https://discord.com/api/webhooks/..."
"""

import json
import logging
import os
from typing import Any, Dict, Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)


class DiscordNotifier:
    """Send notifications to Discord."""

    def __init__(self, config_manager):
        self.config_manager = config_manager
        self.webhook_url = self._resolve_env(
            config_manager.get_config_value("notifications.discord.webhook_url")
        )

    def _resolve_env(self, value: Optional[str]) -> Optional[str]:
        """Resolve environment variable references like ${VAR}."""
        if not value:
            return value
        if value.startswith("${") and value.endswith("}"):
            env_var = value[2:-1]
            return os.environ.get(env_var)
        return value

    def send(self, context: Dict[str, Any]) -> None:
        """Send notification to Discord.

        Args:
            context: Job context with job_id, exit_code, status, etc.

        Raises:
            Exception: If sending fails.
        """
        if not self.webhook_url:
            raise ValueError("Discord webhook_url not configured")

        payload = self._build_message(context)

        data = json.dumps(payload).encode("utf-8")
        request = Request(
            self.webhook_url,
            data=data,
            headers={"Content-Type": "application/json"},
        )

        try:
            with urlopen(request, timeout=30) as response:
                # Discord returns 204 No Content on success
                if response.status not in (200, 204):
                    raise Exception(f"Discord webhook returned {response.status}")
        except HTTPError as e:
            raise Exception(f"Discord webhook failed: {e.code} {e.reason}")
        except URLError as e:
            raise Exception(f"Discord webhook connection failed: {e.reason}")

    def _build_message(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """Build Discord message payload using embeds."""
        job_id = context.get("job_id", "unknown")
        job_name = context.get("job_name", job_id.split(".")[0])
        status = context.get("status", "UNKNOWN")
        exit_code = context.get("exit_code", "?")
        output_tail = context.get("output_tail")

        # Status emoji and color
        if status == "SUCCESS":
            emoji = "✅"
            color = 0x00FF00  # Green
        else:
            emoji = "❌"
            color = 0xFF0000  # Red

        # Build embed
        embed = {
            "title": f"{emoji} Job {status}: {job_name}",
            "color": color,
            "fields": [
                {"name": "Job ID", "value": f"`{job_id}`", "inline": True},
                {"name": "Exit Code", "value": f"`{exit_code}`", "inline": True},
            ],
        }

        # Add output tail if available
        if output_tail:
            # Discord has a 1024 char limit per field
            truncated = output_tail[:1000]
            embed["fields"].append(
                {
                    "name": "Last 10 lines of output",
                    "value": f"```\n{truncated}\n```",
                    "inline": False,
                }
            )

        # Add helpful commands in footer
        embed["footer"] = {
            "text": (
                f"View output: qxub history show {job_id} --output | "
                f"View logs: qxub history show {job_id} --logs"
            )
        }

        return {
            "content": None,  # No plain text, just embed
            "embeds": [embed],
        }
