"""Posting to Slack and Discord without leaking their secrets.

A webhook URL is itself the credential, and a bot token is sent in a
header. Errors from urllib or http.client can quote either, and
``qxub notify send`` runs in a job whose stderr goes to a file other users
may read. So every error raised from here is a fixed description plus
non-secret detail: the exception's class or the HTTP status.
"""

import json
from typing import Any, Dict, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


class NotificationError(Exception):
    """A notification couldn't be sent. The message holds no secrets."""


def check_webhook_url(url: str, key: str) -> None:
    """Refuse a webhook URL that can never work, naming ``key``, not the URL.

    It must be a well-formed ``https://`` URL with no whitespace or control
    characters.
    """
    text = str(url)
    if any(c.isspace() or ord(c) < 32 or ord(c) == 127 for c in text):
        raise NotificationError(
            f"{key} is not a valid https:// URL (it has whitespace or "
            "control characters)"
        )
    try:
        parts = urlsplit(text)
        valid = parts.scheme == "https" and bool(parts.hostname)
    except ValueError:
        valid = False
    if not valid:
        raise NotificationError(f"{key} is not a valid https:// URL")


def post_json(
    url: str,
    payload: Dict[str, Any],
    label: str,
    headers: Optional[Dict[str, str]] = None,
    ok_statuses=(200,),
) -> bytes:
    """POST ``payload`` as JSON to ``url`` and return the response body.

    Raises :class:`NotificationError` naming ``label`` on any failure,
    without the URL or the headers' values.
    """
    data = json.dumps(payload).encode("utf-8")
    try:
        request = Request(
            url,
            data=data,
            headers={"Content-Type": "application/json", **(headers or {})},
        )
        with urlopen(request, timeout=30) as response:
            status = response.status
            body = response.read()
    except HTTPError as e:
        raise NotificationError(f"{label} failed: HTTP {e.code}") from None
    except URLError as e:
        # A socket error's own text (not the URL): "Name or service not known"
        reason = e.reason
        detail = (
            reason.strerror
            if isinstance(reason, OSError) and reason.strerror
            else type(e).__name__
        )
        raise NotificationError(f"{label} connection failed: {detail}") from None
    except Exception as e:  # pylint: disable=broad-except
        # ValueError, http.client.InvalidURL, socket and TLS errors and the
        # like: their text can quote the URL or a header, so only the
        # class is reported.
        raise NotificationError(f"{label} failed ({type(e).__name__})") from None
    if status not in ok_statuses:
        raise NotificationError(f"{label} returned HTTP {status}")
    return body
