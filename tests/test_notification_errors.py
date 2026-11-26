"""Notification errors never quote the webhook URL or the bot token.

``qxub notify send`` runs in a job whose stderr goes to a file other users
may read, and a Slack or Discord webhook URL is itself the credential.
"""

import pytest
from click.testing import CliRunner

from qxub.config import manager as config_mod
from qxub.notifications import send_all_notifications, webhook
from qxub.notifications.webhook import NotificationError, check_webhook_url
from qxub.notify_cli import notify_cli

SECRET = "T0SECRET/B0SECRET/xyzSECRETxyz"
CONTEXT = {"job_id": "1.x", "exit_code": 0, "job_name": "j", "status": "SUCCESS"}


class FakeConfig:
    def __init__(self, **values):
        self.values = values

    def get_config_value(self, key):
        return self.values.get(key)


def all_text(caplog, *texts):
    return "\n".join([*texts, caplog.text])


@pytest.mark.parametrize(
    "url",
    [
        f"hooks.slack.com/services/{SECRET}",  # no scheme
        f"https://hooks.slack.com/services/{SECRET[:8]} {SECRET[8:]}",  # wrapped paste
        f"http://hooks.slack.com/services/{SECRET}",  # not https
        f"https://hooks.slack.com/services/{SECRET}\n",
        f"https:///services/{SECRET}",  # no host
    ],
)
def test_malformed_webhook_is_refused_without_quoting_it(url, caplog, monkeypatch):
    def no_network(*args, **kwargs):
        raise AssertionError("no request should be made")

    monkeypatch.setattr(webhook, "urlopen", no_network)
    for key in ("notifications.slack.webhook_url", "notifications.discord.webhook_url"):
        with pytest.raises(Exception) as excinfo:
            send_all_notifications(CONTEXT, FakeConfig(**{key: url}))
        message = str(excinfo.value)
        assert key in message
        assert "SECRET" not in all_text(caplog, message)


def test_check_webhook_url_accepts_a_real_one():
    check_webhook_url(f"https://hooks.slack.com/services/{SECRET}", "k")
    check_webhook_url(f"https://discord.com/api/webhooks/1/{SECRET}", "k")


@pytest.mark.parametrize(
    "error",
    [
        ValueError(f"unknown url type: https://hooks.slack.com/{SECRET}"),
        __import__("http.client").client.InvalidURL(f"bad url {SECRET}"),
        OSError(f"socket trouble {SECRET}"),
        RuntimeError(f"anything {SECRET}"),
    ],
)
def test_send_errors_name_only_the_class(error, caplog, monkeypatch):
    def failing_urlopen(*args, **kwargs):
        raise error

    monkeypatch.setattr(webhook, "urlopen", failing_urlopen)
    config = FakeConfig(
        **{
            "notifications.slack.webhook_url": f"https://hooks.slack.com/services/{SECRET}",
            "notifications.discord.webhook_url": f"https://discord.com/api/webhooks/{SECRET}",
        }
    )
    with pytest.raises(Exception) as excinfo:
        send_all_notifications(CONTEXT, config)
    message = str(excinfo.value)
    assert type(error).__name__ in message
    assert "SECRET" not in all_text(caplog, message)


def test_bot_token_never_appears_in_errors(caplog, monkeypatch):
    def failing_urlopen(request, *args, **kwargs):
        # http.client rejects a bad header value by quoting it
        raise ValueError(
            f"Invalid header value {request.get_header('Authorization')!r}"
        )

    monkeypatch.setattr(webhook, "urlopen", failing_urlopen)
    from qxub.notifications.slack import SlackNotifier

    notifier = SlackNotifier(
        FakeConfig(
            **{
                "notifications.slack.bot_token": f"xoxb-{SECRET}\n",
                "notifications.slack.user_id": "U1",
            }
        )
    )
    with pytest.raises(NotificationError) as excinfo:
        notifier.send(CONTEXT)
    assert "SECRET" not in all_text(caplog, str(excinfo.value))


def test_notify_cli_output_has_no_secret(monkeypatch):
    config = FakeConfig(
        **{
            "notifications.slack.webhook_url": f"hooks.slack.com/services/{SECRET}",
            "notifications.discord.webhook_url": f"https://discord.com/api/webhooks/{SECRET} x",
        }
    )
    monkeypatch.setattr(config_mod, "config_manager", config)
    runner = CliRunner(mix_stderr=False) if _has_mix_stderr() else CliRunner()
    for args in (["send", "--job-id", "1.x", "--exit-code", "0"], ["test", "all"]):
        result = runner.invoke(notify_cli, args)
        assert result.exit_code == 1
        output = result.output + (getattr(result, "stderr", "") or "")
        assert "webhook_url is not a valid https:// URL" in output
        assert "SECRET" not in output


def _has_mix_stderr() -> bool:
    import inspect

    return "mix_stderr" in inspect.signature(CliRunner.__init__).parameters
