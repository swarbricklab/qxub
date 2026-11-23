"""Tests for `qxub cancel all` parsing of `qstat -u` listings."""

import sys
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).parent.parent))

from click.testing import CliRunner  # noqa: E402

from qxub import cancel_cli  # noqa: E402
from qxub.backends import pbs  # noqa: E402

# PBS Pro `qstat -w -u` output in the layout Gadi produces (leading blank line
# included). The Queue column shows execution queues, e.g. `normal-exec`.
QSTAT_W_U = """
gadi-pbs:
                                                                                                   Req'd  Req'd   Elap
Job ID                         Username        Queue           Jobname         SessID   NDS  TSK   Memory Time  S Time
------------------------------ --------------- --------------- --------------- -------- ---- ----- ------ ----- - -----
153392916.gadi-pbs             abc123          normal-exec     align_sample1      12345    1     4   16gb 04:00 R 00:12
153392917.gadi-pbs             abc123          express-exec    qc                    --    1     1    4gb 01:00 Q   --
153392918[].gadi-pbs           abc123          normal-exec     array_job             --    1     1    2gb 00:30 B   --
153392918[3].gadi-pbs          abc123          normal-exec     array_job          999    1     1    2gb 00:30 R 00:01
153392919.gadi-pbs             abc123          normal-exec     held_job              --    1     2    8gb 02:00 H   --
"""


def _listing(output):
    return [(j.job_id, j.queue, j.state) for j in pbs.parse_qstat_user_listing(output)]


def test_parse_skips_headers_and_reads_state_and_queue():
    jobs = _listing(QSTAT_W_U)
    assert jobs == [
        ("153392916.gadi-pbs", "normal-exec", "R"),
        ("153392917.gadi-pbs", "express-exec", "Q"),
        ("153392918[].gadi-pbs", "normal-exec", "B"),
        ("153392918[3].gadi-pbs", "normal-exec", "R"),
        ("153392919.gadi-pbs", "normal-exec", "H"),
    ]


def test_parse_empty_listing():
    assert _listing("") == []
    assert _listing("\ngadi-pbs:\n") == []


def _run_cancel_all(args):
    completed = mock.Mock(stdout=QSTAT_W_U, returncode=0)
    with (
        mock.patch.object(pbs.subprocess, "run", return_value=completed) as run,
        mock.patch.object(cancel_cli, "cancel_by_id") as cancel_by_id,
    ):
        CliRunner().invoke(
            cancel_cli.cancel_all, ["--user", "abc123", "--force", *args]
        )
    return run, cancel_by_id


def _cancelled_ids(cancel_by_id):
    if not cancel_by_id.called:
        return []
    return list(cancel_by_id.call_args.kwargs["job_ids"])


def test_cancel_all_uses_wide_listing_without_queue_summary_flag():
    run, _ = _run_cancel_all(["--queue", "express"])
    assert run.call_args.args[0] == ["qstat", "-w", "-u", "abc123"]


def test_cancel_all_filters_by_status():
    _, cancel_by_id = _run_cancel_all(["--status", "R"])
    assert _cancelled_ids(cancel_by_id) == [
        "153392916.gadi-pbs",
        "153392918[3].gadi-pbs",
    ]


def test_cancel_all_queue_matches_submission_queue_name():
    # Users pass the queue they submitted to; qstat shows its `-exec` twin.
    _, cancel_by_id = _run_cancel_all(["--queue", "express"])
    assert _cancelled_ids(cancel_by_id) == ["153392917.gadi-pbs"]


def test_cancel_all_queue_matches_execution_queue_name():
    _, cancel_by_id = _run_cancel_all(["--queue", "express-exec"])
    assert _cancelled_ids(cancel_by_id) == ["153392917.gadi-pbs"]


def test_cancel_all_queue_is_not_a_prefix_match():
    # A partial name must not match: `norm` selects neither `normal-exec` nor
    # anything else in the listing.
    _, cancel_by_id = _run_cancel_all(["--queue", "norm"])
    assert _cancelled_ids(cancel_by_id) == []


def test_cancel_all_never_targets_header_tokens():
    _, cancel_by_id = _run_cancel_all([])
    ids = _cancelled_ids(cancel_by_id)
    assert len(ids) == 5
    assert all(i[0].isdigit() for i in ids)
