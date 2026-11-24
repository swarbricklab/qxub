"""ResourceSpec, the one size and duration parser, and PBS ``-l`` rendering."""

import pytest
from test_golden_pbs_submission import dry_run_qsub_argv

from qxub.backends.pbs import parse_pbs_resources, pbs_resource_list
from qxub.resources import (
    ResourceSpec,
    format_duration,
    format_size,
    parse_duration,
    parse_memory_size,
    parse_size,
    parse_walltime,
    size_to_bytes,
    time_to_seconds,
)

GB = 1024**3
MB = 1024**2


@pytest.mark.parametrize(
    "text, expected",
    [
        ("4GB", 4 * GB),
        ("4gb", 4 * GB),
        ("4G", 4 * GB),
        ("4 GB", 4 * GB),
        ("1.5GB", int(1.5 * GB)),
        ("400000M", 400000 * MB),
        ("135580kb", 135580 * 1024),
        ("4294967296b", 4 * GB),
        ("2TB", 2 * 1024**4),
        ("1024", 1024 * MB),  # bare number: MB by default
        (512, 512 * MB),
    ],
)
def test_parse_size(text, expected):
    assert parse_size(text) == expected


def test_parse_size_bare_unit():
    assert parse_size("4096", bare_unit="B") == 4096
    assert parse_size(2, bare_unit="GB") == 2 * GB


@pytest.mark.parametrize("text", ["", "GB", "4XB", "four", "-1GB", "4GB8"])
def test_parse_size_rejects(text):
    with pytest.raises(ValueError):
        parse_size(text)


@pytest.mark.parametrize(
    "text, expected",
    [
        ("01:00:00", 3600),
        ("1:00:00", 3600),
        ("100:30:00", 100 * 3600 + 1800),
        ("2:30", 2 * 3600 + 1800),  # two fields: hours:minutes
        ("2h", 7200),
        ("90m", 5400),
        ("1h30m", 5400),
        ("1H30M", 5400),
        ("45s", 45),
        ("1.5h", 5400),
        ("1d", 86400),
        ("90", 5400),  # bare number: minutes by default
        (90, 5400),
    ],
)
def test_parse_duration(text, expected):
    assert parse_duration(text) == expected


def test_parse_duration_bare_unit():
    assert parse_duration("2", bare_unit="hours") == 7200
    assert parse_duration("30", bare_unit="seconds") == 30


@pytest.mark.parametrize("text", ["", "h", "1:2:3:4", "1:xx", "30m1h", "soon"])
def test_parse_duration_rejects(text):
    with pytest.raises(ValueError):
        parse_duration(text)


def test_formatting():
    assert format_size(4 * GB) == "4GB"
    assert format_size(400000 * MB) == "400000MB"
    assert format_size(1536 * MB) == "1536MB"
    assert format_size(1000) == "1000B"
    assert format_duration(5400) == "01:30:00"
    assert format_duration(100 * 3600) == "100:00:00"


def test_older_parsers_use_the_same_rules():
    """The older helpers keep their units but share the parser."""
    assert parse_memory_size("4G") == 4 * GB
    assert parse_memory_size("1024") == 1024 * MB
    assert parse_memory_size("rubbish") is None
    assert parse_walltime("1h30m") == 1.5  # used to be None
    assert parse_walltime("2") == 2.0  # bare: hours
    assert size_to_bytes("135580kb") == 135580 * 1024
    assert size_to_bytes("") == 0
    assert time_to_seconds("01:45:30") == 6330


def test_from_options_and_merging():
    spec = ResourceSpec.from_options(memory="4GB", walltime="1h", cpus="8")
    assert spec == ResourceSpec(cpus=8, memory=4 * GB, walltime=3600)
    assert not ResourceSpec.from_options(memory=None, cpus=0)
    assert bool(spec)

    base = ResourceSpec(cpus=2, memory=GB, native={"storage": "gdata/ab12"})
    merged = base.updated(ResourceSpec(cpus=4, native={"place": "excl"}))
    assert merged == ResourceSpec(
        cpus=4, memory=GB, native={"storage": "gdata/ab12", "place": "excl"}
    )
    assert base.cpus == 2  # unchanged

    assert merged.as_dict() == {
        "cpus": 4,
        "memory": "1GB",
        "native": {"storage": "gdata/ab12", "place": "excl"},
    }


def test_from_options_rejects_bad_values():
    with pytest.raises(ValueError):
        ResourceSpec.from_options(cpus="lots")
    with pytest.raises(ValueError):
        ResourceSpec.from_options(memory="4XB")


def test_parse_pbs_resources():
    spec = parse_pbs_resources(
        ["ncpus=48,mem=128GB", "walltime=2:00:00", "storage=gdata/ab12+scratch/ab12"]
        + ["jobfs=100GB", "ngpus=1", "ncpus=4"]  # later wins
    )
    assert spec == ResourceSpec(
        cpus=4,
        gpus=1,
        memory=128 * GB,
        walltime=7200,
        disk=100 * GB,
        native={"storage": "gdata/ab12+scratch/ab12"},
    )


def test_parse_pbs_resources_follows_pbs_units():
    """A plain PBS number is bytes or seconds; MM:SS is minutes."""
    spec = parse_pbs_resources(["mem=4096", "walltime=30:00"])
    assert (spec.memory, spec.walltime) == (4096, 1800)
    assert parse_pbs_resources(["walltime=3600"]).walltime == 3600


@pytest.mark.parametrize("item", ["ncpus", "ncpus=", "ncpus=four", "mem=lots"])
def test_parse_pbs_resources_rejects(item):
    with pytest.raises(ValueError):
        parse_pbs_resources([item])


def test_pbs_resource_list_order_and_round_trip():
    spec = ResourceSpec(
        walltime=3600, memory=4 * GB, cpus=2, native={"storage": "gdata/ab12"}
    )
    rendered = pbs_resource_list(spec)
    assert rendered == ["ncpus=2", "mem=4GB", "walltime=01:00:00", "storage=gdata/ab12"]
    assert parse_pbs_resources(rendered) == spec


# -- qxub exec ---------------------------------------------------------------


def _l_values(argv):
    return [argv[i + 1] for i, arg in enumerate(argv) if arg == "-l"]


def test_exec_option_overrides_the_same_l_request(sandbox):
    """--cpus wins over -l ncpus=, and qsub gets one ncpus, not two."""
    result = sandbox.run(
        "exec", "--dry-run", "-P", "ab12", "-l", "ncpus=8,mem=8GB",
        "--cpus", "4", "--", "ls",
    )  # fmt: skip
    assert _l_values(dry_run_qsub_argv(result)) == ["ncpus=4", "mem=8GB"]


def test_exec_accepts_compound_runtime(sandbox):
    result = sandbox.run(
        "exec", "--dry-run", "-P", "ab12", "--runtime", "1h30m", "--mem", "1.5GB",
        "--", "ls",
    )  # fmt: skip
    assert _l_values(dry_run_qsub_argv(result)) == ["mem=1536MB", "walltime=01:30:00"]


def test_exec_rejects_a_bad_resource(sandbox):
    result = sandbox.run("exec", "--dry-run", "-P", "ab12", "--mem", "lots", "--", "ls")
    assert result.returncode == 2
    assert "Invalid resource request" in result.stderr

    result = sandbox.run("exec", "--dry-run", "-P", "ab12", "-l", "ncpus", "--", "ls")
    assert result.returncode == 2
    assert "NAME=VALUE" in result.stderr
