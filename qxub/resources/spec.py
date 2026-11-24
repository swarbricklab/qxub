"""What a job asks for, independent of any scheduler.

:class:`ResourceSpec` holds parsed values (bytes, seconds, counts). The
scheduler backend turns it into its own syntax: ``PbsBackend`` writes
``-l ncpus=…``. Everything that reads a memory size or a duration uses
:func:`parse_size` and :func:`parse_duration`, so ``--mem 4GB``, a queue's
``max_memory`` and a ``qstat`` field are read the same way.
"""

import re
from dataclasses import dataclass, field, fields, replace
from typing import Any, Dict, Optional, Union

Number = Union[int, float]

_SIZE_UNITS = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
_SIZE_RE = re.compile(r"^(\d+(?:\.\d+)?)([KMGT]?)(B?)$")

_DURATION_UNITS = {"seconds": 1, "minutes": 60, "hours": 3600}
_DURATION_PARTS_RE = re.compile(
    r"^(?:(?P<d>\d+(?:\.\d+)?)d)?(?:(?P<h>\d+(?:\.\d+)?)h)?"
    r"(?:(?P<m>\d+(?:\.\d+)?)m)?(?:(?P<s>\d+(?:\.\d+)?)s)?$"
)
_NUMBER_RE = re.compile(r"^\d+(?:\.\d+)?$")


def parse_size(value: Union[str, Number], bare_unit: str = "MB") -> int:
    """A memory or disk size in bytes.

    Accepts ``4GB``, ``4G``, ``4gb``, ``1.5 GB``, ``512mb``, ``135580kb`` and
    ``4294967296b``. Units are binary (1KB = 1024 bytes). A plain number is
    in ``bare_unit``: qxub's ``--mem 4096`` means MB, while PBS reads
    ``mem=4096`` as bytes. Raises ``ValueError`` for anything else.
    """
    multiplier = _size_multiplier(bare_unit)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return int(value * multiplier)
    text = str(value).replace(" ", "").upper()
    match = _SIZE_RE.match(text)
    if not match:
        raise ValueError(f"invalid size: {value!r}")
    number, prefix, byte = match.groups()
    if prefix or byte:
        multiplier = _SIZE_UNITS[prefix]
    return int(float(number) * multiplier)


def _size_multiplier(unit: str) -> int:
    prefix = unit.upper().rstrip("B")
    if prefix not in _SIZE_UNITS:
        raise ValueError(f"unknown size unit: {unit!r}")
    return _SIZE_UNITS[prefix]


def format_size(num_bytes: int) -> str:
    """``num_bytes`` in the largest unit that divides it exactly: ``4GB``."""
    for unit in ("T", "G", "M", "K"):
        if num_bytes and num_bytes % _SIZE_UNITS[unit] == 0:
            return f"{num_bytes // _SIZE_UNITS[unit]}{unit}B"
    return f"{num_bytes}B"


def parse_duration(value: Union[str, Number], bare_unit: str = "minutes") -> int:
    """A duration in seconds.

    Accepts ``HH:MM:SS``, ``HH:MM``, and unit forms such as ``2h``, ``90m``,
    ``1h30m``, ``45s`` and ``1.5h``. A plain number is in ``bare_unit``
    (``"seconds"``, ``"minutes"`` or ``"hours"``): ``--runtime 90`` means
    minutes, while a platform file's ``max_walltime: 48`` means hours.
    Raises ``ValueError`` for anything else.
    """
    if bare_unit not in _DURATION_UNITS:
        raise ValueError(f"unknown duration unit: {bare_unit!r}")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return int(value * _DURATION_UNITS[bare_unit])
    text = str(value).replace(" ", "").lower()
    if not text:
        raise ValueError("empty duration")
    if ":" in text:
        parts = text.split(":")
        if len(parts) not in (2, 3) or not all(p.isdigit() for p in parts):
            raise ValueError(f"invalid duration: {value!r}")
        numbers = [int(p) for p in parts] + [0] * (3 - len(parts))
        hours, minutes, seconds = numbers
        return hours * 3600 + minutes * 60 + seconds
    if _NUMBER_RE.match(text):
        return int(float(text) * _DURATION_UNITS[bare_unit])
    match = _DURATION_PARTS_RE.match(text)
    if not match or not any(match.groupdict().values()):
        raise ValueError(f"invalid duration: {value!r}")
    scale = {"d": 86400, "h": 3600, "m": 60, "s": 1}
    return int(sum(float(n) * scale[k] for k, n in match.groupdict().items() if n))


def format_duration(seconds: int) -> str:
    """``seconds`` as ``HH:MM:SS``, the form PBS and Slurm both accept."""
    seconds = int(seconds)
    return f"{seconds // 3600:02d}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


@dataclass
class ResourceSpec:
    """The resources a job requests.

    Unset fields are ``None``: the scheduler's or queue's default applies.
    ``native`` holds scheduler-specific requests that have no neutral field
    yet, passed through as ``name=value`` (on Gadi, ``storage=gdata/ab12``).
    """

    cpus: Optional[int] = None
    memory: Optional[int] = None  # bytes
    walltime: Optional[int] = None  # seconds
    gpus: Optional[int] = None
    disk: Optional[int] = None  # bytes of node-local scratch (PBS jobfs)
    native: Dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_options(
        cls,
        *,
        cpus: Any = None,
        memory: Any = None,
        walltime: Any = None,
        gpus: Any = None,
        disk: Any = None,
    ) -> "ResourceSpec":
        """Parse qxub's own options (``--mem 4GB --runtime 2h --cpus 8``).

        Empty values are left unset. A plain number is MB for ``memory`` and
        ``disk``, and minutes for ``walltime``.
        """
        return cls(
            cpus=_count(cpus, "cpus") if cpus else None,
            memory=parse_size(memory) if memory else None,
            walltime=parse_duration(walltime) if walltime else None,
            gpus=_count(gpus, "gpus") if gpus else None,
            disk=parse_size(disk) if disk else None,
        )

    def __bool__(self) -> bool:
        return any(getattr(self, f.name) for f in fields(self))

    def updated(self, other: "ResourceSpec") -> "ResourceSpec":
        """A copy with every value set in ``other`` taking precedence."""
        changes = {
            f.name: getattr(other, f.name)
            for f in fields(self)
            if f.name != "native" and getattr(other, f.name) is not None
        }
        return replace(self, native={**self.native, **other.native}, **changes)

    def as_dict(self) -> Dict[str, Any]:
        """The set values, readable: ``{"cpus": 2, "memory": "4GB", …}``."""
        out: Dict[str, Any] = {}
        if self.cpus is not None:
            out["cpus"] = self.cpus
        if self.gpus is not None:
            out["gpus"] = self.gpus
        if self.memory is not None:
            out["memory"] = format_size(self.memory)
        if self.walltime is not None:
            out["walltime"] = format_duration(self.walltime)
        if self.disk is not None:
            out["disk"] = format_size(self.disk)
        if self.native:
            out["native"] = dict(self.native)
        return out


def _count(value: Any, name: str) -> int:
    try:
        number = int(str(value).strip())
    except ValueError:
        raise ValueError(f"invalid {name}: {value!r}") from None
    if number < 0:
        raise ValueError(f"invalid {name}: {value!r}")
    return number
