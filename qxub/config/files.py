"""
Where qxub's config files live, whether to trust them, and writing them safely.

Kept free of heavy imports: both the config manager and the shortcut manager
use it at import time.
"""

import os
import stat
import tempfile
from pathlib import Path
from typing import List, Optional

DEFAULT_XDG_CONFIG_DIRS = "/etc/xdg"


def xdg_config_dirs() -> List[Path]:
    """``$XDG_CONFIG_DIRS`` as absolute directories, in precedence order.

    Empty and relative entries are dropped: they would resolve against the
    working directory, which anyone who can write there controls. Unset,
    empty, or nothing usable left all mean ``/etc/xdg``, as the XDG spec says.
    """
    raw = os.environ.get("XDG_CONFIG_DIRS") or ""
    dirs = [Path(d) for d in raw.split(":") if d and os.path.isabs(d)]
    if not dirs:
        dirs = [Path(DEFAULT_XDG_CONFIG_DIRS)]
    return dirs


def trust_problem(path: Path, require_owner: bool = True) -> Optional[str]:
    """Why ``path`` should not be trusted as config, or None if it is fine.

    Someone other than the invoking user must not be able to change it: it
    may not be group- or world-writable, and with ``require_owner`` it must be
    owned by the invoking user or root. Symlinks are judged by their target.
    """
    try:
        st = os.stat(path)
    except OSError as e:
        return f"cannot stat {path}: {e}"
    if require_owner and st.st_uid not in (os.getuid(), 0):
        return f"{path} is owned by another user (uid {st.st_uid})"
    if st.st_mode & (stat.S_IWGRP | stat.S_IWOTH):
        return f"{path} is writable by group or others"
    return None


def file_trust_problem(path: Path, require_owner: bool = True) -> Optional[str]:
    """:func:`trust_problem` for a config file and the directory holding it."""
    path = Path(os.path.realpath(path))
    return trust_problem(path.parent, require_owner) or trust_problem(
        path, require_owner
    )


def write_private_file(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` readable by its owner only (0600).

    The file is written to a temporary file next to it and renamed into
    place, so readers never see a partial file, a failure leaves the old one
    intact, and an existing looser file is replaced by a 0600 one. A missing
    directory is created owner-only. A symlinked ``path`` keeps its link and
    the target is replaced.
    """
    path = Path(os.path.realpath(path))
    if not path.parent.exists():
        path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.chmod(tmp, 0o600)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
