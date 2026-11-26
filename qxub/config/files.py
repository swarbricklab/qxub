"""
Where qxub's config files live, whether to trust them, and writing them safely.

Kept free of heavy imports: both the config manager and the shortcut manager
use it at import time.
"""

import errno
import os
import stat
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, TextIO, Tuple

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


def _stat_trust_problem(
    st: os.stat_result, path: Path, require_owner: bool = True
) -> Optional[str]:
    """Why a file with stat ``st`` should not be trusted, or None."""
    if require_owner and st.st_uid not in (os.getuid(), 0):
        return f"{path} is owned by another user (uid {st.st_uid})"
    if st.st_mode & (stat.S_IWGRP | stat.S_IWOTH):
        return f"{path} is writable by group or others"
    return None


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
    return _stat_trust_problem(st, path, require_owner)


def open_trusted_dir(path: Path) -> Tuple[Optional[int], Optional[str]]:
    """Open directory ``path`` for :func:`open_trusted_file`, if it is trusted.

    Returns ``(fd, None)``, or ``(None, problem)``. The directory itself must
    not be a symlink, and the checks of :func:`trust_problem` are made on the
    opened directory, so it can't be swapped between the check and the reads.
    """
    if os.path.islink(path):
        return None, f"{path} is a symbolic link"
    try:
        fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    except OSError as e:
        return None, f"cannot open {path}: {e}"
    st = os.fstat(fd)
    problem = (
        None if stat.S_ISDIR(st.st_mode) else f"{path} is not a directory"
    ) or _stat_trust_problem(st, path)
    if problem:
        os.close(fd)
        return None, problem
    return fd, None


def open_trusted_file(
    dir_fd: int, name: str, path: Path
) -> Tuple[Optional[TextIO], Optional[str]]:
    """Open ``name`` in the directory ``dir_fd`` if it is a trusted file.

    Returns ``(file, None)``; ``(None, None)`` if there is no such file; or
    ``(None, problem)``. A symlink is not followed, and the checks of
    :func:`trust_problem` are made on the opened file, which the caller reads
    (``path`` is only for messages).
    """
    try:
        fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=dir_fd)
    except FileNotFoundError:
        return None, None
    except OSError as e:
        if e.errno == errno.ELOOP:
            return None, f"{path} is a symbolic link"
        return None, f"cannot open {path}: {e}"
    st = os.fstat(fd)
    problem = (
        None if stat.S_ISREG(st.st_mode) else f"{path} is not a regular file"
    ) or _stat_trust_problem(st, path)
    if problem:
        os.close(fd)
        return None, problem
    return os.fdopen(fd, "r", encoding="utf-8"), None


def file_trust_problem(path: Path, require_owner: bool = True) -> Optional[str]:
    """:func:`trust_problem` for a config file and the directory holding it."""
    path = Path(os.path.realpath(path))
    return trust_problem(path.parent, require_owner) or trust_problem(
        path, require_owner
    )


# Config directories already checked by ensure_private_config_dir this process.
_PRIVATE_DIRS_CHECKED: set = set()


def ensure_private_config_dir(directory: Path) -> None:
    """Create qxub's user config ``directory`` owner-only, or tighten it once.

    It holds the job database, history and config that may contain tokens.
    A missing directory is created 0700. An existing one that is ours and not
    a symlink loses its group/other permissions (``chmod go-rwx``), with a
    one-line notice on stderr when that changed anything. Checked once per
    directory per process.
    """
    directory = Path(directory)
    key = str(directory)
    if key in _PRIVATE_DIRS_CHECKED:
        return
    _PRIVATE_DIRS_CHECKED.add(key)
    if not os.path.lexists(directory):
        directory.mkdir(mode=0o700, parents=True, exist_ok=True)
        return
    try:
        st = os.lstat(directory)
        if (
            not stat.S_ISDIR(st.st_mode)
            or st.st_uid != os.getuid()
            or not st.st_mode & 0o077
        ):
            return
        os.chmod(directory, stat.S_IMODE(st.st_mode) & ~0o077)
    except OSError:
        return
    sys.stderr.write(
        f"qxub: removed group/other access from {directory} "
        "(it holds your job history and config)\n"
    )


def write_private_file(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` readable by its owner only (0600).

    The file is written to a temporary file next to it and renamed into
    place, so readers never see a partial file, a failure leaves the old one
    intact, and an existing looser file is replaced by a 0600 one. A missing
    directory is created owner-only. A symlinked ``path`` keeps its link and
    the target is replaced.
    """
    ensure_private_config_dir(Path(path).parent)
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
