"""One-JSON-file-per-record storage for history (#129).

History used to be two YAML files that every qxub process rewrote in full.
Concurrent processes lost each other's records, and a reader that caught the
file mid-rewrite could wipe it. Here each record is its own file:

- A new record gets a name no other process can take (``O_CREAT | O_EXCL``).
- Every write goes to a temporary file that is renamed over the record, so a
  reader sees the old or the new content, never a partial file.
- Only the process that created a record updates it.

No locks are needed, so this does not depend on ``flock``/``fcntl`` working
across nodes on NFS or Lustre. Rename within one directory is atomic on both.

Records hold full command lines, so directories are created ``0700`` and
files ``0600`` whatever the umask. A directory left readable by older qxub is
tightened, with the records in it, the first time it is used.
"""

import json
import logging
import os
import stat
import time
import uuid
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Tuple

logger = logging.getLogger(__name__)

SUFFIX = ".json"
# Temporary files start with "." so listings skip them.
_TMP_PREFIX = "."
_STALE_TMP_SECONDS = 3600
_FILE_MODE = 0o600
_DIR_MODE = 0o700


def make_owner_only(path: Path) -> bool:
    """Remove group and other permissions from ``path`` if it is ours.

    Returns True if the mode was changed. Symlinks are left alone.
    """
    try:
        st = os.lstat(path)
    except OSError:
        return False
    if stat.S_ISLNK(st.st_mode) or st.st_uid != os.getuid():
        return False
    if not st.st_mode & 0o077:
        return False
    try:
        os.chmod(path, stat.S_IMODE(st.st_mode) & 0o700)
    except OSError:
        return False
    return True


def ensure_private_dir(directory: Path) -> bool:
    """Create ``directory`` (and missing parents) owner-only, or tighten it.

    Returns True if an existing directory had to be tightened, meaning files
    in it may still be readable by others.
    """
    directory = Path(directory)
    missing = []
    p = directory
    while not p.exists() and p != p.parent:
        missing.append(p)
        p = p.parent
    for p in reversed(missing):
        try:
            p.mkdir(mode=_DIR_MODE)
        except FileExistsError:
            pass
    return make_owner_only(directory)


class RecordStore:
    """A directory of JSON records, keyed by file name."""

    def __init__(self, directory: Path):
        self.directory = Path(directory)
        self._dir_ready = False

    def _path(self, key: str) -> Path:
        return self.directory / f"{key}{SUFFIX}"

    def _ensure_dir(self) -> None:
        """Create the directory owner-only; tighten it and its records once."""
        if self._dir_ready and self.directory.is_dir():
            return
        if ensure_private_dir(self.directory):
            # Written by older qxub with the umask's mode: tighten the records
            for name in os.listdir(self.directory):
                make_owner_only(self.directory / name)
        self._dir_ready = True

    def write(self, key: str, record: Dict[str, Any]) -> None:
        """Atomically create or replace the record ``key``."""
        self._ensure_dir()
        tmp = self.directory / f"{_TMP_PREFIX}{key}.{uuid.uuid4().hex}.tmp"
        try:
            fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, _FILE_MODE)
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump(record, f, default=str, indent=1)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp, self._path(key))
        except BaseException:
            try:
                tmp.unlink()
            except OSError:
                pass
            raise

    def create_unique(self, key: int, record: Dict[str, Any]) -> str:
        """Write ``record`` under the first free integer key >= ``key``.

        The name is claimed with ``O_EXCL`` before the content is written, so
        two processes that start in the same microsecond get different keys.
        Until the content lands the claimed file is empty, and readers skip it.
        """
        self._ensure_dir()
        while True:
            try:
                fd = os.open(self._path(str(key)), os.O_CREAT | os.O_EXCL, _FILE_MODE)
            except FileExistsError:
                key += 1
                continue
            os.close(fd)
            break
        self.write(str(key), record)
        return str(key)

    def write_if_absent(self, key: str, record: Dict[str, Any]) -> bool:
        """Write ``record`` unless ``key`` exists. Returns True if written.

        Two processes can both find it absent; both then write, and since the
        content only depends on ``key`` it does not matter which one lands.
        """
        if self._path(key).exists():
            return False
        self.write(key, record)
        return True

    def read(self, key: str) -> Optional[Dict[str, Any]]:
        """The record ``key``, or None if missing or not (yet) readable."""
        try:
            with open(self._path(key), encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, ValueError):
            return None
        return data if isinstance(data, dict) else None

    def keys(self) -> List[str]:
        """Keys of all records, in no particular order."""
        try:
            names = os.listdir(self.directory)
        except FileNotFoundError:
            return []
        return [
            n[: -len(SUFFIX)]
            for n in names
            if n.endswith(SUFFIX) and not n.startswith(_TMP_PREFIX)
        ]

    def items_newest_first(self) -> Iterator[Tuple[str, Dict[str, Any]]]:
        """(key, record) for integer-keyed records, newest first."""
        keys = sorted((k for k in self.keys() if k.isdigit()), key=int, reverse=True)
        for key in keys:
            record = self.read(key)
            if record is not None:
                yield key, record

    def delete(self, key: str) -> None:
        try:
            self._path(key).unlink()
        except FileNotFoundError:
            pass

    def prune(
        self,
        keep: int,
        slack: int = 0,
        protect: Optional[Callable[[str], bool]] = None,
    ) -> None:
        """Keep the newest ``keep`` integer-keyed records.

        Nothing happens until there are more than ``keep + slack``, so pruning
        runs now and then rather than on every write. Concurrent pruners may
        try to delete the same file; that is harmless.

        A record for which ``protect(key)`` is true is kept even if it is
        among the oldest; newer records are not deleted in its place.
        """
        keys = [k for k in self.keys() if k.isdigit()]
        if len(keys) <= keep + slack:
            return
        for key in sorted(keys, key=int)[: len(keys) - keep]:
            if protect is not None and protect(key):
                continue
            self.delete(key)
        self._remove_stale_tmp()

    def _remove_stale_tmp(self) -> None:
        """Remove temporary files left by processes that died mid-write."""
        cutoff = time.time() - _STALE_TMP_SECONDS
        try:
            names = os.listdir(self.directory)
        except FileNotFoundError:
            return
        for name in names:
            if name.startswith(_TMP_PREFIX) and name.endswith(".tmp"):
                path = self.directory / name
                try:
                    if path.stat().st_mtime < cutoff:
                        path.unlink()
                except OSError:
                    pass

    def clear(self) -> None:
        for key in self.keys():
            self.delete(key)
