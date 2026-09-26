"""Deferred filesystem cleanup, run only after the catalog commits.

Refgenie keeps state in three places that must agree: the SQLite catalog, the
RefgetStore, and plain files on disk. Only the catalog has real transactions, so
it is the participant that decides. Everything else is cleanup that *follows* a
committed decision.

This module is the mechanism for that ordering. A mapper listener that would
otherwise ``shutil.rmtree`` during ``before_delete`` -- inside the flush, before
``COMMIT``, where a rollback restores the row but cannot restore the bytes --
instead records its intent here. The queue is drained by ``after_commit`` and
discarded by ``after_rollback``.

Two rules follow from that:

* **Enqueueing never touches the disk.** The whole point is that nothing
  irreversible happens until the reversible half has succeeded.
* **A failed cleanup is not an error.** By the time the queue runs, the catalog
  already says the data is gone, and that is the answer every reader gets. A
  file left behind is drift for a repair pass to sweep, not a reason to fail a
  command that has already succeeded.

Paths must still be *computed* while the row is present -- the row is where the
path comes from -- so handlers resolve the path as they always did and hand the
resolved path to :func:`enqueue_path_removal`.
"""

import shutil
from pathlib import Path
from collections.abc import Callable
from typing import Literal

from sqlalchemy.orm import Session, object_session

from refgenie.logger import logger

#: Key under which the pending-cleanup queue is stashed on ``Session.info``.
PENDING_KEY = "_pending_cleanups"

RemovalKind = Literal["tree", "file", "empty_dir"]


def _pending(session: Session) -> list:
    return session.info.setdefault(PENDING_KEY, [])


def _resolve_session(target, session: Session | None) -> Session | None:
    if session is not None:
        return session
    return object_session(target)


def enqueue_path_removal(
    target,
    path: Path,
    kind: RemovalKind = "tree",
    *,
    session: Session | None = None,
) -> None:
    """
    Record that ``path`` should be removed once the catalog commits.

    Args:
        target: The ORM object whose deletion motivates the removal. Its session
            is the one the removal is queued on.
        path: The path to remove. Resolve it now, while ``target`` is still
            present -- after the commit the row it was derived from is gone.
        kind: ``"tree"`` for a directory tree (or a symlink to one), ``"file"``
            for a single file or symlink, ``"empty_dir"`` to rmdir a directory
            only if it is empty.
        session: The owning session, when the caller already knows it.
    """
    enqueue_callable(target, lambda: _remove_path(path, kind), session=session)


def enqueue_callable(target, fn: Callable[[], None], *, session: Session | None = None) -> None:
    """
    Record an arbitrary cleanup callable to run once the catalog commits.

    Args:
        target: The ORM object whose write motivates the cleanup.
        fn: A zero-argument callable. It runs outside the transaction, so it
            must not emit SQL.
        session: The owning session, when the caller already knows it.
    """
    resolved = _resolve_session(target, session)
    if resolved is None:
        # No session means no transaction to defer to, so there is nothing this
        # could be ordered after. Run it now rather than dropping it, and say so.
        logger.warning(
            f"No session for {target!r}; running its cleanup immediately rather than "
            f"deferring it to a commit that will never come."
        )
        _run_one(fn)
        return
    _pending(resolved).append(fn)


def _remove_path(path: Path, kind: RemovalKind) -> None:
    """Remove one path. Missing is success: the goal state is 'not there'."""
    if kind == "tree":
        if path.is_symlink():
            path.unlink()
            logger.info(f"Removed symlink: {path}")
        elif path.exists():
            shutil.rmtree(path)
            logger.info(f"Removed directory: {path}")
    elif kind == "file":
        if path.is_file() or path.is_symlink():
            path.unlink()
            logger.info(f"Removed file: {path}")
    elif kind == "empty_dir":
        if path.is_dir() and not any(path.iterdir()):
            path.rmdir()
            logger.info(f"Removed empty directory: {path}")
    else:
        raise ValueError(f"Unknown removal kind: {kind!r}")


def _run_one(fn: Callable[[], None]) -> None:
    try:
        fn()
    except Exception as e:  # noqa: BLE001 - see module docstring
        logger.error(f"Deferred cleanup failed: {e}")


def run_pending_cleanups(session: Session) -> None:
    """Drain and run the queue. Called from ``after_commit``; never raises."""
    pending = session.info.pop(PENDING_KEY, None)
    if not pending:
        return
    logger.debug(f"Running {len(pending)} deferred cleanup(s) after commit")
    for fn in pending:
        _run_one(fn)


def discard_pending_cleanups(session: Session) -> None:
    """Drop the queue without running it. Called from the rollback events."""
    pending = session.info.pop(PENDING_KEY, None)
    if pending:
        logger.debug(f"Discarding {len(pending)} deferred cleanup(s) after rollback")
