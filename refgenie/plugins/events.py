"""The event sink managers record into, and the ``@update_scope`` marker.

Managers never load or call plugins. They emit pull/build events and record
committed changes here; the sink hands them to one dispatcher, which the
``Refgenie`` root supplies (``PluginHost``). A manager built without a root
gets ``NULL_EVENTS`` and behaves exactly as it did before plugins existed.
"""

import functools
import threading
from collections.abc import Callable
from contextlib import contextmanager

from refgenie.plugins.hooks import POST_UPDATE, Change, HookEvent

Dispatch = Callable[[HookEvent], None]


class EventSink:
    """Collects events from managers and hands them to one dispatcher.

    ``emit(event)`` dispatches at once (pre/post pull and build).
    ``record(change)`` buffers until the outermost ``scope()`` on this thread
    exits, then dispatches one ``post_update`` carrying every change.
    With no dispatcher (the default) everything is a no-op.
    """

    def __init__(self, dispatch: Dispatch | None = None):
        self._dispatch = dispatch
        # Depth and pending changes are per thread: the dash runs jobs on
        # worker threads that share one Refgenie, and one job's scope must
        # never swallow another's changes.
        self._local = threading.local()

    def _depth(self) -> int:
        return getattr(self._local, "depth", 0)

    def _pending(self) -> list[Change]:
        pending = getattr(self._local, "pending", None)
        if pending is None:
            pending = self._local.pending = []
        return pending

    def emit(self, event: HookEvent) -> None:
        """Dispatch an event now."""
        if self._dispatch is not None:
            self._dispatch(event)

    def record(self, change: Change) -> None:
        """Record a committed change; it is dispatched when the scope closes."""
        if self._dispatch is None:
            return
        self._pending().append(change)
        if self._depth() == 0:  # recorded outside any scope: flush now
            self._flush()

    @contextmanager
    def scope(self):
        """Group changes so ``post_update`` fires once, when the outermost scope exits."""
        self._local.depth = self._depth() + 1
        try:
            yield
        finally:
            self._local.depth -= 1
            if self._local.depth == 0:
                # Also on exception: committed changes are real.
                self._flush()

    def _flush(self) -> None:
        pending = self._pending()
        changes = tuple(pending)
        pending.clear()
        if changes and self._dispatch is not None:
            self._dispatch(HookEvent(hook=POST_UPDATE, changes=changes))


def update_scope(method):
    """Mark a method that may change local assets, aliases or genomes.

    Nested calls share one scope, so ``post_update`` fires once, when the
    outermost marked call returns or raises. Requires ``self._events`` (an
    ``EventSink``).
    """

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        with self._events.scope():
            return method(self, *args, **kwargs)

    return wrapper


#: Shared no-op sink, the default for managers constructed without a root.
NULL_EVENTS = EventSink()
