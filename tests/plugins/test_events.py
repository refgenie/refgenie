"""EventSink and update_scope: buffering, flushing, and per-thread scopes."""

import threading

import pytest

from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import POST_UPDATE, PRE_PULL, Change, HookEvent


def _sink():
    seen: list[HookEvent] = []
    return EventSink(seen.append), seen


def _change(n: int) -> Change:
    return Change(action="asset_added", asset=f"a{n}")


def test_sink_without_dispatcher_is_a_no_op():
    sink = EventSink()
    sink.emit(HookEvent(hook=PRE_PULL))
    with sink.scope():
        sink.record(_change(1))
    sink.record(_change(2))
    assert sink._pending() == []


def test_null_events_has_no_dispatcher():
    NULL_EVENTS.record(_change(1))
    assert NULL_EVENTS._pending() == []


def test_emit_dispatches_at_once():
    sink, seen = _sink()
    event = HookEvent(hook=PRE_PULL, genome="g")
    sink.emit(event)
    assert seen == [event]


def test_record_outside_a_scope_flushes_at_once():
    sink, seen = _sink()
    sink.record(_change(1))
    assert seen == [HookEvent(hook=POST_UPDATE, changes=(_change(1),))]


def test_nested_scopes_flush_once_in_order():
    sink, seen = _sink()
    with sink.scope():
        sink.record(_change(1))
        with sink.scope():
            sink.record(_change(2))
            with sink.scope():
                sink.record(_change(3))
            assert seen == []
    assert seen == [HookEvent(hook=POST_UPDATE, changes=(_change(1), _change(2), _change(3)))]


def test_a_raising_scope_still_flushes():
    sink, seen = _sink()
    with pytest.raises(ValueError):
        with sink.scope():
            sink.record(_change(1))
            raise ValueError("after the commit")
    assert [e.changes for e in seen] == [(_change(1),)]


def test_a_scope_with_no_records_dispatches_nothing():
    sink, seen = _sink()
    with sink.scope():
        pass
    assert seen == []


def test_scopes_on_two_threads_do_not_share_state():
    sink, seen = _sink()
    inside = threading.Event()
    release = threading.Event()

    def worker():
        with sink.scope():
            sink.record(_change(2))
            inside.set()
            release.wait(5)

    with sink.scope():
        sink.record(_change(1))
        thread = threading.Thread(target=worker)
        thread.start()
        assert inside.wait(5)
    # The main thread's scope closed while the worker's is still open.
    assert [e.changes for e in seen] == [(_change(1),)]
    release.set()
    thread.join(5)
    assert [e.changes for e in seen] == [(_change(1),), (_change(2),)]


def test_update_scope_decorator_uses_self_events():
    sink, seen = _sink()

    class Manager:
        _events = sink

        @update_scope
        def outer(self):
            self.inner()
            sink.record(_change(2))

        @update_scope
        def inner(self):
            sink.record(_change(1))

    Manager().outer()
    assert [e.changes for e in seen] == [(_change(1), _change(2))]
