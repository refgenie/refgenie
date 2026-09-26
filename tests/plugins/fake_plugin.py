"""Plugin functions the plugin tests register through fake entry points."""

import threading

from refgenie.plugins import HOOKS, HookEvent

#: The entry-point target of `record`.
RECORD = "tests.plugins.fake_plugin:record"
BOOM = "tests.plugins.fake_plugin:boom"

#: (hook, event) for every call to `record`, in order.
CALLS: list[tuple[str, HookEvent]] = []
#: The thread each `record` call ran on.
THREADS: list[str] = []
#: What `read_settings` saw in `rg.plugins.settings("fake")`.
SEEN_SETTINGS: list[dict[str, str]] = []


def record(rg, event: HookEvent) -> None:
    CALLS.append((event.hook, event))
    THREADS.append(threading.current_thread().name)


def boom(rg, event: HookEvent) -> None:
    raise RuntimeError("boom")


def interrupt(rg, event: HookEvent) -> None:
    raise KeyboardInterrupt


def reenter(rg, event: HookEvent) -> None:
    """Dispatch again from inside a plugin, as a plugin that pulls would."""
    CALLS.append(("reenter", event))
    rg.dispatch(event)


def read_settings(rg, event: HookEvent) -> None:
    SEEN_SETTINGS.append(rg.plugins.settings("fake"))


def seek_added(rg, event: HookEvent) -> None:
    """Seek every asset a post_update says was added, to prove the tree is there."""
    for change in event.changes:
        if change.action == "asset_added":
            CALLS.append(("seek", rg.asset.seek(change.genome, change.asset_group, change.asset)))


def on_every_hook(name: str = "record", target: str = RECORD) -> dict:
    """An ``install_plugins`` spec registering one function on every hook."""
    return {hook: [(name, target)] for hook in HOOKS}
