"""Find, load, disable and describe installed plugins.

Discovery is lazy: nothing here runs at import. The first hook that fires makes
one ``entry_points()`` scan, which is cached for the life of the process, and
loads only that hook's plugins.
"""

import functools
from collections.abc import Callable
from dataclasses import dataclass
from importlib.metadata import EntryPoint, entry_points
from typing import Literal

from refgenie.config import config
from refgenie.logger import logger
from refgenie.plugins.hooks import ENTRY_POINT_PREFIX, HOOKS

#: (hook, plugin name) -> the error message from a failed load.
_LOAD_ERRORS: dict[tuple[str, str], str] = {}
#: hook -> loaded (name, function) pairs, in name order.
_LOADED: dict[str, list[tuple[str, Callable]]] = {}

STATUS_OK = "ok"
STATUS_DISABLED = "disabled"
STATUS_UNKNOWN_HOOK = "unknown hook (never fires)"


@dataclass(frozen=True)
class PluginInfo:
    """One installed entry point, as ``refgenie plugins`` shows it."""

    hook: str
    name: str
    target: str
    dist: str | None
    version: str | None
    status: str


@functools.cache
def _discover() -> dict[str, list[EntryPoint]]:
    """Every ``refgenie.hooks.*`` entry point, keyed by hook suffix.

    Unknown suffixes (legacy hooks such as ``pre_tag``) are kept so
    ``describe`` can show them; ``load`` never loads them.
    """
    found: dict[str, list[EntryPoint]] = {}
    for ep in entry_points():
        if ep.group.startswith(ENTRY_POINT_PREFIX):
            hook = ep.group[len(ENTRY_POINT_PREFIX) :]
            found.setdefault(hook, []).append(ep)
    for eps in found.values():
        eps.sort(key=lambda ep: ep.name)
    return found


def disabled() -> set[str] | Literal["all"]:
    """The plugin names ``REFGENIE_DISABLE_PLUGINS`` turns off, or ``"all"``."""
    value = (config.disable_plugins or "").strip()
    if value.lower() in ("", "0", "false"):
        return set()
    if value.lower() in ("1", "true", "all"):
        return "all"
    return {name.strip() for name in value.split(",") if name.strip()}


def _is_disabled(name: str, off: set[str] | Literal["all"]) -> bool:
    return off == "all" or name in off


def load(hook: str) -> list[tuple[str, Callable]]:
    """The loaded plugins for one hook, in name order.

    Disabled plugins are skipped. A plugin that fails to load is logged,
    skipped for the rest of the process, and reported by ``describe``.
    """
    if hook in _LOADED:
        return _LOADED[hook]
    off = disabled()
    loaded: list[tuple[str, Callable]] = []
    if off != "all" and hook in HOOKS:
        for ep in _discover().get(hook, []):
            if _is_disabled(ep.name, off):
                continue
            try:
                loaded.append((ep.name, ep.load()))
            except Exception as e:
                logger.warning(
                    f"refgenie plugin '{ep.name}' ({ep.value}) failed to load for hook "
                    f"'{hook}': {e}"
                )
                logger.debug("Plugin load traceback", exc_info=True)
                _LOAD_ERRORS[(hook, ep.name)] = str(e)
    _LOADED[hook] = loaded
    return loaded


def describe() -> list[PluginInfo]:
    """Every installed entry point and whether it loads. Loads known plugins."""
    off = disabled()
    infos = []
    for hook, eps in sorted(_discover().items()):
        if hook in HOOKS:
            load(hook)
        for ep in eps:
            if hook not in HOOKS:
                status = STATUS_UNKNOWN_HOOK
            elif _is_disabled(ep.name, off):
                status = STATUS_DISABLED
            elif (hook, ep.name) in _LOAD_ERRORS:
                status = f"load error: {_LOAD_ERRORS[(hook, ep.name)]}"
            else:
                status = STATUS_OK
            dist = ep.dist
            infos.append(
                PluginInfo(
                    hook=hook,
                    name=ep.name,
                    target=ep.value,
                    dist=dist.name if dist is not None else None,
                    version=dist.version if dist is not None else None,
                    status=status,
                )
            )
    return infos


def clear_cache() -> None:
    """Forget discovery, loaded plugins and load errors. Used by tests."""
    _discover.cache_clear()
    _LOADED.clear()
    _LOAD_ERRORS.clear()
