"""``PluginHost``: the single place refgenie calls plugin code."""

import threading
from typing import Any

from refgenie.logger import logger
from refgenie.plugins import registry
from refgenie.plugins.hooks import HookEvent


class PluginHost:
    """Calls installed plugins for one Refgenie instance. The only caller of plugin code.

    ``rg`` is typed ``Any`` on purpose: this package sits below the root in the
    layering and must not import it. Plugins use it by duck typing.
    """

    def __init__(self, rg: Any):
        self._rg = rg
        self._local = threading.local()

    def __call__(self, event: HookEvent) -> None:
        if getattr(self._local, "dispatching", False):
            # A plugin that pulls or builds must not re-trigger hooks.
            logger.debug(f"Not dispatching {event.hook} from inside a plugin")
            return
        plugins = registry.load(event.hook)
        if not plugins:
            return
        self._local.dispatching = True
        try:
            for name, func in plugins:
                logger.debug(f"Running {event.hook} plugin: {name}")
                try:
                    func(self._rg, event)
                except Exception as e:
                    logger.warning(
                        f"refgenie plugin '{name}' failed in {event.hook}: {e!r}. "
                        f"Set REFGENIE_DISABLE_PLUGINS={name} to skip it."
                    )
                    logger.debug("Plugin traceback", exc_info=True)
        finally:
            self._local.dispatching = False
