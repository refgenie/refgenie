"""PluginSettingsManager - per-plugin settings stored in the database (``rg.plugins``)."""

import re

from sqlalchemy.engine import Engine

from refgenie.managers.base import ResourceManager
from refgenie.managers.configuration import latest_configuration
from refgenie.plugins import registry
from refgenie.plugins.registry import PluginInfo

_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


def _check_name(kind: str, value: str) -> None:
    if not isinstance(value, str) or not _NAME.match(value):
        raise ValueError(
            f"Invalid plugin {kind} '{value}': use letters, digits, '_', '.' or '-' only."
        )


class PluginSettingsManager(ResourceManager):
    """Per-plugin key/value settings on the configuration row in force. ``rg.plugins``.

    Settings live in one JSON column, ``configuration.plugin_settings``: a dict
    keyed by plugin name whose values are flat ``{key: str}`` dicts. A plugin
    reads its own settings inside its hook with ``rg.plugins.settings(name)``.
    Changing settings is not a change to local assets, so it fires no hook.
    """

    def __init__(self, database_engine: Engine, *, enabled: bool):
        super().__init__(database_engine)
        #: Whether this Refgenie instance dispatches hooks at all.
        self.enabled = enabled

    def settings(self, plugin: str) -> dict[str, str]:
        """A copy of one plugin's settings; ``{}`` when none are stored."""
        return dict(self.all_settings().get(plugin, {}))

    def all_settings(self) -> dict[str, dict[str, str]]:
        """A copy of every plugin's settings."""
        with self._database_session as session:
            stored = latest_configuration(session).plugin_settings or {}
            return {name: dict(values) for name, values in stored.items()}

    def set(self, plugin: str, **values: str) -> dict[str, str]:
        """Merge ``values`` into the plugin's settings and return the result.

        Raises:
            ValueError: If the plugin name or a key is invalid.
        """
        _check_name("name", plugin)
        for key in values:
            _check_name("key", key)
        with self._database_session as session:
            cfg = latest_configuration(session)
            old = dict(cfg.plugin_settings or {})
            merged = {**old.get(plugin, {}), **{k: str(v) for k, v in values.items()}}
            # Assign a new dict: a plain JSON column does not see in-place edits.
            cfg.plugin_settings = {**old, plugin: merged}
            session.add(cfg)
            session.commit()
        return dict(merged)

    def unset(self, plugin: str, *keys: str) -> dict[str, str]:
        """Drop the given keys. With no keys, drop the plugin's entry entirely.

        Returns:
            What is left of the plugin's settings.
        """
        _check_name("name", plugin)
        with self._database_session as session:
            cfg = latest_configuration(session)
            old = dict(cfg.plugin_settings or {})
            remaining = {k: v for k, v in old.get(plugin, {}).items() if k not in keys}
            new = {name: values for name, values in old.items() if name != plugin}
            if keys and remaining:
                new[plugin] = remaining
            cfg.plugin_settings = new
            session.add(cfg)
            session.commit()
        return dict(remaining) if keys else {}

    def installed(self) -> list[PluginInfo]:
        """What is installed, and whether each plugin loads. See ``registry.describe``."""
        return registry.describe()
