"""The `plugins` command group: list installed plugins and manage their settings."""

from collections.abc import Callable

from pydantic import BaseModel, Field
from pydantic_settings import CliPositionalArg, CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import fail
from refgenie.logger import logger


class PluginsListModel(BaseModel):
    """plugins list: show installed plugins, their status and settings."""

    pass


class PluginsSetModel(BaseModel):
    """plugins set: store settings for a plugin."""

    plugin: CliPositionalArg[str] = Field(description="Plugin name (its entry-point name).")
    values: CliPositionalArg[list[str]] = Field(description="Settings, each as key=value.")


class PluginsUnsetModel(BaseModel):
    """plugins unset: remove settings (all of them when no key is given)."""

    plugin: CliPositionalArg[str] = Field(description="Plugin name (its entry-point name).")
    keys: CliPositionalArg[list[str]] = Field(
        default=[], description="Keys to remove. Omit to remove every setting."
    )


class PluginsParser(BaseModel):
    """Intermediate parser for plugins subcommands."""

    list: CliSubCommand[PluginsListModel] = Field(
        description="List installed plugins and their settings."
    )
    set: CliSubCommand[PluginsSetModel] = Field(description="Store settings for a plugin.")
    unset: CliSubCommand[PluginsUnsetModel] = Field(description="Remove settings for a plugin.")


def _dispatch_off_reason(refgenie) -> str | None:
    """Why this instance runs no plugins, or None when it runs them."""
    from refgenie.config import config
    from refgenie.plugins import registry

    if registry.disabled() == "all":
        return f"REFGENIE_DISABLE_PLUGINS={config.disable_plugins}"
    if not refgenie.plugins.enabled:
        return "server mode (set REFGENIE_SERVER_PLUGINS=true to enable)"
    return None


def handle_plugins_list(cmd, refgenie) -> None:
    from refgenie.utils.tables import build_table

    installed = refgenie.plugins.installed()
    if installed:
        rprint(
            build_table(
                "Refgenie plugins",
                ["Hook", "Plugin", "Target", "Package", "Status"],
                [
                    [
                        p.hook,
                        p.name,
                        p.target,
                        f"{p.dist} {p.version}" if p.dist else "",
                        p.status,
                    ]
                    for p in installed
                ],
            )
        )
    else:
        rprint("No refgenie plugins installed.")

    names = {p.name for p in installed}
    rows = [
        [name if name in names else f"{name} (not installed)", key, value]
        for name, values in sorted(refgenie.plugins.all_settings().items())
        for key, value in sorted(values.items())
    ]
    if rows:
        rprint(build_table("Plugin settings", ["Plugin", "Key", "Value"], rows))

    if (reason := _dispatch_off_reason(refgenie)) is not None:
        rprint(f"Plugins are not run here: {reason}")


def _print_settings(plugin: str, values: dict[str, str]) -> None:
    if not values:
        rprint(f"No settings stored for '{plugin}'.")
        return
    for key, value in sorted(values.items()):
        rprint(f"{plugin}: {key}={value}")


def handle_plugins_set(cmd, refgenie) -> None:
    pairs: dict[str, str] = {}
    for item in cmd.values:
        key, sep, value = item.partition("=")
        if not sep or not key:
            fail(f"Settings must look like key=value, got '{item}'.")
        pairs[key] = value
    if cmd.plugin not in {p.name for p in refgenie.plugins.installed()}:
        logger.warning(f"No installed plugin is named '{cmd.plugin}'. Storing the settings anyway.")
    try:
        result = refgenie.plugins.set(cmd.plugin, **pairs)
    except ValueError as e:
        fail(str(e))
    _print_settings(cmd.plugin, result)


def handle_plugins_unset(cmd, refgenie) -> None:
    try:
        result = refgenie.plugins.unset(cmd.plugin, *cmd.keys)
    except ValueError as e:
        fail(str(e))
    _print_settings(cmd.plugin, result)


PLUGINS_DISPATCH: dict[type, Callable] = {
    PluginsListModel: handle_plugins_list,
    PluginsSetModel: handle_plugins_set,
    PluginsUnsetModel: handle_plugins_unset,
}


def handle_plugins_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=False)
    if leaf is None:  # bare `refgenie plugins` lists
        handle_plugins_list(PluginsListModel(), refgenie)
        return
    handler = PLUGINS_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown plugins subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
