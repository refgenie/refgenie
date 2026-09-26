# Write a refgenie plugin

A plugin is a small Python package that refgenie calls when something happens: an asset is pulled, a build finishes, or local assets change. Use one to keep another tool in sync with refgenie. For example, [refgenie_nfcore](https://github.com/databio/refgenie_nfcore) rewrites a Nextflow genomes config every time your assets change.

## How it works

Your package registers a function in an entry-point group named `refgenie.hooks.<hook>`. When that hook fires, refgenie calls the function with two arguments:

```python
def my_hook(rg, event):
    ...
```

- `rg` is the `Refgenie` instance that fired the hook. Use it to read paths and settings.
- `event` is a `refgenie.plugins.HookEvent` that says what happened.

The return value is ignored. A plugin cannot stop an operation.

## The hooks

| Hook | When it fires | `HookEvent` fields that are set |
| --- | --- | --- |
| `pre_pull` | before each asset pull, bulk pulls included | `genome`, `asset_group`, `asset` (if the caller named one) |
| `post_pull` | after each asset pull, whether it worked or not | `genome`, `asset_group`, `asset`, `succeeded` |
| `pre_build` | before a build | `genome`, `asset_group`, `asset` (if the caller named one) |
| `post_build` | after a build, whether it worked or not | `genome`, `asset_group`, `asset`, `succeeded` |
| `post_update` | once, after any command that changed local assets, aliases or genomes | `changes` |

`event.genome` is the genome as the caller named it, usually an alias.

`post_update` carries `event.changes`, a tuple of `refgenie.plugins.Change` records. Each has an `action` (`asset_added`, `asset_removed`, `asset_renamed`, `default_changed`, `genome_added`, `genome_removed`, `alias_added` or `alias_removed`) and the fields that apply: `genome` (a digest), `asset_group`, `asset`, `previous`, `digest` and `alias`. A bulk pull of ten assets fires ten `pre_pull`/`post_pull` pairs and one `post_update` listing every change.

## Register the function

In your package's `pyproject.toml`:

```toml
[project]
name = "refgenie_myplugin"
dependencies = ["refgenie>=1.0.0a1"]

[project.entry-points."refgenie.hooks.post_update"]
myplugin = "refgenie_myplugin:on_update"
```

The name on the left (`myplugin`) is your plugin's name. Use the same name for your settings (below) and it is also the name to put in `REFGENIE_DISABLE_PLUGINS`.

## A minimal plugin

```python
import logging

from refgenie.plugins import HookEvent

_LOGGER = logging.getLogger(__name__)


def on_update(rg, event: HookEvent) -> None:
    for change in event.changes:
        _LOGGER.info(f"refgenie changed: {change.action} {change.asset_group}:{change.asset}")
```

## Read asset paths

`rg.paths()` is a read-only view of every local asset path, keyed by genome alias, then asset group, then seek key:

```python
paths = rg.paths()
for genome in paths:            # genome aliases
    for group in paths[genome]:  # asset groups that have a default asset
        for seek_key, path in paths[genome][group].items():
            print(genome, group, seek_key, path)
```

To get one path, use `rg.asset.seek_components` with a registry path:

```python
from refgenie.models import AssetRegistryPathComponents

path = rg.asset.seek_components(
    AssetRegistryPathComponents(genome="hg38", asset_group="fasta")
)
```

## Settings

A plugin keeps its settings in the refgenie database, not in environment variables. The user sets them from the command line:

```console
refgenie plugins set myplugin greeting=hi
refgenie plugins unset myplugin greeting
```

The plugin reads them inside its hook:

```python
greeting = rg.plugins.settings("myplugin").get("greeting", "hello")
```

Values are always strings. Use your entry-point name as the settings name.

## Rules

- **Do not raise to stop an operation.** If your plugin raises, refgenie logs a warning with the plugin's name and carries on. Other plugins still run.
- **Hooks do not fire from inside a plugin.** If your plugin pulls or builds, those operations do not call plugins again.
- **Keep `post_update` fast.** It runs at the end of every command that changes local assets.

## Check that it is installed

```console
refgenie plugins
```

This lists every installed plugin with its hook and status. `ok` means it loads. `load error` shows why it does not. `unknown hook (never fires)` means the entry point uses a hook name refgenie does not have, such as a legacy `pre_tag`.

## Turning plugins off

- `REFGENIE_DISABLE_PLUGINS=1` turns off every plugin. `REFGENIE_DISABLE_PLUGINS=myplugin,other` turns off just those.
- Plugins do not run on a server (`refgenie serve`) unless the server sets `REFGENIE_SERVER_PLUGINS=true`. They do run in `refgenie dash`, which is your own local refgenie.

## Examples

- [refgenie_myplugin](https://github.com/databio/refgenie_myplugin): a template to copy. It logs `post_update`, `pre_pull` and `post_build`, and shows how to read a setting.
- [refgenie_nfcore](https://github.com/databio/refgenie_nfcore): keeps a Nextflow genomes config in sync with refgenie.
