# Data channels

A *data channel* is a web address (or local folder) that publishes [asset classes](asset_class_specification.md) and [recipes](recipe_specification.md). An asset class says which files make up a kind of asset, such as `fasta` or `bwa_index`. A recipe says how to build one. Refgenie reads both from data channels, so anyone can add new kinds of assets without changing refgenie itself.

A new install has no asset classes or recipes. You need a data channel to build assets. Pulling does not need one, because the server sends the asset class with the asset.

!!! tip "Make your own"
    This page shows how to use a channel. To publish one, see [Set up your own data channel](set_up_data_channel.md).

## The refgenie channel

The refgenie team publishes a channel from the [refgenie-registry](https://github.com/refgenie/refgenie-registry) repository, with classes and recipes for common assets such as `fasta`, `bwa_index`, and `bowtie2_index`. Add and sync it:

```bash
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok
```

`data-channel add` takes a name, a type (`http`, `https`, `ftp`, or `local`), and the address of the channel's `index.yaml`. `sync` downloads every asset class and recipe the channel lists and registers them. Run `sync` again later to get new ones. `--exists-ok` skips items you already have; `--exists-overwrite` replaces them.

Recipes run shell commands on your computer when you build, so add only channels you trust.

## Manage channels

```bash
refgenie data-channel list              # all channels
refgenie data-channel show refgenie     # one channel
refgenie data-channel validate refgenie # check that its index can be read
refgenie data-channel remove refgenie
```

## What a channel contains

A channel is a set of YAML files listed in an `index.yaml`:

```yaml
asset_class:
  dir: asset_classes
  files:
    - fasta.yaml
    - bowtie2_index.yaml
recipe:
  dir: recipes
  files:
    - fasta/recipe.yaml
    - bowtie2_index/recipe.yaml
```

Each file path is relative to its section's `dir`, which is relative to the folder that holds `index.yaml`.

## From Python

```python
from refgenie import Refgenie

r = Refgenie()
r.sources.add_channel(
    name="refgenie",
    type="https",
    index_address="https://refgenie.github.io/refgenie-registry/index.yaml",
)
report = r.sources.sync_channel("refgenie", exists_ok=True)
print(report.asset_classes_added, report.recipes_added)
```

## Channels on a refgenie server

A refgenie server can re-serve the channels it is configured with, at `/v4/data_channel/<channel>/index.yaml`. Add that address like any other channel. See [Serve data channels](server/data_channels.md).
