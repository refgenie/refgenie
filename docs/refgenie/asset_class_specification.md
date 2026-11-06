<h1>Asset class specification</h1>

[TOC]

# Overview

An asset class describes a kind of asset, such as `fasta` or `bowtie2_index`. It lists the asset's *seek keys*: the files, folders, or values a user can look up (*seek*) in it. Each asset is an instance of one class, and each [recipe](recipe_specification.md) builds assets of one class.

An asset class is a YAML file with these keys:

- `name`: (REQUIRED) A string identifying the asset class.
- `version`: (REQUIRED) A [semantic version](https://semver.org/) string, like `0.1.0`.
- `description`: (REQUIRED, may be `null`) A string describing the asset class.
- `seek_keys`: (REQUIRED) A dictionary of seek key names and their definitions.
- `serving_modes`: (OPTIONAL) How a server offers assets of this class: a list of `file`, `archive`, or `none`. Defaults to `[file]`.

## Example asset class

```yaml
name: fasta
version: 0.1.0
description: Sequences in the FASTA format, indexed FASTA and chromosome sizes file
serving_modes: [file, archive]
seek_keys:
  fasta:
    value: "{genome}.fa"
    type: file
    description: FASTA file
  fai:
    value: "{genome}.fa.fai"
    type: file
    description: FASTA index file
  chrom_sizes:
    value: "{genome}.chrom.sizes"
    type: file
    description: Chromosome sizes file
```

The [refgenie-registry](https://github.com/refgenie/refgenie-registry/tree/master/asset_classes) repository has many real asset classes.

# Details

## name

Identifies the asset class. Recipes and assets refer to it by this name.

## version

Tells apart versions of a changing asset class. Refgenie rejects versions that are not valid semantic versions. When a recipe names a class without a version, refgenie uses the latest one.

## description

Free text describing the asset class and its seek keys.

## seek_keys

Each seek key has three keys:

- `type`: one of the types below.
- `value`: a template that tells refgenie how to find the file or folder in a built asset.
- `description`: free text describing the seek key.

The type decides how `value` is used:

| Type | How refgenie finds it |
|---|---|
| `file` | Each `{...}` in the template becomes a wildcard (`{genome}.fa.gz` becomes `*.fa.gz`). Exactly one file in the asset folder must match. |
| `prefix` | The common prefix of all files in the asset folder, as for aligner indexes. The template is not used for matching. |
| `directory` | A folder with exactly the name in `value`, inside the asset folder. |
| `string`, `json` | A value, not a path. The building recipe sets it in `custom_seek_keys`. |

After a build, refgenie checks that every path seek key matches, so a missing or extra file is caught.

## serving_modes

Decides what [staging](staging.md) creates and how clients download the asset. See [Serving modes](serving_modes.md).

# Managing asset classes

```bash
refgenie asset-class list
refgenie asset-class add --source /path/to/asset_class.yaml
refgenie asset-class add --source https://example.com/asset_classes/my_asset.yaml
refgenie asset-class show fasta
refgenie asset-class remove my_asset_class
```

`add --force` replaces an existing class with the same name and version. `show` and `remove` take `--asset-class-version` to choose a version; `remove` needs it when more than one version exists. An asset class that is used by an asset or a recipe cannot be removed.

Syncing a [data channel](data_channels.md) registers all of its asset classes at once.
