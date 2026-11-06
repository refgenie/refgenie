<h1>Recipe specification</h1>

[TOC]

# Overview

A recipe tells refgenie how to build an asset of one [asset class](asset_class_specification.md). It is a YAML file with these keys:

- `name`: (REQUIRED) - A string identifying the recipe.
- `version`: (REQUIRED) - A [semantic version](https://semver.org/) string identifying the version of the recipe, like `0.1.0`.
- `output_asset_class`: (REQUIRED) - The asset class of the asset the recipe produces, as `name` or `name:version`.
- `description`: (REQUIRED, may be `null`) - A string describing the recipe.
- `input_assets`, `input_files`, `input_params`: (REQUIRED, may be `null`) - Dictionaries describing the assets, files and parameters the recipe takes as inputs.
- `command_templates`: (REQUIRED) - A list of [Jinja2](https://jinja.palletsprojects.com/) templates, each a shell command with variables written like `{{values.variable}}`.
- `default_asset`: (REQUIRED) - The default name for the built asset; may be a Jinja2 template.
- `docker_image`: (REQUIRED, may be `null`) - A container image that has the recipe's tools.
- `custom_seek_keys`: (OPTIONAL) - A dictionary of extra seek keys for the asset produced by the recipe, each computed by a shell command.
- `inherent`: (OPTIONAL) - A list of file patterns naming which output files define the asset's identity (its content digest).

Keys marked "may be `null`" must be present in the file, but can be set to `null` when they do not apply.

## Example recipe

```yaml
name: my_asset
version: 0.0.1
output_asset_class: my_asset
description: This asset consists of FASTA, JSON and HTML files.
input_files:
  html:
    description: A HTML file to include
input_params:
  cores:
    description: Number of cores to use
    default: 1
input_assets:
  fasta:
    asset_class: fasta
    description: fasta asset for genome
    default: fasta
docker_image: docker.io/myorg/jsonmaker:1.0.0
command_templates:
  - cp {{values.files["html"]}} {{values.output_folder}}/{{values.genome_digest}}.html
  - cp {{values.genome_folder}}/{{values.assets["fasta"].seek_keys_dict["fasta"]}} {{values.output_folder}}/{{values.genome_digest}}.fa
  - jsonMaker --cores {{values.params["cores"]}} --input-file {{values.output_folder}}/{{values.genome_digest}}.fa
custom_seek_keys:
  version: jsonMaker --version
  settings: jsonMaker --settings
default_asset: "{{values.custom_seek_keys.version}}"
```

The [refgenie-registry](https://github.com/refgenie/refgenie-registry/tree/master/recipes) repository has many real recipes.

## Managing Recipes via CLI

### Listing recipes

To see all registered recipes:

```bash
refgenie recipe list
```

### Adding a recipe

You can add a recipe from a local file or from a URL with `--source`:

```bash
# From a local file
refgenie recipe add --source /path/to/recipe.yaml

# From a URL
refgenie recipe add --source https://example.com/recipes/my_recipe.yaml
```

The recipe's `output_asset_class` and the asset class of every input asset must already be registered. If a recipe with the same name and version already exists, use `--force` to overwrite it:

```bash
refgenie recipe add --source /path/to/recipe.yaml --force
```

Syncing a [data channel](data_channels.md) registers all of its recipes and asset classes at once.

### Showing recipe details

To view the full YAML definition of a recipe:

```bash
refgenie recipe show my_recipe
```

If multiple versions of a recipe exist, refgenie shows the latest one. You can specify a particular version:

```bash
refgenie recipe show my_recipe --recipe-version 0.0.2
```

### Removing a recipe

To remove a recipe from your configuration:

```bash
refgenie recipe remove my_recipe
```

If more than one version of the recipe exists, you must say which one to remove:

```bash
refgenie recipe remove my_recipe --recipe-version 0.0.1
```

A recipe that was used to build an existing asset cannot be removed.

### Checking recipe requirements

To see what inputs (files, parameters, and assets) a recipe requires:

```bash
refgenie recipe requirements my_recipe
```

This shows the files, parameters, and assets you must provide to build it. `refgenie build <genome>/<asset_group> -q` shows the same table.

# Details of the recipe components

## name

The recipe name identifies the recipe and should be unique among your recipes. When no recipe is explicitly specified (with `refgenie build --recipe-name`), refgenie uses the recipe whose `name` matches the asset group being built.

## version

The recipe version tells apart versions of a changing recipe. It must be a valid [semantic version](https://semver.org/), like `0.1.0`; refgenie rejects recipes with other version strings. When you do not ask for a version (for example with `refgenie build --recipe-version`), refgenie uses the latest one.

## output_asset_class

The asset class of the asset this recipe builds, for example `bowtie2_index`. Give it as `name` to use the latest registered version of the asset class, or as `name:version` to pin one. The asset class must be registered before the recipe is added.

## description

Free text describing the recipe, such as its outputs and the software it uses. Refgenie uses it as the `description` of the asset that is created by the recipe, unless you provide one with `refgenie build --asset-description`.

## docker_image

The container image the recipe's tools live in. When you build with `refgenie build --container`, the commands run in this image, with Docker or Apptainer, whichever this computer has.

## Inputs

A recipe's inputs are in three top-level keys: `input_assets`, `input_files` and `input_params`. Each is a dictionary keyed by input name, or `null` if the recipe takes no inputs of that kind.

### input_assets

The assets the recipe uses. Refgenie checks that they exist before running the recipe (`refgenie build --pull-parents` pulls missing ones). An asset built from other assets is a _derived asset_.

Each input asset has three required keys:

- `description`: A string briefly describing the required asset
- `asset_class`: A string that identifies the asset class of the asset, as `name` or `name:version`. The asset class must be registered when the recipe is added.
- `default`: The asset group, in the same genome, to use for a build if none is specified. Refgenie uses that group's default asset.

Users can choose a different input asset with `refgenie build --assets fasta=hg38/fasta:default`.

### input_files

The files the recipe uses. Users provide them with `refgenie build --files html=/path/to/file.html`.

Each input file has these keys:

- `description`: A string briefly describing the required file.
- `default`: (optional) A default path. A file with no default must be provided at build time.

### input_params

Text values passed to the recipe, such as the number of cores or simple settings. Users provide them with `refgenie build --params cores=4`.

Each input parameter has these keys:

- `description`: A string briefly describing the required parameter
- `default`: (optional) The default value to use for a build if none is specified. A parameter with no default must be provided at build time.

## command_templates

A list of [Jinja2](https://jinja.palletsprojects.com/) templates, one per shell command, run in order.

The templates can use one object, `values`, with double-brace syntax: `{{values.variable}}`, or `{{values.variable["key"]}}` for dictionaries. Templates are rendered strictly, so referring to a variable that does not exist is an error. `values` holds:

- genome_digest: a `str` representing the digest uniquely identifying the genome
- asset_group_name: a `str` representing the asset group name
- output_folder: a `pathlib.Path` object indicating the output folder; write the asset's files here
- genome_folder: a `pathlib.Path` object indicating the genome folder
- refget_store_path: a `str` path to the RefgetStore that holds the genome sequences
- assets: a `dict[str, Asset]` mapping input asset names to the resolved input assets: specified by the user, or defaults. Use `seek_keys_dict` to get a file path relative to `genome_folder`, like `{{values.genome_folder}}/{{values.assets["fasta"].seek_keys_dict["fasta"]}}`
- params: a `dict[str, Any]` mapping of parameter names to values: specified by the user, or defaults
- files: a `dict[str, Path]` mapping of file names to file paths: specified by the user, or defaults
- custom_seek_keys: a `dict[str, Any]` mapping of custom seek key values, may include software versions, etc.

An example `values` object:

```python
values = {
    "genome_digest": "pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2",
    "asset_group_name": "my_asset",
    "output_folder": "/path/to/genomes/...",
    "genome_folder": "/path/to/genomes",
    "refget_store_path": "/path/to/genomes/.refget_store",
    "assets": {"fasta": Asset(...)},
    "params": {"cores": "2"},
    "files": {"html": "/path/to/my_asset.html"},
    "custom_seek_keys": {"version": "0.0.1", "settings": "a=1, b=2"},
}
```

## default_asset

The name a built asset gets, unless the user gives one, as in `refgenie build hg38/my_asset:my_name`. It can be a plain string or a Jinja2 template that combines other `values` into a name. A good choice is `{{values.custom_seek_keys.version}}`, as in the example above, so assets are named after the version of the tool that built them. If the template renders to an empty string, the build stops and asks you to name the asset yourself. Recipes with no tool version to report can use a literal name, like `default_asset: "default"`.

## custom_seek_keys

A mapping of extra seek key names to shell commands. Refgenie runs each command and stores its output as the value, which is a way to record details such as the version of the tool that built the asset. The results are available to `command_templates` and `default_asset` as `values.custom_seek_keys`, and are stored as seek keys on the built asset. The other `values` are not available to these commands. If the recipe has a `docker_image`, the commands run inside that image during a `--container` build, and whenever refgenie names an asset outside a build (the web build form, a generated Snakefile).

## inherent

An optional, ordered list of gitignore-style file patterns. It says which of the files the recipe produces are part of the asset's identity, so the asset's content digest is computed over only those files. The last matching pattern wins, a `!` prefix marks matching files as not part of the identity, and a file that no pattern matches is not part of it either. If the recipe leaves out `inherent`, every file counts. For example, `inherent: ["*", "!*.log"]` counts every file except log files.
