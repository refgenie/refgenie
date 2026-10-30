<h1>Recipe specification</h1>

[TOC]

# Introduction

# Overview

Let's start with a brief description of a recipe and it's components. Additionally, we'll look at an example recipe definition that will be used throughout this document.

## Overview of recipe components

A refgenie recipe is a `yaml` file that may consist of the following keys:

- `name`: (REQUIRED) - A string identifying the recipe.
- `version`: (REQUIRED) - A [semantic version](https://semver.org/) string identifying the version of the recipe, like `0.1.0`.
- `output_asset_class`: (REQUIRED) - The asset class of the asset the recipe produces, as `name` or `name:version`.
- `description`: (REQUIRED, may be `null`) - A string describing the recipe, which may include a brief description of the outputs of the recipe, software used to generate the outputs, and other relevant information.
- `input_assets`, `input_files`, `input_params`: (REQUIRED, may be `null`) - Dictionaries describing the assets, files and parameters the recipe takes as inputs.
- `command_templates`: (REQUIRED) - A list of strings, each of which is a [Jinja2](https://jinja.palletsprojects.com/en/3.0.x/) command template for a particular command. This basically means each string is a shell command with variables encoded like `{{values.variable}}`.
- `default_asset`: (REQUIRED) - A string, which is the default name to use for the asset produced by the recipe. The string can be a Jinja2 template, which will be evaluated to determine the name.
- `docker_image`: (REQUIRED, may be `null`) - a container registry path, which can be used to run the recipe.
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
docker_image: docker.io/databio/refgenie
command_templates:
  - cp {{values.files["html"]}} {{values.output_folder}}/{{values.genome_digest}}.html
  - cp {{values.genome_folder}}/{{values.assets["fasta"].seek_keys_dict["fasta"]}} {{values.output_folder}}/{{values.genome_digest}}.fa
  - jsonMaker --cores {{values.params["cores"]}} --input-file {{values.output_folder}}/{{values.genome_digest}}.fa
custom_seek_keys:
  version: jsonMaker --version
  settings: jsonMaker --settings
default_asset: "{{values.custom_seek_keys.version}}"
```

This GitHub repository contains numerous recipe examples: [refgenie/recipes](https://github.com/refgenie/recipes/tree/master/recipes)

## Managing Recipes via CLI

Refgenie provides several CLI commands for managing recipes. These commands allow you to list, add, inspect, and remove recipes from your refgenie configuration.

### Listing recipes

To see all registered recipes:

```bash
refgenie recipe list
```

This displays a table of all recipes currently available in your refgenie installation.

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

Most users do not add recipes one by one: syncing a data channel registers its recipes and asset classes for you (see [Set up a data channel](set_up_data_channel.md)).

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

This displays a table showing all required inputs, which is useful before building an asset to understand what you need to provide.

# Details of the recipe components

Let's describe the components of the recipe specification in more detail, i.e. their internal structure and the impact they have on the asset build process.

## name

The recipe name is an arbitrary string that identifies the recipe, therefore it should be unique among your recipes. When no recipe is explicitly specified (with `refgenie build --recipe-name`), refgenie uses the recipe whose `name` matches the asset group being built.

## version

The recipe version identifies the version of the recipe. It is used to distinguish between different versions of an evolving recipe. It must be a valid [semantic version](https://semver.org/), like `0.1.0`; refgenie rejects recipes with other version strings. When you do not ask for a version (for example with `refgenie build --recipe-version`), refgenie uses the latest one.

## output_asset_class

The asset class of the asset this recipe builds, for example `bowtie2_index`. Give it as `name` to use the latest registered version of the asset class, or as `name:version` to pin one. The asset class must be registered before the recipe is added.

## description

The recipe description is a freeform text that describes the recipe, which may include a brief description of the outputs of the recipe, software used to generate the outputs, and other relevant information. Refgenie uses it as the `description` of the asset that is created by the recipe, unless you provide one with `refgenie build --asset-description`.

## docker_image

The container registry path, which can be used to run the recipe. When you build with `refgenie build --docker`, the commands run in this image.

## Inputs

Inputs are the parameters that are used to generate the outputs of the recipe. They are specified in three top-level keys: `input_assets`, `input_files` and `input_params`. Each is a dictionary keyed by input name, or `null` if the recipe takes no inputs of that kind. Let's take a look at each of them.

### input_assets

This section is used to specify the assets that are used by the recipe. Refgenie will make sure the assets specified in this section are already managed by refgenie before running the recipe (`refgenie build --pull-parents` pulls missing ones). Assets that require other assets as inputs are referred to as _derived assets_.

There are 3 required keys for each input asset:

- `description`: A string briefly describing the required asset
- `asset_class`: A string that identifies the asset class of the asset, as `name` or `name:version`. The asset class must be registered when the recipe is added.
- `default`: The asset group, in the same genome, to use for a build if none is specified. Refgenie uses that group's default asset.

Users can choose a different input asset with `refgenie build --assets fasta=hg38/fasta:default`.

### input_files

This section is used to specify the files that are used by the recipe. Users provide them with `refgenie build --files html=/path/to/file.html`.

Each input file has these keys:

- `description`: A string briefly describing the required file.
- `default`: (optional) A default path. A file with no default must be provided at build time.

### input_params

This section is used to specify the parameters that are used by the recipe. Technically, this a way to pass any text data to the recipe from the command line, like max number of cores to use, simple settings, etc. Users provide them with `refgenie build --params cores=4`.

Each input parameter has these keys:

- `description`: A string briefly describing the required parameter
- `default`: (optional) The default value to use for a build if none is specified. A parameter with no default must be provided at build time.

## command_templates

The command templates are the most critical part of the recipe. They are a list of strings, each of which is a [Jinja2](https://jinja.palletsprojects.com/en/3.0.x/) command template for a particular command. The command templates are evaluated in the order they are specified in the list.

Within the `command_templates`, you have access to a single object named `values`, which holds variables from several sources. You can access them in the command template using the double-brace Jinja2 template language syntax: `{{values.variable}}`, or `{{values.variable["key"]}}` for dictionaries. Templates are rendered strictly, so referring to a variable that does not exist is an error. Below is a list of the variables available to the command templates:

- genome_digest: a `str` representing the digest uniquely identifying the genome
- asset_group_name: a `str` representing the asset group name
- output_folder: a `pathlib.Path` object indicating the output folder; write the asset's files here
- genome_folder: a `pathlib.Path` object indicating the genome folder
- refget_store_path: a `str` path to the RefgetStore that holds the genome sequences
- assets: a `dict[str, Asset]` mapping input asset names to the resolved input assets: specified by the user, or defaults. Use `seek_keys_dict` to get a file path relative to `genome_folder`, like `{{values.genome_folder}}/{{values.assets["fasta"].seek_keys_dict["fasta"]}}`
- params: a `dict[str, Any]` mapping of parameter names to values: specified by the user, or defaults
- files: a `dict[str, Path]` mapping of file names to file paths: specified by the user, or defaults
- custom_seek_keys: a `dict[str, Any]` mapping of custom seek key values, may include software versions, etc.

This is a dict-like representation of an example `values` object available to the command templates:

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

This section is used to determine the name that the asset will get when it is built. The string can be a Jinja2 template, which will be evaluated to determine the name. Users can override this value by specifying an asset name when building, as in `refgenie build hg38/my_asset:my_name`.
The `default_asset` value can be a plain string or a Jinja2 template that combines other `values` into a name. It is reasonable to use `{{values.custom_seek_keys.version}}`, as presented in the example recipe at the top of this page, so assets are named after the version of the tool that built them. If the template renders to an empty string, the build stops and asks you to name the asset yourself. Recipes with no tool version to report can use a literal name, like `default_asset: "default"`.

## custom_seek_keys

This section is an object where keys are the names of custom seek keys and values are shell commands. Refgenie runs each command and uses its output as the value. This is a way to record any computing environment specific information, like the version of a software package used to produce the asset. The results are available to `command_templates` and `default_asset` as `values.custom_seek_keys`, and are stored as seek keys on the built asset. Please note that the other `values` are not available to the commands in this section. If the recipe has a `docker_image`, the commands run inside that image during a `--docker` build.

## inherent

This optional section is an ordered list of gitignore-style file patterns. It says which of the files the recipe produces are part of the asset's identity, so the asset's content digest is computed over only those files. The last matching pattern wins, a `!` prefix marks matching files as not part of the identity, and a file that no pattern matches is not part of it either. If the recipe leaves out `inherent`, every file counts. For example, `inherent: ["*", "!*.log"]` counts every file except log files.
