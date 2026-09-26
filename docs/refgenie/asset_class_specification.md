<h1>Asset class specification</h1>

[TOC]

# Introduction

Asset classes are data models that define assets. They specify a list of asset attributes of interest, which correspond to information that a user may want to "look up" (or *seek*) for that asset. An asset would be considered an instance of a class

# Motivation

Recipes used to be tightly coupled to assets. For instance, to get a fasta asset, you must build the fasta recipe. This means a namespace could only have one fasta asset (_e.g._ hg38/fasta), and this asset was used for all other assets that require that type of asset as input.

Generally speaking, there was no way to have more than one asset of the same type under one namespace, which should be allowed so that recipes can require multiple inputs of the same type. To do this, the tight 1:1 coupling between recipes and assets has been relaxed by introducing the concept of _asset classes_.

# Overview

Before we start, let's take a moment to understand what an asset class is and what it's components are. Additionally, we'll look at an example asset class definition that will be used throughout this document.

## Overview of asset class components

A refgenie asset class is a `yaml` file that may consist of the following keys:

- `name`: (REQUIRED) - A string identifying the asset class.
- `version`: (REQUIRED) - A string identifying the version of the asset class.
- `description`: (REQUIRED) - A string describing the asset class, which may include a brief description of the files the class requires and other relevant information.
- `seek_keys`: (REQUIRED) - A dictionary of key names and seek key schemas that include template values that will be used to seek for files in the asset based on the asset class definition.
- `parents`: (OPTIONAL) - A list of asset class names that are the parents of this asset class.

## Example asset class

```yaml
name: my_asset
version: 0.0.1
description: A collection of a FASTA, JSON and HTML file.
seek_keys:
  fasta:
    value: "{genome}.fa.gz"
    type: file
    description: A FASTA file
  html:
    value: "{genome}.html"
    type: file
    description: An HTML file
  json:
    value: "{genome}.json"
    type: file
    description: A JSON file
parents: []
```

This GitHub repository contains numerous asset class examples: [refgenie/recipes](https://github.com/refgenie/recipes/tree/master/asset_classes)

# Asset class inheritance

Asset classes can inherit from other asset classes. This allows for asset classes to extend one or more other asset classes.

Each child asset class will _update_ the seek keys or values of the parent asset classes. Therefore every child asset class always has the same seek keys as its parent asset classes and possibly more.

Let's consider the following example:

**Parent asset class:**

```yaml
name: my_asset_parent
version: 0.0.1
description: A collection of a FASTA and JSON file.
seek_keys:
  fasta:
    value: "{genome}.fa.gz"
    type: file
    description: A FASTA file
  json:
    value: "{genome}.json"
    type: file
    description: A JSON file
parents: []
```

**Child asset class:**

```yaml
name: my_asset_child
version: 0.0.1
description: A collection of a FASTA, JSON and HTML file.
seek_keys:
  fasta:
    value: "{genome}.fa.gz"
    type: file
    description: A FASTA file
  json:
    value: "{genome}_child.json"
    type: file
    description: A JSON file
  html:
    value: "{genome}.html"
    type: file
    description: An HTML file
parents:
  - my_asset_parent
```

**The effective seek keys of the child asset class are:**

- fasta: `{genome}.fa.gz` -- seek key and value are inherited unchanged
- json: `{genome}_child.json` -- seek key value is altered
- html: `{genome}.html` -- seek key and value are specific to the child asset class

# Details of the asset class components

## name

The asset class name is an arbitrary string that identifies the asset class, therefore it should be unique among your asset classes.

## version

The asset class version is an arbitrary string that identifies the version of the asset class. It is used to distinguish between different versions of an evolving asset class.

## description

The asset class description is a freeform text that describes the asset class, which may include a brief description of the seek keys defined by the the asset class and other relevant information.

## seek_keys

The seek keys are the most important part of an asset class definition. They define the files or values that can be looked up in the asset.
The seek keys are a dictionary of key names and seek key schemas that include template values that will be used to locate files in the asset based on the asset class definition.

Each seek key schema must include:

- **value**: a template string that will be populated based on the files existing in the asset directory and the seek key type:
  - For `file` type: The template (e.g., `{genome}.fa.gz`) is converted to a glob pattern by replacing template variables with wildcards (e.g., `*.fa.gz`). The first and only file matching this pattern in the asset directory is used as the value. If no file or multiple files match, an error is raised.
  - For `prefix` type: The value is determined by finding the common prefix among all files in the asset directory (excluding certain files, such as build stats). The template is not used for globbing, but the prefix is computed from the actual files present.
  - For `directory` type: The template is used as a literal directory name. The existence of a directory with this name is checked in the asset directory. The value is the directory name if it exists.
- **type**: a string that specifies the type of the seek key. The type must be either a [jsonschema basic type](http://json-schema.org/understanding-json-schema/reference/type.html#type-specific-keywords) or a special type (`file`, `directory` or `prefix`) is used. The special types are used to specify that the value is a file, a directory or a prefix of a file. This way the file existence can be checked after the asset is built.
- **description**: a freeform text that describes the seek key

## parents

The `parents` key is a list of asset class names and/or definition sources (absolute file path or URL) of the parents of this asset class. **If parents are specified**, then either these parent asset classes must be already managed by refgenie, or the definition files must be present in the same directory as the child asset class definition file.

## Managing Asset Classes via CLI

Refgenie provides CLI commands to manage asset classes. All commands are available under the `asset_class` subcommand.

### Listing asset classes

To see all registered asset classes:

```bash
refgenie1 asset_class list
```

This displays a table of all asset classes currently managed by refgenie.

### Adding an asset class

To register a new asset class from a YAML file:

```bash
refgenie1 asset_class add /path/to/asset_class.yaml
```

You can also add an asset class from a URL:

```bash
refgenie1 asset_class add https://example.com/asset_classes/my_asset.yaml
```

If an asset class with the same name already exists, use `--force` to overwrite it:

```bash
refgenie1 asset_class add /path/to/asset_class.yaml --force
```

### Showing asset class details

To view the YAML definition of a registered asset class:

```bash
refgenie1 asset_class show fasta
```

You can show multiple asset classes at once:

```bash
refgenie1 asset_class show fasta bowtie2_index
```

### Removing an asset class

To remove a registered asset class:

```bash
refgenie1 asset_class remove my_asset_class
```

You can remove multiple asset classes at once:

```bash
refgenie1 asset_class remove asset_class1 asset_class2
```
