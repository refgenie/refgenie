<h1>Asset class specification</h1>

[TOC]

# Introduction

Recipes used to be tightly coupled to assets. For instance, to get a fasta asset, you must build the fasta recipe. This means a namespace could only have one fasta asset (_e.g._ hg38/fasta), and this asset was used for all other assets that require that type of asset as input.

Generally speaking, there was no way to have more than one asset of the same type under one namespace, which should be allowed so that recipes can require multiple inputs of the same type. To do this, the tight 1:1 coupling between recipes and assets has been relaxed by introducing the concept of _asset classes_.

# Overview

Before we start, let's take a moment to understand what an asset class is and what it's components are. Additionally, we'll look at an example asset class definition that will be used throughout this document.

## Overview of asset class components

A refgenie asset class my consist of the following keys:

- `name`: (REQUIRED) - A string identifying the asset class.
- `version`: (REQUIRED) - A string identifying the version of the asset class.
- `description`: (REQUIRED) - A string describing the asset class, which may include a brief description of the files the class requires and other relevant information.
- `seek_keys`: (REQUIRED) - A dictionary of key names and template values that will be used to seek for files in the asset based on the asset class definition.
- `parents`: (OPTIONAL) - A list of asset class names that are the parents of this asset class.

## Example asset class

```yaml
name: my_asset
version: 0.0.1
description: A collection of FASTA, JSON and HTML files.
seek_keys:
  fasta: "{genome}.fa.gz"
  html: "{genome}.html"
  json: "{genome}.json"
parents: []
```

# Asset class inheritance

Asset classes can inherit from other asset classes. This allows for asset classes to extend one or more other asset classes.

Each child asset class will _update_ the seek keys or values of the parent asset classes. Therefore every child asset class always has the same seek keys as its parent asset classes and possibly more.

Let's consider the following example:

**Parent asset class:**

```yaml
name: my_asset_parent
version: 0.0.1
description: A collection of FASTA and JSON files.
seek_keys:
  fasta: "{genome}.fa.gz"
  json: "{genome}.json"
parents: []
```

**Child asset class:**

```yaml
name: my_asset_child
version: 0.0.1
description: A collection of FASTA, JSON and HTML files.
seek_keys:
  html: "{genome}.html"
  json: "{genome}_child.json"
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

The seek keys are the most important part of an asset class definition. They define the files that will be sought for in the asset based on the asset class definition.
The seek keys are a dictionary of key names and template values that will be used to seek for files in the asset based on the asset class definition. This is the list of namespaces that are avilable to use in the seek keys:

- genome (`str`): The genome digest.
- asset (`str`): The asset name.
- params (`Dict[str, Dict[str, str]]`): The _resolved_ input parameters.
- files (`Dict[str, Dict[str, str]]`): The _resolved_ input files.
- assets (`Dict[str, Dict[str, str]]`): The _resolved_ input assets.
- custom*properties (`Dict[str, Dict[str, Any]]`): The \_resolved* custom properties.

## parents

The `parents` key is a list of asset class names and/or definition sources (absolute file path or URL) of the parents of this asset class.

**If names are listed** then these parent asset classes must be either already managed by refgenie or the definition files must be present in the same directory as the child asset class definition file.
