<h1>Recipe specification</h1>

[TOC]

# Introduction

# Overview

Let's start with a brief description of a recipe and it's components. Additionally, we'll look at an example recipe definition that will be used throughout this document.

## Overview of recipe components

A refgenie recipe is a `yaml` file that may consist of the following keys:

- `name`: (REQUIRED) - A string identifying the recipe.
- `version`: (REQUIRED) - A string identifying the version of the recipe.
- `description`: (REQUIRED) - A string describing the recipe, which may include a brief description of the outputs of the recipe, software used to generate the outputs, and other relevant information.
- `inputs`: (REQUIRED) - A nested dictionary of input parameters for the recipe, which includes 3 keys: `files`, `assets` and `params`.
- `command_template_list`: (REQUIRED) - A list of strings, each of which is a [Jinja2](https://jinja.palletsprojects.com/en/3.0.x/) command template for a particular command. This basically means each string is a shell command with allowed variables encoded like `{variable}`.
- `default_tag`: (REQUIRED) - A string, which is the default tag to use for the asset produced by the recipe. The string can be a Jinja2 template, which will be evaluated to determine the tag.
- `container`: (RECOMMENDED) - a container registry path, which can be used to run the recipe.
- `custom_properties`: (OPTIONAL) - A dictionary of custom properties of the asset produced by the recipe.
- `test`: (RECOMMENDED) - A dictionary of test settings.

## Example recipe

```yaml
name: my_asset1
version: 0.0.1
output_asset_class: my_asset
description: This asset consists of FASTA, JSON and HTML files.
inputs:
  files:
    html:
      description: A HTML file to include
  params:
    cores:
      description: Number of cores to use
      default: 1
  assets:
    fasta:
      asset_class: fasta
      description: gzipped fasta file
      default: fasta
container: docker.io/databio/refgenie
command_template_list:
  - cp {{files.html}} {{asset_outfolder}}/{{genome}}.html
  - cp {{assets.fasta}} {{asset_outfolder}}/{{genome}}.fa.gz
  - jsonMaker --cores {{params.cores}} --input-file {{assets.fasta}}
custom_properties:
  version: jsonMaker --version
  settings: jsonMaker --settings
default_tag: "{{custom_properties.version}}"
```

This GitHub repository contains numerous recipe examples: [refgenie/recipes](https://github.com/refgenie/recipes/tree/master/recipes)

# Details of the recipe components

Let's describe the components of the recipe specification in more detail, i.e. their internal structure and the impact they have on the asset build process.

## name

The recipe name is an arbitrary string that identifies the recipe, therefore it should be unique among your recipes. Refgenie tries to match the `name` attribute to the requested asset build, when no recipe is explicitly specified.

## version

The recipe version is an arbitrary string that identifies the version of the recipe. It is used to distinguish between different versions of an evolving recipe.

## description

The recipe description is a freeform text that describes the recipe, which may include a brief description of the outputs of the recipe, software used to generate the outputs, and other relevant information. Refgenie attaches a `description` attribute to the asset that is created by the recipe.

## container

The container registry path, which can be used to run the recipe.

## inputs

Inputs are the parameters that are used to generate the outputs of the recipe. The inputs are specified as a nested dictionary of input parameters, which include: `files`, `assets` and `params`. Let's take a look at each of them.

### assets

This subsection is used to specify the assets that are used by the recipe. Refgnie will make sure the assets specified in this section are already managed by refgenie before running the recipe. Assets that require other assets as inputs are referred to as _derived assets_.

There are 3 required keys in the `assets` subsection:

- `description`: A string briefly describing the required asset
- `asset_class`: A string that identifies the asset class of the asset. Asset class mismatches will result in errors before the recipes are run.
- `default`: A string that identifies the default asset to use for a build if none are specified.

### files

This subsection is used to specify the files that are used by the recipe.

There is only one required key in the `files` subsection:

- `description`: A string briefly describing the required file.

### params

This subsection is used to specify the parameters that are used by the recipe. Technically, this a way to pass any text data to the recipe from the command line, like max number of cores to use, simple settings, etc.

There are 2 required keys in the `params` subsection:

- `description`: A string briefly describing the required parameter
- `default`: A string that sepecifies the default parameter to use for a build if none are specified.

## command_template_list

The command template list is the most critical part of the recipe. It is a list of strings, each of which is a [Jinja2](https://jinja.palletsprojects.com/en/3.0.x/) command template for a particular command. The command templates are evaluated in the order they are specified in the list.

Within the `command_template_list`, you have access to variables from several sources. These variables are divided into namespaces depending on the variable source. You can access the values of these variables in the command template using the double-brace Jinja2 template language syntax: `{{namespace.variable}}`. Below is a list of the namespaces available to the command templates:

- genome (`str`): The genome digest.
- asset (`str`): The asset name.
- params (`Dict[str, Dict[str, str]]`): The _resolved_ input parameters.
- files (`Dict[str, Dict[str, str]]`): The _resolved_ input files.
- assets (`Dict[str, Dict[str, str]]`): The _resolved_ input assets.
- custom*properties (`Dict[str, Dict[str, Any]]`): The \_resolved* custom properties.
- test (`Dict[str, Dict[str, str]]`): The recipe testing setup.

This is a dict-like representation of an example namespaces object available to the command templates:

```python
namespaces = {
    "genome": "6c5f19c9c2850e62cc3f89b04047fa05eee911662bd77905",
    "asset": "my_asset1",
    "assets": {
        "fasta": "/path/to/genomes/fasta.fa.gz"
    },
    "params": {
        "cores": "2"
    },
    "files": {
        "html": "/path/to/my_asset1.html"
    },
    "custom_properties": {
        "version": "0.0.1",
        "settings": "a=1, b=2"
    }
}
```

## default_tag

This section is used to determine the string that the asset will be tagged with when it is built. The string can be a Jinja2 template, which will be evaluated to determine the tag. Users can override this value by specifying a tag when building the asset.
The `default_tag` value can be a string or a Jinja2 template, that combines the values of other namespaces into a tag. It is reasonable to use `{{custom_properties.version}}` as the default tag, as presented in the example recipe at the top of this page.

## custom_properties

This section is an object where keys are the names of custom properties and values are shell commands that will be executed in the asset building requirement. This is a way to record any computing enviroment specific information, like the version of a software package used to produce the asset. Please note that the build namespaces are not available to the commands in this section.

## test

This section is used to specify the test setup for the recipe. The test setup is a dictionary of test parameters, which include:

- `genome`: the genome digest to use for derived asset tests. Any input assets will be sourced from this genome namespace.
- `inputs`: the input files to the recipes. The values must be URLs to remote files, so that the recipes are portable. The keys in this section must match the keys in the `inputs.files` section or the recipe. Currently `inputs.params` and `inputs.assets` cannot be changed, default values from the top-level `inputs` section are used.
- `outputs`: the outputs to test the results of running the recipe commands against. There are two ways the outputs can be tested: 1) by comparing the file checksum (`checksum`), or 2) by comparing the value (`value`). The keys in this section must match the `seek_keys` in the asset class definition.

```yaml
test:
  genome: 6c5f19c9c2850e62cc3f89b04047fa05eee911662bd77905
  inputs:
    fasta: /path/to/genome/fasta.fa
  outputs:
    html:
      md5sum: d41d8cd98f00b204e9800998ecf8427e
      # OR
      value: 100
```
