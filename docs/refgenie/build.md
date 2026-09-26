<h1>Build</h1>
[TOC]

!!! tip "See also"
    For a complete end-to-end walkthrough covering build, stage, serve, and push, see the [Building and Serving Tutorial](building_tutorial.md).

## Introduction
Once you've [installed refgenie](install.md), you can use `refgenie pull` to [download pre-built assets](pull.md) without installing any additional software. However, you may need to use the `build` function for genomes or assets that are not available on the server. You can build assets for any genome for which you can provide the required inputs.

Building assets is a bit more complicated than pulling them. You'll need to set up 2 things: 1) any *inputs* required by the asset build recipe; 2) Any *software* required by the recipe. Below, we'll walk you through each of these requirements, but first, how can you tell *what* refgenie can build in the first place?

## What assets can refgenie build?

Refgenie now features a flexible asset class and recipe system, allowing users and tool developers to define, share, and manage their own asset types and build methods. Instead of being limited to a fixed set of internally-defined assets, you can now introduce new asset classes and recipes via external YAML or JSON files—without modifying refgenie's source code. This enables multiple assets of the same class to coexist, supports community-driven innovation, and decouples asset types from recipes for greater workflow flexibility.

**Data Channels**: Refgenie's extensibility is powered by **data channels**—external sources that provide definitions for asset classes and recipes. Data channels enable the community to publish, distribute, and synchronize new types of reference genome assets and build instructions across the refgenie ecosystem. Through data channels, organizations and tool developers can share their custom asset types with the broader community without requiring changes to the core refgenie codebase.

For a detailed overview of this system and how to define or use custom asset types and recipes, see the [Flexible Asset Types documentation](flexible_asset_types.md). To learn more about how data channels work and how to create or use them, see the [Data Channels documentation](data_channels.md)</h1>

## Recipes require inputs

Each asset requires some inputs, which we classify as _assets_, _files_ or _parameters_.

| **input class** |    **recipe name**   | **command line argument** |    **input type**   |
|:---------------:|:--------------------:|:-------------------------:|:-------------------:|
| assets          |   `required_assets`  |         `--assets`        | asset registry path |
| files           |   `required_files`   |         `--files`         |    file/dir path    |
| parameters      | `required_parametes` |       `--params`          |        string       |


> To view required inputs for an asset, add an `-q` or `--requirements` flag to the `refgenie build` command:

```console
$ refgenie build hg38/bowtie2_index -q

                                                                               Recipes
┏━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Name          ┃ Version ┃ Output asset class ┃ Input asset classes                ┃ Input files ┃ Input params                        ┃ Docker image               ┃
┡━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ bowtie2_index │ 0.0.1   │ bowtie2_index      │ • fasta (fasta asset for genome)   │ None        │ • threads (Number of threads to     │ docker.io/databio/refgenie │
│               │         │                    │ default=fasta                      │             │ use) default=1                      │                            │
└───────────────┴─────────┴────────────────────┴────────────────────────────────────┴─────────────┴─────────────────────────────────────┴────────────────────────────┘─────────────────────────────┴─────────────┴─────────────────────────────────────┴────────────────────────────┘
```

Notice how 'fasta' appears under `assets` and not under `files` or `params`. This means to build a bowtie2 index, you do *not* provide a fasta file as an *argument*; instead, you *must already have a fasta asset managed by `refgenie`*. One advantage of this is that it allows refgenie to keep a record of how you've built your assets, so refgenie can remember the link between this bowtie2 asset and the fasta asset, which turns out to be very useful for maintaining provenance of your assets. It also makes it easier to build derived assets like this, because you don't actually have to pass any additional arguments to build them.

So, you'll need to build the `fasta` asset for `hg38` genome before building `bowtie2_index`, but once you have that, building this asset is as simple as typing:

```console
refgenie build hg38/bowtie2_index
```

For many of the built-in recipes, a pre-existing `fasta` asset is the only requirement and refgenie will use the correct one by default. However, if you wish to build an asset with a different asset as an input, refgenie provides full flexibility. For instance, you can use `fasta:other_tag` (non-default tag) or even `hg38_cdna/fasta` (`fasta` asset from a different namespace) by adding `--assets` argument to the `refgenie build` command, like so:

```console
refgenie build hg38/bowtie2_index --assets fasta=hg38_cdna/fasta
```

This will build a `bowtie2_index` asset in `hg38` namespace but based on a transcriptome `fasta`.

Next, here's an example of an asset that requires an argument, but not a pre-existing asset:

```console
$ refgenie build hg38/refgene_anno -q

┏━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━┓
┃ Name         ┃ Version ┃ Output asset class ┃ Input asset classes ┃ Input files                                                  ┃ Input params ┃ Docker image     ┃
┡━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━┩
│ refgene_anno │ 0.0.1   │ refgene_anno       │ None                │ • refgene (gzipped RefGene database annotation file)         │ None         │ databio/refgenie │
│              │         │                    │                     │ default=None                                                 │              │                  │
└──────────────┴─────────┴────────────────────┴─────────────────────┴──────────────────────────────────────────────────────────────┴──────────────┴──────────────────┘
```

You'll need to provide this recipe with a `refgene` file, like this:

```console
refgenie build hg38/refgene_anno --files refgene=REFGENE_FILE.txt.gz
```

You can see the [example build output](build_output.md).

## Recipes require software

If you want to build assets, you'll need to get the software required by the asset you want to build. You have three choices to get that software: you can either install it natively, use a docker image, or use a bulker manifest.

### Install build software natively

`Refgenie` expects to find in your `PATH` any tools needed for building a desired asset. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie). We discourage this approach because it makes the assets dependent on your particular uncontrolled environment, which is not ideal. As a result, we don't have great documentation for what is required if you want to use this native approach. As we develop a custom asset system, we're planning to revamp this to provide more detailed way to see what requirements are for a specific recipe.

### Build assets with docker

If you don't want to install all the software needed to build all these assets (and I don't blame you), then you can just use `docker`. Each of our recipes knows about a `docker image` that has everything it needs. If you have `docker` installed, you should be able to simply run `refgenie build` with the `-d` flag. For example:

```console
refgenie build -d genome/asset ...
```

This tells `refgenie` to execute the building in a `docker container` requested by the particular asset recipe you specify. `Docker` will automatically pull the image it needs when you call this. If you like, you can build the `docker container` yourself like this:

```console
git clone https://github.com/databio/refgenie.git
cd refgenie/containers
make refgenie
```

or pull it directly from [dockerhub](https://hub.docker.com/r/databio/refgenie) like this:

```console
docker pull databio/refgenie
```

### Build assets with bulker

For an even more seamless integration of containers with `refgenie`, learn about [bulker](http://bulker.io), our multi-container environment manager. Here, you'd just need to do this:

```console
pip install bulker

# Next, configure bulker according to your local compute environment

bulker load databio/refgenie:0.7.0
bulker activate databio/refgenie:0.7.0
refgenie build ...
```

Bulker works on both singularity and docker systems. The bulker docs also contain a [more complete tutorial of using bulker and refgenie together](http://bulker.databio.org/en/latest/refgenie_tutorial/).

## Staging after building

Once you have built an asset, you can **stage** it to make it available for serving or pushing to cloud storage. You can do this in a separate step:

```console
refgenie build hg38/bowtie2_index
refgenie stage stage hg38/bowtie2_index
```

Or combine building and staging in one command with the `--stage` flag:

```console
refgenie build hg38/bowtie2_index --stage
```

You can also create push intent records at the same time using `--push-to`:

```console
refgenie build hg38/bowtie2_index --stage --push-to my-s3-remote
```

The `--push-to` flag accepts one or more remote names or IDs. After the build and stage complete, refgenie creates records indicating that these assets should be pushed to the specified remotes. You then run `refgenie push` to execute the actual upload.

See [Stage assets for serving](staging.md) for details on the staging system and [Configure remote storage](remotes.md) for setting up remotes.

## Versioning the assets

Refgenie supports tags to facilitate management of multiple "versions" of the same asset. Simply add a `:your_asset_name` appendix to the asset registry path in the `refgenie build` command and the created asset will be tagged:

```console
refgenie build hg38/bowtie2_index:my_asset_name
```
