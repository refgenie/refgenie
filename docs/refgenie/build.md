# Build assets

`refgenie pull` [downloads pre-built assets](pull.md). Use `refgenie build` for genomes or assets the server does not have. To build, you need the recipe's inputs and the software it runs.

!!! tip "See also"
    The [Building and serving tutorial](building_tutorial.md) walks through build, stage, serve, and push from start to finish.

## What can refgenie build?

Refgenie builds anything it has a recipe for. Recipes and asset classes come from [data channels](data_channels.md); a new install has none. Add the refgenie channel first:

```console
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok
refgenie recipe list
```

To write your own, see the [recipe specification](recipe_specification.md) and [asset class specification](asset_class_specification.md).

## Recipes need inputs

A recipe takes three kinds of input:

| Input | Recipe key | Command-line flag | Value |
|---|---|---|---|
| assets | `input_assets` | `--assets` | asset registry path |
| files | `input_files` | `--files` | file or folder path |
| parameters | `input_params` | `--params` | string |

To see what a recipe needs, add `-q` (`--requirements`):

```console
$ refgenie build hg38/bowtie2_index -q
```

The table shows that `bowtie2_index` needs a `fasta` *asset*, not a FASTA *file*. You must already have a `fasta` asset for `hg38`. Refgenie uses it automatically and records the link, so you can always tell which FASTA an index was built from:

```console
refgenie build hg38/bowtie2_index
```

To use a different input asset, such as another asset name or another genome's `fasta`, pass `--assets`:

```console
refgenie build hg38/bowtie2_index --assets fasta=hg38_cdna/fasta
```

Add `--pull-parents` to pull a missing input asset from your server before building.

Other recipes take files or parameters:

```console
refgenie build hg38/refgene_anno --files refgene=REFGENE_FILE.txt.gz
refgenie build hg38/bowtie2_index --params threads=8
```

## Recipes need software

Each recipe runs command-line tools. Get them one of three ways.

### Run in the recipe's container

Each recipe names a container image in `docker_image` that has every tool it needs. With Docker or Apptainer installed, add `--container`:

```console
refgenie build hg38/bowtie2_index --container
```

Refgenie runs each of the recipe's commands inside that image. It uses Docker when the Docker daemon answers, and Apptainer otherwise (most clusters have Apptainer). Set `REFGENIE_CONTAINER_ENGINE` to `docker` or `apptainer` to choose. The image is pulled the first time it is needed. On a cluster, point `REFGENIE_APPTAINER_CACHEDIR` (or Apptainer's own `APPTAINER_CACHEDIR`) at scratch space so the image cache does not fill your home folder. Use `--volumes` to make extra folders visible inside the container.

### Install the tools yourself

Without `--container`, every tool the recipe needs must be on your `PATH`; see the recipe's `docker_image` for what it uses (`refgenie recipe show <name>`).

### Use bulker

[Bulker](https://bulker.io) puts containerized tools on your `PATH` as if they were installed, so a build without `--container` can still run each tool in a container.

## Name the asset

Each build creates an asset inside an asset group. Recipes usually name the asset after the tool version, such as `2.5.4`. To choose the name yourself, add it after a colon:

```console
refgenie build hg38/bowtie2_index:my_asset_name
```

The asset you build becomes the group's default. See [Asset groups and naming](naming_assets.md).

## Stage after building

To serve an asset or push it to cloud storage, [stage](staging.md) it. Stage in a separate step:

```console
refgenie stage add hg38/bowtie2_index
```

or while building:

```console
refgenie build hg38/bowtie2_index --stage
```

Add `--push-to` with one or more remote names to mark the staged asset for upload; `refgenie push` then uploads it:

```console
refgenie build hg38/bowtie2_index --stage --push-to my-s3-remote
refgenie push
```

See [Configure remote storage](remotes.md).
