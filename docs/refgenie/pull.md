# Download pre-built reference genome assets

## Introduction

Use the `refgenie` command-line interface to download and organize genome assets. You do this by simply running `refgenie` from the command line.

The `listr` command *lists remote assets* to see what's available:

```console
refgenie listr
```

The `pull` *downloads* the specific asset of your choice:

```console
refgenie pull GENOME/ASSET
```

Where `GENOME` refers to a genome key (*e.g.* hg38) and `ASSET` refers to one or more specific asset keys (*e.g.* bowtie2_index). For example:

```console
refgenie pull hg38/bowtie2_index
```

## What a pull sets up for you

A pull needs nothing set up beyond `refgenie init`. You do not need to add a data channel first.

- **Asset classes.** Every asset belongs to an asset class (such as `fasta` or `bwa_index`). If you do not have the asset's class yet, refgenie downloads its definition from the server you are pulling from and registers it. You only need a [data channel](data_channels.md) to *build* assets yourself.
- **Parent assets.** Some assets are made from others. A `bwa_index` is built from a `fasta`, so the `fasta` is its parent. If a parent is missing, refgenie pulls it first, from the same server, and then pulls the asset you asked for. The log says so:

    ```
    Pulling parent fasta:default (digest a0e94e...) required by bwa_index:0.7.19
    ```

    Parents of parents are pulled the same way. The size flags below (`--large`, `--skip-large`, `--size-cutoff`, `--batch`) apply to parent pulls too. A parent that was pulled stays, even if the asset you asked for then fails.

## Pull multiple assets

You can pull several assets in one command by listing them as positional arguments:

```console
refgenie pull hg38/bowtie2_index hg38/bwa_index mm10/fasta
```

## Batch pulling with `--genome`

To pull all available assets for a specific genome, use the `--genome` flag combined with `--all`:

```console
refgenie pull --genome hg38 --all
```

This downloads every asset available on the server for the `hg38` genome. You can also pull all assets across all genomes:

```console
refgenie pull --all --all-genomes
```

## Controlling large asset downloads

Some assets are very large. Refgenie provides flags to control how large assets are handled:

| Flag | Description |
|---|---|
| `--skip-large` | Skip assets larger than the size cutoff without prompting |
| `--large` | Download large assets without prompting |
| `--size-cutoff` | Size threshold in bytes (assets above this are "large") |
| `--batch` | Non-interactive mode: apply the large/skip-large policy without prompting |

For example, to pull all assets but skip anything over 1 GB:

```console
refgenie pull --genome hg38 --all --skip-large --size-cutoff 1073741824 --batch
```

## Initialize genome on pull

If you are pulling an asset for a genome that has not been initialized locally, use the `--init` flag to automatically initialize the genome:

```console
refgenie pull hg38/fasta --init
```

## Force re-download

To re-download an asset that already exists locally, use the `--force` flag:

```console
refgenie pull hg38/fasta --force
```

## Skipping specific asset classes or recipes

You can exclude certain asset classes or recipes from batch pulls:

```console
refgenie pull --genome hg38 --all --skip-asset-class star_index
refgenie pull --genome hg38 --all --skip-recipe star_index
```

## Error handling

If a pull fails (e.g., network error, asset not found on server), refgenie logs a warning and continues with the remaining assets. Check the log output for details about any failures.

Common issues:

- **"No local digest for genome alias"**: The genome alias is not yet known locally. Refgenie will attempt to resolve it with the server automatically, or you can use `--init`.
- **"No server subscriptions found"**: Your configuration subscribes to no server (it was created with `--no-default-server`, or you unsubscribed). Run `refgenie subscribe -s https://api.refgenie.org`, or use the **Servers** page of `refgenie dash`.
- **Server unreachable**: Verify your server subscription with `refgenie config get` and check your network connection.
- **Asset not found**: Use `refgenie listr` to verify the asset exists on the server.
- **"Could not pull parent ..."**: The asset needs a parent that you do not have, and the parent could not be pulled. The message names the parent (genome digest, asset group, asset name and digest) and says why its pull failed.
- **"Asset class '...' not found"**: You do not have the asset's class, and the server could not send it (older servers may not). Add and sync a data channel to get it: `refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml && refgenie data-channel sync refgenie --exists-ok`.

To see more details, consult the usage docs by running `refgenie pull --help`.

## Downloading manually

You can also browse and download pre-built `refgenie` assemblies manually at [refgenomes.databio.org](http://refgenomes.databio.org).
