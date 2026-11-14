# Download pre-built reference genome assets

`refgenie listr` lists the assets on the servers you subscribe to. A new configuration subscribes to [api.refgenie.org](https://api.refgenie.org).

```console
refgenie listr
```

`refgenie pull` downloads one or more assets, named by [registry path](asset_registry_paths.md):

```console
refgenie pull hg38/bowtie2_index
refgenie pull hg38/bowtie2_index hg38/bwa_index mm10/fasta
```

## What a pull sets up for you

A pull needs nothing set up beyond `refgenie init`.

- **Asset classes.** If you do not have the asset's class (such as `fasta` or `bwa_index`), refgenie downloads it from the server and registers it. You need a [data channel](data_channels.md) only to *build* assets yourself.
- **Parent assets.** Some assets are made from others: a `bowtie2_index` is built from a `fasta`. A pull brings only the asset you name, and that asset works on its own. Pull the parent as well if you need it (`refgenie pull hg38/fasta`). The exception is an asset that keeps links to its parent's files in its own folder, such as `bismark_bt2_index`. Refgenie pulls that parent first, from the same server, and creates the links.

## Pull many assets at once

Pull every asset of one or more genomes:

```console
refgenie pull -g hg38 --all
refgenie pull -g hg38,mm10 --all
```

Pull one asset group across several genomes, or across every genome on your servers:

```console
refgenie pull -g hg38,mm10 --asset fasta
refgenie pull --all-genomes --asset fasta
```

To copy everything on your servers, use `refgenie mirror`.

## Large assets

Refgenie asks before it downloads an archive bigger than the size cutoff (10 GB by default).

| Flag | What it does |
|---|---|
| `--size-cutoff` | The cutoff, in GB |
| `--skip-large` | Skip archives over the cutoff without asking |
| `--pull-large` | Download archives over the cutoff without asking |
| `--batch` | Never ask; download everything (same as `--pull-large`) |

For example, pull all of hg38 but skip anything over 1 GB:

```console
refgenie pull -g hg38 --all --skip-large --size-cutoff 1
```

## Other options

- `--init` registers genomes (aliases and metadata) without downloading any asset files. Name them with `-g` or `--genome-digest`: `refgenie pull -g hg38 --init`.
- `--force` (`-f`) replaces an asset you already have and skips confirmation prompts.

Run `refgenie pull --help` for the full list.

## Errors

When one asset in a list fails, refgenie logs the error and goes on with the rest.

- **"No server subscriptions found"**: you subscribe to no server (the configuration was made with `--no-default-server`, or you unsubscribed). Run `refgenie subscribe -s https://api.refgenie.org`, or use the **Servers** page of `refgenie dash`.
- **Server unreachable**: check your subscriptions with `refgenie config get`, and your network connection.
- **Asset not found**: run `refgenie listr` to see what the server has.
- **"Could not pull parent ..."**: the asset needs a parent you do not have, and the parent could not be pulled. The message names the parent and says why its pull failed.
- **"Asset class '...' not found"**: the server could not send the asset's class. Add and sync a data channel to get it:

    ```console
    refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
    refgenie data-channel sync refgenie --exists-ok
    ```
