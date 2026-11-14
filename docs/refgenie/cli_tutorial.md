# CLI tutorial

This tutorial covers common refgenie tasks on the command line: downloading assets, getting paths to them, managing aliases, and building your own. The examples use `rCRSd`, the small human mitochondrial genome, so downloads take seconds.

## Prerequisites

- Python 3.11 or higher
- An internet connection

## Step 1: Install

```bash
pip install "refgenie>=1.0.0a2"
refgenie --version
```

## Step 2: Initialize

```bash
refgenie init
```

This creates `~/.refgenie/` with a SQLite database, a `genomes/` folder for assets, and a database config file. It also subscribes you to the public server at `https://api.refgenie.org`.

!!! tip "Custom location"
    Set `REFGENIE_HOME_PATH` to use a different base directory:
    ```bash
    export REFGENIE_HOME_PATH=~/my_genomes
    refgenie init
    ```

## Step 3: Check your configuration

```bash
refgenie config get
```

This shows your genome folder, stage folder, database config path, and subscribed servers. See [Configure](configuration.md) for details.

## Step 4: Manage server subscriptions

A refgenie server hosts pre-built assets. You are already subscribed to `https://api.refgenie.org`. To add another server, such as one your lab runs:

```bash
refgenie subscribe -s https://my-lab-server.example.org
```

To remove it:

```bash
refgenie unsubscribe -s https://my-lab-server.example.org
```

## Step 5: List remote assets

```bash
refgenie listr
```

```
                         Refgenie assets. Source: https://api.refgenie.org
┏━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Aliases          ┃ Genome digest                    ┃ Asset group        ┃ Asset   ┃
┡━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ rCRSd            │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ bowtie2_index      │ 2.5.5   │
│ rCRSd            │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ bwa_index          │ 0.7.19  │
│ rCRSd            │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ fasta              │ default │
...
```

Each row is one asset. The *asset group* is the kind of asset (`bowtie2_index`), and the *asset* is a named build within the group (`2.5.5`). To list the assets for one genome, name it: `refgenie listr -g rCRSd`. It need not be on your machine yet.

## Step 6: Pull an asset

```bash
refgenie pull rCRSd/fasta
```

Refgenie looks up the `rCRSd` alias on the server to find its genome digest, downloads the asset, and registers it in your local database. Pull several at once by listing more registry paths:

```bash
refgenie pull rCRSd/bowtie2_index rCRSd/bwa_index
```

## Step 7: List local assets

```bash
refgenie list
```

```
                    Refgenie assets. Source: local
┏━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Aliases ┃ Genome digest                    ┃ Asset group   ┃ Asset   ┃
┡━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ rCRSd   │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ bowtie2_index │ 2.5.5   │
│ rCRSd   │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ bwa_index     │ 0.7.19  │
│ rCRSd   │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ fasta         │ default │
└─────────┴──────────────────────────────────┴───────────────┴─────────┘
```

## Step 8: Get asset paths

`seek` returns the local path to an asset:

```bash
refgenie seek rCRSd/fasta
```

```
/home/user/.refgenie/genomes/alias/rCRSd/fasta/default/rCRSd.fa
```

Some assets contain several files. The `fasta` asset holds the FASTA file, its index, and chromosome sizes. Pick a file with a *seek key*:

```bash
refgenie seek rCRSd/fasta.fai          # FASTA index
refgenie seek rCRSd/fasta.chrom_sizes  # chromosome sizes
```

For aligner indexes, `seek` returns the index prefix that the aligner expects:

```bash
bowtie2 -x $(refgenie seek rCRSd/bowtie2_index) -U reads.fq
```

Use the bowtie2 version that built the index: `refgenie seek rCRSd/bowtie2_index.version` prints it (`2.5.5`). Install that version with conda (`conda install -c bioconda bowtie2=2.5.5`), or, without conda, download the prebuilt binary from the [bowtie2 GitHub releases](https://github.com/BenLangmead/bowtie2/releases).

Using `seek` in scripts, instead of fixed paths, makes them work on any computer.

### Remote paths

`seekr` returns the URL of an asset file on a server, without downloading it:

```bash
refgenie seekr rCRSd/fasta
```

## Step 9: Manage genome aliases

Refgenie identifies each genome by a digest computed from its sequences. Aliases like `rCRSd` or `hg38` are names that point to these digests. When you pull by alias, refgenie gets the matching digest from the server.

```bash
# List aliases
refgenie alias get

# Add aliases to a genome (comma-separate several)
refgenie alias set -a mito,human_mt -d jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP

# Remove an alias
refgenie alias remove -a human_mt
```

See [Manage aliases](alias.md).

## Step 10: Build a custom asset

For genomes or assets that no server has, build them yourself. Building needs recipes, which come from a [data channel](data_channels.md):

```bash
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok
```

Register your genome from its FASTA file. This also builds its `fasta` asset:

```bash
refgenie genome init --fasta /path/to/my_genome.fa.gz --name my_genome --description "My custom genome"
```

Then build other assets from it (this one needs `bowtie2` installed):

```bash
refgenie build my_genome/bowtie2_index
```

See [Build assets](build.md), and the [Building and serving tutorial](building_tutorial.md) for staging, serving, and pushing assets to the cloud.

## Quick reference

| Command | Description |
|---------|-------------|
| `refgenie init` | Initialize configuration |
| `refgenie config get` | View current configuration |
| `refgenie subscribe -s URL` | Add a server |
| `refgenie unsubscribe -s URL` | Remove a server |
| `refgenie listr` | List remote assets |
| `refgenie pull GENOME/ASSET` | Download an asset |
| `refgenie list` | List local assets |
| `refgenie seek GENOME/ASSET` | Get local asset path |
| `refgenie seekr GENOME/ASSET` | Get remote asset URL |
| `refgenie alias get` | View genome aliases |
| `refgenie alias set -a ALIAS -d DIGEST` | Set a genome alias |
| `refgenie genome init --fasta FILE --name NAME` | Register a genome from a FASTA file |
| `refgenie build GENOME/ASSET` | Build an asset |

## Next steps

- [Configure](configuration.md): other locations and PostgreSQL
- [Build assets](build.md): building custom assets and using containers
- [Data channels](data_channels.md): adding asset types and recipes
- [Retrieve paths to assets](seek.md): more on `seek`
- [Download pre-built assets](pull.md): more on `pull`
