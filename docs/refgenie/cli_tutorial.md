# CLI Tutorial

This tutorial walks you through common refgenie operations, from installation to managing reference genome assets. By the end, you will be able to download, build, and retrieve paths to reference genome resources using the command line.

## Prerequisites

- Python 3.9 or higher
- pip package manager
- Internet connection (for downloading remote assets)

## Step 1: Installation

!!! warning "Transitionary Period"
    During the current transition period, install `refgenie1` and use the `refgenie1` command. Once the transition is complete, the package will be available as `refgenie` on PyPI.

Download the latest wheel from [GitHub releases](https://github.com/refgenie/refgenie1/releases) and install:

```bash
pip install refgenie1-*.whl
```

Verify the installation:

```bash
refgenie1 --version
```

## Step 2: Initialize Refgenie

Before using refgenie, initialize the configuration and database:

```bash
refgenie1 init
```

This creates a configuration directory at `~/.refgenie/` with:

- A SQLite database for tracking assets
- A `genomes/` folder for storing downloaded assets
- A database configuration file

!!! tip "Custom Location"
    Set the `REFGENIE_HOME_PATH` environment variable to use a different base directory:
    ```bash
    export REFGENIE_HOME_PATH=~/my_genomes
    refgenie1 init
    ```

## Step 3: Check Your Configuration

View the current refgenie configuration:

```bash
refgenie1 config get
```

Expected output:

```
                                               Refgenie configuration
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ genome_folder                       ┃ version ┃ genome_archive_folder                ┃ servers ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ /home/user/.refgenie/genomes        │ 1       │ /home/user/.refgenie/archives        │ []      │
└─────────────────────────────────────┴─────────┴──────────────────────────────────────┴─────────┘

Environment-based configuration

log_level: LogLevel.INFO
genome_folder: /home/user/.refgenie/genomes
genome_archive_folder: /home/user/.refgenie/archives
database_config_path: /home/user/.refgenie/refgenie_db_config.yaml
```

## Step 4: Subscribe to a Server

To download pre-built assets, subscribe to a remote refgenie server:

```bash
refgenie1 subscribe http://refgenomes.databio.org
```

This adds the public refgenie server to your configuration. You can subscribe to multiple servers.

!!! note "What is a refgenie server?"
    A refgenie server hosts pre-built reference genome assets. The public server at `refgenomes.databio.org` provides common genomes like hg38, mm10, and more. Organizations can also run private servers.

To unsubscribe from a server:

```bash
refgenie1 unsubscribe -s http://refgenomes.databio.org
```

## Step 5: List Remote Assets

Browse what assets are available on subscribed servers:

```bash
refgenie1 listr
```

Expected output:

```
                 Refgenie assets. Source: http://refgenomes.databio.org
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Genome digest                                    ┃ Asset group             ┃ Asset   ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 │ fasta                   │ default │
│ 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 │ bowtie2_index           │ default │
│ 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 │ bwa_index               │ default │
│ ...                                              │ ...                     │ ...     │
└──────────────────────────────────────────────────┴─────────────────────────┴─────────┘
```

!!! tip "Filtering Results"
    Filter by genome to see only relevant assets:
    ```bash
    refgenie1 listr -g hg38
    ```

## Step 6: Pull an Asset

Download a pre-built asset from the server:

```bash
refgenie1 pull hg38/fasta
```

This downloads the human reference genome (hg38) FASTA file and related indexes. Refgenie automatically:

1. Resolves the `hg38` alias to its unique genome digest
2. Downloads the asset archive
3. Extracts and organizes the files
4. Registers the asset in the local database

Expected output:

```
INFO     Setting 'hg38' identity with server: http://refgenomes.databio.org
INFO     Determined digest for hg38: 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4
INFO     Set genome alias: hg38
hg38/fasta:default ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% • 833.4/833.4 MB
```

!!! tip "Pull Multiple Assets"
    Download several assets at once:
    ```bash
    refgenie1 pull --genome hg38 fasta bowtie2_index bwa_index
    ```

## Step 7: List Local Assets

View assets you have downloaded or built locally:

```bash
refgenie1 list
```

Expected output:

```
                                      Local refgenie assets
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━┓
┃ Genome                                           ┃ Asset         ┃ Tag      ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━┩
│ hg38                                             │ fasta         │ default  │
└──────────────────────────────────────────────────┴───────────────┴──────────┘
```

## Step 8: Seek Asset Paths

Retrieve the local path to a downloaded asset using `seek`:

```bash
refgenie1 seek hg38/fasta
```

Expected output:

```
/home/user/.refgenie/genomes/alias/hg38/fasta/default
```

This returns the absolute path to the asset folder, making your scripts portable across different computing environments.

### Using Seek Keys

Some assets contain multiple files. The `fasta` asset, for example, includes a FASTA file, an index (.fai), and chromosome sizes. Access specific files with seek keys:

```bash
# Get the FASTA file path
refgenie1 seek hg38/fasta

# Get the FASTA index (.fai) path
refgenie1 seek hg38/fasta.fai

# Get the chromosome sizes file path
refgenie1 seek hg38/fasta.chrom_sizes
```

!!! tip "Portable Scripts"
    Instead of hardcoding paths in your scripts, use refgenie seek:
    ```bash
    # Before (not portable)
    bwa mem /path/to/genomes/hg38/bwa_index reads.fq

    # After (portable)
    bwa mem $(refgenie1 seek hg38/bwa_index) reads.fq
    ```

### Remote Seek

Get paths to assets on remote servers without downloading them:

```bash
refgenie1 seekr hg38/fasta
```

This is useful for cloud workflows where you want to stream data directly from S3 or other remote storage.

## Step 9: Set Genome Aliases

Refgenie uses sequence-derived digests to uniquely identify genomes. Aliases like `hg38` or `mm10` are human-friendly names that map to these digests.

### View Existing Aliases

```bash
refgenie1 alias get
```

### Set a New Alias

If you know the genome digest, you can set an alias manually:

```bash
refgenie1 alias set --aliases hg38 --digest 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4
```

You can also set multiple aliases for the same genome:

```bash
refgenie1 alias set --aliases GRCh38 human_reference --digest 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4
```

!!! note "Automatic Alias Resolution"
    When you pull an asset using a name like `hg38`, refgenie automatically queries the server to resolve the alias to its digest. This ensures you always get the correct genome, even if different servers use different naming conventions.

### Remove an Alias

```bash
refgenie1 alias remove --aliases my_old_alias
```

## Step 10: Build a Custom Asset

For genomes or assets not available on remote servers, you can build them locally. Here is a brief example of building a FASTA asset:

```bash
refgenie1 build my_genome/fasta --files fasta=/path/to/my_genome.fa.gz --genome-description "My custom genome"
```

Once you have a FASTA asset, you can build derived assets like indexes:

```bash
refgenie1 build my_genome/bowtie2_index
```

!!! info "Full build/stage/push workflow"
    For the complete workflow — including how to stage assets for serving, serve them locally, build in bulk with Snakemake, and push to cloud remotes — see the [Building and Serving Tutorial](building_tutorial.md).

    For reference documentation on building, see [Build](build.md).

## Quick Reference

| Command | Description |
|---------|-------------|
| `refgenie1 init` | Initialize configuration |
| `refgenie1 config get` | View current configuration |
| `refgenie1 subscribe URL` | Add a remote server |
| `refgenie1 unsubscribe -s URL` | Remove a remote server |
| `refgenie1 listr` | List remote assets |
| `refgenie1 listr -g GENOME` | List remote assets for a genome |
| `refgenie1 pull GENOME/ASSET` | Download an asset |
| `refgenie1 list` | List local assets |
| `refgenie1 seek GENOME/ASSET` | Get local asset path |
| `refgenie1 seekr GENOME/ASSET` | Get remote asset path |
| `refgenie1 alias get` | View genome aliases |
| `refgenie1 alias set -a ALIAS -d DIGEST` | Set a genome alias |
| `refgenie1 build GENOME/ASSET` | Build an asset |

## Next Steps

- [Configuration](configuration.md) - Advanced configuration options, PostgreSQL setup
- [Build](build.md) - Building custom assets and using Docker
- [Data Channels](data_channels.md) - Adding custom asset types and recipes
- [Seek](seek.md) - More on retrieving asset paths
- [Pull](pull.md) - More on downloading assets
