# Refgenie CLI

The installed command is `refgenie`. This document walks through a first session
and then documents every command and flag.

## Installation

Install from [PyPI](https://pypi.org/project/refgenie/). Because 1.0 is still a
pre-release, `pip` must be told to accept it:

```bash
pip install --pre refgenie
```

The base install provides the CLI and the Python API. `refgenie dash`,
`refgenie serve`, `refgenie-mcp`, and snakemake-based bulk building each require
an extra (`dash`, `server`, `mcp`, `snakemake`); see the
[README installation section](../README.md#installation) for the extras table.
Running a command without its extra prints a message naming the extra to install.

Check the installed version:

```bash
$ refgenie --version
refgenie 1.0.0a1
```

## Getting started

Refgenie stores asset metadata in a database — SQLite by default, PostgreSQL
optionally. With no configuration, refgenie creates a `refgenie` SQLite file
under `~/.refgenie` (or `$REFGENIE_HOME_PATH`). Configuration is initialized
automatically on first use, but running `refgenie init` explicitly is the
supported starting point: it also creates the genome folder and the stage folder.

```bash
$ refgenie init
INFO     Database configuration file created at /home/user/.refgenie/refgenie_db_config.yaml.
INFO     Genome folder ready: /home/user/.refgenie/genomes
INFO     Genome stage folder ready: /home/user/.refgenie/archives
INFO     Initialized refgenie backend: 'sqlite:////home/user/.refgenie/refgenie'
```

No asset classes or recipes — not even `fasta` — are registered by the package.
They come from a data channel, which you must register and sync before building
anything.

### Order of operations

1. `refgenie init` — initialize the config, database, and folders.
2. `refgenie data-channel add ...` then `refgenie data-channel sync ...` —
   register asset classes and recipes.
3. `refgenie genome init ...` — register a genome. This auto-builds the `fasta`
   asset, but only once a `fasta` recipe is registered by step 2. Running
   `genome init` before syncing a data channel skips the auto-build and prints a
   message; `--no-build` skips the attempt entirely.
4. `refgenie build ...` — build other assets (e.g. `bwa_index`) whose recipes
   are registered.

### Register a data channel

A data channel is an index of asset class and recipe definitions. The canonical
channel is published at
[refgenie-registry](https://refgenie.github.io/refgenie-registry).

```bash
$ refgenie data-channel add my-fav-channel https https://refgenie.github.io/refgenie-registry/index.yaml
INFO     Added data channel: my-fav-channel
```

The three positional arguments are the channel name, its type, and the address
of its `index.yaml`.

List registered channels:

```bash
$ refgenie data-channel list
                                 Data Channels
┏━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┓
┃ Name           ┃ Type  ┃ Index Address               ┃ Description ┃ Credentials set ┃
┡━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━┩
│ my-fav-channel │ https │ https://refgenie.github.io… │             │ False           │
└────────────────┴───────┴─────────────────────────────┴─────────────┴─────────────────┘
```

### Sync asset classes and recipes

```bash
$ refgenie data-channel sync my-fav-channel --exists-ok
INFO     Registered 'fasta' recipe
...
INFO     Successfully synced from channel 'my-fav-channel'
```

`--exists-ok` skips items already present instead of erroring. Use
`--exists-overwrite` to replace conflicting items.

Recipes contain shell commands, and building an asset runs them. Sync only
channels you trust.

Inspect what arrived:

```bash
refgenie recipe list
refgenie asset-class list
```

### Register a genome and build its fasta asset

`genome init` computes the genome's sequence-collection digest and, when a
`fasta` recipe is available, builds the `fasta` asset in the same step.

```bash
$ refgenie genome init --fasta rCRSd.fa --name rCRSd --species 'Homo sapiens' \
    --description 'human mitochondrial genome'
INFO     Asset 'rCRSd/fasta:default' build succeeded
INFO     Added: 'rCRSd/fasta:default'
INFO     Fasta asset built successfully for rCRSd
```

The genome is now listed, and its assets have resolvable paths:

```bash
$ refgenie list
                    Refgenie assets. Source: local
┏━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Aliases ┃ Genome digest                    ┃ Asset group ┃ Asset   ┃
┡━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━┩
│ rCRSd   │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ fasta       │ default │
└─────────┴──────────────────────────────────┴─────────────┴─────────┘

$ refgenie seek rCRSd/fasta
/home/user/.refgenie/genomes/alias/rCRSd/fasta/default/rCRSd.fa
```

To build the `fasta` asset separately — for instance after a `genome init
--no-build` — call `build` with the source file:

```bash
refgenie build rCRSd/fasta --files fasta=rCRSd.fa
```

### Build another asset

Once the genome exists and the recipe is registered, other assets build from it:

```bash
refgenie build rCRSd/bwa_index
```

Recipes run real tools. If `bwa` is not on `PATH`, either install it or pass
`-d/--docker` to run the recipe in its container. `-q/--requirements` prints
what a recipe needs without building:

```bash
refgenie build rCRSd/bwa_index -q
```

### Retrieve sequence

```bash
$ refgenie getseq -g rCRSd -l 'rCRSd:0-30'
GATCACAGGTCTATCACCCTATTAACCACT
```

## Registry paths

Most asset commands take one or more *registry paths*:

```
<genome>/<asset>
<genome>/<asset>:<tag>
<genome>/<asset>.<seek_key>:<tag>
```

`<genome>` is an alias (e.g. `hg38`) or a sequence-collection digest. Omitting
`:<tag>` uses the default tag. The `.seek_key` suffix selects a specific file
within an asset (e.g. `hg38/fasta.fai`).

## Command reference

`refgenie --help` groups commands as follows.

| Group | Commands |
| --- | --- |
| Asset management | `list`, `asset`, `seek`, `add`, `remove`, `rename`, `id`, `build`, `populate` |
| Remote operations | `listr`, `seekr`, `pull`, `push`, `mirror`, `populater`, `compare` |
| Genome management | `genome` |
| Server | `serve`, `dash`, `subscribe`, `unsubscribe`, `catalog-export` |
| Configuration | `init`, `purge`, `config`, `alias`, `recipe`, `asset-class`, `stage` |
| Asset definitions | `data-channel`, `generate`, `remote` |
| Sequences | `getseq` |

Boolean flags shown as `-f, --force` also accept the negated form
(`--no-force`). Command names use hyphens (`asset-class`, `data-channel`,
`catalog-export`); `data_channel` is additionally accepted as an underscore
alias for `data-channel`.

### Asset management

#### `list`

List available local assets.

```
refgenie list [-g GENOME ...]
```

| Option | Description |
| --- | --- |
| `-g, --genome` | One or more genomes to restrict the listing to. |

#### `asset list`

Alias group. `refgenie asset list` is equivalent to `refgenie list` and takes
the same `-g, --genome` option.

#### `seek`

Print the local path of an asset.

```
refgenie seek ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `-e, --check-exists` | Check the returned path for existence on disk. |
| `--abs` | Return the digest-addressed content path under `data/` instead of the human-readable alias path. |

#### `add`

Register an asset that already exists on disk.

```
refgenie add ASSET-REGISTRY-PATHS [...] -p PATH -c ASSET_CLASS
```

| Option | Description |
| --- | --- |
| `-p, --path` | Relative local path to the asset. **Required.** |
| `-c, --asset-class` | Name of the asset's asset class. **Required.** |
| `-d, --description` | Description of the asset. |
| `-k, --seek-keys` | Non-path seek key values, as `name=value`. Repeat for multiple keys. |

#### `remove`

Remove a local asset.

```
refgenie remove ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `-f, --force` | Do not prompt before removing. |
| `-a, --aliases` | Also remove the genome alias if this was the genome's last asset. |

#### `rename`

Rename an asset.

```
refgenie rename ASSET-REGISTRY-PATHS [...] -n NEW_ASSET_NAME
```

| Option | Description |
| --- | --- |
| `-n, --new-asset-name` | New name for the asset. **Required.** |

#### `id`

Return a digest. A genome name yields the genome digest; a registry path yields
the asset digest.

```
refgenie id REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `-v, --verbose` | Show detailed genome metadata (sequence count, total length, source). |
| `--validate-store` | Verify that the genome's RefgetStore exists and is valid. |
| `--remote` | Query subscribed seqcolapi servers when the genome is not found locally. |
| `--info` | Given a digest, show its aliases and metadata. |

#### `build`

Build genome assets.

```
refgenie build ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `--asset-description` | Asset-level description (e.g. `built with version 0.3.2`). |
| `--recipe-name` | Recipe to use. |
| `--recipe-version` | Recipe version to use. |
| `-d, --docker` | Run all commands in the refgenie docker container. |
| `--pull-parents` | Automatically pull a required parent asset that was not provided. |
| `-q, --requirements` | Show the build requirements for the asset and exit. |
| `--stage` | Stage the asset after building. Requires the genome stage folder to be set. |
| `--push-to` | Remote names/IDs to create push intent records for after staging. |
| `--pipeline-kwargs` | Extra arguments for the build pipeline, as `arg_name=arg_val`. |
| `--assets` | Override the genome, asset, and tag of parents, e.g. `fasta=hg38/fasta:default`. |
| `--files` | Paths to required input files, e.g. `fasta=/path/to/file.fa.gz`. |
| `--params` | Required parameter values, e.g. `param1=value1`. |
| `--volumes` | Additional folders to mount as volumes when using docker. |

#### `populate`

Replace refgenie registry paths with local paths. Reads the file given by `-f`,
or stdin when `-f` is omitted.

```
refgenie populate [-f FILE]
```

| Option | Description |
| --- | --- |
| `-f, --file` | File containing registry paths to populate. |

### Remote operations

#### `listr`

List assets available on subscribed servers.

```
refgenie listr [-g GENOME ...] [-s URL ...]
```

| Option | Description |
| --- | --- |
| `-g, --genome` | One or more genomes to restrict the listing to. |
| `-s, --genome-server` | One or more server URLs to use for this call only; not persisted to config. |
| `-p, --append-server` | Append the provided servers to the configured list rather than replacing it. |

#### `seekr`

Print the remote path of an asset.

```
refgenie seekr ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `-s, --genome-server` | One or more server URLs to use for this call only; not persisted. |
| `-p, --append-server` | Append the provided servers to the configured list. |

#### `pull`

Download assets from subscribed servers. Registry paths may be given
positionally or with `--asset-registry-paths`.

```
refgenie pull ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `--asset-registry-paths` | Registry paths to pull (equivalent to the positional form). |
| `-g, --genome` | Reference assembly ID, e.g. `mm10`. |
| `--all` | Pull all assets for the specified genome(s). |
| `--all-genomes` | Apply the operation to all genomes available on subscribed servers. |
| `--asset` | Pull one asset type across the specified genomes, e.g. `--asset fasta`. |
| `--init` | Register genome(s) locally (aliases, metadata) without downloading asset files. |
| `--skip-large` | Do not pull archives over the size cutoff. |
| `--pull-large` | Pull all archives regardless of size. |
| `--size-cutoff` | Maximum archive size, in GB, to pull without confirmation. Default `10`. |
| `--batch` | Batch mode: pull all archives regardless of size. |
| `-f, --force` | Skip confirmation prompts for multi-asset operations. |

#### `push`

Upload staged assets to cloud remotes.

```
refgenie push [-r REMOTE] [-g GENOME]
```

| Option | Description |
| --- | --- |
| `-r, --remote` | Push only to this remote (by name or id). Default: all remotes with unpushed assets. |
| `-g, --genome` | Push only assets for this genome. |
| `-n, --dry-run` | Show what would be pushed without executing. |
| `--strategy` | `per_asset` (upload each asset) or `folder_sync` (sync the whole genome stage folder). Default `per_asset`. |

#### `mirror`

Mirror all assets from all genomes on subscribed servers.

```
refgenie mirror
```

| Option | Description |
| --- | --- |
| `--skip-large` | Do not pull archives over the size cutoff. |
| `--pull-large` | Pull all archives regardless of size. |
| `--size-cutoff` | Maximum archive size, in GB, to pull without confirmation. Default `10`. |
| `--batch` | Batch mode: pull all archives regardless of size. |
| `-f, --force` | Skip the confirmation prompt. |

#### `populater`

Replace refgenie registry paths with remote paths. Reads the file given by `-f`,
or stdin when `-f` is omitted.

```
refgenie populater [-f FILE]
```

| Option | Description |
| --- | --- |
| `-f, --file` | File containing registry paths to populate. |
| `-s, --genome-server` | One or more server URLs to use for this call only; not persisted. |
| `-p, --append-server` | Append the provided servers to the configured list. |

#### `compare`

Compare two genomes for compatibility.

```
refgenie compare GENOME1 GENOME2
```

### Genome management

#### `genome init`

Initialize a genome from a FASTA file, a refgenie server, or a RefgetStore.
When initialized from a FASTA file it also builds the `fasta` asset (`fa`,
`fai`, `chrom.sizes`).

```
refgenie genome init -n NAME [--fasta PATH | --server URL | --store URL]
```

| Option | Description |
| --- | --- |
| `-n, --name` | One or more alias names for the genome. **Required.** |
| `--fasta` | Path to a local FASTA file. |
| `--server` | URL of a refgenie server to initialize from. |
| `--store` | URL of a RefgetStore to initialize from; requires `--digest` or `--namespace`. |
| `--namespace` | Namespace for alias lookup when using `--store`. |
| `--digest` | Seqcol digest of the genome. |
| `-d, --description` | Genome description, e.g. `Human genome build 38`. |
| `-s, --species` | Species name, e.g. `Homo sapiens`. |
| `--fhr` | Path to an FHR `.fhr.json` metadata file to apply after init. When given it is authoritative for description/species and writes the RefgetStore sidecar. |
| `-f, --force` | Allow re-initialization of an existing genome (adds new aliases). |
| `--build` / `--no-build` | Build the fasta asset after initialization. Default on; runs only when a `fasta` recipe is registered. |

#### `genome set-metadata`

Apply FHR metadata to an already-registered genome, with no rebuild. Updates the
genome's description and species and the RefgetStore sidecar.

```
refgenie genome set-metadata (-n NAME | --digest DIGEST) --fhr PATH
```

| Option | Description |
| --- | --- |
| `-n, --name` | Genome alias to update. |
| `--digest` | Genome seqcol digest to update (alternative to `--name`). |
| `--fhr` | Path to the FHR `.fhr.json` metadata file to apply. **Required.** |

#### `genome list`

List all genomes, with digests, aliases, source, species, and description.

```
refgenie genome list
```

#### `genome remove`

Remove a genome and all its assets.

```
refgenie genome remove --genome NAME [...]
```

| Option | Description |
| --- | --- |
| `--genome` | Genome digest(s) or alias(es) to remove. **Required.** |
| `-f, --force` | Do not prompt before removing. |

#### `genome browse`

Browse genomes available on a refgenie server or RefgetStore.

```
refgenie genome browse [--server-url URL]
```

| Option | Description |
| --- | --- |
| `--server-url` | URL of a refgenie server or RefgetStore. Defaults to subscribed server(s). |
| `--page` | Page number for paginated results. Default `0`. |
| `--page-size` | Number of results per page. Default `20`. |

#### `genome sync`

Bulk-register all genomes from subscribed servers or a remote source.

```
refgenie genome sync [--server-url URL]
```

| Option | Description |
| --- | --- |
| `--server-url` | URL of a remote source to sync from. Defaults to all subscribed server(s). |
| `--page-size` | Number of collections to request per page. Default `1000`. |

### Server

#### `serve`

Start the production refgenie server. Requires the `server` extra.

```
refgenie serve [-p PORT]
```

| Option | Description |
| --- | --- |
| `-p, --port` | Port to run the server on. Default `8000`. |
| `-r, --reload` | Enable auto-reload on code changes (for development). |

#### `dash`

Start the local refgenie web UI. Requires the `dash` extra.

```
refgenie dash [-p PORT] [-b {off,read,full}]
```

| Option | Description |
| --- | --- |
| `-p, --port` | Port to run the dashboard on. Default `8080`. |
| `-b, --bridge` | Localhost-bridge mode for this run, overriding `$REFGENIE_BRIDGE_MODE`: `off` = no cross-origin access, `read` = allowlisted public origins may read, `full` = additionally allows cross-origin pull. |

#### `subscribe`

Add refgenieserver URLs to the config.

```
refgenie subscribe -s URL [...]
```

| Option | Description |
| --- | --- |
| `-s, --genome-server` | One or more URLs to add to the subscription list. |
| `-r, --reset` | Overwrite the current list of server URLs. |

#### `unsubscribe`

Remove refgenieserver URLs from the config.

```
refgenie unsubscribe -s URL [...]
```

| Option | Description |
| --- | --- |
| `-s, --genome-server` | One or more URLs to remove from the subscription list. |

#### `catalog-export`

Export a publish catalog covering pushed assets only, for a server to import.

```
refgenie catalog-export [--dest PATH] [--https-prefix URL]
```

| Option | Description |
| --- | --- |
| `--dest` | Path to write the publish-catalog SQLite artifact to. |
| `--https-prefix` | Public https base URL mirroring the stage folder that `refgenie push` uploaded to, e.g. `https://<bucket>.s3.amazonaws.com/assets`. Download links are served from here. |

### Configuration

#### `init`

Initialize the refgenie configuration, database, and folders.

```
refgenie init
```

| Option | Description |
| --- | --- |
| `-f, --genome-folder` | Absolute path to the parent folder for refgenie-managed assets. |
| `-a, --genome-stage-folder` | Absolute path to the parent stage folder for refgenie-managed assets; used by refgenieserver. |
| `-v, --config-version` | Config version to initialize the config file with. |

#### `purge`

Purge the genome configuration.

```
refgenie purge
```

| Option | Description |
| --- | --- |
| `-f, --force` | Do not prompt before purging. |

#### `config get` / `config set`

```
refgenie config get
```

`config get` displays the current configuration, including the database
connection and the environment-derived settings. `config set` is not yet
implemented.

#### `alias get` / `alias set` / `alias remove`

```
refgenie alias get [-a ALIAS ...] [-g DIGEST ...]
refgenie alias set -a ALIAS [...] [-d DIGEST]
refgenie alias remove -a ALIAS [...]
```

`alias get` options (mutually exclusive):

| Option | Description |
| --- | --- |
| `-a, --aliases` | Aliases to get the digests for. |
| `-g, --genome-digests` | Genome digests to get the aliases for. |

`alias set` options:

| Option | Description |
| --- | --- |
| `-a, --aliases` | Aliases to set. **Required.** |
| `-d, --digest` | Digest to set the aliases on. |
| `-r, --reset` | Remove all aliases before setting the new ones. |
| `-f, --force` | Force the action even if the genome does not exist. |

`alias remove` options:

| Option | Description |
| --- | --- |
| `-a, --aliases` | Aliases to remove. **Required.** |

#### `recipe`

```
refgenie recipe list
refgenie recipe show RECIPE-NAME [--recipe-version VERSION]
refgenie recipe requirements RECIPE-NAME [--recipe-version VERSION]
refgenie recipe add --source PATH_OR_URL [-f]
refgenie recipe remove RECIPE-NAME [--recipe-version VERSION]
```

| Subcommand | Description |
| --- | --- |
| `list` | List local recipes. |
| `show` | Display a recipe. |
| `requirements` | Show a recipe's requirements. |
| `add` | Add a recipe from a path or URL. `--source` is **required**; `-f, --force` overwrites. |
| `remove` | Remove a recipe. |

#### `asset-class`

```
refgenie asset-class list
refgenie asset-class show ASSET-CLASS-NAME [--asset-class-version VERSION]
refgenie asset-class add --source PATH_OR_URL [-f]
refgenie asset-class remove ASSET-CLASS-NAME [--asset-class-version VERSION]
```

| Subcommand | Description |
| --- | --- |
| `list` | List local asset classes. |
| `show` | Display an asset class. |
| `add` | Add an asset class from a path or URL. `--source` is **required**; `-f, --force` forces the action. |
| `remove` | Remove an asset class. |

#### `stage`

Manage staged assets — the archive area that `refgenie push` uploads from.

```
refgenie stage add ASSET-REGISTRY-PATHS [...]
refgenie stage remove ASSET-REGISTRY-PATHS [...]
refgenie stage list
```

| Subcommand | Description |
| --- | --- |
| `add` | Stage an asset. |
| `remove` | Unstage an asset. |
| `list` | List staged assets, with digest, name, mode, and size. |

### Asset definitions

#### `data-channel`

```
refgenie data-channel add NAME TYPE INDEX-ADDRESS [-d DESCRIPTION]
refgenie data-channel list
refgenie data-channel show NAME
refgenie data-channel validate NAME
refgenie data-channel sync NAME [--exists-ok | --exists-overwrite]
refgenie data-channel remove NAME
```

`add` positional arguments: the channel name, its type, and the address of its
index YAML file.

| `add` option | Description |
| --- | --- |
| `-d, --description` | Description of the data channel. |
| `--username` | Username for authentication. |
| `--password` | Password for authentication. |
| `--token` | Authentication token. |

| `sync` option | Description |
| --- | --- |
| `--exists-ok` | Skip existing assets/recipes without error. |
| `--exists-overwrite` | Delete conflicting items before adding. |

`--exists-ok` and `--exists-overwrite` are mutually exclusive.

#### `generate snakefile`

Generate a Snakemake file from the refgenie configuration. Requires the
`snakemake` extra to run the result.

```
refgenie generate snakefile -o OUTPUT_PATH [-s TEMPLATE_PATH]
```

| Option | Description |
| --- | --- |
| `-o, --output-path` | Path to save the generated Snakefile. **Required.** |
| `-s, --snakefile-template-path` | Path to the Snakefile template. |

#### `remote`

Configure the cloud destinations that `refgenie push` uploads to. A remote is
identified by its type.

```
refgenie remote add --type {s3,http,https} --prefix PREFIX --description TEXT [--push-command CMD]
refgenie remote list
refgenie remote status [-r REMOTE]
refgenie remote remove --type {s3,http,https}
```

| `add` option | Description |
| --- | --- |
| `--type` | Type of the remote: `s3`, `http`, or `https`. **Required.** |
| `--prefix` | Prefix/identifier for the remote. **Required.** |
| `--description` | Description of the remote. **Required.** |
| `--push-command` | Shell command template for pushing assets. Placeholders: `{local_path}`, `{relative_path}`, `{prefix}`, `{genome_stage_folder}`. Example: `aws s3 cp {local_path} s3://bucket/{relative_path}`. |

| `status` option | Description |
| --- | --- |
| `-r, --remote` | Show status for only this remote (by name or id). |

| `remove` option | Description |
| --- | --- |
| `--type` | Type of the remote to remove. **Required.** |

### Sequences

#### `getseq`

Retrieve a sequence region from a genome. Coordinates are 0-based and half-open.

```
refgenie getseq -g GENOME -l LOCUS
```

| Option | Description |
| --- | --- |
| `-g, --genome` | Reference assembly ID, e.g. `mm10`. **Required.** |
| `-l, --locus` | Coordinates of the desired sequence, e.g. `chr1:50000-50200`. **Required.** |

## Related documentation

- [README](../README.md) — installation, extras, environment variables, and data channels
- [Migrating from legacy refgenie](migration-from-legacy.md) — 0.x to 1.0 command and API changes
- [Build tutorial](refgenie_cli_build_tutorial.md) — a longer end-to-end build walkthrough
