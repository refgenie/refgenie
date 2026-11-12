# Refgenie CLI

The installed command is `refgenie`. This page walks through a first session,
then lists every command and flag.

## Installation

Install from [PyPI](https://pypi.org/project/refgenie/). Refgenie 1.0 is a
pre-release, so ask for it by version:

```bash
pip install "refgenie>=1.0.0a1"
```

The base install provides the CLI and the Python API. `refgenie dash`,
`refgenie serve`, `refgenie-mcp`, and snakemake-based bulk building each require
an extra (`dash`, `server`, `mcp`, `snakemake`); see [Install](install.md).
Running a command without its extra prints a message naming the extra to install.

Check the installed version:

```bash
$ refgenie --version
refgenie 1.0.0a1
```

## Getting started

Refgenie stores asset metadata in a database: SQLite by default, or PostgreSQL.
With no configuration, refgenie creates a `refgenie` SQLite file under
`~/.refgenie` (or `$REFGENIE_HOME_PATH`). Start with `refgenie init`, which also
creates the genome folder and the stage folder.

```bash
$ refgenie init
INFO     Database configuration file created at /home/user/.refgenie/refgenie_db_config.yaml. It uses SQLite by default; edit it only to switch to PostgreSQL.
INFO     Genome folder ready: /home/user/.refgenie/genomes
INFO     Genome stage folder ready: /home/user/.refgenie/archives
INFO     Initialized refgenie backend: 'sqlite:////home/user/.refgenie/refgenie'
INFO     Subscribed to the default server: https://api.refgenie.org
```

Asset classes and recipes, including `fasta`, come from a data channel, so
register and sync one before building anything.

### Order of operations

1. `refgenie init` sets up the config, database, and folders.
2. `refgenie data-channel add ...` then `refgenie data-channel sync ...`
   register asset classes and recipes.
3. `refgenie genome init ...` registers a genome and builds its `fasta` asset,
   once step 2 has registered a `fasta` recipe. `--no-build` skips the build.
4. `refgenie build ...` builds other assets (e.g. `bwa_index`).

### Register a data channel

A data channel is an index of asset class and recipe definitions. The main
channel is [refgenie-registry](https://refgenie.github.io/refgenie-registry).

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

To build the `fasta` asset separately, for instance after `genome init
--no-build`, call `build` with the source file:

```bash
refgenie build rCRSd/fasta --files fasta=rCRSd.fa
```

### Build another asset

Once the genome exists and the recipe is registered, other assets build from it:

```bash
refgenie build rCRSd/bwa_index
```

Recipes run real tools. If `bwa` is not on `PATH`, either install it or pass
`--container` to run the recipe in its container image, with Docker or
Apptainer. `-q/--requirements` prints what a recipe needs without building:

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

`<genome>` is an alias (e.g. `hg38`). Omitting `:<tag>` uses the default tag.
The `.seek_key` suffix selects a specific file within an asset (e.g.
`hg38/fasta.fai`).

To name a genome by its sequence-collection digest instead, leave the genome out
of the path and pass `--genome-digest`:

```
refgenie seek fasta --genome-digest DIGEST
```

This reaches a genome that has no alias. A digest written into the path is read
as an alias, and naming the genome both ways at once is an error. `build` takes
an alias only, because the build folder is named after it.

## Command reference

Results go to stdout; messages and logs go to stderr. Logs show warnings and
errors by default; set `REFGENIE_LOG_LEVEL=INFO` for more.

`refgenie --help` groups commands as follows.

| Group | Commands |
| --- | --- |
| Asset management | `list`, `asset`, `seek`, `add`, `remove`, `rename`, `id`, `build`, `populate` |
| Remote operations | `listr`, `seekr`, `pull`, `push`, `mirror`, `populater`, `compare` |
| Genome management | `genome` |
| Server | `serve`, `dash`, `subscribe`, `unsubscribe`, `catalog-export`, `catalog-prune` |
| Configuration | `init`, `purge`, `config`, `plugins`, `alias`, `recipe`, `asset-class`, `stage` |
| Asset definitions | `data-channel`, `generate`, `remote` |
| Federation | `store` |
| Sequences | `getseq` |

Give several values to a list option comma-separated (`-a hg38,GRCh38`) or by
repeating the option (`-a hg38 -a GRCh38`). Values separated by spaces are an
error. Registry paths are positional and take spaces:
`refgenie pull hg38/fasta hg38/bwa_index`.

Boolean flags shown as `-f, --force` also accept the negated form
(`--no-force`). Command names use hyphens (`asset-class`, `data-channel`,
`catalog-export`).

### Asset management

#### `list`

List available local assets.

```
refgenie list [-g ALIAS[,ALIAS...]] [--genome-digest DIGEST[,DIGEST...]]
```

| Option | Description |
| --- | --- |
| `-g, --genome` | One or more genome aliases to restrict the listing to. |
| `--genome-digest` | One or more genome digests to restrict the listing to. |
| `--json` | Print JSON instead of a table. |

With `--json`, `list` prints one entry per source, here `local`. The `assets`
strings are `group:name`, so `<alias>/<asset>` is a registry path:

```json
[{"source": "local", "genomes": [{"digest": "AHyVecap0miRdrz6ZffizkZQUOJJBO2c", "aliases": ["athaliana"], "assets": ["fasta:default", "bowtie2_index:default"]}]}]
```

#### `asset list`

Alias group. `refgenie asset list` is equivalent to `refgenie list` and takes
the same `-g, --genome` and `--genome-digest` options.

#### `seek`

Print the local path of an asset.

```
refgenie seek ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `-e, --check-exists` | Check the returned path for existence on disk. |
| `--abs` | Return the digest-addressed content path under `data/` instead of the human-readable alias path. |
| `--genome-digest` | Name the genome by digest instead of by the alias in the registry path; leave the genome out of the path (e.g. `seek fasta --genome-digest DIGEST`). Reaches a genome with no alias. |

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
| `--genome-digest` | Name the genome by digest instead of by the alias in the registry path; leave the genome out of the path (e.g. `seek fasta --genome-digest DIGEST`). Reaches a genome with no alias. |

#### `remove`

Remove a local asset.

```
refgenie remove ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `-f, --force` | Do not prompt before removing. |
| `-a, --aliases` | Also remove the genome alias if this was the genome's last asset. |
| `--genome-digest` | Name the genome by digest instead of by the alias in the registry path; leave the genome out of the path (e.g. `seek fasta --genome-digest DIGEST`). Reaches a genome with no alias. |

#### `rename`

Rename an asset.

```
refgenie rename ASSET-REGISTRY-PATHS [...] -n NEW_ASSET_NAME
```

| Option | Description |
| --- | --- |
| `-n, --new-asset-name` | New name for the asset. **Required.** |
| `--genome-digest` | Name the genome by digest instead of by the alias in the registry path; leave the genome out of the path (e.g. `seek fasta --genome-digest DIGEST`). Reaches a genome with no alias. |

#### `id`

Return a digest. A genome alias yields the genome digest; a registry path yields
the asset digest. A positional name is always an alias; name a genome by digest
with `--genome-digest`.

```
refgenie id REGISTRY-PATHS [...]
refgenie id --genome-digest DIGEST [--remote]
```

| Option | Description |
| --- | --- |
| `--genome-digest` | Name the genome by digest instead of by alias. Alone, prints the digest if the genome is known; with a registry path, leave the genome out of it (e.g. `id fasta --genome-digest DIGEST`). |
| `-v, --verbose` | Show detailed genome metadata (sequence count, total length, source). |
| `--validate-store` | Verify that the genome's RefgetStore exists and is valid. |
| `--remote` | With `--genome-digest`, query subscribed servers when the genome is not found locally. |
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
| `--container` | Run the recipe's commands in its container image, if it has one, with Docker or Apptainer (see `REFGENIE_CONTAINER_ENGINE`). |
| `--pull-parents` | Automatically pull a required parent asset that was not provided. |
| `-q, --requirements` | Show the build requirements for the asset and exit. |
| `--stage` | Stage the asset after building. Requires the genome stage folder to be set. |
| `--push-to` | Remotes, by name or id, to create push intent records for after staging. An unknown remote is an error. |
| `--pipeline-kwargs` | Extra arguments for the build pipeline, as `arg_name=arg_val`. |
| `--assets` | Override the genome, asset, and tag of parents, e.g. `fasta=hg38/fasta:default`. |
| `--files` | Paths to required input files, e.g. `fasta=/path/to/file.fa.gz`. |
| `--params` | Required parameter values, e.g. `param1=value1`. |
| `--volumes` | With `--container`, also make these folders visible inside the container. |

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
refgenie listr [-g ALIAS[,ALIAS...]] [--genome-digest DIGEST[,DIGEST...]] [-s URL[,URL...]]
```

| Option | Description |
| --- | --- |
| `-g, --genome` | One or more genome aliases to restrict the listing to. A name not registered here is looked up on the servers. |
| `--genome-digest` | One or more genome digests to restrict the listing to. The genomes need not exist locally. |
| `-s, --genome-server` | One or more server URLs to use for this call only; not persisted to config. |
| `-p, --append-server` | Append the provided servers to the configured list rather than replacing it. |
| `--json` | Print JSON instead of a table: a list with one entry per server, in the shape `list --json` uses, with the server URL as `source`. |

#### `seekr`

Print the remote path of an asset.

```
refgenie seekr ASSET-REGISTRY-PATHS [...]
```

| Option | Description |
| --- | --- |
| `-s, --genome-server` | One or more server URLs to use for this call only; not persisted. |
| `-p, --append-server` | Append the provided servers to the configured list. |
| `--genome-digest` | Name the genome by digest instead of by the alias in the registry path; leave the genome out of the path (e.g. `seek fasta --genome-digest DIGEST`). Reaches a genome with no alias. |

#### `pull`

Download assets from subscribed servers. Name assets by registry path.

```
refgenie pull GENOME/ASSET [...]
```

| Option | Description |
| --- | --- |
| `-g, --genome` | Genome alias(es), e.g. `mm10`; comma-separate several for `--all` or `--asset`. |
| `--genome-digest` | Genome digest(s), in place of aliases. With registry paths, one digest, and leave the genome out of the path (`pull fasta --genome-digest DIGEST`). With `--all` or `--asset`, comma-separate several. |
| `--all` | Pull all assets for the specified genome(s). |
| `--all-genomes` | Apply the operation to all genomes available on subscribed servers. |
| `--asset` | Pull one asset type across the specified genomes, e.g. `--asset fasta`. |
| `--init` | Register genome(s), named with `-g` or `--genome-digest`, without downloading asset files. |
| `--skip-large` | Do not pull archives over the size cutoff. |
| `--pull-large` | Pull all archives regardless of size. |
| `--size-cutoff` | Maximum archive size, in GB, to pull without confirmation. Default `10`. |
| `--batch` | Batch mode: pull all archives regardless of size. |
| `-f, --force` | Skip confirmation prompts for multi-asset operations. |

#### `push`

Upload staged assets to cloud remotes.

```
refgenie push [-r REMOTE] [-g ALIAS | --genome-digest DIGEST]
```

| Option | Description |
| --- | --- |
| `-r, --remote` | Push only to this remote (by name or id). Default: all remotes with unpushed assets. |
| `-g, --genome` | Push only assets for the genome with this alias. |
| `--genome-digest` | Push only assets for the genome with this digest. |
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

Compare two genomes for compatibility. A genome not registered here is
compared on the subscribed server that has it.

```
refgenie compare ALIAS1 ALIAS2
refgenie compare ALIAS1 --genome-digest DIGEST2
refgenie compare --genome-digest DIGEST1,DIGEST2
```

| Option | Description |
| --- | --- |
| `--genome-digest` | Genome digest(s) to compare, in place of aliases. Aliases and digests together must name exactly two genomes. |

### Genome management

#### `genome init`

Initialize a genome from a FASTA file, a refgenie server, or a RefgetStore,
and build its `fasta` asset (`fa`, `fai`, `chrom.sizes`). For a genome from a
server or store, the build first downloads its sequences.

```
refgenie genome init -n NAME [--fasta PATH | --server URL | --store URL]
```

| Option | Description |
| --- | --- |
| `-n, --name` | One or more alias names for the genome. **Required.** |
| `--fasta` | Path to a local FASTA file. |
| `--server` | URL of a refgenie server to initialize from. For a RefgetStore URL, use `--store`. |
| `--store` | URL of a RefgetStore to initialize from; requires `--digest` or `--namespace`. |
| `--namespace` | Namespace for alias lookup when using `--store`. |
| `--digest` | Seqcol digest of the genome. |
| `-d, --description` | Genome description, e.g. `Human genome build 38`. |
| `-s, --species` | Species name, e.g. `Homo sapiens`. |
| `--fhr` | Path to an FHR `.fhr.json` metadata file to apply after init. When given it is authoritative for description and species (species comes from the record's `taxon.name`) and writes the RefgetStore sidecar. |
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
| `--fhr` | Path to the FHR `.fhr.json` metadata file to apply. **Required.** It is authoritative for description and species (species comes from the record's `taxon.name`). |

#### `genome list`

List all genomes, with digests, aliases, source, species, and description.

```
refgenie genome list
```

#### `genome remove`

Remove a genome and all its assets.

```
refgenie genome remove (--genome ALIAS[,...] | --genome-digest DIGEST[,...])
```

| Option | Description |
| --- | --- |
| `--genome` | Genome alias(es) to remove. |
| `--genome-digest` | Genome digest(s) to remove. One of the two is required. |
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

Add refgenie server URLs to the subscription list.

```
refgenie subscribe -s URL[,URL...]
```

| Option | Description |
| --- | --- |
| `-s, --genome-server` | One or more URLs to add to the subscription list. |
| `-r, --reset` | Overwrite the current list of server URLs. |

#### `unsubscribe`

Remove refgenie server URLs from the subscription list.

```
refgenie unsubscribe -s URL[,URL...]
```

| Option | Description |
| --- | --- |
| `-s, --genome-server` | One or more URLs to remove from the subscription list. |

#### `catalog-export`

Export a publish catalog covering pushed assets only, for a server to import.

```
refgenie catalog-export [--dest PATH] [--https-prefix URL]
    [--catalog-store-name NAME --catalog-store-url URL [--catalog-store-priority N]]
```

| Option | Description |
| --- | --- |
| `--dest` | Path to write the publish-catalog SQLite artifact to. |
| `--https-prefix` | Public https base URL mirroring the stage folder that `refgenie push` uploaded to, e.g. `https://<bucket>.s3.amazonaws.com/assets`. Download links are served from here. |
| `--catalog-store-name` | Name of the public store that holds the genomes this node built itself, e.g. `refgenie-main`. Those genomes belong to no registered store; the catalog lists them under this one. |
| `--catalog-store-url` | Public root URL of that store. Required with `--catalog-store-name`. |
| `--catalog-store-priority` | Priority of that store. Defaults to `30`. |

The catalog carries the [store list](#store): every enabled `remote` store in
the database, plus the catalog store if one is named. `on_disk` stores are left
out, since their URL is a path on this machine. A disabled store is left out
too, along with its genomes and their aliases, assets and download records; the
export logs how many genomes it omitted for each such store. Enable the store
again (`refgenie store apply`) and the next export carries them. A server that
imports the catalog adds or updates the listed stores and disables any store
the catalog does not list. The export fails if a genome belongs to any other
store it does not publish (an `on_disk` store, or one missing from the store
table), because the server could not import that genome.

#### `catalog-prune`

Delete what a server holds but its publish catalog no longer has. A catalog
import never deletes assets, genomes or stores, so run this on the server
after it imports the catalog, with `--dry-run` first. It removes:

- assets of genomes in the catalog that the catalog no longer lists, such as
  the old copy of an asset rebuilt on the build node, with their names, links
  and staged rows;
- genomes the catalog no longer carries that have no assets, with their
  aliases (such as the genomes of a retired store). A dropped genome that
  still has assets is listed and kept;
- disabled stores left with no genome.

A catalog with no genomes deletes no genomes.

```
refgenie catalog-prune [--source PATH_OR_URL] [--dry-run]
```

| Option | Description |
| --- | --- |
| `--source` | The publish catalog the server last imported (path or http(s) URL). Defaults to `$REFGENIE_CATALOG_URL`. |
| `--dry-run` | List what would be deleted without deleting it. |

### Configuration

#### `init`

Initialize the refgenie configuration, database, and folders.

```
refgenie init
```

| Option | Description |
| --- | --- |
| `-f, --genome-folder` | Absolute path to the parent folder for refgenie-managed assets. |
| `-a, --genome-stage-folder` | Absolute path to the stage folder, where staged assets wait for `refgenie push`. |
| `-v, --config-version` | Config version to initialize the config file with. |
| `--default-server` / `--no-default-server` | Subscribe a new configuration to `https://api.refgenie.org` (default on; init makes no network call). `--no-default-server`, or `REFGENIE_SUBSCRIBE_DEFAULT_SERVER=false`, starts with no subscriptions. |

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

`config get` prints the current configuration, including the database
connection and the environment-derived settings, as YAML. `config set` is not yet
implemented.

#### `plugins list` / `plugins set` / `plugins unset`

```
refgenie plugins [list]
refgenie plugins set PLUGIN KEY=VALUE [KEY=VALUE ...]
refgenie plugins unset PLUGIN [KEY ...]
```

`plugins list` (also bare `refgenie plugins`) shows every installed plugin:
its hook, entry-point name, target, package and status (`ok`, `disabled`,
`load error: ...`, or `unknown hook (never fires)` for an entry point in a
`refgenie.hooks.*` group refgenie does not fire). Below that it shows the stored
plugin settings; settings for a plugin that is not installed are marked
`(not installed)`. When this refgenie runs no plugins it says why:
`REFGENIE_DISABLE_PLUGINS`, or server mode.

`plugins set` stores settings for a plugin in the database, merging them with
what is already there. Each setting is `key=value`, split on the first `=`.
Settings may be stored before the plugin is installed; refgenie warns and
stores them anyway. `plugins unset` removes the given keys, or every setting
for the plugin when no key is given.

```bash
refgenie plugins set nfcore config_path=/abs/path/nf.config
refgenie plugins unset nfcore config_path
```

#### `alias get` / `alias set` / `alias remove`

```
refgenie alias get [-a ALIAS[,ALIAS...]] [-g DIGEST[,DIGEST...]]
refgenie alias set -a ALIAS[,ALIAS...] [-d DIGEST]
refgenie alias remove -a ALIAS[,ALIAS...]
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

Manage staged assets, the area that `refgenie push` uploads from.

```
refgenie stage add ASSET-REGISTRY-PATHS [...] [--push-to REMOTE[,REMOTE...]]
refgenie stage remove ASSET-REGISTRY-PATHS [...]
refgenie stage list
```

| Subcommand | Description |
| --- | --- |
| `add` | Stage an asset. `--push-to REMOTE` also marks it for `refgenie push`. |
| `remove` | Unstage an asset. |
| `list` | List staged assets, with digest, name, mode, and size. |

`stage add` and `stage remove` also take `--genome-digest DIGEST` to name the
genome by digest instead of by the alias in the registry path.

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
refgenie generate snakefile -o OUTPUT_PATH [-s TEMPLATE_PATH] [--container]
```

| Option | Description |
| --- | --- |
| `-o, --output-path` | Path to save the generated Snakefile. **Required.** |
| `-s, --snakefile-template-path` | Path to the Snakefile template. |
| `--container` | Make every build rule pass `--container`, so each recipe runs in its image. |

#### `remote`

Configure the cloud destinations that `refgenie push` uploads to. Every command
that takes a remote (`remote remove`, `remote status -r`, `push -r`,
`build --push-to`, `stage add --push-to`) accepts its name or its numeric id. Several remotes can share
a type; names are unique.

```
refgenie remote add --type {s3,http,https} --prefix PREFIX --name NAME [--push-command CMD]
refgenie remote list
refgenie remote status [-r REMOTE]
refgenie remote remove REMOTE
```

| `add` option | Description |
| --- | --- |
| `--type` | Type of the remote: `s3`, `http`, or `https`. **Required.** |
| `--prefix` | Prefix/identifier for the remote. **Required.** |
| `--name` | Name of the remote. Must be unique and not all digits, since digits read as an id. **Required.** |
| `--push-command` | Shell command template for pushing assets. Placeholders: `{local_path}`, `{relative_path}`, `{prefix}`, `{genome_stage_folder}`. Example: `aws s3 cp {local_path} s3://bucket/{relative_path}`. |

| `status` option | Description |
| --- | --- |
| `-r, --remote` | Show status for only this remote (by name or id). |

| `remove` argument | Description |
| --- | --- |
| `REMOTE` | Name or id of the remote to remove. **Required.** |

### Federation

#### `store`

Manage the registry of refget stores a server serves genomes from. When two
stores map the same alias to different genomes, the store with the lower
priority number wins.

```
refgenie store add --name NAME -u URL [-t {remote,on_disk}] [-p PRIORITY] [--description TEXT]
refgenie store apply -f FILE
refgenie store sync [--name NAME] [--page-size N]
refgenie store list
refgenie store remove --name NAME
refgenie store conflicts
```

| Subcommand | Description |
| --- | --- |
| `add` | Register a store. `--name` (unique; `local` is reserved) and `-u, --url` (store root URL, or a local `.refget_store` path with `-t on_disk`) are **required**. `-p, --priority` defaults to `100`. |
| `apply` | Make the store list match a JSON file. `-f, --file` is **required**: a JSON array of stores, each `{"name", "url", "type"?, "priority"?, "description"?, "enabled"?}`. Listed stores are added or updated; stores not listed are disabled, never deleted. Prints what was added, updated and disabled. |
| `sync` | Load a store's collections and aliases into the database. With no `--name`, syncs every enabled store. `--page-size` defaults to `100`. A store's aliases are replaced as a set: names it no longer offers are removed, and names with '/' are skipped. A store that offers under half the aliases it holds, or a sync with failures, leaves its aliases unchanged. |
| `list` | List registered stores. |
| `remove` | Remove a store from the registry. `--name` is **required**. |
| `conflicts` | Show aliases that different stores map to different genomes. |

### Sequences

#### `getseq`

Retrieve a sequence region from a genome. Coordinates are 0-based and half-open
unless you pass `--one-based`. It reads only the bases you ask for and caches
nothing locally.

```
refgenie getseq (-g ALIAS | --genome-digest DIGEST) -l LOCUS [--one-based] [--fasta]
```

| Option | Description |
| --- | --- |
| `-g, --genome` | Genome alias, e.g. `mm10`. A genome not registered here is first registered from the subscribed servers (metadata only, like `pull --init`). |
| `--genome-digest` | Genome digest, in place of an alias. One of `-g` or `--genome-digest` is required. |
| `-l, --locus` | Coordinates of the desired sequence, e.g. `chr1:50000-50200`. **Required.** |
| `--one-based` | Read the locus as 1-based, closed coordinates (samtools style). |
| `--fasta` | Print a FASTA record: a header line with the locus, then 60 bases per line. |

```bash
refgenie getseq -g hg38 -l chrM:1-60 --one-based   # same bases as chrM:0-60
refgenie getseq -g hg38 -l chrM:0-120 --fasta > chrM_start.fa
```

## Related documentation

- [Install](install.md): installation and extras
- [README](https://github.com/refgenie/refgenie/blob/master/README.md#environment-variables): environment variables
- [Migrating from legacy refgenie](migration-from-legacy.md): 0.x to 1.0 command and API changes
