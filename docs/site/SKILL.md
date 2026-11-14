---
name: refgenie
description: >-
  Manage reference genome assets with refgenie: download, build, and locate
  genome indexes and sequence files (FASTA, bwa/bowtie2/salmon/hisat2 indexes,
  chrom.sizes, and more), initialize genomes into a content-addressable
  RefgetStore, retrieve sequences by coordinate, and verify genome identity
  with GA4GH sequence-collection digests. Use whenever a task involves
  reference genomes, genome indexes, chromosome sizes, obtaining a portable
  path to a genome asset for a pipeline, or the refgenie CLI, Python API, or
  MCP server.
homepage: https://refgenie.org
---

# Refgenie

Refgenie manages reference genome resources. Its command line and Python API
**download** pre-built genome "assets" (like the indexes aligners need),
**build** assets for custom genomes, and **seek** the local path to any asset,
so pipelines can switch genomes without code changes:

- `refgenie pull hg38/bwa_index` downloads a pre-built asset from a server.
- `refgenie build mygenome/bowtie2_index` builds an asset locally.
- `refgenie seek hg38/salmon_index` prints the local path to an asset.

Refgenie stores its catalog in a database (SQLite by default, PostgreSQL for
large servers) and identifies genomes by **GA4GH digests computed from their
sequences**, so you can check that two genomes are the same.

Full documentation: <https://docs.refgenie.org/>. To use a refgenie server over
HTTP with nothing installed, read <https://api.refgenie.org/SKILL.md>.

---

## Install

Refgenie 1.0 is a pre-release, so ask for it by version:

```bash
pip install "refgenie>=1.0.0a2"
refgenie --version
```

Extras: `mcp` (the MCP server), `dash` (the local web dashboard), `server`
(`refgenie serve`, includes `dash` and `mcp`), and `snakemake`.

---

## The MCP server (read-only)

To *answer questions* about a user's genomes and assets, use the read-only MCP
server instead of shelling out:

```bash
pip install "refgenie[mcp]>=1.0.0a2"
claude mcp add refgenie "$(which refgenie-mcp)"
```

`"$(which refgenie-mcp)"` saves the full path, so the server starts even when
refgenie's virtual environment is not active; run the line with it active.
`refgenie-mcp` reads the same `REFGENIE_*` variables as the CLI.
`REFGENIE_HOME_PATH` picks the refgenie folder (default `~/.refgenie`). To use
another folder, put `-e` after the server name and `--` before the command:
`claude mcp add refgenie -e REFGENIE_HOME_PATH=/data/refgenie -- "$(which refgenie-mcp)"`.
If you set `REFGENIE_DB_CONFIG_PATH` for the CLI, pass it the same way. Other
clients take the full path as `command` and the variables in the server's `env`
block. The public server also offers the same tools over HTTP, with nothing
installed:

```bash
claude mcp add --transport http refgenie-public https://api.refgenie.org/mcp
```

| Tool | Purpose |
|------|---------|
| `catalog_summary` | counts of genomes, asset groups and assets, and the species that have assets |
| `list_genomes` | genomes, most assets first; filters `has_assets`, `taxon_id` |
| `search_genomes` | substring search over digest, alias, species names, other genome names, description, accession |
| `get_genome` | one genome by `alias` or `digest`, with its FHR metadata |
| `list_assets` | assets, optionally of one genome and/or `asset_class` |
| `get_asset` | one asset by digest: seek keys, parents, children, local path, `pull_command` |
| `list_asset_classes` | asset classes and their seek key names |
| `list_recipes` | recipes, optionally by `output_asset_class` |
| `get_recipe` | one recipe in full: inputs, command templates, docker image |
| `lookup_digest` | what a digest is: a genome, an asset, or nothing |
| `get_genome_sequences` | a genome's sequences (name, length, digest), paged; filter with `name_contains` |
| `get_sequence` | the bases of one region (`name`, `start`, `end`; 0-based, half-open; at most 100,000 bases) |
| `compare_genomes` | seqcol comparison of two genomes, with a one-sentence summary and a map of renamed sequences |

A tool that takes a genome has separate `alias` and `digest` parameters; give
exactly one. List tools page with a `token` (pass back
`pagination.next_page_token`).

To change anything (pull, build, init, alias edits), use the CLI below.

---

## Core concepts

- **Genome digest**: the GA4GH sequence-collection digest, computed from the
  sequences, their names and their order. The same named sequences in the same
  order give the same digest anywhere. Two genomes with the same sequences under
  different names (`chrI` vs `I`) have different digests; `compare` shows they
  match.
- **Alias**: a human name like `hg38` or `GRCh38` that points to a digest.
  Each refgenie instance keeps its own aliases. **An alias is not a digest.**
- **Asset group**: a named slot for one kind of asset on one genome, like
  `hg38`'s `bwa_index`. A group holds one or more assets and has a default.
- **Asset**: one built artifact in a group, with a name (`default` for
  `fasta`; the tool version for an index, like `bowtie2_index:2.5.5`)
  and its own asset digest.
- **Asset class**: the type of an asset, which defines its seek keys.
- **Seek key**: a named file or value inside an asset. The `fasta` class has
  seek keys `fasta` (the `.fa`), `fai` (the `.fa.fai`), and `chrom_sizes` (the
  `.chrom.sizes`).
- **Registry path**: the CLI's way to name an asset,
  `genome/asset_group.seek_key:asset`, e.g. `hg38/fasta.fai:default`. Only
  `genome/asset_group` is required; the seek key defaults to the asset's
  default seek key and the asset to the group's default asset. So
  `refgenie seek hg38/fasta` gives the `.fa` and `refgenie seek hg38/fasta.fai`
  gives its index.
- **RefgetStore**: a content-addressable sequence store
  (`~/.refgenie/genomes/.refget_store/`) where each sequence is keyed by its
  GA4GH digest. `genome init` fills it; `getseq` reads subsequences from it
  without a FASTA on disk.
- **Config location**: `~/.refgenie/` (SQLite database and `genomes/` folder).
  Change it with `export REFGENIE_HOME_PATH=~/my_genomes`.

---

## Download workflow (use assets from a server)

```bash
# 1. One-time setup: create the config and local database.
#    This also subscribes to the public server, https://api.refgenie.org.
refgenie init

# 2. Optional: subscribe to another server
refgenie subscribe -s https://my.server.org

# 3. Browse what's available (subscribed servers; no download)
refgenie listr
refgenie listr -g hg38

# 4. Download assets (resolves the alias, fetches, extracts, registers)
refgenie pull hg38/fasta
refgenie pull hg38/bowtie2_index hg38/bwa_index   # several at once

# 5. See what you have locally
refgenie list

# 6. Get the local path to an asset
refgenie seek hg38/bwa_index
refgenie seek hg38/fasta.fai          # a specific file, by seek key
```

Use `seek` in pipelines so paths are never hardcoded:

```bash
bwa mem $(refgenie seek hg38/bwa_index) reads.fq > aligned.sam
```

`pull` downloads only the asset you name, except that a few assets that need
their parent's files beside them (such as `bismark_bt2_index`) also pull that
parent. `seek` and `seekr` print only the path or URL. Set
`REFGENIE_LOG_LEVEL=INFO` to see progress, and add `--json` to `list` or `listr`
for machine-readable output.

`seekr` gives the path to a file on the server (e.g. on S3) without
downloading it, for cloud workflows:

```bash
refgenie seekr hg38/fasta.fai
```

---

## Build workflow (custom or unavailable genomes)

Register the genome from a FASTA (`.fa`, `.fa.gz`, `.fasta`, `.fasta.gz`), then
build assets from it. `genome init` loads the sequences into the RefgetStore,
computes the seqcol digest, and builds the `fasta` asset (`.fa`, `.fa.fai`,
`.chrom.sizes`); add `--no-build` to skip that build.

```bash
refgenie genome init \
  --fasta /path/to/mygenome.fa.gz \
  --name mygenome \
  --description "My custom genome" \
  --species "Homo sapiens"

refgenie genome list                    # check it is registered

refgenie build mygenome/bowtie2_index   # build a derived index

refgenie seek mygenome/fasta            # the .fa
refgenie seek mygenome/fasta.fai        # the index
refgenie seek mygenome/fasta.chrom_sizes
```

---

## Recipes

### The aligner that matches an index

Refgenie gives you the index, not the aligner. The asset name is the tool
version that built it.

```bash
refgenie seek yeast_s288c/bowtie2_index.version   # 2.5.5
conda install -c bioconda bowtie2=2.5.5           # or any installer; match the version
bowtie2 -x "$(refgenie seek yeast_s288c/bowtie2_index)" -U reads.fq -S out.sam
```

### Chromosome sizes without the genome

```bash
refgenie pull hg38/fasta_index                 # about 12 KB: the .fai and .chrom.sizes
refgenie seek hg38/fasta_index.chrom_sizes
```

Use `fasta_index.chrom_sizes`. It is the same file as `fasta.chrom_sizes`,
without the 3 GB FASTA. Over MCP: page through `get_genome_sequences`.

### A region of a genome you have not downloaded

```bash
refgenie getseq -g hg38 -l chr17:43044294-43044394   # first 100 bp of BRCA1; 0-based, half-open
refgenie getseq -g hg38 -l chr17:43044295-43044394 --one-based   # same 100 bp, 1-based inclusive
refgenie getseq -g hg38 -l chr17:43044294-43044394 --fasta > brca1_start.fa
```

The first call registers hg38 from your server; no FASTA is downloaded. Over
MCP: `get_sequence(alias="hg38", name="chr17", start=43044294, end=43044394)`.

### Check whether two genomes match

```bash
refgenie compare hg38 hg38_noalt_decoy   # works for genomes you have not pulled
```

Same `sequences` and `names`: the same genome. Same `sequences`, different
`names`: renamed. Shared `name_length_pairs`: the same coordinates.

Not sure which hg38 to pull? See the [choosing hg38 page](https://docs.refgenie.org/refgenie/choosing_hg38/).

---

## Request a genome or an asset

- A genome refgenie does not list: open
  https://github.com/refgenie/refgenie-registry/issues/new?template=genome_proposal.yml
- A listed genome missing an asset: open
  https://github.com/refgenie/refgenie-registry/issues/new?template=asset_request.yml&genome=<alias>&genome_digest=<digest>&title=Asset+request:+<alias>
- Need it now: build it yourself (see "Build workflow").

---

## Python API

The `Refgenie` object gives the same operations as the CLI, against the
configured database:

```python
from refgenie import Refgenie

rgc = Refgenie()
digest = rgc.alias.resolve("hg38")  # seek takes a genome digest
rgc.asset.seek(digest, "bowtie2_index")  # local path to an asset
```

---

## Aliases

```bash
refgenie alias get                                    # list aliases
refgenie alias set -a hg38,GRCh38 -d <SEQCOL_DIGEST>  # several aliases: comma-separated
refgenie alias remove -a old_name
```

`pull` resolves aliases to digests automatically.

---

## Serving assets

```bash
refgenie stage add hg38/fasta   # mark built assets for serving
refgenie push                   # upload staged assets to a cloud remote (e.g. S3)
refgenie serve                  # start the HTTP server (REST API, web UI, /mcp)
refgenie dash                   # start the local dashboard on port 8080
```

`serve` needs `pip install "refgenie[server]>=1.0.0a2"`; `dash` needs
`"refgenie[dash]"`. See <https://docs.refgenie.org/refgenie/building_tutorial/>
for the full build, stage, push, and serve steps.

---

## Command reference (top level)

**Asset management:** `list`, `asset`, `seek`, `add`, `remove`, `rename`, `id`,
`build`, `populate`
**Remote operations:** `listr`, `seekr`, `pull`, `push`, `mirror`, `populater`,
`compare`
**Genome management:** `genome` (`init`, `list`, `remove`, `browse`, `sync`,
`set-metadata`)
**Server:** `serve`, `dash`, `subscribe`, `unsubscribe`, `catalog-export`, `catalog-prune`
**Configuration:** `init`, `purge`, `config`, `plugins`, `alias`, `recipe`,
`asset-class`, `stage`
**Asset definitions:** `data-channel`, `generate`, `remote`
**Federation:** `store`
**Sequences:** `getseq`

Every command takes `--help`. Start with `refgenie --help`.

---

## Where to read more

- CLI tutorial: <https://docs.refgenie.org/refgenie/cli_tutorial/>
- Genomes and the RefgetStore: <https://docs.refgenie.org/refgenie/genome_tutorial/>
- Connecting AI assistants (MCP): <https://docs.refgenie.org/refgenie/mcp/>
- Building and serving: <https://docs.refgenie.org/refgenie/building_tutorial/>
- Refget Python package (sequences and seqcol): <https://docs.refgenie.org/refget/>
- Using a refgenie server over HTTP, with nothing installed: <https://api.refgenie.org/SKILL.md>
- Source: <https://github.com/refgenie/refgenie>
