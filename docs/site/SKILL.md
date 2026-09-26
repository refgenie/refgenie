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

Refgenie manages the storage, access, and transfer of reference genome
resources. It gives you command-line and Python interfaces to **download**
pre-built genome "assets" (like the indexes bioinformatics tools need) and to
**build** those assets for custom genomes. It also provides a standard folder
structure so software can swap from one genome to another without changing
code.

Think of it as "GitHub for reference genomes":

- `refgenie pull hg38/bwa_index` — download a pre-built asset from a server.
- `refgenie build mygenome/bowtie2_index` — build an asset locally.
- `refgenie seek hg38/salmon_index` — get the local path to an asset, so
  pipelines stay portable across machines.

Refgenie is backed by a database (SQLite by default, PostgreSQL for scale) and
identifies genomes by **sequence-derived GA4GH digests**, not just arbitrary
names like "hg38", so genome identity is verifiable.

Full documentation: <https://docs.refgenie.org>

---

## Command naming during the transition

Refgenie has been rewritten as a single unified package (the old `refgenie` +
`refgenconf` + `refgenieserver` are end-of-life). During the current
transition, the package on PyPI is published as **`refgenie1`** and the command
may be **`refgenie1`**. This document writes every command as `refgenie` for
future-proofing.

**If `refgenie` is not found, substitute `refgenie1`** (both the pip package
and the command). To find out which you have:

```bash
refgenie --version || refgenie1 --version
```

---

## Fastest path for an AI assistant: the MCP server

If your goal is to *answer questions* about a user's genomes and assets (not
change anything), use the built-in **read-only MCP server** instead of shelling
out. It never modifies data.

Register it once (Claude Code):

```bash
claude mcp add refgenie refgenie-mcp
```

Or in Claude Desktop / any stdio MCP client, point the client at the command
`refgenie-mcp` (no arguments — it reads the user's default refgenie config).

Read-only tools it exposes:

| Tool | Purpose |
|------|---------|
| `list_genomes` | All genomes with aliases, species, description |
| `search_genomes` | Search genomes by species/alias/description substring |
| `get_genome` | Detailed genome info by alias or digest, incl. asset groups |
| `list_asset_classes` | Registered asset classes, seek keys, serving modes |
| `list_recipes` | Registered recipes, output asset class, inputs |
| `list_assets` | Assets, optionally filtered by genome and/or asset class |
| `get_asset` | Asset detail by digest: seek keys, parents, children |
| `lookup_digest` | Universal digest lookup (genome first, then asset) |
| `get_genome_metadata` | Seqcol metadata: #sequences, total length, source |
| `get_genome_sequences` | Sequence-level names, lengths, sequence digests |
| `compare_genomes` | Seqcol comparison of two genomes |

Genomes may be referenced by alias (`hg38`) or digest; the server resolves
aliases automatically.

For anything that *changes state* (pull, build, init, alias edits), use the CLI
below.

---

## Core concepts

- **Asset registry path**: `genome/asset_class:tag`, e.g. `hg38/bwa_index` or
  `hg38/fasta:default`. The tag defaults to `default` and is usually omitted.
- **Seek key**: a named file within an asset, addressed as
  `genome/asset.seek_key`. The `fasta` asset has seek keys `fasta` (the `.fa`),
  `fai` (the `.fa.fai` index), and `chrom_sizes` (the `.chrom.sizes` file):
  `refgenie seek hg38/fasta.fai`.
- **Alias vs digest**: users say `hg38`; refgenie tracks the genome by its
  GA4GH sequence-collection digest and maps aliases to it. Two genomes built
  from the same sequences get the same digest, anywhere.
- **RefgetStore**: a content-addressable store (`~/.refgenie/genomes/.refget_store/`)
  where every sequence is keyed by its GA4GH SHA-512/24 digest. Filled by
  `genome init`; enables dedup, verification, and direct subsequence retrieval
  with `getseq` — without a FASTA on disk.
- **Config location**: defaults to `~/.refgenie/` (SQLite DB + `genomes/`
  folder). Override the base with `export REFGENIE_HOME_PATH=~/my_genomes`.

---

## Download workflow (use assets from a server)

```bash
# 1. One-time setup: create config + local SQLite database
refgenie init

# 2. Subscribe to a server that hosts pre-built assets
refgenie subscribe -s https://api.refgenie.org

# 3. Browse what's available (filter with -g GENOME)
refgenie listr
refgenie listr -g hg38

# 4. Download an asset (resolves the alias, fetches, extracts, registers it)
refgenie pull hg38/fasta
refgenie pull --genome hg38 fasta bowtie2_index bwa_index   # several at once

# 5. See what you have locally
refgenie list

# 6. Get a portable local path to an asset (the key value of refgenie)
refgenie seek hg38/bwa_index
refgenie seek hg38/fasta.fai          # a specific file via its seek key
```

Use `seek` in pipelines so paths never get hardcoded:

```bash
bwa mem $(refgenie seek hg38/bwa_index) reads.fq > aligned.sam
```

**Remote mode** — get a path to a file on the server (e.g. S3) *without*
downloading it, for cloud workflows:

```bash
refgenie seekr hg38/fasta.fai
```

---

## Build workflow (custom or unavailable genomes)

Register a genome, then build assets from it. `genome init` loads sequences
into the RefgetStore; `build` exports/derives asset files on disk.

```bash
# Register a genome from a local FASTA (.fa/.fa.gz/.fasta/.fasta.gz).
# This loads sequences into the RefgetStore and computes the seqcol digest.
refgenie genome init \
  --fasta /path/to/mygenome.fa.gz \
  --name mygenome \
  --description "My custom genome" \
  --species "Homo sapiens"

refgenie genome list                  # verify registration

# Build the FASTA asset: exports {digest}.fa, .fa.fai (computed natively —
# no samtools needed), and .chrom.sizes
refgenie build mygenome/fasta

# Build a derived index from the fasta asset
refgenie build mygenome/bowtie2_index

# Seek the built files
refgenie seek mygenome/fasta            # the .fa
refgenie seek mygenome/fasta.fai        # the index
refgenie seek mygenome/fasta.chrom_sizes
```

You can also build the FASTA asset directly with per-file inputs:

```bash
refgenie build mygenome/fasta --files fasta=/path/to/mygenome.fa.gz \
  --genome-description "My custom genome"
```

---

## Retrieve sequences directly (no FASTA needed)

`getseq` pulls sequence data straight from the RefgetStore. Coordinates are
**0-based, half-open**.

```bash
refgenie getseq -g hg38 -l chr1              # whole chromosome
refgenie getseq -g hg38 -l chr1:0-1000       # a subsequence
refgenie getseq -g hg38 -l chr1:50000        # from a position to the end
```

If the genome was initialized from a remote source, `getseq` fetches on demand
and caches locally.

---

## Compare genomes

```bash
refgenie compare hg38 hg38_custom
```

Reports which sequences are shared vs unique, and whether two genomes differ
only in sequence names/order or in actual sequence content (via GA4GH seqcol
comparison).

---

## Python API

The `Refgenie` object mirrors the CLI against the configured database:

```python
from refgenie import Refgenie

rgc = Refgenie()
rgc.asset.seek("hg38", "bowtie2_index")     # local path to an asset
```

---

## Aliases

```bash
refgenie alias get                                        # list aliases
refgenie alias set --aliases hg38 GRCh38 --digest <SEQCOL_DIGEST>
refgenie alias remove --aliases old_name
```

Aliases are resolved to the correct digest automatically when you `pull`.

---

## Serving assets (advanced)

Refgenie includes a built-in server — no separate package needed.

```bash
refgenie stage ...      # prepare built assets for serving (file or archive mode)
refgenie push ...       # push staged assets to a cloud remote (e.g. S3)
refgenie serve          # start the HTTP server (REST API + web interface)
refgenie dash           # start the local dashboard UI for browsing assets
```

Install extras when needed: `pip install "refgenie1[dash]"` or
`pip install "refgenie1[server]"`.

See <https://docs.refgenie.org/refgenie/building_tutorial> for the full
build → stage → serve → push lifecycle.

---

## Command reference (top level)

**Asset management:** `list`, `asset`, `seek`, `add`, `remove`, `rename`, `id`,
`build`, `populate`
**Remote operations:** `listr`, `seekr`, `pull`, `push`, `mirror`, `populater`,
`compare`
**Genome management:** `genome` (`init`, `list`, `remove`, `browse`, `sync`,
`set-metadata`)
**Server:** `serve`, `dash`, `subscribe`, `unsubscribe`, `catalog-export`
**Configuration:** `init`, `purge`, `config`, `alias`, `recipe`, `asset-class`,
`stage`
**Asset definitions:** `data-channel`, `generate`, `remote`
**Sequences:** `getseq`

Every command supports `--help`. Start with `refgenie --help`.

---

## Where to read more

- Getting started / CLI tutorial: <https://docs.refgenie.org/refgenie/cli_tutorial>
- Genomes and the RefgetStore: <https://docs.refgenie.org/refgenie/genome_tutorial>
- Connecting AI assistants (MCP): <https://docs.refgenie.org/refgenie/mcp>
- Building & serving: <https://docs.refgenie.org/refgenie/building_tutorial>
- Refget Python package (sequences & seqcol): <https://docs.refgenie.org/refget>
- Source: <https://github.com/refgenie/refgenie1>
