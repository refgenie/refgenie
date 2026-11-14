# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0a2] - 2026-11-14

### Added

- The genome page has Sequences, Extract, Download and Compare tabs: browse a
  genome's sequences, extract regions or a BED file of regions, download FASTA
  in the browser, and compare with another genome or your own FASTA file.
- Collection pages (`REFGENIE_COLLECTIONS`) show a branded slice of the
  catalog, such as the VGP assemblies.
- `GET /v4/genomes/{digest}/sequences` and `GET /v4/genomes/{digest}/compare/{other}`.
- `list --json` and `listr --json`; `getseq --fasta` and `--zero-based`.
- MCP: a `get_sequence` tool; `compare_genomes` returns a one-sentence summary and
  a map of renamed sequences; `refgenie-mcp --help`; and the instructions and
  `catalog_summary` say which server this is.
- `refgenie compare` and `GET /v4/genomes/{digest}/compare/{other}` give the same
  one-sentence summary and rename map as MCP; `compare --json` adds the full
  GA4GH comparison.
- `refgenie init` says where assets go and which server it subscribed to; an
  empty `list` says what to run next; `seek` suggests close matches (such as
  `hg38/fasta_index.chrom_sizes` for `hg38/chrom_sizes`) and the `pull` command.

### Changed

- `getseq` reads a text region as 1-based and inclusive, like samtools, UCSC
  and IGV, and the website's Extract tab: `chr17:43,044,295-43,044,394` is 100
  bases (commas are fine). `--zero-based` reads it as 0-based and half-open.
  Separate start and end numbers (the HTTP API, MCP, Python) stay 0-based and
  half-open, like GA4GH refget. getseq says on stderr how it read the region.
- `list` and `listr` name each genome once, by its shortest plain alias;
  `--json` has every alias and the digest.
- MCP `get_asset` over stdio gives `seek_command` (the asset is already here);
  on a server it gives `pull_command`.
- `REFGENIE_LOG_LEVEL=INFO` no longer prints a line for every HTTP request;
  `DEBUG` does.
- refgenie.org's install line follows the server's version.
- refgenie uses httpx2, the maintained successor of httpx, for HTTP.
- `refgenie pull` brings a parent asset only when the asset colocates its files
  (such as `bismark_bt2_index`). `refgenie pull hg38/fasta_index` now downloads
  12 KB, not the FASTA.
- A genome is read only from the store its `store_name` names; stores are no
  longer probed. A genome whose store the server does not serve is left out of
  every listing and answers 404.
- A publish catalog with a genome that names no store, or a store the catalog
  does not list, is refused at import.
- Compare returns the same GA4GH comparison from the server, the CLI and MCP.
- MCP `get_asset` builds its pull command from the shortest plain alias, gives
  absolute URLs, and lists each file in `files`; `get_genome_sequences` is paged
  and takes `name_contains`.
- `listr -g` and `compare` find genomes on the subscribed servers.
- `getseq` registers an unregistered genome first (metadata only).
- The default log level is WARNING (set `REFGENIE_LOG_LEVEL=INFO` for the old
  output). Results and hints print on stderr without the logger.
- `config get` prints YAML.
- `getseq` reads byte ranges and no longer caches sequences.
- `pull --init` accepts `--genome-digest`, and `--all-genomes` registers
  alias-less genomes by digest.
- `genome init --server` takes only server URLs (use `--store` for stores).
- `refgenie store sync` replaces each store's aliases as a set: a name the
  store dropped is removed, a moved name follows its collection, and names
  with '/' are skipped. A store that suddenly offers under half its aliases
  keeps them. Only a genome's owning store writes its description and other
  FHR columns.
- A catalog import rewrites the aliases of each genome the catalog carries,
  skipping names with '/', and refuses a catalog that would remove more than
  half of them. `refgenie catalog-prune` also deletes genomes the catalog no
  longer carries that have no assets, with their aliases, and disabled stores
  left empty.

### Fixed

- Pulling a genome by a second alias links all its assets under that alias, so
  `seek` no longer prints a path that does not exist.
- `SKILL.md` is served as `text/markdown` on macOS and on Python before 3.13.
- `genome init --server`, `alias set` and `id --remote` read a genome from the
  store that holds it, not only from the server's first store.
- Pulling a child asset (bwa_index, fasta_index, ...) of hg38, celegans and similar genomes no longer fails and deletes the genome when the server lists an alias containing '/'. A parent lands under the child's alias and adds none of the server's other aliases.
- `pull --init` now fails with a clear error, instead of saving an alias that cannot be resolved, when a genome's sequence collection is not on the server (#299).
- Quadratic `getseq` output time, whole-chromosome downloads for a few bases,
  truncated tables off a terminal, the download progress bar on stdout, and
  `build` exiting silently on an unknown alias.

### Removed

- `refgenie config set`, which was never implemented. `refgenie config` alone
  prints the configuration.
- `refgenie pull --asset-registry-paths`; pass registry paths as positional
  arguments.
- MCP: the `get_genome_metadata` tool (`get_genome_sequences` pages carry
  `n_sequences` and `total_length`), and `get_asset`'s `downloads` and `local_path`.

## [1.0.0a1] - 2026-11-06

First public alpha of refgenie 1.0, a rewrite that replaces the `refgenie` 0.x,
`refgenconf`, and `refgenieserver` packages with one package. Install it with
`pip install "refgenie>=1.0.0a1"`. It does not read 0.x configuration files or genome
folders; the [upgrade guide](./docs/refgenie/upgrade.md) migrates an existing
install, and [migration details](./docs/refgenie/migration-from-legacy.md) lists
what changed for pipelines and Python code.

### Added

- A database (SQLite by default, or PostgreSQL) holds genomes, assets, aliases,
  and settings.
- Genomes are identified by GA4GH sequence collection (seqcol) digests, with
  human-readable aliases such as `hg38` on top. Sequences live in a refget
  store, and `refgenie getseq` reads a sequence or region straight from it.
- Assets are stored by content digest. `refgenie seek` returns a readable path
  under `alias/`; `refgenie seek --abs` returns the digest path.
- Data channels: asset classes and recipes come from channels you subscribe to,
  and you can publish your own (`refgenie data-channel`, `recipe`,
  `asset-class`).
- Recipes can run in a container image under Docker or Apptainer
  (`refgenie build --container`, `REFGENIE_CONTAINER_ENGINE`).
- Asset classes choose how they are served: individual files, a `.tgz`
  archive, or metadata only. `refgenie stage`, `push`, and `remote` copy built
  assets to cloud storage for serving.
- A built-in server (`refgenie serve`) with a JSON API at `/v4` and GA4GH DRS
  and seqcol endpoints. One server can serve several refget stores
  (`refgenie store`). The web UI for the public server is
  at https://refgenie.org.
- `refgenie dash`: a local web dashboard for your genomes and assets, which
  runs pulls and builds as background jobs.
- An MCP server for AI assistants: `refgenie-mcp` over stdio, and `/mcp` on the
  HTTP server.
- Plugins: packages can hook into pull, build, and update events
  (`refgenie plugins`).
- Integrations for looper (`refgenie.integrations.looper.populate`) and
  Snakemake (`refgenie.integrations.snakemake`, `snakemake` extra).
- Optional extras: `dash`, `server`, `mcp`, and `snakemake`.

### Changed

- The default server is https://api.refgenie.org.
- `$REFGENIE` is no longer read. Point `$REFGENIE_DB_CONFIG_PATH` at the
  database settings file instead; see the README for every environment
  variable.
- Asset classes declare typed seek keys (`file`, `directory`, `prefix`,
  `string`, `json`). Aligner indexes return their prefix directly, so
  `hg38/bowtie2_index` replaces `hg38/bowtie2_index.dir` plus a hand-built
  suffix.
- The looper pre-submit hook is `refgenie.integrations.looper.populate`
  (was `refgenconf.looper_refgenie_populate`).
- In Python, use the `Refgenie` class in place of refgenconf's `RGC`.
  `seek` returns a `str`.
- Logs and errors go to stderr; stdout carries only command output.

### Removed

- The `genome_config.yaml` file. Settings live in the database.
- The implicit `.dir` seek key on every asset.
- The separate `refgenconf` and `refgenieserver` packages.

[Unreleased]: https://github.com/refgenie/refgenie/compare/v1.0.0a2...HEAD
[1.0.0a2]: https://github.com/refgenie/refgenie/compare/v1.0.0a1...v1.0.0a2
[1.0.0a1]: https://github.com/refgenie/refgenie/releases/tag/v1.0.0a1
