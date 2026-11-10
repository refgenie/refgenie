# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- The genome page has Sequences, Extract, Download and Compare tabs: browse a
  genome's sequences, extract regions or a BED file of regions, download FASTA
  in the browser, and compare with another genome or your own FASTA file.
- Collection pages (`REFGENIE_COLLECTIONS`) show a branded slice of the
  catalog, such as the VGP assemblies.
- `GET /v4/genomes/{digest}/sequences` and `GET /v4/genomes/{digest}/compare/{other}`.

### Changed

- A genome is read only from the store its `store_name` names; stores are no
  longer probed. A genome whose store the server does not serve is left out of
  every listing and answers 404.
- A publish catalog with a genome that names no store, or a store the catalog
  does not list, is refused at import.
- Compare returns the same GA4GH comparison from the server, the CLI and MCP.

### Fixed

- `pull --init` now fails with a clear error, instead of saving an alias that cannot be resolved, when a genome's sequence collection is not on the server (#299).

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

[Unreleased]: https://github.com/refgenie/refgenie/compare/v1.0.0a1...HEAD
[1.0.0a1]: https://github.com/refgenie/refgenie/releases/tag/v1.0.0a1
