# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Importing refgenie no longer changes Python's logging setup. Looper, which
  imports the refgenie hook, stopped printing its own progress after the first
  sample. Only the `refgenie` and `refgenie-mcp` commands and the server set up
  logging now.

## [1.0.0a3] - 2026-11-17

### Added

- A docs guide to using refgenie genomes in R: BiocRefgetStore turns any
  refgenie genome into a Bioconductor BSgenome package.
- `refgenie catalog-export --catalog-store-priority` sets the catalog store's
  priority.

### Changed

- A server's catalog import also prunes: genomes, assets and stores the catalog
  dropped leave on the next import, with no `refgenie catalog-prune` run. An
  import refuses a prune that would delete more than half the server's genomes
  or assets.
- refgenie installs next to current Snakemake 9.x: the sqlmodel minimum is
  0.0.37 again, and the `snakemake` extra takes any 9.x release. The `dash`
  extra needs fastapi 0.162.4 or later.
- A genome is named by the registry's name for it (`hg38_initial`), not a
  shorter synonym (`hg38-ddbj`), on the website, in compare results and in MCP
  pull targets.
- Website: a genome page with no assets opens on its Sequences tab; the Compare
  box picks a genome on Enter and says when none is chosen; comparing your own
  FASTA lists odd renames such as `chrM` to `Mito` first.
- Website: hg38, mm39 and yeast searches and genome pages link the docs guide
  on which one to use; the alpha notice shows the exact version to pin; the
  no-JavaScript text gives the install command.
- Website: the asset page breadcrumb names the genome, and small asset files
  download under the genome's name (`hg38.chrom.sizes`).
- Website: plainer words, with menu hover hints, no "asset groups" count and a
  simpler About intro.
- Docs: the install page says plain `pip install refgenie` gives the old 0.13,
  a table maps 0.x commands to 1.0, and the sidebar is shorter.

### Fixed

- `refgenie populater` fills `refgenie://genome/asset` with the default asset,
  like `populate`.
- `populate` and `populater` exit nonzero when a `refgenie://` path is left
  unfilled.
- An empty `$REFGENIE` no longer triggers the legacy-variable warning.
- `seek` errors name the genome as you typed it, not its digest. Python
  `asset.seek` says when it was given an alias instead of a digest.
- Old 0.x commands such as `refgenie seek -g hg38 -a fasta` print a hint with
  the 1.0 form, `refgenie seek hg38/fasta`.
- The `compare` sentence lists renames that break the pattern (such as
  `chrM->Mito`), warns when a name holds different sequences in the two
  genomes, and says "1 sequence", not "1 sequences".
- CLI hints on stderr print after the result they describe, even when output is
  piped.

## [1.0.0a2] - 2026-11-14

### Added

- Genome pages have Sequences, Extract, Download and Compare tabs: browse the
  sequences (with alphabet charts), extract regions or a BED file of regions,
  download all or chosen sequences as FASTA in the browser (large files in
  Firefox and Safari too), and compare with another genome or your own FASTA
  file. The Details tab finds other genomes that share each digest, and the
  Download and Extract tabs show the same steps in Python, R, Rust and JavaScript.
- Collection pages (`REFGENIE_COLLECTIONS`) show a branded slice of the catalog,
  such as the VGP assemblies, with a species tree and a genome summary pop-up.
- Genome search ranks exact alias matches first, then word starts, then the rest.
- `refgenie compare`, `GET /v4/genomes/{digest}/compare/{other}` and the MCP
  `compare_genomes` tool give one sentence on how two genomes differ and a map of
  renamed sequences; `compare --json` adds the GA4GH comparison.
- `GET /v4/genomes/{digest}/sequences`; MCP `get_sequence`; `refgenie-mcp --help`.
- `list --json`, `listr --json`, `getseq --fasta` and `getseq --zero-based`.
- `init` says what it set up; an empty `list` says what to run next; `seek`
  suggests close matches or the `pull` command.
- A server picks up stores added to its publish catalog without a restart.

### Changed

- `getseq` reads a text region as 1-based and inclusive, like samtools and the
  website (commas are fine), says on stderr how it read it, reads byte ranges
  instead of whole sequences, and registers an unregistered genome first.
  Start/end numbers (HTTP API, MCP, Python) stay 0-based.
- `refgenie pull` brings a parent asset only when the asset colocates its files:
  `refgenie pull hg38/fasta_index` downloads 12 KB, not the FASTA.
- Output is quiet by default (log level WARNING; HTTP requests show only at
  DEBUG); results and hints go to stderr. `config get` prints YAML.
- `list` and `listr` name each genome once, by its shortest plain alias.
- `listr -g` and `compare` find genomes on the subscribed servers;
  `pull --init` accepts `--genome-digest`; `genome init --server` takes only
  server URLs (use `--store` for stores).
- MCP `get_asset` lists each file and gives a `seek_command` over stdio or a
  `pull_command` on a server; `get_genome_sequences` is paged.
- Servers read each genome only from the store its catalog row names, and keep
  catalog and store-sync aliases in step, refusing changes that would drop more
  than half of them; `refgenie catalog-prune` also removes emptied genomes and
  stores.
- refgenie uses httpx2, the maintained successor of httpx.

### Fixed

- Pulling a child asset of hg38 and other genomes whose server aliases contain
  '/' no longer fails and deletes the genome.
- Pulling a genome by a second alias links all its assets under that alias, so
  `seek` no longer prints a path that does not exist.
- `genome init --server`, `alias set` and `id --remote` find genomes in every
  store, not only the first.
- `pull --init` fails clearly when a genome's sequence collection is not on the
  server (#299).
- Slow `getseq` output, whole-chromosome downloads for a few bases, truncated
  tables off a terminal, and `build` exiting silently on an unknown alias.
- `SKILL.md` is served as `text/markdown` on macOS and Python before 3.13.

### Removed

- `REFGENIE_STORES`: a server reads its stores from the publish catalog.
- `refgenie config set`, which was never implemented.
- `refgenie pull --asset-registry-paths`; pass registry paths as arguments.
- MCP `get_genome_metadata` and `get_asset`'s `downloads` and `local_path`.

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

[Unreleased]: https://github.com/refgenie/refgenie/compare/v1.0.0a3...HEAD
[1.0.0a3]: https://github.com/refgenie/refgenie/compare/v1.0.0a2...v1.0.0a3
[1.0.0a2]: https://github.com/refgenie/refgenie/compare/v1.0.0a1...v1.0.0a2
[1.0.0a1]: https://github.com/refgenie/refgenie/releases/tag/v1.0.0a1
